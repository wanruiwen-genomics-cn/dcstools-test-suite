#!/usr/bin/env python3
"""Audit an unpacked dcs-tools bundle for packaging defects before release.

This catches the class of bug where the tarball looks complete but the binaries
inside do not actually run against the libraries shipped with them.

  rpath       every ELF file must use $ORIGIN-relative RPATH. An absolute path
              left over from the build machine is a release blocker: DT_RPATH
              outranks LD_LIBRARY_PATH, so bin/dcs cannot override it. On a
              customer machine the stale path does not exist and the loader
              falls back, so the bundle appears fine; on a build machine that
              shares the filesystem the tool silently loads somebody's
              development build instead. That makes internal release testing
              validate an artefact that is not the one being shipped.
  resolution  every NEEDED library must resolve either inside the bundle or to
              a small allowlist of genuine system libraries. Anything else is
              the same leak, observed from the loader's side.
  symbols     GLIBCXX versions required by each file against what the bundled
              libstdc++ provides, so "version `GLIBCXX_3.4.29' not found" is
              caught at package time rather than at first run.
  layout      the files bin/dcs expects are present.

Exit codes: 0 pass, 1 packaging defect, 2 error.

Example:
  bundle_audit.py --dcs-home /path/to/dcs-tools --json reports/bundle.json
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lib import report as rp
from lib import tools

# Libraries that legitimately come from the host rather than the bundle.
SYSTEM_LIBS = {
    "linux-vdso.so.1", "ld-linux-x86-64.so.2", "libc.so.6", "libm.so.6",
    "libdl.so.2", "librt.so.1", "libpthread.so.0", "libresolv.so.2",
    "libnsl.so.1", "libutil.so.1", "libgcc_s.so.1",
}

# Standard locations a host library may legitimately come from. Anything
# resolving outside both the bundle and these is a build-tree leak.
SYSTEM_DIRS = ("/lib/", "/usr/lib/", "/lib64/", "/usr/lib64/")

ELF_MAGIC = b"\x7fELF"
SYMVER = re.compile(rb"GLIBCXX_(\d+)\.(\d+)\.(\d+)")
GLIBC_SYMVER = re.compile(rb"GLIBC_(\d+)\.(\d+)(?:\.(\d+))?")
RPATH_LINE = re.compile(r"\((RPATH|RUNPATH)\).*\[(.*)\]")
NEEDED_LINE = re.compile(r"\(NEEDED\).*\[(.*)\]")
LDD_LINE = re.compile(r"^\s*(\S+)\s+=>\s+(.*?)(?:\s+\(0x[0-9a-f]+\))?$")

SCAN_DIRS = ("libexec", "lib", "bin")


def is_elf(path):
    try:
        with open(path, "rb") as handle:
            return handle.read(4) == ELF_MAGIC
    except OSError:
        return False


def iter_elf_files(root):
    for directory in SCAN_DIRS:
        base = os.path.join(root, directory)
        if not os.path.isdir(base):
            continue
        for dirpath, _dirs, names in os.walk(base):
            for name in sorted(names):
                path = os.path.join(dirpath, name)
                if os.path.isfile(path) and not os.path.islink(path) and is_elf(path):
                    yield path


def read_dynamic(readelf_bin, path):
    """Return (rpath_entries, rpath_tag, needed_sonames)."""
    try:
        proc = tools.run([readelf_bin, "-d", path], check=False)
    except OSError:
        return [], None, []
    entries, tag, needed = [], None, []
    for line in proc.stdout.splitlines():
        match = RPATH_LINE.search(line)
        if match:
            tag = match.group(1)
            entries = [e for e in match.group(2).split(":") if e]
        match = NEEDED_LINE.search(line)
        if match:
            needed.append(match.group(1))
    return entries, tag, needed


def resolve_libs(path, bundle_lib):
    """Ask the loader where each NEEDED library actually comes from."""
    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = bundle_lib + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    try:
        proc = subprocess.run(
            ["ldd", path],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            universal_newlines=True,
            env=env,
            timeout=60,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    resolved = {}
    for line in proc.stdout.splitlines():
        match = LDD_LINE.match(line.rstrip())
        if match:
            soname, target = match.group(1), match.group(2).strip()
            resolved[soname] = None if target == "not found" else target
    return resolved


def max_version(pattern, blob):
    best = None
    for match in pattern.finditer(blob):
        parts = tuple(int(g) for g in match.groups() if g is not None)
        if best is None or parts > best:
            best = parts
    return best


def fmt(version):
    return ".".join(str(p) for p in version) if version else "none"


def check_rpath(report, root, files, readelf_bin):
    offenders = {}
    scanned = 0
    for path in files:
        entries, tag, _needed = read_dynamic(readelf_bin, path)
        if not entries:
            continue
        scanned += 1
        leaked = [e for e in entries if not e.startswith("$ORIGIN")]
        if leaked:
            offenders[os.path.relpath(path, root)] = {"tag": tag, "absolute": leaked}
    report.check(
        "rpath.absolute",
        not offenders,
        "%d file(s) with an RPATH, all $ORIGIN-relative" % scanned,
        "%d file(s) carry build-machine paths in %s -- these override "
        "LD_LIBRARY_PATH and make the bundle load libraries from outside itself"
        % (len(offenders), "/".join(sorted({v["tag"] for v in offenders.values()}))),
        offenders,
    )
    return offenders


def check_resolution(report, root, files, allow_outside):
    bundle_lib = os.path.join(root, "lib")
    real_root = os.path.realpath(root)
    leaked, host_provided, missing, unchecked = {}, {}, {}, []
    for path in files:
        resolved = resolve_libs(path, bundle_lib)
        if resolved is None:
            unchecked.append(os.path.relpath(path, root))
            continue
        for soname, target in sorted(resolved.items()):
            if soname in SYSTEM_LIBS or soname in allow_outside:
                continue
            relative = os.path.relpath(path, root)
            if target is None:
                missing.setdefault(relative, []).append(soname)
            elif os.path.realpath(target).startswith(real_root):
                continue
            elif target.startswith(SYSTEM_DIRS):
                # Not shipped, but supplied by any normal host. Worth listing as
                # the bundle's external dependency set, not a defect.
                host_provided.setdefault(soname, set()).add(os.path.dirname(target))
            else:
                leaked.setdefault(relative, []).append("%s => %s" % (soname, target))

    report.check(
        "libs.buildtree",
        not leaked,
        "no binary loads a library from a build tree",
        "%d file(s) load libraries from a developer build tree instead of the "
        "bundle -- internal testing would not be testing the released artefact"
        % len(leaked),
        leaked,
    )
    report.add(
        "libs.host",
        rp.INFO,
        "%d librar(ies) come from the host: %s"
        % (len(host_provided), ", ".join(sorted(host_provided))),
        {name: sorted(dirs) for name, dirs in host_provided.items()},
    )
    report.check(
        "libs.missing",
        not missing,
        "no unresolved libraries",
        "%d file(s) have unresolved libraries" % len(missing),
        missing,
    )
    if unchecked:
        report.add("libs.unchecked", rp.INFO,
                   "%d file(s) could not be probed with ldd" % len(unchecked),
                   {"files": unchecked})


def check_symbols(report, root, files, readelf_bin):
    libstdcxx = os.path.join(root, "lib", "libstdc++.so.6")
    if not os.path.isfile(libstdcxx):
        report.add("symbols.glibcxx", rp.SKIP, "no libstdc++ bundled")
        return
    with open(libstdcxx, "rb") as handle:
        provided = max_version(SYMVER, handle.read())

    offenders = {}
    glibc_required = None
    for path in files:
        _entries, _tag, needed = read_dynamic(readelf_bin, path)
        with open(path, "rb") as handle:
            blob = handle.read()
        glibc = max_version(GLIBC_SYMVER, blob)
        if glibc and (glibc_required is None or glibc > glibc_required):
            glibc_required = glibc
        if "libstdc++.so.6" not in needed or os.path.realpath(path) == os.path.realpath(libstdcxx):
            continue
        required = max_version(SYMVER, blob)
        if required and provided and required > provided:
            offenders[os.path.relpath(path, root)] = {
                "requires": "GLIBCXX_" + fmt(required),
                "bundle provides": "GLIBCXX_" + fmt(provided),
            }

    report.check(
        "symbols.glibcxx",
        not offenders,
        "all C++ binaries fit the bundled libstdc++ (GLIBCXX_%s)" % fmt(provided),
        "%d file(s) need a newer libstdc++ than the bundle ships" % len(offenders),
        offenders or {"bundle provides": "GLIBCXX_" + fmt(provided)},
    )
    report.add("symbols.glibc", rp.INFO,
               "bundle needs glibc >= %s on the target host" % fmt(glibc_required))


def check_layout(report, root):
    missing = [
        rel for rel in ("bin/dcs", "libexec")
        if not os.path.exists(os.path.join(root, rel))
    ]
    has_lib = os.path.isdir(os.path.join(root, "lib"))
    has_tarball = os.path.isfile(os.path.join(root, "lib.tar.gz"))
    if not has_lib and not has_tarball:
        missing.append("lib/ or lib.tar.gz")
    report.check(
        "layout",
        not missing,
        "bundle layout complete",
        "missing from the bundle: %s" % ", ".join(missing),
        {"lib_dir": has_lib, "lib_tarball": has_tarball},
    )


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--dcs-home", default=os.environ.get("DCS_HOME"),
                        help="unpacked dcs-tools directory (default: $DCS_HOME)")
    parser.add_argument("--allow-outside", action="append", default=[], metavar="SONAME",
                        help="extra library allowed to come from the host (repeatable)")
    parser.add_argument("--readelf-bin", help="path to readelf (default: PATH)")
    parser.add_argument("--skip-ldd", action="store_true",
                        help="skip the loader resolution check")
    return rp.add_common_args(parser)


def main(argv=None):
    args = build_parser().parse_args(argv)
    report = rp.Report("bundle_audit", dcs_home=args.dcs_home or "unset")

    if not args.dcs_home or not os.path.isdir(args.dcs_home):
        report.add("layout", rp.ERROR,
                   "set --dcs-home (or $DCS_HOME) to an unpacked dcs-tools directory")
        return report.finish(args.json, args.verbose, args.warn_is_failure)

    root = os.path.abspath(args.dcs_home)
    try:
        readelf_bin = tools.resolve("readelf", "READELF", args.readelf_bin)
    except tools.ToolNotFound as exc:
        report.add("tool.readelf", rp.ERROR, str(exc))
        return report.finish(args.json, args.verbose, args.warn_is_failure)

    check_layout(report, root)
    files = list(iter_elf_files(root))
    if not files:
        report.add("scan", rp.ERROR, "no ELF files found under %s" % ", ".join(SCAN_DIRS))
        return report.finish(args.json, args.verbose, args.warn_is_failure)
    report.add("scan", rp.INFO, "%d ELF file(s) under %s" % (len(files), ", ".join(SCAN_DIRS)))

    try:
        check_rpath(report, root, files, readelf_bin)
        if not args.skip_ldd:
            check_resolution(report, root, files, set(args.allow_outside))
        check_symbols(report, root, files, readelf_bin)
    except KeyboardInterrupt:
        report.add("run", rp.ERROR, "interrupted")

    return report.finish(args.json, args.verbose, args.warn_is_failure)


if __name__ == "__main__":
    sys.exit(main())
