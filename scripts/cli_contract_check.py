#!/usr/bin/env python3
"""Verify that every flag the pipeline scripts pass to `dcs <tool>` still exists.

A flag removed from a binary while a pipeline script still passes it is invisible
until someone runs that exact code path, and the failure surfaces at a customer
site rather than in the build. This check is cheap enough to run on every commit.

It does three things:

  dispatcher  every tool advertised by `dcs` with no arguments answers --help,
              so a missing or unlaunchable libexec binary is caught immediately.
  contract    shell scripts are parsed for `dcs <tool> ...` invocations -- both
              literal flags and flags held in variables such as
              idxFlag="--use-low-memory-index" -- and each flag is looked up in
              that tool's usage text.
  coverage    flags a tool accepts that no script or case file ever passes.
              This is the parameter matrix gap list, not a failure.

Optionally pins versions: --expect-version aligner=1.1.3 fails the run if the
packaged binary reports something else, which catches a release tarball built
from stale artefacts.

Exit codes: 0 pass, 1 contract violation, 2 error.

Examples:
  cli_contract_check.py --dcs-home /path/to/dcs-tools
  cli_contract_check.py --dcs-home /path/to/dcs-tools --scan-dir cases \\
      --expect-version aligner=1.1.3 --expect-version variantCaller=2.5.0 \\
      --json reports/cli_contract.json
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lib import report as rp

# Flags eaten by the dcs wrapper before the tool ever sees them.
WRAPPER_FLAGS = {
    "Mutect2": {"--java-options"},
    "FilterMutectCalls": {"--java-options"},
}

# Tools whose --help needs a JVM; skipped automatically when java is absent.
JAVA_TOOLS = {"Mutect2", "FilterMutectCalls", "jointcaller"}

# A tool that cannot launch prints an error where its usage should be. Treating
# that as "this binary accepts no flags" would report every flag as a contract
# violation, so these are detected and reported as environment errors instead.
LAUNCH_ERRORS = (
    "error while loading shared libraries",
    "version `GLIBC",
    "not found (required by",
    "cannot open shared object file",
    "undefined symbol",
    "is not JDK",
    "JAVA_HOME",
    "command not found",
    "Permission denied",
    "Segmentation fault",
)

LONG_FLAG = re.compile(r"--[A-Za-z][A-Za-z0-9_-]*")
SHORT_FLAG = re.compile(r"(?<![\w-])-[A-Za-z](?![\w-])")
ASSIGNMENT = re.compile(r"^\s*(?:export\s+|local\s+)?([A-Za-z_][A-Za-z0-9_]*)=(.*)$")
INVOCATION = re.compile(r"(?:^|[\s;&|(`])(?:\S*/)?dcs\s+([A-Za-z0-9_][A-Za-z0-9_-]*)(.*)$")
VAR_REF = re.compile(r"\$\{?([A-Za-z_][A-Za-z0-9_]*)\}?")
VERSION = re.compile(r"[Vv]ersion[:\s]+v?(\d+\.\d+\.\d+)")
ADVERTISED = re.compile(r"^\s{2}([A-Za-z][A-Za-z0-9_-]+)\s{2,}\S")


def extract_flags(text):
    return set(LONG_FLAG.findall(text)) | set(SHORT_FLAG.findall(text))


def logical_lines(text):
    """Collapse backslash continuations so an invocation is one line."""
    return text.replace("\\\n", " ").splitlines()


def parse_script(path):
    """Return {tool: {flags}} for every `dcs <tool>` invocation in the script.

    Flags stashed in variables are resolved, because fastq2vcf.sh builds most of
    its command lines that way (idxFlag, bedParam, ploidyParam, recalTable, ...).
    """
    with open(path, errors="replace") as handle:
        text = handle.read()

    var_flags = {}
    for line in logical_lines(text):
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        match = ASSIGNMENT.match(line)
        if match:
            name, value = match.group(1), match.group(2)
            flags = extract_flags(value)
            if flags:
                var_flags.setdefault(name, set()).update(flags)

    used = {}
    invocations = 0
    for line in logical_lines(text):
        stripped = line.strip()
        if stripped.startswith("#") or stripped.startswith("echo"):
            continue
        for match in INVOCATION.finditer(line):
            tool, rest = match.group(1), match.group(2)
            if tool.startswith("-"):
                continue
            invocations += 1
            flags = extract_flags(rest)
            for var in VAR_REF.findall(rest):
                flags |= var_flags.get(var, set())
            used.setdefault(tool, set()).update(flags)
    return used, invocations


def scan_for_flags(paths, suffixes=(".sh", ".yaml", ".yml", ".json", ".txt", ".md")):
    """Union of every flag mentioned anywhere under the given paths."""
    mentioned = set()
    for root_path in paths:
        if os.path.isfile(root_path):
            candidates = [root_path]
        else:
            candidates = [
                os.path.join(dirpath, name)
                for dirpath, _dirs, names in os.walk(root_path)
                for name in names
                if name.endswith(suffixes)
            ]
        for candidate in candidates:
            try:
                with open(candidate, errors="replace") as handle:
                    mentioned |= extract_flags(handle.read())
            except OSError:
                continue
    return mentioned


def looks_like_usage(text):
    """Tell real usage text apart from a binary that failed to start."""
    if not text or len(text.strip()) < 40:
        return False, "no output"
    for pattern in LAUNCH_ERRORS:
        if pattern.lower() in text.lower():
            for line in text.splitlines():
                if pattern.lower() in line.lower():
                    return False, line.strip()[:200]
    if len(set(LONG_FLAG.findall(text))) < 2 and "sage" not in text and "ptions" not in text:
        return False, "output does not look like usage text"
    return True, None


def run_help(dcs_bin, tool, timeout, help_cache=None):
    """Return (usage_text, ok, reason). Only real usage is cached."""
    cache_file = os.path.join(help_cache, "%s.help.txt" % tool) if help_cache else None
    if cache_file and os.path.isfile(cache_file):
        with open(cache_file, errors="replace") as handle:
            return handle.read(), True, None

    best, best_reason = "", "no output"
    for argv in ([dcs_bin, tool, "--help"], [dcs_bin, tool, "-h"], [dcs_bin, tool]):
        try:
            proc = subprocess.run(
                argv,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            best, best_reason = best, "timed out after %gs" % timeout
            continue
        except OSError as exc:
            best, best_reason = best, str(exc)
            continue
        # Usage often goes to stderr and exits non-zero, so judge on content.
        text = proc.stdout or ""
        ok, reason = looks_like_usage(text)
        if ok:
            if cache_file:
                os.makedirs(help_cache, exist_ok=True)
                with open(cache_file, "w") as handle:
                    handle.write(text)
            return text, True, None
        if len(text.strip()) > len(best.strip()):
            best, best_reason = text, reason
    return best, False, best_reason


def advertised_tools(dcs_bin, timeout):
    try:
        proc = subprocess.run(
            [dcs_bin],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            timeout=timeout,
        )
    except (subprocess.TimeoutExpired, OSError):
        return []
    tools_found = []
    for line in proc.stdout.splitlines():
        match = ADVERTISED.match(line)
        if match and not line.strip().startswith("$"):
            tools_found.append(match.group(1))
    return tools_found


def check_dispatcher(report, dcs_bin, tools_list, helps, skipped):
    dead = {
        tool: helps[tool][2]
        for tool in tools_list
        if tool not in skipped and tool in helps and not helps[tool][1]
    }
    report.check(
        "dispatcher.tools",
        not dead,
        "%d advertised tool(s) launch and print usage%s"
        % (len(tools_list) - len(skipped), " (%d skipped)" % len(skipped) if skipped else ""),
        "%d advertised tool(s) failed to launch: %s" % (len(dead), ", ".join(sorted(dead))),
        dead or {"advertised": tools_list, "skipped": sorted(skipped)},
    )


def check_contracts(report, used_by_tool, helps, skipped):
    total_violations = 0
    for tool in sorted(used_by_tool):
        flags = used_by_tool[tool]
        if tool in skipped:
            report.add("contract.%s" % tool, rp.SKIP, "skipped by request")
            continue
        usage, ok, reason = helps.get(tool, ("", False, "not probed"))
        if not ok:
            # Cannot distinguish "flag removed" from "binary never started", so
            # this is an environment error rather than a contract violation.
            report.add(
                "contract.%s" % tool,
                rp.ERROR,
                "cannot verify %d flag(s): %s" % (len(flags), reason),
                {"used": sorted(flags), "output": (usage or "").strip()[:400]},
            )
            continue
        accepted = extract_flags(usage) | WRAPPER_FLAGS.get(tool, set())
        missing = sorted(flags - accepted)
        total_violations += len(missing)
        report.check(
            "contract.%s" % tool,
            not missing,
            "%d flag(s) used, all accepted" % len(flags),
            "%d flag(s) are passed but not accepted by the binary: %s"
            % (len(missing), ", ".join(missing)),
            {"used": sorted(flags), "not_accepted": missing},
        )
    return total_violations


def check_coverage(report, helps, mentioned, skipped):
    gaps = {}
    for tool, (usage, ok, _reason) in sorted(helps.items()):
        if not ok or tool in skipped:
            continue
        accepted = {f for f in extract_flags(usage) if f.startswith("--")}
        untested = sorted(accepted - mentioned - {"--help"})
        if untested:
            gaps[tool] = untested
    report.add(
        "coverage.flags",
        rp.INFO,
        "%d flag(s) across %d tool(s) are never exercised by any script or case"
        % (sum(len(v) for v in gaps.values()), len(gaps)),
        gaps,
    )


def check_versions(report, helps, expectations):
    found = {}
    for tool, (usage, ok, _reason) in sorted(helps.items()):
        if not ok:
            continue
        match = VERSION.search(usage)
        if match:
            found[tool] = match.group(1)
    if not expectations:
        report.add("versions", rp.INFO, "reported versions: %s"
                   % ", ".join("%s=%s" % item for item in sorted(found.items())), found)
        return
    mismatched = {
        tool: {"expected": expected, "reported": found.get(tool, "unknown")}
        for tool, expected in expectations.items()
        if found.get(tool) != expected
    }
    report.check(
        "versions",
        not mismatched,
        "%d pinned version(s) match the packaged binaries" % len(expectations),
        "%d version mismatch(es) -- the release tarball may contain stale binaries"
        % len(mismatched),
        mismatched or found,
    )


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--dcs-home", default=os.environ.get("DCS_HOME"),
                        help="unpacked dcs-tools directory (default: $DCS_HOME)")
    parser.add_argument("--script", action="append", default=[], metavar="PATH",
                        help="shell script to parse (default: $DCS_HOME/libexec/fastq2vcf.sh)")
    parser.add_argument("--scan-dir", action="append", default=[], metavar="PATH",
                        help="also scan here when computing flag coverage (repeatable)")
    parser.add_argument("--skip-tool", action="append", default=[], metavar="NAME",
                        help="do not invoke this tool (repeatable)")
    parser.add_argument("--help-cache", metavar="DIR",
                        help="read/write usage text here so the check can run without the bundle")
    parser.add_argument("--expect-version", action="append", default=[], metavar="TOOL=VER",
                        help="fail if the tool reports a different version (repeatable)")
    parser.add_argument("--timeout", type=float, default=60.0,
                        help="seconds to wait for each --help (default 60)")
    return rp.add_common_args(parser)


def main(argv=None):
    args = build_parser().parse_args(argv)

    report = rp.Report("cli_contract_check", dcs_home=args.dcs_home or "unset")

    if not args.dcs_home or not os.path.isdir(args.dcs_home):
        report.add("tool.dcs", rp.ERROR,
                   "set --dcs-home (or $DCS_HOME) to an unpacked dcs-tools directory")
        return report.finish(args.json, args.verbose, args.warn_is_failure)

    dcs_bin = os.path.join(args.dcs_home, "bin", "dcs")
    if not os.access(dcs_bin, os.X_OK):
        report.add("tool.dcs", rp.ERROR, "%s is missing or not executable" % dcs_bin)
        return report.finish(args.json, args.verbose, args.warn_is_failure)
    report.add("tool.dcs", rp.INFO, dcs_bin)

    scripts = args.script or [os.path.join(args.dcs_home, "libexec", "fastq2vcf.sh")]
    used_by_tool = {}
    for script in scripts:
        if not os.path.isfile(script):
            report.add("script.%s" % os.path.basename(script), rp.ERROR, "not found: %s" % script)
            continue
        used, invocations = parse_script(script)
        for tool, flags in used.items():
            used_by_tool.setdefault(tool, set()).update(flags)
        report.check(
            "script.%s" % os.path.basename(script),
            invocations > 0,
            "%d dcs invocation(s) across %d tool(s)" % (invocations, len(used)),
            "no dcs invocation found -- the parser or the script layout changed",
            {"tools": sorted(used)},
            fail_status=rp.WARN,
        )

    skipped = set(args.skip_tool)
    if not shutil.which("java") and not os.environ.get("JAVA_HOME"):
        skipped |= JAVA_TOOLS
        report.add("tool.java", rp.INFO, "java not available; skipping %s"
                   % ", ".join(sorted(JAVA_TOOLS)))

    tools_list = advertised_tools(dcs_bin, args.timeout)
    wanted = sorted(set(tools_list) | set(used_by_tool))
    helps = {
        tool: run_help(dcs_bin, tool, args.timeout, args.help_cache)
        for tool in wanted
        if tool not in skipped
    }

    if tools_list:
        check_dispatcher(report, dcs_bin, tools_list, helps, skipped)
    else:
        report.add("dispatcher.tools", rp.ERROR, "could not read the tool list from `dcs`")

    check_contracts(report, used_by_tool, helps, skipped)

    mentioned = scan_for_flags(scripts + args.scan_dir)
    check_coverage(report, helps, mentioned, skipped)

    expectations = {}
    for item in args.expect_version:
        if "=" in item:
            tool, version = item.split("=", 1)
            expectations[tool.strip()] = version.strip()
    check_versions(report, helps, expectations)

    return report.finish(args.json, args.verbose, args.warn_is_failure)


if __name__ == "__main__":
    sys.exit(main())
