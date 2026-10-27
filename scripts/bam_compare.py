#!/usr/bin/env python3
"""Compare a test BAM against a golden BAM with timestamp noise normalised away.

Why not plain md5: a BAM carries @PG lines with the full command line and the
tool version, so two byte-identical runs of different builds never match. This
script strips that noise and compares what actually matters.

What it does:

  header    @SQ (reference contigs), @HD sort order and @RG are compared
            exactly; @PG is skipped unless --compare-pg is given.
  records   every alignment is canonicalised (optional tags sorted, noisy tags
            dropped) and folded into two digests -- one order-sensitive, one
            order-insensitive. Matching unordered digests with differing
            ordered digests means the same reads came out in a different order,
            which is what multi-threaded sort/markdup instability looks like.
  stats     flagstat / idxstats / samtools stats metrics compared against the
            tolerances in thresholds, so a run that is not bit-identical can
            still be judged on mapped rate, duplicate rate, insert size, etc.

Exit codes: 0 pass, 1 fail, 2 error.

Examples:
  bam_compare.py -b golden/wgs_human_small.bam -t out/sample.align.bam
  bam_compare.py -b a.bam -t b.bam --region chr20 --json reports/bam.json
  bam_compare.py -b a.bam -t b.bam --skip-stats --max-records 2000000
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import threading
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lib import report as rp
from lib import tools
from lib.digest import StreamDigest

# Expose the status constants at module level so the lock-aware wrappers in
# compare_stats can use them without a global lookup.
PASS = rp.PASS
WARN = rp.WARN
FAIL = rp.FAIL
ERROR = rp.ERROR
SKIP = rp.SKIP
INFO = rp.INFO

DEFAULTS = {
    # Record-level tags dropped before hashing. PG points at an @PG ID that
    # changes with the build; it says nothing about the alignment itself.
    "ignore_tags": ["PG"],
    # samtools flagstat / idxstats / stats: exact match (no relative tolerance).
    # Diffs are reported with full baseline/test tables in the JSON detail.
    "flagstat_rel_tol": 0.0,
    "idxstats_rel_tol": 0.0,
    "rate_abs_tol": 0.0,
    "insert_size_rel_tol": 0.0,
    "error_rate_rel_tol": 0.0,
    "file_size_rel_tol": 0.05,
    # Aligner output is expected to ship an index; turn this off with
    # --no-require-index when comparing intermediate BAMs.
    "require_index": True,
    # How many differing records to print when the ordered digests disagree.
    "max_reported_diffs": 5,
    # Auto-sampling: when a BAM exceeds auto_sample_threshold_gb, automatically
    # limit the record digest to auto_sample_records reads.  A 30x WGS BAM
    # carries ~1 billion reads; Python-side canonicalisation runs at ~150k
    # reads/s, so a full digest takes ~2 hours.  Sampling 5M reads still
    # covers every aligner code path and finishes in ~30 s.  The flagstat /
    # idxstats / stats checks still run over the whole file, so a regression
    # that shifts mapping rates or insert sizes is still caught.
    "auto_sample_threshold_gb": 10,
    "auto_sample_records": 5_000_000,
}

SN_KEYS = [
    "raw total sequences",
    "reads mapped",
    "reads mapped and paired",
    "reads properly paired",
    "reads duplicated",
    "reads MQ0",
    "reads unmapped",
    "average quality",
    "insert size average",
    "insert size standard deviation",
    "error rate",
    "mismatches",
    "bases mapped (cigar)",
]

RATE_KEYS = {
    "mapped_rate": ("reads mapped", "raw total sequences"),
    "properly_paired_rate": ("reads properly paired", "raw total sequences"),
    "duplicate_rate": ("reads duplicated", "raw total sequences"),
    "mq0_rate": ("reads MQ0", "raw total sequences"),
}

FLOAT_SN_KEYS = {
    "average quality",
    "insert size average",
    "insert size standard deviation",
    "error rate",
}


# --------------------------------------------------------------------------- header


def read_header(samtools_bin, path):
    proc = tools.run([samtools_bin, "view", "-H", path])
    groups = {}
    for line in proc.stdout.splitlines():
        if not line.startswith("@"):
            continue
        groups.setdefault(line[:3], []).append(line)
    return groups


def parse_header_line(line):
    fields = line.rstrip("\n").split("\t")
    parsed = {}
    for field in fields[1:]:
        if ":" in field:
            key, value = field.split(":", 1)
            parsed[key] = value
    return parsed


def compare_sq(report, base_header, test_header):
    base = [(p.get("SN"), p.get("LN")) for p in map(parse_header_line, base_header.get("@SQ", []))]
    test = [(p.get("SN"), p.get("LN")) for p in map(parse_header_line, test_header.get("@SQ", []))]
    if base == test:
        report.add(
            "header.sq",
            rp.PASS,
            "reference identical: %d contigs" % len(base),
            {"contigs": len(base)},
        )
        return True

    base_names = [name for name, _ in base]
    test_names = [name for name, _ in test]
    only_base = [n for n in base_names if n not in set(test_names)]
    only_test = [n for n in test_names if n not in set(base_names)]
    length_diffs = [
        "%s: %s vs %s" % (name, base_len, dict(test).get(name))
        for name, base_len in base
        if name in dict(test) and dict(test)[name] != base_len
    ]
    report.add(
        "header.sq",
        rp.FAIL,
        "reference mismatch -- the two BAMs were aligned against different references",
        {
            "baseline_contigs": len(base),
            "test_contigs": len(test),
            "only_in_baseline": only_base[:10],
            "only_in_test": only_test[:10],
            "length_mismatch": length_diffs[:10],
            "order_changed": bool(not only_base and not only_test and base_names != test_names),
        },
    )
    return False


def compare_hd(report, base_header, test_header):
    def sort_order(header):
        lines = header.get("@HD", [])
        return parse_header_line(lines[0]).get("SO", "unknown") if lines else "absent"

    base, test = sort_order(base_header), sort_order(test_header)
    report.check(
        "header.hd",
        base == test,
        "sort order: %s" % base,
        "sort order differs: baseline=%s test=%s" % (base, test),
        {"baseline": base, "test": test},
    )


def compare_rg(report, base_header, test_header):
    base = sorted(base_header.get("@RG", []))
    test = sorted(test_header.get("@RG", []))
    if base == test:
        report.add(
            "header.rg",
            rp.PASS,
            "%d read group(s) identical" % len(base) if base
            else "no read groups in either file",
        )
        return
    report.add(
        "header.rg",
        rp.WARN,
        "read group differs (check the --read-group passed to the aligner)",
        {"only_in_baseline": [l for l in base if l not in test][:5],
         "only_in_test": [l for l in test if l not in base][:5]},
    )


def compare_pg(report, base_header, test_header):
    def strip(lines):
        # CL holds absolute paths and VN the build version; neither is a result.
        return sorted(
            "\t".join(f for f in line.split("\t") if not f.startswith(("CL:", "VN:")))
            for line in lines
        )

    base, test = strip(base_header.get("@PG", [])), strip(test_header.get("@PG", []))
    report.check(
        "header.pg",
        base == test,
        "@PG chain identical (ignoring CL/VN)",
        "@PG chain differs -- the pipeline steps are not the same",
        {"baseline": base[:6], "test": test[:6]},
    )


# --------------------------------------------------------------------------- records


def make_canonicaliser(ignore_tags):
    ignored = set(ignore_tags)

    def canonicalise(line):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            return line.rstrip("\n")
        optional = [tag for tag in fields[11:] if tag[:2] not in ignored]
        optional.sort()
        return "\t".join(fields[:11] + optional)

    return canonicalise


def view_cmd(samtools_bin, path, region, threads):
    cmd = [samtools_bin, "view"]
    if threads > 1:
        cmd += ["-@", str(threads)]
    cmd.append(path)
    if region:
        cmd.append(region)
    return cmd


def digest_bam(samtools_bin, path, region, threads, canonicalise, max_records):
    digest = StreamDigest()
    for line in tools.stream_lines(view_cmd(samtools_bin, path, region, threads), max_records):
        digest.update(canonicalise(line))
    return digest


def first_differences(samtools_bin, base, test, region, threads, canonicalise, limit, max_records):
    """Walk both record streams in lockstep and collect the first differences."""
    base_iter = tools.stream_lines(view_cmd(samtools_bin, base, region, threads), max_records)
    test_iter = tools.stream_lines(view_cmd(samtools_bin, test, region, threads), max_records)
    diffs = []
    index = 0
    while len(diffs) < limit:
        base_line = next(base_iter, None)
        test_line = next(test_iter, None)
        if base_line is None and test_line is None:
            break
        index += 1
        left = canonicalise(base_line) if base_line is not None else "<end of file>"
        right = canonicalise(test_line) if test_line is not None else "<end of file>"
        if left != right:
            diffs.append({"record": index, "baseline": left[:300], "test": right[:300]})
    base_iter.close()
    test_iter.close()
    return diffs


def compare_records(report, args, config, samtools_bin, lock=None):
    def _add(*a, **kw):
        if lock:
            with lock:
                return report.add(*a, **kw)
        return report.add(*a, **kw)

    canonicalise = make_canonicaliser(config["ignore_tags"])
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [
            pool.submit(
                digest_bam,
                samtools_bin,
                path,
                args.region,
                args.threads,
                canonicalise,
                args.max_records,
            )
            for path in (args.baseline, args.test)
        ]
        base_digest, test_digest = (future.result() for future in futures)

    detail = {"baseline": base_digest.as_dict(), "test": test_digest.as_dict()}
    if args.max_records:
        detail["note"] = "limited to the first %d records" % args.max_records

    if base_digest.count == 0 and test_digest.count == 0:
        # Two empty streams hash identically. Passing here would turn a typo in
        # --region, or an empty BAM, into a green run.
        _add(
            "records.digest",
            rp.ERROR,
            "no records on either side -- nothing was compared"
            + (" (does --region %s match a contig in the header?)" % args.region
               if args.region else ""),
            detail,
        )
        return

    if base_digest.ordered == test_digest.ordered:
        _add(
            "records.digest",
            rp.PASS,
            "%d records identical (same content, same order)" % base_digest.count,
            detail,
        )
        return

    if base_digest.unordered == test_digest.unordered and base_digest.count == test_digest.count:
        _add(
            "records.digest",
            rp.WARN,
            "same %d records but emitted in a different order -- suspect non-deterministic "
            "sort/markdup; rerun with --threads 1 on both sides to confirm" % base_digest.count,
            detail,
        )
        return

    detail["record_count_delta"] = test_digest.count - base_digest.count
    try:
        detail["first_differences"] = first_differences(
            samtools_bin,
            args.baseline,
            args.test,
            args.region,
            args.threads,
            canonicalise,
            config["max_reported_diffs"],
            args.max_records,
        )
    except tools.CommandFailed as exc:
        detail["first_differences"] = "unavailable: %s" % exc

    _add(
        "records.digest",
        rp.FAIL,
        "alignment records differ (baseline %d records, test %d)"
        % (base_digest.count, test_digest.count),
        detail,
    )


# --------------------------------------------------------------------------- stats


def parse_flagstat(samtools_bin, path):
    try:
        proc = tools.run([samtools_bin, "flagstat", "-O", "tsv", path])
        metrics = {}
        for line in proc.stdout.splitlines():
            parts = line.split("\t")
            if len(parts) >= 3:
                metrics[parts[2].strip()] = int(parts[0])
        if metrics:
            return metrics
    except tools.CommandFailed:
        pass
    # Older samtools without -O tsv.
    proc = tools.run([samtools_bin, "flagstat", path])
    metrics = {}
    for line in proc.stdout.splitlines():
        match = re.match(r"^(\d+) \+ (\d+) ([^(]+)", line)
        if match:
            metrics[match.group(3).strip()] = int(match.group(1))
    return metrics


def parse_idxstats(samtools_bin, path):
    proc = tools.run([samtools_bin, "idxstats", path])
    metrics = {}
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 4:
            metrics[parts[0]] = (int(parts[2]), int(parts[3]))
    return metrics


def parse_stats(samtools_bin, path, threads, region):
    cmd = [samtools_bin, "stats"]
    if threads > 1:
        cmd += ["-@", str(threads)]
    cmd.append(path)
    if region:
        cmd.append(region)
    proc = tools.run(cmd)
    metrics = {}
    for line in proc.stdout.splitlines():
        if not line.startswith("SN\t"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        key = parts[1].rstrip(":")
        if key not in SN_KEYS:
            continue
        raw = parts[2].strip()
        try:
            metrics[key] = float(raw) if key in FLOAT_SN_KEYS else int(raw)
        except ValueError:
            continue
    return metrics


def compare_metrics(report, name, base, test, tolerance, absolute=False, label=""):
    """Compare two flat metric dicts, reporting every out-of-tolerance key."""
    if not base and not test:
        report.add(name, rp.SKIP, "no %s metrics available" % (label or name))
        return
    offenders = {}
    for key in sorted(set(base) | set(test)):
        left, right = base.get(key), test.get(key)
        if left is None or right is None:
            offenders[key] = {"baseline": left, "test": right, "note": "missing on one side"}
            continue
        delta = abs(left - right) if absolute else rp.rel_diff(left, right)
        if delta > tolerance:
            offenders[key] = {
                "baseline": left,
                "test": right,
                ("abs_diff" if absolute else "rel_diff"): round(delta, 6),
            }
    # Always keep baseline/test values so the HTML report can render a
    # side-by-side metrics table even when every key is within tolerance.
    detail = {
        "tolerance": tolerance,
        "compared": len(set(base) | set(test)),
        "baseline": base,
        "test": test,
    }
    if offenders:
        detail["offenders"] = offenders
    report.check(
        name,
        not offenders,
        "%d %s metrics within tolerance (%s%s)"
        % (len(set(base) | set(test)), label or name, tolerance, "" if absolute else " relative"),
        "%d %s metric(s) out of tolerance" % (len(offenders), label or name),
        detail,
    )


def derive_rates(metrics):
    rates = {}
    for name, (numerator, denominator) in RATE_KEYS.items():
        total = metrics.get(denominator)
        value = metrics.get(numerator)
        if total:
            rates[name] = value / total
    return rates


def compare_stats(report, args, config, samtools_bin, lock=None):
    def _add(name, status, message, detail=None):
        if lock:
            with lock:
                return report.add(name, status, message, detail)
        return report.add(name, status, message, detail)

    def _check(name, ok, ok_message, fail_message, detail=None, fail_status=FAIL):
        if lock:
            with lock:
                return report.check(name, ok, ok_message, fail_message, detail, fail_status)
        return report.check(name, ok, ok_message, fail_message, detail, fail_status)

    # Wrap report with lock-aware versions for the metric comparison helpers.
    class _LockedReport:
        def add(self, name, status, message, detail=None):
            return _add(name, status, message, detail)
        def check(self, name, ok, ok_message, fail_message, detail=None, fail_status=FAIL):
            return _check(name, ok, ok_message, fail_message, detail, fail_status)

    locked_report = _LockedReport()
    threads = max(1, int(args.threads or 16))

    # Run baseline and test samtools invocations in parallel (2 workers), each
    # with -@ threads. Wall time is roughly one side of the pair, not the sum.
    def _flagstat_pair():
        with ThreadPoolExecutor(max_workers=2) as pool:
            fb = pool.submit(parse_flagstat, samtools_bin, args.baseline)
            ft = pool.submit(parse_flagstat, samtools_bin, args.test)
            return fb.result(), ft.result()

    def _idxstats_pair():
        with ThreadPoolExecutor(max_workers=2) as pool:
            fb = pool.submit(parse_idxstats, samtools_bin, args.baseline)
            ft = pool.submit(parse_idxstats, samtools_bin, args.test)
            return fb.result(), ft.result()

    def _stats_pair():
        with ThreadPoolExecutor(max_workers=2) as pool:
            fb = pool.submit(parse_stats, samtools_bin, args.baseline, threads, args.region)
            ft = pool.submit(parse_stats, samtools_bin, args.test, threads, args.region)
            return fb.result(), ft.result()

    scope = " (whole file; --region only applies to samtools stats)" if args.region else ""
    try:
        base_flag, test_flag = _flagstat_pair()
    except tools.CommandFailed as exc:
        _add("stats.flagstat", rp.ERROR, str(exc))
    else:
        compare_metrics(
            locked_report,
            "stats.flagstat",
            base_flag,
            test_flag,
            config["flagstat_rel_tol"],
            label="flagstat" + scope,
        )

    if has_index(args.baseline) and has_index(args.test):
        try:
            base_idx, test_idx = _idxstats_pair()
        except tools.CommandFailed as exc:
            _add("stats.idxstats", rp.ERROR, str(exc))
        else:
            base_flat = {"%s:%s" % (k, kind): v[i]
                         for k, v in base_idx.items()
                         for i, kind in enumerate(("mapped", "unmapped"))}
            test_flat = {"%s:%s" % (k, kind): v[i]
                         for k, v in test_idx.items()
                         for i, kind in enumerate(("mapped", "unmapped"))}
            compare_metrics(
                locked_report,
                "stats.idxstats",
                base_flat,
                test_flat,
                config["idxstats_rel_tol"],
                label="per-contig",
            )
    else:
        _add("stats.idxstats", rp.SKIP, "needs a .bai/.csi index on both files")

    try:
        base_sn, test_sn = _stats_pair()
    except tools.CommandFailed as exc:
        _add("stats.summary", rp.ERROR, str(exc))
        return

    insert_keys = {"insert size average", "insert size standard deviation"}
    compare_metrics(
        locked_report,
        "stats.insert_size",
        {k: v for k, v in base_sn.items() if k in insert_keys},
        {k: v for k, v in test_sn.items() if k in insert_keys},
        config["insert_size_rel_tol"],
        label="insert size",
    )
    compare_metrics(
        locked_report,
        "stats.error_rate",
        {k: v for k, v in base_sn.items() if k in ("error rate", "average quality", "mismatches")},
        {k: v for k, v in test_sn.items() if k in ("error rate", "average quality", "mismatches")},
        config["error_rate_rel_tol"],
        label="error rate",
    )
    compare_metrics(
        locked_report,
        "stats.rates",
        derive_rates(base_sn),
        derive_rates(test_sn),
        config["rate_abs_tol"],
        absolute=True,
        label="mapping rate",
    )


# --------------------------------------------------------------------------- files


def has_index(path):
    return any(
        os.path.isfile(candidate)
        for candidate in (path + ".bai", path + ".csi", os.path.splitext(path)[0] + ".bai")
    )


def check_files(report, args, config):
    missing = [p for p in (args.baseline, args.test) if not os.path.isfile(p)]
    if missing:
        report.add("file.exists", rp.ERROR, "missing input: %s" % ", ".join(missing))
        return False
    report.add("file.exists", rp.PASS, "both inputs present")

    if config["require_index"]:
        report.check(
            "file.index",
            has_index(args.test),
            "index present next to the test BAM",
            "no .bai/.csi next to the test BAM -- the aligner should emit one",
            {"test": args.test},
        )
    else:
        report.add("file.index", rp.INFO,
                   "present" if has_index(args.test) else "absent (not required)")

    base_size = os.path.getsize(args.baseline)
    test_size = os.path.getsize(args.test)
    delta = rp.rel_diff(base_size, test_size)
    report.check(
        "file.size",
        delta <= config["file_size_rel_tol"],
        "size within %.1f%% (baseline %d, test %d)" % (config["file_size_rel_tol"] * 100, base_size, test_size),
        "size differs by %.2f%% (baseline %d, test %d)" % (delta * 100, base_size, test_size),
        {"baseline_bytes": base_size, "test_bytes": test_size, "rel_diff": round(delta, 6)},
        fail_status=rp.WARN,
    )
    return True


# --------------------------------------------------------------------------- main


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-b", "--baseline", required=True, help="golden BAM")
    parser.add_argument("-t", "--test", required=True, help="BAM produced by the build under test")
    parser.add_argument("--samtools-bin", help="path to samtools (default: $SAMTOOLS, env/tools.conf, PATH)")
    parser.add_argument("--region", help="restrict the record comparison, e.g. chr20:1-10000000 (needs an index)")
    parser.add_argument("--max-records", type=int, help="stop after N records (smoke tier)")
    parser.add_argument("--ignore-tag", action="append", default=[], metavar="XX",
                        help="drop this optional tag before hashing (repeatable)")
    parser.add_argument("--compare-pg", action="store_true", help="also compare the @PG chain")
    parser.add_argument("--no-require-index", action="store_true",
                        help="do not fail when the test BAM has no .bai/.csi")
    parser.add_argument("--skip-stats", action="store_true", help="digest only, skip flagstat/idxstats/stats")
    parser.add_argument("--threads", type=int, default=16,
                        help="threads handed to samtools stats -@; baseline and "
                             "test also run in parallel (default 16)")
    parser.add_argument("--auto-sample-threshold", type=float, default=None, metavar="GB",
                        help="auto-cap the digest to --auto-sample-records when a BAM exceeds this size (default: %(default)s)")
    parser.add_argument("--auto-sample-records", type=int, default=None, metavar="N",
                        help="record cap used by auto-sampling (default: %(default)s)")
    parser.add_argument("--no-auto-sample", action="store_true",
                        help="disable auto-sampling even for very large BAMs")
    return rp.add_common_args(parser)


def main(argv=None):
    args = build_parser().parse_args(argv)
    config = rp.load_config(DEFAULTS, args.config)
    if args.ignore_tag:
        config["ignore_tags"] = sorted(set(config["ignore_tags"]) | set(args.ignore_tag))
    if args.no_require_index:
        config["require_index"] = False
    if args.auto_sample_threshold is not None:
        config["auto_sample_threshold_gb"] = args.auto_sample_threshold
    if args.auto_sample_records is not None:
        config["auto_sample_records"] = args.auto_sample_records
    if args.no_auto_sample:
        config["auto_sample_threshold_gb"] = float("inf")

    report = rp.Report(
        "bam_compare",
        baseline=args.baseline,
        test=args.test,
        region=args.region or "whole file",
        ignored=",".join(config["ignore_tags"]) or "none",
    )

    try:
        samtools_bin = tools.samtools(args.samtools_bin)
    except tools.ToolNotFound as exc:
        report.add("tool.samtools", rp.ERROR, str(exc))
        return report.finish(args.json, args.verbose, args.warn_is_failure)
    report.add("tool.samtools", rp.INFO, samtools_bin)

    try:
        if not check_files(report, args, config):
            return report.finish(args.json, args.verbose, args.warn_is_failure)

        base_header = read_header(samtools_bin, args.baseline)
        test_header = read_header(samtools_bin, args.test)
        same_reference = compare_sq(report, base_header, test_header)
        compare_hd(report, base_header, test_header)
        compare_rg(report, base_header, test_header)
        if args.compare_pg:
            compare_pg(report, base_header, test_header)

        if not same_reference:
            report.add("records.digest", rp.SKIP, "skipped: the two BAMs use different references")
        else:
            # Auto-sample large BAMs: if neither --max-records nor --region was
            # given and the BAM exceeds the threshold, cap the digest.
            if not args.max_records and not args.region:
                base_gb = os.path.getsize(args.baseline) / (1024 ** 3)
                test_gb = os.path.getsize(args.test) / (1024 ** 3)
                threshold = config["auto_sample_threshold_gb"]
                if max(base_gb, test_gb) >= threshold:
                    args.max_records = config["auto_sample_records"]
                    report.add("records.sampling", rp.INFO,
                                "auto-sampling %d records (BAM %.1f GB >= %d GB threshold); "
                                "flagstat/idxstats/stats still cover the whole file"
                                % (args.max_records, max(base_gb, test_gb), threshold))

            # Run the record digest and the stats checks concurrently.  The
            # digest is CPU-bound in Python (GIL), but it spends most of its
            # wall time in samtools subprocess I/O, which releases the GIL;
            # the stats checks are pure subprocess I/O.  A lock serialises
            # writes to the shared report.
            report_lock = threading.Lock()
            if not args.skip_stats:
                def run_stats():
                    compare_stats(report, args, config, samtools_bin, report_lock)
                def run_digest():
                    compare_records(report, args, config, samtools_bin, report_lock)
                with ThreadPoolExecutor(max_workers=2) as pool:
                    f_stats = pool.submit(run_stats)
                    f_digest = pool.submit(run_digest)
                    f_stats.result()
                    f_digest.result()
            else:
                compare_records(report, args, config, samtools_bin, None)
    except tools.CommandFailed as exc:
        report.add("run", rp.ERROR, str(exc))
    except KeyboardInterrupt:
        report.add("run", rp.ERROR, "interrupted")

    return report.finish(args.json, args.verbose, args.warn_is_failure)


if __name__ == "__main__":
    sys.exit(main())
