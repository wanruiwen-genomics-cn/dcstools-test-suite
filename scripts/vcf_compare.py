#!/usr/bin/env python3
"""Compare a test VCF/gVCF against a baseline, site by site and genotype by genotype.

Two modes:

  --mode regression  (default) the baseline is a golden file produced by a
                     known-good dcstools build. Any difference is suspicious,
                     so the thresholds default to "must be identical".
  --mode truth       the baseline is an external truth set (GIAB NA12878/HG002)
                     or another caller such as GATK. The body digest is then
                     informational and the verdict comes from
                     precision/recall/F1 against the configured minimums.

What it does:

  header       contig list, sample names and INFO/FORMAT definitions. Noise
               lines (##fileDate, ##source, ##*Command*) are ignored.
  records      full-body digest of every VCF record (INFO + FORMAT + sample
               columns kept verbatim). Regression mode requires byte-identical
               records (ordered digest match).
  concordance  streaming merge join on (contig, pos, ref, alt) giving shared /
               baseline-only / test-only counts split by variant type, plus
               genotype concordance on the shared sites. Constant memory, so it
               handles a WGS VCF.
  stats        bcftools stats summary metrics (Ti/Tv, counts per type).

Both inputs must be coordinate sorted with the same contig order; that is what
the header check verifies before the join runs.

Exit codes: 0 pass, 1 fail, 2 error.

Examples:
  vcf_compare.py -b golden/wgs_small.vcf.gz -t out/sample.genotyper.vcf.gz
  vcf_compare.py -b giab_hg002.vcf.gz -t out/sample.vcf.gz --mode truth \\
                 --region chr20 --norm-ref /ref/hg38.fa --json reports/vcf.json
"""

from __future__ import annotations

import argparse
import collections
import hashlib
import os
import shutil
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from lib import report as rp
from lib import tools
from lib.digest import StreamDigest

DEFAULTS = {
    # Strict regression: keep every INFO/FORMAT field and compare the full VCF
    # body as emitted by bcftools view (byte-identical after a light
    # normalisation that only sorts INFO key order for stable hashing -- set
    # ignore_info empty and sort_info False for true column identity).
    "ignore_info": [],
    "ignore_qual": False,
    "sort_info": False,
    # Exact match: any record-count or site drift fails.
    "count_rel_tol": 0.0,
    "count_abs_tol": 0,
    "site_jaccard_min": 1.0,
    "gt_concordance_min": 1.0,
    "tstv_rel_tol": 0.0,
    "precision_min": None,
    "recall_min": None,
    "f1_min": None,
    # How many differing records to embed in records.digest detail for the
    # HTML report (locus + truncated baseline/test lines).
    "max_reported_diffs": 20,
    # Somatic only. 0 keeps the WGS rule (any INFO/FORMAT byte difference
    # fails). Above 0, a record may differ in QUAL/INFO/FORMAT counts and
    # still pass when chrom, pos, ref, alt, FILTER and GT match, as long as
    # the fraction of such records stays within this relative tolerance.
    "annotation_rel_tol": 0.0,
}

TRUTH_DEFAULTS = {
    "count_rel_tol": 0.05,
    "site_jaccard_min": None,
    "gt_concordance_min": 0.99,
    "tstv_rel_tol": 0.05,
    "precision_min": 0.99,
    "recall_min": 0.99,
    "f1_min": 0.99,
}

NOISE_PREFIXES = ("##fileDate", "##source", "##reference", "##bcftools", "##GATKCommandLine",
                  "##DCSCommandLine", "##filedate")

QUERY_FORMAT = "%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%FILTER[\\t%GT]\\n"

STATS_KEYS = ("number of SNPs:", "number of indels:", "number of MNPs:",
              "number of others:", "number of multiallelic sites:", "number of records:")


# --------------------------------------------------------------------------- model

Record = collections.namedtuple("Record", "key contig pos ref alt qual filt gts vtype")


def classify(ref, alt):
    alts = [a for a in alt.split(",") if a not in (".", "*")]
    if not alts:
        return "other"
    if any(a.startswith("<") or "[" in a or "]" in a for a in alts):
        return "symbolic"
    if all(len(a) == len(ref) for a in alts):
        return "snp" if len(ref) == 1 else "mnp"
    return "indel"


def normalise_gt(raw, strict_phase=False):
    if strict_phase:
        return raw
    if "|" not in raw and "/" not in raw:
        return raw
    separator = "|" if "|" in raw else "/"
    alleles = raw.split(separator)
    return "/".join(sorted(alleles))


# --------------------------------------------------------------------------- header


def read_header(bcftools_bin, path):
    return tools.run([bcftools_bin, "view", "-h", path]).stdout.splitlines()


def header_contigs(lines):
    contigs = []
    for line in lines:
        if line.startswith("##contig=<"):
            body = line[len("##contig=<"):].rstrip(">")
            fields = dict(
                part.split("=", 1) for part in body.split(",") if "=" in part
            )
            contigs.append((fields.get("ID"), fields.get("length")))
    return contigs


def header_samples(lines):
    for line in lines:
        if line.startswith("#CHROM"):
            return line.split("\t")[9:]
    return []


def header_definitions(lines):
    defs = set()
    for line in lines:
        for kind in ("##INFO=<ID=", "##FORMAT=<ID=", "##FILTER=<ID="):
            if line.startswith(kind):
                defs.add(kind.strip("#<=") + ":" + line[len(kind):].split(",", 1)[0].rstrip(">"))
    return defs


def compare_header(report, base_lines, test_lines):
    base_contigs = header_contigs(base_lines)
    test_contigs = header_contigs(test_lines)
    same_contigs = base_contigs == test_contigs
    if same_contigs:
        report.add("header.contigs", rp.PASS, "%d contigs identical and in the same order"
                   % len(base_contigs), {"contigs": len(base_contigs)})
    else:
        base_names = [c[0] for c in base_contigs]
        test_names = [c[0] for c in test_contigs]
        report.add(
            "header.contigs",
            rp.FAIL,
            "contig list differs -- different reference, or a different contig order",
            {
                "baseline_contigs": len(base_contigs),
                "test_contigs": len(test_contigs),
                "only_in_baseline": [n for n in base_names if n not in set(test_names)][:10],
                "only_in_test": [n for n in test_names if n not in set(base_names)][:10],
                "order_only": sorted(base_names) == sorted(test_names),
            },
        )

    base_samples = header_samples(base_lines)
    test_samples = header_samples(test_lines)
    same_samples = base_samples == test_samples
    report.check(
        "header.samples",
        same_samples,
        "%d sample(s): %s" % (len(base_samples), ", ".join(base_samples) or "none"),
        "sample names/order differ",
        {"baseline": base_samples, "test": test_samples},
    )

    base_defs = header_definitions(base_lines)
    test_defs = header_definitions(test_lines)
    report.check(
        "header.definitions",
        base_defs == test_defs,
        "%d INFO/FORMAT/FILTER definitions identical" % len(base_defs),
        "INFO/FORMAT/FILTER definitions differ",
        {
            "only_in_baseline": sorted(base_defs - test_defs)[:15],
            "only_in_test": sorted(test_defs - base_defs)[:15],
        },
        fail_status=rp.WARN,
    )
    return same_contigs, same_samples


# --------------------------------------------------------------------------- records


def view_cmd(bcftools_bin, path, region, samples, header=False):
    cmd = [bcftools_bin, "view", "-H" if not header else "-h"]
    if region:
        cmd += ["-r", region]
    if samples:
        cmd += ["-s", samples]
    cmd.append(path)
    return cmd


def make_canonicaliser(ignore_info, ignore_qual, sort_info=False):
    """Return a line transformer used before hashing.

    Default (ignore_info=[], sort_info=False, ignore_qual=False) is a pure
    identity on the VCF body columns -- INFO and FORMAT/sample fields are
    compared as emitted, so two files must be byte-identical at the record
    level (modulo trailing newline) to get matching digests.
    """
    ignored = set(ignore_info or [])

    def canonicalise(line):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            return line.rstrip("\n")
        if ignore_qual:
            fields[5] = "."
        if ignored or sort_info:
            info = [
                item for item in fields[7].split(";")
                if item and item.split("=", 1)[0] not in ignored
            ]
            if sort_info:
                info.sort()
            fields[7] = ";".join(info)
        # FORMAT + sample columns (fields[8:]) kept verbatim.
        return "\t".join(fields)

    return canonicalise


def digest_vcf(bcftools_bin, path, region, samples, canonicalise):
    digest = StreamDigest()
    for line in tools.stream_lines(view_cmd(bcftools_bin, path, region, samples)):
        digest.update(canonicalise(line))
    return digest


def file_sha256(path, chunk=8 * 1024 * 1024):
    """Hash the on-disk bytes (gzip container). Used as a fast identity shortcut."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def _record_locus(line):
    """Pull CHROM:POS REF>ALT from a VCF body line for diff captions."""
    fields = line.split("\t")
    if len(fields) < 5:
        return "?"
    return "%s:%s %s>%s" % (fields[0], fields[1], fields[3], fields[4])


def call_key(line):
    """Identity that ignores annotation jitter: site, alleles, FILTER, GT.

    QUAL, INFO and the non-GT FORMAT fields (AD, DP, AF, TLOD, ...) are not
    part of the key. Used by the somatic tolerance path only.
    """
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 8:
        return None
    chrom, pos, ref, alt, filt = fields[0], fields[1], fields[3], fields[4], fields[6]
    gts = ()
    if len(fields) >= 10:
        keys = fields[8].split(":")
        try:
            gi = keys.index("GT")
        except ValueError:
            gi = None
        parsed = []
        for sample in fields[9:]:
            if gi is None:
                parsed.append(".")
                continue
            vals = sample.split(":")
            raw = vals[gi] if gi < len(vals) else "."
            parsed.append(normalise_gt(raw))
        gts = tuple(parsed)
    return (chrom, pos, ref, alt, filt, gts)


def lockstep_digest_and_diffs(bcftools_bin, base, test, region, samples, canonicalise, limit,
                              classify_annotations=False):
    """One pass over both bodies: digest both streams and collect first diffs.

    Wall time is roughly one full read of the slower side (both bcftools view
    processes run concurrently via the two iterators' underlying pipes).

    When classify_annotations is set, also count mismatched lines whose
    call_key still matches (annotation-only) versus lines that move a site,
    allele, FILTER or GT (structural).
    """
    base_iter = tools.stream_lines(view_cmd(bcftools_bin, base, region, samples))
    test_iter = tools.stream_lines(view_cmd(bcftools_bin, test, region, samples))
    base_digest = StreamDigest()
    test_digest = StreamDigest()
    diffs = []
    annotation_only = 0
    structural = 0
    index = 0
    try:
        while True:
            base_line = next(base_iter, None)
            test_line = next(test_iter, None)
            if base_line is None and test_line is None:
                break
            index += 1
            left = canonicalise(base_line) if base_line is not None else None
            right = canonicalise(test_line) if test_line is not None else None
            if left is not None:
                base_digest.update(left)
            if right is not None:
                test_digest.update(right)
            if left == right:
                continue
            if classify_annotations:
                if (base_line is not None and test_line is not None
                        and call_key(base_line) == call_key(test_line)):
                    annotation_only += 1
                else:
                    structural += 1
            if len(diffs) < limit:
                locus = _record_locus(left if left is not None else right)
                diffs.append({
                    "record": index,
                    "locus": locus,
                    "baseline": (left if left is not None else "<end of file>")[:500],
                    "test": (right if right is not None else "<end of file>")[:500],
                })
    finally:
        base_iter.close()
        test_iter.close()
    return base_digest, test_digest, diffs, annotation_only, structural


def compare_records(report, args, config, bcftools_bin):
    canonicalise = make_canonicaliser(
        config.get("ignore_info") or [],
        config.get("ignore_qual", False),
        config.get("sort_info", False),
    )
    detail = {}
    diffs = []

    # Fast path: identical on-disk bytes => identical decoded records. Skip the
    # multi-GB bcftools view stream (common when golden was copied from work).
    if (not args.region and not args.samples
            and os.path.getsize(args.baseline) == os.path.getsize(args.test)):
        try:
            b_hash = file_sha256(args.baseline)
            t_hash = file_sha256(args.test)
        except OSError as exc:
            b_hash = t_hash = None
            detail["file_hash_error"] = str(exc)
        else:
            detail["file_sha256"] = {"baseline": b_hash, "test": t_hash}
            if b_hash == t_hash:
                report.add(
                    "records.digest",
                    rp.PASS,
                    "on-disk files identical (sha256 %s…); skipped body stream "
                    "(full INFO+FORMAT implied)" % b_hash[:12],
                    detail,
                )
                return 0, 0

    # Single lockstep pass: parallel I/O via two bcftools pipes + digests + diffs.
    ann_tol = float(config.get("annotation_rel_tol") or 0)
    try:
        base, test, diffs, annotation_only, structural = lockstep_digest_and_diffs(
            bcftools_bin,
            args.baseline,
            args.test,
            args.region,
            args.samples,
            canonicalise,
            int(config.get("max_reported_diffs") or 20),
            classify_annotations=ann_tol > 0,
        )
    except tools.CommandFailed as exc:
        report.add("records.digest", rp.ERROR, str(exc))
        return 0, 0

    detail["baseline"] = base.as_dict()
    detail["test"] = test.as_dict()

    if base.count == 0 and test.count == 0:
        # Two empty streams hash identically. Passing here would turn a typo in
        # --region, or an empty output file, into a green run.
        report.add(
            "records.digest",
            rp.ERROR,
            "no records on either side -- nothing was compared"
            + (" (does --region %s match a contig in the header?)" % args.region
               if args.region else ""),
            detail,
        )
        return base.count, test.count

    if base.ordered == test.ordered:
        report.add("records.digest", rp.PASS,
                   "%d records identical (full INFO+FORMAT)" % base.count, detail)
    elif (ann_tol > 0 and structural == 0 and base.count == test.count
          and base.count and annotation_only / float(base.count) <= ann_tol):
        frac = annotation_only / float(base.count)
        detail["annotation_only"] = annotation_only
        detail["structural"] = 0
        detail["annotation_fraction"] = frac
        detail["annotation_rel_tol"] = ann_tol
        detail["record_count_delta"] = 0
        detail["first_differences"] = diffs
        report.add(
            "records.digest", rp.PASS,
            "%d records; %d differ in INFO/FORMAT only (sites, alleles, FILTER "
            "and GT identical; %.4f%% <= %.2f%% somatic tolerance)"
            % (base.count, annotation_only, 100.0 * frac, 100.0 * ann_tol),
            detail,
        )
    elif args.mode == "truth":
        report.add("records.digest", rp.INFO,
                   "baseline %d records, test %d (differences expected in truth mode)"
                   % (base.count, test.count), detail)
    elif base.unordered == test.unordered and base.count == test.count:
        report.add("records.digest", rp.WARN,
                   "same %d records in a different order" % base.count, detail)
    else:
        abs_tol = int(config.get("count_abs_tol") or 0)
        rel = rp.rel_diff(base.count, test.count)
        within = (abs(base.count - test.count) <= abs_tol) or (rel <= config["count_rel_tol"])
        detail["record_count_delta"] = test.count - base.count
        detail["first_differences"] = diffs
        if within and (config["count_rel_tol"] > 0 or abs_tol > 0):
            report.add("records.digest", rp.WARN,
                       "records differ within tolerance (baseline %d, test %d, "
                       "rel %.6f)" % (base.count, test.count, rel), detail)
        else:
            report.add("records.digest", rp.FAIL,
                       "records differ (baseline %d, test %d; full INFO+FORMAT compare)"
                       % (base.count, test.count), detail)
    return base.count, test.count


# --------------------------------------------------------------------------- concordance


def iter_records(bcftools_bin, path, region, samples, rank, strict_phase, label):
    """Yield Records in coordinate order, refusing to continue on unsorted input."""
    cmd = [bcftools_bin, "query", "-f", QUERY_FORMAT]
    if region:
        cmd += ["-r", region]
    if samples:
        cmd += ["-s", samples]
    cmd.append(path)

    previous = None
    for line in tools.stream_lines(cmd):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            continue
        contig, pos, ref, alt, qual, filt = fields[:6]
        key = (rank.get(contig, len(rank)), contig, int(pos), ref, alt)
        if previous is not None and key < previous:
            raise ValueError(
                "%s is not coordinate sorted (or its contig order disagrees with the "
                "header) at %s:%s" % (label, contig, pos)
            )
        previous = key
        yield Record(
            key, contig, int(pos), ref, alt, qual, filt,
            tuple(normalise_gt(g, strict_phase) for g in fields[6:]),
            classify(ref, alt),
        )


def merge_join(left, right):
    a = next(left, None)
    b = next(right, None)
    while a is not None or b is not None:
        if b is None or (a is not None and a.key < b.key):
            yield a, None
            a = next(left, None)
        elif a is None or b.key < a.key:
            yield None, b
            b = next(right, None)
        else:
            yield a, b
            a = next(left, None)
            b = next(right, None)


def compare_concordance(report, args, config, bcftools_bin, base_lines, compare_gt):
    rank = {name: index for index, (name, _) in enumerate(header_contigs(base_lines))}
    samples = header_samples(base_lines)

    shared = collections.Counter()
    only_base = collections.Counter()
    only_test = collections.Counter()
    gt_match = collections.Counter()
    gt_total = collections.Counter()
    filter_base = collections.Counter()
    filter_test = collections.Counter()
    examples = []

    left = iter_records(bcftools_bin, args.baseline, args.region, args.samples, rank,
                        args.strict_phase, "baseline")
    right = iter_records(bcftools_bin, args.test, args.region, args.samples, rank,
                         args.strict_phase, "test")

    try:
        for a, b in merge_join(left, right):
            if a is not None:
                filter_base[a.filt] += 1
            if b is not None:
                filter_test[b.filt] += 1
            if a is not None and b is not None:
                shared[a.vtype] += 1
                if compare_gt and a.gts and len(a.gts) == len(b.gts):
                    for index, (left_gt, right_gt) in enumerate(zip(a.gts, b.gts)):
                        name = samples[index] if index < len(samples) else "sample%d" % index
                        gt_total[name] += 1
                        if left_gt == right_gt:
                            gt_match[name] += 1
                        elif len(examples) < config["max_reported_diffs"]:
                            examples.append("%s:%d %s>%s %s: %s vs %s"
                                            % (a.contig, a.pos, a.ref, a.alt, name,
                                               left_gt, right_gt))
            elif a is not None:
                only_base[a.vtype] += 1
                if len(examples) < config["max_reported_diffs"]:
                    examples.append("baseline only %s:%d %s>%s" % (a.contig, a.pos, a.ref, a.alt))
            else:
                only_test[b.vtype] += 1
                if len(examples) < config["max_reported_diffs"]:
                    examples.append("test only %s:%d %s>%s" % (b.contig, b.pos, b.ref, b.alt))
    except ValueError as exc:
        report.add("concordance.sites", rp.ERROR, str(exc))
        return

    total_shared = sum(shared.values())
    total_base_only = sum(only_base.values())
    total_test_only = sum(only_test.values())
    union = total_shared + total_base_only + total_test_only
    if union == 0:
        report.add(
            "concordance.sites",
            rp.ERROR,
            "no variants on either side -- nothing was compared"
            + (" (does --region %s match a contig in the header?)" % args.region
               if args.region else ""),
        )
        return
    jaccard = total_shared / union
    precision = total_shared / (total_shared + total_test_only) if (total_shared + total_test_only) else 0.0
    recall = total_shared / (total_shared + total_base_only) if (total_shared + total_base_only) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

    detail = {
        "shared": total_shared,
        "baseline_only": total_base_only,
        "test_only": total_test_only,
        "jaccard": round(jaccard, 6),
        "precision": round(precision, 6),
        "recall": round(recall, 6),
        "f1": round(f1, 6),
        "by_type_shared": dict(shared),
        "by_type_baseline_only": dict(only_base),
        "by_type_test_only": dict(only_test),
        "examples": examples,
    }

    failures = []
    for name, value in (("site_jaccard_min", jaccard), ("precision_min", precision),
                        ("recall_min", recall), ("f1_min", f1)):
        minimum = config.get(name)
        if minimum is not None and value < minimum:
            failures.append("%s=%.6f < %s" % (name[:-4], value, minimum))

    report.check(
        "concordance.sites",
        not failures,
        "shared %d, baseline-only %d, test-only %d (jaccard %.6f, F1 %.6f)"
        % (total_shared, total_base_only, total_test_only, jaccard, f1),
        "site concordance below threshold: %s" % "; ".join(failures),
        detail,
    )

    if compare_gt and sum(gt_total.values()):
        per_sample = {
            name: round(gt_match[name] / total, 6) for name, total in gt_total.items() if total
        }
        overall = sum(gt_match.values()) / sum(gt_total.values())
        minimum = config.get("gt_concordance_min")
        report.check(
            "concordance.genotypes",
            minimum is None or overall >= minimum,
            "genotype concordance %.6f over %d shared genotypes" % (overall, sum(gt_total.values())),
            "genotype concordance %.6f below %s" % (overall, minimum),
            {"overall": round(overall, 6), "per_sample": per_sample,
             "mismatched": sum(gt_total.values()) - sum(gt_match.values())},
        )
    else:
        report.add("concordance.genotypes", rp.SKIP, "no comparable genotype columns")

    compare_counter(report, "counts.filter", filter_base, filter_test, config["count_rel_tol"],
                    "FILTER", config.get("count_abs_tol", 0))


def compare_counter(report, name, base, test, tolerance, label, abs_tol=0):
    offenders = {}
    for key in sorted(set(base) | set(test)):
        left, right = base.get(key, 0), test.get(key, 0)
        # a small absolute drift is noise for low-count categories (a single
        # record flipping FILTER in a category of ~450 reads is a 0.2% rel
        # diff but not a regression); gate on both relative and absolute.
        if rp.rel_diff(left, right) > tolerance and abs(left - right) > abs_tol:
            offenders[key] = {"baseline": left, "test": right,
                              "rel_diff": round(rp.rel_diff(left, right), 6)}
    report.check(
        name,
        not offenders,
        "%s breakdown within tolerance (%d categories)" % (label, len(set(base) | set(test))),
        "%s breakdown differs in %d category(ies)" % (label, len(offenders)),
        offenders or {"baseline": dict(base), "test": dict(test)},
    )


# --------------------------------------------------------------------------- stats


def parse_stats(bcftools_bin, path, region, samples):
    cmd = [bcftools_bin, "stats"]
    if region:
        cmd += ["-r", region]
    if samples:
        cmd += ["-s", samples]
    cmd.append(path)
    metrics = {}
    for line in tools.run(cmd).stdout.splitlines():
        if line.startswith("SN\t"):
            parts = line.split("\t")
            if len(parts) >= 4 and parts[2] in STATS_KEYS:
                try:
                    metrics[parts[2].rstrip(":")] = int(parts[3])
                except ValueError:
                    pass
        elif line.startswith("TSTV\t"):
            parts = line.split("\t")
            if len(parts) >= 5:
                try:
                    metrics["ts/tv"] = float(parts[4])
                except ValueError:
                    pass
    return metrics


def compare_stats(report, args, config, bcftools_bin):
    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            fb = pool.submit(parse_stats, bcftools_bin, args.baseline,
                             args.region, args.samples)
            ft = pool.submit(parse_stats, bcftools_bin, args.test,
                             args.region, args.samples)
            base, test = fb.result(), ft.result()
    except tools.CommandFailed as exc:
        report.add("stats.summary", rp.ERROR, str(exc))
        return
    _compare_stats_maps(report, config, base, test)


def _compare_stats_maps(report, config, base, test):
    offenders = {}
    for key in sorted(set(base) | set(test)):
        left, right = base.get(key, 0), test.get(key, 0)
        tolerance = config["tstv_rel_tol"] if key == "ts/tv" else config["count_rel_tol"]
        delta = rp.rel_diff(left, right)
        if delta > tolerance:
            offenders[key] = {"baseline": left, "test": right, "rel_diff": round(delta, 6)}
    report.check(
        "stats.summary",
        not offenders,
        "bcftools stats metrics within tolerance (%s)"
        % ", ".join("%s=%s" % (k, base[k]) for k in sorted(base)),
        "%d bcftools stats metric(s) out of tolerance" % len(offenders),
        offenders or {"baseline": base, "test": test},
    )


# --------------------------------------------------------------------------- normalisation


def normalise_inputs(bcftools_bin, paths, reference, workdir):
    """Split multi-allelics and left-align so the two callers use one representation."""
    normalised = []
    for index, path in enumerate(paths):
        out = os.path.join(workdir, "norm%d.vcf.gz" % index)
        tools.run([bcftools_bin, "norm", "-m", "-any", "-f", reference,
                   "-Oz", "-o", out, path])
        tools.run([bcftools_bin, "index", "-f", "-t", out])
        normalised.append(out)
    return normalised


# --------------------------------------------------------------------------- main


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-b", "--baseline", required=True, help="golden VCF, or the truth set")
    parser.add_argument("-t", "--test", required=True, help="VCF produced by the build under test")
    parser.add_argument("--mode", choices=("regression", "truth"), default="regression",
                        help="regression: must match the golden; truth: judge on precision/recall")
    parser.add_argument("--bcftools-bin", help="path to bcftools (default: $BCFTOOLS, env/tools.conf, PATH)")
    parser.add_argument("--region", help="restrict to a region, e.g. chr20 (needs a .tbi/.csi index)")
    parser.add_argument("--samples", help="comma separated sample subset passed to bcftools -s")
    parser.add_argument("--norm-ref", metavar="FASTA",
                        help="split multi-allelics and left-align both inputs first")
    parser.add_argument("--ignore-info", action="append", default=[], metavar="KEY",
                        help="drop this INFO key before hashing (repeatable)")
    parser.add_argument("--ignore-qual", action="store_true", help="blank QUAL before hashing")
    parser.add_argument(
        "--annotation-rel-tol", type=float, default=None, metavar="FRAC",
        help="somatic tolerance: PASS when at most this fraction of records "
             "differ, and every difference is confined to QUAL/INFO/FORMAT "
             "annotations (chrom, pos, ref, alt, FILTER and GT must still "
             "match). Unset or 0 keeps the strict WGS comparison",
    )
    parser.add_argument("--strict-phase", action="store_true",
                        help="treat 0|1 and 1|0 as different genotypes")
    parser.add_argument("--skip-stats", action="store_true", help="skip the bcftools stats check")
    parser.add_argument("--skip-concordance", action="store_true",
                        help="skip site/genotype concordance (also auto-skipped in "
                             "regression mode when records.digest already PASS/FAIL)")
    return rp.add_common_args(parser)


def main(argv=None):
    args = build_parser().parse_args(argv)
    defaults = dict(DEFAULTS)
    if args.mode == "truth":
        defaults.update(TRUTH_DEFAULTS)
    config = rp.load_config(defaults, args.config)
    if args.ignore_info:
        config["ignore_info"] = sorted(set(config["ignore_info"]) | set(args.ignore_info))
    if args.ignore_qual:
        config["ignore_qual"] = True
    if args.annotation_rel_tol is not None:
        config["annotation_rel_tol"] = args.annotation_rel_tol

    report = rp.Report(
        "vcf_compare",
        baseline=args.baseline,
        test=args.test,
        mode=args.mode,
        region=args.region or "whole file",
    )

    try:
        bcftools_bin = tools.bcftools(args.bcftools_bin)
    except tools.ToolNotFound as exc:
        report.add("tool.bcftools", rp.ERROR, str(exc))
        return report.finish(args.json, args.verbose, args.warn_is_failure)
    report.add("tool.bcftools", rp.INFO, bcftools_bin)

    missing = [p for p in (args.baseline, args.test) if not os.path.isfile(p)]
    if missing:
        report.add("file.exists", rp.ERROR, "missing input: %s" % ", ".join(missing))
        return report.finish(args.json, args.verbose, args.warn_is_failure)
    report.add("file.exists", rp.PASS, "both inputs present")

    workdir = None
    try:
        if args.norm_ref:
            workdir = tempfile.mkdtemp(prefix="vcf_compare_")
            args.baseline, args.test = normalise_inputs(
                bcftools_bin, [args.baseline, args.test], args.norm_ref, workdir
            )
            report.add("records.normalised", rp.INFO,
                       "both inputs split with bcftools norm -m -any")

        base_lines = read_header(bcftools_bin, args.baseline)
        test_lines = read_header(bcftools_bin, args.test)
        same_contigs, same_samples = compare_header(report, base_lines, test_lines)

        compare_records(report, args, config, bcftools_bin)

        # Concordance re-streams the whole file. Only skip it when the body
        # digest already FAILed (first_differences are enough); on PASS keep
        # concordance so the HTML report still has Jaccard / GT rows like 30x.
        digest_status = next(
            (c["status"] for c in report.checks if c["name"] == "records.digest"),
            None,
        )
        skip_conc = (
            args.skip_concordance
            or (args.mode == "regression" and digest_status == rp.FAIL)
        )
        if skip_conc:
            reason = "records.digest already decisive (%s)" % (digest_status or "?")
            if args.skip_concordance:
                reason = "--skip-concordance"
            report.add("concordance.sites", rp.SKIP,
                       "skipped: %s" % reason)
        elif same_contigs:
            compare_concordance(report, args, config, bcftools_bin, base_lines, same_samples)
        else:
            report.add("concordance.sites", rp.SKIP,
                       "skipped: the merge join needs a shared contig order")

        if not args.skip_stats:
            compare_stats(report, args, config, bcftools_bin)
    except tools.CommandFailed as exc:
        report.add("run", rp.ERROR, str(exc))
    except KeyboardInterrupt:
        report.add("run", rp.ERROR, "interrupted")
    finally:
        if workdir:
            shutil.rmtree(workdir, ignore_errors=True)

    return report.finish(args.json, args.verbose, args.warn_is_failure)


if __name__ == "__main__":
    sys.exit(main())
