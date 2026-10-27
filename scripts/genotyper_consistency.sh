#!/bin/bash
# genotyper_consistency.sh -- run genotyper N times on the SAME gVCF and check
# whether the output VCFs are identical.
#
# This isolates genotyper run-to-run determinism from the rest of the pipeline:
# the input gVCF is fixed, so any difference between runs is produced inside
# genotyper itself (thread scheduling, floating-point reduction order, ...).
# The QD diffs seen in the 20260929_human_retest run came from two FULL
# pipeline runs (fresh align + bqsr + variantCaller each time); this script
# answers whether genotyper alone can produce such drift:
#   - outputs identical  -> genotyper is deterministic on this input; the QD
#                           drift originates upstream (variantCaller / bqsr /
#                           aligner), most likely the known BQSR race.
#   - outputs differ     -> genotyper itself is nondeterministic; the differing
#                           field histogram shows whether it is QD-only again.
#
# Usage:
#   scripts/genotyper_consistency.sh -i sample.g.vcf.gz [options] [-- extra genotyper args]
#
#   -i|--input PATH   input gVCF (.g.vcf.gz). Default: first existing of
#                       $WORK_ROOT/wgs_human_30x/sample.variantCaller.g.vcf.gz
#                       $WORK_ROOT/wgs_human_30x_ert/sample.variantCaller.g.vcf.gz
#                       $WORK_ROOT/wgs_human_small/*.variantCaller.g.vcf.gz
#   -n|--runs N       how many times to run genotyper (default 2)
#   -o|--outdir DIR   where to put runK.vcf.gz + reports
#                       (default $WORK_ROOT/genotyper_consistency/<timestamp>)
#   -s FLOAT          -s stand_call_conf passed to genotyper (default 30.0)
#   --threads N       genotyper threads (default $THREADS from env/cases.conf)
#
# Env (picked up from env/cases.conf when run inside the suite):
#   DCS_HOME   release bundle under test; point it at another bundle to test
#              a different genotyper build.
#   PROXY_HOST license proxy (must already be running).
#
# Exit: 0 all runs identical, 1 outputs differ, 2 setup/run error.

set -uo pipefail

SUITE_ROOT=$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/.." && pwd)
[ -f "$SUITE_ROOT/env/cases.conf" ] && . "$SUITE_ROOT/env/cases.conf"
export DCS_HOME PROXY_HOST
DCS="$DCS_HOME/bin/dcs"

INPUT=""
RUNS=2
OUT=""
STAND_CONF=30.0
THREADS=${THREADS:-16}
EXTRA=()
while [ $# -gt 0 ]; do
    case "$1" in
        -i|--input)  INPUT=$2; shift 2 ;;
        -n|--runs)   RUNS=$2; shift 2 ;;
        -o|--outdir) OUT=$2; shift 2 ;;
        -s)          STAND_CONF=$2; shift 2 ;;
        --threads)   THREADS=$2; shift 2 ;;
        --)          shift; EXTRA+=("$@"); break ;;
        -h|--help)   sed -n '2,33p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown argument: $1 (see --help)" >&2; exit 2 ;;
    esac
done

if [ "$RUNS" -lt 2 ]; then
    echo "error: -n must be >= 2" >&2
    exit 2
fi

if [ -z "$INPUT" ]; then
    for cand in "$WORK_ROOT"/wgs_human_30x/sample.variantCaller.g.vcf.gz \
                "$WORK_ROOT"/wgs_human_30x_ert/sample.variantCaller.g.vcf.gz \
                "$WORK_ROOT"/wgs_human_small/*.variantCaller.g.vcf.gz; do
        if [ -f "$cand" ]; then INPUT=$cand; break; fi
    done
fi
if [ -z "$INPUT" ] || [ ! -f "$INPUT" ]; then
    echo "error: input gVCF not found; pass -i <file.g.vcf.gz>" >&2
    exit 2
fi
if [ ! -x "$DCS" ]; then
    echo "error: dcs not executable: $DCS (set DCS_HOME)" >&2
    exit 2
fi
if ! "$DCS" lickit --ping "$PROXY_HOST" >/dev/null 2>&1; then
    echo "error: no license proxy at $PROXY_HOST (start it with: dcs licproxy --license <file>)" >&2
    exit 2
fi

OUT=${OUT:-"${WORK_ROOT:-$SUITE_ROOT/reports}/genotyper_consistency/$(date +%Y%m%d_%H%M%S)"}
mkdir -p "$OUT"

echo "input    : $INPUT"
echo "dcs      : $DCS"
echo "outdir   : $OUT"
echo "runs     : $RUNS (genotyper --threads $THREADS -s $STAND_CONF ${EXTRA[*]+"${EXTRA[*]}"})"
echo

# -- run genotyper N times on the same input ----------------------------------
for k in $(seq 1 "$RUNS"); do
    out_vcf="$OUT/run$k.vcf.gz"
    echo ">>> run $k -> $out_vcf"
    start=$(date +%s)
    "$DCS" genotyper \
        --threads "$THREADS" \
        --input "$INPUT" \
        --output "$out_vcf" \
        -s "$STAND_CONF" ${EXTRA[@]+"${EXTRA[@]}"} \
        > "$OUT/run$k.log" 2>&1
    rc=$?
    if [ "$rc" -ne 0 ]; then
        echo "error: genotyper run $k exited $rc (see $OUT/run$k.log)" >&2
        exit 2
    fi
    echo "    ok in $(( $(date +%s) - start ))s"
done
echo

# -- field-level diff across runs (no external deps) ---------------------------
python3 - "$OUT" "$RUNS" <<'PY'
import gzip, sys
from collections import Counter

outdir, runs = sys.argv[1], int(sys.argv[2])

def parse(path):
    recs, order = {}, []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            key = (cols[0], cols[1], cols[3], cols[4])
            order.append(key)
            recs[key] = cols
    return recs, order

data = {k: parse(f"{outdir}/run{k}.vcf.gz") for k in range(1, runs + 1)}
base_recs, base_order = data[1]

print("records per run: " + ", ".join(f"run{k}={len(data[k][0])}" for k in data))
same_keys = all(set(data[k][0]) == set(base_recs) for k in data)
same_order = all(data[k][1] == base_order for k in data)
print(f"locus sets identical : {same_keys}")
print(f"record order identical: {same_order}")

field_hist = Counter()
qd_only, other = [], []
for key in base_order:
    lines = {k: data[k][0].get(key) for k in data}
    lines = {k: v for k, v in lines.items() if v is not None}
    if len(lines) < runs:
        continue  # locus missing in some run; already flagged by same_keys
    if len({"\t".join(c) for c in lines.values()}) == 1:
        continue  # identical record in every run
    infos = {k: dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
             for k, c in lines.items()}
    all_f = set().union(*[set(d) for d in infos.values()])
    diff_f = {f for f in all_f if len({d.get(f) for d in infos.values()}) > 1}
    ncols = min(len(c) for c in lines.values())
    fmt_diff = any(len({lines[k][i] for k in lines}) > 1 for i in range(8, ncols))
    for f in sorted(diff_f):
        field_hist[f] += 1
    if fmt_diff:
        field_hist["<FORMAT/sample>"] += 1
    entry = {
        "locus": f"{key[0]}:{key[1]} {key[2]}>{key[3]}",
        "fields": sorted(diff_f) + (["<FORMAT/sample>"] if fmt_diff else []),
        "qual": {k: lines[k][5] for k in sorted(lines)},
        "qd":   {k: infos[k].get("QD") for k in sorted(lines)},
        "dp":   {k: infos[k].get("DP") for k in sorted(lines)},
        "ad":   {k: (lines[k][9].split(":")[1] if len(lines[k]) > 9 and ":" in lines[k][9] else None)
                 for k in sorted(lines)},
    }
    (qd_only if diff_f == {"QD"} and not fmt_diff else other).append(entry)

total = len(qd_only) + len(other)
print(f"differing loci: {total}  (QD-only: {len(qd_only)}, other: {len(other)})")
if field_hist:
    print("differing-field histogram: " + ", ".join(f"{f} x{n}" for f, n in field_hist.most_common()))
for e in (qd_only + other)[:20]:
    qd = " vs ".join(f"run{k}={v}" for k, v in e["qd"].items())
    print(f"  {e['locus']}  fields={e['fields']}  QD: {qd}")
    print(f"      QUAL={e['qual']}  DP={e['dp']}  AD={e['ad']}")
sys.exit(1 if total else 0)
PY
diff_rc=$?
echo

# -- rigorous full-record comparison via the suite checker ---------------------
compare_rc=0
for k in $(seq 2 "$RUNS"); do
    json="$OUT/compare_run1_vs_run$k.json"
    python3 "$SUITE_ROOT/scripts/vcf_compare.py" \
        -b "$OUT/run1.vcf.gz" -t "$OUT/run$k.vcf.gz" --json "$json" >>"$OUT/compare.log" 2>&1
    rc=$?
    case "$rc" in
        0) echo "vcf_compare run1 vs run$k: IDENTICAL" ;;
        1) echo "vcf_compare run1 vs run$k: DIFFERS (detail: $json)"; compare_rc=1 ;;
        *) echo "vcf_compare run1 vs run$k: checker error (see $OUT/compare.log)" ;;
    esac
done
echo

if [ "$diff_rc" = 0 ] && [ "$compare_rc" = 0 ]; then
    echo "=> CONSISTENT: $RUNS genotyper runs (--threads $THREADS) on the same gVCF produced identical VCFs."
    echo "   (deterministic at this thread count; try --threads 32 to probe the"
    echo "    fixTooHighQD shared-RNG jitter -- see env/cases.conf GENOTYPER_THREADS)"
    exit 0
else
    echo "=> INCONSISTENT: genotyper produced different VCFs from the SAME gVCF (--threads $THREADS)."
    echo "   products kept under $OUT for inspection."
    exit 1
fi
