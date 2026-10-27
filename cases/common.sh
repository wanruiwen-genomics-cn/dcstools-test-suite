#!/bin/bash
# Shared runtime for the case scripts in this directory.
#
# A case script describes one scenario: which data goes in, which dcs commands
# run, and which products are compared against the golden baseline. Everything
# else -- preflight, resource accounting, golden handling, reporting, exit codes
# -- lives here so the cases stay readable.
#
# Usage inside a case:
#
#   source "$(dirname "$0")/common.sh"
#   case_begin wgs_human_small
#   need_file FA "$DATA_ROOT/ref/chr20.fa" "reference"
#   preflight_done
#   run_step align "$DCS" aligner --threads "$THREADS" ...
#   expect_file "$WORK/sample.align.bam"
#   check_bam align_bam "$WORK/sample.align.bam"
#   case_end
#
# Exit codes match the python checkers: 0 pass, 1 a check failed, 2 the case
# could not run (missing data, missing license, tool crash).

set -uo pipefail

SUITE_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SCRIPTS_DIR="$SUITE_ROOT/scripts"

# --update-golden installs this run's products as the new baseline instead of
# comparing against it. Deliberately explicit: silently refreshing a golden on
# failure is how a regression suite stops catching regressions.
UPDATE_GOLDEN=0
KEEP_WORK=0
COMPARE_ONLY=0
for arg in "$@"; do
    case "$arg" in
        --update-golden) UPDATE_GOLDEN=1 ;;
        --keep-work)     KEEP_WORK=1 ;;
        --compare-only)  COMPARE_ONLY=1; KEEP_WORK=1 ;;
        -h|--help)
            sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'
            exit 0 ;;
    esac
done

[ -f "$SUITE_ROOT/env/cases.conf" ] && . "$SUITE_ROOT/env/cases.conf"

# The tools read PROXY_HOST and DCS_HOME from the environment, not from their
# arguments. Exporting them is what makes the license check in the preflight
# test the same mechanism the pipeline will use.
export DCS_HOME PROXY_HOST

DCS="$DCS_HOME/bin/dcs"

_missing_data=()
_checks=()
_steps=()
_fail_count=0
_error_count=0
# A case that bailed out before doing any work is not a case that passed. Kept
# separate from the failure counters because it must not turn the run red
# either -- "no JDK 8 on this host" is a fact about the host, not the build.
_skipped=0

_now() { date +%s.%N; }
_stamp() { date +%Y%m%d_%H%M%S; }

_log() { echo "$*" | tee -a "$LOG" ; }
_rule() { _log "------------------------------------------------------------"; }

_record_check() {
    # name | status | message [detail_json]
    # detail_json is an optional 4th arg: a JSON string embedded in the pipe-
    # delimited record.  Empty or absent means no detail.
    local name=$1 status=$2 message=$3 detail=${4:-}
    _checks+=("$name|$status|$message|$detail")
    case "$status" in
        FAIL) _fail_count=$((_fail_count + 1)) ;;
        ERROR) _error_count=$((_error_count + 1)) ;;
    esac
    printf '  [%-5s] %-22s %s\n' "$status" "$name" "$message" | tee -a "$LOG"
}

case_begin() {
    CASE_NAME=$1
    # run.sh exports SUITE_RUN_ID so every case of one suite run lands under the
    # same stamp and the runner can find each report without parsing its output.
    RUN_ID=${SUITE_RUN_ID:-$(_stamp)}
    WORK="$WORK_ROOT/$CASE_NAME"
    GOLDEN="$GOLDEN_ROOT/$CASE_NAME"
    REPORT_DIR="$REPORT_ROOT/$CASE_NAME/$RUN_ID"
    mkdir -p "$WORK" "$REPORT_DIR"
    # Start from a clean work directory. Products from a previous run (e.g. a
    # golden-baseline run with --update-golden) would otherwise collide with
    # the current run: SeqArc refuses to overwrite an existing .arc file, and
    # stale BAMs mask alignment regressions. The reports/ and golden/ trees
    # preserve previous results; work/ is ephemeral.
    # --compare-only skips the cleanup and the pipeline: it assumes the work
    # directory already holds the products from a previous run and jumps
    # straight to the golden comparison. This is the recovery path when a
    # checker crashed after a long pipeline run and you do not want to rerun
    # the 2-hour alignment just to redo a 30-second digest.
    if [ "$COMPARE_ONLY" = 0 ] && [ "${RESUME_VC:-0}" = 0 ]; then
        rm -rf "$WORK"/* "$WORK"/.* 2>/dev/null || true
    fi
    LOG="$REPORT_DIR/run.log"
    : > "$LOG"
    CASE_START=$(_now)

    echo
    echo "=== case: $CASE_NAME  ($RUN_ID) ==="
    _log "dcs      : $DCS"
    _log "data     : $DATA_ROOT"
    _log "work     : $WORK"
    _log "golden   : $GOLDEN"
    _log "threads  : $THREADS"
    [ "$UPDATE_GOLDEN" = 1 ] && _log "mode     : UPDATE GOLDEN"
    _rule

    if [ ! -x "$DCS" ]; then
        _record_check env.dcs ERROR "not executable: $DCS"
        case_end
    fi
}

# Collect every missing input before giving up, so the operator fixes the data
# layout in one pass instead of one file per run.
need_file() {
    local var=$1 path=$2 desc=${3:-}
    if [ -e "$path" ]; then
        printf -v "$var" '%s' "$path"
        export "$var"
    else
        _missing_data+=("$path  ($desc)")
    fi
}

need_dir() { need_file "$@"; }

preflight_done() {
    if [ ${#_missing_data[@]} -gt 0 ]; then
        _record_check data.present ERROR "${#_missing_data[@]} input(s) missing"
        for item in "${_missing_data[@]}"; do _log "        $item"; done
        _log "        set DATA_ROOT in env/cases.conf, or drop the data in place"
        case_end
    fi
    _record_check data.present PASS "all inputs present"

    # The license proxy gates every tool; failing here gives a one-line answer
    # instead of an opaque crash three steps into the pipeline.
    if "$DCS" lickit --ping "$PROXY_HOST" >>"$LOG" 2>&1; then
        _record_check license.proxy PASS "proxy reachable at $PROXY_HOST"
    else
        _record_check license.proxy ERROR "no license proxy at $PROXY_HOST (start it with: dcs licproxy --license <file>)"
        case_end
    fi
}

# note <name> <message> -- record a measurement that is reported but not gated.
note() { _record_check "$1" INFO "$2"; }

skip_case() {
    _skipped=1
    _record_check "$1" SKIP "$2"
    _log ""
    _log "case skipped: $2"
    _write_report
    exit 0
}

# run_step <label> <command...>
# Times the command, records peak RSS, and aborts the case on failure -- there
# is no point comparing products of a pipeline that did not finish.
run_step() {
    local label=$1; shift
    local timefile="$REPORT_DIR/$label.time"
    local outfile="$REPORT_DIR/$label.out"
    local start rc wall rss threads
    # Pick the thread count out of the command line so the steps table can show
    # it per step. Falls back to $THREADS when the tool has no --threads flag.
    threads=$THREADS
    local prev=""
    for a in "$@"; do
        case "$prev" in
            --threads|-threads|--thread) threads="$a" ;;
        esac
        case "$a" in
            --threads=*|-threads=*|--thread=*) threads="${a#*=}" ;;
        esac
        prev="$a"
    done
    _log ""
    _log ">>> step: $label"
    _log "    $*"
    # --compare-only: skip the pipeline step, assume the product already exists.
    if [ "$COMPARE_ONLY" = 1 ]; then
        _log "    [skipped: --compare-only]"
        _steps+=("$label|0|0|0|$threads")
        return 0
    fi
    start=$(_now)
    # Each step gets its own output file so a failure report shows that step's
    # output rather than the tail of everything that ran before it.
    if [ -x /usr/bin/time ]; then
        /usr/bin/time -v -o "$timefile" "$@" >"$outfile" 2>&1
        rc=$?
    else
        "$@" >"$outfile" 2>&1
        rc=$?
    fi
    cat "$outfile" >> "$LOG"
    wall=$(awk -v s="$start" -v e="$(_now)" 'BEGIN{printf "%.1f", e-s}')
    rss=$(awk -F': ' '/Maximum resident set size/{printf "%.0f", $2/1024}' "$timefile" 2>/dev/null)
    _steps+=("$label|$rc|$wall|${rss:-0}|$threads")
    if [ "$rc" -ne 0 ]; then
        _record_check "step.$label" ERROR "exit code $rc after ${wall}s (see $outfile)"
        grep -viE '^\s*(\[(memory|cpu|file)|Model:|Cores:|Arch:|Total:|Free:|Available:|Soft Limit:|Hard Limit:|-{10,})' \
            "$outfile" | tail -15 | sed 's/^/        /'
        case_end
    fi
    _record_check "step.$label" PASS "ok in ${wall}s, peak RSS ${rss:-?} MB"
}

# expect_file <path> [min_bytes]
expect_file() {
    local path=$1 min=${2:-1} size name
    name="product.$(basename "$path")"
    if [ "$COMPARE_ONLY" = 1 ]; then
        [ -f "$path" ] && _record_check "$name" PASS "exists (compare-only)" \
                     || _record_check "$name" FAIL "missing (compare-only needs prior run)"
        return
    fi
    if [ ! -f "$path" ]; then
        _record_check "$name" FAIL "not produced"
        return
    fi
    size=$(stat -L -c%s "$path")
    if [ "$size" -lt "$min" ]; then
        _record_check "$name" FAIL "only $size bytes (expected >= $min)"
        return
    fi
    # Keep text products (bamStat / vcfStat / recal.table) next to the BAM/VCF
    # golden so later reports can size-compare and diff their contents.
    if [ "$UPDATE_GOLDEN" = 1 ]; then
        mkdir -p "$GOLDEN"
        cp -f "$path" "$GOLDEN/$(basename "$path")"
    fi
    _record_check "$name" PASS "$size bytes"
}

# check_report <report_dir>
#
# The customer-facing report/ directory (HTML + SVG + xls). On --update-golden
# it is snapshotted next to the BAM/VCF goldens; on a normal run every file is
# compared by name and md5 so a silent chart/stat regression is caught.
check_report() {
    local product=$1
    local name=product.report
    local golden="$GOLDEN/report"

    if [ ! -d "$product" ]; then
        _record_check "$name" FAIL "report directory not produced"
        return
    fi

    local n_test
    n_test=$(find "$product" -type f | wc -l)

    if [ "$UPDATE_GOLDEN" = 1 ]; then
        mkdir -p "$GOLDEN"
        rm -rf "$golden"
        cp -a "$product" "$golden"
        _write_provenance
        _record_check "$name" INFO "golden updated: $golden ($n_test files)"
        return
    fi

    if [ ! -d "$golden" ]; then
        _record_check "$name" SKIP "no golden report yet (run once with --update-golden)" \
            "{\"test_files\":$n_test}"
        return
    fi

    local n_base
    n_base=$(find "$golden" -type f | wc -l)
    local detail json_rc
    detail=$(python3 - "$product" "$golden" <<'PY'
import hashlib, json, os, sys
test_dir, base_dir = sys.argv[1], sys.argv[2]

def walk(root):
    out = {}
    for dirpath, _, files in os.walk(root):
        for f in files:
            full = os.path.join(dirpath, f)
            rel = os.path.relpath(full, root)
            h = hashlib.md5()
            with open(full, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(chunk)
            out[rel.replace("\\", "/")] = {
                "md5": h.hexdigest(),
                "size": os.path.getsize(full),
            }
    return out

test, base = walk(test_dir), walk(base_dir)
only_test = sorted(set(test) - set(base))
only_base = sorted(set(base) - set(test))
shared = sorted(set(test) & set(base))
mismatch = []
for rel in shared:
    if test[rel]["md5"] != base[rel]["md5"]:
        mismatch.append({
            "file": rel,
            "test_size": test[rel]["size"],
            "baseline_size": base[rel]["size"],
            "test_md5": test[rel]["md5"],
            "baseline_md5": base[rel]["md5"],
        })
print(json.dumps({
    "test_files": len(test),
    "baseline_files": len(base),
    "only_in_test": only_test[:20],
    "only_in_baseline": only_base[:20],
    "mismatched": mismatch[:30],
    "mismatch_count": len(mismatch),
}, ensure_ascii=False))
sys.exit(1 if (only_test or only_base or mismatch) else 0)
PY
)
    json_rc=$?
    if [ "$json_rc" -eq 0 ]; then
        _record_check "$name" PASS \
            "$n_test files identical to golden ($n_base files)" "$detail"
    else
        local n_mis
        n_mis=$(python3 -c "import json,sys; print(json.load(sys.stdin).get('mismatch_count',0))" <<<"$detail" 2>/dev/null || echo "?")
        _record_check "$name" FAIL \
            "differs from golden (test $n_test / golden $n_base files, $n_mis content mismatch)" \
            "$detail"
    fi
}

_golden_path() { echo "$GOLDEN/$(basename "$1")"; }

# Shared by check_bam/check_vcf: handle the update and no-baseline-yet paths.
# Returns 0 when an actual comparison should run.
_golden_ready() {
    local name=$1 product=$2 golden
    golden=$(_golden_path "$product")
    if [ ! -f "$product" ]; then
        _record_check "$name" FAIL "product missing, nothing to compare"
        return 1
    fi
    if [ "$UPDATE_GOLDEN" = 1 ]; then
        mkdir -p "$GOLDEN"
        cp -f "$product" "$golden"
        for suffix in .bai .csi .tbi; do
            [ -f "$product$suffix" ] && cp -f "$product$suffix" "$golden$suffix"
        done
        _write_provenance
        _record_check "$name" INFO "golden updated: $golden"
        return 1
    fi
    if [ ! -f "$golden" ]; then
        _record_check "$name" SKIP "no golden yet (run once with --update-golden)"
        return 1
    fi
    return 0
}

_run_checker() {
    local name=$1 script=$2; shift 2
    local json="$REPORT_DIR/$name.json" rc
    python3 "$SCRIPTS_DIR/$script" "$@" --json "$json" >>"$LOG" 2>&1
    rc=$?
    local summary
    summary=$(python3 - "$json" <<'PY' 2>/dev/null
import json, sys
try:
    data = json.load(open(sys.argv[1]))
except Exception:
    sys.exit(0)
bad = [c for c in data["checks"] if c["status"] in ("FAIL", "ERROR")]
if bad:
    print("; ".join("%s: %s" % (c["name"], c["message"]) for c in bad[:3]))
else:
    print("; ".join(c["message"] for c in data["checks"]
                    if c["name"].endswith("digest") or "concordance" in c["name"]))
PY
)
    case $rc in
        0) _record_check "$name" PASS "${summary:-matches golden}" ;;
        1) _record_check "$name" FAIL "${summary:-differs from golden} (detail: $json)" ;;
        *) _record_check "$name" ERROR "${summary:-checker could not run} (detail: $json)" ;;
    esac
}

# check_bam <name> <product> [extra bam_compare.py args...]
check_bam() {
    local name=$1 product=$2; shift 2
    _golden_ready "$name" "$product" || return 0
    # samtools stats uses 16 threads; baseline/test also run in parallel inside
    # bam_compare.py. Cap at $THREADS when the host is smaller.
    local bam_threads=16
    if [ "${THREADS:-16}" -lt 16 ] 2>/dev/null; then
        bam_threads=$THREADS
    fi
    _run_checker "$name" bam_compare.py -b "$(_golden_path "$product")" -t "$product" \
        --threads "$bam_threads" "$@"
}

# check_vcf <name> <product> [extra vcf_compare.py args...]
check_vcf() {
    local name=$1 product=$2; shift 2
    _golden_ready "$name" "$product" || return 0
    _run_checker "$name" vcf_compare.py -b "$(_golden_path "$product")" -t "$product" "$@"
}

_samtools() {
    local s=${SAMTOOLS:-}
    if [ -z "$s" ] && [ -f "$SUITE_ROOT/env/tools.conf" ]; then
        s=$(sed -n 's/^[[:space:]]*SAMTOOLS[[:space:]]*=[[:space:]]*//p' \
            "$SUITE_ROOT/env/tools.conf" | tr -d '"'\''' | head -1)
    fi
    [ -z "$s" ] && s=$(command -v samtools 2>/dev/null)
    [ -n "$s" ] && [ -x "$s" ] && { echo "$s"; return 0; }
    return 1
}

# check_base_quality <name> <bam> [threshold]
#
# variantCaller drops every base scoring below -Q (default 18) before it builds
# assembly graphs. Data that sits entirely under that floor still aligns, still
# runs to completion and still writes a full-length gVCF -- the gVCF is simply
# almost all <NON_REF>, and every step reports success. wgsim is the usual
# source: it emits a flat Phred 17, exactly one below the default, which is
# enough to reduce a 3x chr20 callset to a handful of QUAL~0.006 noise sites.
# Comparing such a run against a golden built the same way passes forever.
check_base_quality() {
    local name=$1 bam=$2 threshold=${3:-18} st stats maxq frac
    st=$(_samtools) || { _record_check "$name" SKIP "samtools not found (see env/tools.conf)"; return; }
    stats=$("$st" view "$bam" 2>/dev/null | head -2000 | awk -v T="$threshold" '
        BEGIN { for (i = 33; i < 127; i++) ord[sprintf("%c", i)] = i }
        {
            n = length($11)
            for (i = 1; i <= n; i++) {
                v = ord[substr($11, i, 1)] - 33
                total++
                if (v >= T) ok++
                if (v > max) max = v
            }
        }
        END { printf "%d %.4f", max, (total ? ok / total : -1) }')
    read -r maxq frac <<< "$stats"
    local pct
    pct=$(awk -v f="$frac" 'BEGIN{printf "%.1f%%", f*100}')
    local bam_base
    bam_base=$(basename "$bam")
    local detail="{\"bam\":\"${bam_base}\",\"threshold\":${threshold},\"max_phred\":${maxq},\"qualified_ratio\":${frac},\"sampled_reads\":2000}"
    if [ "$frac" = "-1.0000" ]; then
        _record_check "$name" ERROR "could not read base qualities from $(basename "$bam")" "$detail"
    elif [ "$maxq" -lt "$threshold" ]; then
        # Simulated data (wgsim etc.) often emits a flat Phred just under the
        # default -Q floor -- a property of the simulator, not a build defect.
        # Warn rather than fail so a small simulated dataset does not turn the
        # whole case red; the near-empty callset is still caught by vcf.nonempty.
        _record_check "$name" WARN \
            "max Phred $maxq is below -Q $threshold (simulated data is often low-quality); variantCaller may emit a near-empty callset" \
            "$detail"
    elif awk -v f="$frac" 'BEGIN { exit !(f < 0.5) }'; then
        _record_check "$name" WARN \
            "only $pct of bases reach -Q $threshold (max Phred $maxq)" "$detail"
    else
        _record_check "$name" PASS \
            "$pct of bases reach -Q $threshold (max Phred $maxq)" "$detail"
    fi
}

# check_cmd <name> <message> <command...> -- generic pass/fail assertion.
check_cmd() {
    local name=$1 message=$2; shift 2
    if "$@" >>"$LOG" 2>&1; then
        _record_check "$name" PASS "$message"
    else
        _record_check "$name" FAIL "$message"
    fi
}

# Collect "tool=version" lines for the four main tools. Used by the golden
# provenance writer (PROVENANCE.txt).
_collect_versions() {
    local tool ver
    for tool in aligner bqsr variantCaller genotyper; do
        ver=$("$DCS" "$tool" --help 2>&1 | grep -oiE 'version:? *v?[0-9]+\.[0-9]+\.[0-9]+' | head -1)
        ver=${ver#*[Vv]ersion}          # strip the "version"/"Version" label
        ver=$(printf '%s' "$ver" | tr -d ' :')
        printf '%s=%s\n' "$tool" "$ver"
    done
}

_write_provenance() {
    {
        echo "case:      $CASE_NAME"
        echo "generated: $(date -Is)"
        echo "host:      $(hostname)"
        echo "dcs_home:  $DCS_HOME"
        echo "data_root: $DATA_ROOT"
        echo "threads:   $THREADS"
        echo "versions:"
        local line tool ver
        while IFS='=' read -r tool ver; do
            [ -n "$tool" ] && printf '  %-14s version: %s\n' "$tool" "$ver"
        done < <(_collect_versions)
    } > "$GOLDEN/PROVENANCE.txt"
}

_join_lines() {
    # Newline-join an array that may be empty, which "${arr[*]}" cannot do
    # safely under set -u.
    [ "$#" -eq 0 ] && return 0
    printf '%s\n' "$@"
}

_write_report() {
    local status="PASS"
    [ "$_skipped" = 1 ] && status="SKIP"
    [ "$_fail_count" -gt 0 ] && status="FAIL"
    [ "$_error_count" -gt 0 ] && status="ERROR"
    python3 - "$REPORT_DIR/case.json" "$CASE_NAME" "$RUN_ID" "$status" \
        "$(awk -v s="$CASE_START" -v e="$(_now)" 'BEGIN{printf "%.1f", e-s}')" \
        "$(_join_lines ${_checks[@]+"${_checks[@]}"})" \
        "$(_join_lines ${_steps[@]+"${_steps[@]}"})" \
        "$GOLDEN/case.json" <<'PY'
import json, os, sys
out, name, run_id, status, duration, checks_raw, steps_raw, golden_path = sys.argv[1:9]
checks = []
for line in checks_raw.splitlines():
    if not line:
        continue
    parts = line.split("|", 3)
    c = dict(zip(("name", "status", "message"), parts[:3]))
    if len(parts) > 3 and parts[3]:
        try:
            c["detail"] = json.loads(parts[3])
        except (ValueError, TypeError):
            c["detail"] = {"raw": parts[3]}
    else:
        c["detail"] = {}
    checks.append(c)
steps = [dict(zip(("name", "exit_code", "wall_s", "peak_rss_mb", "threads"), line.split("|")))
         for line in steps_raw.splitlines() if line]
for step in steps:
    step["exit_code"] = int(step["exit_code"])
    step["wall_s"] = float(step["wall_s"])
    step["peak_rss_mb"] = int(step["peak_rss_mb"])
    try:
        step["threads"] = int(step.get("threads", 0)) or 0
    except (ValueError, TypeError):
        step["threads"] = 0
# Embed the golden baseline's step timings so the HTML report can render a
# side-by-side test-vs-golden table without knowing the golden path. The
# golden case.json is written by --update-golden; if absent, golden_steps is
# empty and the renderer just shows the test side.
golden_steps = []
gp = golden_path
try:
    if gp and os.path.isfile(gp):
        gd = json.load(open(gp))
        golden_steps = gd.get("steps", [])
except Exception:
    pass
json.dump({"tool": "case", "schema_version": 1, "case": name, "run_id": run_id,
           "status": status, "duration_s": float(duration),
           "checks": checks, "steps": steps, "golden_steps": golden_steps},
          open(out, "w"), indent=2)
PY
}

case_end() {
    local status="PASS" rc=0
    [ "$_skipped" = 1 ] && status="SKIP"
    [ "$_fail_count" -gt 0 ] && { status="FAIL"; rc=1; }
    [ "$_error_count" -gt 0 ] && { status="ERROR"; rc=2; }

    _rule
    if [ ${#_steps[@]} -gt 0 ]; then
        _log "resource usage:"
        for row in "${_steps[@]}"; do
            IFS='|' read -r label _rc wall rss _threads <<< "$row"
            _log "$(printf '  %-14s %8ss  %6s MB  %s thr' "$label" "$wall" "$rss" "${_threads:-?}")"
        done
    fi
    _log ""
    _log "=> $CASE_NAME: $status   (report: $REPORT_DIR)"
    _write_report
    # On --update-golden, also snapshot this run's case.json (steps + threads +
    # checks) into golden/, so later runs can render a side-by-side test-vs-
    # golden steps table. Golden products are already copied by _golden_ready;
    # this only adds the report metadata.
    if [ "$UPDATE_GOLDEN" = 1 ] && [ -f "$REPORT_DIR/case.json" ]; then
        mkdir -p "$GOLDEN"
        cp -f "$REPORT_DIR/case.json" "$GOLDEN/case.json"
        # When refreshing golden, the baseline IS this run, so the side-by-side
        # steps table should show this run's timings on both sides. _write_report
        # ran before the golden copy, so golden_steps was empty; backfill it from
        # the steps we just wrote so a single update-golden run renders a complete
        # test-vs-golden table (a later normal run reads the real golden instead).
        python3 - "$REPORT_DIR/case.json" <<'PY'
import json, sys
p = sys.argv[1]
d = json.load(open(p))
d["golden_steps"] = d.get("steps", [])
json.dump(d, open(p, "w"), indent=2)
PY
    fi
    [ "$KEEP_WORK" = 1 ] || true
    exit $rc
}
