#!/bin/bash
# Run the release validation suite and reduce it to one verdict.
#
#   ./run.sh                              the smoke tier
#   ./run.sh --tier full                  every case
#   ./run.sh --tier heavy                 only the large/ploidy cases
#   ./run.sh --tier static                packaging and CLI checks only, no data
#   ./run.sh wgs_human_small mutect2_tn   just these items
#   ./run.sh --tier full --update-golden  install this run as the baseline
#   ./run.sh --list
#
# Items run in the listed order and a failure does not stop the run: the point
# of a suite is to come back with the whole picture, not the first problem.
# Use --fail-fast when bisecting.
#
# Exit codes: 0 everything passed, 1 something failed, 2 something could not run.

set -uo pipefail

SUITE_ROOT=$(cd "$(dirname "$(readlink -f "$0")")" && pwd)
CASES_DIR="$SUITE_ROOT/cases"
SCRIPTS_DIR="$SUITE_ROOT/scripts"

[ -f "$SUITE_ROOT/env/cases.conf" ] && . "$SUITE_ROOT/env/cases.conf"
export DCS_HOME PROXY_HOST

# Tiers, in execution order. An item is "<kind>:<name>":
#   check:<stem>  a python checker in scripts/
#   case:<stem>   a case script in cases/
#
# Cheap items come first so an obviously broken bundle is reported in seconds
# rather than after the alignment runs. The static pair needs no test data at
# all, which is also what makes it the right gate for a build pipeline.
TIER_static="check:bundle_audit check:cli_contract_check"
# The heavy cases need real, large datasets and run for minutes to hours:
# 30x human WGS and the potato polyploid (jointcaller/Spark). Cycas and wheat
# stay in cases/ but are left out of heavy/full until their data and goldens
# are ready. They are excluded from the smoke tier so a developer running
# `make smoke` is not blocked on data they may not have, but `make heavy` and
# `make full` run the rest as a pre-release gate.
TIER_heavy_cases="case:wgs_human_30x case:wgs_human_30x_ert case:wgs_human_30x_minibwa \
case:wgs_rice case:wgs_rice_ert case:wgs_rice_minibwa \
case:wgs_mouse case:wgs_mouse_ert case:wgs_mouse_minibwa \
case:wgs_potato case:wgs_potato_ert case:wgs_potato_minibwa \
case:potato_polyploid"
TIER_smoke="$TIER_static case:wgs_human_small"
TIER_heavy="$TIER_static $TIER_heavy_cases"
TIER_full="$TIER_static case:mutect2_tn case:wgs_human_small case:wes_human \
           case:jointcall case:seqarc_roundtrip $TIER_heavy_cases"

TIER=smoke
FAIL_FAST=0
DO_LIST=0
PASSTHRU=()
SELECTED=()

while [ $# -gt 0 ]; do
    case "$1" in
        --tier)   TIER=${2:-}; shift 2 ;;
        --tier=*) TIER=${1#*=}; shift ;;
        --list)   DO_LIST=1; shift ;;
        --fail-fast) FAIL_FAST=1; shift ;;
        --update-golden|--keep-work) PASSTHRU+=("$1"); shift ;;
        -h|--help)
            sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'
            exit 0 ;;
        -*) echo "unknown option: $1" >&2; exit 2 ;;
        *)  SELECTED+=("$1"); shift ;;
    esac
done

if [ "$DO_LIST" = 1 ]; then
    echo "tiers:"
    for tier in static smoke full heavy; do
        var="TIER_$tier"
        printf '  %-8s %s\n' "$tier" "$(echo ${!var} | sed 's/check://g; s/case://g')"
    done
    echo "cases:"
    for file in "$CASES_DIR"/*.sh; do
        name=$(basename "$file" .sh)
        [ "$name" = common ] && continue
        printf '  %s\n' "$name"
    done
    exit 0
fi

# An explicit item list beats the tier; resolve each name to its kind so the
# caller does not have to remember which of the two directories it lives in.
if [ ${#SELECTED[@]} -gt 0 ]; then
    items=()
    for name in "${SELECTED[@]}"; do
        if [ -f "$CASES_DIR/$name.sh" ]; then
            items+=("case:$name")
        elif [ -f "$SCRIPTS_DIR/$name.py" ]; then
            items+=("check:$name")
        else
            echo "unknown item: $name  (try ./run.sh --list)" >&2
            exit 2
        fi
    done
else
    var="TIER_$TIER"
    items=(${!var-})
    if [ ${#items[@]} -eq 0 ]; then
        echo "unknown tier: $TIER  (have: static smoke full heavy)" >&2
        exit 2
    fi
fi

export SUITE_RUN_ID=${SUITE_RUN_ID:-$(date +%Y%m%d_%H%M%S)}
SUITE_DIR="$REPORT_ROOT/suite/$SUITE_RUN_ID"
mkdir -p "$SUITE_DIR"

echo "suite run $SUITE_RUN_ID   tier=$TIER   items=${#items[@]}"
echo "dcs_home  $DCS_HOME"
echo "data_root $DATA_ROOT"
echo "reports   $SUITE_DIR"

# The checkers take different arguments; keeping that mapping here rather than
# in the tier table means adding a checker does not mean touching three places.
_check_args() {
    case "$1" in
        bundle_audit)
            echo "--dcs-home|$DCS_HOME" ;;
        cli_contract_check)
            # Scan the cases too: they are the scripts whose flags must keep
            # working, and they are where a removed flag would bite first.
            echo "--dcs-home|$DCS_HOME|--scan-dir|$CASES_DIR" ;;
    esac
}

records=()
worst=0

for item in "${items[@]}"; do
    kind=${item%%:*}
    name=${item#*:}
    echo
    echo "### $name"
    start=$(date +%s.%N)

    if [ "$kind" = check ]; then
        json="$SUITE_DIR/$name.json"
        IFS='|' read -r -a extra <<< "$(_check_args "$name")"
        python3 "$SCRIPTS_DIR/$name.py" ${extra[@]+"${extra[@]}"} \
            --json "$json" 2>&1 | tee "$SUITE_DIR/$name.out"
        rc=${PIPESTATUS[0]}
    else
        json="$REPORT_ROOT/$name/$SUITE_RUN_ID/case.json"
        bash "$CASES_DIR/$name.sh" ${PASSTHRU[@]+"${PASSTHRU[@]}"}
        rc=$?
    fi

    wall=$(awk -v s="$start" -v e="$(date +%s.%N)" 'BEGIN{printf "%.1f", e-s}')
    records+=("--item" "$name|$kind|$rc|$wall|$json")
    [ "$rc" -gt "$worst" ] && worst=$rc

    if [ "$rc" -ne 0 ] && [ "$FAIL_FAST" = 1 ]; then
        echo "stopping: $name exited $rc and --fail-fast is set" >&2
        break
    fi
done

python3 "$SCRIPTS_DIR/suite_report.py" \
    --run-id "$SUITE_RUN_ID" --tier "$TIER" \
    --out "$SUITE_DIR/suite.json" "${records[@]}"

# Generate the HTML report from the JSON we just wrote. Best-effort: a failure
# here does not change the exit code, which is owned by suite_report.py.
if [ -f "$SUITE_DIR/suite.json" ]; then
    python3 "$SCRIPTS_DIR/html_report.py" "$SUITE_DIR/suite.json" \
        -o "$SUITE_DIR/report.html" 2>/dev/null
    # Refresh the top-level entry page so reports/index.html always reflects
    # the latest state of every case across all runs. Also best-effort.
    python3 "$SCRIPTS_DIR/index_report.py" \
        --reports-dir "$REPORT_ROOT" -o "$REPORT_ROOT/index.html" 2>/dev/null || true
fi

exit $?
