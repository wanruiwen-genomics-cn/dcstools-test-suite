#!/bin/bash
# Full release validation: regenerate golden with the baseline build, then
# run the release-candidate build against that golden.
# Long-running (~10h). Logs to logs/release_validation.log.
set -uo pipefail
cd "$(dirname "$(readlink -f "$0")")"

export JAVA_HOME=/data/disk2/private/wanruiwen/software/jdk-17.0.18+8
export PATH="$JAVA_HOME/bin:$PATH"
export PROXY_HOST=127.0.0.1:8909
export WGS30X_DATA=/zdswhst2/ST_BIOINTEL/P24Z32400N0004/wanruiwen/dcstools_data/HG001_30x
export WES_DATA=/zdswhst2/ST_BIOINTEL/P24Z32400N0004/wanruiwen/dcstools_data/HG001_30x
# JAVA8_HOME, SPARK_BIN, DATA_ROOT, GOLDEN_DCS_HOME, THREADS come from env/cases.conf

BASELINE=/data/disk2/private/wanruiwen/dcstools/dcs-tools
RELEASE=/data/disk2/private/wanruiwen/dcstools/dcs-archive/dcs-tools
# Cases with data on this machine (cycas/potato/wheat have no golden/data here).
ITEMS="bundle_audit cli_contract_check mutect2_tn wgs_human_small wes_human jointcall seqarc_roundtrip wgs_human_30x"

LOG=logs/release_validation.log
mkdir -p logs
: > "$LOG"

log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

log "===== RELEASE VALIDATION START ====="
log "baseline  (golden build) : $BASELINE"
log "release   (candidate)     : $RELEASE"
log "items: $ITEMS"
log "disk free: $(df -h /data/disk2 | awk 'NR==2{print $4}')"

log "===== STAGE 1: UPDATE GOLDEN with baseline build ====="
DCS_HOME="$BASELINE" ./run.sh --update-golden $ITEMS 2>&1 | tee -a "$LOG"
rc1=${PIPESTATUS[0]}
log "STAGE 1 rc=$rc1"
log "disk free after stage 1: $(df -h /data/disk2 | awk 'NR==2{print $4}')"

log "===== STAGE 2: TEST with release build ====="
DCS_HOME="$RELEASE" ./run.sh $ITEMS 2>&1 | tee -a "$LOG"
rc2=${PIPESTATUS[0]}
log "STAGE 2 rc=$rc2"

log "===== DONE  stage1=$rc1  stage2=$rc2 ====="
