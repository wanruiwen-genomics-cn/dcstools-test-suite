#!/bin/bash
# Verify the aligner's --useSharedMemIndex actually shares the index.
#
# This case does NOT compare against a golden baseline. It only checks that the
# shared-memory index mechanism engages, which requires concurrency to be
# meaningful (a single process always "creates"; sharing only shows up when a
# second process "attaches"). It launches two aligners on the same small index
# at the same time and asserts:
#   1. a /dev/shm/aligner_idx_* object appears and has a plausible size, and
#   2. the logs show exactly one "created shared memory" and >=1 "attached".
# Both are reported as checks; the shm name/size/role counts are recorded.

source "$(dirname "$(readlink -f "$0")")/common.sh"

case_begin aligner_shared_mem

need_file FA  "$DATA_ROOT/ref/chr20.fa"        "reference fasta"
need_file FQ1 "$DATA_ROOT/fastq/simu.r1.fq.gz" "read 1"
need_file FQ2 "$DATA_ROOT/fastq/simu.r2.fq.gz" "read 2"
preflight_done

# The shared-memory index needs the classic BWA index sidecars present.
for file in "$FA".0123 "$FA".amb "$FA".ann "$FA".bwt "$FA".pac "$FA".sa; do
    if [ ! -f "$file" ]; then
        run_step build_index "$DCS" aligner --threads "$THREADS" --build-index "$FA"
        break
    fi
done

group="@RG\\tID:shmtest\\tLB:shmtest\\tSM:shmtest\\tPL:COMPLETE\\tCN:DCS"

# Snapshot existing aligner shm objects so we only count ones this run creates.
before=$(mktemp)
ls -1 /dev/shm/aligner_idx_* 2>/dev/null | sort > "$before" || true

mkdir -p "$WORK/p1/qc" "$WORK/p2/qc"
log1="$WORK/p1/align.log"
log2="$WORK/p2/align.log"

_log ""
_log ">>> launching 2 concurrent aligners with --useSharedMemIndex"

# Run the two aligners directly (not via run_step) so they overlap in time;
# `dcs aligner` execs the aligner, so pids are the aligner processes.
"$DCS" aligner --threads 4 --aln-index "$FA" --fq1 "$FQ1" --fq2 "$FQ2" \
    --read-group "$group" --no-filter --fastq-qc-dir "$WORK/p1/qc" \
    --bam-out "$WORK/p1/out.bam" --useSharedMemIndex > "$log1" 2>&1 &
pid1=$!
"$DCS" aligner --threads 4 --aln-index "$FA" --fq1 "$FQ1" --fq2 "$FQ2" \
    --read-group "$group" --no-filter --fastq-qc-dir "$WORK/p2/qc" \
    --bam-out "$WORK/p2/out.bam" --useSharedMemIndex > "$log2" 2>&1 &
pid2=$!

# Poll /dev/shm while both are alive to catch the shm object at its peak.
shm_seen=""
shm_size=0
for _ in $(seq 1 60); do
    cur=$(ls -1 /dev/shm/aligner_idx_* 2>/dev/null | sort)
    new=$(comm -13 "$before" <(printf '%s\n' "$cur"))
    if [ -n "$new" ]; then
        shm_seen=$(printf '%s\n' "$new" | head -1)
        shm_size=$(stat -c %s "$shm_seen" 2>/dev/null || echo 0)
        break
    fi
    # Stop early if both already exited without creating shm.
    if ! kill -0 "$pid1" 2>/dev/null && ! kill -0 "$pid2" 2>/dev/null; then
        break
    fi
    sleep 1
done

wait "$pid1"; rc1=$?
wait "$pid2"; rc2=$?
rm -f "$before"

# Give the last process a moment to unlink, then check for leftovers.
sleep 2
leftover=$(ls -1 /dev/shm/aligner_idx_* 2>/dev/null | wc -l)

created=$(grep -c "created shared memory" "$log1" "$log2" 2>/dev/null | awk -F: '{s+=$NF} END{print s+0}')
attached=$(grep -c "attached to existing shared memory" "$log1" "$log2" 2>/dev/null | awk -F: '{s+=$NF} END{print s+0}')

_log "aligner rc: p1=$rc1 p2=$rc2"
_log "shm object: ${shm_seen:-<none>} size=${shm_size}"
_log "roles: created=$created attached=$attached leftover=$leftover"

# Both aligners must finish OK.
if [ "$rc1" -eq 0 ] && [ "$rc2" -eq 0 ]; then
    _record_check shm.align_exit PASS "both aligners exited 0"
else
    _record_check shm.align_exit FAIL "aligner exit codes p1=$rc1 p2=$rc2"
fi

# A shared-memory object must have appeared with non-trivial size.
if [ -n "$shm_seen" ] && [ "$shm_size" -gt 1000000 ]; then
    _record_check shm.object_present PASS \
        "created $(basename "$shm_seen") ($shm_size bytes)" \
        "{\"shm\":\"$(basename "$shm_seen")\",\"bytes\":$shm_size}"
else
    _record_check shm.object_present FAIL \
        "no /dev/shm/aligner_idx_* object observed (size=$shm_size)"
fi

# Exactly one creator, and at least one attacher -> real sharing.
if [ "$created" -eq 1 ] && [ "$attached" -ge 1 ]; then
    _record_check shm.roles PASS "created=$created attached=$attached (index shared)" \
        "{\"created\":$created,\"attached\":$attached}"
else
    _record_check shm.roles FAIL "created=$created attached=$attached (expected created=1, attached>=1)" \
        "{\"created\":$created,\"attached\":$attached}"
fi

# No orphan left behind after both processes exited.
if [ "$leftover" -eq 0 ]; then
    _record_check shm.cleanup PASS "no leftover /dev/shm/aligner_idx_* after exit"
else
    _record_check shm.cleanup WARN "$leftover leftover shm object(s) after exit"
fi

case_end
