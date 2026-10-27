#!/bin/bash
# Joint genotyping of many gVCFs through dcs jointcaller (DPGT on Spark).
#
# Two things make this case different from the rest of the suite. It needs an
# external environment -- JDK 8 and spark-submit -- which is itself worth
# asserting, since a wrong JDK is the most common reason this step fails. And it
# is the path a polyploid sample takes: fastq2vcf.sh routes ploidy > 2 to
# jointcaller instead of genotyper, so this code is on the critical path for the
# potato dataset even for a single sample.
#
#   cases/jointcall.sh [--update-golden]
#   JAVA8_HOME=/path/to/jdk8 SPARK_BIN=/path/to/spark/bin cases/jointcall.sh

source "$(dirname "$(readlink -f "$0")")/common.sh"

case_begin jointcall

region=${REGION:-chr20:1000000-2000000}
memory=${JOINTCALL_MEMORY:-8g}

need_dir  GVCF_DIR "$DATA_ROOT/gvcf"     "gvcf directory"
need_file FA       "$DATA_ROOT/ref/chr20.fa" "reference fasta"

# Environment first: skipping with a clear reason beats a Java stack trace.
if [ -n "$JAVA8_HOME" ]; then
    export JAVA_HOME="$JAVA8_HOME"
fi
[ -n "$SPARK_BIN" ] && export PATH="$SPARK_BIN:$PATH"

if [ -z "${JAVA_HOME:-}" ]; then
    skip_case env.java "JAVA_HOME unset; set JAVA8_HOME in env/cases.conf"
fi
java_version=$("$JAVA_HOME/bin/java" -version 2>&1 | head -1)
case "$java_version" in
    *'"1.8'*|*'"8.'*) ;;
    *) skip_case env.java "jointcaller requires JDK 8, found: $java_version" ;;
esac
command -v spark-submit >/dev/null 2>&1 || \
    skip_case env.spark "spark-submit not on PATH; set SPARK_BIN in env/cases.conf"

preflight_done
note env.java "$java_version"

# Step 1: the input is a file of gVCF paths, one per line.
GVCF_FOFN="$WORK/gvcf.fofn"
ls "$GVCF_DIR"/*.g.vcf.gz > "$GVCF_FOFN" 2>/dev/null || ls "$GVCF_DIR"/*.gz > "$GVCF_FOFN"
sample_count=$(wc -l < "$GVCF_FOFN")
note input.samples "$sample_count gVCF(s) joint called over $region"

check_cmd input.fofn "gvcf list is not empty" test "$sample_count" -gt 0

# Every input needs its index, otherwise DPGT fails deep inside a Spark task
# where the message is almost unreadable.
check_cmd input.indexed "every input gVCF has a .tbi/.csi" \
    bash -c 'while read -r f; do [ -f "$f.tbi" ] || [ -f "$f.csi" ] || exit 1; done < "$1"' _ "$GVCF_FOFN"

# Step 2: joint calling.
OUT_DIR="$WORK/result"
run_step jointcall "$DCS" jointcaller \
    --input "$GVCF_FOFN" \
    --output "$OUT_DIR" \
    --reference "$FA" \
    --jobs "$THREADS" \
    --memory "$memory" \
    -l "$region"

# --- verification -----------------------------------------------------------

# DPGT names its output by region, so find it rather than assuming.
joint_vcf=$(find "$OUT_DIR" -name "*.vcf.gz" | sort | head -1)
if [ -z "$joint_vcf" ]; then
    _record_check product.joint_vcf FAIL "no vcf.gz under $OUT_DIR"
    case_end
fi
_record_check product.joint_vcf PASS "$(basename "$joint_vcf"), $(stat -c%s "$joint_vcf") bytes"

# The joint VCF must carry every input sample; a silently dropped sample is the
# classic joint-calling regression and record counts alone will not reveal it.
check_cmd joint.sample_count "output has all $sample_count samples" \
    bash -c 'test "$(zcat "$1" | grep -m1 "^#CHROM" | cut -f10- | wc -w)" -eq "$2"' \
    _ "$joint_vcf" "$sample_count"

check_cmd joint.nonempty "at least one variant emitted" \
    bash -c 'test "$(zcat "$1" | grep -vc "^#")" -gt 0' _ "$joint_vcf"

cp -f "$joint_vcf" "$WORK/jointcall.vcf.gz"
[ -f "$joint_vcf.tbi" ] && cp -f "$joint_vcf.tbi" "$WORK/jointcall.vcf.gz.tbi"
check_vcf joint_vcf "$WORK/jointcall.vcf.gz"

case_end
