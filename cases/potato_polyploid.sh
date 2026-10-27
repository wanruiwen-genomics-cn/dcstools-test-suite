#!/bin/bash
# Polyploid variant calling on potato (土豆), an autotetraploid (4n).
#
# Every other case in the suite is diploid: variantCaller defaults to ploidy 2
# and genotyper calls diploid genotypes. This one exercises the path
# fastq2vcf.sh takes when ploidy > 2: variantCaller gets --ploidy 4, and
# genotyping goes through jointcaller (DPGT on Spark) rather than genotyper,
# because genotyper cannot emit polyploid GTs. A regression that breaks
# --ploidy, or that silently falls back to diploid, would pass every other
# case and fail only here.
#
# The case needs JDK 8 and spark-submit for the jointcaller step. It skips
# with a clear reason when those are missing, the same way jointcall does.
#
#   cases/potato_polyploid.sh
#   cases/potato_polyploid.sh --update-golden
#   POTATO_DATA=/data/potato POTATO_SAMPLE=tuber03 PLOIDY=4 cases/potato_polyploid.sh
#   JAVA8_HOME=/path/to/jdk8 SPARK_BIN=/path/to/spark/bin cases/potato_polyploid.sh

source "$(dirname "$(readlink -f "$0")")/common.sh"

case_begin potato_polyploid

data=${POTATO_DATA:-$DATA_ROOT/potato}
sample=${POTATO_SAMPLE:-potato}
ploidy=${PLOIDY:-4}
group="@RG\\tID:${sample}\\tLB:${sample}\\tSM:${sample}\\tPL:COMPLETE\\tCN:DCS"

need_file FA  "$data/ref/potato.fa"              "potato reference fasta"
need_file FQ1 "$data/fastq/${sample}.r1.fq.gz"   "read 1"
need_file FQ2 "$data/fastq/${sample}.r2.fq.gz"   "read 2"
# Known sites are optional.
known_sites=()
for vcf in "$data"/known_sites/*.vcf.gz; do
    [ -f "$vcf" ] && known_sites+=("$vcf")
done

# Environment for the jointcaller step: same preflight as jointcall.
if [ -n "$JAVA8_HOME" ]; then
    export JAVA_HOME="$JAVA8_HOME"
fi
[ -n "$SPARK_BIN" ] && export PATH="$SPARK_BIN:$PATH"

if [ -z "${JAVA_HOME:-}" ]; then
    skip_case env.java "JAVA_HOME unset; set JAVA8_HOME in env/cases.conf for polyploid genotyping"
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
note data.ploidy "$ploidy"

# Step 0: build the alignment index if this dataset has not been indexed yet.
for file in "$FA".0123 "$FA".amb "$FA".ann "$FA".bwt "$FA".index1 "$FA".index2 "$FA".pac "$FA".sa; do
    if [ ! -f "$file" ]; then
        run_step build_index "$DCS" aligner --threads "$THREADS" --build-index "$FA"
        break
    fi
done

mkdir -p "$WORK/${sample}.qcdir"

# Step 1: alignment, sorting and duplicate marking in one pass.
run_step align "$DCS" aligner \
    --threads "$THREADS" \
    --aln-index "$FA" \
    --fq1 "$FQ1" --fq2 "$FQ2" \
    --read-group "$group" \
    --mark-split-hits-secondary \
    --use-soft-clipping-for-sup-align \
    --fastq-qc-dir "$WORK/${sample}.qcdir" \
    --bam-out "$WORK/${sample}.align.bam"

check_base_quality data.base_quality "$WORK/${sample}.align.bam"

# Step 2: BQSR (only if known sites were provided).
recal_arg=""
if [ ${#known_sites[@]} -gt 0 ]; then
    ks_args=()
    for vcf in "${known_sites[@]}"; do
        ks_args+=(--known-sites "$vcf")
    done
    run_step bqsr "$DCS" bqsr \
        --threads "$THREADS" \
        --in-bam "$WORK/${sample}.align.bam" \
        --reference "$FA" \
        "${ks_args[@]}" \
        --recal-table "$WORK/${sample}.recal.table"
    recal_arg="--recal-table $WORK/${sample}.recal.table"
fi

# Step 3: variant calling with --ploidy. This is the flag under test. The
# gVCF that comes out carries polyploid likelihoods; genotyper would collapse
# them to diploid, so the next step uses jointcaller instead.
run_step variant_call "$DCS" variantCaller \
    --threads "$THREADS" \
    --input "$WORK/${sample}.align.bam" \
    --output "$WORK/${sample}.variantCaller.g.vcf.gz" \
    --reference "$FA" \
    --ploidy "$ploidy" \
    $recal_arg

# Step 4: genotyping via jointcaller (DPGT on Spark). fastq2vcf.sh routes
# ploidy > 2 here because genotyper cannot emit polyploid GTs. The input is a
# FOFN with one gVCF line; the same interface as the multi-sample jointcall
# case, but here it is a single sample exercising the polyploid genotype path.
GVCF_FOFN="$WORK/gvcf.fofn"
echo "$WORK/${sample}.variantCaller.g.vcf.gz" > "$GVCF_FOFN"

OUT_DIR="$WORK/result"
run_step genotype "$DCS" jointcaller \
    --input "$GVCF_FOFN" \
    --output "$OUT_DIR" \
    --reference "$FA" \
    --jobs "$THREADS" \
    --memory "${JOINTCALL_MEMORY:-4g}"

# --- verification -----------------------------------------------------------

# DPGT names its output by region; find it rather than assuming a fixed name.
joint_vcf=$(find "$OUT_DIR" -name "*.vcf.gz" | sort | head -1)
if [ -z "$joint_vcf" ]; then
    _record_check product.joint_vcf FAIL "no vcf.gz under $OUT_DIR"
    case_end
fi
_record_check product.joint_vcf PASS "$(basename "$joint_vcf"), $(stat -c%s "$joint_vcf") bytes"

expect_file "$WORK/${sample}.align.bam"                  1000
expect_file "$WORK/${sample}.align.bam.bai"              100
expect_file "$WORK/${sample}.variantCaller.g.vcf.gz"    1000

check_cmd joint.nonempty "at least one variant emitted" \
    bash -c 'test "$(zcat "$1" | grep -vc "^#")" -gt 0' _ "$joint_vcf"

# The assertion that makes this case distinct from every other in the suite:
# the genotypes must be polyploid. A GT field with exactly $ploidy components
# (e.g. 0/0/0/1 for tetraploid) proves --ploidy was honoured end to end. A
# release that silently drops to diploid would emit 0/0, 0/1, 1/1 -- three
# components or fewer -- and this check catches it.
check_cmd genotype.ploidy "at least one genotype has $ploidy allele components" \
    bash -c 'zcat "$1" | grep -v "^#" | cut -f10 | tr "," "\n" | grep -E "^[0-9](/[0-9]){$((ploidy-1))}$" | head -1 | grep -q .' \
    _ "$joint_vcf"

cp -f "$joint_vcf" "$WORK/potato_polyploid.vcf.gz"
[ -f "$joint_vcf.tbi" ] && cp -f "$joint_vcf.tbi" "$WORK/potato_polyploid.vcf.gz.tbi"

check_bam align_bam "$WORK/${sample}.align.bam"
check_vcf gvcf      "$WORK/${sample}.variantCaller.g.vcf.gz"
check_vcf joint_vcf "$WORK/potato_polyploid.vcf.gz"

case_end
