#!/bin/bash
# Germline WGS on a small human dataset: the full fastq -> vcf path.
#
# Covers aligner (default index), bqsr with four known-sites files, variantCaller
# in gVCF mode, genotyper, and the stat/report tools. This is the case that runs
# on every release; the large genomes and the polyploid case reuse the same shape
# with different inputs.
#
#   cases/wgs_human_small.sh                  compare against golden
#   cases/wgs_human_small.sh --update-golden  install this run as the baseline

source "$(dirname "$(readlink -f "$0")")/common.sh"

case_begin wgs_human_small

sample=simu
group="@RG\\tID:${sample}\\tLB:${sample}\\tSM:${sample}\\tPL:COMPLETE\\tCN:DCS"

need_file FA  "$DATA_ROOT/ref/chr20.fa"                   "reference fasta"
need_file FQ1 "$DATA_ROOT/fastq/simu.r1.fq.gz"            "read 1"
need_file FQ2 "$DATA_ROOT/fastq/simu.r2.fq.gz"            "read 2"
need_file KN1 "$DATA_ROOT/known_sites/1000G_phase1.snps.high_confidence.hg38.chr20.vcf.gz" "known sites"
need_file KN2 "$DATA_ROOT/known_sites/dbsnp_151.hg38.chr20.vcf.gz"                         "dbsnp"
need_file KN3 "$DATA_ROOT/known_sites/hapmap_3.3.hg38.chr20.vcf.gz"                        "known sites"
need_file KN4 "$DATA_ROOT/known_sites/Mills_and_1000G_gold_standard.indels.hg38.chr20.vcf.gz" "known sites"
preflight_done

if [ "${RESUME_VC:-0}" = 0 ]; then
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

# Qualities under variantCaller's -Q floor survive alignment and bqsr untouched
# -- recalibration does not lift them back over the threshold -- and then yield a
# full-length gVCF that is almost entirely <NON_REF> with every step exiting 0.
check_base_quality data.base_quality "$WORK/${sample}.align.bam"

# Step 2: base quality score recalibration.
run_step bqsr "$DCS" bqsr \
    --threads "$THREADS" \
    --in-bam "$WORK/${sample}.align.bam" \
    --reference "$FA" \
    --known-sites "$KN1" --known-sites "$KN2" \
    --known-sites "$KN3" --known-sites "$KN4" \
    --recal-table "$WORK/${sample}.recal.table"
else
    # RESUME_VC: reuse existing align BAM + recal table, only re-run vc + downstream.
    [ -f "$WORK/${sample}.recal.table" ] || { echo "RESUME_VC: missing recal.table" >&2; exit 1; }
fi

# Step 3: variant calling into a gVCF.
run_step variant_call "$DCS" variantCaller \
    --threads "$THREADS" \
    --input "$WORK/${sample}.align.bam" \
    --output "$WORK/${sample}.variantCaller.g.vcf.gz" \
    --reference "$FA" \
    --recal-table "$WORK/${sample}.recal.table"

# Step 4: single sample genotyping. Pinned to GENOTYPER_THREADS (default 1, see
# env/cases.conf) because genotyper's QD jitter RNG is a shared singleton and
# multi-threaded runs assign the jittered values in scheduling order.
run_step genotype "$DCS" genotyper \
    --threads "$GENOTYPER_THREADS" \
    --input "$WORK/${sample}.variantCaller.g.vcf.gz" \
    --output "$WORK/${sample}.genotyper.vcf.gz" \
    --dbsnp "$KN2" -s 30.0

# Step 5: statistics and the customer-facing report.
# bam-stat requires --source: the reference for WGS, the target BED for WES.
# Note that quick_start_wgs.sh (and fastq2vcf.sh on the no-BED path) hands
# bamStat.txt to `report` without ever creating it, so generate it here.
run_step bam_stat "$DCS" bam-stat \
    --bam "$WORK/${sample}.align.bam" \
    --source "$FA" \
    --output "$WORK/bamStat.txt"
run_step vcf_stat bash -c "\"$DCS\" vcf-stat --vcf \"$WORK/${sample}.genotyper.vcf.gz\" > \"$WORK/vcfStat.txt\""
run_step report "$DCS" report \
    --sample "$sample" \
    --filter "$WORK/${sample}.qcdir" \
    --bam_stat "$WORK/bamStat.txt" \
    --vcf_stat "$WORK/vcfStat.txt" \
    --output "$WORK/report"

# --- verification -----------------------------------------------------------

expect_file "$WORK/${sample}.align.bam"              1000
expect_file "$WORK/${sample}.align.bam.bai"          100
expect_file "$WORK/${sample}.recal.table"            100
expect_file "$WORK/${sample}.variantCaller.g.vcf.gz" 1000
expect_file "$WORK/${sample}.genotyper.vcf.gz"       1000
expect_file "$WORK/${sample}.genotyper.vcf.gz.tbi"   50
expect_file "$WORK/bamStat.txt"
expect_file "$WORK/vcfStat.txt"
check_report "$WORK/report"

check_bam align_bam  "$WORK/${sample}.align.bam"
check_vcf gvcf       "$WORK/${sample}.variantCaller.g.vcf.gz"
check_vcf vcf        "$WORK/${sample}.genotyper.vcf.gz"

case_end
