#!/bin/bash
# Germline WGS at realistic 30x coverage on human data.
#
# wgs_human_small uses 3x wgsim data whose base qualities sit at Phred 17 --
# exactly one below variantCaller's default -Q floor -- so the pipeline runs
# green end to end and emits a near-empty gVCF. That case guards the plumbing.
# This one guards the biology: at 30x with real Illumina qualities the caller
# must actually find variants. A release that silently drops the callset to
# noise would pass wgs_human_small and fail here.
#
# BQSR is optional: pass known sites and the case recalibrates; leave them
# unset and it skips straight to calling, which is the path a non-model
# organism or a re-analysis without known sites takes.
#
#   cases/wgs_human_30x.sh
#   cases/wgs_human_30x.sh --update-golden
#   WGS30X_DATA=/data/wgs30x cases/wgs_human_30x.sh

source "$(dirname "$(readlink -f "$0")")/common.sh"

case_begin wgs_human_30x

# Per-dataset subdirectory so several genomes can live under one DATA_ROOT
# without colliding on ref/ or fastq/. Override with WGS30X_DATA for a custom
# layout outside DATA_ROOT.
data=${WGS30X_DATA:-$DATA_ROOT/wgs30x}
sample=${WGS30X_SAMPLE:-sample}
group="@RG\\tID:${sample}\\tLB:${sample}\\tSM:${sample}\\tPL:COMPLETE\\tCN:DCS"

need_file FA  "$data/ref/hg38.fa"          "reference fasta"
need_file FQ1 "$data/fastq/${sample}.r1.fq.gz" "read 1"
need_file FQ2 "$data/fastq/${sample}.r2.fq.gz" "read 2"
# Known sites are optional -- the case runs BQSR only if at least one is present.
known_sites=()
for vcf in "$data"/known_sites/*.vcf.gz; do
    [ -f "$vcf" ] && known_sites+=("$vcf")
done

preflight_done

# Step 0: build the alignment index if this dataset has not been indexed yet.
# A 30x human run is usually whole-genome, so the index is large and slow to
# build; checking first avoids rebuilding it on every run. The .bwt.2bit.64
# file is the ERT index the aligner loads even in default mode; .index1/.index2
# are the legacy large index and may be absent when only the ERT index was built.
for file in "$FA".0123 "$FA".amb "$FA".ann "$FA".bwt "$FA".pac "$FA".sa "$FA".bwt.2bit.64; do
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

# The guard that wgs_human_small cannot provide: at 30x with real qualities
# the caller has enough signal to work. If this fails the data is not what the
# case was built for, and every downstream check is meaningless.
check_base_quality data.base_quality "$WORK/${sample}.align.bam"

# Step 2: base quality score recalibration (only if known sites were provided).
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

# Step 3: variant calling into a gVCF.
run_step variant_call "$DCS" variantCaller \
    --threads "$THREADS" \
    --input "$WORK/${sample}.align.bam" \
    --output "$WORK/${sample}.variantCaller.g.vcf.gz" \
    --reference "$FA" \
    $recal_arg

# Step 4: single sample genotyping. Pinned to GENOTYPER_THREADS (default 1, see
# env/cases.conf) because genotyper's QD jitter RNG is a shared singleton and
# multi-threaded runs assign the jittered values in scheduling order.
run_step genotype "$DCS" genotyper \
    --threads "$GENOTYPER_THREADS" \
    --input "$WORK/${sample}.variantCaller.g.vcf.gz" \
    --output "$WORK/${sample}.genotyper.vcf.gz" \
    -s 30.0

# Step 5: statistics and the customer-facing report.
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

expect_file "$WORK/${sample}.align.bam"                  1000
expect_file "$WORK/${sample}.align.bam.bai"              100
expect_file "$WORK/${sample}.variantCaller.g.vcf.gz"    1000
expect_file "$WORK/${sample}.genotyper.vcf.gz"           1000
expect_file "$WORK/${sample}.genotyper.vcf.gz.tbi"        50
expect_file "$WORK/bamStat.txt"
expect_file "$WORK/vcfStat.txt"
check_report "$WORK/report"

# Variant count vs golden: display-only comparison (no absolute threshold).
# Real callset regressions are caught by check_vcf / vcf_compare.
variant_count=$(zcat "$WORK/${sample}.genotyper.vcf.gz" 2>/dev/null | grep -vc '^#' || echo 0)
golden_vcf="$GOLDEN/${sample}.genotyper.vcf.gz"
if [ -f "$golden_vcf" ]; then
    golden_count=$(zcat "$golden_vcf" 2>/dev/null | grep -vc '^#' || echo 0)
    delta=$((variant_count - golden_count))
    detail="{\"test\":$variant_count,\"baseline\":$golden_count,\"delta\":$delta}"
    _record_check variant.count PASS \
        "test $variant_count / golden $golden_count (Δ $delta)" "$detail"
else
    _record_check variant.count INFO \
        "$variant_count variants (no golden VCF yet; run --update-golden)" \
        "{\"test\":$variant_count}"
fi

check_bam align_bam  "$WORK/${sample}.align.bam"
check_vcf gvcf       "$WORK/${sample}.variantCaller.g.vcf.gz"
check_vcf vcf        "$WORK/${sample}.genotyper.vcf.gz"

case_end
