#!/bin/bash
# Germline WGS 30x on human data, alignment with the minibwa index.
#
# Same biology as wgs_human_30x (HG001 30x FASTQ + calling), but the aligner
# step is forced onto --use-minibwa so a regression in the minibwa (.l2b/.mbw)
# path is caught independently of the default / ERT path.
#
# Reference defaults to the cluster minibwa bundle (hg38.fa + .l2b/.mbw).
# FASTQ still comes from the HG001_30x dataset.
#
#   cases/wgs_human_30x_minibwa.sh
#   cases/wgs_human_30x_minibwa.sh --update-golden
#   WGS30X_DATA=/data/HG001_30x WGS30X_MINIBWA_REF=/ref/hg38.fa cases/wgs_human_30x_minibwa.sh

source "$(dirname "$(readlink -f "$0")")/common.sh"

case_begin wgs_human_30x_minibwa

data=${WGS30X_DATA:-$DATA_ROOT/wgs30x}
sample=${WGS30X_SAMPLE:-sample}
group="@RG\\tID:${sample}\\tLB:${sample}\\tSM:${sample}\\tPL:COMPLETE\\tCN:DCS"

# Minibwa reference. Override with WGS30X_MINIBWA_REF if the bundle moves.
need_file FA "${WGS30X_MINIBWA_REF:-/zdswhst2/ST_BIOINTEL/P24Z32400N0004/wanruiwen/ref_minibwa/hg38.fa}" \
    "minibwa reference fasta"
need_file FQ1 "$data/fastq/${sample}.r1.fq.gz" "read 1"
need_file FQ2 "$data/fastq/${sample}.r2.fq.gz" "read 2"

# Prefer known sites from the 30x dataset; fall back to a sibling of the ref.
known_sites=()
for dir in "$data/known_sites" "$(dirname "$FA")/known_sites"; do
    [ -d "$dir" ] || continue
    for vcf in "$dir"/*.vcf.gz; do
        [ -f "$vcf" ] && known_sites+=("$vcf")
    done
    [ ${#known_sites[@]} -gt 0 ] && break
done

preflight_done

# Build the minibwa index only when .l2b / .mbw are missing. Pass
# --use-minibwa on build so the right sidecars are produced.
if [ ! -f "$FA.l2b" ] || [ ! -f "$FA.mbw" ]; then
    run_step build_index "$DCS" aligner --threads "$THREADS" \
        --use-minibwa --build-index "$FA"
fi
check_cmd index.minibwa "minibwa index present ($FA.l2b / $FA.mbw)" \
    test -f "$FA.l2b" -a -f "$FA.mbw"

mkdir -p "$WORK/${sample}.qcdir"

# Alignment under test: force the minibwa index path.
run_step align "$DCS" aligner \
    --threads "$THREADS" \
    --aln-index "$FA" \
    --use-minibwa \
    --fq1 "$FQ1" --fq2 "$FQ2" \
    --read-group "$group" \
    --mark-split-hits-secondary \
    --use-soft-clipping-for-sup-align \
    --fastq-qc-dir "$WORK/${sample}.qcdir" \
    --bam-out "$WORK/${sample}.align.bam"

check_base_quality data.base_quality "$WORK/${sample}.align.bam"

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

run_step variant_call "$DCS" variantCaller \
    --threads "$THREADS" \
    --input "$WORK/${sample}.align.bam" \
    --output "$WORK/${sample}.variantCaller.g.vcf.gz" \
    --reference "$FA" \
    $recal_arg

# Single-sample genotyping. Pinned to GENOTYPER_THREADS (default 1) because
# genotyper's QD jitter RNG is a shared singleton; multi-threaded runs are
# nondeterministic and would not match the golden VCF.
run_step genotype "$DCS" genotyper \
    --threads "$GENOTYPER_THREADS" \
    --input "$WORK/${sample}.variantCaller.g.vcf.gz" \
    --output "$WORK/${sample}.genotyper.vcf.gz" \
    -s 30.0

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

expect_file "$WORK/${sample}.align.bam"                  1000
expect_file "$WORK/${sample}.align.bam.bai"              100
expect_file "$WORK/${sample}.variantCaller.g.vcf.gz"    1000
expect_file "$WORK/${sample}.genotyper.vcf.gz"           1000
expect_file "$WORK/${sample}.genotyper.vcf.gz.tbi"        50
expect_file "$WORK/bamStat.txt"
expect_file "$WORK/vcfStat.txt"
check_report "$WORK/report"

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
