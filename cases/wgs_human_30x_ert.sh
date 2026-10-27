#!/bin/bash
# Germline WGS 30x on human data, alignment with the ERT / high-memory index.
#
# Same biology as wgs_human_30x (HG001 30x FASTQ + calling), but the aligner
# step is forced onto --use-high-memory-index so a regression in the ERT
# index path is caught independently of the default / minibwa path.
#
# Reference defaults to the cluster ERT bundle (hg38.fa + .index1/.index2 or
# .bwt.2bit.64). FASTQ still comes from the HG001_30x dataset.
#
#   cases/wgs_human_30x_ert.sh
#   cases/wgs_human_30x_ert.sh --update-golden
#   WGS30X_DATA=/data/HG001_30x WGS30X_ERT_REF=/ref/hg38.fa cases/wgs_human_30x_ert.sh

source "$(dirname "$(readlink -f "$0")")/common.sh"

case_begin wgs_human_30x_ert

data=${WGS30X_DATA:-$DATA_ROOT/wgs30x}
sample=${WGS30X_SAMPLE:-sample}
group="@RG\\tID:${sample}\\tLB:${sample}\\tSM:${sample}\\tPL:COMPLETE\\tCN:DCS"

# ERT / high-memory reference. Override with WGS30X_ERT_REF if the bundle moves.
need_file FA "${WGS30X_ERT_REF:-/zdswhst2/ST_BIOINTEL/P24Z32400N0004/wanruiwen/ref_ert/hg38.fa}" \
    "ERT high-memory reference fasta"
need_file FQ1 "$data/fastq/${sample}.r1.fq.gz" "read 1"
need_file FQ2 "$data/fastq/${sample}.r2.fq.gz" "read 2"

# Prefer known sites next to the ERT ref; fall back to the 30x dataset layout.
known_sites=()
for dir in "$(dirname "$FA")/known_sites" "$data/known_sites"; do
    [ -d "$dir" ] || continue
    for vcf in "$dir"/*.vcf.gz; do
        [ -f "$vcf" ] && known_sites+=("$vcf")
    done
    [ ${#known_sites[@]} -gt 0 ] && break
done

preflight_done

# Build the ERT index only when neither flavour of high-memory sidecar is present.
# Human ERT bundles on this cluster ship .index1/.index2; wheat-style bundles
# ship .bwt.2bit.64. Accept either so we do not rebuild a multi-hour index.
hm_index=""
if [ -f "$FA.bwt.2bit.64" ]; then
    hm_index="$FA.bwt.2bit.64"
elif [ -f "$FA.index1" ]; then
    hm_index="$FA.index1"
fi
if [ -z "$hm_index" ]; then
    run_step build_index "$DCS" aligner --threads "$THREADS" --build-index "$FA"
    if [ -f "$FA.bwt.2bit.64" ]; then
        hm_index="$FA.bwt.2bit.64"
    elif [ -f "$FA.index1" ]; then
        hm_index="$FA.index1"
    fi
fi
check_cmd index.high_memory "ERT high-memory index present ($hm_index)" \
    test -n "$hm_index" -a -f "$hm_index"

mkdir -p "$WORK/${sample}.qcdir"

# Alignment under test: force the high-memory / ERT index path.
run_step align "$DCS" aligner \
    --threads "$THREADS" \
    --aln-index "$FA" \
    --use-high-memory-index \
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

# Pinned to GENOTYPER_THREADS (default 1, see env/cases.conf): genotyper's QD
# jitter RNG is a shared singleton; multi-threaded runs are nondeterministic.
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
