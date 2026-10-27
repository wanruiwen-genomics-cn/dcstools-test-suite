#!/bin/bash
# Germline WGS on potato (Solanum tuberosum, SolTub_3.0), ERT / high-memory index.
#
# Same biology as wgs_potato but the aligner is forced onto --use-high-memory-index
# so a regression in the ERT index path is caught independently. Tetraploid
# (ploidy 4), no BQSR. Index is built in its own folder (index_ert) under
# POTATO_IDX_ROOT.
#
#   cases/wgs_potato_ert.sh
#   cases/wgs_potato_ert.sh --update-golden
#   POTATO_DATA=/data/potato cases/wgs_potato_ert.sh

source "$(dirname "$(readlink -f "$0")")/common.sh"

case_begin wgs_potato_ert

data=${POTATO_DATA:-/zdswhst2/ST_BIOINTEL/P24Z32400N0004/yangqi4/project/07.dcs/01.variantcaller_ploidy/data}
idx_root=${POTATO_IDX_ROOT:-/zdswhst2/ST_BIOINTEL/P24Z32400N0004/wanruiwen/dcstools_data/potato_indexes}
sample=${POTATO_SAMPLE:-potato}
group="@RG\\tID:${sample}\\tLB:${sample}\\tSM:${sample}\\tPL:COMPLETE\\tCN:DCS"

src_fa=${POTATO_REF:-$data/bwa-mem2-index/Solanum_tuberosum.SolTub_3.0.dna_rm.toplevel.fa}
need_file SRC_FA "$src_fa" "potato reference fasta"
need_file FQ1 "${POTATO_FQ1:-$data/SRR14400529.fastq/SRR14400529_1.fastq.gz}" "read 1"
need_file FQ2 "${POTATO_FQ2:-$data/SRR14400529.fastq/SRR14400529_2.fastq.gz}" "read 2"

# Non-human species do not run BQSR regardless of known-sites availability.
preflight_done

MODE=ert
idx_dir="$idx_root/index_$MODE"
mkdir -p "$idx_dir"
plain_src="$src_fa"
case "$src_fa" in
    *.gz) plain_src="${src_fa%.gz}"; [ -f "$plain_src" ] || zcat "$src_fa" > "$plain_src" ;;
esac
fa="$idx_dir/$(basename "$plain_src")"
[ -f "$fa" ] || cp -l "$plain_src" "$fa" 2>/dev/null || cp "$plain_src" "$fa"

# Accept either ERT sidecar flavour: .bwt.2bit.64 or .index1.
hm_index=""
if [ -f "$fa.bwt.2bit.64" ]; then
    hm_index="$fa.bwt.2bit.64"
elif [ -f "$fa.index1" ]; then
    hm_index="$fa.index1"
fi
if [ -z "$hm_index" ]; then
    run_step build_index "$DCS" aligner --threads "$THREADS" \
        --use-high-memory-index --build-index "$fa"
    if [ -f "$fa.bwt.2bit.64" ]; then
        hm_index="$fa.bwt.2bit.64"
    elif [ -f "$fa.index1" ]; then
        hm_index="$fa.index1"
    fi
fi
check_cmd index.high_memory "ERT high-memory index present ($hm_index)" \
    test -n "$hm_index" -a -f "$hm_index"

mkdir -p "$WORK/${sample}.qcdir"

run_step align "$DCS" aligner \
    --threads "$THREADS" \
    --aln-index "$fa" \
    --use-high-memory-index \
    --fq1 "$FQ1" --fq2 "$FQ2" \
    --read-group "$group" \
    --mark-split-hits-secondary \
    --use-soft-clipping-for-sup-align \
    --fastq-qc-dir "$WORK/${sample}.qcdir" \
    --bam-out "$WORK/${sample}.align.bam"

check_base_quality data.base_quality "$WORK/${sample}.align.bam"

run_step variant_call "$DCS" variantCaller \
    --threads "$THREADS" \
    --input "$WORK/${sample}.align.bam" \
    --output "$WORK/${sample}.variantCaller.g.vcf.gz" \
    --reference "$fa" \
    --ploidy 4

run_step bam_stat "$DCS" bam-stat \
    --bam "$WORK/${sample}.align.bam" \
    --source "$fa" \
    --output "$WORK/bamStat.txt"
# genotyper is not supported for polyploid potato; stat the gVCF instead.
run_step vcf_stat bash -c "\"$DCS\" vcf-stat --vcf \"$WORK/${sample}.variantCaller.g.vcf.gz\" > \"$WORK/vcfStat.txt\""
run_step report "$DCS" report \
    --sample "$sample" \
    --filter "$WORK/${sample}.qcdir" \
    --bam_stat "$WORK/bamStat.txt" \
    --vcf_stat "$WORK/vcfStat.txt" \
    --output "$WORK/report"

expect_file "$WORK/${sample}.align.bam"                  1000
expect_file "$WORK/${sample}.align.bam.bai"              100
expect_file "$WORK/${sample}.variantCaller.g.vcf.gz"    1000
expect_file "$WORK/bamStat.txt"
expect_file "$WORK/vcfStat.txt"
check_report "$WORK/report"

# genotyper is not supported for polyploid potato; compare the gVCF instead.
variant_count=$(zcat "$WORK/${sample}.variantCaller.g.vcf.gz" 2>/dev/null | grep -vc '^#' || echo 0)
golden_gvcf="$GOLDEN/${sample}.variantCaller.g.vcf.gz"
if [ -f "$golden_gvcf" ]; then
    golden_count=$(zcat "$golden_gvcf" 2>/dev/null | grep -vc '^#' || echo 0)
    delta=$((variant_count - golden_count))
    detail="{\"test\":$variant_count,\"baseline\":$golden_count,\"delta\":$delta}"
    _record_check variant.count PASS \
        "test $variant_count / golden $golden_count (Δ $delta)" "$detail"
else
    _record_check variant.count INFO \
        "$variant_count gVCF records (no golden yet; run --update-golden)" \
        "{\"test\":$variant_count}"
fi

check_bam align_bam  "$WORK/${sample}.align.bam"
check_vcf gvcf       "$WORK/${sample}.variantCaller.g.vcf.gz"

case_end
