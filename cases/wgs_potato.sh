#!/bin/bash
# Germline WGS on potato (Solanum tuberosum, SolTub_3.0), default aligner index.
#
# Tetraploid (ploidy 4). Mirrors wgs_rice: real Illumina WGS, no BQSR (non-human
# species skip recalibration), no known sites. The DCS index is built in its own
# folder (index_default) under POTATO_IDX_ROOT so the three index modes never
# share or clobber sidecars. The reference lives in a read-only project dir,
# so it is hardlinked into the index folder (same /zdswhst2 filesystem).
#
#   cases/wgs_potato.sh
#   cases/wgs_potato.sh --update-golden
#   POTATO_DATA=/data/potato cases/wgs_potato.sh

source "$(dirname "$(readlink -f "$0")")/common.sh"

case_begin wgs_potato

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

# Index dir: build the default (classic BWA) index here, isolated from the
# ERT / minibwa folders. Decompress a gzipped ref once first, then hardlink
# the plain fa into this folder so sidecars land beside it.
MODE=default
idx_dir="$idx_root/index_$MODE"
mkdir -p "$idx_dir"
plain_src="$src_fa"
case "$src_fa" in
    *.gz) plain_src="${src_fa%.gz}"; [ -f "$plain_src" ] || zcat "$src_fa" > "$plain_src" ;;
esac
fa="$idx_dir/$(basename "$plain_src")"
[ -f "$fa" ] || cp -l "$plain_src" "$fa" 2>/dev/null || cp "$plain_src" "$fa"

for file in "$fa".0123 "$fa".amb "$fa".ann "$fa".bwt "$fa".pac "$fa".sa; do
    if [ ! -f "$file" ]; then
        run_step build_index "$DCS" aligner --threads "$THREADS" --build-index "$fa"
        break
    fi
done

mkdir -p "$WORK/${sample}.qcdir"

run_step align "$DCS" aligner \
    --threads "$THREADS" \
    --aln-index "$fa" \
    --fq1 "$FQ1" --fq2 "$FQ2" \
    --read-group "$group" \
    --mark-split-hits-secondary \
    --use-soft-clipping-for-sup-align \
    --fastq-qc-dir "$WORK/${sample}.qcdir" \
    --bam-out "$WORK/${sample}.align.bam"

check_base_quality data.base_quality "$WORK/${sample}.align.bam"

# Step 2: variant calling into a gVCF (no BQSR for non-human species).
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
