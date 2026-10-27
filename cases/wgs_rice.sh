#!/bin/bash
# Germline WGS on rice (Oryza sativa, IRGSP-1.0), default aligner index.
#
# Mirrors wgs_human_30x but on rice: real Illumina WGS, diploid calling, no
# BQSR (non-human species skip recalibration). The index is built in its own folder
# (index_default) under the dataset so the three index modes (default / ERT /
# minibwa) never share or clobber sidecars.
#
#   cases/wgs_rice.sh
#   cases/wgs_rice.sh --update-golden
#   RICE_DATA=/data/Oryza_sativa cases/wgs_rice.sh

source "$(dirname "$(readlink -f "$0")")/common.sh"

case_begin wgs_rice

data=${RICE_DATA:-/zdswhst2/ST_BIOINTEL/P24Z32400N0004/wanruiwen/dcstools_data/Oryza_sativa}
sample=${RICE_SAMPLE:-rice}
group="@RG\\tID:${sample}\\tLB:${sample}\\tSM:${sample}\\tPL:COMPLETE\\tCN:DCS"

src_fa=${RICE_REF:-$data/ref/Oryza_sativa.IRGSP-1.0.dna_rm.toplevel.fa}
need_file SRC_FA "$src_fa" "rice reference fasta"
need_file FQ1 "${RICE_FQ1:-$data/SRR32814349_1.fastq}" "read 1"
need_file FQ2 "${RICE_FQ2:-$data/SRR32814349_2.fastq}" "read 2"

# Non-human species do not run BQSR regardless of known-sites availability.
preflight_done

# Index dir: build the default (classic BWA) index here, isolated from the
# ERT / minibwa folders. Decompress a gzipped ref once into ref/ first, then
# hardlink the plain fa into this folder so sidecars land beside it.
MODE=default
idx_dir="$data/index_$MODE"
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
    --ploidy 2

# Pinned to GENOTYPER_THREADS (default 1, see env/cases.conf): genotyper's QD
# jitter RNG is a shared singleton; multi-threaded runs are nondeterministic.
run_step genotype "$DCS" genotyper \
    --threads "$GENOTYPER_THREADS" \
    --input "$WORK/${sample}.variantCaller.g.vcf.gz" \
    --output "$WORK/${sample}.genotyper.vcf.gz" \
    -s 30.0

run_step bam_stat "$DCS" bam-stat \
    --bam "$WORK/${sample}.align.bam" \
    --source "$fa" \
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
