#!/bin/bash
# Germline WES: same tools as the WGS case, but driven by a target BED.
#
# The point of keeping this separate from wgs_human_small is the BED path. A
# target file turns on -L/--interval in variantCaller and changes bam-stat into
# its on-target mode, and contig-name or interval bugs only show up here
# (commit 6c46f1d "修复bed区间contig id bug" was exactly this class of defect).
#
# This case runs on the HG001 30x dataset (real Illumina qualities, so the
# caller actually finds variants) restricted to a simulated exome panel in
# $WES_DATA/target.bed. The BED is a handful of ~400 kb regions spread across
# chr1-22 + chrX (~8 Mb total) -- enough at 30x to produce a meaningful callset
# while staying much cheaper than a whole-genome run. Override the data root
# with WES_DATA and the panel with WES_BED.
#
#   cases/wes_human.sh [--update-golden]

source "$(dirname "$(readlink -f "$0")")/common.sh"

case_begin wes_human

data=${WES_DATA:-$DATA_ROOT/HG001_30x}
bed=${WES_BED:-$data/target.bed}
sample=${WES_SAMPLE:-sample}
group="@RG\\tID:${sample}\\tLB:${sample}\\tSM:${sample}\\tPL:COMPLETE\\tCN:DCS"
# Padding around each target; exercised so the flag does not silently rot.
padding=100

need_file FA  "$data/ref/hg38.fa"            "reference fasta"
need_file FQ1 "$data/fastq/${sample}.r1.fq.gz" "read 1"
need_file FQ2 "$data/fastq/${sample}.r2.fq.gz" "read 2"
need_file BED "$bed"                          "target regions (simulated exome panel)"
# Known sites are optional -- BQSR runs only if at least one is present.
known_sites=()
for vcf in "$data"/known_sites/*.vcf.gz; do
    [ -f "$vcf" ] && known_sites+=("$vcf")
done
preflight_done

# Build the alignment index once per dataset; no-op if already present.
for file in "$FA".0123 "$FA".amb "$FA".ann "$FA".bwt "$FA".pac "$FA".sa "$FA".bwt.2bit.64; do
    if [ ! -f "$file" ]; then
        run_step build_index "$DCS" aligner --threads "$THREADS" --build-index "$FA"
        break
    fi
done

mkdir -p "$WORK/${sample}.qcdir"

# The BED contig names must match the reference, otherwise variantCaller quietly
# calls nothing. Assert it before spending the alignment time.
check_cmd data.bed_contigs "BED contig names match the reference" \
    bash -c 'cut -f1 "$1" | sort -u | while read -r c; do grep -qx "$c" <(cut -f1 "$2") || exit 1; done' \
    _ "$BED" "$FA.fai"

run_step align "$DCS" aligner \
    --threads "$THREADS" \
    --aln-index "$FA" \
    --fq1 "$FQ1" --fq2 "$FQ2" \
    --read-group "$group" \
    --mark-split-hits-secondary \
    --use-soft-clipping-for-sup-align \
    --fastq-qc-dir "$WORK/${sample}.qcdir" \
    --bam-out "$WORK/${sample}.align.bam"

# Real 30x Illumina data passes the -Q floor; the check is a guard, not a gate.
check_base_quality data.base_quality "$WORK/${sample}.align.bam"

# BQSR (only if known sites were provided).
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

# Targeted calling: -L restricts to the BED, -P pads each interval.
run_step variant_call "$DCS" variantCaller \
    --threads "$THREADS" \
    --input "$WORK/${sample}.align.bam" \
    --output "$WORK/${sample}.variantCaller.g.vcf.gz" \
    --reference "$FA" \
    -L "$BED" \
    -P "$padding" \
    $recal_arg

# Pinned to GENOTYPER_THREADS (default 1, see env/cases.conf): genotyper's QD
# jitter RNG is a shared singleton; multi-threaded runs are nondeterministic.
run_step genotype "$DCS" genotyper \
    --threads "$GENOTYPER_THREADS" \
    --input "$WORK/${sample}.variantCaller.g.vcf.gz" \
    --output "$WORK/${sample}.genotyper.vcf.gz" \
    -s 30.0

# On-target statistics: bam-stat takes the BED as --source in WES mode.
run_step bam_stat "$DCS" bam-stat \
    --bam "$WORK/${sample}.align.bam" \
    --source "$BED" \
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
expect_file "$WORK/${sample}.variantCaller.g.vcf.gz" 1000
expect_file "$WORK/${sample}.genotyper.vcf.gz"       1000
expect_file "$WORK/${sample}.genotyper.vcf.gz.tbi"   50
expect_file "$WORK/bamStat.txt"
expect_file "$WORK/vcfStat.txt"
check_report "$WORK/report"

# A targeted run that produced zero variants means the BED never matched;
# the golden comparison would happily pass two empty files, so assert it here.
check_cmd vcf.nonempty "genotyper emitted at least one variant" \
    bash -c 'test "$(zcat "$1" | grep -vc "^#")" -gt 0' _ "$WORK/${sample}.genotyper.vcf.gz"

check_bam align_bam "$WORK/${sample}.align.bam"
check_vcf gvcf      "$WORK/${sample}.variantCaller.g.vcf.gz"
check_vcf vcf       "$WORK/${sample}.genotyper.vcf.gz"

case_end
