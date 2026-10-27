#!/bin/bash
# Large-genome alignment on bread wheat (Triticum aestivum / 小麦).
#
# The ~14.5 Gb hexaploid genome is the originally planned large-plant case
# (wheat_large_genome in the suite outline). Like cycas_large_genome it
# exercises --use-high-memory-index / the ERT .bwt.2bit.64 index, which is
# what a production run on a genome this size uses. Hexaploid genotyping
# itself is covered by potato_polyploid; this case stays on the large-genome
# alignment and diploid call path so a regression in the ERT index is not
# mixed up with a ploidy failure.
#
# BQSR is optional: wheat known-sites files are typically absent. The case
# recalibrates when they are present and skips straight to calling when they
# are not.
#
#   cases/wheat_large_genome.sh
#   cases/wheat_large_genome.sh --update-golden
#   WHEAT_DATA=/data/Triticum_aestivum WHEAT_SAMPLE=wheat cases/wheat_large_genome.sh
#   WHEAT_FQ1=/path/r1.fq.gz WHEAT_FQ2=/path/r2.fq.gz cases/wheat_large_genome.sh

source "$(dirname "$(readlink -f "$0")")/common.sh"

case_begin wheat_large_genome

data=${WHEAT_DATA:-$DATA_ROOT/wheat}
sample=${WHEAT_SAMPLE:-wheat}
group="@RG\\tID:${sample}\\tLB:${sample}\\tSM:${sample}\\tPL:COMPLETE\\tCN:DCS"

# Prefer the RefSeq v2 fasta name used on the cluster; fall back to wheat.fa
# if someone has installed a standardised layout under DATA_ROOT/wheat.
if [ -n "${WHEAT_REF:-}" ]; then
    ref_fa=$WHEAT_REF
elif [ -f "$data/ref/Triticum_aestivum_refseqv2.fa" ]; then
    ref_fa="$data/ref/Triticum_aestivum_refseqv2.fa"
else
    ref_fa="$data/ref/wheat.fa"
fi
need_file FA  "$ref_fa" "wheat reference fasta"

# FASTQs live in a different tree on this cluster than the reference. WHEAT_FQ1
# / WHEAT_FQ2 override the standard fastq/<sample>.r{1,2}.fq.gz layout.
need_file FQ1 "${WHEAT_FQ1:-$data/fastq/${sample}.r1.fq.gz}" "read 1"
need_file FQ2 "${WHEAT_FQ2:-$data/fastq/${sample}.r2.fq.gz}" "read 2"
# Known sites are optional for a non-model organism.
known_sites=()
for vcf in "$data"/known_sites/*.vcf.gz; do
    [ -f "$vcf" ] && known_sites+=("$vcf")
done

preflight_done

# Step 0: build the alignment index if this dataset has not been indexed yet.
# Wheat on the cluster already has the ERT high-memory files (.bwt.2bit.64,
# .0123) and typically no legacy .bwt/.sa. Checking only the ERT set avoids
# a multi-hour rebuild of a 14 Gb index just because .sa is absent.
for file in "$FA".0123 "$FA".amb "$FA".ann "$FA".pac "$FA".bwt.2bit.64; do
    if [ ! -f "$file" ]; then
        run_step build_index "$DCS" aligner --threads "$THREADS" --build-index "$FA"
        break
    fi
done

# Assert the high-memory index file is present before attempting to use it.
# --use-high-memory-index silently falls back or crashes depending on the
# build; failing here with a clear message is better than a mysterious OOM
# three minutes into the alignment.
hm_index="$FA.bwt.2bit.64"
check_cmd index.high_memory "ERT index file present ($hm_index)" \
    test -f "$hm_index"

mkdir -p "$WORK/${sample}.qcdir"

# Step 1: alignment with the high-memory (ERT) index. This is the path under
# test. --use-high-memory-index loads the .bwt.2bit.64 file for faster seeding
# at the cost of more RAM, which is the trade-off a large genome demands.
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

# Step 3: variant calling into a gVCF.
run_step variant_call "$DCS" variantCaller \
    --threads "$THREADS" \
    --input "$WORK/${sample}.align.bam" \
    --output "$WORK/${sample}.variantCaller.g.vcf.gz" \
    --reference "$FA" \
    $recal_arg

# Step 4: single sample genotyping.
# Pinned to GENOTYPER_THREADS (default 1, see env/cases.conf): genotyper's QD
# jitter RNG is a shared singleton; multi-threaded runs are nondeterministic.
run_step genotype "$DCS" genotyper \
    --threads "$GENOTYPER_THREADS" \
    --input "$WORK/${sample}.variantCaller.g.vcf.gz" \
    --output "$WORK/${sample}.genotyper.vcf.gz" \
    -s 30.0

# Step 5: statistics.
run_step bam_stat "$DCS" bam-stat \
    --bam "$WORK/${sample}.align.bam" \
    --source "$FA" \
    --output "$WORK/bamStat.txt"
run_step vcf_stat bash -c "\"$DCS\" vcf-stat --vcf \"$WORK/${sample}.genotyper.vcf.gz\" > \"$WORK/vcfStat.txt\""

# --- verification -----------------------------------------------------------

expect_file "$WORK/${sample}.align.bam"                  1000
expect_file "$WORK/${sample}.align.bam.bai"              100
expect_file "$WORK/${sample}.variantCaller.g.vcf.gz"    1000
expect_file "$WORK/${sample}.genotyper.vcf.gz"           1000
expect_file "$WORK/${sample}.genotyper.vcf.gz.tbi"        50
expect_file "$WORK/bamStat.txt"
expect_file "$WORK/vcfStat.txt"

# A large-genome run that produced zero variants means the reference contig
# names did not match the read names, or the caller silently dropped everything.
# The golden comparison would pass two empty files, so assert nonempty here.
check_cmd vcf.nonempty "genotyper emitted at least one variant" \
    bash -c 'test "$(zcat "$1" | grep -vc "^#")" -gt 0' _ "$WORK/${sample}.genotyper.vcf.gz"

check_bam align_bam  "$WORK/${sample}.align.bam"
check_vcf gvcf       "$WORK/${sample}.variantCaller.g.vcf.gz"
check_vcf vcf        "$WORK/${sample}.genotyper.vcf.gz"

case_end
