#!/bin/bash
# Large-genome alignment on Cycas (苏铁), a gymnosperm with a ~10-20 Gb genome.
#
# The small human case and the 30x case both use the default minibwa index.
# This one exercises the other code path: --use-high-memory-index, the ERT
# index, which is what a production run on a large plant genome uses because
# the default index is too slow or too large to fit. A regression that breaks
# the ERT path would pass every other case in the suite.
#
# BQSR is optional: Cycas is a non-model organism and may not have known-sites
# files. The case recalibrates when they are present and skips straight to
# calling when they are not.
#
#   cases/cycas_large_genome.sh
#   cases/cycas_large_genome.sh --update-golden
#   CYCAS_DATA=/data/cycas CYCAS_SAMPLE=leaf01 cases/cycas_large_genome.sh

source "$(dirname "$(readlink -f "$0")")/common.sh"

case_begin cycas_large_genome

data=${CYCAS_DATA:-$DATA_ROOT/cycas}
sample=${CYCAS_SAMPLE:-cycas}
group="@RG\\tID:${sample}\\tLB:${sample}\\tSM:${sample}\\tPL:COMPLETE\\tCN:DCS"

need_file FA  "$data/ref/cycas.fa"              "cycas reference fasta"
need_file FQ1 "$data/fastq/${sample}.r1.fq.gz"  "read 1"
need_file FQ2 "$data/fastq/${sample}.r2.fq.gz"  "read 2"
# Known sites are optional for a non-model organism.
known_sites=()
for vcf in "$data"/known_sites/*.vcf.gz; do
    [ -f "$vcf" ] && known_sites+=("$vcf")
done

preflight_done

# Step 0: build the alignment index if this dataset has not been indexed yet.
# The ERT high-memory index is the .bwt.2bit.64 file; --build-index produces
# both the default and the high-memory index, so one build covers both paths.
for file in "$FA".0123 "$FA".amb "$FA".ann "$FA".bwt "$FA".index1 "$FA".index2 "$FA".pac "$FA".sa; do
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
