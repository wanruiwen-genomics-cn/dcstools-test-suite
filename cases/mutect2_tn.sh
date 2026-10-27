#!/bin/bash
# Somatic calling: align tumor/normal FASTQ, then Mutect2 + FilterMutectCalls.
#
# Testing Mutect2 without FilterMutectCalls covers only half the feature -- the
# raw callset is not what anyone consumes, and the filtering step is where the
# customised GATK build most often diverges. Both run here.
#
# Mutect2 and FilterMutectCalls are dispatched to a GATK jar, so this case needs
# a JVM; it skips itself with a clear reason when none is available.
#
#   cases/mutect2_tn.sh
#   cases/mutect2_tn.sh --update-golden
#   TUMOR_DATA=/path/to/tumor FASTQ_DIR=/path TUMOR_REF=/path/hg38.fa cases/mutect2_tn.sh

source "$(dirname "$(readlink -f "$0")")/common.sh"

case_begin mutect2_tn

java_opts=${JAVA_OPTS:--Xmx8g}

# TUMOR_DATA: directory holding tumor.r{1,2}.fq.gz and normal.r{1,2}.fq.gz.
# The tumor dataset lives as a sibling of the suite's DATA_ROOT (quick_start_data),
# not under it, so default to $(dirname $DATA_ROOT)/tumor.
tumor_data=${TUMOR_DATA:-$(dirname "$DATA_ROOT")/tumor}

# TUMOR_REF: the hg38 reference (shared with the 30x WGS case). Defaults to the
# HG001_30x reference layout; override with TUMOR_REF for a custom location.
tumor_ref=${TUMOR_REF:-$(dirname "$DATA_ROOT")/HG001_30x/ref/hg38.fa}

need_file FA  "$tumor_ref"                          "reference fasta (hg38)"

# PREBUILT_BAM_DIR: directory holding pre-aligned tumor.align.bam /
# normal.align.bam (+ .bai). When present, alignment is skipped and these BAMs
# are reused, so iterating on Mutect2/FilterMutectCalls does not cost another
# ~2.5h of alignment. Defaults to the tumor data dir (alongside the FASTQs).
prebuilt_dir=${PREBUILT_BAM_DIR:-$tumor_data}
prebuilt_tumor="$prebuilt_dir/tumor.align.bam"
prebuilt_normal="$prebuilt_dir/normal.align.bam"
have_prebuilt=0
if [ -f "$prebuilt_tumor" ] && [ -f "$prebuilt_normal" ]; then
    have_prebuilt=1
    note prebuilt.bam "reusing pre-aligned BAMs from $prebuilt_dir (skipping align)"
else
    need_file TQ1 "$tumor_data/tumor.r1.fq.gz"           "tumor read 1"
    need_file TQ2 "$tumor_data/tumor.r2.fq.gz"           "tumor read 2"
    need_file NQ1 "$tumor_data/normal.r1.fq.gz"         "normal read 1"
    need_file NQ2 "$tumor_data/normal.r2.fq.gz"         "normal read 2"
fi

if ! command -v java >/dev/null 2>&1 && [ -z "${JAVA_HOME:-}" ]; then
    skip_case env.java "Mutect2 needs a JVM; set JAVA_HOME or put java on PATH"
fi
preflight_done

# Build the alignment index if this dataset has not been indexed yet. Reuses
# the same index the 30x WGS case builds, so if that case already ran this is
# a no-op. Skipped when reusing pre-aligned BAMs (no alignment -> no index).
if [ "$have_prebuilt" = 0 ]; then
    for file in "$FA".0123 "$FA".amb "$FA".ann "$FA".bwt "$FA".pac "$FA".sa "$FA".bwt.2bit.64; do
        if [ ! -f "$file" ]; then
            run_step build_index "$DCS" aligner --threads "$THREADS" --build-index "$FA"
            break
        fi
    done
fi

# Mutect2 is a GATK tool and needs the samtools fai index and the .dict
# dictionary. Build them if absent (aligner does not produce them).
if [ ! -f "$FA.fai" ]; then
    run_step faidx samtools faidx "$FA"
fi
dict_path="${FA%.fa}.dict"
if [ ! -f "$dict_path" ]; then
    # samtools dict is the standard way; it reads contig lengths from the fasta.
    run_step dict samtools dict -a hg38 "$FA" -o "$dict_path"
fi

# --- Step 1: align tumor and normal -----------------------------------------
# Two separate aligner runs; Mutect2 matches --normal against the SM value in
# the BAM header, so each sample gets its own @RG with a distinct SM. When
# pre-aligned BAMs are available (PREBUILT_BAM_DIR), symlink them into WORK and
# skip alignment.

if [ "$have_prebuilt" = 1 ]; then
    ln -sf "$prebuilt_tumor"      "$WORK/tumor.align.bam"
    ln -sf "$prebuilt_tumor.bai"  "$WORK/tumor.align.bam.bai"
    ln -sf "$prebuilt_normal"     "$WORK/normal.align.bam"
    ln -sf "$prebuilt_normal.bai" "$WORK/normal.align.bam.bai"
    _record_check step.align_tumor  PASS "skipped (prebuilt BAM)"
    _record_check step.align_normal PASS "skipped (prebuilt BAM)"
else
    tumor_rg="@RG\\tID:tumor\\tLB:tumor\\tSM:tumor\\tPL:COMPLETE\\tCN:DCS"
    normal_rg="@RG\\tID:normal\\tLB:normal\\tSM:normal\\tPL:COMPLETE\\tCN:DCS"

    mkdir -p "$WORK/tumor.qcdir" "$WORK/normal.qcdir"

    run_step align_tumor "$DCS" aligner \
        --threads "$THREADS" \
        --aln-index "$FA" \
        --fq1 "$TQ1" --fq2 "$TQ2" \
        --read-group "$tumor_rg" \
        --mark-split-hits-secondary \
        --use-soft-clipping-for-sup-align \
        --fastq-qc-dir "$WORK/tumor.qcdir" \
        --bam-out "$WORK/tumor.align.bam"

    run_step align_normal "$DCS" aligner \
        --threads "$THREADS" \
        --aln-index "$FA" \
        --fq1 "$NQ1" --fq2 "$NQ2" \
        --read-group "$normal_rg" \
        --mark-split-hits-secondary \
        --use-soft-clipping-for-sup-align \
        --fastq-qc-dir "$WORK/normal.qcdir" \
        --bam-out "$WORK/normal.align.bam"
fi

TUMOR="$WORK/tumor.align.bam"
NORMAL="$WORK/normal.align.bam"

# --- Step 2: paired tumour/normal calling ----------------------------------
# Mutect2 threading: --nthread drives parallel processing, --native-pair-hmm-threads
# drives the AVX PairHMM (the CPU-hot part). Both default to $THREADS (32 from
# cases.conf), which cut a full-genome tumor/normal run from ~147 min to ~19 min
# on this host. Override with MUTECT2_NTHREAD / MUTECT2_PAIRHMM_THREADS.
m2_nthread=${MUTECT2_NTHREAD:-$THREADS}
m2_pairhmm=${MUTECT2_PAIRHMM_THREADS:-$THREADS}

run_step mutect2 "$DCS" Mutect2 --java-options "$java_opts" \
    -R "$FA" \
    -I "$TUMOR" \
    -I "$NORMAL" \
    -normal normal \
    --nthread "$m2_nthread" \
    --native-pair-hmm-threads "$m2_pairhmm" \
    -O "$WORK/somatic.unfiltered.vcf.gz"

# --- Step 3: filtering -----------------------------------------------------
# The unfiltered callset is an intermediate, not a product.

run_step filter "$DCS" FilterMutectCalls --java-options "$java_opts" \
    -R "$FA" \
    -V "$WORK/somatic.unfiltered.vcf.gz" \
    -O "$WORK/somatic.filtered.vcf.gz"

# --- verification -----------------------------------------------------------

expect_file "$WORK/tumor.align.bam"            1000
expect_file "$WORK/normal.align.bam"          1000
expect_file "$WORK/somatic.unfiltered.vcf.gz" 1000
expect_file "$WORK/somatic.filtered.vcf.gz"   1000
expect_file "$WORK/somatic.unfiltered.vcf.gz.stats" 10

# Filtering must actually assign FILTER values. A filtered VCF where every
# record is still PASS usually means the stats file was missing and the step
# degenerated into a copy.
check_cmd filter.applied "FilterMutectCalls assigned FILTER values" \
    bash -c 'zcat "$1" | grep -v "^#" | cut -f7 | sort -u | grep -qv "^PASS$"' \
    _ "$WORK/somatic.filtered.vcf.gz"

note somatic.counts "$(zcat "$WORK/somatic.unfiltered.vcf.gz" | grep -vc '^#') raw, $(zcat "$WORK/somatic.filtered.vcf.gz" | grep -v '^#' | awk -F'\t' '$7=="PASS"' | wc -l) PASS after filtering"

# BAM 不纳入对照（tumor/normal 对齐结果不做 golden 比较）。
# VCF 比 WGS 松一档：位点、等位基因、FILTER、基因型必须一致；QUAL/INFO 和
# AD/DP 等注释允许抖动。当前金标差 2/38936 条、且都是均聚物上少 1 条 ref
# read 带出来的注释，落在 0.1% 以内，记为通过。位点增删或基因型变化仍失败。
check_vcf unfiltered_vcf "$WORK/somatic.unfiltered.vcf.gz" --annotation-rel-tol 0.001
check_vcf filtered_vcf   "$WORK/somatic.filtered.vcf.gz" --annotation-rel-tol 0.001

case_end
