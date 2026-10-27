#!/bin/bash
# SeqArc compression: index -> compress -> decompress, verified as lossless.
#
# The only acceptance criterion that matters for a compressor is that the data
# survives the round trip. "The script ran without an error" says nothing: a
# compressor that drops or corrupts reads exits 0 just as happily. So this case
# compares the decompressed FASTQ against the original byte for byte and treats
# any difference as a release blocker.
#
#   cases/seqarc_roundtrip.sh

source "$(dirname "$(readlink -f "$0")")/common.sh"

case_begin seqarc_roundtrip

need_file FA  "$DATA_ROOT/ref/chr20.fa"        "reference fasta"
need_file FQ1 "$DATA_ROOT/fastq/simu.r1.fq.gz" "read 1"
preflight_done

ARCREF="$WORK/arcref"
mkdir -p "$ARCREF"

# Step 1: build the SeqArc reference index used by the referenced mode.
run_step index "$DCS" SeqArc index -r "$FA" -o "$ARCREF"

# `index` writes only the sidecar files (<name>.fa.mini_*, <name>.fa.similar_*),
# but `compress -r <dir>/<name>.fa` also wants the FASTA itself at that path.
# quick_start_seqarc.sh relies on it already being there; link it so the case
# works against a clean data directory.
[ -e "$ARCREF/$(basename "$FA")" ] || ln -sf "$(readlink -f "$FA")" "$ARCREF/$(basename "$FA")"
note seqarc.index_size "$(du -sh "$ARCREF" | cut -f1) of index for $(basename "$FA")"

# Step 2: compress into a dedicated directory. SeqArc derives the output name
# from the input rather than using -o verbatim (-o simu.r1 yields
# simu.r1.fq.arc), so the products are located by search afterwards instead of
# being assumed -- quick_start_seqarc.sh assumes "<prefix>.arc" and would break.
ARCDIR="$WORK/arc"
mkdir -p "$ARCDIR"
run_step compress "$DCS" SeqArc compress \
    -t "$THREADS" \
    -r "$ARCREF/$(basename "$FA")" \
    -i "$FQ1" \
    -o "$ARCDIR/$(basename "${FQ1%.gz}")"

ARC=$(find "$ARCDIR" -name "*.arc" | sort | head -1)
if [ -z "$ARC" ]; then
    _record_check product.arc FAIL "compress produced no .arc under $ARCDIR"
    case_end
fi
_record_check product.arc PASS "$(basename "$ARC"), $(stat -c%s "$ARC") bytes"

# Step 3: decompress with the same reference, again into its own directory.
OUTDIR="$WORK/restored"
mkdir -p "$OUTDIR"
run_step decompress "$DCS" SeqArc decompress \
    -t "$THREADS" \
    -r "$ARCREF/$(basename "$FA")" \
    -i "$ARC" \
    -o "$OUTDIR/$(basename "${FQ1%.gz}")"

RESTORED=$(find "$OUTDIR" -name "*.fq" -o -name "*.fastq" | sort | head -1)
if [ -z "$RESTORED" ]; then
    _record_check product.restored FAIL "decompress produced no FASTQ under $OUTDIR"
    case_end
fi
_record_check product.restored PASS "$(basename "$RESTORED"), $(stat -c%s "$RESTORED") bytes"

# --- verification -----------------------------------------------------------

# The acceptance criterion: the decompressed reads must equal the input exactly.
#
# Both sides are normalised to a plain stream first. The input is gzipped, and
# `decompress` also emits BGZF even though -o asked for a name ending in .fq --
# comparing the files as-is would report a false loss. Streaming also avoids
# materialising a second copy of a WGS FASTQ on disk.
_plain() {
    if [ "$(head -c2 "$1" | od -An -tx1 | tr -d ' \n')" = "1f8b" ]; then
        gzip -cd "$1"
    else
        cat "$1"
    fi
}

if [ "$(head -c2 "$RESTORED" | od -An -tx1 | tr -d ' \n')" = "1f8b" ]; then
    note seqarc.output_format "decompress emitted gzip/BGZF despite the .fq output name"
fi

if cmp -s <(_plain "$FQ1") <(_plain "$RESTORED"); then
    _record_check seqarc.lossless PASS "decompressed FASTQ is byte-identical to the input"
else
    _record_check seqarc.lossless FAIL "round trip is not lossless"
    _log "    first difference:"
    cmp <(_plain "$FQ1") <(_plain "$RESTORED") 2>&1 | sed 's/^/        /' | head -3 | tee -a "$LOG"
    _log "    reads: original=$(( $(_plain "$FQ1" | wc -l) / 4 )) decompressed=$(( $(_plain "$RESTORED" | wc -l) / 4 ))"
fi

# Report the ratio so a release that silently stops compressing is visible.
raw=$(_plain "$FQ1" | wc -c)
arc=$(stat -c%s "$ARC")
gz=$(stat -c%s "$FQ1")
note seqarc.ratio "$(awk -v r="$raw" -v a="$arc" -v g="$gz" \
    'BEGIN{printf "raw %.1f MB -> arc %.1f MB (%.2fx, gzip was %.2fx)", r/1048576, a/1048576, r/a, r/g}')"

case_end
