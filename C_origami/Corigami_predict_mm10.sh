#!/usr/bin/env bash
# corigami_predict_mm10.sh
# Run C.Origami prediction over chosen chromosomes / windows on mm10.

set -euo pipefail

### --- USER CONFIG (edit once; reuse everywhere) ---
SAMPLE=${SAMPLE:-sample}
OUT=${OUT:-"$PWD/Predict_TADs_${SAMPLE}"}
CS=${CS:-/home/opt/maxatac/data/mm10/mm10.chrom.sizes}
SEQ_DIR=${SEQ_DIR:-/projects/Reference/mm10_perChromosome/dna_sequence}
MODEL=${MODEL:-/home/Reference/corigami_data/model_weights/corigami_base.ckpt}
ATAC=${ATAC:-"$PWD/out_maxatac_prepare/${SAMPLE}_IS_slop20_RP20M.bw"}
CTCF=${CTCF:-"$PWD/pred_${SAMPLE}_CTCF/${SAMPLE}_CTCF.bw"}  # maxATAC predict output

# Windowing (2,097,152 bp default). STEP=WIN for butt-joined tiles; set STEP=1048576 for 50% overlap.
WIN=${WIN:-2097152}
STEP=${STEP:-2097152}

mkdir -p "$OUT"

### --- CLI HELP ---
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<USAGE
Usage:
  $(basename "$0")                # default: chr1 only
  $(basename "$0") 3              # chr1..chr3
  $(basename "$0") chr1 chrX      # explicit list
Env overrides (optional):
  SAMPLE OUT CS SEQ_DIR MODEL ATAC CTCF WIN STEP
USAGE
  exit 0
fi

### --- choose chromosomes from CLI args ---
DEFAULT_CHRS="chr1"
if (( $# == 0 )); then
  CHRS="$DEFAULT_CHRS"
elif [[ "$1" =~ ^[0-9]+$ ]]; then
  N="$1"
  CHRS=""
  for ((i=1; i<=N; i++)); do CHRS+=" chr${i}"; done
  CHRS="${CHRS# }"
else
  CHRS="$*"
fi
echo "[info] Chromosomes: $CHRS"

### --- sanity checks ---
for f in "$MODEL" "$ATAC" "$CTCF" "$CS"; do
  [[ -s "$f" ]] || { echo "[error] Missing file: $f" >&2; exit 1; }
done

command -v corigami-predict >/dev/null || { echo "[error] corigami-predict not found in PATH"; exit 1; }

### --- run ---
for CHR in $CHRS; do
  LEN=$(awk -v c="$CHR" '$1==c{print $2}' "$CS")
  if [[ -z "${LEN:-}" ]]; then
    echo "[warn] Skip $CHR (not in $CS)"; continue
  fi
  if [[ ! -s "$SEQ_DIR/${CHR}.fa.gz" ]]; then
    echo "[warn] Skip $CHR (missing ${SEQ_DIR}/${CHR}.fa.gz)"; continue
  fi

  echo "[info] $CHR length: $LEN; WIN=$WIN STEP=$STEP"

  # main tiles
  for ((START=0; START <= LEN-WIN; START+=STEP)); do
    echo "[run] $CHR:$START"
    corigami-predict \
      --out      "$OUT" \
      --celltype "$SAMPLE" \
      --chr      "$CHR" \
      --start    "$START" \
      --model    "$MODEL" \
      --seq      "$SEQ_DIR" \
      --ctcf     "$CTCF" \
      --atac     "$ATAC"
  done

  # tail tile to cover chromosome end if STEP doesn’t divide LEN (optional; keeps last 2Mb anchored to end)
  TAIL=$(( LEN - WIN ))
  if (( TAIL > 0 )) && (( (LEN-1) % STEP != 0 )); then
    echo "[run] tail $CHR:$TAIL"
    corigami-predict \
      --out      "$OUT" \
      --celltype "$SAMPLE" \
      --chr      "$CHR" \
      --start    "$TAIL" \
      --model    "$MODEL" \
      --seq      "$SEQ_DIR" \
      --ctcf     "$CTCF" \
      --atac     "$ATAC"
  fi
done

echo "[done] Output in: $OUT"
