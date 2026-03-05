#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") -s SAMPLE -b INPUT.bam -o OUT_DIR [-t THREADS] [-C chr1,chr2,...] [--no-dedup] [--loglevel LEVEL]

Runs: maxatac prepare
  -s, --sample      Sample name (used for -n)
  -b, --bam         Input BAM (chr-prefixed contigs)
  -o, --outdir      Output directory (created if missing)
  -t, --threads     Threads for maxATAC (default: 8)
  -C, --chroms      Chromosomes (comma OR space separated). Default: chr1..chr19
      --no-dedup    Do NOT pass --skip_deduplication (dedup will run)
      --loglevel    maxATAC loglevel (default: debug)
  -h, --help        Show help
EOF
}

SAMPLE=""
BAM=""
OUT=""
THREADS=8
CHROMS_DEFAULT="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19"
CHROMS="$CHROMS_DEFAULT"
SKIP_DEDUP=1
LOGLEVEL="debug"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--sample)   SAMPLE="$2"; shift 2;;
    -b|--bam)      BAM="$2"; shift 2;;
    -o|--outdir)   OUT="$2"; shift 2;;
    -t|--threads)  THREADS="${2}"; shift 2;;
    -C|--chroms)   CHROMS="$2"; shift 2;;
    --no-dedup)    SKIP_DEDUP=0; shift 1;;
    --loglevel)    LOGLEVEL="$2"; shift 2;;
    -h|--help)     usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

[[ -z "${SAMPLE}" ]] && { echo "ERROR: -s/--sample required"; exit 1; }
[[ -z "${BAM}"    ]] && { echo "ERROR: -b/--bam required"; exit 1; }
[[ -z "${OUT}"    ]] && { echo "ERROR: -o/--outdir required"; exit 1; }
[[ -f "${BAM}"    ]] || { echo "ERROR: BAM not found: ${BAM}"; exit 1; }

# Make output dir
mkdir -p "${OUT}"

# Chromosomes: allow comma or space separated input; split to array for maxATAC
CHROMS="${CHROMS//,/ }"
read -r -a CHR_ARR <<< "${CHROMS}"

# Optional dedup flag
DEDUP_FLAG="--skip_deduplication"
[[ "${SKIP_DEDUP}" -eq 0 ]] && DEDUP_FLAG=""

echo "[maxATAC-prepare] sample=${SAMPLE}"
echo "[maxATAC-prepare] bam=${BAM}"
echo "[maxATAC-prepare] out=${OUT}"
echo "[maxATAC-prepare] threads=${THREADS}"
echo "[maxATAC-prepare] chroms=${CHROMS}"
echo "[maxATAC-prepare] dedup_flag='${DEDUP_FLAG}'"
echo "[maxATAC-prepare] loglevel=${LOGLEVEL}"

maxatac prepare \
  -i "${BAM}" \
  -o "${OUT}" \
  -n "${SAMPLE}" \
  ${DEDUP_FLAG} \
  -t "${THREADS}" \
  --loglevel "${LOGLEVEL}" \
  --chromosomes "${CHR_ARR[@]}"

# Sentinel for Nextflow (optional)
touch "${OUT}/.done_prepare"

