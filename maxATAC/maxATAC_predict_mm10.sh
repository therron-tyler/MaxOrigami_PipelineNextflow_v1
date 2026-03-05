#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") \
  -s SAMPLE \
  -c CHROM \
  -i PREP_DIR \
  -o OUT_DIR \
  [-t THREADS] \
  [-d DATA_DIR] \
  [--signal PATH.bw] \
  [--tf TFNAME ...] \
  [--sequence PATH.2bit] \
  [--chrom-sizes PATH.sizes] \
  [--blacklist PATH.bw] \
  [--loglevel LEVEL]

Runs: maxatac predict for a single chromosome (per-chrom fanout by Nextflow).

Required:
  -s, --sample        Sample name
  -c, --chrom         One chromosome (e.g., chr7)
  -i, --prepdir       Directory from maxATAC prepare (contains signal .bw)
  -o, --outdir        Output directory to write predictions

Optional:
  -t, --threads       Threads (default: 4)
  -d, --data          Dir with default resources (sequence/chrom.sizes/blacklist)
      --signal        Explicit path to signal .bw (overrides auto-detect)
      --tf            TF name (repeatable). Default: CTCF if none given
      --sequence      Genome 2bit (overrides -d default)
      --chrom-sizes   Chrom sizes (overrides -d default)
      --blacklist     Blacklist bigWig (overrides -d default)
      --loglevel      debug|info|warning|error (default: debug)

Notes:
- If --signal is not provided, the script tries PREP_DIR/\$SAMPLE*.bw then PREP_DIR/*.bw.
- If -d is provided and no explicit resources are given, it uses:
    \$DATA_DIR/hg38.2bit
    \$DATA_DIR/hg38.chrom.sizes
    \$DATA_DIR/hg38_maxatac_blacklist.bw
  (Rename to mm10 files as appropriate for your install.)
EOF
}

SAMPLE=""
CHROM=""
PREP=""
OUT=""
THREADS=4
DATA_DIR=""
SIGNAL=""
LOGLEVEL="debug"
TFs=()
SEQ=""
CHROMSIZES=""
BLACKLIST=""

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--sample)       SAMPLE="$2"; shift 2;;
    -c|--chrom)        CHROM="$2"; shift 2;;
    -i|--prepdir)      PREP="$2"; shift 2;;
    -o|--outdir)       OUT="$2"; shift 2;;
    -t|--threads)      THREADS="$2"; shift 2;;
    -d|--data)         DATA_DIR="$2"; shift 2;;
    --signal)          SIGNAL="$2"; shift 2;;
    --tf)              TFs+=("$2"); shift 2;;
    --sequence)        SEQ="$2"; shift 2;;
    --chrom-sizes)     CHROMSIZES="$2"; shift 2;;
    --blacklist)       BLACKLIST="$2"; shift 2;;
    --loglevel)        LOGLEVEL="$2"; shift 2;;
    -h|--help)         usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

[[ -z "${SAMPLE}" ]] && { echo "ERROR: -s/--sample required"; exit 1; }
[[ -z "${CHROM}"  ]] && { echo "ERROR: -c/--chrom required"; exit 1; }
[[ -z "${PREP}"   ]] && { echo "ERROR: -i/--prepdir required"; exit 1; }
[[ -z "${OUT}"    ]] && { echo "ERROR: -o/--outdir required"; exit 1; }
[[ -d "${PREP}"   ]] || { echo "ERROR: prepdir not found: ${PREP}"; exit 1; }
mkdir -p "${OUT}"

# Resolve defaults from data dir if given (you can rename hg38 -> mm10 in your install)
if [[ -n "${DATA_DIR}" ]]; then
  : "${SEQ:=${DATA_DIR}/hg38.2bit}"
  : "${CHROMSIZES:=${DATA_DIR}/hg38.chrom.sizes}"
  : "${BLACKLIST:=${DATA_DIR}/hg38_maxatac_blacklist.bw}"
fi
[[ -n "${SEQ}"        ]] || { echo "ERROR: --sequence or -d DATA_DIR required"; exit 1; }
[[ -n "${CHROMSIZES}" ]] || { echo "ERROR: --chrom-sizes or -d DATA_DIR required"; exit 1; }
[[ -n "${BLACKLIST}"  ]] || { echo "ERROR: --blacklist or -d DATA_DIR required"; exit 1; }

# CTCF
TFs=( "CTCF" )


# Find signal bigWig if not explicitly provided
if [[ -z "${SIGNAL}" ]]; then
  candidates=()
  shopt -s nullglob
  for pat in "${PREP}/${SAMPLE}"*.bw "${PREP}"/*.bw; do
    [[ -f "$pat" ]] && candidates+=("$pat")
  done
  shopt -u nullglob
  if [[ ${#candidates[@]} -eq 0 ]]; then
    echo "ERROR: no .bw signal found in ${PREP}; pass --signal" >&2
    exit 1
  fi
  # pick largest by size (robust to spaces)
  max_size=-1; SIGNAL=""
  for f in "${candidates[@]}"; do
    sz=$(stat -c%s "$f" 2>/dev/null || stat -f%z "$f")
    if (( sz > max_size )); then max_size=$sz; SIGNAL="$f"; fi
  done
fi
[[ -f "${SIGNAL}" ]] || { echo "ERROR: signal not found: ${SIGNAL}"; exit 1; }

echo "[maxATAC-predict] sample=${SAMPLE} chrom=${CHROM}"
echo "[maxATAC-predict] prep=${PREP}"
echo "[maxATAC-predict] signal=${SIGNAL}"
echo "[maxATAC-predict] outdir=${OUT}"
echo "[maxATAC-predict] seq=${SEQ}"
echo "[maxATAC-predict] sizes=${CHROMSIZES}"
echo "[maxATAC-predict] blacklist=${BLACKLIST}"
echo "[maxATAC-predict] TFs=${TFs[*]}"
echo "[maxATAC-predict] threads=${THREADS} loglevel=${LOGLEVEL}"

for TF in "${TFs[@]}"; do
  maxatac predict \
    -n "${SAMPLE}_${TF}" \
    -tf "${TF}" \
    --signal "${SIGNAL}" \
    --sequence "${SEQ}" \
    --chrom_sizes "${CHROMSIZES}" \
    --blacklist "${BLACKLIST}" \
    --chromosomes "${CHROM}" \
    --output "${OUT}" \
    --loglevel "${LOGLEVEL}"
done