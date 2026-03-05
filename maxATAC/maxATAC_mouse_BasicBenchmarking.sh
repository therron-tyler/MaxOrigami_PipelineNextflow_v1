#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage:
  $(basename "$0") -s SAMPLE -l PRED_DIRS.txt -o OUT_DIR \
                   -a ATAC.bw -b BLACKLIST.bed -c CHROM.sizes \
                   [-f TF_NAME] [-C chr1,chr2,...] [-t THRESHOLD]

Computes enrichment: fraction of TF peaks inside ATAC regions vs random expectation.

Required:
  -s, --sample        Sample name
  -l, --list          Text file: one per-chrom maxATAC prediction directory per line
  -o, --outdir        Output directory
  -a, --atac          ATAC signal bigWig (from maxATAC prepare)
  -b, --blacklist     BED of blacklisted regions
  -c, --chrom-sizes   Chrom sizes file (e.g., mm10.chrom.sizes)

Optional:
  -f, --tf            TF name (default: CTCF) -> expects <SAMPLE>_<TF>_peaks.bed in each pred dir
  -C, --chroms        Comma or space separated chroms (default: all in chrom.sizes)
  -t, --threshold     ATAC min-max value threshold (default: 0.30)

Outputs:
  OUT_DIR/ATAC.narrowPeak
  OUT_DIR/TF_peaks.filtered.bed
  OUT_DIR/benchmark.txt
EOF
}

SAMPLE=""
LIST=""
OUT=""
ATAC=""
BL=""
CS=""
TF="CTCF"
CHROMS=""
THRESH="0.30"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--sample)        SAMPLE="$2"; shift 2;;
    -l|--list)          LIST="$2"; shift 2;;
    -o|--outdir)        OUT="$2"; shift 2;;
    -a|--atac)          ATAC="$2"; shift 2;;
    -b|--blacklist)     BL="$2"; shift 2;;
    -c|--chrom-sizes)   CS="$2"; shift 2;;
    -f|--tf)            TF="$2"; shift 2;;
    -C|--chroms)        CHROMS="$2"; shift 2;;
    -t|--threshold)     THRESH="$2"; shift 2;;
    -h|--help)          usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

[[ -n "$SAMPLE" && -n "$LIST" && -n "$OUT" && -n "$ATAC" && -n "$BL" && -n "$CS" ]] || { usage; exit 1; }
[[ -s "$ATAC"  ]] || { echo "[error] ATAC bw not found: $ATAC" >&2; exit 1; }
[[ -s "$BL"    ]] || { echo "[error] blacklist not found: $BL" >&2; exit 1; }
[[ -s "$CS"    ]] || { echo "[error] chrom.sizes not found: $CS" >&2; exit 1; }
mkdir -p "$OUT"

# Build chrom filter regex (or take all)
CHR_RE=".*"
if [[ -n "$CHROMS" ]]; then
  CHROMS="${CHROMS//,/ }"
  # Turn "chr1 chr2" into "^(chr1|chr2)$"
  RE="^("
  for c in $CHROMS; do RE+="${c}|"; done
  CHR_RE="${RE%|})$"
fi

# 1) Derive ATAC narrowPeak across requested chroms
#    Requires: bigWigToBedGraph, bedtools
bigWigToBedGraph "$ATAC" stdout \
| awk -v re="$CHR_RE" -v t="$THRESH" 'BEGIN{OFS="\t"} $1 ~ re && $4>=t {print $1,$2,$3}' \
| bedtools sort -i - \
| bedtools merge -i - \
| bedtools intersect -v -a - -b "$BL" \
> "$OUT/ATAC.narrowPeak"

# 2) Collect TF peak beds from each per-chrom prediction dir; filter to chroms
#    Expected filename in each dir: <SAMPLE>_<TF>_peaks.bed
TMP_TF="$OUT/TF_peaks.all.bed"
: > "$TMP_TF"
while IFS= read -r d; do
  [[ -z "$d" ]] && continue
  BED="${d%/}/${SAMPLE}_${TF}_peaks.bed"
  [[ -s "$BED" ]] || { echo "[warn] missing bed in $d: $(basename "$BED")"; continue; }
  cat "$BED" >> "$TMP_TF"
done < "$LIST"

[[ -s "$TMP_TF" ]] || { echo "[error] no TF peaks collected (TF=$TF)"; exit 1; }

awk -v re="$CHR_RE" '$1 ~ re' "$TMP_TF" > "$OUT/TF_peaks.filtered.bed"

# 3) Stats
a_in_b=$(bedtools intersect -u -a "$OUT/TF_peaks.filtered.bed" -b "$OUT/ATAC.narrowPeak" | wc -l | awk '{print $1}')
a_all=$(wc -l < "$OUT/TF_peaks.filtered.bed")

# genome length over chosen chroms
L=$(awk -v re="$CHR_RE" '$1 ~ re {sum+=$2} END{print (sum?sum:0)}' "$CS")
# ATAC coverage over chosen chroms
covered=$(awk '{sum+=($3-$2)} END{print (sum?sum:0)}' "$OUT/ATAC.narrowPeak")

# fractions and pretty-printed percentages
obs=$(awk -v a="$a_in_b" -v b="$a_all"     'BEGIN{ if(b==0) print 0; else print a/b }')
exp=$(awk -v c="$covered" -v L="$L"        'BEGIN{ if(L==0) print 0; else print c/L }')
cov_pct=$(awk -v c="$covered" -v L="$L"    'BEGIN{ printf("%.2f", (L==0?0:100.0*c/L)) }')
obs_pct=$(awk -v o="$obs"                  'BEGIN{ printf("%.2f", 100.0*o) }')
exp_pct=$(awk -v e="$exp"                  'BEGIN{ printf("%.2f", 100.0*e) }')
enrich_line=$(awk -v o="$obs" -v e="$exp"  'BEGIN{
  if(e==0) printf("Observed: %.2f%%  Expected: 0%%  Enrichment: n/a\n", 100*o);
  else     printf("Observed: %.2f%%  Expected: %.2f%%  Enrichment: %.2fx\n", 100*o, 100*e, o/e);
}')

{
  printf "Sample: %s\n" "$SAMPLE"
  printf "TF: %s\n" "$TF"
  printf "Chromosomes: %s\n" "${CHROMS:-ALL_FROM_${CS}}"
  printf "ATAC threshold: %s\n" "$THRESH"
  printf "TF peaks in ATAC: %s / %s\n" "$a_in_b" "$a_all"
  printf "ATAC coverage: %s%%\n" "$cov_pct"
  printf "%s" "$enrich_line"
} > "$OUT/benchmark.txt"

echo "[done] Wrote: $OUT/benchmark.txt"
