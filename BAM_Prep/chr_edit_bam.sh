#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") -i INPUT.bam -o OUTPUT.chr.bam [-t THREADS] [--no-index]

Adds 'chr' prefixes to @SQ SN: fields in a BAM header (M/MT -> chrM).
Leaves scaffolds/haplotypes already containing 'chr' or non-standard names unchanged.

  -i, --in         Input BAM
  -o, --out        Output BAM (recommended: <sample>.chr.bam)
  -t, --threads    Threads for samtools index (default: 4)
      --no-index   Skip indexing the output BAM
  -h, --help       Show this help
EOF
}

IN=""
OUT=""
THREADS=4
DO_INDEX=1

# --- parse args ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--in)      IN="$2"; shift 2;;
    -o|--out)     OUT="$2"; shift 2;;
    -t|--threads) THREADS="$2"; shift 2;;
    --no-index)   DO_INDEX=0; shift 1;;
    -h|--help)    usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

[[ -z "${IN}"  ]] && { echo "ERROR: -i/--in required" >&2; usage; exit 1; }
[[ -z "${OUT}" ]] && { echo "ERROR: -o/--out required" >&2; usage; exit 1; }
[[ -f "${IN}"  ]] || { echo "ERROR: input not found: ${IN}" >&2; exit 1; }

# Ensure output dir exists; temp files live next to OUT for HPC friendliness
OUTDIR="$(dirname -- "${OUT}")"
mkdir -p "${OUTDIR}"
TMP_HDR_ORIG="$(mktemp -p "${OUTDIR}" .hdr.orig.XXXXXX)"
TMP_HDR_EDIT="$(mktemp -p "${OUTDIR}" .hdr.edit.XXXXXX)"
TMP_OUT="${OUT}.tmp"

cleanup() { rm -f "${TMP_HDR_ORIG}" "${TMP_HDR_EDIT}" "${TMP_OUT}"; }
trap cleanup EXIT

# 1) extract header
samtools view -H "${IN}" > "${TMP_HDR_ORIG}"

# 2) rewrite @SQ SN: fields → add chr, handle M/MT
awk '
BEGIN{OFS="\t"}
/^@SQ/{
  for(i=1;i<=NF;i++){
    if($i ~ /^SN:/){
      split($i,a,":"); sn=a[2];
      # already has chr -> leave unchanged
      if(sn !~ /^chr/){
        if(sn=="MT" || sn=="M"){ sn="chrM" }
        else if(sn ~ /^(X|Y|[0-9]+)$/){ sn="chr" sn }
        # other contigs (GL/KI/JH/…) left as-is
      }
      $i="SN:" sn
    }
  }
  print; next
}
{print}
' "${TMP_HDR_ORIG}" > "${TMP_HDR_EDIT}"

# 3) reheader → new BAM (write to tmp then mv for atomicity)
samtools reheader "${TMP_HDR_EDIT}" "${IN}" > "${TMP_OUT}"
mv -f "${TMP_OUT}" "${OUT}"

# 4) index (optional; Nextflow also checks and indexes if missing)
if [[ "${DO_INDEX}" -eq 1 ]]; then
  samtools index -@ "${THREADS}" "${OUT}"
fi

# 5) quick sanity
echo "Wrote: ${OUT}"
echo "First few contigs in ${OUT}:"
samtools idxstats "${OUT}" | awk 'NR<=10{print $1}'
