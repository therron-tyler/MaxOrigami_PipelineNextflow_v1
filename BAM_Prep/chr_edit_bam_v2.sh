#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") -i INPUT.bam -o OUTPUT.chr.bam [-t THREADS] [--no-index]
                        [--downsample N] [--mapq Q] [--exclude-flags INT] [--seed INT]

Adds 'chr' prefixes to @SQ SN: fields in a BAM header (M/MT -> chrM).
Optionally downsamples the BAM to approximately N eligible alignments using samtools.

Downsampling:
  --downsample N        Target number of eligible alignments to KEEP (default: 0 = disabled)
  --mapq Q              MAPQ filter for counting/downsampling (default: 30)
  --exclude-flags INT   Exclude flags for counting/downsampling (default: 0xF00)
                        (0x100 secondary + 0x200 QCfail + 0x400 dup + 0x800 supplementary)
  --seed INT            Seed for samtools -s (default: 42)

Other:
  -t, --threads         Threads for sort/index (default: 4)
  --no-index            Skip indexing the output BAM
  -h, --help            Show this help
EOF
}

IN=""
OUT=""
THREADS=4
DO_INDEX=1

DOWNSAMPLE=0
MAPQ=30
EXCL_FLAGS=0xF00
SEED=42

# --- parse args ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--in)           IN="$2"; shift 2;;
    -o|--out)          OUT="$2"; shift 2;;
    -t|--threads)      THREADS="$2"; shift 2;;
    --no-index)        DO_INDEX=0; shift 1;;
    --downsample)      DOWNSAMPLE="$2"; shift 2;;
    --mapq)            MAPQ="$2"; shift 2;;
    --exclude-flags)   EXCL_FLAGS="$2"; shift 2;;
    --seed)            SEED="$2"; shift 2;;
    -h|--help)         usage; exit 0;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1;;
  esac
done

[[ -z "${IN}"  ]] && { echo "ERROR: -i/--in required" >&2; usage; exit 1; }
[[ -z "${OUT}" ]] && { echo "ERROR: -o/--out required" >&2; usage; exit 1; }
[[ -f "${IN}"  ]] || { echo "ERROR: input not found: ${IN}" >&2; exit 1; }

OUTDIR="$(dirname -- "${OUT}")"
mkdir -p "${OUTDIR}"

TMP_HDR_ORIG="$(mktemp -p "${OUTDIR}" .hdr.orig.XXXXXX)"
TMP_HDR_EDIT="$(mktemp -p "${OUTDIR}" .hdr.edit.XXXXXX)"
TMP_REHDR="${OUT}.rehdr.tmp.bam"
TMP_DS="${OUT}.ds.tmp.bam"
TMP_SORT="${OUT}.sort.tmp.bam"

cleanup() { rm -f "${TMP_HDR_ORIG}" "${TMP_HDR_EDIT}" "${TMP_REHDR}" "${TMP_DS}" "${TMP_SORT}"; }
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
      if(sn !~ /^chr/){
        if(sn=="MT" || sn=="M"){ sn="chrM" }
        else if(sn ~ /^(X|Y|[0-9]+)$/){ sn="chr" sn }
      }
      $i="SN:" sn
    }
  }
  print; next
}
{print}
' "${TMP_HDR_ORIG}" > "${TMP_HDR_EDIT}"

# 3) reheader → tmp BAM
samtools reheader "${TMP_HDR_EDIT}" "${IN}" > "${TMP_REHDR}"

# 4) optional downsample
if [[ "${DOWNSAMPLE}" -gt 0 ]]; then
  echo "[downsample] target=${DOWNSAMPLE} mapq=${MAPQ} excl_flags=${EXCL_FLAGS} seed=${SEED}"

  # Count eligible alignments (same filters used for sampling)
  n_total=$(samtools view -c -q "${MAPQ}" -F "${EXCL_FLAGS}" "${TMP_REHDR}")
  echo "[downsample] eligible_alignments=${n_total}"

  if [[ "${n_total}" -le "${DOWNSAMPLE}" ]]; then
    echo "[downsample] eligible <= target; skipping downsample"
    cp "${TMP_REHDR}" "${TMP_DS}"
  else
    # probability = target / total
    prob=$(python - <<PY
n=${n_total}
t=${DOWNSAMPLE}
print(t/float(n))
PY
)
    # samtools -s expects SEED.FRAC (frac ~ prob); format with 6 digits
    sarg=$(python - <<PY
seed=${SEED}
p=float("${prob}")
print(f"{seed}.{int(p*1e6):06d}")
PY
)
    echo "[downsample] samtools -s ${sarg} (prob=${prob})"

    samtools view -@ "${THREADS}" -b -q "${MAPQ}" -F "${EXCL_FLAGS}" -s "${sarg}" "${TMP_REHDR}" > "${TMP_DS}"
    # NOTE: this outputs only alignments passing filters; if you want to keep ALL reads but subsample,
    # remove -q/-F here and only use them for counting.
  fi
else
  cp "${TMP_REHDR}" "${TMP_DS}"
fi

# 5) sort to ensure indexable + consistent output
samtools sort -@ "${THREADS}" -o "${TMP_SORT}" "${TMP_DS}"
mv -f "${TMP_SORT}" "${OUT}"

# 6) index (optional)
if [[ "${DO_INDEX}" -eq 1 ]]; then
  samtools index -@ "${THREADS}" "${OUT}"
fi

# 7) quick sanity
echo "Wrote: ${OUT}"
echo "First few contigs in ${OUT}:"
samtools idxstats "${OUT}" | awk 'NR<=10{print $1}'
