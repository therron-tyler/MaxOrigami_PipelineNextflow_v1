#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-.}"          # search root (default: current dir)
MODE="${2:-dryrun}"     # 'dryrun' (default) or 'apply'

# Find .../maxatac/predict/<chr>/predict_<chr>
find "$ROOT" -type d -name 'predict_*' -path '*/maxatac/predict/*/predict_*' -print0 |
while IFS= read -r -d '' pred_dir; do
  chr_dir="$(dirname "$pred_dir")"          # .../maxatac/predict/<chr>
  chr_base="$(basename "$chr_dir")"         # <chr>  e.g., chr9
  leaf_base="$(basename "$pred_dir")"       # predict_<chr> e.g., predict_chr9
  expected="predict_${chr_base}"

  # Guard: only flatten if names match (prevents accidents)
  if [[ "$leaf_base" != "$expected" ]]; then
    echo "SKIP (name mismatch): $pred_dir  (expected leaf '$expected')"
    continue
  fi

  # What are we moving?
  mapfile -d '' items < <(find "$pred_dir" -mindepth 1 -maxdepth 1 -print0)
  if ((${#items[@]} == 0)); then
    echo "EMPTY predict dir, nothing to move: $pred_dir"
    # Try to remove empty dir anyway
    rmdir "$pred_dir" 2>/dev/null || true
    continue
  fi

  echo
  echo "Will move ${#items[@]} item(s):"
  printf ' - %s\n' "${items[@]}"
  echo "From: $pred_dir"
  echo "  To: $chr_dir"

  if [[ "$MODE" != "apply" ]]; then
    echo "Mode=dryrun (no changes). To apply: $0 \"$ROOT\" apply"
    continue
  fi

  # Merge-safe move (preserves perms/timestamps; merges directories)
  rsync -a "$pred_dir"/ "$chr_dir"/

  # Remove source directory after successful copy
  rm -rf "$pred_dir"

  echo "Flattened: $leaf_base -> $chr_dir"
done

echo
echo "Done. Mode: $MODE"
