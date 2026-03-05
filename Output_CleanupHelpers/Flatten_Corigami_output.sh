#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-.}"          # search root (default: current dir)
MODE="${2:-dryrun}"     # 'dryrun' (default) or 'apply'

# Find: */corigami/origami_chr*/<sample>/prediction
find "$ROOT" -type d -name prediction -path '*/corigami/origami_chr*/*/prediction' -print0 |
while IFS= read -r -d '' pred_dir; do
  sample_dir="$(dirname "$pred_dir")"       # .../origami_chrN/<sample>
  base_dir="$(dirname "$sample_dir")"       # .../origami_chrN

  # Simpler/accurate sanity check for your tree
  case "$base_dir" in
    *"/corigami/origami_chr"*) ;;           # OK: matches corigami/origami_chr*
    *) echo "SKIP (unexpected path): $pred_dir"; continue;;
  esac

  # Items directly under prediction/
  mapfile -d '' items < <(find "$pred_dir" -mindepth 1 -maxdepth 1 -print0)

  if ((${#items[@]} == 0)); then
    echo "EMPTY prediction dir, nothing to move: $pred_dir"
    continue
  fi

  echo
  echo "Will move ${#items[@]} item(s):"
  printf ' - %s\n' "${items[@]}"
  echo "From: $pred_dir"
  echo "  To: $base_dir"

  if [[ "$MODE" != "apply" ]]; then
    echo "Mode=dryrun (no changes). To apply: $0 \"$ROOT\" apply"
    continue
  fi

  # Merge-safe move: copies contents up, preserves perms, merges dirs
  mkdir -p "$base_dir"
  rsync -a "$pred_dir"/ "$base_dir"/

  # Remove source and clean empty wrapper
  rm -rf "$pred_dir"
  rmdir "$sample_dir" 2>/dev/null || true

  echo "Moved and cleaned: $pred_dir -> $base_dir"
done

echo
echo "Done. Mode: $MODE"
