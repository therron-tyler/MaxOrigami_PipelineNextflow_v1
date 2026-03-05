#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-.}"          # where to search; default is current dir
MODE="${2:-dryrun}"     # 'dryrun' (default) or 'delete'

# Collect matches safely (handles spaces/newlines)
mapfile -d '' MATCHES < <(
  find "$ROOT" -type d -path '*/maxatac/predict/*/prepare' -print0
)

if ((${#MATCHES[@]} == 0)); then
  echo "No nested prepare dirs found under: $ROOT"
  exit 0
fi

echo "Found ${#MATCHES[@]} nested prepare dirs under: $ROOT"
printf ' - %s\n' "${MATCHES[@]}"

if [[ "$MODE" != "delete" ]]; then
  echo
  echo "Dry run only. To actually delete, run:"
  echo "  $0 \"$ROOT\" delete"
  exit 0
fi

echo
echo "Deleting..."
# Delete each directory (null-delimited for safety)
printf '%s\0' "${MATCHES[@]}" | xargs -0 rm -rf
echo "Done."
