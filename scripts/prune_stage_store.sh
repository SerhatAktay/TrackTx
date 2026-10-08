#!/usr/bin/env bash
# Remove superseded stage-store directories.
#
# Every change to a module, script or param gives that stage (and everything
# downstream) a new stamp, and the old directory stays so a past run can still be
# restored. This keeps only the newest directory per <module>/<label> and lists the
# rest. Dry run unless --apply is given.
#
# Usage: scripts/prune_stage_store.sh <output_dir> [--apply]
set -euo pipefail

out="${1:?usage: $0 <output_dir> [--apply]}"
apply="${2:-}"
store="$out/.store"
[[ -d "$store" ]] || { echo "No stage store at $store" >&2; exit 1; }

freed=0
for mod in "$store"/*/; do
  [[ "$(basename "$mod")" == _manifests ]] && continue
  # newest first; the first dir seen for each label is kept
  seen=""   # newline-separated labels (macOS ships bash 3.2: no associative arrays)
  while IFS= read -r d; do
    label="$(basename "$d" | sed -E 's/-[0-9a-f]{16}$//')"
    if grep -qxF -- "$label" <<< "$seen"; then
      kb=$(du -sk "$d" | cut -f1)
      freed=$((freed + kb))
      echo "old: $d (${kb} KB)"
      if [[ "$apply" == "--apply" ]]; then
        rm -rf "$d"
        rm -f "$store/_manifests/$(basename "$mod")/$(basename "$d").txt"
      fi
    else
      seen+="$label"$'\n'
    fi
  done < <(ls -dt "$mod"*/ 2>/dev/null | sed 's:/$::')
done

echo "$((freed / 1024)) MB in superseded directories$([[ "$apply" == "--apply" ]] && echo ' (removed)' || echo ' (dry run, pass --apply to remove)')"
