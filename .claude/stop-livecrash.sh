#!/usr/bin/env bash
# Stop hook: re-injects the "fix live-FVS crashes" GOAL + doctrine until every crash class is root-caused,
# minimally patched, documented, and the crashing stands run clean. Off-switch: touch docs/FVS_LIVECRASH_COMPLETE
set -euo pipefail
DONE=/workspace/FVSjl/docs/FVS_LIVECRASH_COMPLETE
GOAL=/workspace/FVSjl/docs/FVS_LIVECRASH_FIX_GOAL.md
[ -f "$DONE" ] && exit 0
[ -f "$GOAL" ] || exit 0
{
  echo "=== ACTIVE GOAL REMINDER (Stop hook) — fix live-FVS crashes on real FIA (minimal upstream patches) ==="
  cat "$GOAL"
  echo ""
  echo "=== Working log: docs/FVS_LIVECRASH_AUDIT.md ; patches: docs/patches/ ; bug registry: docs/FVS_SOURCE_BUGS.md ==="
  echo "=== DISCIPLINE: instrument .f -> build -> run -> RESTORE source + clean .o -> verify oracle pristine. FVSjl runs these clean = correct side. ==="
  echo "=== When all 6 crash classes are patched + documented + the 12 stands run clean, run: touch $DONE ==="
} >&2
exit 2
