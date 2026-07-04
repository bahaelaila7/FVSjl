#!/usr/bin/env bash
# Stop hook: re-injects the test-tolerance closure GOAL + the doctrine until every tolerance is
# either BIT-EXACT or PROVEN-IRREDUCIBLE-ULP (cornered to the exact Float32 op).
# Off-switch: touch /workspace/FVSjl/docs/TOLERANCE_COMPLETE
set -euo pipefail
DONE=/workspace/FVSjl/docs/TOLERANCE_COMPLETE
GOAL=/workspace/FVSjl/docs/TOLERANCE_GOAL.md
[ -f "$DONE" ] && exit 0
[ -f "$GOAL" ] || exit 0
# exit 2 + stderr => Claude Code blocks the stop and feeds stderr back as a reminder
{
  echo "=== ACTIVE GOAL REMINDER (Stop hook) — test-tolerance closure campaign ==="
  cat "$GOAL"
  echo ""
  echo "=== Working checklist: docs/TOLERANCE_AUDIT.md (tick each ⬜ → BIT-EXACT or ULP+root) ==="
  echo "=== When the checklist is 100% closed, run: touch $DONE ==="
} >&2
exit 2
