#!/usr/bin/env bash
# Stop hook: re-injects the FVSjl modernization GOAL + doctrine until all FOUR pillars are met —
# (1) 100% drop-in correctness/coverage for SN+NE+CS+LS at the bit-exact-or-cornered floor,
# (2) allocation-free memory path, (3) SoA + massively-parallel (serial==parallel bit-identical),
# (4) idiomatic/maintainable/type-stable. Every change must keep the tolerance-campaign bit-exact floor.
# Off-switch: touch /workspace/FVSjl/docs/MODERNIZATION_COMPLETE
set -euo pipefail
DONE=/workspace/FVSjl/docs/MODERNIZATION_COMPLETE
GOAL=/workspace/FVSjl/docs/MODERNIZATION_GOAL.md
[ -f "$DONE" ] && exit 0
[ -f "$GOAL" ] || exit 0
# exit 2 + stderr => Claude Code blocks the stop and feeds stderr back as a reminder
{
  echo "=== ACTIVE GOAL REMINDER (Stop hook) — FVSjl modernization campaign ==="
  cat "$GOAL"
  echo ""
  echo "=== Working checklist: docs/MODERNIZATION_AUDIT.md (tick each slice: bit-exact re-verified + metric) ==="
  echo "=== CORRECTNESS FLOOR (never regress): suite green + tolerance-campaign bit-exact-or-cornered state ==="
  echo "=== When ALL FOUR pillars are met and every slice is bit-exact-verified, run: touch $DONE ==="
} >&2
exit 2
