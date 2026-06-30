#!/usr/bin/env bash
# Stop hook: re-injects the divergence-fix-campaign GOAL + doctrine until complete.
# Off-switch: touch /workspace/FVSjl/docs/DIVERGENCE_COMPLETE
set -euo pipefail
DONE=/workspace/FVSjl/docs/DIVERGENCE_COMPLETE
GOAL=/workspace/FVSjl/docs/DIVERGENCE_GOAL.md
[ -f "$DONE" ] && exit 0
[ -f "$GOAL" ] || exit 0
{
  echo "=== ACTIVE GOAL REMINDER (Stop hook) — non-ULP divergence fix campaign ==="
  cat "$GOAL"
  echo "=== Ledger: docs/DIVERGENCE_FIX_CAMPAIGN.md · off-switch: touch $DONE ==="
} >&2
exit 2
