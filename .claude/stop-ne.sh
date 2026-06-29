#!/usr/bin/env bash
# Stop hook: re-injects the NE-port GOAL + the doctrine until the NE port is complete.
# Off-switch: touch /workspace/FVSjl/docs/NE_COMPLETE
set -euo pipefail
DONE=/workspace/FVSjl/docs/NE_COMPLETE
GOAL=/workspace/FVSjl/docs/NE_GOAL.md
[ -f "$DONE" ] && exit 0
[ -f "$GOAL" ] || exit 0
# exit 2 + stderr => Claude Code blocks the stop and feeds stderr back as a reminder
{
  echo "=== ACTIVE GOAL REMINDER (Stop hook) — NE variant port ==="
  cat "$GOAL"
  echo "=== If the NE port is genuinely complete, run: touch $DONE ==="
} >&2
exit 2
