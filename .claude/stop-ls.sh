#!/usr/bin/env bash
# Stop hook: re-injects the LS-port GOAL + the doctrine until the LS port is complete.
# Off-switch: touch /workspace/FVSjl/docs/LS_COMPLETE
set -euo pipefail
DONE=/workspace/FVSjl/docs/LS_COMPLETE
GOAL=/workspace/FVSjl/docs/LS_GOAL.md
[ -f "$DONE" ] && exit 0
[ -f "$GOAL" ] || exit 0
# exit 2 + stderr => Claude Code blocks the stop and feeds stderr back as a reminder
{
  echo "=== ACTIVE GOAL REMINDER (Stop hook) — LS (Lake States) variant port ==="
  cat "$GOAL"
  echo "=== If the LS port is genuinely complete, run: touch $DONE ==="
} >&2
exit 2
