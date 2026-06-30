#!/usr/bin/env bash
# Stop hook: re-injects the CS-port GOAL + the doctrine until the CS port is complete.
# Off-switch: touch /workspace/FVSjl/docs/CS_COMPLETE
set -euo pipefail
DONE=/workspace/FVSjl/docs/CS_COMPLETE
GOAL=/workspace/FVSjl/docs/CS_GOAL.md
[ -f "$DONE" ] && exit 0
[ -f "$GOAL" ] || exit 0
# exit 2 + stderr => Claude Code blocks the stop and feeds stderr back as a reminder
{
  echo "=== ACTIVE GOAL REMINDER (Stop hook) — CS variant port ==="
  cat "$GOAL"
  echo "=== If the CS port is genuinely complete, run: touch $DONE ==="
} >&2
exit 2
