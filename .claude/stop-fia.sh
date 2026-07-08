#!/usr/bin/env bash
# Stop hook: re-injects the FVSjl FIA/FVS behaviour-compatibility GOAL + doctrine until all FOUR pillars
# are met — (1) scaled stratified real-FIA sample per variant, (2) multi-cycle projection bit-exact-or-
# cornered vs live FVS, (3) management-scenario compatibility, (4) divergence taxonomy fully cornered —
# without regressing the closed campaigns' floor (38527 pass / 143 broken / 0 fail).
# Off-switch: touch /workspace/FVSjl/docs/FIA_FVS_COMPAT_COMPLETE
set -euo pipefail
DONE=/workspace/FVSjl/docs/FIA_FVS_COMPAT_COMPLETE
GOAL=/workspace/FVSjl/docs/FIA_FVS_COMPAT_GOAL.md
[ -f "$DONE" ] && exit 0
[ -f "$GOAL" ] || exit 0
# exit 2 + stderr => Claude Code blocks the stop and feeds stderr back as a reminder
{
  echo "=== ACTIVE GOAL REMINDER (Stop hook) — FVSjl FIA/FVS behaviour-compatibility campaign ==="
  cat "$GOAL"
  echo ""
  echo "=== Working checklist: docs/FIA_FVS_COMPAT_AUDIT.md (tick each slice: plots covered, pass rate, divergences cornered/fixed) ==="
  echo "=== CORRECTNESS FLOOR (never regress): suite 38527/143/0 + tolerance-campaign bit-exact-or-cornered state ==="
  echo "=== When all four pillars are met and every divergence is bit-exact-or-cornered, run: touch $DONE ==="
} >&2
exit 2
