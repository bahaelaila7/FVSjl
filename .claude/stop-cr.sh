#!/usr/bin/env bash
DONE=/workspace/FVSjl/docs/CR_VARIANT_PORT_COMPLETE
GOAL=/workspace/FVSjl/docs/CR_VARIANT_PORT_GOAL.md
[ -f "$DONE" ] && exit 0
echo "=== ACTIVE GOAL (Stop hook) — port the CR (Central Rockies) variant ===" >&2
cat "$GOAL" >&2
echo "=== DISCIPLINE: validate bit-exact vs LIVE FVScr per chunk; MEASURE don't infer; per-record treelist diff INVALID after tripling (use pre-split window / stand .sum); reuse the shared engine. Off: touch docs/CR_VARIANT_PORT_COMPLETE ===" >&2
exit 2
