#!/usr/bin/env bash
DONE=/workspace/FVSjl/docs/FVS_LARGEST_DIV_COMPLETE
GOAL=/workspace/FVSjl/docs/FVS_LARGEST_DIV_GOAL.md
[ -f "$DONE" ] && exit 0
echo "=== ACTIVE GOAL (Stop hook) — verify LARGEST FIA divergences irreducible ===" >&2
cat "$GOAL" >&2
echo "=== DISCIPLINE: per-record treelist diff INVALID after tripling; use stand-level .sum or pre-split window. Off: touch docs/FVS_LARGEST_DIV_COMPLETE ===" >&2
exit 2
