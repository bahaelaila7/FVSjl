#!/usr/bin/env bash
# Active campaign: WESTERN FVS variant cluster + extensions rollout (NOT a single variant).
# Goal: docs/WESTERN_ROLLOUT_GOAL.md ; off-switch: touch docs/WESTERN_ROLLOUT_COMPLETE.
DONE=/workspace/FVSjl/docs/WESTERN_ROLLOUT_COMPLETE
GOAL=/workspace/FVSjl/docs/WESTERN_ROLLOUT_GOAL.md
[ -f "$DONE" ] && exit 0
echo "=== ACTIVE GOAL (Stop hook) — Western FVS variant cluster + extensions rollout ===" >&2
cat "$GOAL" >&2
echo "=== DISCIPLINE: validate bit-exact vs LIVE FVS per chunk; MEASURE don't infer; per-record treelist diff INVALID after tripling (use .sum / pre-split); reuse the shared engine. Do NOT stop at one variant. Remaining: CI[#142]/EM[#137]/BM[#140]/IE-AUTOES[#143]/Climate-FVS. Off: touch docs/WESTERN_ROLLOUT_COMPLETE ===" >&2
exit 2
