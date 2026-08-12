#!/usr/bin/env bash
# Active campaign: WESTERN FVS variant cluster + extensions rollout (NOT a single variant).
# Goal: docs/WESTERN_ROLLOUT_GOAL.md ; off-switch: touch docs/WESTERN_ROLLOUT_COMPLETE.
DONE=/workspace/FVSjl/docs/WESTERN_ROLLOUT_COMPLETE
GOAL=/workspace/FVSjl/docs/WESTERN_ROLLOUT_GOAL.md
[ -f "$DONE" ] && exit 0
echo "=== ACTIVE GOAL (Stop hook) — Western FVS variant cluster + extensions rollout ===" >&2
cat "$GOAL" >&2
echo "=== DISCIPLINE: validate bit-exact vs LIVE FVS per chunk; MEASURE don't infer; per-record treelist diff INVALID after tripling (use .sum / pre-split); reuse the shared engine. Do NOT stop at one variant. DONE-or-cornered: #140(BM), #137(EM self-thin), EM-sub-inch-DG(crown-dub), #191(TT aspen), #143(IE pvref1 + EM AUTOES re-measured cornered +4-10% TPA, BA/TopHt match), #194(CI esgent birth-cycle fixed, transition cornered), CI-vol(bit-exact @cyc0), Climate-FVS. WHOLE WESTERN CLUSTER growth+vol bit-exact-or-cornered. #196: BC metric-DATABASE units bug FOUND+FIXED (42eb555, cyc0 bit-exact TPA/SDI/TopHt/QMD); REMAINING = YSM multi-cycle under-mortalization, attribution NOT isolated (candidates: NEWSPRED spatial DM unported+jl-mistletoe-inert, OR BC BAMAX self-thin under-kill) — oracle A/B BLOCKED (FVSbc_clean SIGSEGV dbstreesin.f:57 on DATABASE); next = inline-TREEDATA FVSbc with/without MISTOE. base mistoe.f done N-Rockies+CR, NEWSPRED+BC-DM unported. + BC V2/non-ICH. Off (USER's call): touch docs/WESTERN_ROLLOUT_COMPLETE ===" >&2
exit 2
