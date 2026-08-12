#!/usr/bin/env bash
# Active campaign: WESTERN FVS variant cluster + extensions rollout (NOT a single variant).
# Goal: docs/WESTERN_ROLLOUT_GOAL.md ; off-switch: touch docs/WESTERN_ROLLOUT_COMPLETE.
DONE=/workspace/FVSjl/docs/WESTERN_ROLLOUT_COMPLETE
GOAL=/workspace/FVSjl/docs/WESTERN_ROLLOUT_GOAL.md
[ -f "$DONE" ] && exit 0
echo "=== ACTIVE GOAL (Stop hook) — Western FVS variant cluster + extensions rollout ===" >&2
cat "$GOAL" >&2
echo "=== DISCIPLINE: validate bit-exact vs LIVE FVS per chunk; MEASURE don't infer; per-record treelist diff INVALID after tripling (use .sum / pre-split); reuse the shared engine. Do NOT stop at one variant. DONE-or-cornered: #140(BM), #137(EM self-thin), EM-sub-inch-DG(crown-dub), #191(TT aspen), #143(IE pvref1 + EM AUTOES re-measured cornered +4-10% TPA, BA/TopHt match), #194(CI esgent birth-cycle fixed, transition cornered), CI-vol(bit-exact @cyc0), Climate-FVS. WHOLE WESTERN CLUSTER growth+vol bit-exact-or-cornered. #196: BC metric-DATABASE units bug FOUND+FIXED (42eb555, cyc0 bit-exact TPA/SDI/TopHt/QMD); REMAINING = YSM multi-cycle under-mortalization ISOLATED to MISSING dwarf mistletoe: DM-free V3 control mrun/all_BC.key shows jl's BC V3 self-thin MATCHES oracle (TPA 2089→1253 vs 2087→1292), so YSM extra mort = the unported NEWSPRED spatial DM (canada/newmist ~50 rtns incl DMAUTO) + BC not in DM dispatch. FIX = port newmist NEWSPRED + wire BC (a genuine extension port). base mistoe.f done N-Rockies+CR. + BC V2/non-ICH. Off (USER's call): touch docs/WESTERN_ROLLOUT_COMPLETE ===" >&2
exit 2
