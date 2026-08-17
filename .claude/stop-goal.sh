#!/usr/bin/env bash
# Active campaign: FVS → FVSjl port. Western variant cluster + all major extensions are DONE.
# Goal: docs/WESTERN_ROLLOUT_GOAL.md ; off-switch: touch docs/WESTERN_ROLLOUT_COMPLETE.
DONE=/workspace/FVSjl/docs/WESTERN_ROLLOUT_COMPLETE
GOAL=/workspace/FVSjl/docs/WESTERN_ROLLOUT_GOAL.md
[ -f "$DONE" ] && exit 0
echo "=== ACTIVE GOAL (Stop hook) — FVS→FVSjl port: post-WRD backlog ===" >&2
cat "$GOAL" >&2
echo "=== DISCIPLINE (unchanged): validate BIT-EXACT-OR-CORNERED vs the LIVE relinked FVS oracle per chunk; MEASURE don't infer (instrument the Fortran, g16 single-.o-swap dump-replay, instrumented .sum byte-identical before trusting); per-record treelist diff INVALID after tripling (use .sum aggregates / pre-split window); reuse the shared engine; never FFI an RNG; commit ONLY validated chunks; verify-by-run on the MERGED tree after every merge; block-on-signal via a foreground-blocking Bash wait (NOT idle stop-hook cycling). ===" >&2
echo "=== STATUS 2026-08-17 (kt-variant-port HEAD 8139b42d): ALL ~22 US variants + BC growth+volume DONE; FFE/ECON/Climate-FVS/base+spatial-mistletoe(NEWSPRED #196)/ORGANON/DBS/SVS-data-path DONE; WRD (rd/) COMPLETE + live-validated end-to-end across all 15 base-rd/ variants (RD suite 1159/1159, rd−ctrl .sum DELTA cornered vs each FVS<v>_clean); OP/SVS(visualization) loose-ends CLOSED. #196 BC NEWSPRED in-tree + closed-cornered (DM matches oracle's minor ~113 TPH; YSM residual = 2 BC DB-reader bugs FIXED, growth #206 straddle). ===" >&2
echo "=== NEXT (USER ROADMAP 2026-08-17, work in ORDER): (1) DMNTRD (dmntrd.f 202ln — NEWSPRED cyc≥2 tripled-offspring infection remap; gated off newspred.jl:672; closes #196). (2) INSECT/PATHOGEN family (~267 f, keyword-no-ops in FVSjl = feature gap): pilot dfb(29 DF Beetle) end-to-end → dftm(31)/wpbr(43)/wwpb(53)/lpmpb(54)/wsbwe(57). (3) COVER understory covr(11)+vcovr(3)+pg(14). (4) ON (canada/on) Ontario/Penner variant. Off-switch (USER's call): touch docs/WESTERN_ROLLOUT_COMPLETE ===" >&2
exit 2
