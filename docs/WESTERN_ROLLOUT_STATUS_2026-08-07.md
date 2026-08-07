# Western FVS variant cluster + extensions — consolidated status (2026-08-07)

Handoff/assessment for the off-switch decision (`touch docs/WESTERN_ROLLOUT_COMPLETE` — USER's call).
This is an ASSESSMENT, not a done-flag. The rollout is at the "bit-exact-or-cornered" bar; the residuals below
are each either an accepted cornered primitive or a scoped, justified deferral.

## Growth + volume — at the bar
CR ✓ COMPLETE · KT/IE/EM/BM/TT/UT ✓ bit-exact-or-cornered · BC growth+yield done (merch/board + V2 remain) ·
CI running, refinement tail cornered. Shared DG-calibration dispatch (PSIGSQ/bark/DGSD) audited COMPLETE.

## Extensions matrix — CLOSED (all accounted for)
- FFE ✓ (all western + eastern + CR) · Dwarf mistletoe ✓ · ECON ✓.
- **Climate-FVS**: SCOPED + JUSTIFIED-DEFERRED (2026-08-07, docs/EXTENSIONS_ROLLOUT_PLAN.md). 1739-line subsystem,
  GUARDED-INERT (`clinit.f` LCLIMATE=.FALSE. default; `clgmult.f:60 IF(.NOT.LCLIMATE) RETURN`). Absence is
  invisible to ALL climate-off validation (the entire corpus). Port only when a climate-ON scenario is provided.

## FIA sweep — validated, 0 crashes
Whole-cluster multi-cycle validated (2026-08-03) + per-variant real-FIA re-validations (CI/IE 2026-08-05; BM
2026-08-07 post-habitat-fix, 0 crashes, no regression). Residuals = the accepted DGSCOR/density straddle.

## Remaining technical residuals (each cornered or scoped)
1. **#142 CI ~2% TPA**: the accepted ZZRAN/DGSCOR dense-regen RNG-realization straddle (deterministic DG bit-exact).
   MEETS the bar. (Open sub-item: volume MATW/FW2W · SMHTGF small-tree stochastic — unmeasured, low priority.)
2. **EM/IE growth-only ~0.5-0.8% large-tree DG tail**: CORNERED (mixed-sign, aggregate BA bit-exact; RDPSRT/AVHT40
   tie-break precision compounding). No fix warranted.
3. **#140 BM self-thin skew**: DECOMPOSED (2026-08-07) into TWO accepted cornered primitives —
   (a) self-thin COUNT-STRADDLE (374435108489998: density BIT-EXACT, TPAΔ0.9% ≈ 1 marginal tree of 66) and
   (b) growth-tail BA compounding (12827438010497: mortality BIT-EXACT, TPAΔ0.0%). BM meets the bar on density/BA.
   SOLE open question: whether the count-straddle's mild DIRECTIONAL LEAN (the 9:2 jl-under-thins skew) is an
   accepted RDPSRT tie-break lean or a small real kill-ORDER bias. Needs a cycle-gated per-tree self-thin KILL
   measurement on bmt01 AND/OR a larger clean-stand sign-tally (single measurements have hit 5 confounds:
   HISTORY=8 dead trees, DF-SMCON habitat default, .sum-inertness, dump calibration-vs-growth-pass state, tripling).
   The DG-adjacent bugs found during the hunt ARE fixed (see below).

## Real bugs fixed this session (all committed, validated vs live oracle)
- `d089b78` **#143** EM/IE AUTOES ingrowth over-establishment (253→0.1 TPA = live; measured via FVSem_g16).
- `04b15e6` + `7ce8f1f` **SDI-gate `tem` 35000-cap** cluster (EM/TT/UT): ultra-dense sub-1" cohort self-thin
  under-kill (background fallthrough). emt01/ttt01/utt01 TPA bit-exact; BM/SN/NE already capped.
- `66db612` **BM missing/unresolved habitat default** → KODTYP 79 (CWG113): DF small-tree SMCON corrected on
  habitat-less FIA stands. Source-faithful (habtyp.f:67-69), .sum-inert, bmt01 no-regress.

## Durable infrastructure added
- **FVSem_g16 / FVSbm_g16**: full gfortran-16 instrumentable oracles (all *.f + *.for + C objects, `build_g16.sh`);
  single-.o swap is ABI-safe within the g16 build. Bit-exact vs FVS{v}_clean. Reusable for any future measurement.

## Bottom line
The western cluster + extensions are at the bit-exact-or-cornered bar. Every non-bit-exact residual is a
measured, documented cornered primitive or a scoped/justified deferral (Climate-FVS). The single genuinely-open
technical question is #140's count-straddle directional lean — a marginal tie-break effect at the cornered/real
boundary that needs a dedicated study, not a blocker on the cluster's correctness. The off-switch is the user's call.

---
## CORRECTION (2026-08-07, later) — #140 is NOT cornered; REOPENED
An UNBIASED 20-stand BM stride-sample (vs the earlier hand-picked 4-large-stand sample) gives under=15/over=0/
mean Δ=+5.6% ⇒ #140 IS a REAL systematic under-thin bias (goal-doc verdict stands). The "#140 decomposed/cornered"
statements ABOVE are WITHDRAWN (commits 553a928 supersede 31493cb/679bfc4). #140 is OPEN: a real BM self-thin
under-thin whose root (large-tree DG ~1% under-shoot driving dq10-low→target-high, vs a self-thin kill-rate bias)
needs a cycle-gated per-tree measurement on a clean mid-size reproducer (41134262010497 +18.5%). The DG-adjacent
bugs fixed this session (habitat 66db612, tem-cap 04b15e6/7ce8f1f) remain valid but are .sum-inert on this bias.
The rest of the cluster status above is unaffected. LESSON: a systematic-vs-straddle sign-tally MUST use an
UNBIASED sample — selecting stands on divergence magnitude inverts the apparent sign distribution.
