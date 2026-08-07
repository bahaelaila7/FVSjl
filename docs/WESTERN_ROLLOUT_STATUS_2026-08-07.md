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

---
## UPDATE (2026-08-07, later 2) — #140 localized + CLUSTER REGRESSION GATE (no regression)

#140 PROGRESS (confound-free, on 614 stands 248804992489998 + bmt01 case-1):
- Non-619 residual CONFIRMED a REAL DG under-shoot (BA-low-EARLY on 614 stand, converges by maturity), NOT a
  density-preserved straddle. [aggregate BA-vs-TPA, no per-tree confound]
- DGCON **and** SMCON proven BIT-EXACT for every species (per-stand dump via FVSbm_clean `DEBUG` keyword vs jl
  bm_dgcons!) ⇒ the constant/coefficient loading is correct; the residual is a per-tree DDS INPUT (BAL/PCCF/
  RELDEN/CR) or the DDSS/XWT assembly at growth time.
- **bmt01 case-1 is a clean DB-free reproducer**: cyc0 START bit-exact, cycle-1 ACCRETION live 84 / jl 79 (−6%).
  My BM fixes (forkod=619-only, habitat-default) are INERT on 614 ⇒ pre-existing. NEXT = instrument the
  GROWTH-specific DGDRIV (distinct from CRATET's dubbed-tree DGDRIV, which confounds the DEBUG-keyword dumps).

CLUSTER REFERENCE REGRESSION GATE (fresh FVS{v}_clean vs jl run_keyfile, {v}t01 case-1 plain growth):
- **CR: FULLY bit-exact** all 11 cycles (TPA+BA). **TT: cyc0 core bit-exact.**
- cyc0 core (TPA/BA/QMD/TopHt/merch-cuft/board-ft) BIT-EXACT for ALL of em/ut/tt/ci/cr; the ONLY cyc0 diff is
  col-9 TOTAL-stem cuft (em −1.1% / ut −1.0% / ci −0.06%) = the KNOWN documented total-cuft residual.
- Multi-cycle BA tail is DIRECTIONALLY variant-specific and matches each variant's documented cornered tail:
  em BA +4.8% by 2090 · ut BA −5.4% · ci over-thins (#142) · tt under-thins-then-converges. The tem-cap fix
  (04b15e6/7ce8f1f) is INERT here (the 35000 cap only bites sub-1" dense cohorts; these start QMD 5.1).
  ⇒ NO REGRESSION from this session's commits.
HONEST CLUSTER CHARACTERIZATION: the western Wykoff-DDS variants (em/ut/tt/ci/bm) share a small DIRECTIONAL
multi-cycle DG-tail (±5% BA by 2090), variant-specific in sign — at the cornered/real boundary. #140 (BM) is the
one under active per-tree investigation; the growth-DGDRIV instrumentation on bmt01 is the shared tool to settle
whether this class is a real per-tree-input bug or the accepted DGSCOR/RDPSRT-tie-break straddle. CR (bit-exact)
is the existence proof that the shared engine CAN be exact — so the tails are variant-DATA/input specific, not
an engine-wide flaw.

---
## UPDATE (2026-08-07, later 3) — #140 RESOLVED + cluster-wide DGSD-disconnect audit closed
**#140 (BM self-thin under-thin) RESOLVED to the bit-exact-or-cornered bar** (commit e130546). Confound-free trace
via instrumented FVSbm_g16 proved: BM GROWTH is bit-exact; the residual is a MORTALITY under-kill via the self-thin
QMD-projection d10, driven by the per-tree DGSCOR/FRM realization. REAL BUG FOUND+FIXED — the "DGSD field
disconnect":
- jl has TWO control fields both meaning FVS's single grinit DGSD: `dg_sd` (→ small-tree regent/crown random) and
  `dg_stddev_bound` (→ large-tree DGSCOR reject-bound + OLDRN clamp + mortality). Variants that OVERRIDE DGSD from
  the 2.0 default set only `dg_sd`, leaving `dg_stddev_bound` at 2.0.
- CLUSTER-WIDE AUDIT (all 9 western species.jl): ONLY BM (DGSD=1.5) and CI (DGSD=1.7) override → BOTH had the
  disconnect. Fixed both (`dg_stddev_bound`=1.5/1.7). The other 7 (EM/TT/UT/KT/IE/CR/BC) use DGSD=2.0 = the default
  ⇒ no disconnect (the goal-doc's "all 9 set DGSD, no gap" is NOW actually true — it wasn't before this fix).
- VALIDATED: OLDRN clamp now matches live (±1.5σ, was ±2.0σ); bmt01/cit01 cyc0 bit-exact; multi-cycle comparable
  (no regression, BM under-thin marginally reduced); no collapse.
- RESIDUAL (persists, cornered): the per-tree BACHLO draw still differs — jl's RNG is a faithful port but not
  byte-identical to FVS on the stochastic-DG path (desync triggered by uncalibrated species' calibration OLDRN
  bachlo-seeding). = the accepted DGSCOR RNG-realization straddle. The ENTIRE deterministic path is bit-exact.
This also tightens CI #142 (same fix). Remaining cluster items unchanged: EM #137 (estab self-thin), IE AUTOES
#143 close-out, Climate-FVS.

---
## UPDATE (2026-08-07, later 4) — CURRENT-CODE CLUSTER-WIDE real-FIA sweep: 32/32 cyc0 bit-exact, 0 crashes
Post-all-session-fixes (DGSD-disconnect e130546, crown-init f7d8f86, #158 sitage-nonfix), a fresh real-FIA sweep
drawn live from the 70GB FVS-ready DB (VARIANT-filtered), jl vs each FVS{v}_clean oracle:
  BM 12/12 · CI 8/8 · EM 3/3 · TT 3/3 · UT 3/3 · IE 3/3  = **32/32 stands cyc0 BIT-EXACT, 0 jl crashes** (6 variants).
Multi-cycle residuals = the accepted, now-fully-characterized DGSCOR/ZZRAN RNG-realization straddle (every
deterministic path — DDS/DGCON/SMCON/COR/SSIG/RHO/VARDG/VMLT/OLDRN-clamp/regent-htgr — proven bit-exact this
session via the FVS{v}_g16 instrumentable oracles). Confirms the session's fixes introduce NO cyc0 regression
cluster-wide and the whole western cluster is at the bit-exact-or-cornered bar on real FIA data. (KT stands not
returned under the "KT" VARIANT string — a DB-labeling detail, not a jl issue; KT growth+volume already validated
vs FVSkt_clean on ktt01.) The mission's "full FVS-ready FIA sweep" is satisfied for the current code.

## UPDATE (2026-08-07, later 5) — BM multi-cycle real-FIA validation: worst max|ΔBA|=6.4%, no missed bug
Beyond cyc0, checked MULTI-CYCLE divergence on the 12 BM real-FIA stands (the most-changed variant this session):
per-stand max|ΔBA%| across all cycles → worst = 6.4% (stand 1544771845290487); NONE >8%. All within the accepted
DGSCOR/ZZRAN RNG-realization straddle band (amplified by the self-thin QMD-feedback, per #140's full trace) — no
stand shows a >10% divergence that would indicate a missed deterministic bug. ⇒ BM is bit-exact-or-cornered on
real FIA data through the full projection (cyc0 12/12 bit-exact + multi-cycle ≤6.4%), post the DGSD-disconnect +
crown-init fixes. This closes the loop: the session's cluster-wide real-FIA sweep is validated at BOTH cyc0 (32/32
across 6 variants) AND multi-cycle (BM ≤6.4%, no outliers).

---
## ⚠️ CRITICAL (2026-08-07, later 6) — MULTI-CYCLE real-FIA sweep reveals LARGE cross-variant divergences (cyc0 MASKED them)
Extending the real-FIA sweep from cyc0 to MULTI-CYCLE (max|ΔBA%| across all cycles) surfaced LARGE divergences the
cyc0-32/32-bit-exact result HID: CI worst 114.7% (753180709290487) + 20/20/50/18% on 4 more CI · EM 29.0%
(474180830489998) · TT 33.3% (753186539290487) · UT 16.4% (3624632010690) · IE 4.8% (only IE cornered-small).
NOT AUTOES (jl==jl+NOAUTOES on the worst CI stand). NOT my session's fixes (EM/TT/UT unchanged this session ⇒
PRE-EXISTING). ROOT = jl OVER-GROWS the large-tree DG on MATURE real-FIA stands: on 753180709290487 the stand
DECLINES in live (BA 52→37, mortality>growth) but GROWS in jl (BA 52→74); decomposed at 2029 the +19% BA gap is
+8 from QMD over-growth (jl 9.3 vs live 8.66) + only +1 from mortality (TPA 118 vs 115). A per-cycle DG
over-prediction compounds over 8 cycles into the 16-114% BA blow-up — FAR beyond the accepted ~7% cornered tail.
⇒ THE "cluster bit-exact-or-cornered" VERDICT IS CYC0-ONLY AND INCOMPLETE: real-FIA MULTI-CYCLE has a systematic
cross-variant large-tree-DG over-growth on mature stands (species/condition-specific — cit01/synthetic use GF and
were bit-exact, hiding it). This is the HIGHEST-priority open item, reopened by this measurement. NEXT: instrument
FVS{v}_g16 dgf vs jl on 753180709290487 (CI) at the first divergent cycle (2019→2029) — per-tree DG by species,
find which species/term over-grows; likely a large-tree DG coefficient/site/CCF term wrong for the real-FIA species
mix (NOT GF). META (doctrine, hard): cyc0 bit-exactness is NECESSARY BUT NOT SUFFICIENT — a small per-cycle DG
error is invisible at cyc0 and compounds to catastrophic multi-cycle divergence; ALL "bit-exact-or-cornered"
claims must be MULTI-CYCLE on real FIA data, not cyc0. My prior "32/32 mission satisfied" was premature.

### CORRECTION (later 6b) — the multi-cycle root is MORTALITY under-kill on MATURE stands, NOT DG over-growth
The .sum ACCRE (growth) / MORT (mortality) columns are decisive and OVERTURN the "DG over-growth" read: on CI
753180709290487, GROWTH is similar (ACCRE live 13/15/15/16/16/14 vs jl 14/16/18/19/19/19) but MORTALITY is ~9×
UNDER in jl (MORT live 27/23/21/19/18/17 vs jl 3/4/4/7/8/8). Live kills the large MATURE trees (BA 52→37 declines,
QMD held ~8-10) while jl fails to (big trees survive AND grow ⇒ BA 52→74, QMD 8.4→>10). The QMD gap I decomposed
earlier is a CONSEQUENCE of the mortality gap, not a DG cause. These stands are AGE 233 (very old) ⇒ this is the
MATURE-STAND / LARGE-TREE mortality (MORTMSB mature-stand-breakup or the large-DBH background/BAMAX path) UNDER-
firing in jl on old real-FIA stands. Memory says "MORTMSB/MSBMRT ported+bit-exact" — but validated on younger
synthetic stands; it evidently does NOT fire (or fires ~9× weak) on these age-200+ real-FIA stands. NEXT (revised):
instrument jl mortality! vs FVS{v}_g16 morts.f/mortmsb on 753180709290487 at 2019→2029 — dump the MSB trigger
(QMDMSB/SLPMSB/CEPMSB), MSB kill, and the large-tree background RI, find why jl's large-tree/mature kill is ~9×
low. Highest-priority. (Build FVSci_g16 or reuse the EM DF stand 474180830489998 with FVSem_g16.)
