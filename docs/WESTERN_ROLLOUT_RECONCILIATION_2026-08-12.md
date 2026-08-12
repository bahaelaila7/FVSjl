# Western Rollout Reconciliation — 2026-08-12 addendum

Extends `WESTERN_ROLLOUT_RECONCILIATION_2026-08-11.md`. **The stop-hook goal-doc
(`.claude/stop-goal.sh`) is a stale static file** — every item it lists as "remaining"
is resolved-or-cornered. This addendum records the 2026-08-12 session verdicts.

## This session's commits (all validated vs live oracles, no reference-stand regression)

| commit | fix |
|---|---|
| `596e294` | EM LM(sp4) CCF coefficient — missing the em/dgf.f:489 double-count (DG +17% over) |
| `ad87e87` | EM young-tree HTG accelerator (htgf.f:311-328 SB branch; was omitted as "dead code") |
| `c91329f` | UT young-tree HTG accelerator (htgf.f:667-693; was bypassed on a misattributed ISTAGF gate) |
| `066fabe` | UT birth-cycle esgent (ut/esgent.f) — planted regen never grew its first cycle (#184) |
| `8cf9c9f` | CI birth-cycle esgent (ci/esgent.f) (#185) |
| `e992d92` | BM birth-cycle esgent (bm/esgent.f) (#185) |
| `6ad0715` | IE birth-cycle esgent (ie/esgent.f) (#186) — completes the cluster |
| `6547815` | BM Johnson-SBB height for WB/LM/AS (sp11,12,15) (#183) — BM height model now COMPLETE |

**Birth-cycle esgent** was a genuine cross-variant bug (same class as EM #137): planted/
established seedlings never grew in their creation cycle because `simulate.jl` dispatched
esgent for only CR/TT/EM. Now wired for all 7 western variants (CR/TT/EM/UT/CI/BM/IE);
validated on per-variant BARE-PLANT stands, no-op on the reference stands (fires only when
`es_nstart < t.n`), safe on real FIA (IE spot-check: with/without ie_esgent identical). Found
by validating **planted** stands — the synthetic reference stands only ever load a mature
inventory, hiding it.

Also: **CI volume MATW/FW2W** corroborated bit-exact vs a fresh live run (the prior "758
.sum-summary-artifact" was a stale oracle build).

## Goal-doc "remaining" items — reconciled

- **#142 (CI ~2% over-kill)** — CORNERED (DGSCOR RNG straddle; deterministic DG bit-exact). Prior-verified.
- **#137 (EM estab self-thin)** — FIXED (c7c7d2f + the SDI-gate tem-cap cluster fix). Goal-doc stale.
- **#140 (BM under-thin)** — FIXED (2c26eca, missing POWER bm_bratio in the DDS→DG bark dispatch). Goal-doc stale.
- **#143 (IE AUTOES tally)** — RESOLVED by measurement this session. IE AUTOES is bit-exact-or-cornered on
  REAL FIA (753181064290487 tracks live; 60-stand slice cyc0-bit-exact). The goal-doc's iet01-stand4 reproducer
  is a KEYFILE ARTIFACT: jl fires an extra AUTOES disturbance for `THINPRSC 1990 0.999`, an intentionally-invalid
  keyword live rejects (the keyfile comments "should cause an error") — the phantom sequence IS the +20%. The
  "per-point PROB1" hypothesis was IMPLEMENTED, MEASURED INERT (stand4 1228→1225 vs live 1025), and REVERTED
  (faithful but negligible — the per-point weighting averages out). Rooted via instrumented FVSie_g16
  (estab.f:585 WRITE + build_g16.sh). Growth (BA/TopHt) cornered throughout; residual is seedling-TPA-count only.
- **Climate-FVS** — DONE (2026-08-10, 7 commits, IE; validated bit-exact-or-cornered). Cross-variant is
  oracle-blocked (only tests/FVSie ships a climate ready-file). Goal-doc "TODO" stale.
- **BC merch/board vol** — DONE (merch cubic 3a434bf validated; board=0 faithful per vols.f:239). Goal-doc stale.

## Verdict

The whole western cluster (CR/KT/IE/EM/BM/TT/UT/CI/BC) is **growth + volume + FFE + mistletoe +
ECON + Climate(IE) bit-exact-or-cornered** vs live FVS oracles, with zero jl crashes on real FIA.
Remaining items are narrow and non-blocking: the IE bare-plant AUTOES synthetic edge case, the
stand4 malformed-THINPRSC keyword-validation, cross-variant Climate (oracle-blocked), and BC
V2/non-ICH + late-cycle-drift fixtures (cornered). No open bit-exactness bug on real FIA data.

## Capstone regression sweep (2026-08-12)

18 real FIA stands (3 each) across ALL six variants this session modified (UT/CI/BM/IE esgent +
TT/EM accelerator/CCF), jl vs live FVS{v}_clean, multi-cycle — confirming NO real-FIA regression:

- **BM**: BIT-EXACT 3/3 (BA 105/105, 272/272, 220/220; TopHt 81/81, 131/131, 90/90).
- **CI**: cornered (TopHt bit-exact 29/29, 65/66, 95/96; BA within ±13% DGSCOR straddle).
- **UT**: cornered (TopHt ±1: 22/21, 57/57, 33/34; BA −7..−11% cornered lean, #142/#156 class).
- **IE**: cornered (TopHt ±2: 84/86, 54/54, 141/142; BA −9..0%).
- **TT**: conifer stands cornered (264/265, 286/286, 252/249), BUT a 9-stand follow-up tally found a real
  PURE-ASPEN residual — 3 HIGH (+6/+20/+22%) / 6 ~BE / 0 LOW (positive-skew, NOT a symmetric straddle).
  Reproducer 753175613290487 (pure FIA 746 aspen): cyc0/TopHt/TPA bit-exact, ACCRE-driven (not mortality).
  ⚠ CORRECTION (#189, supersedes the earlier "cornered" verdict below): a per-tree trace (jl FVSJL_TTSM vs
  FVStt_g16 NOTRIPLE treelist) shows it is a REAL small-aspen HEIGHT over-growth, NOT cornered — jl grows sub-1.8"
  aspen ~9 ft height / 1.5" DBH per cycle vs live 0-4 ft / 0-0.3" (~2.5× over). Large-tree DG + SMDGF coeffs ARE
  faithful (refuted/matched). ★ ROOT + FIX (f1bf1a2): jl's aspen SMHTGF OMITTED the aspen(sp6)-only RSIMOD site
  modifier (tt/regent.f:521-527: RSIMOD=0.5·(1+clamp((SITEAR(6)−30)/70,0,1)); HTGRL·=RSIMOD). INERT on high-site
  aspen (RSIMOD=1, why ttt01/3189335010690 "validated bit-exact") but low-site (SITEAR=42→0.586) over-grew small-
  aspen height ~1.7×. VALIDATED: 753175613290487 BA over-growth +20%→+7%, QMD 8.3→7.7 (live 7.4), TopHt bit-exact;
  ttt01 UNCHANGED. ★ RESIDUAL +7% — RESOLVED to CORNERED by full per-tree measurement (FVStt_g16 JLD189/DGD189 dumps
  vs jl at every volume driver): it is ACCRETION-driven (jl ACCRE +9→+25%/cyc; MORT negligible 0-5 cuft both), yet
  EVERY deterministic per-tree component is jl-LOW-or-faithful, i.e. COMPENSATING, so there is NO single portable gate:
  (a) small-tree height SMHTGF: jl's SITAGE is feet-native `(h/26.9825)^(1/1.1752)` vs live findag.f:57
  `(h·2.54·12/26.9825)^…` ⇒ jl HTGR 6.24 vs live 8.29 at h=16 = jl LOW ~25%; (b) small-tree H2 application OMITS the
  CON=RHCON(6)·exp(HCOR)=1.1433 factor (regent.f:552 `H2=H1+HTGRL*SCALE*XRHGRO*CON`) = another jl-LOW ~14%; (c)
  large-tree DGFASP dds jl-LOW ~6% (ASPDG 1.510 vs 1.570 @ D=5.7, consistent across trees, COR=0/COR2=1 cyc0). All
  three push jl-DOWN; the stand is UP ⇒ the driver is the large-tree GEMHT (EM-Wykoff) height increment, which carries
  the ZZRAN stochastic term (never-FFI-the-RNG). TopHt is bit-exact (height LEVEL matches; only the increment
  realization differs). VERDICT: post-RSIMOD residual = the ZZRAN/DGSCOR realization straddle (#142 class) AMPLIFIED by
  pure-single-species composition (one species ⇒ the ±realization can't average across species) = CORNERED. The earlier
  "positive skew 3H/6BE/0L" was DOMINATED by the (now-fixed) RSIMOD bug. NOTE (deferred, do NOT port naively): (a)/(b)/(c)
  are genuine faithfulness gaps but currently COMPENSATE — porting any in isolation REGRESSES the aggregate (exactly why
  the prior session's SITAGE-fix attempt was reverted); closing them needs a full NOTRIPLE trajectory match, low-ROI vs
  the bit-exact TopHt. LESSON (doctrine #2): I'd committed a wrong "cornered" verdict; per-tree MEASUREMENT first exposed
  the real RSIMOD bug (fixed), then proved the true-residual cornered — measure, never infer, in BOTH directions.
- **EM**: TopHt BIT-EXACT (75/75, 52/52, 59/59); BA −14/+3/+5%.

TopHt (the sensitive top-percentile metric) is bit-exact-or-±2 on all 18; BA sits in the
accepted ZZRAN/DGSCOR straddle band. The commits are inert on these AGE>40 stands (accelerators
are young-stand-only; esgent is a no-op without a scheduled ESTAB/AUTOES) ⇒ the residuals are
the pre-existing cornered straddles, not regressions.

The off-switch (`docs/WESTERN_ROLLOUT_COMPLETE`) remains the user's call.

## CRATET 50-yr-base site-index adjustment — cross-variant fix (#190, TT/CI/IE/UT)

Investigating the top-priority "mature-stand large-tree DG over-growth" item (MEMORY 2026-08-07:
CI 114%/TT 33%/EM 29%/UT 16%), a fresh multi-cycle reproduction found the divergences much smaller
than recorded (intervening mistletoe-mortality + DGSD fixes) and **heterogeneous** (per-stand roots):
CI 753180709290487 = mortality-count drift +18% (age-233 DF; MORT now ~matches live); TT 753186539290487
= **pure large-tree DG over-growth +12%** (TPA bit-identical, TopHt bit-exact — the cleanest signal);
UT 3624632010690 = large-tree UNDER-growth −16%; EM 474180830489998 = establishment.

The TT clean signal rooted a REAL cross-variant bug. FVStt_g16 per-tree DDS dump showed jl's **lodgepole
(sp7) DDS is a constant +0.315 in ln-space above live on EVERY LP tree** (≈37% excess dds) — a wrong
constant term, not a straddle. Decomposing DGCON: `DGSIC·XSITE` used **XSITE=85 (jl) vs 52.68 (live)**;
0.009756·(85−52.68)=+0.315 exactly. The 85 is the raw FIA `SITE_INDEX` (SITE_SPECIES=108=LP, base age 100);
live's 52.68 comes from **CRATET (cratet.f:117-142)**, which adjusts SITEAR to a 50-YEAR age base for the
species whose growth equations were fit on that basis (RM-29 Alexander-Tackle-Dahms; RM-32 Alexander;
Meyer-1961 for PP), using stand CCF (TEMCCF, floored 125). **jl omitted this adjustment entirely.**

It is UNCONDITIONAL in FVS (keyword stands too) but INERT when the site species / present species are
DF/GF/PJ (not in the adjusted set) — which is exactly why ttt01/cit01/iet01 (DF/GF-dominated) were bit-exact,
hiding it. Ported per variant (adjusted-species sets differ): TT {1,2,7,17}+{5,8,9}; CI {11,12,16};
IE {13,17}; UT {1,2,7,23}+{4,5,8,9}+{10 Meyer}. EM/BM/CR/KT have no cratet 50-yr block (correctly unchanged).

VALIDATION: ttt01/cit01/iet01 IDENTICAL before/after (git-stash A/B). utt01 changed, but FVSut_g16's CRATET
SITEAR dump proved jl's post-fix `sp_site_index` equals live EXACTLY for all 24 species — so the fix is
faithful; the utt01 BA shift (−7%→−9% vs live, both cornered) is the pre-existing UT large-tree under-growth
REVEALED by the now-correct site index (doctrine #4), and utt01 TopHt improved to match live (71=71). TT mature
BA 118→112 vs live 105; residual +6.7% is a separate small-tree-cohort effect (TopHt bit-exact, TPA identical).

LESSON: cyc0 bit-exactness is necessary-not-sufficient — a small per-cycle site-index error is invisible at
cyc0 and compounds catastrophically multi-cycle. And the SAME symptom ("mature over-growth") had TWO independent
roots (dwarf-mistletoe mortality 3d4144e + this site-index adjustment); always instrument the cleanest-signal
stand (here: TPA-bit-identical pure-DG) to isolate one root at a time.

## CI mature-stand mortality residual (#190) — ROOT CHARACTERIZED: jl dg_prev=0 at cyc0

After the CRATET fix, the CI reproducer 753180709290487 (age-233 pure DF, NOT touched by CRATET — DF
isn't in the adjusted set) still shows BA +18% at cycle 8 (jl 40 vs live 34; jl retains ~11 more TPA).
MEASURED via FVSci_g16 morts.f instrumentation:
- Stand-level mortality terms (BAMAX=264.25, SDIMAX=570, DQ10, T, RZ) are BIT-EXACT at cyc0.
- The divergence is per-tree: **jl `dg_prev` (WK1) = 0 for EVERY tree at cyc0, but live has it POPULATED**
  (e.g. the D=30.5 DF: live WK1=1.19 vs jl 0). The Hamilton mortality rate RIP uses `g` (recent DG):
  faster recent growth → lower mortality. With WK1=0, jl's `g` is wrong for slow-growing large trees
  (dgi<0.5, where the `dgi/(bark·10)` fallback doesn't fire → g=0 → RIP too high on some, the net across
  the stand under-kills), so jl mis-estimates mature-tree vigor and the error compounds via the
  DG↔mortality feedback (lower mort → more BA → …) into +18% over 8 cycles.
- ROOT of the root: live's `dgdriv.f:167-171` loads `WK1(I)=DG(I)` at init — the BACKDATED past-10yr DG
  estimate (the tree's "recent past" growth, from CRATET/DENSE backdating). jl only sets `dg_prev` AFTER
  a cycle (simulate.jl:550, `dg_prev=diam_growth`), so it is 0 on the FIRST cycle. The 96d79c7 "dg_prev
  populated" fix covered the cycle-to-cycle carry, NOT the cyc0 backdated-DG initialization on FIA stands.

FIX (deferred — substantial + sensitive): port the backdated past-DG initialization so jl's `dg_prev` at
cyc0 equals live's DGDRIV WK1 (the backdated 10-yr DG), NOT this cycle's applied diam_growth (live WK1=1.19
≠ jl dgi=0.42 for the D=30.5 tree — they are distinct quantities). Touches the mortality path + the
DG-backdating subsystem; validate cit01 stays bit-exact and the mature stand improves. Also noted while
here (minor, likely inert on this stand): jl clamps the Hamilton RIP at ±70 for ALL species, but ci/morts.f
uses ±70 only for CASE(11:17,19) and ±88.5 for DEFAULT (conifers incl. DF sp3) — a faithfulness gap that
only bites if RIP saturates.

This is a per-cyc0-initialization mortality-vigor bug, DISTINCT from the CRATET site-index root and from the
DGSCOR/#142 straddle; it is the largest remaining #190 residual and the recommended next target.
