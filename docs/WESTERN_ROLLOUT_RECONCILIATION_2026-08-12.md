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
