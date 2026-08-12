# ACTIVE GOAL — Western FVS variant cluster + extensions rollout (FVSjl)

## Mission (user-cemented 2026-08-03; re-anchored 2026-08-05)
"Work unattended until ALL extensions for the variants implemented so far are ported and validated. Do a full
FVS-ready FIA sweep as well. Do NOT stop." Port + validate the WESTERN FVS variant cluster and all extensions,
bit-exact-or-cornered vs live FVS oracles, chunk by chunk. Branch: `kt-variant-port`.
DO NOT narrow scope to a single variant — CR is DONE; the goal is the whole cluster + extensions + FIA sweep.

## Variant status — growth+volume ports (oracle = live FVS relinked from bin/FVS{v}_buildDir/*.o)
- **CR** (Central Rockies) ★★ COMPLETE (2026-08-05: 3 bugs fixed — backdated-density dub / forkod imodty /
  strict site-species — DB sweep 1/40→39/40; every residual bit-exact or measured accepted primitive).
- **KT / IE / EM / BM / TT / UT** ★★ growth+volume bit-exact-or-cornered.
- **BC** (British Columbia) — growth+yield complete; remaining: merch/board vol, V2/non-ICH.
- **CI** (Central Idaho) ◐ IN PROGRESS — 9/9 chunks running, cit01 1990 bit-exact; refinement tail OPEN.

## Systematic DG-calibration dispatch audit — COMPLETE (2026-08-05)
Swept every shared COR-shrinkage/bark dispatch for CI-class missing-variant branches:
- **PSIGSQ** (COR Bayes-shrinkage prior variance): was missing CI (fixed 3b9aa35) AND IE (fixed 96cde22, was on the
  SN 0.0898 scalar default). Dispatch NOW COMPLETE: NE/CR/KT/EM/TT/UT/BM/BC/CI/IE all wired. Both fixes source-
  verified faithful, `.sum`-inert on cit01/iet01 (no regression), correct for other-species stands.
- **Bark**: CR/TT/BM/BC/CI have POWER/special bratio wired; EM/IE/UT/KT barks are linear-encodable (c.bark_a/bark_b)
  and the shared bark_ratio clamps [0.80,0.99]==live bratio.f ⇒ their linear fallback is FAITHFUL. CI's POWER
  ci_bratio was the ONLY genuine missing-bark branch. No further bark bugs.
- **DGSD**: all 9 western variants set it explicitly from their grinit.f (2.0 except BM 1.5, CI 1.7). No gap.
⇒ No remaining missing-branch bugs in the shared DG-calibration path cluster-wide.

## Extensions matrix — ALL DONE-OR-CORNERED ✓ (2026-08-12)
- **FFE**: ALL western validated-cornered ✓ (+ eastern + CR).  **Dwarf mistletoe**: ALL western DONE ✓.
  **ECON**: DONE ✓.  **Climate-FVS**: ✓ DONE (~95%, FAITHFUL) — 2026-08-12 line-by-line re-assessment: the
  CLIMDATA reader + clgmult(growth) + clmorts(viability + SPMORT2 transfer-distance DMORT) + clmaxden + clim_autoestb
  are ALL ported, WIRED, and cycle-0 bit-exact vs FVSie_clean; `apply_climate_mort!` matches clmorts.f:205-230 line
  for line. The old "TODO/inert" label was STALE (SPMORT2 was already ported, contradicting a stale in-code comment).
  Residual = a multi-cycle climate-modified self-thin realization (BA cornered, TPA straddle) — same accepted class.
- **FIA sweep**: whole-cluster multi-cycle validated (2026-08-03); residuals = ZZRAN/DGSCOR dense-regen straddle.
  ★ 2026-08-05 POST-FIX real-FIA re-validation (docs/WESTERN_FIA_VALIDATION_2026-08-05.md, stands drawn live from
  the 70GB FVS-ready DB by VARIANT): **CI** 25-stand slice → 6 treed, 0 jl crashes, cyc0 5/6 bit-exact, remaining
  = Δ1-NINT + the accepted DGSCOR/density compounding tail; no regression from the bark/CI_PSIGSQ fixes. **IE**
  60-stand slice → 9+ treed, cyc0 8/9 bit-exact, 0 crashes (monotone partial; no regression from IE_PSIGSQ).
  Both variants bit-exact-or-cornered with ZERO jl crashes on real FIA data. Harness: extract_sample.jl + the
  generalized scratchpad/fia_sweep_check.jl (any variant, reusable cluster-wide).

## ★ 2026-08-12 SESSION UPDATE — #140 RESOLVED, #137 self-thin EXONERATED, Climate-FVS DONE
Two of the four listed "remaining" items are now RESOLVED-or-reframed by end-to-end measurement (via the new
scoped-DEBUG capability that unblocks BM/EM/IE live instrumentation past the fvsvol volume-DEBUG crash — the
DEBUG keyword needs a NON-BLANK field 2 to read a routine onto DBSTK; bare DEBUG=ALLSUB and crashes):
- **#140 BM — RESOLVED, CORNERED.** The "consistent under-thin bias" framing in item 3 below is SUPERSEDED. Full
  chain measured on a dense self-thinner (22960873010497): cyc1 self-thin BIT-EXACT (sdimax/d10/tn10 all match);
  cyc2 divergence traces to the `bm_dubscr` crown-dubbing bachlo N(0,sd) RNG draw (deterministic cr_arg=1.1765 vs
  live 1.1742 = BIT-EXACT; only the random draw byte-differs), which cascades crown→vigor→sub-inch HTGR→breast-
  height crossing→QMD-projection→self-thin (hyper-sensitive). Entire deterministic chain FAITHFUL. = ZZRAN/DGSCOR
  accepted-RNG-primitive class; sign varies by stand (NOT a fixed under-thin bias). docs/…RECONCILIATION_2026-08-12.
- **#137 EM — self-thin EXONERATED (faithful); root reframed to EM sub-inch DG.** On the em_dense reproducer the EM
  self-thin is FAITHFUL (jl tn10=t85d10=29750=live); the divergence is UPSTREAM — jl's EM sub-inch seedlings never
  accumulate DBH (QMD frozen at the 0.3 DIA0 floor while live climbs 0.3→0.7). MEASURED LEAD (2026-08-12): the
  em_dense seedlings carry crown_pct=0 in jl, so `_em_smhtgf`'s beta2·cr height term vanishes ⇒ height crawls
  (h 1.5→2.4 over 2 cycles, never crosses 4.5ft) ⇒ SMDGF never assigns DBH. NEXT: confirm live's crown for these
  seedlings (jl-cr-dub vs live) — if live dubs cr>0 it's a crown-init bug; NOTE em_dense is SYNTHETIC (40000 TPA),
  so the real-FIA EM priority is #143. Self-thin needs NO further work.
- **Climate-FVS — DONE** (see Extensions matrix above; ~95%, faithful, cyc0 bit-exact).

## REMAINING WORK — genuinely open (task-tracker #142/#143/#191/#194 + EM sub-inch DG)
1. **CI refinement tail [#142]**: cit01 jl OVER-KILLS TPA ~2%. RELIABLE STATE (2026-08-05, after FOUR wrong
   root-causes corrected by measurement — backdated-density/GF-COR/bark/deferred-ZZRAN all refuted): CI deterministic
   DG is BIT-EXACT (GF DDS jl==live), serial-corr is ACTIVE (real-run c.sigma[4]=0.26, NOT deferred), COR applied
   (c.dg_cor[4]=0.05693=live). ⇒ the ~2% over-kill is the DGSCOR RNG-realization = the accepted "ZZRAN/DGSCOR
   dense-regen straddle" (cornered; straddles ~0 across stands per the 2026-08-03 FIA-sweep memo) → MEETS the bar.
   ★ HARD LESSON: `each_stand` returns PRE-calibration state (sigma/cor=0) — measure calibration-dependent quantities
   in the REAL run only. 2 real adjacent bugs FIXED (faithful, .sum-inert cit01): 0fa9677 bark branch, 3b9aa35
   CI_PSIGSQ branch. STILL open (unmeasured): volume MATW/FW2W · SMHTGF small-tree stochastic. Oracle FVSci_clean.
   ★ 2026-08-05 SETTLED: the EM/IE growth-only ~7%-BA-by-2090 compounding OVER-GROWTH tail is CORNERED, not a bug.
   Full-precision cyc0-DG test (live EM D@ICYC=2 vs jl exact d2000, NOTRIPLE): per-tree DG diffs are real ~0.5-0.8%
   (large-tree) but MIXED-SIGN and mostly-cancelling (aggregate BA bit-exact) = the accepted RDPSRT/AVHT40 BA-
   percentile/crown-ratio tie-break precision compounding. No fix warranted. (docs/EM_VARIANT_PORT_AUDIT.md)
2. **EM/IE AUTOES establishment [#143]** — the genuine real-FIA establishment priority. AUTOES over-establishes on
   bare/establishment FIA stands (IE fixed earlier — d089b78, jl 253→0.1=live; EM multi-cycle AUTOES still open,
   jl+NOAUTOES bit-exact w/ live). Validate vs stand4_booktpa on UNMODIFIED FVSie_clean. See fvsjl-em/ie memories.
3. **EM sub-inch small-tree DG** (#137 follow-on) — jl's EM sub-inch seedlings under-grow (QMD frozen) on dense
   cohorts; measured lead = crown_pct=0 killing `_em_smhtgf`'s beta2·cr term (regent.jl:340-345). Confirm vs live
   crown-dub; the MIRROR of BM #140 (BM over-grew sub-inch; EM under-grows). Real-FIA impact still to be scoped.
4. **TT aspen bug PAIR [#191]** — sub-1" regent subcycle over-growth + DGFASP-RMSQD under-growth (entangled; the
   #158 single-step-suppressed-model port, entangled with DGFASP crown-init). Land together.
5. **CI bare-establishment regen [#194]** — ci_esgent over-establishment (+282%), same class as the fixed #193.
6. **CI volume tail** — MATW/FW2W merch volume + SMHTGF small-tree stochastic (unmeasured), on cit01. Oracle FVSci_clean.
   (#140 BM, #137 EM self-thin, and Climate-FVS are RESOLVED/exonerated/done — see the SESSION UPDATE above.)

## DOCTRINE (hard-won — carry from the FIA campaign)
1. Validate vs LIVE FVS oracle, bit-exact per chunk. 2. MEASURE, don't infer — instrument the Fortran.
3. Per-record treelist INVALID after tripling — use .sum aggregates / pre-split window. 4. Port faithfully then
validate; a regression on a faithful chunk = examine the oracle. 5. Reuse the shared engine — only add
variant-specific equations + data. 6. Document every chunk verdict in docs/{VARIANT}_VARIANT_PORT_AUDIT.md.

## Off-switch
`touch docs/WESTERN_ROLLOUT_COMPLETE` (USER's call). Per-variant done-flags: docs/{V}_VARIANT_PORT_COMPLETE.
Charters: docs/EXTENSIONS_ROLLOUT_PLAN.md. Memory: fvsjl-ci-variant-port, fvsjl-extensions-rollout,
fvsjl-{em,bm,ie,ut,tt,kt,bc}-variant-port. CR sub-goal retired → docs/CR_VARIANT_PORT_COMPLETE.
