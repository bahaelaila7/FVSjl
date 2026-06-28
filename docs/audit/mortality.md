# Bandaid Audit — Southern Mortality (MORTS / VARMRT / SDICHK / SDICAL)

File audited: `src/variants/southern/mortality.jl`
FVS sources checked: `sn/morts.f`, `sn/varmrt.f`, `bin/FVSsn_buildDir/sdical.f`,
`bin/FVSsn_buildDir/sdichk.f`.

## Verdict

The **core spine is faithful**. The Pretzsch self-thinning line solve (`_pretzsch_tn10`
vs morts.f:334–468), the QMD-convergence iteration, the background/density RIP merge
(morts.f:504–525), VARMRT geometric distribution (varmrt.f verbatim incl. PEFF, VARADJ,
MINSTP, npass/short_v logic), the species-sorted SDI-sum order, the linear `(DG/BARK)·(FINT/5)`
trajectory at all four call sites, SIZCAP (morts.f:691–694), BAMAX iteration (morts.f:704–771),
and the line-reset test (morts.f:245) all trace to specific FVS lines and match. The
`t.crown_ratio`-as-PCT field naming (FVS PCT = basal-area percentile, ARRAYS.F77:201) is
intentional and correct (standstats.jl:152). No bandaids found in the *tested* path.

All flags below are **GAPs**: faithful on the snt01/5-yr/SN-no-climate path but silently
divergent from FVS source on a real but untested/unreached case.

---

## GAP 1 — `sdimax < 5` does background-only, but FVS kills the whole stand

- jl: `mortality!`, lines 276–277. `if sdimax < 5f0; _varmrt!(... bg_tokill)`.
- FVS `sn/morts.f:343–346`: `IF(SDIMAX .LT. 5) THEN; TN10=0.; GO TO 271`. Comment (340–342):
  "ASSUME CLIMATE HAS CHANGED ENOUGH THAT THE SITE WILL NO LONGER SUPPORT TREES, AND KILL ALL
  EXISTING TREES." With TN10=0, RN=1 (line 472) and SUMTRE=T−TN10=T (line 558), so VARMRT is
  handed **all TPA**, not the Hamilton background total.
- jl instead passes `bg_tokill` (the per-tree Hamilton sum), so only background mortality is
  applied where FVS removes every tree.
- Impact: currently **unreachable** — `stand_sdimax` never calls CLMAXDEN, so jl `sdimax` is the
  raw BA-weighted SDIDEF (~300–600) and never < 5. Wrong only if the Climate extension is later
  ported. Faithful today; will under-kill catastrophically on a climate-kill stand.

## GAP 2 — TPAMRT (next-cycle line-reset basis) is threshold-filtered; FVS sums all trees

- jl: `mortality!`, lines 369–373. `surv` sums `t.tpa[i]-killed[i]` only for `t.dbh[i] >= dthresh`,
  then `s.density.tpa_mort = surv`.
- FVS `sn/morts.f:710–722` (and the BAMAX re-loop 743–758) computes `TNEW = Σ(PROB−WK2)` over
  `DO 36 I=1,ITRN` with **no `IF(D.LT.DBHSTAGE/DBHZEIDE)` guard**, then `TPAMRT=TNEW` (line 772).
- The reset test next cycle (morts.f:245, `ABS(T−TPAMRT).GT.1`) compares the threshold-filtered
  `T` against this **unfiltered** TPAMRT. jl makes both sides filtered.
- Impact: identical when every tree exceeds threshold (snt01). In a stand carrying sub-threshold
  (< DBHSTAGE) regen, FVS's TPAMRT includes those stems while next-cycle `T` excludes them →
  FVS tends to reset the self-thinning line every cycle; jl will not. Diverges the persisted
  SLPMRT/CEPMRT trajectory on regen/seedling stands.

## GAP 3 — SDICHK uses the 0.3-floored QMD for the decision and reset; FVS uses the unfloored QMD

- jl: `sdi_max_check!`, lines 71–78. A single `dq0` is floored (`dq0 < 0.3f0 && (dq0 = 0.3f0)`)
  and then used for both `temmax` (the over-density test) and `const_v2`/`tem2` (the SDImax reset).
- FVS `sdichk.f`: keeps **two** diameters. `DQ0` is floored to 0.3 but feeds only `TMD0`→`UPLIM`
  (the cosmetic warning). The *decision* uses `TEMMAX = CONST*(TEMD0**-1.605)` and the *reset*
  uses `EXP(ALOG(TEMTPA+1)+1.605*ALOG(TEMD0))`, where `TEMD0 = RMSQD` (or `DR016`) is **never
  floored**.
- Impact: equal when QMD ≥ 0.3 (all normal/snt01 stands). For an over-dense stand with mean
  diameter < 0.3 in (dense sub-inch regen), jl uses 0.3 where FVS uses the true smaller QMD, so
  the reset `tem2` SDImax differs. Wired live via simulate.jl:28.

## GAP 4 — user-set BAMAX (LBAMAX) ignored in `stand_sdimax`/`bamax`

- jl: `stand_sdimax` (35–44) always returns the BA-weighted SDIDEF average, and `mortality!`
  line 346 recomputes `bamax = sdimax*0.5454154*pmsdiu` from it.
- FVS `sdical.f:203–209`: `IF(.NOT.LBAMAX) BAMAX = XMAX*0.5454154*PMSDIU` **ELSE** invert the
  user BAMAX into `XMAX = BAMAX/(0.5454154*PMSDIU)`. When the BAMAX keyword is set, FVS drives the
  whole density limit from the user value, not the BA-weighted SDIDEF.
- Impact: no effect unless the BAMAX keyword is used (not in snt01). If a user supplies BAMAX, jl's
  SDImax (and hence tn10, RN, and the BAMAX cap) diverge.

## GAP 5 — MSB / alternate-mortality (MSBMRT) is entirely unported

- jl: no MSB code anywhere in `mortality!`.
- FVS `sn/morts.f:619–680`: when `SLPMSB ≠ 0` and `D10 > QMDMSB`, FVS computes the mature-stand
  boundary, kills `TMORE = TN − T85MSB` additional TPA via `CALL MSBMRT`, and flags recalibration.
- Impact: no-op unless the alternate-mortality (mature stand boundary) keyword is set. Silent
  under-mortality on stands that schedule it; out of scope for the default SN suite.

## GAP 6 — T>35000 cap applied before vs after the DQ0/D10 computation (negligible)

- jl: `mortality!` line 249 caps `tt > 35000f0 && (tt = 35000f0)` **before** computing `dia0`/`d10`
  (lines 250–251).
- FVS `sn/morts.f`: computes `DQ0/DR0/DQ10/DR10` at lines 253–267 with the **uncapped** T, then
  caps `IF(T .GT. 35000.) T=35000.` at line 349.
- Impact: only matters at > 35000 TPA (absurdly dense; not realistic for SN). Listed for
  completeness — the start-of-cycle mean diameters would use a slightly different denominator at
  that extreme.
