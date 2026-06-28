# Density module audit — standstats.jl + crown_width.jl

Scope: `src/engine/standstats.jl` (NOTRE expansion, summary stats, CCF, SDI, PTBAL, PCT)
and `src/engine/crown_width.jl` (CWCALC crown-width library).

FVS sources checked: `base/dense.f`, `base/notre.f`, `base/pctile.f`,
`bin/FVSsn_buildDir/ccfcal.f`, `bin/FVSsn_buildDir/sdical.f`, `bin/FVSsn_buildDir/ptbal.f`,
`bin/FVSsn_buildDir/cwcalc.f`, `bin/FVSsn_buildDir/grinit.f`,
plus the data CSVs `data/southern/crown_width_{equations,species}.csv`.

## Verdict

The module is overwhelmingly FAITHFUL. Every non-trivial constant, ordering choice, threshold,
and cap traces to a specific FVS line. One item could not be fully confirmed within scope and is
flagged UNVERIFIED.

### Verified-faithful highlights (not individually flagged)

- `BA_PER_TREE = 0.005454154` matches `dense.f:190` (BATREE) and `ptbal.f:148`. (The truncated
  `0.0054542` in `sdical.f:123` is a *different* routine — the SDImax BA weighting, not used here.)
- `notre!` expansion (fp/vp/fp2, BAF≤0 branch, BRK split, 9e-25 floor, ×GROSPC) matches
  `notre.f:43-69` line-for-line.
- `stand_top_height` (largest-DBH 40 TPA, last tree prorated) matches `dense.f:285-298`.
- `point_basal_area!` (descending-DBH per-point accumulation, PTBALT = BA-in-larger, scale
  PI/GROSPC, PTBAA = point BA) matches `ptbal.f` Western/Southern default branch (SN is not
  CS/LS/NE/ON, so it takes `CASE DEFAULT` at `ptbal.f:69-154`).
- `stand_pct!` (cumulative D²·TPA from smallest up, ÷total·100) matches `pctile.f:49-71`
  driven by `dense.f:273` `CALL PCTILE(ITRN,IND,WK5,...)` with `WK5 = D·D·P` (`dense.f:187`).
- `stand_ccf`: `0.001803·CW²·P` for D>0.1 else `0.001·P`, CR fixed at 90, IWHO=1 — matches
  `ccfcal.f:54-63` exactly.
- `stand_sdi` Zeide form `Σ TPA·(D/10)^1.605` over D≥dbh_zeide matches `sdical.f:326`; the SN
  default `LZEIDE=.TRUE.`, `DBHZEIDE=0`, `DBHSTAGE=0` are confirmed at `grinit.f:129,262,263`.
- `stand_sdi_reineke` A/B/SDIC Stage form matches `sdical.f:281-283` (and the whole-stand reduction
  `SPROB·A + B·SDSQ` = `Σ(A+B·D²)·TPA`, `sdical.f:327`).
- `hopkins_index` matches `cwcalc.f:91-96` term-for-term.
- `_cw_eval` four families match `cwcalc.f`: bechtold MIND=5 floor + inner `min(D,dbh_cap)`
  (the `IF(D.LT.30/40/50/24/18)` caps) + optional outer `max_cw`; smith `(a+b·Dcm+c·Dcm²)·3.28084`
  with Dcm=D·2.54 and OMIND=3; ek `a+b·D^power` OMIND=3; braggm (eq 76102) `a+b·D^power` with
  MIND=5 floor + cap 52 (`cwcalc.f:1954-1960`). Spot-checked dbh_cap/max_cw rows in the CSV
  (09701/11101/26101/31601/12501/12801/11001/12101/11005/09104) — all match the Fortran caps.
- Final clamp `[0.5, 99.9]` and unknown-species → 0.5 match `cwcalc.f:2456-2460` (CW initialized 0,
  clamped up to 0.5). The 23 pure-`bragg` "02" equations are dead data in both FVS and the CSV
  (never selected by the SN species map) — confirmed none are referenced in
  `crown_width_species.csv`.

decisionsReviewed ≈ 19.

---

## FLAG 1 — `notre!` omits the FINT/FINTM inflation FVS applies to dead records

- **jl symbol / line:** `notre!`, `src/engine/standstats.jl:54-62` (loop over `1:(t.n + t.ndead)`).
- **Claim / comment:** line 54-55 — *"expand live records and the dead partition (n+1:n+ndead)
  **alike** — dead trees carry their expanded TPA into the backdated calibration BA."* The dead
  partition is expanded with the **same** fp/vp/fp2 factors as the live records.
- **FVS source checked:** `base/notre.f:99-125`. After processing the projectable records,
  NOTRE loops the non-projectable / recent-dead records (`I1=IREC2 … I2=MAXTRE`) but FIRST
  multiplies the expansion factors by FINT/FINTM:
  `VP=VP*(FINT/FINTM)`, `FP=FP*(FINT/FINTM)`, `FP2=FP2*(FINT/FINTM)` (`notre.f:122-124`).
  So FVS does **not** treat dead trees "alike" — their stored PROB is inflated by FINT/FINTM.
  FVS later *de-inflates* when consuming them, e.g. `sdical.f:139` `DPROB = PROB(II)/(FINT/FINTM)`.
- **Severity:** UNVERIFIED.
- **Faithfulness impact:** When FINT==FINTM (the common case; FINTM defaults to 5, and many SN
  scenarios have a 5-yr first cycle) the factor is identity and `notre!` is exactly faithful.
  When FINT≠FINTM, jl stores the *true* (de-inflated) dead TPA and uses that value consistently
  downstream (no FINT/FINTM division was found anywhere in `src/` except the GROWTH-keyword
  default), whereas FVS stores the inflated value and de-inflates per consumer. The two
  conventions are equivalent **iff** every FVS de-inflation exactly cancels its inflation for all
  consumers of the dead partition. That cancellation is plausible and the natural-process path is
  reported bit-exact vs the oracle, but I could not confirm it within this module's scope.
  **To verify:** (1) confirm no jl consumer of the `n+1:n+ndead` partition applies a FINT/FINTM
  factor, and (2) confirm a live-FVS run with FINT≠FINTM (e.g. INVYEAR offset giving a non-5-yr
  first cycle plus a recent-mortality tree list) produces the same backdated calibration BA / SDI
  as jl. If both hold, downgrade to faithful.
