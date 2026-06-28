# Bandaid audit — Height growth (large + small tree) + regen height

Module: `src/variants/southern/height_growth.jl`, `src/variants/southern/small_tree_growth.jl`
FVS source checked: `sn/htgf.f`, `sn/htcalc.f`, `sn/regent.f`, `sn/cratet.f`, `base/gradd.f`,
`base/update.f`, `bin/FVSsn_buildDir/blkdat.f`, `bin/FVSsn_buildDir/initre.f` (NOHTDREG),
plus `_htdbh_dbh` in `src/engine/volume.jl`.

## Summary

`height_growth.jl` (HTGF/HTCALC) is **clean** — every coefficient, modifier and clamp traces
to `sn/htgf.f` / `sn/htcalc.f` exactly (yellow-poplar coeffs htcalc.f:150-156; HGMDCR Hoerl
htgf.f:223; the five HGMDRH Chapman-Richards factors htgf.f:235-242; HTGMOD clamp htgf.f:256-257;
the `HTMAX-HTI<=1 ⇒ 0.10·XHT·SCALE·EXP(HTCON)` branch htgf.f:192-193; SIZCAP(,4) cap htgf.f:286-289;
`scale = fint/5` = FINT/YR with YR=5 from blkdat.f:61). No flags.

`small_tree_growth.jl` (REGENT) is faithful on the dominant increment. In particular the
**`scale2 = 1f0` hardcode is NOT a bandaid**: FVS regent.f:135 sets `SCALE2 = YR/FNT (= 5/FINT)`
which *shrinks* the DDS, and GRADD then *re-expands* by `SCALE = FINT/YR` (gradd.f:79-87); the two
cancel to the plain FINT-period increment. FVSjl folds GRADD's expansion into the large-tree DDS
(`dds*(sfint/5f0)`, diameter_growth.jl:622) and correspondingly uses the identity (`scale2=1`,
`sqrt((d·B)²+dg·(2·B·d+dg))−B·d ≡ dg`) for the regen path — net result equals FVS for all FINT.
Verified `base/update.f:115` adds `DG/BRATIO` directly (no second expansion), so there is no
double-count. Three lower-severity gaps remain, all only reachable off the snt01/5-yr-cycle path.

## Flags

### 1. GAP (low) — budwidth DIAM-floor + dg_bound applied at the wrong period scale for FINT≠5
- jl: `small_tree_growth.jl:112-113`
  ```
  (d + dg) < regent_diam[sp] && (dg = regent_diam[sp] - d)
  dg = dg_bound(dlo_v, dhi_v, sp, d, dg, sizcap)
  ```
- FVS: `sn/regent.f:359-370` sequences **scale2-shrink (line 359-360) → DIAM budwidth floor
  (363-365) → DGBND (370)**, i.e. the floor and bound act on the *YR-shrunk* (≈5-yr) DG; GRADD
  (gradd.f:79-87) then re-expands by `FINT/YR`. FVSjl, having folded GRADD in, applies the floor and
  `dg_bound` to the **full-FINT** DG instead.
- Impact: identical at FINT=5 (scale2=1, shrink/expand both =1). For a 10-yr cycle (timeint10/s5/s9)
  the budwidth floor `DIAM(ISPC)` and the DGBND clamp bind at a different effective magnitude than
  FVS (FVS floors the 5-yr DG then ×2 the DDS; FVSjl floors the 10-yr DG once). Sub-inch, and only
  when the floor/bound actually binds (very slow-growing sub-3" regen). Faithful otherwise.

### 2. GAP (low) — `HK ≤ 4.5` micro-DBH bump dropped
- jl: `small_tree_growth.jl:35` (`_regent_dg`): `hk <= 4.5f0 && return 0f0`
- FVS: `sn/regent.f:285-287` does more than zero the increment — it assigns
  `DBH(K) = D + 0.001*HK` (a budwidth-scaled stem bump) in addition to `DG(K)=0.0`. FVSjl returns
  `dg=0` and leaves dbh unchanged, silently dropping the `+0.001·HK` bump.
- Impact: only for records whose post-growth height stays ≤4.5 ft, i.e. dbh≈0 seedlings carried in
  the live list (dbh>0 implies height>4.5 by definition of breast height, so unreachable for normal
  small trees). The dropped term is ≤0.0045 in. Negligible but a real, un-ported source rule. The
  jl comment ("DBH bump path; not hit for snt01") acknowledges it.

### 3. GAP (low) — height–diameter calibration branch (LHTDRG / NOHTDREG) not modeled
- jl: header comment `small_tree_growth.jl:18` ("LHTDRG=false ⇒ the HTDBH-inverse branch"); the
  Wykoff `DKK=(BX/(ALOG(HK-4.5)-AX))-1` form is never used — `_regent_dg` always calls
  `_htdbh_dbh` (the inventory inverse). `state.jl:282` initializes `ht_drag_sp` all-false and nothing
  ever sets it true.
- FVS: `sn/regent.f:294-327` keeps the Wykoff DKK/DK (lines 300-305) when
  `LHTDRG(ISPC) .AND. IABFLG(ISPC).NE.1`; otherwise HTDBH (inventory) overwrites them. `LHTDRG`
  defaults FALSE (`sn/grinit.f:104`) and is only set TRUE by the **NOHTDREG/HTDREG keyword**
  (`bin/FVSsn_buildDir/initre.f:2619/2640/2663`, field-2>0).
- Impact: faithful for the default and for snt01 (LHTDRG false ⇒ inventory inverse, exactly what
  FVSjl does). It is a GAP only if a user supplies the HTDREG calibration keyword for a species that
  also re-fits IABFLG≠1 in CRATET — then FVSjl would use the inventory inverse where FVS uses the
  Wykoff calibrated inverse. The keyword path is unported; a self-consistent scoping decision rather
  than an output-matching hack.

## Faithful (reviewed, not flagged)
HTGF: yellow-poplar special coeffs; PCOM/`eco_unit[1]=='M'` montane test; HTCALC age/incr/height
formulae (htcalc.f:166-192); `htcalc_age` defensive ratio clamp (guarded by `HTMAX-HTI>1`, never
engages on a valid path); RELHT≤1.5 clamp; HGMDCR/HGMDRH; WTCR=0.25 weighting; HTGMOD∈[0.1,2];
`EXP(HTCON)` calibration via HTCONS; uniform RHM/RHXS/RHK scalars (htgf.f:128/133/138); SIZCAP cap on
both the over-max and main branches.
REGENT: `scale=fint/5`; `scale2=1` (GRADD-fold cancellation, verified); `dgmx=5·scale`;
`con=RHCON·exp(HCOR)` with RHCON=RCOR2-or-1 (regent.f:583-585); XWT blend; BACHLO rejection loop +
`±0.1·HTGR` random effect under DGSD≥1; per-record size cap; `_regent_dg` negative-DK fallback
(`htg·0.2·bark`), DG≥0.1 / ≤DGMX clamps; per-stand `bark_ratio`; species sort + tripling stash order.
