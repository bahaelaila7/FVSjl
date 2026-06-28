# Bandaid audit — FFE carbon + biomass (Jenkins) + consumption

Module: FFE carbon / biomass / consumption / crown biomass
Files audited:
- `src/engine/fire/carbon.jl`
- `src/engine/fire/biomass.jl`
- `src/engine/fire/consumption.jl`
- `src/engine/fire/crown_biomass.jl`

FVS sources checked: `bin/FVSsn_buildDir/fmcbio.f`, `fmcrbout.f`, `fmdout.f`, `fmcons.f`,
`fmcrowe.f`, `fmsvol.f` (FMSVL2 entry), `fmscut.f`, `fmchrvout.f`, `fminit.f`, `fmsdit.f`,
`fmmain.f`, `fmcrow.f`, `fmvinit.f`, `METRIC.F77`, `FMPARM.F77`.

Overall: `biomass.jl` and `crown_biomass.jl` are faithful, bit-for-bit transcriptions of FMCBIO
and FMCROWE (every coefficient, the <2.5cm/<1in scaling, the SG=V2T/2000 rescale, the small-tree
`MAX(X,MCF)` bole, the UMBTW cones, and the per-form size-class assembly were verified against
source). The Jenkins live-carbon pools, FAPROP harvested-wood-product fate distribution, down-wood
volume/cover coefficients, and the 0.5 / 0.37 carbon fractions are all faithful. Five concerns below.

---

## FLAG 1 — BANDAID: inventory crown snapshot creates cycle-1 crown-lift that FVS gates off

- jl symbol: `write_carbon_report` calling `snapshot_ffe_oldcrown!(stand)` at inventory,
  `carbon.jl:388` (mirrored in the main path at `io/summary.jl:101`).
- Claim/comment (carbon.jl:385-388): "Snapshot the INVENTORY crown as the OLD state ... so the FIRST
  cycle's crown-lift has a valid OLDCRW. Without this the 1st cycle's crown-lift is skipped
  (ffe_oldht=0), losing ~1.9 t/ac of fine down-wood — the DDW sizes-1-3 gap."
- FVS source checked: `fmsdit.f:93` `IF (ICYC.GT.1) THEN` wraps the ENTIRE OLDCRW crown-lift-fall
  computation (the `X = ((NEWBOT-OLDBOT)/OLDCRL)/CYCLEN` scaling that FMCADD later adds to the CWD
  pools, `fmsdit.f:103-119`). On the first growth cycle (ICYC=1) the fall material is **exactly
  zero** regardless of whether OLDCRW was populated. FMMAIN does call FMOLDC at the end of the
  inventory pass (`fmmain.f:268`), so the comment's premise ("FVS calls FMOLDC in the inventory
  FMMAIN") is literally true — but it is *irrelevant*, because the ICYC>1 gate suppresses the
  cycle-1 fall anyway.
- Self-contradiction: the same function's later comment (carbon.jl:402-404) asserts the crown-lift
  is "no-op on the first grow (ffe_old* unset ⇒ zero), matching FVS's ICYC>1 gate." Line 388
  populates `ffe_old*` at inventory, so that premise is false and the first grown cycle's
  `compute_crown_lift!` (carbon.jl:406) returns ~1.9 t/ac of down-wood that FVS does not add.
- Severity: BANDAID. Justification is pure output-matching ("the DDW sizes-1-3 gap"); FVS mandates
  cycle-1 crown-lift = 0. This is the same pattern the audit brief cites as a confirmed bandaid.
- Faithfulness impact: first grown cycle's DDW (and total stand carbon) overstated by the crown-lift
  fall FVS gates off; compounds forward as that down-wood decays through later cycles.

---

## FLAG 2 — GAP: Stand Carbon total omits belowground-dead (root) carbon, which FVS includes by default

- jl symbol: `stand_carbon_report`, `carbon.jl:321-323`
  (`total = above + below + sd + dw + ff + sh`); comment at 298-301 cites "fmcrbout.f:178".
- FVS source checked: `fmcrbout.f:178-183`:
  ```
  V(9) = V(1)+V(3)+V(5)+V(6)+V(7)+V(8)
  IF (LDCAY) THEN
     V(9) = V(9) + V(4)        ! V(4) = BELOWGROUND DEAD (decaying roots)
  ELSE
     V(4) = -1.0
  ENDIF
  ```
  `LDCAY` is set from `CRDCAY > 0` (`fmcrbout.f:60-64`), and `fminit.f:918` sets the SN default
  `CRDCAY = 0.0425` (>0). So by default `LDCAY` is TRUE and the reported total **includes** the
  dead-root pool V(4). The jl computes `belowground_dead` (carbon.jl:316-317, 321) as a column but
  never adds it to `total`. The cited line (:178) is real, but the implementation stops there and
  ignores the very next conditional (:179-180), which fires under the default config.
- Severity: GAP. Faithful only when `CRDCAY = 0` (root decay disabled); under the SN default it
  understates Total Stand Carbon by the dead-root pool whenever any tree has died (`BIOROOT > 0`).
- Faithfulness impact: Total Stand Carbon column low by `belowground_dead_carbon` after the first
  mortality/harvest; grows as `bioroot` accumulates.

---

## FLAG 3 — GAP: fire carbon release applies 0.5 to consumed litter/duff (FVS uses 0.37) and omits live/crown consumption

- jl symbol: `apply_fire_consumption!`, `consumption.jl:38-51` — returns `released * 0.5f0` where
  `released` sums consumed biomass over all 11 surface classes, **including class 10 (litter) and
  class 11 (duff)**.
- FVS source checked: `fmdout.f:286-287` `BIOCON(1) = BURNED(3,10) + BURNED(3,11)` (consumed litter +
  duff), `BIOCON(2) = TOTCON - BIOCON(1)`; then `fmcrbout.f:151`
  `V(11) = BIOCON(1) * 0.37 + BIOCON(2) * 0.50`. The consumed forest-floor (litter+duff) carbon is
  released at **0.37**, not 0.5 (consistent with the forest-floor pool fraction). The jl releases it
  at 0.5, overstating the litter/duff release by (0.50−0.37)/0.37 ≈ 35%.
- Additionally: FVS `TOTCON` includes `BURNLV(1)+BURNLV(2)+BURNCR` (consumed live herb/shrub and
  crown; `fmdout.f:269`), so V(11) carries that burned carbon at 0.5. `apply_fire_consumption!` only
  touches `fs.cwd` surface pools — no `flive`, no crown — so the release also omits live/crown
  consumption. (The file docstring scopes itself to "natural-unpiled-fuels," so the live/crown
  omission is partly disclosed; the 0.37-vs-0.5 litter/duff fraction is not.)
- Severity: GAP. The "Carbon Released from Fire" value (used in the BurnReport via `fmburn.jl:100`)
  is biased high on the forest-floor component and low on the live/crown component.
- Faithfulness impact: incorrect per-fire carbon-release figure; magnitude depends on duff/litter
  load vs live fuel load.

---

## FLAG 4 — GAP: consumption fractions drop the PSBURN scaling and the empty-1-3"-pool special case

- jl symbol: `fire_consumption_fractions`, `consumption.jl:23-29`; applied in
  `apply_fire_consumption!` and `fmburn.jl:100` with no PSBURN argument.
- FVS source checked:
  - `fmcons.f:204-208` scales every `PRBURN(1,I)` (and the live `PLVBRN`) by `PSBURN/100`, where
    PSBURN is the percentage of the stand actually burned. The jl applies the raw fractions, i.e.
    assumes PSBURN = 100 (whole stand burns). For a partial burn the jl over-consumes the pools and
    over-releases carbon.
  - `fmcons.f:147-160`: the <1" classes (PRBURN(1,1), PRBURN(1,2)) are 0.9 only when
    `BURNZ(1,3) > 0`; when the 1-3" pool is empty they burn 1.0. The jl hardcodes 0.9 for classes
    1-2 unconditionally.
- Severity: GAP (edge/partial-burn cases). Faithful for a full-stand burn (PSBURN=100) with a
  non-empty 1-3" pool; otherwise the consumed loadings and the carbon release diverge.
- Faithfulness impact: over-consumption / over-release on partial-area burns; minor over/under on
  the <1" classes when the 1-3" pool is empty.

---

## FLAG 5 — GAP: FFE-method live carbon uses gross cubic volume where FVS uses merch (MAX(X,MCF))

- jl symbols: `ffe_live_carbon`, `carbon.jl:132` (`stem = _fm_cuft(s, sp, d, h) * v2t[sp]`, default
  `merch=false` ⇒ returns `v[1]` gross TCF) and the identical `ffe_fuel_loadings` stem,
  `carbon.jl:101`.
- FVS source checked: `fmdout.f:243-258` (FFE live `TOTLIV`) and `fmcrbout.f:123-130` (FFE V(2) merch)
  both call `FMSVL2(JS,D,H,X=-1,VT,...,LMERCH=.FALSE.,...)`. The FMSVL2 entry in `fmsvol.f:148-155`
  for SN returns:
  ```
  VOL2HT = MAX(X,MCF)        ! X = 0.005454154*H ; MCF = MERCH cubic
  IF (LMERCH) VOL2HT = SCF
  ```
  i.e. with `LMERCH=.FALSE.` for SN the stem volume is the **merch** cubic `MAX(X,MCF)`, NOT the
  gross TCF. The same routine's own small-tree-bole branch in `crown_biomass` (carbon's sibling
  `crown_biomass.jl:118`) correctly uses `merch=true` = `MAX(X,MCF)`, so the gross choice in
  `ffe_live_carbon`/`ffe_fuel_loadings` is an inconsistency, not a different model.
- Severity: GAP. This affects only the non-default FFE carbon method (`carbon_method == 0`; default
  is Jenkins = 1, `state.jl:326`) and the FVS_Fuels DBS standing-live column — paths the file itself
  marks "not yet end-to-end validated" (carbon.jl:120-122). For trees above merch size, gross > merch,
  so the FFE-method aboveground/merch live carbon and the live fuel loadings are biased high by the
  unmerchantable top + stump volume.
- Faithfulness impact: FFE-method (CARBCALC=0) Aboveground/Merch live carbon and the standing-live
  fuel loadings overstated; Jenkins default path unaffected.

---

## Verified faithful (counted, not individually flagged)

- `jenkins_biomass` (biomass.jl): all 10 aboveground coeff pairs, the softwood/hardwood merch+root
  coeffs, `KGtoTI = TMtoTI/1000 = 1.102311/1000`, `INtoCM=2.54`, the <2.5cm 2.5-cm-scaling branch,
  and the `D >= DBHMIN` merch gate all match `fmcbio.f:50-109` + `METRIC.F77`.
- `crown_biomass` (crown_biomass.jl): SG = V2T (post `/2000` rescale, `fmvinit.f:1094`, passed as
  `fmcrow.f:183`), the species totabv-coef selection, foliage/bark/wood fractions, the small-tree
  linear scaling, the small-tree bole `MAX(X,MCF)` via FMSVL2, the DOBF/HTF/LILPCE geometry, the
  UMBTW cone/cylinder weights, all clamps, and the four per-form XV assemblies match `fmcrowe.f`.
- Down-wood volume (`ffe_down_wood`): `CWDDEN` 18.72/24.96 and the cover power-law coefficients match
  `fmdout.f:311-376`; `BIODDW = sum(cwd[1:9])` matches `SMALL2+LARGE2` (`fmdout.f:281`).
- Carbon fractions: 0.5 for woody/live/snag/DDW and 0.37 for forest floor match `fmcrbout.f:155-161`;
  `_TONAC_TO_MTHA` matches `TItoTM/ACRtoHA`.
- HWP: `accrue_harvest_carbon!` group (BIOGRP>5), product (DBH>CDBRK), CDBRK 9/11 (`fminit.f:919-920`)
  match `fmscut.f:147-151`; `harvested_carbon_report` FAPROP fate distribution, emissions residual
  `(1-ΣFAPROP)`, and `stored=V1+V2`, `removed=V3+V4+V5` match `fmchrvout.f:83-104`.
- `fire_consumption_fractions` DIARED/large-class/duff/litter forms match `fmcons.f:168-197`.
