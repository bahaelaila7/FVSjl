# Faithfulness audit — FFE fire behavior (fuel models + Rothermel + burn + live fuels)

Module: FFE fuel models / Rothermel / burn / live fuels.
Files audited:
- src/engine/fire/fuel_model.jl
- src/engine/fire/rothermel.jl
- src/engine/fire/fmburn.jl
- src/engine/fire/fmcba.jl
- src/engine/fire/fire_effects.jl
- src/engine/fire/fuel_moisture.jl

FVS source checked (all in bin/FVSsn_buildDir/): fmfint.f, fmcfmd.f, fmdyn.f, fmburn.f,
fmcba.f, fmeff.f, fmbrkt.f, fmmois.f, fmgfmv.f, fmcfmd2.f, fminit.f, fmtret.f.

Overall: the Rothermel core, the fuel-moisture/wind tables, the standard-model SAV defaults,
the FMCFMD candidate selection, the FMDYN inverse-distance geometry, the FMEFF mortality
logistics + bark thickness, scorch height, and the FMCBA dead-fuel default split are all
FAITHFUL to the cited FVS source. Five items are flagged below (one BANDAID, four GAPs).

---

## FLAG 1 — BANDAID: weighted flame length is summed per-model instead of recomputed from the weighted Byram intensity

- jl symbol/line: `fmburn!` — fmburn.jl:58-65 (`flame_raw += r.flame * w` then `flame = flame_raw * flmult`);
  same pattern in `potential_fire`/`scenario` — fmburn.jl:299-303 (`flame += r.flame * w`).
- Claim/comment (fmburn.jl:55-56): "integrate Rothermel over them, summing the weighted flame & Byram."
- FVS source checked — fmfint.f:507-541. Inside the per-model loop FVS *does* accumulate a
  weighted flame (`FLAME = FLAME + (0.45*(BYRAMT/60.0)**0.46)*FWT(INB)`, line 508), but after the
  loop it **discards that accumulation and recomputes flame from the weighted Byram**:
  ```
  C     COMPUTE FLAME LENGTH AS A FUNCTION OF BYRAM RATHER THAN USING
  C     THE WEIGHTED AVARAGE AS COMPUTED IN THE BLOCK ABOVE (NLC 21 AUG 2003)
        FLAME = (0.45 * (BYRAM/ 60.0) ** 0.46)            ! fmfint.f:541, BYRAM = Σ BYRAMT·FWT
  ```
  i.e. final `FLAME = 0.45·(Σ wᵢ·byramᵢ / 60)^0.46`, NOT `Σ wᵢ·0.45·(byramᵢ/60)^0.46`.
  fmburn.f:436/439 then takes that flame directly (`OLDFL = FLAME`) for the SN path.
- The jl `byram` accumulation (fmburn.jl:63, used for scorch) is correct, but the jl `flame`
  is the pre-2003 weighted-sum-of-flames that fmfint.f:538-541 explicitly removed. Because
  `x^0.46` is concave, the jl flame is biased LOW relative to FVS whenever >1 model is weighted
  (Jensen). This is the source of the known, *accepted* residual recorded in memory
  ("flame 3.84 / scorch 15.65 vs Fortran 3.9 / 15.9") — accepting an output-mismatched value
  with no source basis is the bandaid signature.
- Faithfulness impact: REAL. The flame value feeds `fire_tree_mortality` for SN mortality
  groups 1–5 (oaks/hickory/red-maple/black-gum: `charht = flame·0.7`, fire_effects.jl:94) and is
  the reported flame length; a low-biased flame under-kills those species and mis-states the
  PotFire/Burn report. Fix: drop `flame_raw`, set `flame = 0.45·(byram/60)^0.46` from the summed
  `byram` (then apply `flmult`), per fmfint.f:541.

---

## FLAG 2 — GAP: scorch height not re-derived from the multiplied Byram when FLAMEADJ flame multiplier ≠ 1

- jl symbol/line: `fmburn!` — fmburn.jl:65-67 (`flame = flame_raw * flmult` … `sch = scorch_height(byram, …)`).
- Claim/comment: `flmult` documented as the FLAMEADJ flame-length multiplier (fmburn.jl:31).
- FVS source checked — fmburn.f:441-472. When the multiplier changes the flame
  (`IF (FLMULT .NE. 1.0) FLAME = OLDFL*FLMULT`, line 443), FVS *inverts* the flame→Byram relation
  before computing scorch: `IF (FLAME .NE. OLDFL) BYRAM = 60.0*((FLAME/0.45)**(1.0/0.46))`
  (lines 459-464), then `SCH` uses that revised BYRAM (lines 470-472).
- jl multiplies `flame` but computes `sch` from the *original* unmultiplied `byram`, so with
  `flmult ≠ 1` the scorch height (and hence crown-scorch mortality) does not track the adjusted
  flame the way FVS does.
- Faithfulness impact: NONE on the default/tested path (FLAMEADJ absent ⇒ `flmult = 1`, both
  agree). Silent divergence only when the FLAMEADJ multiplier keyword is used. Low severity.

---

## FLAG 3 — GAP: per-tree fire-mortality multiplier FMORTMLT is not applied

- jl symbol/line: `fmburn!` — fmburn.jl:82-85 (mortality computed and clamped with no per-record
  multiplier); likewise `potential_fire`/`scenario` fmburn.jl:311-314.
- FVS source checked — fmeff.f:340: `PMORT = PMORT*FMORTMLT(I)` is applied to every tree record
  after the logistic and before the `>1`/`<0` clamps (fmeff.f:342-343).
- jl has no `FMORTMLT(I)` term anywhere in the kill loop.
- Faithfulness impact: NONE on the default/tested path (`FMORTMLT` defaults to 1.0 for every
  record). Diverges only when a keyword sets a non-unit per-tree mortality multiplier. Low severity.

---

## FLAG 4 — GAP: FMDYN's "truncate to 4 models + reweight" final step is omitted

- jl symbol/line: `_fmdyn` — fuel_model.jl:159-293. Output is the up-to-`MXFMOD`(=5) `(model,weight)`
  list assembled at lines 277-291; there is no final cap-to-4 renormalization, and no descending sort.
- FVS source checked — fmdyn.f:374-393. After the collinear-split, FVS sorts by weight (RDPSRT,
  line 364) and then **truncates to 4 models and renormalizes their weights to sum to 1**:
  ```
  XWT = Σ_{I=1,4} FWT(I);  IF (XWT>1e-6) { FWT(1:4)/=XWT; FMOD/FWT(5:MXFMOD)=0 }   ! lines 381-393
  ```
- The jl can emit 5 entries when one collinear bracket (the shared 5,15 litter iso-line carrying
  candidates 5+8+9) coincides with ≥2 other brackets (e.g. iffeft 1–3 with 4 < SMALL ≤ 6, which
  also turns on EQWT(5)); FVS would drop the lowest-weight 5th model and renormalize the other four,
  jl keeps all five un-renormalized. The missing descending sort is harmless (the burn loop only
  uses the weighted sums), but the missing 4-cap/reweight leaves the weights slightly off in that case.
- Faithfulness impact: narrow. Only fires when >4 weighted models survive collinear expansion;
  for the common SN candidate sets (≤4 distinct models) jl and FVS agree exactly. Low severity.

---

## FLAG 5 — GAP: FMCBA dead-fuel "soft" column (FUELSOFT / photo-code) is not representable

- jl symbol/line: `fmcba!` — fmcba.jl:62-77. `stfuel = ffe_dead_fuel_loading(...)` is a 1-D vector,
  and the loading is added only to `fs.cwd[isz, 2, idc]` (the J=2 = hard column).
- Claim/comment (fmcba.jl:60-61): "STFUEL's 'soft' column is 0 here, so only the dead (J=2) pool is
  populated."
- FVS source checked — fmcba.f:284-285 (default `STFUEL(ISZ,2)=FUINI(...); STFUEL(ISZ,1)=0`) and
  the split loop fmcba.f:378-396 which runs **`DO J=1,2`** and adds `STFUEL(ISZ,J)` into
  `CWD(1,ISZ,J,IDC)` for BOTH columns. The soft column `STFUEL(*,1)` is non-zero when the FUELSOFT
  keyword (fmcba.f:354-362) or the photo-code path (`STFUEL(I,1)=FOTOVALS(I)`, fmcba.f:298) is used.
- The comment is correct for the DEFAULT case (soft column = 0, so only the hard J=2 pool loads),
  and jl is faithful there. But `ffe_dead_fuel_loading` returning a single 1-D vector means jl
  structurally cannot carry the soft column, so FUELSOFT / photo-code soft loadings are silently
  dropped.
- Faithfulness impact: NONE on the default/tested path (validated bit-exact for DDW carbon per
  memory). Silent divergence only under FUELSOFT or photo-code keywords. Low severity.

---

### Items reviewed and found FAITHFUL (not individually flagged)

fuel_moisture table (fmmois.f:52-98); fire_wind_reduction ALGSLP over CANCLS/CORFAC
(fmburn.f:390, fmvinit.f:36-44); standard_fuel_model SAV defaults 109/30/1500/1500 and the
dead-herb=live-herb SAV rule (fminit.f:167-170, fmgfmv.f:79); select_fuel_models candidate
logic for iffeft 1–9 incl. the redcedar Σht/Σtpa average (fmcfmd.f:131-206); _fmdyn collinear
EQMOD/rescale + signed-distance neighbor search + inverse-distance weights (fmdyn.f:108-330);
_small_large_fuel size classification (fmtret.f:378-388); Rothermel constants RHOP/TMIN/SILFRE
and the full reaction-intensity / propagating-flux / wind-slope / Byram / flame body
(fmfint.f:91-535); live moisture-of-extinction (fmfint.f:350-358); scorch_height (fmburn.f:470-472);
crown_volume_scorched (fmeff.f:170-186); fire_mortality_group (fmeff.f:207-221); MORTB0/1/2
coefficients (fmeff.f:87-89); the Regelbrugge-Smith and Reinhardt logistics incl. the MNMORT
overflow guard (fmeff.f:188-252); fire_bark_thickness B1 table + shortleaf Harmon quadratic
(fmbrkt.f:21-129); fire_mortality_adjust no-op for SN (fmeff.f:196,278,304 — all LS/ON/NE-gated);
the DBH≤1 & CSV>50 ⇒ 1.0 universal rule (fmeff.f:330); crown-fire share of curkil (fmeff.f:543-544);
fire-killed BA constant 0.005454154; FMCBA cover-type/percent-cover (fmcba.f); build_dynamic_fuel_model
woody/herb loads + depth + moisture-of-extinction (fmcfmd2.f:562-587) and its dead-herb split
(fmgfmv.f:88-97).
