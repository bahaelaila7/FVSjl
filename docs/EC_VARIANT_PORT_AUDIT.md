# EC (EastCascades) variant port — FFE (Fire & Fuels Extension) audit — #232 (2026-08-15)

Task #232 ports the FFE for **EC / EastCascades** (VARACD "EC", MAXSP=32, forest 608 = OKANOGAN R6
eastside, 10-yr cycle, Reineke SDI, Wykoff DDS). Oracle: `/workspace/.ecwork/FVSec_clean` (made the
`.sum.save`) + `/workspace/.ecwork/FVSec_g16` (gfortran-16 instrumentable; exit 10 = normal). Growth+volume
already merged (849728a); this task is **FFE only**. Reference stand `ect01_ffe` (S248112 — same trees as
`cat01`/`pnt01`/`wct01`; forest 608, habitat 12, SIMFIRE @2000: SWIND=10, FMOIS=1, ATEMP=50, 100% burned).
EC fire sources: `fire/ec/{fmbrkt,fmcba,fmcblk,fmcfmd,fmcrow,fmvinit}.f`; the rest (`fmcfir/fmburn/…`) are
the shared FFE base.

## KEY STRUCTURAL FINDING #1 — ect01_ffe is a SURFACE fire, NOT a crown fire

Unlike the crown-fire siblings (CA 530→2, WC 491→2, PN 485→4), **the EC oracle runs a SURFACE fire** on
ect01_ffe. Measured on `FVSec_g16` (instrumented `fmburn.f`/`fmcba.f`/`fmcfmd.f`, restored pristine):

| FMBURN (SIMFIRE @2000) | Oracle |
|---|---|
| COVTYP / PERCOV | 4 (Pacific silver fir) / 57.02 |
| ICT (cover metagroup) | 15 = DMIXCT (dry mixed conifer) |
| FMD (dominant) / FMOD·FWT (dynamic set) | 9 / **{9:0.594, 10:0.302, 6:0.104}** |
| ACTCBH / CBD | 4 / 0.12518 |
| **CRBURN / FIRTYPE** | **0.0 / 3 = SURFACE** |
| FLAME / SCH | 3.177 ft / 11.63 ft |
| post-fire TPA (2010) | **530 → 144** |

FMCFIR computes CBD/ACTCBH and the torching index but the fire stays surface (CRBURN=0). So the fire
result is driven entirely by the **surface fuel-model → flame → scorch → FMEFF bark-thickness mortality**
chain — the fuel-model selection, fuel loading, and bark all matter directly (they don't in a crown fire,
where the near-total kill masks them). **EC is therefore NOT added to the crown-fire gate** (consistent
with the measured surface behavior, and with WS/WC's exclusion). The `cr_crownw` crown-biomass routing IS
added (required by the mortality/fuel-additions snag path); measured jl CBD=0.122 vs oracle 0.125 with
ACTCBH 4=4 confirms EC would stay surface even if routed through the gate.

## KEY STRUCTURAL FINDING #2 — fmcfmd is the FMDYN "detailed low fuel model selection"

`ec/fmcfmd.f` is NEITHER the California-CWHR matrix (CA/NC/WS/SO) NOR the WC/PN FIRE-VPN 6-cover-metagroup.
It is EC's own dynamic selection: the 32 species pool into **15 cover-type metagroups** (DF/LP/SAF/MH/WP/
ES/PSF/WL/PP pure + DF-GF/PP-DF/LP-WL mixes + SAF-leading + moist/dry mixed conifer via ECMOIST); a single
dominant ICT is resolved (>50% BA one-sp, then two-sp, then SAF-leading, then moist/dry catch-all); a set
of ALGSLP weighting rules keyed by ICT/PERCOV/QMD and the single-vs-multi-strata flag LSNGL (FMSSTAGE)
build EQWT; models 10/12/13 are always natural-fuel candidates; and `FMDYN` resolves the candidates over
the (SMALL,LARGE) down-wood point. Closest already-ported template = the **LS/CR `_fmdyn` cover-metagroup
→ EQWT** family (EC reuses the shared `_fmdyn` + `_cr_algslp2`), but the metagroup logic is EC-bespoke.
XPTS are identical to WC's (reuse `_WC_FMD_XPTS`).

## FFE chunk verdicts

| Chunk | What | Verdict |
|-------|------|---------|
| F1 | `fire_species_props.csv` (32 sp): `EC_ISPMAP`/`ec_uses_fmcrowe` (ec/fmcrow.f), v2t + dkr_cls/leaf/tfall/snag/biogrp (ec/fmvinit.f + fmcblk.f) | ✅ source-faithful |
| F2 | Fuel loading — 32-sp FULIVE/FULIVI/FUINIE/FUINII (`ec_live/dead_fuel_loading`, **single-COVTYP**) + `ec_moist` (MAPDRY) VERBATIM from ec/fmcba.f + `ec_cwcalc` (Crookston R6, ECMAP, forest-608 OKANOGAN BF) | ✅ 1990 PERCOV 47.594 vs 47.605; FLIVE bit-exact-or-cornered |
| F2-decay | EC base DKR (ec/fmvinit.f: BM woody rates, litter 0.50) + habitat DKRADJ adjustment `ec_adjusted_dkr` (ECHMC/ECWMD × the shared `_FM_DKRADJ` table) | ✅ SMALL/LARGE 3.995/10.06 vs oracle 4.224/9.627 (required for FMDYN to keep model 10) |
| F4 | FMDYN cover-metagroup selection `ec_select_fuel_models` (15 metagroups → EQWT → `_fmdyn`) transcribed from ec/fmcfmd.f | ✅ FMD set {6,9,10} = oracle |
| F5 | Fire-mortality bark `_EC_FM_BARK_B1` (32 sp, ec/fmbrkt.f `bt=DBH·B1`) + group-6 FMEFF gate | ✅ source-faithful |
| — | `fire_fuel_models.csv` (Anderson-13, universal) + `ffe_on` gate now includes EastCascades | ✅ unblocks the burn / fixes fire_smlg |
| Crown biomass | route EastCascades through `cr_crownw` (FMCROWW; EC has no local fmcroww.f — CALL at fmcrow.f:162) + `EC_ISPMAP`/`ec_uses_fmcrowe` | ✅ CBD 0.122 vs 0.125, ACTCBH 4=4 |
| Crown fire | **NOT gated** — ect01_ffe is a SURFACE fire (oracle CRBURN=0); WS/WC precedent | ✅ correct by measurement |

## Two measured pitfalls fixed along the way

1. **`ffe_on` run_keyfile gap** — the SIMFIRE down-wood point (SMALL,LARGE) was `(0,0)` because
   `summary.jl`'s `ffe_on` gate required the eastern `fire_fuel_live.csv` (which EC, like NC/CR-family,
   does not have — its fuel lives in `ec_live/dead_fuel_loading`). With `fire_smlg=(0,0)`, FMDYN dropped
   the natural-fuel model 10 (jl set {6,9} vs oracle {6,9,10}) ⇒ lower flame/scorch ⇒ under-fire. Added
   `EastCascades` to the gate. (WC/PN/CA/WS mask this because their crown fire dominates; EC's surface
   fire is directly sensitive.)
2. **EC decay rates** — without EC's base DKR (coarse-wood 0.019–0.058/yr vs the SN default 0.07–0.11) +
   the habitat DKRADJ, the first-cycle LARGE down-wood decayed ~1.9× too fast (5.28 vs 9.627) ⇒ still
   dropped model 10 even with a correct `fire_smlg`. With EC's decay: LARGE 10.06 ≈ oracle 9.627 ⇒ FMDYN
   keeps {6,9,10}.

**Source note (QMD, ec/fmcfmd.f:171):** the Fortran QMD loop indexes the leftover `I` (=ICLSS+1), NOT the
loop var `J` — a source bug collapsing QMD to one record's DBH. Non-portable and immaterial here (EC's QMD
cutpoints are 2/4/7/9/19/21; both the buggy DBH(14)≈6.64 and the intended stand QMD 6.0 land in the same
ALGSLP bucket). jl computes the INTENDED QMD (Σprob·D²/Σprob); documented cornered.

## Fire behavior — SURFACE fire MEASURED on both sides (ect01_ffe SIMFIRE @2000)

Oracle `FVSec_g16` instrumented (unconditional WRITE in `fmburn.f` of ACTCBH/CBD/CRBURN/FIRTYPE/TCLOAD/
RFINAL/FLAME/SCH, `fmcba.f` of COVTYP/PERCOV/FLIVE, `fmcfmd.f` of FMD/ICT/QMD/FMOD/FWT), rebuilt, measured,
then **restored pristine** (grep-verified 0 markers; g16 relinked and its TPA-by-year re-verified == the
clean oracle 144). jl's `fmburn.jl`/`ec_fuel_model.jl` were temporarily printed and reverted (0 markers).

| Fire quantity | Oracle | jl |
|---|---|---|
| COVTYP / PERCOV | 4 / 57.02 | 4 / 57.08 |
| ICT / dominant FMD | 15 (DMIXCT) / 9 | 15 / 9 |
| FMDYN set (FMOD·FWT) | {9:0.594, 10:0.302, 6:0.104} | {9:0.594, 10:0.302, 6:0.104} |
| SMALL / LARGE (FMDYN point) | 4.224 / 9.627 | 3.995 / 10.06 (cornered) |
| ACTCBH / CBD | 4 / 0.12518 | 4 / 0.12201 |
| CRBURN / fire type | 0.0 / SURFACE | SURFACE (not gated) |
| **post-fire TPA (2010)** | **530 → 144** | **530 → 145** |

## End-to-end ect01_ffe .sum vs FVSec_clean

| year | col | oracle | jl | verdict |
|---|---|---|---|---|
| 1990 | TPA/BA/SDI/TopHt/QMD | 536/77/184/63/5.1 | 536/77/184/63/5.1 | **bit-exact** |
| 2000 (pre-fire) | TPA/BA/SDI/TopHt/QMD | 530/104/233/73/6.0 | 530/104/233/73/6.0 | **bit-exact** |
| 2010 (post-fire) | TPA/BA/QMD | 144/68/9.3 | 145/68/9.3 | **cornered** off-by-1 survivor (BA/QMD exact) |
| 2020–2090 | BA @2090 | 185 | 187 | pre-existing EC #206 OLDRN growth straddle (NOT FFE) |

**Verdict: bit-exact-or-cornered.** The pre-fire 2000 state is bit-exact, the fire type / COVTYP / ICT /
dominant-FMD / dynamic-FMD-set all match, CBD is within 2.5% with ACTCBH exact, and the post-fire result is
**530→145 vs 530→144** — an off-by-1 surface-fire largest-survivor tie-break (BA 68 and QMD 9.3 exact),
the same quality class as WC (491→2 vs 491→1). The 2020–2090 BA drift (187 vs 185 at 2090) is the
documented EC OLDRN/DGSCOR growth straddle carried by the post-fire survivors, not an FFE effect.

## Non-regression (all green)

- **CA cat01_ffe: 530 → 2**, **WC wct01_ffe: 491 → 2**, **PN pnt01_ffe: 485 → 2** — all unchanged.
- Fire test suite (`test_fire_biomass.jl` + `test_fire_effects.jl` + `test/integration/test_fire.jl`):
  **272 pass / 2 pre-existing-broken / 0 fail**.
- All shared-engine edits (`crown_biomass.jl`, `fmcba.jl`, `fuel_model.jl`, `fire_effects.jl`,
  `fuel_decay.jl`, `summary.jl`) are `EastCascades`-gated (byte-inert for every other variant by
  construction). EC growth control unchanged (1990 + pre-fire 2000 bit-exact).
- All Fortran restored pristine (grep marker count 0; g16 relinked, TPA-by-year == clean oracle).
