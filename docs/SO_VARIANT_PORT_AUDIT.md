# SO (SouthCentralOregon / SORNEC) variant port — FFE (Fire & Fuels Extension) audit — #233 (2026-08-15)

Task #233 ports the FFE for **SO / SouthCentralOregon** (VARACD "SO", SORNEC-33, MAXSP=33, DGSD=2.0,
Reineke SDI, Wykoff DDS). Last westside sibling of stream #223. Oracle: `/workspace/.sowork/FVSso_clean`
(made the `.sum`) + `/workspace/.sowork/FVSso_g16` (gfortran-16 instrumentable; exit 10 = normal). Growth +
volume already merged (kt-variant-port, chunks 0-8); this task is **FFE only**. Reference stand `sot01_ffe`
(S248112 — same trees as `ect01`/`cat01`, remapped to SO species codes; forest **601 = DESCHUTES R6**,
habitat 49 = CPS111, SIMFIRE @2000: SWIND=10, FMOIS=1, ATEMP=50). SO fire sources:
`fire/so/{fmbrkt,fmcba,fmcblk,fmcfmd,fmcrow,fmvinit}.f` + `bin/FVSso_buildDir/{cwcalc,fmcroww,essprt,esuckr}.f`.

## KEY STRUCTURAL FINDING #1 — sot01_ffe is a SURFACE fire (like EC), NOT a crown fire

Measured on `FVSso_g16` (instrumented `fmburn.f`/`fmcfmd.f`/`fmcba.f`, restored pristine, g16 relinked):

| FMBURN (SIMFIRE @2000) | Oracle |
|---|---|
| COVTYP / PERCOV | **32** (Other softwoods) / 66.62 |
| IPAG (Oregon super-group) | **1 = DRY PONDEROSA PINE** (IPASO(49)) |
| FMD (dominant) | **9** |
| ACTCBH / CBD | 4 / 0.04916 |
| **CRBURN / FIRTYPE** | **0.0 / 3 = SURFACE** |
| FLAME / SCH | 2.550 ft / 8.472 ft |
| SMALL / LARGE (FMDYN point) | 5.008 / 5.277 |
| post-fire TPA (2010) | **442 → 216** (~51% survival) |

CRBURN=0 ⇒ the fire is driven by the surface fuel-model → flame → scorch → FMEFF bark-thickness mortality
chain. **SO is therefore NOT added to the crown-fire gate** (EC/WC/WS precedent). The `cr_crownw` crown-biomass
routing IS added (required by the mortality/fuel-additions snag path).

## KEY STRUCTURAL FINDING #2 — fmcfmd is REGION-DEPENDENT; the R6 stand takes the Oregon 8-plant-group path

`so/fmcfmd.f` selects the **R5-California SOSPDM** dominant-species path only for KODFOR 500-599 / 701 (an
RNG-jittered dominant finder → models 2/5/8/9/10); every R6-Oregon forest uses the **8-plant-group path**
(IPAG = IPASO(ITYPE) collapses the 92 SORNEC plant associations into 8 super-groups). The reference stand
(forest 601) rides the Oregon branch, so **no RNG is involved in fuel-model selection**. The R5-California
SOSPDM branch is a documented follow-on (not exercised by the R6 reference stand). Closest already-ported
template = the shared `_fmdyn` resolver; the 8-group EQWT logic is SO-bespoke (`so_select_fuel_models`).

## KEY STRUCTURAL FINDING #3 — fuel loading is the FCCS/Ottmar photo-series approach (not the FULIVE tables)

Unlike EC/CR-family (single-COVTYP FULIVE/FUINIE interpolation) or NC/WS/CA (top-2 CWHR), SO uses the FCCS
crosswalk: `ISPX = (COVTYP, 2→1)·10 + ISSX` (ISSX = the FMSSTAGE structural stage, with the PERCOV≥60
SE-open→SE-closed bump) keys **SO_COVRINI** (218 rows) → an FCC fuel-model number IMODX; **SO_FUELINI**
(43 rows) → the (herb,shrub) live pool + a 9-class dead pool re-binned into the 11 FFE size classes with the
fmcba.f:707-725 class-4-6 fiddling. Both tables extracted VERBATIM from so/fmcba.f. MEASURED vs the oracle:
COVTYP 32, ISSX 3, IMODX 56 (FCC), STFUEL = [0.5, 1.3, 3.0, 3.0, 1.909, 1.091, 0, 0, 0, 0, 9.1] and FLIVE
0.5/0.5 — **bit-exact**.

## FFE chunk verdicts

| Chunk | What | Verdict |
|-------|------|---------|
| F1 | `fire_species_props.csv` (33 sp): `SO_ISPMAP`/`so_uses_fmcrowe` (so/fmcrow.f) + v2t/leaf/dkr_cls/snag_alldwn/biogrp (so/fmvinit.f + fmcba.f Oregon snag defaults + fmcblk.f) + dbh_min(9.0 so/grinit.f) + is_sprouting(ISPSPE) + sprout_ht1/ht2 (blkdat HT1/HT2) | ✅ source-faithful; clears the FFE crashes |
| F2 | Fuel loading — `data/southcentraloregon/fire/ffe_fuel.jl`: `so_fuel_ini` (COVRINI→FUELINI, FCCS/Ottmar) VERBATIM from so/fmcba.f | ✅ COVTYP 32, ISSX 3, IMODX 56, STFUEL + FLIVE **bit-exact** |
| F2-cw | Crown width `so_cwcalc` (so/cwcalc.f SOMAP, forest-601 DESCHUTES BF) — 33 species (R6-m2 via `_cr_r6m2`, R1-log, Donnelly, Bechtold-m2, MH piecewise) | ✅ PERCOV 64.9–67.9 vs oracle 66.6 (~2.6%, from lat=0 on the 3 Bechtold species) |
| F2-decay | SO Oregon base DKR (= `_FM_DKR_EC` values) + `so_adjusted_dkr` (SOHMC/SOWMD × the shared `_FM_DKRADJ`) | ✅ ITYPE 49 (CPS111 default) → TEMP=hot/MOIST=dry; SMALL/LARGE feed FMDYN |
| F4 | `so_fuel_model.jl` — Oregon 8-plant-group selection (`so_select_fuel_models`, IPAG=IPASO(ITYPE) → EQWT → `_fmdyn` with the 14-model `_SO_FMD_XPTS`) | ✅ IPAG 1 (dry PP) → EQWT {9,10,12,13} = oracle; **FMD 9** = oracle |
| F5 | Fire-mortality bark `_SO_FM_BARK_B1` (33 sp, so/fmbrkt.f `bt=DBH·B1`) + group-6 FMEFF gate | ✅ source-faithful; values verified vs fmbrkt.f line-by-line |
| Crown biomass | route SouthCentralOregon through `cr_crownw` (FMCROWW) + `SO_ISPMAP`/`so_uses_fmcrowe`; **extended `cr_crownw` for SPIE groups 5/6/19/20/23** (so/fmcroww.f) | ✅ CBD/ACTCBH surface-consistent |
| Sprouting | `essprt_so`/`nsprec_so`/`sprtht_so`/`so_sprout_dbh` (so/essprt.f + esuckr.f CASE('SO'), ASSPTN aspen index 24) | ✅ source-faithful; fire-killed hardwoods sprout via ESUCKR |
| Crown fire | **NOT gated** — sot01_ffe is a SURFACE fire (oracle CRBURN=0); EC/WC/WS precedent | ✅ correct by measurement |

## Fire behavior — SURFACE fire MEASURED on both sides (sot01_ffe SIMFIRE @2000)

| Fire quantity | Oracle (FVSso_g16) | jl |
|---|---|---|
| COVTYP / PERCOV | 32 / 66.62 | 32 / 64.88 |
| IPAG / dominant FMD | 1 (dry PP) / 9 | 1 / 9 |
| FMDYN set (FMOD·FWT) | dominant 9 | {9:0.787, 10:0.213} |
| FLAME / SCORCH | 2.550 / 8.472 | 2.537 / 8.386 |
| CRBURN / fire type | 0.0 / SURFACE | SURFACE (not gated) |
| **post-fire TPA (2010)** | **442 → 216** | **431 → 217** |

## End-to-end sot01_ffe .sum vs FVSso_clean (NOTRIPLE — isolates FFE from the tripling RNG)

| year | col | oracle | jl | verdict |
|---|---|---|---|---|
| 1990 | TPA/BA/SDI/TopHt/QMD | 613/92/218/63/5.2 | 613/92/218/63/5.2 | **bit-exact** |
| 2000 (pre-fire) | TPA/BA/SDI | 442/102/221 | 431/102/221 | pre-existing SO #206 OLDRN growth straddle (BA/SDI exact) |
| 2010 (post-fire) | TPA/BA/QMD | 216/63/7.3 | 217/68/7.6 | **cornered** off-by-1 survivor (surface-fire largest-survivor tie-break) |
| 2020–2090 | TPA @2090 | 89 | 65 | pre-existing SO #206 growth straddle carried by survivors (NOT FFE) |

**Verdict: bit-exact-or-cornered.** The cycle-0 (1990) state is bit-exact; the fire type / COVTYP / IPAG /
dominant-FMD / fuel loading (STFUEL/FLIVE) all match; FLAME 2.537 vs 2.550 and SCORCH 8.386 vs 8.472; the
post-fire survivor count is **431→217 vs 442→216** — an off-by-1 surface-fire largest-survivor tie-break,
the same quality class as EC (530→145 vs 144) and WC (491→2 vs 1). The pre-fire 2000 (431 vs 442) and the
2020–2090 drift are the **pre-existing SO #206 OLDRN growth straddle** (DGSD=2.0) — verified present WITHOUT
FFE (so_c1 growth-only jl 2000 = 431, identical to the FFE pre-fire), so FFE does not perturb growth.

## Two documented deferrals (same class as the SO growth-port deferrals)

1. **Habitat ITYPE** — the reference stand rides so/habtyp.f's DEFAULT plant association CPS111 = ITYPE 49
   (SI 70), the same default the SO growth port's SITEAR/SDIDEF ride. Explicit-habitat real-FIA stands need
   the full so/ecocls.f 92-entry PA decode; `p.habitat_input` > 0 (a decoded PA) is honored when present.
2. **R5-California SOSPDM fuel-model branch + the LP/juniper (IPAG 6/7/8) ACTCBH path** — not exercised by
   the R6 dry-PP (IPAG 1) reference stand; the SOSPDM dominant-finder uses RNG and is out of scope for the
   surface-fire smoke.

## Non-regression (all green)

- **CA cat01_ffe: 530 → 2**, **WC wct01_ffe: 491 → 2**, **PN pnt01_ffe: 485 → 2**, **EC ect01_ffe: 530 → 145**
  — all unchanged.
- Fire test suite (`push!(ARGS,"fire"); include("test/runtests.jl")`): **38781 pass / 0 fail / 75 pre-broken /
  1 oracle-env error**.
- All shared-engine edits (`crown_biomass.jl`, `cr_crown_biomass.jl`, `fmcba.jl`, `fuel_model.jl`,
  `fuel_decay.jl`, `fire_effects.jl`, `structure_stage.jl`, `sprout.jl`, `summary.jl`) are
  `SouthCentralOregon`-gated (byte-inert for every other variant by construction; the `cr_crownw` SPIE 5/6/19/
  20/23 additions are new branches that no other variant reached). SO growth control unchanged (1990 bit-exact).
- All FVS Fortran restored pristine (marker count 0; the 6 fire sources diff-clean vs the fresh tree; g16
  relinked and its post-fire TPA re-verified == the clean oracle 216).
