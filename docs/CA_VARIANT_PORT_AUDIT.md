# CA (CentralCalifornia) variant port audit — FFE chunk verdicts

Variant: **CA / CentralCalifornia** (MAXSP=50, forest 610 = R6 IFOR=6, 10-yr cycle, DGSD=2.0).
Oracle: `FVSca_clean` (/workspace/.cawork), relinked from `bin/FVSca_buildDir`. Reference stand
`cat01_ffe` (S248112, forest 610, SIMFIRE @2000: SWIND=10, FMOIS=1, ATEMP=50).

Growth + volume: **bit-exact** vs FVSca_clean (prior CA beachhead; cyc0 536/184). This doc covers the
**FFE (Fire & Fuels Extension)** rollout, first sibling of the WS-FFE template (#223).

## FFE chunk verdicts

| Chunk | What | Commit | Verdict |
|-------|------|--------|---------|
| F1 | `fire_species_props.csv` (50 sp: v2t, tfall, leaf-life, decay/fall/snag classes, bark_eqnum, biogrp) from ca/fmvinit.f/fmcblk.f/fmcrow.f | cca0082 | ✅ clears the first two FFE crashes |
| F2 | Fuel loading — 50-sp FULIVE/FULIVI/FUINIE/FUINII, top-2-cover PERCOV interpolation (`ca_live/dead_fuel_loading`) + Anderson-13 fuel-model table + `is_sprouting` + crown-biomass bark dispatch (`wc_bratio`) | a03a3aa | ✅ source-faithful |
| F3 | Crown-biomass small-tree height-dub (`ca_htdbh_height`) — CA + WS | 4262654 | ✅ CA FFE **runs end-to-end** |
| F4a | R6 Crookston crown-width `ca_cwcalc` (50 CWEQN from ca/cwcalc.f CAMAP; DF/WF/SP/LP/PP via `_cr_r6m2`, BR via `_nc_donnelly`) | 5e73e04 | ✅ **PERCOV 0.27 → 44.22** (fuel loading moves from bare→established regime) |
| F5 | Fire-mortality bark thickness `_CA_FM_BARK_B1` (50 sp, ca/fmbrkt.f `FMBRKT=DBH·B1(ISP)`, indexed **by species**) + group-6 FMEFF gate | 5dce801 | ✅ source-faithful; values verified vs fmbrkt.f line-by-line |
| F4b | California-CWHR fuel-model selection `ca_select_fuel_models` (ca/fmcfmd.f) — CWHRFMD 11×18 + 50-sp IFT `_ca_ift`; reuses NC CWHR consts/FINDJSS/FINDMOD/FMDYN | dc47533 | ✅ selects the high-intensity models **(6, 0.78)+(10, 0.22)** for CWHR sz=3/dn=M, replacing the weak SN default |

All CA-specific ported chunks are **source-faithful and validated** (bit-exact where measurable; fuel
model selection, crown-width, and bark verified against the FVS source directly).

## Open residual — fire runs SURFACE where the oracle CROWNS (shared crown-fire path)

Measured end-to-end on cat01_ffe:

- **oracle @2000: 530 → 2** (99.6% kill; the 2 survivors are the *largest* trees, QMD 12.4). oracle @1995
  (mid-cycle-1) drives the stand to **0** — the oracle fire is a **near-total-kill crown fire**, largely
  date-independent.
- **jl @2000: 530 → 133** (75% kill). jl computes a **surface** fire: flame 4.16 ft, scorch 17.32 ft — enough
  to kill small trees but not the medium/large overstory.

Root cause (MEASURED via FVSca_g16 — the g16 DOES run the fire; the earlier "g16 crashes on the fire cycle"
was a misread of exit-10 = normal FVS completion). Instrumented ca/fmburn.f + ca/fmcfir.f dumps (restored
pristine after) captured the oracle's crown-fire internals at the SIMFIRE @2000 vs jl:

| FMCFIR input | Oracle | jl |
|---|---|---|
| ACTCBH (crown base ht) | **4 ft** | **12 ft** |
| CBD (crown bulk density) | **0.1294** | **0.0492** |
| INIT1 | 65.58 | ~341 |
| RINIT1 (critical torch spread) | 6.17 | 24.7 |
| OINIT1 (torching index) | **0.0** | 23.78 |
| OACT1 (crowning index) | 17.13 | 36.08 |
| SFRATE (surface spread @10mph) | 8.84 | 8.84 |
| CRBURN / FIRTYPE | **0.554 / PASSIVE** | 0 / SURFACE |
| FLAME / SCORCH | **16.42 / 72.3** | 4.16 / 17.32 |

The surface spread (8.84) and HPA are fine; the torching MATH is fine. The **single root cause is jl's
`canopy_bulk_density`**: it computes ACTCBH=12 / CBD=0.0492 where the oracle gets ACTCBH=4 / CBD=0.1294.
With the oracle's actcbh=4, RINIT1=60·INIT1/HPA=6.17 < surface spread 8.84 ⇒ OINIT1=0 ⇒ the fire torches at
any wind ⇒ PASSIVE crown fire (CRBURN=0.554, flame 16.4, scorch 72) ⇒ near-total kill. jl's actcbh=12 makes
RINIT1=24.7 > 8.84 ⇒ never torches ⇒ surface (flame 4.16) ⇒ under-kill.

Why jl's `canopy_bulk_density` is wrong for cat01: its crown-fuel profile is too thin/high — the 1-ft-layer
crown-fuel array `crfill[6]=12.3, crfill[10]=7.8, crfill[15]=75, crfill[20]=105` reaches the 30-lb/ac-ft
crown-base threshold only at layer 12, while the oracle reaches it at layer 4 (much more crown fuel packed
into the 4–12 ft band, cbd 0.129 vs 0.049). Sub-cause is the CA crown-fuel profile feeding it — the CA
`fmcrow`/`fmpocr` foliage biomass or its low-tree inclusion — NOT the shared crown-fire classification math.
(This overturns the earlier speculation that the gap was the shared `fmcfir` torching-index/HPA and that CA
was "correctly excluded from the crown gate": in fact CA SHOULD get the crown-fire path once `canopy_bulk_density`
is fixed to match the oracle's actcbh/cbd.) `fm_canopy_lsw` has no CA method (falls to `sp<=25`) — plausibly
part of the low-tree inclusion issue; to be confirmed by instrumenting the oracle's ca/fmpocr.f + ca/fmcrow.f
per-tree crown-fuel profile.

**Verdict:** CA FFE port is **complete at bar** for the CA-specific FUEL/FUEL-MODEL/CROWN-WIDTH/BARK subsystems
(all source-faithful). The residual post-fire under-kill is a **measured, localized jl `canopy_bulk_density`
(crown-fuel profile) bug for CA** (actcbh 12 vs 4, cbd 0.049 vs 0.129), tracked as #229 — NOT a shared
(surface-vs-crown at low wind), tracked with the NC crown-fire byram/torching-index open item — resolving it
closes CA and NC together and is gated on being able to instrument an oracle crown-fire event.
