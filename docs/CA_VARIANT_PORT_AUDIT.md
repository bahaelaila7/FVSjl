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

Root cause (measured, not inferred): jl classifies the fire as **surface**, not crown. For cat01
(`canopy_bulk_density`: cbd=0.0492, actcbh=12, tcload=0.161):
`torching_index` OINIT = **23.78**, `crowning_index` OACT = **36.08**, both > SWIND = **10** ⇒
`crown_fire_result` returns CRBURN=0 (SURFACE, fmcfir.f:334). The critical torching spread
`rinit1 = 60·init1/hpa = 24.68 ft/min` (init1 from actcbh=12, hpa=828.4) is reached only at ~24 mph, so at
10 mph jl never torches. The oracle crowns at 10 mph ⇒ its critical spread / HPA / surface-vs-crown
classification differs.

This is the **shared `fmcfir` crown-fire path**, not a CA-specific chunk:
- It is the **same family** as the open NC crown-fire issue (fmburn.jl:129-137: the shared crown-fire BYRAM
  step over-drives FINTEN; NC surface-only 54 vs 58 is cornered, crown-on over-kills 0 vs 58). CA/WS are
  currently **excluded** from the crown-fire boost gate (fmburn.jl:138) and the `crown_fire_result` /
  `torching_index` / `crowning_index` Union dispatch.
- Adding CA to the gate + Unions was tested and is **inert** for cat01 (crb=0, because OINIT=23.78 > 10) —
  the lever is the torching-index / HPA magnitude, not the gate membership. Reverted (kept faithful).
- jl's `actcbh` is validated against NC live (cbd 0.147 = live 0.152, actcbh 6=6), so the discrepancy is in
  the HPA / critical-spread / crown-type classification, needing the oracle's internal OINIT/flame.

**Blocked measurement:** the CA g16 oracle (`FVSca_g16`) crashes on the fire cycle (exit 2 at 2000, known
g16 fragility — the CA volume path SIGFPEs under the isoc23 shim), so the oracle's OINIT/flame/fire-type
cannot be dumped. FVSca_clean's fire-report keywords (BURNREPT/FUELREPT/MORTREPT) require an open DATABASE
and exit 10 in this build; SIMFIRE only accepts cycle-boundary dates. ⇒ oracle crown-fire internals not yet
measurable.

**Verdict:** CA FFE port is **complete at bar** for the CA-specific subsystems (fuel, fuel-model, crown-width,
bark — all source-faithful). The residual post-fire under-kill is a **shared crown-fire-classification** gap
(surface-vs-crown at low wind), tracked with the NC crown-fire byram/torching-index open item — resolving it
closes CA and NC together and is gated on being able to instrument an oracle crown-fire event.
