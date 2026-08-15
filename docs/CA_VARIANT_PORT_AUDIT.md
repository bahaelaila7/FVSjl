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

## ✅ #229 RESOLVED (a4a0815) — CA crown-fire now BIT-EXACT (530→2)

The crown-fire under-kill below is **FIXED**. Root cause (measured, table retained for the record):
`crown_biomass` omitted CentralCalifornia from the `cr_crownw` (FMCROWW) dispatch → CA fell to the Jenkins
FMCROWE eastern path → conifer crown biomass 6–20× too low → `canopy_bulk_density` gave actcbh=12/cbd=0.049
(vs oracle 4/0.129) → RINIT1=24.7 > surface spread 8.84 → never torched → surface fire.

**Fix** (two commits): `2c45d0c` (part 1) added `CA_ISPMAP`/`ca_uses_fmcrowe` from ca/fmcrow.f;
`a4a0815` (part 2) wired CentralCalifornia into the `cr_crownw` gate + spie/spils selectors in
`crown_biomass.jl`, and into the crown-fire gate (fmburn.jl:138) + the 3 crown-fire dispatch Unions
(crowning_index/torching_index/crown_fire_result → shared `crown_fire_result`, mirroring WestSierra).

**VALIDATED:** `cat01_ffe` 2000: 530 → 2010: **2** — BIT-EXACT vs FVSca_clean (the 2 survivors are the
largest trees, growing through 2090). Non-regression: fire suite **272 pass / 2 pre-existing-broken / 0 fail**;
all edits are CentralCalifornia-gated (inert for CR/NC/WS/BM/eastern by construction). ⇒ **CA FFE is now
complete-at-bar end-to-end** (fuel loading, fuel-model selection, crown-width, bark, AND crown fire).

NOTE: WS and WC remain **excluded** from the crown-fire gate — they carry the same latent surface-vs-crown
gap on crowning stands (their crown biomass is correct via `cr_crownw`, but they don't enter the crown
classification). That is a separate follow-up (each variant's own verdict), not chased here.

## Open residual — fire runs SURFACE where the oracle CROWNS (shared crown-fire path) — ⬆ RESOLVED, see above

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

SUB-CAUSE PINNED (instrumented oracle ca/fmpocr.f + ca/fmcrow.f, restored pristine). The FMPOCR algorithm is
byte-for-byte identical to jl's (CBHCUT=30, 13-ft CBD running mean, 3-ft crown-base mean; CA has no LBHPP case
so it uses the uniform crown-fill — same as jl). The ONLY difference is the per-tree crown biomass `CROWNW`
feeding the profile — jl's `crown_biomass` is **6–20× too low on foliage and ≈0 on the 1-hr woody**:

| tree (D",H,ICR) | jl fol / w1 | oracle CROWNW0 / CROWNW1 |
|---|---|---|
| SP12 4.85/86.5/39 | 6.28 / 0.01 | 39.5 / 33.25 |
| SP4  2.97/54.9/50 | 2.25 / 0.01 | 48.0 / 19.7 |
| SP7  2.09/33.5/28 | 1.23 / 0.01 | 19.0 / 9.88 |

ROOT: `crown_biomass` (src/engine/fire/crown_biomass.jl:103-104) routes CR/BM/Klamath/WestSierra through
`cr_crownw` (the FMCROWW western engine) but **omits CentralCalifornia** — so CA falls through to the Jenkins
FMCROWE eastern-hardwood path (via `ls_spi`), which gives tiny conifer foliage and ~0 woody. `ca/fmcroww.f`
is BYTE-IDENTICAL (md5 f764dce1) to CR/WC/WS's, so CA reuses `cr_crownw` directly. **FIX** (turnkey, deferred
until the WC-FFE agent frees crown_biomass.jl, which it is editing for WC's own routing): add a `CA_ISPMAP[50]`
(`7,20,7,4,4,4,3,6,24,14, 11,11,11,11,15,15,15,13,13,11, 16,18,19,7,11,17,17,21,17,21, 21,21,17,5,44,23,10,17,
56,29, 46,17,60,41,17,64,17,17,21,19` from ca/fmcrow.f) + `ca_uses_fmcrowe(sp)= sp∈{35,39,40,41,43,44,45,46}`
+ route CentralCalifornia through `cr_crownw` (mirror WestSierra). Then add CA to the crown-fire gate/Unions
(fmburn.jl:138/359/379/412) and validate cat01_ffe 530→2; regress NC (actcbh 6=6)/WS/CR. NOTE the F2 fuel
loading was unaffected because it uses the FULIVE/FUINI tables (`ca_live/dead_fuel_loading`), not `crown_biomass`
— which is why the cyc0 PERCOV validation didn't catch this. This overturns the earlier "shared fmcfir
classification, CA correctly excluded from the crown gate" framing.

**Verdict:** CA FFE port is **complete at bar** for the CA-specific FUEL/FUEL-MODEL/CROWN-WIDTH/BARK subsystems
(all source-faithful). The residual post-fire under-kill is a **measured, localized jl `canopy_bulk_density`
(crown-fuel profile) bug for CA** (actcbh 12 vs 4, cbd 0.049 vs 0.129), tracked as #229 — NOT a shared
(surface-vs-crown at low wind), tracked with the NC crown-fire byram/torching-index open item — resolving it
closes CA and NC together and is gated on being able to instrument an oracle crown-fire event.
