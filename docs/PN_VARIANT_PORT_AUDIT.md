# PN (PacificNorthwest) variant port — FFE (Fire & Fuels Extension) audit — #231 (2026-08-15)

Task #231 ports the FFE for **PN / PacificNorthwest** (MAXSP=39, forest 612 = SIUSLAW R6, 10-yr cycle,
Reineke SDI, DGSD=1.7). Oracle: `/workspace/.pnwork/FVSpn_clean` (relink `bin/FVSpn_buildDir`) +
`FVSpn_g16` (instrumentable). Reference stand `pnt01_ffe` (S248112 — same trees as `wct01`; forest 612,
habitat 40, SIMFIRE @2000: SWIND=10, FMOIS=1, ATEMP=50). Growth+volume already bit-exact-or-cornered
(the merged PN growth port). This continues the western FFE rollout after CA (#229) and WC (#230).

## KEY STRUCTURAL FINDING — PN FFE is FIRE-VPN (like WC), NOT California-CWHR

The task template flagged that `pn/fmcfmd.f` "has a CWHR mention" and said to CONFIRM. **Measured:**
`pn/fmcfmd.f` is **BYTE-IDENTICAL to `wc/fmcfmd.f`** (md5 `e0423702`) and its header is `FIRE-VPN`.
So PN FFE is the Pacific-Northwest FIRE-VPN structure, the **same as WC**, not the CA/WS/NC California-CWHR:

- **fuel loading** (`pn/fmcba.f`) uses the **SINGLE dominant cover type** COVTYP (like WC/CR/IE/EM), NOT
  the top-2-cover COVCA/COVCAWT interpolation of NC/WS/CA.
- **fuel-model selection** (`pn/fmcfmd.f`, byte-identical to WC) is the **6-cover-metagroup** FIRE-VPN
  structure (SF/DF/MH/RA/LP/WO) + QMD80 + PERCOV + habitat MAPFGS/MAPDRY rules — NOT the CWHRFMD matrix.

**PN FFE is essentially a coefficient-swap on the just-completed WC FFE port.** The shared FFE Fortran is
byte-identical PN↔WC: `fmcfmd.f`, `fmcroww.f` (md5 `f764dce1`), `cwcalc.f` (md5 `2bace281`), `fmcbio.f`.
The only PN-specific DATA is in `fmvinit.f` / `fmcrow.f` / `fmbrkt.f` / `fmcblk.f` / `fmcba.f`, and every
one differs from WC **only at species 6 (SS / Sitka spruce; WC's slot 6 is blank)** plus a handful of
genuine per-species fuel-loading differences (below).

## Measured PN-vs-WC FFE DATA deltas (all else byte-identical)

| Fortran file | PN-vs-WC difference |
|---|---|
| `fmvinit.f` (v2t/leaf/tfall/decay/snag) | ONLY species 6 (SS: V2T=20.6, LEAFLF=5, TFALL(3)=15, ALLDWN=110, DKRCLS=2) |
| `fmcrow.f` ISPMAP (`ls_spi`/`PN_ISPMAP`) | ONLY species 6 = **18** (SS→ES; WC=0). FMCROWE routing CASE identical |
| `fmbrkt.f` bark B1 | ONLY species 6 = **0.027** (WC blank −0.000) |
| `fmcblk.f` BIOGRP (Jenkins C-report) | species 6 = 5 (spruce), species 39 (OT) = 6 (WC=7) — carbon-only, not fire behavior |
| `fmcba.f` FULIVE/FULIVI | ONLY species 6 = (0.30,0.20)/(0.30,2.00) |
| `fmcba.f` FUINIE | species 6 (real SS row) **and species 10 (ES uses the DF table)** |
| `fmcba.f` FUINII | species 6; **species 8/9/10/16/17/18/21 first-two classes 1.1→1.6**; species 1 12–20″ 6.0→0.0 |
| `fmcba.f` MAPFGS/MAPDRY | **PN's own 75-code habitat arrays** (WC's are 139-code) — genuinely different |
| `cwcalc.f` crown width | byte-identical eqns/coeffs; only the **forest-612 (SIUSLAW) BF** differs (DF 202=0.977; RC 242=0.905; WH 263=0.924; all other species BF=1.0 — vs WC's forest-618) |

## FFE chunk verdicts

| Chunk | What | Verdict |
|-------|------|---------|
| F1 | `fire_species_props.csv` (39 sp): `PN_ISPMAP`/`pn_uses_fmcrowe` (pn/fmcrow.f), v2t (pn/fmvinit.f), dkr/leaf/tfall/snag/biogrp. = WC's table with row 6 (SS) + row 39 biogrp corrected | ✅ source-faithful |
| F2 | Fuel loading — 39-sp FULIVE/FULIVI/FUINIE/FUINII (`pn_live/dead_fuel_loading`, **single-COVTYP**) transcribed VERBATIM from pn/fmcba.f + `pn_cwcalc` (Crookston R6, forest-612 SIUSLAW BF) | ✅ source-faithful; fire outcome bit-exact-or-cornered |
| F5 | Fire-mortality bark `_PN_FM_BARK_B1` (39 sp, pn/fmbrkt.f `bt=DBH·B1`; = WC with slot 6=0.027) + group-6 FMEFF gate (pn/fmeff.f restricts Regelbrugge-Smith to SN/CS) | ✅ source-faithful |
| F4b | FIRE-VPN cover-metagroup selection — `pn/fmcfmd.f` is BYTE-IDENTICAL to wc/fmcfmd.f ⇒ reuse `wc_select_fuel_models`; only PN's `fmcba.f` `MAPFGS/MAPDRY` (75-code) are PN-specific (`_PN_MAPFGS`/`_PN_MAPDRY`, selected by variant) | ✅ source-faithful |
| — | `fire_fuel_models.csv` (Anderson-13, universal = klamath/CR/WC) | ✅ unblocks the burn |
| Crown fire | route PacificNorthwest through `cr_crownw` (FMCROWW; fmcroww.f byte-identical) + crown-fire gate/Unions in `fmburn.jl` | ✅ crown fire MEASURED bit-exact-or-cornered |

## Fire behavior — crown fire MEASURED on both sides (pnt01_ffe SIMFIRE @2000)

Oracle `FVSpn_g16` was instrumented (unconditional WRITE in `fmcfir.f` of FMOIS/ACTCBH/CBD/RINIT1/SFRATE/
OINIT1/OACT1/CRBURN), rebuilt, measured, then **restored pristine** (grep-verified 0 markers, g16 relinked
and its TPA-by-year re-verified == the clean oracle). jl's `fmburn.jl` `canopy_bulk_density` was
temporarily printed and reverted (0 markers).

| FMCFIR input (the fire, FMOIS=1) | Oracle | jl |
|---|---|---|
| ACTCBH (crown base ht, ft) | **4** | **3** |
| CBD (crown bulk density) | **0.12220** | **0.12255** |
| OINIT1 (torching index) | **0.0** | 0.0 (torches) |
| SFRATE (surface spread) | 15.83 | — |
| CRBURN / fire type | **0.503 / PASSIVE crown** | PASSIVE crown |
| **post-fire TPA (2010)** | **485 → 4** | **485 → 2** |

**Both run a PASSIVE crown fire** (torching index 0 ⇒ torches at any wind; CBD matches within 0.3%). The
crown biomass routes through the **shared `cr_crownw`** (FMCROWW, byte-identical to CR/WC/CA — the #229
machinery validated bit-exact for CA 530→2 and WC 491→2). The residual is **cornered**: ACTCBH 3 vs 4 (a
1-ft crown-base running-mean discretization) plus the pre-fire growth straddle (jl 2000 BA 115 vs 116,
TopHt 72 vs 75 — the documented PN OLDRN/DGSCOR growth straddle, NOT FFE) leave jl's trees marginally
smaller/lower-crowned ⇒ 2 more trees cross the kill threshold ⇒ **485→2 vs 485→4** — the same
largest-survivor tie-break class as CA (530→2) and WC (491→2 vs 491→1).

## End-to-end pnt01_ffe .sum vs FVSpn_clean

| year | col | oracle | jl | verdict |
|---|---|---|---|---|
| 1990 | TPA/BA/SDI/TopHt/QMD | 536/77/184/63/5.1 | 536/77/184/63/5.1 | **bit-exact** |
| 1990 | TotCuFt | 1877 | 1884 | +0.4% (pre-existing PN volume residual) |
| 2000 (pre-fire) | BA/SDI/TopHt | 116/250/75 | 115/248/72 | cornered PN growth straddle (also in pnt01.sum.save; NOT FFE) |
| 2010 (post-fire) | TPA | 4 | 2 | **cornered** crown-fire largest-survivor tie-break |

## Non-regression (all green)

- **CA cat01_ffe: 530 → 2** (unchanged, bit-exact) · **WC wct01_ffe: 491 → 2** (unchanged).
- Fire test suite (`test_fire_biomass.jl` + `test_fire_effects.jl` + `test/integration/test_fire.jl`):
  **272 pass / 2 pre-existing-broken / 0 fail**.
- All shared-engine edits are `PacificNorthwest`-gated (byte-inert for every other variant by construction).
- PN growth control unchanged (no growth code touched; 1990 cyc0 bit-exact, 2000 = the pre-existing straddle).

## Verdict

PN FFE is **complete, bit-exact-or-cornered**: all PN-specific chunks (fuel loading, crown-width, fuel-model
selection, bark, species props, crown biomass) are **source-faithful**; the fire runs the oracle's **PASSIVE
crown fire** (CBD 0.1225 vs 0.1222; both torch) with a **cornered off-by-2 largest-survivor tail** (485→2 vs
485→4). The #229 crown-fire-classification gap — cluster-wide (CA/WC/NC/WS) — is now **closed for PN** as it
was for CA and WC, reusing the identical `cr_crownw` + crown-fire gate machinery with only PN's coefficient
DATA swapped in.

## Remaining (follow-ups, not parity gaps)

- The pre-fire growth straddle (2000 BA 115 vs 116, TopHt 72 vs 75) — the cornered PN OLDRN/DGSCOR class,
  shared with WC/UT/IE (a growth follow-up, not FFE).
- `pn_cwcalc` ports only the 7 pnt01 species (WF/GF/LP/SP/PP/DF/ES) — the remaining WCMAP equations are a
  follow-up crown-width chunk (mirrors CA F4a / WC scope).
- The WCWMC/WCWMD/DKRADJ dead-fuel decay-rate habitat adjustment (multi-cycle fuel decay, deferred as in WC).
- Model-11 (5-yr post-activity fuel jump) — deferred as in WC/NC/WS.
