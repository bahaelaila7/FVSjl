# WS (WestSierra) Variant Port — Audit

Branch: `kt-variant-port`. Oracle: `FVSws_clean` (production relink) + `FVSws_g16` (gfortran-16
instrumentable rebuild), `/workspace/.wswork`. MAXSP=43, Zeide SDI (LZEIDE, no reset), DGSD=2.0, 10-yr
period, seed 55329, R5 California forest 500 (wst01 = forest 511). Doctrine: bit-exact vs LIVE oracle per
chunk; MEASURE don't infer; per-record treelist INVALID after tripling (use `.sum`).

## Growth + volume (chunks 0–8) — BIT-EXACT-OR-CORNERED
Validated in prior sessions (see memory `fvsjl-ws-westsierra-scoping`). cyc0 TPA/BA/QMD/TopHt + MCuFt/BdFt
bit-exact vs `FVSws_g16` controlled-input; multi-cycle within the #206 DGSD=2.0 OLDRN serial-corr straddle.
Notable fixes: `ws_sitset!` BAMAX/PMSDIU; `ws_forkod!` forest-dependent TLAT (ILAT band); CRATET missing-height
dubbing (Wykoff native + HTDBH MODE=0 surrogate); TPA over-mortalization. `:v2t` cut-biomass (thinning path)
is a cluster-wide follow-on.

## FFE (Fire & Fuels Extension) — BIT-EXACT-OR-CORNERED end-to-end

WS FFE was unported (the "runs end-to-end" claim was optimistic — the FFE fuel init crashed on an empty table).
Ported chunk by chunk, each measured against the live oracle:

| chunk | commit | what | verdict |
|---|---|---|---|
| F0 | 8e620e9 | crown-biomass species map (WS_ISPMAP, fmcrow) + `:dbh_min` | advances FFE path |
| F1 | a81a790 | `fire_species_props.csv` (12 cols × 43 sp, ws/fmvinit.f) | FFE runs (no crash) |
| F2 | a585e3f | initial fuel loading (FULIVE/FULIVI/FUINIE/FUINII top-2 cover, ws/fmcba.f) + **R5CRWD** crown-width (ws/r5crwd.f — WS branches to R5CRWD, NOT the WSMAP Crookston eqns) | **cyc0 fuel BIT-EXACT** |
| F3 | c205fd6 | fire path: standard Anderson-13 `fire_fuel_models.csv` + merch standards (ws/grinit.f) + bark dispatch → ws_bratio | fire path runs |
| F4 | c0a2c5f | **fmcfmd California-CWHR fuel-model selection** (ws/fmcfmd.f + cwhr.f): CWHRFMD 12×18 (models 25/26) + 43-sp IFT + CWHR consts; parameterized nc_cwhr → shared `_ca_cwhr` | models 8+11 = oracle |
| F5 | d62ed1c | fire-mortality **bark thickness** (`_WS_FM_BARK_B1`, ws/fmbrkt.f) + group-6 gate (ws/fmeff.f:196) | over-kill RESOLVED |

### cyc0 fuel loading — BIT-EXACT (measured vs instrumented FVSws_g16, restored pristine)
PERCOV jl 41.99 = oracle 41.986; COVCA 1(SP)/3(WF), COVCAWT 0.6153/0.3847 = oracle exactly; DUFF 15.27=15.272,
LITT 0.49, dead 0-3 3.14, 3-6 5.49, 6-12 5.97, HERB 0.246, SHRUB 0.563 — every ALL-FUELS cyc0 bucket matches.

### Fire behaviour + mortality (wst01 + SIMFIRE 2000) — validated vs FVS_BurnReport DBS + `.sum`
The over-kill was chased to ground by measurement (not assumption): the wrong SN-default fuel models
(6+10, flame 4.20) → fmcfmd port (8+11) → bark-thickness fix. Post-fire 2010 TPA **73 → 112 → 209** (oracle 218).

Full trajectory (jl / oracle): 2010 TPA 209/218 BA 205/194; 2020 205/214 270/250; 2040 181/203 372/355;
2090 59/83 330/354. Fire-year mortality tracks within ~4%; the late-cycle TPA drift with BA tracking is the
OLDRN self-thin straddle signature (#206 class) — cornered.

### Cross-variant note (westside FFE template)
The F2–F5 gaps (missing fuel-model table / merch standards / bark dispatch / fmcfmd / bark thickness / group-6)
are the **westside-FFE-port template** — the sibling fresh westside ports (NC/WC/PN/EC/CA/SO) will need the
same chunks when their FFE is brought up. The bark-intercept crash + bark-thickness over-kill are the exact
NC/CR/BM bug class.

## Open (low priority, not reference-parity blockers)
- **#219** FFE multi-cycle DCYMLT (ws/fmcba.f:543-580 Dunning-site-idx decay mult) — ~1-ton DUFF drift.
- Fuel-model weight resolution: jl production picks FM8@100% vs oracle 8@0.78+11@0.22 (fire-basis sm/lg →
  `_fmdyn` weights); contributes the residual ~4-6% in the fire mortality. LOW.
- `:v2t` cut-biomass (thinning path) — cluster-wide gap.
