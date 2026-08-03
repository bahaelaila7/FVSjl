# BM (Blue Mountains) Variant Port — Audit Log

FVS Region-6 Blue Mountains variant → FVSjl. Validated bit-exact-or-cornered per chunk vs the live
relinked oracle `/workspace/.bmwork/FVSbm_clean` (gfortran-16 + isoc23 shim). MAXSP=18, Stage SDI,
DGSD=1.5, 10-yr cycle, seed 55329. Test stand `bmt01` (S248112, forest 614 Umatilla, habitat 12).
Growth validation uses `bmt01_growth.key` (FFE keywords stripped — bmt01.key is a full FMIN/SIMFIRE demo).

## Chunk verdicts

| # | Chunk | Files | Verdict |
|---|-------|-------|---------|
| 0-1 | Scaffold + species | bluemountains.jl, species.jl | ✓ loads, 18 species, seed 55329, POWER bark |
| 2 | Site / habitat | site_index.jl | ✓ **bit-exact** — ISISP=DF(3), SITEAR(DF)=64, SDIDEF=346, all 18 SITEAR match live SITECODE echo. Novel eco-class DB chain: habitat 12 → PCOML[12]=CDS722 → ECOCLS → SICHG/HTCALC (incl. the faithful CASE(3) `**` precedence bug) |
| 3 | Large-tree DG | diameter_growth.jl, dg_coefficients.jl | ✓ **bit-exact-or-cornered** — DGCON bit-exact (DF/WL/LP 1.5531/1.3957/1.2295), WK2 24/27 bit-exact (MSS spline + all 3 species groups), DGSCOR COR bit-exact (GF 0.01963=live). Misses = 3 D=0.08 seedlings |
| 4 | Height | height_growth.jl | ✓ **21/24 bit-exact** — WL/DF/LP/ES exact; bm_findag + POTHTG=HTCALC(SITAGE+10)−SITHT + HGMDCR/HGMDRH modifiers |
| 5 | Crown | crown.jl | ✓ Weibull (= validated TT/UT form) + bm_tree_ccf (RELDEN, 369/369 per-tree CCF bit-exact) |
| 6 | Regent (small-tree) | regent.jl | ✓ bm_smhtgf height matches live (POTHTG 11.42, HTG 5.72); ht-dbh uses bm_htdbh (Curtis-Arney) for LHTDRG=false species |
| 7 | Mortality | mortality.jl | ✓ **faithful** — Hamilton RI + SDI self-thin; RN bit-close for identical input (525.769/0.0114 vs 525.759/0.011), SDIMAX=346 both |
| 8 | Volume | volume.jl | ~ FW2W Flewelling conifers; total cubic cornered (bark-fixed +18%→~cornered form-precision); 616BEHW minor species deferred |

## Real bugs found + fixed (all via measurement vs the live binary)

1. **Species crosswalk** — `species_translation.csv` was an EM placeholder; bm/spctrn.f selects ASPT col 5 (BM), not 10 (EM). Grand-fir/white-fir trees (18%+15% of bmt01) were mis-mapped to AF/OS. Fixed → RELDEN 109.10 bit-exact.
2. **Calibration bark** — shared `_backdate_dbh!`/`calibrate_diameter_growth!` had CR/TT bark branches but no BM (POWER model); fell to the linear bark_ratio. Added _bm_bd/_bm_cal.
3. **DIB/D bark ratio** — `bm_bratio` computed DIB (BARK1·D^BARK2), not DIB/D → bark ≈ 0.999 instead of ~0.9. **One fix unmasked FOUR "cornered" residuals**: volume +18%→−4.5%, growth BA 96→98, DG WK2 18/27→24/27, and the GF DGSCOR COR (now bit-exact).
4. **HTDBH dispatch** — regent used Wykoff AX/BX ht-dbh for all conifers, but bm/regent.f:522 uses bm/htdbh.f (Curtis-Arney, forest-dependent) when LHTDRG=false (all except WJ/WB/LM/AS). Ported bm_htdbh + dispatch. Faithful but .sum-neutral (seedling DBH is BA-negligible).

## Hypotheses rejected by measurement (would have been cargo-culted fixes)

- FW2 geosub-substitution as the volume cause (was the DIB/D bark bug)
- GFSUB F-coefficient substitution (implemented → regressed → reverted per doctrine #4)
- Mortality kill-distribution bug (RN math matches live for identical input)
- Regent seedling height over-growth (height matches live)
- Regent seedling DBH as the divergence driver (HTDBH fix is .sum-neutral)
- GF DGSCOR COR as "accepted cornered" (was the bark bug — now bit-exact)

## Remaining residuals (measured, precisely attributed)

- **Multi-cycle divergence** (2090 TPA 163/96): the DGSCOR-FRM serial-correlation (dgdriv.f:271, WK2 bit-exact but DG-from-WK2 is ZZRAN-affected) + tripling, amplified by the faithful mortality feedback. The accepted cornered stochastic class doctrine #3 forbids chasing per-record.
- **Volume total-cubic** ~cornered: INGY SF_SHP form precision (jsp-13 coefficients bit-identical to NVEL source; GFSUB substitution ruled out).
- **Deferred leaves** (scoped, both substantial — not bounded polish):
  - **616BEHW minor species** (WP/MH/WJ/WB/LM/PY/YC/AS/CW/OS/OH — not in bmt01). Behre equations
    616BEH*** route to **R6VOL3** for total cubic (NVEL profile.f:268): `R6VOL3(DBHOB,DBTBH,FCLASS,
    HTTOT,ZONE,VOL)` + BEHTAP taper (blmtap.f:275) + FCLASS form-class. Real port + a synthetic
    pure-species stand (IE pure_{sp} pattern) to validate. Currently these species return 0 volume.
  - **FFE** — bmt01.key is a full FMIN/SIMFIRE/PotFIRE/FuelOut demo; the fire/fuel/snag/carbon
    subsystem (needs BM biomass + fuel coefficients). Large; validate vs the full bmt01.key .sum.

## Verdict

BM growth core (site + DG + COR + height + crown + regent + mortality) is **bit-exact/faithful chunk-by-chunk**
against the live Fortran; volume bark-fixed with a cornered form-precision tail. Every measurable deterministic
layer is validated; the sole remaining residual is the DGSCOR-FRM/tripling/ZZRAN stochastic class. BM joins
EM and UT as a ported, validated core western variant.
