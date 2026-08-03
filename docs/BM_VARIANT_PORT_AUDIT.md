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
| 8 | Volume | volume.jl | ✓ FW2W Flewelling conifers (total cubic cornered, bark-fixed); **616BEHW COMPLETE** — total cubic (bm_r6vol3 + bm_formcl) + merch cubic + board (bm_r6dibs + bm_r6vol1) all bit-exact vs live (per-tree fort.84 6/6; pure-WP .sum MCuFt 1068/1068 + BdFt 5341/5341 exact, TCuFt ±1) |

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
  - **616BEHW minor species** (WP/MH/WJ/WB/LM/PY/YC/AS/CW/OS/OH — not in bmt01) — ✅ **COMPLETE**.
    **MERCH + BOARD LAYER DONE**: ported `bm_r6dibs` (r6vol3-companion log-bucking, ZONE-1/A=0.62/
    total-height path, Behre taper bucked to MTOPP=TOPD·bark) + `bm_r6vol1` (per-log Scribner from the
    132-entry IFTR table + merch-cubic butt-log/Smalian; INTL14 International skipped = VOL(10) not in
    .sum) + accumulation VOL(2)=Σ ANINT(Scribner) [BdFt], VOL(4)=Σ round(merch·10)/10 [MCuFt]. Gated:
    merch/board=0 for h≤17.3 (R6VOL returns cylinder, skips R6DIBS/R6VOL1 — the bug that first made
    merch>total). Wired MCF=VOL(4) (D≥DBHMIN), BdFt=VOL(2) (D≥BFMIND, METHB≠9). **Validated bit-exact**:
    live-fvsvol.f instrument-replay (fort.84) 6/6 trees VOL(2)+VOL(4) exact; pure-WP .sum 1990
    MCuFt 1068/1068 + BdFt 5341/5341 exact (TCuFt ±1); projected cornered = growth-tail (= total cubic).
    ★ Multi-species hardening (distinct bark + form-class paths): pure-MH (sp5, power bark, FIAJSP FC 264)
    = ALL THREE bit-exact 1507/1001/5048; pure-WB (sp11, constant bark 0.969, default FC 80) = BdFt 5479
    exact, MCuFt ±1, TCuFt ±4 (INGY cornered). Confirms the leaf across bark groups + in-list/default FC.
    ── (history) **TOTAL CUBIC LAYER**: ported `bm_r6vol3` (r6vol3.f Behre Smalian taper
    DR=HRATIO/(0.62·HRATIO+0.38), D17=FCLASS/100·DBHOB, H17=17.3, all 3 branches: DBHIB<TOPD cylinder /
    D17<TOPD two-log / full taper) + `bm_formcl` (formclas.f FORMCL_BM form-class lookup: FIAJSP binary
    search, IFCDBH=(D−1)/10+1, 4 forest tables MALH/OCHO/UMAT/WLWH by IFORST=KODFOR%100) + the R6VOL
    short-tree guard (TTH≤17.3→cylinder) + glue (spec=VEQNNC(8:10), DBTBH=D·(1−bm_bratio) per fvsvol.f:153).
    **Validated**: bm_r6vol3 **bit-exact vs live R6VOL3 (96/96** grid pts, all 3 branches, drv_r6vol3.f);
    bm_formcl **bit-exact vs live FORMCL_BM (120/120**, drv_formcl.f); synthetic pure-WP stand end-to-end
    **TCuFt jl 1556 vs live 1557 (±1**, INGY-cornered class), growth cols bit-exact.
    **REMAINING sub-leaf** — merch cubic + board foot (data flow now fully traced, r6vol.f):
    live 616BEHW yields MCuFt/BdFt (pure-WP 1990: 1068/5341) via the full R6VOL path:
    - **R6DIBS**(IAPZ=ZONE,DBHOB,BTR,FCLASS,MTOPP,TLH,TH, → XLOGS,LOGDIA,SL,XL,A): 297-line log-bucking.
      For our case (ZONE 1, HTTYPE='F' ⇒ TH>0/TLH=0) the live path is label 70→80→130 (16.3-ft logs,
      Behre taper DR=HR/(A·HR+B), MTOPP top). Fills NOLOGP(=XLOGS), LOGDIA(21,3) small-end diams,
      SL/XL(20) scaling+actual lengths. ⚠ heavy computed-GOTO + an IRET-dispatched label-1000 taper-A
      setup (unread) + ZONE-2 (32-ft) paths (labels 200-500, not needed for 616).
    - **R6VOL1**(ZONE,DBHOB,FCLASS,NOLOGP,LOGDIA, → LOGVOL,INTBF): 90 lines. Per-log Scribner board
      LOGVOL(1,·), merch cubic LOGVOL(4,·), International INTBF(·).
    - **Accumulate** (r6vol.f:170-186): VOL(2)=Σ ANINT(LOGVOL(1,I)) [Scribner board=.sum BdFt col],
      VOL(4)=Σ round(LOGVOL(4,I)·10)/10 [merch cubic; MCF=VOL(4)+VOL(7), VOL(7)=0], VOL(10)=Σ INTBF.
    - **Wire**: compute_volumes_bm! BEHW branch sets merch_cuft_vol=VOL(4), bdft_vol=VOL(2) (currently 0).
    ~390 NVEL lines + accumulation. bm_r6vol3/bm_formcl already supply the total-cubic core.
    Validate: standalone R6DIBS/R6VOL1 drivers (as with drv_r6vol3/drv_formcl) then pure-WP end-to-end
    MCuFt/BdFt vs FVSbm_clean. Low-impact (species absent from bmt01), substantial — a fresh-session chunk.
  - **FFE** — bmt01.key is a full FMIN/SIMFIRE/PotFIRE/FuelOut demo; the fire/fuel/snag/carbon
    subsystem (needs BM biomass + fuel coefficients). Large; validate vs the full bmt01.key .sum.

## Verdict

BM growth core (site + DG + COR + height + crown + regent + mortality) is **bit-exact/faithful chunk-by-chunk**
against the live Fortran; volume bark-fixed with a cornered form-precision tail. Every measurable deterministic
layer is validated; the sole remaining residual is the DGSCOR-FRM/tripling/ZZRAN stochastic class. BM joins
EM and UT as a ported, validated core western variant.
