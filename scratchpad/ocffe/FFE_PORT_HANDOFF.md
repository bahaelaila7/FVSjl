# OC/OP FFE port — turnkey handoff (2026-08-20)

USER "do all user-gated" ([[feedback-do-all-user-gated]]). OC/OP FFE is TRACTABLE (live oracle exists) —
unlike WWPB PPE (oracle-blocked). Like NC FFE, an unported extension: the ORGANON variants ported
growth+vol only, so their coef.species lacks the FFE per-species columns jl's shared FFE reads.

## ORACLE TARGET (captured)
scratchpad/ocffe/ocffe_oracle.sum — the oct01 stand-4 "FFE" .sum (FVSoc_clean, full oct01.key via
--keywordfile=; stand 4 = the FFE-TEST stand: FMIn/SNAGINIT/FLAMEADJ/SIMFIRE 2003/PotFIRE/FuelOut).
Fire effects VISIBLE: 1993 residual BA 162 vs 184 (SNAGINIT+THINDBH), and the 2003 SIMFIRE kills
TPA 261→72. The .out carries ALL FUELS REPORT + POTENTIAL FIRE REPORT + FUEL CONSUMPTION + SNAGINIT/
SNAGOUT tables — the full validation surface.
Run: /workspace/.ocwork/FVSoc_clean --keywordfile=oct01_full.key (from a dir with oct01_full.tre);
FFE stand = the -999 "FFE" block. FVSoc links FFE (120 fmcba/fmcros refs) ⇒ live oracle.

## GAP (measured)
jl crashes in fmcba!:27 on `length(coef_col(coef,:dbh_min))`. OC/OP coef.species has 4 cols (ORGANON
growth only); IE has 41. jl's shared FFE reads: :dbh_min (crown-biomass merch gate + MAXSP), :v2t
(crown-biomass V2T, rescaled /2000 at fmvinit.f:1094), the species→crown-biomass-group map (SPILS =
ISPMAP for CR, SO_ISPMAP for SO, :ls_spi for LS), + fuel-model & snag coeffs.

## COEFFICIENT SOURCES (OC buildDir)
  • fmvinit.f — V2T crown-biomass per species (DATA/init) + fuel/snag init.
  • fmcrow.f / fmcrowe.f / fmcroww.f — crown-biomass models + the ISPMAP (fmcrow.f:163 SPILS=ISPMAP(sp)).
  • dbh_min — OC has sp_dbh_min in CONTROL (=7, BLM merch) usable as the :dbh_min FFE merch gate.

## PORT PLAN (chunk by chunk, doctrine bit-exact-or-cornered)
  1. Extract OC/OP V2T(NSP) + ISPMAP + dbh_min → new data/oregoncoast (+ olympic) FFE columns, wire into
     coef.species. Unblock fmcba!. Validate CCF/CANOPY COVER first (report-only, fmcba→fmcbcoe).
  2. Fuel loadings (surface fuel by FFE forest type + snag pools) → ALL FUELS REPORT bit-exact.
  3. Fire behavior (SIMFIRE/PotFIRE → fmburn/fmcfir) → POTENTIAL FIRE + the 2003 fire-mortality .sum row.
  4. Snag dynamics (SNAGINIT/SNAGOUT) → the standing-wood tables.
  Mind the westside-FFE bark_intercept crash [[fvsjl-westside-ffe-bark-intercept-crash]] (POWER-bark ports
  need :bark_intercept). Same crown-biomass/CCF-first order that landed CA/EC/WS/SO FFE.
