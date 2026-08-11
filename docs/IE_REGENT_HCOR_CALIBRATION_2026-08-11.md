# IE regent small-tree HEIGHT calibration (HCOR) ported (#171) — 2026-08-11

## Problem
On dense real-FIA stands the IE small-tree regent OVER-grew: stand 373781950489998 (subalpine)
.sum ACC +39% (jl 57 vs live 41 @2015→25), BA +14% by 2025, growing to +10% by 2065 — flagged #171.

## Root cause (MEASURED via FVSie_g16, layered instrumentation)
- large-tree DG BIT-EXACT (AF sp9 + LL sp14) ⇒ the excess is the small-tree (<3", ~55% of TPA) regent cohort.
- the jl regent DG FORMULA is faithful (regent.jl matches regent.f:938-985) ⇒ over-prediction is in the regent
  HEIGHT growth (htg ~3.3× high), which feeds the height-implied diameter DK.
- per-subcycle HTGRL=CON+RHLH·ln(H)+RHCCF·RDJ+RHBAL·BAL diverges ONLY in CON: jl 0.497 vs live −0.648.
- RHCON bit-exact (jl=live=0.368). CON=RHCON+HCOR ⇒ jl HCOR=+0.129 (leaked large-tree diameter COR) vs live
  HCOR=−1.016 from the REGENT's OWN calibration.
- The regent HCOR is `htg_cor_init` attenuated by the shared dgdriv.f:188-194 formula in
  calibrate_diameter_growth! (`htg_cor_small = dg_cor_goal + cormlt_h·(htg_cor_init − dg_cor_goal)`). IE never
  set `htg_cor_init` (no IE branch in the regent regression) ⇒ it was 0 ⇒ the diameter COR (dg_cor_goal) leaked
  into the regent height CON. The Southern (ht-curve) and Teton (POTHTG) branches DO set it; IE (NIVAR) did not.

## Fix
`ie_regent_hcor_init!` (src/variants/inlandempire/regent.jl) ports ie/regent.f:1138-1337 for NIVAR species:
for each sub-5" tree with a measured height increment, accumulate predicted EDH=HK−H (HK grown over the
subcycles by the NIVAR model with HCOR=0, using the BACKDATED live+dead stand density) and measured
TERM=HTG·SCALE3 (SCALE3=REGYR/FINTH); CORNEW=Σ(TERM·P)/Σ(EDH·P), trapped [0.0821,12.1825];
htg_cor_init[sp]=ln(CORNEW). Called from calibrate_diameter_growth! (IE-guarded, before the shared attenuation).
Inert when fewer than NCALHT(5) sub-5" trees carry a measured HTG ⇒ htg_cor_init stays 0 (iet01).

## Validation (vs FVSie_clean / FVSie_g16)
- Per-tree calibration intermediates matched vs FVSie_g16 anchor: SNY bit-exact (0.9375); measured HTG + backdated
  H exact; EDH within ~3.5% (jl SNX 4.40 vs live 4.25) ⇒ HCOR_init jl −1.546 vs live −1.512 (minor residual, a
  small density/PCT-backdating precision gap — accepted, the stand-level result meets the bar).
- Stand 373781950489998: BA +14%→+2% (2045 EXACT 120=120; 2065 137 vs 136) — #171 RESOLVED to cornered.
- iet01 BYTE-IDENTICAL (calibration doesn't fire — no measured small-tree HTG).
- 4 other IE stands (larch/juniper/dense-seedling): byte-identical to pre-fix ⇒ calibration inert there, 0 crashes.
- IE-guarded (`s.variant isa InlandEmpire`) ⇒ the other 7 western variants untouched by construction.

## Residual / follow-ups
- The ~3.5% EDH precision gap (HCOR_init −1.546 vs −1.512) — likely a backdated-PCT or density-projection subtlety;
  the stand tracks bit-exact-or-cornered regardless. Investigate only if a stand fails the bar.
- Other IE stands' pre-existing residuals are SEPARATE: larch DG/mortality (31366571010690 −9% BA), juniper
  AUTOES/mortality (39447478010690), and #174 dense-seedling under-growth (calibration doesn't fire on sub-1").
