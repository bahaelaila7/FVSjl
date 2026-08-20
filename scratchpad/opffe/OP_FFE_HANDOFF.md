# OP (Olympic) FFE port — turnkey handoff (2026-08-20)

OP FFE is a CLEAN MIRROR of the completed OC FFE port ([[fvsjl-oc-organon-blm-volume]], commits 342b1afc→
216bd02a). Same keyfile structure (SNAGINIT 10/11/50/40/2/50, SIMFIRE 2003, FLAMEADJ, PotFIRE), just
forest **708** (BLM Salem→606 Mt Hood for cwcalc BF), STDINFO elev 7.0, and Olympic's **39 NWO ORGANON
species** (MAXSP=39, not OC's 50). Do it chunk-by-chunk exactly as OC.

## HARNESS (ready)
- Oracle: `/workspace/.opwork/FVSop_clean --keywordfile=opt01_full.key` (opt01_full.{key,tre} = the
  tests/FVSop/opt01.{key,tre} copies). FFE stand = the `-999 "FFE"` .sum block.
- Target: `scratchpad/opffe/opffe_oracle.sum` (13 cycles; 2003 SIMFIRE kills TPA 280→75, mort 328).
- jl harness: `scratchpad/opffe/opffe_full.{key,tre}` (SELF-CONTAINED, TREEFMT inlined — mirrors ocffe_full).
- Current jl state: CRASHES at fmcba! with the legible "FFE not ported for Olympic" guard (needs coeffs).
- Relink recipe: `/workspace/.opwork/build_g16.sh` (or a relink_op.sh mirroring relink_oc.sh) for instrumented
  single-.o dumps (fmburn/fmdyn/fmcba/fmmois) — same as the OC measurement vehicle.

## CHUNKS (mirror OC; extraction sources = op/*.f in bin/FVSop_buildDir)
Each is the OP-39 analogue of the committed OC fix:
  1. **crown-biomass + bark + snag species-props** (like 342b1afc): make an `op/ffe_coefficients.jl` +
     `data/olympic/fire_species_props.csv` (39 sp). Sources: op/fmvinit.f (V2T + LEAFLF + DKRCLS + FALLX +
     ALLDWN, snag_decayx=999), op/fmcrow.f ISPMAP, op/blkdat.f ISPSPE (is_sprouting), op/fmcblk.f BIOGRP
     (Jenkins bio_group — NOT ISPMAP, the OC gotcha), dbh_min (OP BLM merch). Add `_OP_FM_BARK_B1` +
     Olympic branch in fire_bark_thickness. Inject into coefficients(::Olympic). Add Olympic to the
     crown_biomass.jl bark/ls_spi dispatch via op_bratio + OP_FFE_ISPMAP.
     PARTIAL DATA already pulled (VERIFY before use — the awk was rough):
       bark B1(39): .047 .048 .046 .041 .039 .027 .045 .022 .081 .036 .028 .068 .072 .035 .063 .063 .081
         .035 .040 .040 .024 .026 .060 .045 .045 .044 .044 .029 .025 .050 .030 .030 .025 .062 .038 .062
         .041 .0?? .044  (sp38 read 0.000 — re-extract; likely a continuation artifact)
       BIOGRP(partial): 3,3,3,3,3,5,3,1,1,5,10,4,4,4,4,4,2,1,1,3,3,20,7,6,6,8,8,…(cut off — re-extract full 39)
       ISPMAP(40 read, want 39): 4,4,4,1,4,18,4,8,20,18,11,15,15,15,13,3,19,7,6,24,5,23,10,17,17,41,17,17,
         16,1,14,11,7,56,57,61,64,0,41,6 (last值 stray — re-extract exactly 39)
  2. **fire_fuel_models.csv** — the standard Anderson-13 (variant-independent; copy from data/pacificnorthwest).
  3. **fmcba fuel loading + op_cwcalc** (like 3c4a1132): extract op/fmcba.f FUINIE/FUINII/FULIVE/FULIVI (39 sp)
     → data/olympic/fire_fuel_covtype_*.csv (NON-reserved names!). op/cwcalc.f OPMAP → op_cwcalc (forest 708→
     606 Mt Hood BF). Check: OP is R6 westside — op/fmcfmd.f likely == wc/pn (FIRE-VPN), NOT ca-CWHR; op/fmmois.f
     likely == the western/NC dry table. VERIFY by diff vs FVSwc/FVSpn/FVSnc buildDir.
  4. **fmcfmd** (like 440549bc): `diff op/fmcfmd.f {wc,pn,nc,ca}/fmcfmd.f` → route Olympic to the matching
     jl selector (wc_select_fuel_models if == wc/pn; else ca/nc). Add Olympic to the select_fuel_models dispatch.
  5. **snag bole = BLM total cubic** (like 58b6db58): add Olympic branch (op_tree_cuft — the OP analogue of
     oc_tree_cuft; OP uses BLMVOL too) to ALL THREE snag paths (ffe_add_snaginit!, ffe_seed_input_snags!,
     mortality.jl). The R8-Clark else → 0 → Jenkins over-book is the same trap.
  6. **fuel moisture** (like 216bd02a): `diff op/fmmois.f nc/fmmois.f` → fm_mois_table(::Olympic) = the matching
     table (probably _FM_MOIS_NC or a new _FM_MOIS_OP). The wet SN default is WRONG (byram ~25% low).
  7. **ffe_on** (like 7d0fcc30): add Olympic to the summary.jl:200 allowlist (its live-fuel CSV will be
     fire_fuel_covtype_live.csv → ffe_fuel_live empty, same as OC).
  8. **is_sprouting/BIOGRP/DKRT** as needed (op/fmcba.f DKRADJ — check OP forest 708 R6 branch + ITYPE habitat).
Validate each chunk by RUNNING FVSop_clean; end-to-end .sum vs opffe_oracle.sum. Expect the same VERDICT as
OC: fire behavior bit-exact, mortality residual = cornered OP small-tree/growth drift (OP multi-cycle is
already bit-exact-or-cornered per 9594d0a7). Gate: multicycle 339/11.

## KEY OC LESSONS TO CARRY
- bio_group/biogrp = Jenkins BIOGRP (1-10), NOT the crown-biomass ISPMAP (the jenkins_biomass BoundsError).
- Name the OP fuel CSVs fire_fuel_covtype_* (NOT the reserved fire_fuel_{dead,live,models}.csv).
- ffe_on gates the ENTIRE per-cycle fuel machinery — add Olympic or the fire samples an empty pool.
- The snag bole must be BLM TOTAL cubic (op/fmsvol.f LMERCH=.FALSE.), else R8-Clark→0→Jenkins ~3× over-book.
- The fuel-moisture table drives byram ~25% — the wet SN default gives a cool fire.
