# Top-level test entry. Run with: julia --project -e 'using Pkg; Pkg.test()'
#
# Tests are organized as:
#   unit/        — pure-kernel + struct tests (fast, no oracle)
#   integration/ — full-run parity vs Oracle A (snt01/snt02/sndb/sntest)
#
# As each chunk lands, its tests are added here. Integration tests are guarded so
# the suite stays green while the engine is still being built (they activate once
# the corresponding chunk is done).

using Test
using FVSjl

@testset "FVSjl" begin
    include("unit/test_core.jl")           # C0: state, rng, units, variant
    include("unit/test_species.jl")        # C2: SN species tables + blkdat defaults
    include("unit/test_crown_width.jl")    # C2: CSV-driven crown-width library
    include("unit/test_forest_type.jl")    # C3: FORTYP/STKVAL FIA forest-type
    include("unit/test_diameter_growth.jl")# C3: DGF/DGCONS equation core
    include("unit/test_growth.jl")         # C3: DG calibration + HTGF + cycle loop
    include("unit/test_mortality.jl")      # C4: MORTS density (Pretzsch) + SDICAL
    include("unit/test_root_disease.jl")   # WRD Chunk −1: RDIN reader + rdinit defaults + inert seam vs live FVSkt
    include("unit/test_dfb.jl")            # DFB: DFBIND/DFBDBH/DFBER/DFBPRB + DFBRAN/BACHLO/DFBMOD/DFBMRT bit-exact vs relinked FVSie_dfb g16 + gated DFBDRV mortality seam
    include("unit/test_dftm.jl")           # DFTM chunk 0: TMRANN RNG + TMOTPR + dftmin.f keyword reader bit-exact vs pristine dftm/*.f + INERT seam
    include("unit/test_wpbr.jl")           # WPBR chunk 0: BRANN RNG + brin.f keyword reader + BRINIT defaults bit-exact vs pristine wpbr/*.f + INERT seam
    include("unit/test_wwpb.jl")           # WWPB beachhead: BMRANN RNG + bmin.f BMIN output-block reader + BMINIT defaults bit-exact vs pristine wwpb/*.f + INERT seam (PPE outbreak driver absent)
    include("unit/test_ppe_landscape.jl")  # PPE: reconstructed area-weighted landscape harness (PPEXCM PTSTV1) — behavior-faithful (source absent), self-consistent area-weighting over validated per-stand run_keyfile + INERT seam
    include("unit/test_ppe_composite_oracle.jl") # PPE CMADDS/CMPRT2 composite-yield aggregation — oracle-baked (historical FVSppe EC 4-stand run): _ppe_aggregate reproduces the COMPOSITE YIELD STATISTICS table bit-exact (area-weighted TREES/vol, area·period-weighted ACC/MOR)
    include("unit/test_climate_spcalib.jl") # Climate-FVS SPCALIB first-cycle presence-calibration (clmorts.f chunk C): PP ViabMort series bit-exact vs oracle FVS_Climate (present low-viability species no longer over-die early)
    include("unit/test_climate_dbs.jl")     # FVS_Climate DBS table (dbsclsum.f, CLIMREDB): climate_report + write_dbs_climate! — Viability/ViabMort bit-exact vs FVSie_clean, BA/TPA/etc cornered on the OLDRN straddle
    include("unit/test_dm_dbs.jl")          # FVS_DM_Stnd_Sum/Spp_Sum DBS tables (dbsmis.f, MISRPTS): mistletoe_report + write_dbs_dm_* — cyc0 bit-exact vs FVSie_clean .out (all columns incl DM mortality)
    include("unit/test_bm_dbs.jl")          # FVS_BM_Main DBS table (dbsbmmain.f, PPBMMAIN): wwpb_main_report + write_dbs_bm_main! — aggregation bit-exact vs bmout.f gfortran-16 golden (WWPB reconstruction bar; no live oracle)
    include("unit/test_rd_dbs.jl")          # FVS_RD_Sum DBS table (dbsrd.f DBSRD1, RDSUM+RRDOUT): rd_sum_report (rdpr.f aggregation) — LIVE-oracle vs FVSkt_clean, ALL 20 cols (Ave_Pct_Root_Inf=PRINF + New_Inf_Prp_Ins/Exp/Tot=CORINF/EXPINF) BIT-EXACT at inventory+first projected cycle, cornered later (#206 WRD growth+spread straddle)
    include("unit/test_rd_det_dbs.jl")      # FVS_RD_Det DBS table (dbsrd.f DBSRD2, RDDETAIL+RRDOUT): rd_det_report (rddout.f RDPSRT→PCTILE→RDDST per-species percentile DBH) — LIVE-oracle vs FVSkt_clean, 1990 inventory row BIT-EXACT (102/102) + aggregates cornered (#206 WRD growth straddle; DBH-percentile cols amplify it)
    include("unit/test_canprofile_dbs.jl")  # FVS_CanProfile DBS table (dbsfmcanpr.f, CANFPROF): canopy_crfill + write_dbs_canprofile! — cyc0 profile bit-exact vs FVScr_clean (PRE-growth collection), multi-cycle cornered
    include("unit/test_strclass_dbs.jl")    # FVS_StrClass DBS table (dbsstrclass.f, STRCLSDB): structure_report + write_dbs_strclass! — cyc0 all-43-col bit-exact vs FVSoc_clean (+OC oc_cwcalc _ss_strata fix), later cycles cornered
    include("unit/test_calibstats_dbs.jl")  # FVS_CalibStats DBS table (dbscalib.f, CALBSTDB): captured DGSCOR calib stats + write_dbs_calibstats! — all 9 cols bit-exact vs FVScr_clean (WF/ES)
    include("unit/test_oc_treelist_dg.jl")  # FVS_TreeList DG/HtG for OC: ORGANON clobbered diam_growth/ht_growth (0 on projected cycles); shared apply-loop now skips OC + keeps the increment (.sum-inert)
    include("unit/test_op_treelist_dg.jl")  # FVS_TreeList DG/HtG for OP: sibling ORGANON (cooperating driver) — confirms OP has NO OC-style clobber, live-tree DG/HtG bit-exact vs FVSop_clean
    include("unit/test_oc_snag_recruit.jl") # OC FFE snag pool: mortality!(::OregonCoast) was a no-op ⇒ ongoing ORGANON mortality never booked as snags; now books t.mort_pa (1998 SnagSum bit-exact vs FVSoc_clean 92.99)
    include("unit/test_ws_crwidth.jl")      # WS FVS_TreeList CrWidth via R5CRWD (ws_r5crwd) — R5 skips forest-BF ⇒ per-tree bit-exact (29/29 vs FVSws_clean), the 8th variant wired past the eastern 0.5 default
    include("unit/test_ca_crwidth.jl")      # CA FVS_TreeList CrWidth via ca_cwcalc + folded forest-610 BF (SP/LP/PP) — 29/29 bit-exact vs FVSca_clean, INERT on cat01_ffe .sum; 9th variant wired
    include("unit/test_bm_crwidth.jl")     # BM FVS_TreeList CrWidth via bm_cwcalc→cr_cwcalc + folded forest-614 BF (forest_bf flag; FFE path BF-free) — 29/29 vs FVSbm_clean, 10th variant
    include("unit/test_so_crwidth.jl")     # SO FVS_TreeList CrWidth: so_cwcalc forest-601 BF (present) + Hopkins-index fix (so_grinit! lat/lon 42/121 default) — 33/33 vs FVSso_clean, 11th variant
    include("unit/test_nc_crwidth.jl")     # NC/Klamath FVS_TreeList CrWidth via R5CRWD reuse (nc_r5crwd→ws_r5crwd, NC→WS FIA map) — forest 505 R5, 29/29 vs FVSnc_clean, 12th variant
    include("unit/test_ie_crwidth.jl")     # IE FVS_TreeList CrWidth via national cwcalc.f IEMAP dispatch (reuses EM forms + 8 IE codes) — 4399/4399 vs FVSie_clean, 13th variant
    include("unit/test_kt_crwidth.jl")     # KT FVS_TreeList CrWidth via national cwcalc.f KTMAP (= IEMAP[1:11]) → ie_cwcalc — 6132/6132 vs FVSkt_clean, 14th variant
    include("unit/test_ci_crwidth.jl")     # CI FVS_TreeList CrWidth via national cwcalc.f CIMAP (+4 codes 26305/01905/06405/47502) — 4000/4000 vs FVSci_clean, 15th variant
    include("unit/test_tt_crwidth.jl")     # TT FVS_TreeList CrWidth via national cwcalc.f TTMAP (+4 codes 20205/09305/10805/31206) — 4000/4000 vs FVStt_clean, 16th variant
    include("unit/test_ut_crwidth.jl")     # UT FVS_TreeList CrWidth via national cwcalc.f UTMAP (+3 codes 01505/81402/10201) — 4000/4000 vs FVSut_clean, 17th variant
    include("unit/test_ak_crwidth.jl")     # AK FVS_TreeList CrWidth via national cwcalc.f AKMAP (+ '08' form _cw08 + R10 codes) — 3200/3200 vs FVSak_clean, 18th variant
    include("unit/test_bc_crwidth.jl")     # BC FVS_TreeList CrWidth via national cwcalc.f BCMAP (0 new codes; metric) — log-forms bit-exact vs FVSbc_clean, 19th (final western) variant
    include("unit/test_lpmpb.jl")          # LPMPB: COLDBH/COLIND/COLMOD/COLMRT/MPBER Cole rate-of-loss core + MPRANN seed 55329 bit-exact vs relinked FVSie_lpmpb g16 + gated MPBCUP mortality seam
    include("unit/test_lpopdy_chain.jl")   # LPMPB LPOPDY: BETIN/GARBEL/SURFCE/MPBMOD epidemic chain bit-exact vs FVSie_lpmpb (golden fixtures)
    include("unit/test_lpmpb_damage.jl")   # LPMPB INVMORT: treelist MPB damage-code GREINF, cycle-1 mortality delta vs FVSie_lpmpb
    include("unit/test_wsbwe.jl")          # WSBWE beachhead: BWERAN RNG seed 55329 bit-exact vs pristine wsbwe/bweran.f + WSBW keyword reader (keywds.f opt 8) + INERT seam
    include("unit/test_wsbwe_gendefol.jl") # WSBWE GENDEFOL/BUDLITE: ported bwelit.f core (wsbwe_bwelit!) dump-replay bit-exact vs FVSem_wsbwe on the synthetic-weather harness (1990 pulse / 1991 tail / 1992 crash)
    include("unit/test_cover.jl")          # COVER beachhead: CVCW crown-area (CRAREA=Σ CRWDTH²·PROB·0.785398) dump-replay bit-exact vs FVSem_g16 (report-only extension)
    include("unit/test_cvbcal.jl")         # COVER shrub CALIBRATION (cvbcal.f): BHTCF/BPCCF by-layer (SHRBLAYR) + by-species (SHRUBHT/SHRUBPC) correction factors + apply, Float32-hex dump-replay bit-exact vs FVSem_g16
    include("unit/test_ppe_sort.jl")       # PPE PPBASE C11SRT/C26SRT/CH8SRT character index QuickerSort (master stand ordering) — bit-exact vs gfortran-16 golden (recovered PPE source)
    include("unit/test_ppe_search.jl")     # PPE PPBASE C26BSR/CH8BSR/SPBSRX keyed binary search — bit-exact vs gfortran-16 golden
    include("unit/test_ppe_hxindx.jl")     # PPE PPBASE HXINDX hexagonal neighbor indexing — bit-exact vs gfortran-16 golden (Float32-sqrt boundaries)
    include("unit/test_ppe_add1.jl")       # PPE PPBASE ADD1 internal-stand-number digit increment — bit-exact vs gfortran-16 golden
    include("unit/test_ontario.jl")        # ON (Ontario) beachhead: Penner large-tree DGF (annual diameter-increment, on/dgf.f) dump-replay bit-exact vs FVSon_g16
    include("unit/test_ontario_dgf_wired.jl")# ON: engine-WIRED dgf! (stand-context assembly + on_bratio + leftover-BARK DDS) WK2 bit-exact vs production FVSon_g16 on ont01 cyc0
    include("unit/test_ontario_growth_wired.jl")# ON: runnable growth — species table+translation (reader) + coefficients(::Ontario) standalone + DDS→DG on_bratio branch (d_ib/DDS/WKI) bit-exact vs FVSon_wkidump
    include("unit/test_ontario_htg.jl")    # ON: large-tree height growth (htgf.f/htont.f Penner diameter-height) — shipped height_growth!(::Ontario) per-tree HTG + HTONT bit-exact vs instrumented FVSon_g16 on ont01
    include("unit/test_ontario_volume.jl") # ON: per-tree volume (vols.f/varvol.f METHC=8 + volont.f ZAK/HONER + Mowraski cull) GTV/GMV/NMV + on_tree_age dump-replay bit-exact vs FVSon_g16
    include("unit/test_ontario_sum_classification.jl")# ON: cyc0 .sum row bit-exact vs FVSon_g16 through FORTYP/size/stock (metric stkval; row tail 125 11)
    include("unit/test_ontario_allspecies_dgf.jl")   # ON full-port: Penner large-tree DGF (on_penner_dds) bit-exact across ALL 72 species vs FVSon_g16 (72-tree ont_all stand)
    include("unit/test_ontario_allspecies_htg.jl")   # ON full-port: Penner diameter-height (on_htont) bit-exact across ALL 72 species vs FVSon_htdump
    include("unit/test_ontario_allspecies_vol.jl")   # ON full-port: cubic volume (on_zakvol/on_honer) bit-exact across ALL 72 species vs FVSon_voldump
    include("unit/test_ontario_allspecies_cw.jl")    # ON full-port: open-grown crown width (on_open_crown_width) bit-exact across ALL 72 species vs FVSon_cwdump (fixes CCF)
    include("unit/test_ontario_allspecies_mow.jl")   # ON full-port: Mowraski net-merch cull (on_mowraski) bit-exact across ALL 72 species vs FVSon_mowdump
    include("unit/test_ontario_missing_height.jl")  # ON full-port: missing/broken-top height dub crash fix (_on_htdbh_height) + bit-exact volume vs FVSon_g16
    include("unit/test_ontario_database_reader.jl") # ON DATABASE reader: metric cm->in/m->ft conversion (was BC-only gate; ON DB was ~2.5x off)
    include("unit/test_dvee_volume.jl")    # D35: R9 Gevorkiantz '900DVEE' volume vs live
    include("unit/test_ie_estock.jl")      # #143: IE AUTOES ESTOCK P(stocking) vs live FVSie
    include("unit/test_oc_organon_setup.jl")# C2: OC ORGANON PREPARE calibration (TMPCAL) vs live FVSoc
    include("unit/test_op_organon_nwo.jl") # OP ORGANON NWO engine (DG/HG/CR/MORT + PREPARE ACALIB) vs live FVSop
    include("unit/test_op_native_growth.jl") # OP FVS-native large-tree DGF/HTGF (non-ORGANON species) vs live FVSop
    include("unit/test_op_site_and_volume.jl") # OP site-index fan + BLM Behre volume (chunk 2) vs live FVSop
    include("unit/test_op_multicycle_sum.jl") # OP end-to-end multi-cycle .sum (growth+total-cubic bit-exact, merch cornered) vs relinked FVSop_clean
    include("unit/test_oc_regent_smtree.jl") # OC REGENT small-tree height (smhtgf.f) — seedling +7ft bug fixed, vs live FVSoc_clean
    include("integration/test_opt01_cyc0.jl") # OP chunk 3: opt01 cyc0 .sum (density+volume) bit-exact vs live FVSop
    include("integration/test_svs_chunk0.jl") # SVS chunk 0: cyc0 .svs object list bit-exact vs live FVSkt
    include("integration/test_svs_multicycle.jl") # SVS multi-cycle: end-of-projection + begin-cycle pictures vs live FVSkt
    include("integration/test_svs_snag.jl") # SVS mortality→snag: standing-dead aging + display bit-exact vs live FVSkt
    include("integration/test_treedata.jl")# C1: .tre parser vs Oracle A
    include("integration/test_keyword.jl") # C1: keyword lexer vs Oracle A
    include("integration/test_io_formats.jl")# C1b: CSV/format-agnostic round-trips
    include("integration/test_dbs_summary.jl")# C6: DBS FVS_Summary SQLite table vs Fortran
    include("integration/test_dbs_treelist.jl")# C6: DBS FVS_TreeList per-tree table
    include("integration/test_dbs_compute.jl") # C6: DBS FVS_Compute event-monitor vars table
    include("integration/test_dbs_invref.jl")  # C6: DBS FVS_InvReference per-species reference table
    include("integration/test_dbs_cutlist.jl") # C6: DBS FVS_CutList removed-record table
    include("integration/test_init.jl")    # C2: keyword dispatch + tree loading
    include("integration/test_snt01.jl")
include("integration/test_net01.jl")   # C5: .sum cycle-0 bit-exact + cycle-1 tracking
include("integration/test_cst01.jl")   # CS: cst01 cycle-0 stand columns bit-exact (GROSPC<1 path)
include("integration/test_lst01.jl")   # LS: lst01 cycle-0 stand columns bit-exact (all 6, vs live FVSls)
include("integration/test_lst01_ffe.jl")  # LS: FFE fire behavior (fmcfmd model 10 + fmmois) + fire mortality vs live FVSls
include("integration/test_lst01_estab.jl") # LS: BARE-plant establishment stand BIT-EXACT vs live FVSls (RAN-window fix)
include("integration/test_lst01_fire_sprout.jl") # LS: post-fire stump sprouting (fmkill.f:80 fire-kill→ESTUMP→ESUCKR) vs live FVSls
include("integration/test_allspecies.jl")# CS+NE+SN+LS: all-species coverage vs live (per-species coefficient rows; caught the BW crown-width gap)
include("integration/test_canonical_multistand.jl")# SN/NE/CS/LS: full reference multi-stand scenarios (thin/shelterwood/fire/plant) vs live
include("integration/test_ls_sitesweep.jl")# LS: site-productivity sweep — deterministic growth bit-exact across site indices vs live
include("integration/test_r8clark_special.jl")# D7: R8 Clark COEFFSO%DIB17 cypress/green-ash merch volume vs live FVSsn
include("integration/test_simfire_schedule.jl")# D9: SIMFIRE date/cycle default + multiple-fire scheduling vs live FVSsn
include("integration/test_compress_tripling.jl")# COMPRESS still triples its own cycle (merch volume) vs live FVSsn
include("integration/test_growth_fint.jl")# D2: GROWTH FINT!=5 measurement period (first-cycle serial-corr old) vs live FVSsn
include("integration/test_r8_intl_board.jl")# D11: R8 International 1/4in board feet for GW-JF/Ouachita/Ozark forests vs live FVSsn
    include("integration/test_multicycle.jl")# C3/C4/C5: multi-cycle regression vs oracle golden
    include("integration/test_cuts_coverage.jl")# C3: CUTS keyword coverage + gap tracker (decision flow)
    include("integration/test_regen_coverage.jl")# C4: regen/ESTAB coverage + gap tracker (bare stands)
    include("integration/test_fortbragg_coverage.jl")# C5: Fort Bragg (forest 701) KODFOR remap → nonzero volume
    include("integration/test_multistand.jl")# C2/C8: multi-stand driver (each_stand) — TREFMT persist + default INTREE
    include("integration/test_multistand_sum.jl")# C5: multi-stand .sum parity vs Fortran (all 5 stands, state-carry guard)
    include("integration/test_hcor_calib.jl")# C3: REGENT small-tree HCOR calibration vs Fortran
    include("integration/test_multipliers.jl")# C3: growth/mortality keyword multipliers (BAIMULT/HTGMULT/MORTMULT) vs Fortran
    include("integration/test_treeszcp.jl")   # C3: per-species size cap (TREESZCP/SIZCAP) vs Fortran
    include("integration/test_fix_scalers.jl") # C3: FIXDG/FIXHTG one-shot growth scalers vs Fortran
    include("integration/test_htgstp.jl")      # C4: HTGSTOP/TOPKILL top-damage events vs Fortran
    include("integration/test_fixmort.jl")     # C4: FIXMORT forced-mortality override vs Fortran
    include("integration/test_crnmult.jl")     # C3: CRNMULT crown-ratio-change multiplier vs Fortran
    include("integration/test_spgroup.jl")     # C2: SPGROUP species groups (ISPCC<0 refs) vs Fortran
    include("integration/test_tripling.jl")    # C3: NOTRIPLE/NUMTRIP tripling control (ICL4) vs Fortran
    include("integration/test_timeint.jl")     # C3: TIMEINT cycle calendar (period scaling) vs Fortran
    include("integration/test_compute.jl")     # C4: COMPUTE event-monitor user variables vs Fortran
    include("integration/test_volume_override.jl")# C5: VOLUME merch-standard override (DBHMIN gate) vs Fortran
    include("integration/test_mcdefect.jl")     # C5: MCDEFECT/BFDEFECT defect curves vs Fortran
    include("integration/test_pertree_defect.jl")# C5: per-tree DEFECT (damage-code) input vs Fortran
    include("integration/test_voleqnum.jl")     # C5: VOLEQNUM cubic volume-equation override vs Fortran
    include("integration/test_bfvolume.jl")     # C5: BFVOLUME board-foot override + Region-8 vs Fortran
    include("integration/test_minharv.jl")      # C3: MINHARV minimum-harvest cancel gate vs Fortran
    include("integration/test_spleave.jl")      # C3: SPLEAVE/LEAVESP leave-species-during-thin vs Fortran
    include("integration/test_fertiliz.jl")     # C3: FERTILIZE/FFERT fertilizer growth response vs Fortran
    include("integration/test_tcondmlt.jl")     # C3: TCONDMLT tree-condition cut weight vs Fortran
    include("integration/test_bc_dmntrd.jl")    # #196: BC NEWSPRED DMNTRD crown-third remap + tripling DMR carry
    include("integration/test_tfixarea.jl")     # C2: TFIXAREA fixed-plot-area expansion vs Fortran
    include("integration/test_cuteff.jl")       # C3: CUTEFF default cutting efficiency vs Fortran
    include("integration/test_managed.jl")      # C3: MANAGED → DGF planted/managed growth term vs Fortran
    include("integration/test_bamax.jl")         # C3: BAMAX → SDImax self-thinning cap vs Fortran
    include("integration/test_sdimax.jl")     # C3/C4: SDIMAX per-species SDImax + PMSDIL/PMSDIU vs Fortran
    include("integration/test_rannseed.jl")   # RNG: RANNSEED reseed of the main stochastic stream vs Fortran
    include("integration/test_compress.jl")    # COMPRESS keyword recognition + scheduling (algorithm = chunk plan)
    include("integration/test_structure_stage.jl") # SSTAGE structural-stage class (1-6) vs Fortran
    include("integration/test_carbon.jl")          # Stand Carbon Report Jenkins live pools vs Fortran
    include("integration/test_fire.jl")            # FFE fire stand (SIMFIRE) end-to-end vs Fortran
    include("integration/test_growth.jl")       # GROWTH keyword recognition + param capture
    include("integration/test_dgstdev.jl")     # DGSCOR: DGSTDEV DGSD bound on stochastic DG variation vs Fortran
    include("integration/test_nocalib.jl")     # NOCALIB disable DG self-calibration (LDGCAL) vs Fortran
    include("integration/test_serlcorr.jl")    # DGSCOR: SERLCORR ARMA(1,1) phi/theta vs Fortran
    include("integration/test_resetage.jl")    # RESETAGE rebase stand age (resage.f) vs Fortran
    include("integration/test_sdicalc.jl")     # SDICALC SDI method (Zeide/Reineke) + thresholds vs Fortran
    include("integration/test_ccadj.jl")       # CCADJ crown-competition adj: recognized .sum-inert no-op (SN)
    include("integration/test_cycleat.jl")      # CYCLEAT extra cycle boundary (non-uniform IY schedule) vs Fortran
    include("integration/test_readcor.jl")      # READCOR/REUSCOR growth-constant corrections (COR2/HCOR2/RCOR2) vs Fortran
    include("integration/test_setsite.jl")      # SETSITE scheduled mid-run site-index change (act 120) vs Fortran
    include("integration/test_nohtdreg.jl")     # NOHTDREG HT-DBH (LHTDRG) calibration control: suppress no-op / invoke warn
    include("integration/test_mortmsb.jl")       # MORTMSB/MSBMRT alternate mature-stand-breakup mortality vs live FVS
    include("integration/test_sprout_table.jl")  # SPROUT per-species/DBH-range multiplier table (esuckr act 450) vs live FVS
    include("integration/test_estab_pccf.jl")     # regen crown ratio uses stand CCF (PCCF), not 0 — vs live FVS
    include("integration/test_estab_rng_d10.jl")  # D10: establishment :estab RNG stream (RAN window + WK6/NTALLY draws) vs live FVS
    include("integration/test_thinprsc_fragment_d14.jl")  # D14: THINPRSC residual≤0.0005 whole-tree deletion (cuts.f:1632) vs live FVS
    include("integration/test_fire_rng_restore_d15.jl")  # D15: FMEFF RANNGET/RANNPUT RNG save-restore (post-fire growth) vs live FVS
    include("integration/test_mcfdln.jl")       # C5: MCFDLN/BFFDLN form-model coefs (no Fortran oracle — FPE)
    include("unit/test_sprout.jl")              # ESUCKR-B: NSPREC/SPRTHT/ESSPRT sprout sub-routines + Wykoff DBH + cut-log
    include("integration/test_sprout_regen.jl") # ESUCKR-C/D: stump-sprout regen generation loop vs live Fortran
    include("integration/test_thindbh_cycledate.jl") # cuts: blank-date THINDBH = cycle-number date (initre.f:1189)
    include("unit/test_fire_biomass.jl")        # FFE-F1/F2/F3: biomass, crown fuels, surface fuels
    include("unit/test_fire_effects.jl")        # FFE-F6: fire-caused mortality (FMEFF/FMBRKT)
    include("unit/test_rothermel.jl")           # FFE-F5: Rothermel surface fire behavior (FMFINT)
    include("unit/test_fmburn.jl")              # FFE-F5b: fire event driver (FMBURN/FMEFF) → kill TPA
    include("unit/test_carbon.jl")              # FFE-F8: standing live-tree carbon pools (FMCRBOUT)
    include("unit/test_fuel_decay.jl")          # FFE-F3: surface-fuel decay (FMCWD)
    include("unit/test_crown_lift.jl")          # FFE-F3: crown-lift rate X for down-wood additions (FMSDIT)
    include("unit/test_estab_specmult_htadj.jl")  # ESTAB-packet SPECMULT (XESMLT occ mult) + HTADJ (height adj) — parsing + occ-scaling + HTADJ +5.0 vs FVSci_clean
    include("unit/test_estab_minplots.jl")        # ESTAB-packet MINPLOTS (esin.f MINREP → DUPNPT plot replication) — parse+clamp + bit-exact DUPNPT-scaled ESAVE draw chain vs live FVSie
    include("unit/test_estab_mechprep.jl")        # ESTAB-packet MECHPREP/BURNPREP site prep (esetpr.f + estab.f:382-399 WK6 IPPREP sampler) — kernel + live IPPREP vs FVSie_estabdump
    include("unit/test_snag.jl")                # FFE-F7: snag falldown + decay dynamics (FMSFALL)
    include("unit/test_consumption.jl")         # FFE-F7/F8: fire fuel consumption + carbon release (FMCONS)
    include("unit/test_econ.jl")                # C8: ECON economic-analysis core (eccalc.f)
    include("integration/test_longrun.jl")# C4: COMCUP zero-PROB record deletion (long unthinned run)
    include("integration/test_event_monitor.jl")# C4: event monitor (IF/THEN/ENDIF) evaluator + firing
    include("integration/test_keyword_coverage.jl")# C8: SN keyword-coverage drop-in (37 scenarios) vs live FVSsn + yaml==key
    include("integration/test_translate.jl")# C8: YML/CSV<->.key/.tre conversion tool round-trips (translate_io)
    include("integration/test_parallel.jl") # pillar-3: multi-stand parallel run == serial, bit-identical, all variants
    include("integration/test_allocation.jl") # pillar-2: grow_cycle! per-cycle allocation floor guard (net01 NE)
    include("integration/test_kwcov_variants.jl") # pillar-1: NE/CS/LS keyword-isolation coverage vs live Fortran
    include("integration/test_fia_reader.jl") # FIA "FVS-ready" DATABASE/DSNIN reader: cycle-0 bit-exact vs live
    # include("integration/test_sndb.jl")   # enabled at C6
    # include("integration/test_snt02.jl")  # enabled at C8
end
