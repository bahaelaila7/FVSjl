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
    include("integration/test_treedata.jl")# C1: .tre parser vs Oracle A
    include("integration/test_keyword.jl") # C1: keyword lexer vs Oracle A
    include("integration/test_io_formats.jl")# C1b: CSV/format-agnostic round-trips
    include("integration/test_dbs_summary.jl")# C6: DBS FVS_Summary SQLite table vs Fortran
    include("integration/test_dbs_treelist.jl")# C6: DBS FVS_TreeList per-tree table
    include("integration/test_dbs_compute.jl") # C6: DBS FVS_Compute event-monitor vars table
    include("integration/test_dbs_invref.jl")  # C6: DBS FVS_InvReference per-species reference table
    include("integration/test_dbs_cutlist.jl") # C6: DBS FVS_CutList removed-record table
    include("integration/test_init.jl")    # C2: keyword dispatch + tree loading
    include("integration/test_snt01.jl")   # C5: .sum cycle-0 bit-exact + cycle-1 tracking
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
    include("unit/test_snag.jl")                # FFE-F7: snag falldown + decay dynamics (FMSFALL)
    include("unit/test_consumption.jl")         # FFE-F7/F8: fire fuel consumption + carbon release (FMCONS)
    include("unit/test_econ.jl")                # C8: ECON economic-analysis core (eccalc.f)
    include("integration/test_longrun.jl")# C4: COMCUP zero-PROB record deletion (long unthinned run)
    include("integration/test_event_monitor.jl")# C4: event monitor (IF/THEN/ENDIF) evaluator + firing
    # include("integration/test_sndb.jl")   # enabled at C6
    # include("integration/test_snt02.jl")  # enabled at C8
end
