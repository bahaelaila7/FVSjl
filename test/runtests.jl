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
    include("integration/test_longrun.jl")# C4: COMCUP zero-PROB record deletion (long unthinned run)
    include("integration/test_event_monitor.jl")# C4: event monitor (IF/THEN/ENDIF) evaluator + firing
    # include("integration/test_sndb.jl")   # enabled at C6
    # include("integration/test_snt02.jl")  # enabled at C8
end
