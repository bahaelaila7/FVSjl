# test_nocalib.jl — NOCALIB keyword (option 56) → disable DG self-calibration vs Fortran.
#
# NOCALIB clears LDGCAL for a species (0/all, −N group, code) so the DG self-calibration
# COR fit is skipped — the species uses its uncalibrated diameter growth (dgdriv.f:567).
# Before the fix the keyword was unrecognized AND control.dg_calib_sp (LDGCAL) was declared
# but dead (defaulted all-false, never read). Now it defaults all-TRUE and gates the COR fit
# in calibrate_diameter_growth!. (FVSjl does no large-tree HT self-calibration, so the SN
# LHTCAL side is naturally inert.) Checks: SETS the flags; FIRES (≠ calibrated); FORTRAN exact.

using Test, FVSjl

const _NC_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_nc_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]
_nc_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]

@testset "NOCALIB → disable DG self-calibration" begin
    # 1. SETS — default is all-on (LDGCAL); NOCALIB 0 clears every species
    s = FVSjl.StandState(FVSjl.Southern())
    @test all(s.control.dg_calib_sp)
    FVSjl.kw_nocalib!(s, FVSjl.KeywordRecord("NOCALIB ", "", ["0", fill("", 11)...],
                      zeros(Float32, 12), [true, falses(11)...], 12, FVSjl.KW_OK, 0))
    @test !any(s.control.dg_calib_sp)
    # a single species code clears only that species
    s2 = FVSjl.StandState(FVSjl.Southern())
    FVSjl.kw_nocalib!(s2, FVSjl.KeywordRecord("NOCALIB ", "", ["LP", fill("", 11)...],
                      zeros(Float32, 12), [true, falses(11)...], 12, FVSjl.KW_OK, 0))
    @test count(!, s2.control.dg_calib_sp) == 1

    key = joinpath(_NC_DIR, "nocalib.key"); sav = joinpath(_NC_DIR, "nocalib.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "nocalib scenario not available"
    else
        jl = _nc_rows(FVSjl.run_keyfile(key; faithful = true))
        # 2. FIRES — NOCALIB changes the projection vs the calibrated stand
        offkey = joinpath(_NC_DIR, "_nocalib_off.key")
        cp(joinpath(_NC_DIR, "nocalib.tre"), joinpath(_NC_DIR, "_nocalib_off.tre"); force = true)
        open(offkey, "w") do io
            for l in eachline(key); startswith(strip(l), "NOCALIB") || println(io, l); end
        end
        try
            off = _nc_rows(FVSjl.run_keyfile(offkey; faithful = true))
            @test any(jl[i][10] != off[i][10] for i in eachindex(jl))
        finally
            rm(offkey; force = true); rm(joinpath(_NC_DIR, "_nocalib_off.tre"); force = true)
        end
        # 3. FORTRAN — uncalibrated DG matches live Fortran
        ft = _nc_base(sav)
        @test length(jl) == length(ft)
        for (j, f) in zip(jl, ft)
            for c in (3, 4, 5, 6, 7); @test parse(Int, j[c]) == parse(Int, f[c]); end
            @test abs(parse(Float64, j[8]) - parse(Float64, f[8])) <= 0.05
        end
    end
end
