# test_sdicalc.jl — SDICALC keyword (option 1400) → SDI method (Zeide/Reineke) + thresholds vs Fortran.
#
# SDICALC sets the SHARED LZEIDE flag (control.zeide_sdi) + the DBHZEIDE/DBHSTAGE min-DBH
# thresholds. LZEIDE drives BOTH the reported .sum SDI column (stand_sdi) AND the SDImax
# self-thinning mortality (mortality.jl) — so a method change moves the multi-cycle TPA, not
# just cycle-0's SDI. Before the fix the keyword was unrecognized and stand_sdi was Zeide-only.
# Checks: SETS the flag/thresholds; FIRES (≠ default); FORTRAN bit-exact MULTI-CYCLE.

using Test, FVSjl

const _SC_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_sc_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]
_sc_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]

@testset "SDICALC → SDI method + thresholds" begin
    # 1. SETS — field 3 < 1 (blank) ⇒ Reineke; ≥ 1 ⇒ Zeide; thresholds from fields 1/2
    s = FVSjl.StandState(FVSjl.Southern())
    FVSjl.kw_sdicalc!(s, FVSjl.KeywordRecord("SDICALC ", "", ["2.0", "3.0", "", fill("", 9)...],
                      Float32[2, 3, 0, zeros(Float32, 9)...], [true, true, false, falses(9)...], 12, FVSjl.KW_OK, 0))
    @test s.control.zeide_sdi == false && s.control.dbh_stage == 2f0 && s.control.dbh_zeide == 3f0
    FVSjl.kw_sdicalc!(s, FVSjl.KeywordRecord("SDICALC ", "", ["", "", "1", fill("", 9)...],
                      Float32[0, 0, 1, zeros(Float32, 9)...], [false, false, true, falses(9)...], 12, FVSjl.KW_OK, 0))
    @test s.control.zeide_sdi == true

    key = joinpath(_SC_DIR, "sdicalc.key"); sav = joinpath(_SC_DIR, "sdicalc.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "sdicalc scenario not available"
    else
        jl = _sc_rows(FVSjl.run_keyfile(key; faithful = true))
        # 2. FIRES — SDICALC→Reineke changes the SDI column AND (via mortality) the TPA
        offkey = joinpath(_SC_DIR, "_sdicalc_off.key")
        cp(joinpath(_SC_DIR, "sdicalc.tre"), joinpath(_SC_DIR, "_sdicalc_off.tre"); force = true)
        open(offkey, "w") do io
            for l in eachline(key); startswith(strip(l), "SDICALC") || println(io, l); end
        end
        try
            off = _sc_rows(FVSjl.run_keyfile(offkey; faithful = true))
            @test any(jl[i][5] != off[i][5] for i in eachindex(jl))   # SDI column differs
            @test any(jl[i][3] != off[i][3] for i in eachindex(jl))   # TPA differs (mortality follows)
        finally
            rm(offkey; force = true); rm(joinpath(_SC_DIR, "_sdicalc_off.tre"); force = true)
        end
        # 3. FORTRAN — multi-cycle TPA/BA/SDI/CCF/TopHt/QMD match live Fortran
        ft = _sc_base(sav)
        @test length(jl) == length(ft)
        for (j, f) in zip(jl, ft)
            for c in (3, 4, 5, 6, 7); @test parse(Int, j[c]) == parse(Int, f[c]); end
            @test parse(Float64, j[8]) == parse(Float64, f[8])
        end
    end
end
