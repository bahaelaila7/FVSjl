# test_resetage.jl — RESETAGE keyword (option 92, resage.f act 443) → rebase stand age vs Fortran.
#
# RESETAGE makes the stand age equal field 2 at the activity's date (IAGE = age−IDT+IY(1)).
# RESAGE runs after DISPLY, so the reset year's own .sum row keeps the old age and the rebase
# shows the FOLLOWING row. It affects only the AGE + MAI columns. No scheduler is needed — it
# is a pure function of the row year, applied in summary_row. Checks: SETS the fields; the AGE
# column rebases (≠ no-RESETAGE); FORTRAN bit-exact (AGE/MAI match).

using Test, FVSjl

const _RA_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_ra_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]
_ra_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]

@testset "RESETAGE → rebase stand age" begin
    # 1. SETS — resolves the date to a year + stores the age (default none = -1)
    s = FVSjl.StandState(FVSjl.Southern())
    @test s.control.age_reset_year == Int32(-1)
    s.control.cycle_year[1] = Int32(1990); s.control.year = 5f0
    FVSjl.kw_resetage!(s, FVSjl.KeywordRecord("RESETAGE", "", ["2000", "30", fill("", 10)...],
                       Float32[2000, 30, zeros(Float32, 10)...], [true, true, falses(10)...], 12, FVSjl.KW_OK, 0))
    @test s.control.age_reset_year == Int32(2000) && s.control.age_reset_age == Int32(30)

    key = joinpath(_RA_DIR, "resetage.key"); sav = joinpath(_RA_DIR, "resetage.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "resetage scenario not available"
    else
        jl = _ra_rows(FVSjl.run_keyfile(key; faithful = true))
        # 2. FIRES — the AGE column (col 2) rebases after the reset year
        offkey = joinpath(_RA_DIR, "_resetage_off.key")
        cp(joinpath(_RA_DIR, "resetage.tre"), joinpath(_RA_DIR, "_resetage_off.tre"); force = true)
        open(offkey, "w") do io
            for l in eachline(key); startswith(strip(l), "RESETAGE") || println(io, l); end
        end
        try
            off = _ra_rows(FVSjl.run_keyfile(offkey; faithful = true))
            @test any(jl[i][2] != off[i][2] for i in eachindex(jl))   # AGE differs
            # the reset year's row keeps the old age; the next row is rebased
            i2000 = findfirst(r -> r[1] == "2000", jl)
            @test i2000 !== nothing && jl[i2000][2] == off[i2000][2]   # reset-row age unchanged
            @test jl[i2000 + 1][2] == "35"                            # next row = 30 + 5
        finally
            rm(offkey; force = true); rm(joinpath(_RA_DIR, "_resetage_off.tre"); force = true)
        end
        # 3. FORTRAN — AGE (col 2) + MAI (col 25) match live Fortran
        ft = _ra_base(sav)
        @test length(jl) == length(ft)
        for (j, f) in zip(jl, ft)
            @test parse(Int, j[2]) == parse(Int, f[2])              # AGE
            @test parse(Float64, j[25]) == parse(Float64, f[25])    # MAI — BIT-EXACT (measured Δ=0 all rows; the
                                                                    # rendered MAI = total cuft/age renders identically; was ≤0.2 padding)
        end
    end
end
