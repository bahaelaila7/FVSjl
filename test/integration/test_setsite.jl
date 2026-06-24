# test_setsite.jl — SETSITE keyword (opt 138, act 120, initre.f:13800 + grincr.f:89-200).
#
# SETSITE schedules a MID-RUN site change: at the matching cycle it resets the per-species site
# index (SITEAR, directly or as a % change, clamped ≥ 1) and optionally BAMAX / SDImax, then
# recomputes the site-dependent DG constants (RCON = `dgcons!`). The new site index feeds both the
# diameter- and height-growth models from that cycle on. Checks:
#   1. UNIT — species_selector decodes the field (0/all, −group, index); the keyword schedules an
#      act-120 ScheduledActivity with the right params;
#   2. FORTRAN — a SETSITE that raises all-species site index 60→80 at 2000 matches live Fortran on
#      every structural column (volume within ±Scribner Float32 noise) and changes the projection.

using Test, FVSjl

const _SS_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_ss_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]
_ss_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                                y !== nothing && 1900 < y < 2100)]

@testset "SETSITE → scheduled mid-run site change" begin
    # 1. UNIT — selector decode + scheduling.
    s = FVSjl.StandState(FVSjl.Southern())
    @test FVSjl.species_selector(s, "")  == 0      # blank ⇒ all
    @test FVSjl.species_selector(s, "0") == 0      # 0 ⇒ all
    @test FVSjl.species_selector(s, "-3") == -3    # group −3
    # schedule a SETSITE at 2000, all species (field 4 blank), site index 80 (field 5), flag 0
    fields = fill("", 12); fields[1] = "2000"; fields[5] = "80"; fields[6] = "0"
    values = zeros(Float32, 12); values[1] = 2000f0; values[5] = 80f0; values[6] = 0f0
    present = falses(12); present[1] = true; present[5] = true; present[6] = true
    rec = FVSjl.KeywordRecord("SETSITE", "", fields, values, present, 7, FVSjl.KW_OK, 0)
    FVSjl.kw_setsite!(s, rec)
    @test length(s.control.schedule) == 1
    a = s.control.schedule[1]
    @test a.icflag == 120 && a.year == 2000
    @test a.params[4] == 80f0 && a.params[5] == 0f0    # site index = 80, flag = direct

    # 2. FORTRAN — the setsite scenario (SI 60→80 @ 2000) vs live Fortran.
    key = joinpath(_SS_DIR, "setsite.key")
    sav = joinpath(_SS_DIR, "setsite.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "setsite scenario not available"
    else
        jl = _ss_rows(FVSjl.run_keyfile(key; faithful = true))
        ft = _ss_base(sav)
        @test length(jl) == length(ft)
        for (j, f) in zip(jl, ft)
            @test j[1] == f[1]                                   # YEAR
            for col in (3, 4, 5, 6, 7, 8)                        # TPA/BA/SDI/CCF/TopHt/QMD
                @test j[col] == f[col]
            end
            for col in (9, 10, 11, 12)                           # TCuFt/MCuFt/SCuFt/BdFt
                @test abs(parse(Int, j[col]) - parse(Int, f[col])) <= 2
            end
        end
        # the site boost actually changed the projection (vs the same key with SETSITE removed)
        offkey = joinpath(_SS_DIR, "_setsite_off.key")
        cp(joinpath(_SS_DIR, "setsite.tre"), joinpath(_SS_DIR, "_setsite_off.tre"); force = true)
        open(offkey, "w") do io
            for l in eachline(key); startswith(strip(l), "SETSITE") || println(io, l); end
        end
        try
            off = _ss_rows(FVSjl.run_keyfile(offkey; faithful = true))
            # post-2000 stand volume is higher than the unchanged-site baseline
            i2005 = findfirst(r -> r[1] == "2005", jl)
            @test i2005 !== nothing && parse(Int, jl[i2005][9]) > parse(Int, off[i2005][9])
        finally
            rm(offkey; force = true); rm(joinpath(_SS_DIR, "_setsite_off.tre"); force = true)
        end
    end
end
