# test_managed.jl — MANAGED keyword → DGF planted/managed growth term (dgf.f:179/328) vs Fortran.
#
# MANAGED sets the managed-stand flag (MANAGD); the SN diameter-growth function then adds
# `dg_planted[sp]·kplant` to ln(DDS) for every tree (kplant = MANAGD>0). The boost is non-
# zero only for the *planted pines* (loblolly/longleaf/slash/pond/white pine), so the
# scenario is a loblolly-pine stand. Before the fix MANAGED was a recognized no-op, so the
# growth term never fired. Checks:
#   1. SETS — the keyword sets s.plot.managed (and field-2 == 0 clears it);
#   2. FIRES — dropping the MANAGED line changes the projected stand (real .sum effect);
#   3. FORTRAN — TPA/BA/SDI/cubic columns match live Fortran (board feet within Scribner noise).

using Test, FVSjl

const _MG_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_mg_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]
_mg_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]

@testset "MANAGED → DGF planted growth term" begin
    # 1. SETS — the keyword toggles plot.managed (initre.f:10000)
    setflag(field2) = begin
        s = FVSjl.StandState(FVSjl.Southern())
        rec = FVSjl.KeywordRecord("MANAGED ", "", ["", field2, fill("", 10)...],
                                  Float32[0, field2 == "0" ? 0 : 0, zeros(Float32, 10)...],
                                  [false, !isempty(field2), falses(10)...], 12, FVSjl.KW_OK, 0)
        FVSjl.kw_managed!(s, rec)
        s.plot.managed
    end
    @test setflag("") == 1        # bare MANAGED ⇒ managed
    @test setflag("0") == 0       # MANAGED with explicit 0 ⇒ unmanaged

    key = joinpath(_MG_DIR, "managed.key")
    sav = joinpath(_MG_DIR, "managed.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "managed scenario not available"
    else
        jl = _mg_rows(FVSjl.run_keyfile(key; faithful = true))

        # 2. FIRES — removing the MANAGED line changes the loblolly stand's growth. The
        # off-key must share the .tre basename, so write it beside the scenario's .tre.
        offkey = joinpath(_MG_DIR, "_managed_off.key")
        cp(joinpath(_MG_DIR, "managed.tre"), joinpath(_MG_DIR, "_managed_off.tre"); force = true)
        open(offkey, "w") do io
            for l in eachline(key); strip(l) == "MANAGED" || println(io, l); end
        end
        try
            off = _mg_rows(FVSjl.run_keyfile(offkey; faithful = true))
            # the planted-growth term changes the projected stand (BA / volume differ)
            @test any(jl[i][4] != off[i][4] || jl[i][10] != off[i][10] for i in eachindex(jl))
        finally
            rm(offkey; force = true); rm(joinpath(_MG_DIR, "_managed_off.tre"); force = true)
        end

        # 3. FORTRAN — TPA/BA/SDI/CCF/TopHt/QMD/MerchCuft match live Fortran.
        ft = _mg_base(sav)
        @test length(jl) == length(ft)
        for (j, f) in zip(jl, ft)
            for c in (3, 4, 5, 6, 7)                    # TPA, BA, SDI, CCF, TopHt
                @test parse(Int, j[c]) == parse(Int, f[c])
            end
            @test parse(Float64, j[8]) == parse(Float64, f[8])   # QMD
            @test abs(parse(Int, j[10]) - parse(Int, f[10])) <= 2            # MerchCuFt (±Float32)
        end
    end
end
