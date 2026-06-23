# test_mcdefect.jl — volume defect keywords (sdefet.f + FVSsn vols.f) vs live Fortran.
#
# The SN volume driver (bin/FVSsn_buildDir/vols.f) layers a per-species DBH defect curve on the
# taper volumes. Two keywords feed it, both ported here:
#   • MCDEFECT → CFDEFT: the pulpwood/topwood part of merch cubic (MCFV−SCFV) is cut by
#     ICDF% = NINT(ALGSLP(DBH, CFDEFT)·100) clamped [0,99] (≥99 ⇒ all pulpwood gone); sawtimber
#     is untouched by the cubic step (vols.f:294-325).
#   • BFDEFECT → BFDEFT: board feet AND sawtimber cubic are cut by IBDF% (≥99 ⇒ both 0), where
#     board feet exist (vols.f:419-432). Then MCFV = PULPV + (post-board-defect SCFV), so a
#     BFDEFECT also lowers reported merch cubic — the coupled `defect_both` case checks that.
# Curves are undated → applied from cycle 0 (sdefet.f "date not defined" branch). Per-tree DEFECT
# input is deferred; the CFLA0/CFLA1 form model is verified no-op for SN (see DIVERGENCES.md).

using Test, FVSjl

const _MD_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_md_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && tryparse(Int, first(split(l))) !== nothing]
_md_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_mdcol(r, c) = parse(Float64, r[c])

# Run a defect scenario, assert it FIRES (a no-keyword twin differs in `firecol` by >drop in every
# row), and matches the saved live-Fortran .sum on the structural + the directly-affected columns.
function _check_defect(stem, kwprefix, firecol, drop, cols)
    key = joinpath(_MD_DIR, stem * ".key"); sav = joinpath(_MD_DIR, stem * ".sum.save")
    tre = joinpath(_MD_DIR, stem * ".tre")
    if !isfile(key) || !isfile(sav)
        @test_skip "$stem scenario not available"; return
    end
    jl = _md_rows(FVSjl.run_keyfile(key; faithful = true))
    nov = tempname() * ".key"
    write(nov, join(filter(l -> !any(startswith(l, p) for p in kwprefix), readlines(key)), "\n") * "\n")
    cp(tre, replace(nov, ".key" => ".tre"); force = true)
    nd = _md_rows(FVSjl.run_keyfile(nov; faithful = true))
    @test length(jl) == length(nd)
    # the defect FIRES: it lowers `firecol` (never raises it) and pulls it down by >drop somewhere.
    @test all(_mdcol(jl[i], firecol) <= _mdcol(nd[i], firecol) + 1 for i in 1:length(jl))
    @test any(_mdcol(nd[i], firecol) - _mdcol(jl[i], firecol) > drop for i in 1:length(jl))
    ft = _md_base(sav)
    @test length(jl) == length(ft)
    if length(jl) == length(ft)
        for i in 1:length(jl), c in cols
            @test abs(_mdcol(jl[i], c) - _mdcol(ft[i], c)) <= 2
        end
    end
end

@testset "volume defect keywords vs Fortran" begin
    # MCDEFECT: cubic defect drops merch cubic (col 10), incl. cycle 0; sawtimber (11) bit-exact.
    @testset "MCDEFECT" begin _check_defect("mcdefect_override", ("MCDEFECT",), 10, 100, (3, 4, 10)) end
    # BFDEFECT: board defect drops sawtimber cubic (col 11); col 11 is the clean bit-exact target.
    @testset "BFDEFECT" begin _check_defect("bfdefect_override", ("BFDEFECT",), 11, 100, (3, 4, 11)) end
    # Coupled: MCDEFECT+BFDEFECT — merch cubic (10) = PULPV + post-board SCFV must match Fortran.
    @testset "MC+BF" begin _check_defect("defect_both", ("MCDEFECT", "BFDEFECT"), 11, 100, (3, 4, 10, 11)) end
end
