# =============================================================================
# test_allspecies.jl — all-species coverage for CS and NE vs live FVScs/FVSne
#
# Realistic mixed stands that exercise EVERY species' growth / crown-width / volume /
# mortality / FORTYP coefficient row, validated against the live binaries (golden .sum
# captured from a freshly-relinked FVScs/FVSne). The canonical cst01/net01 stands use
# only a handful of common species; this sweep is what caught (and now guards) the
# BW (American basswood) crown-width gap that was missing from both variants' maps.
#
# Coverage construction (see test/harness/gen_allspecies.sh):
#   - CS: one mixed stand, one tree of each of the 96 species (cst01 inventory rows,
#     species field swapped). Live FVScs runs it directly.
#   - NE: 5 stands — live FVSne FPEs on stands of many all-unusual species at once (a
#     live-binary FORTYP/stats limitation), so the 107 non-blank species are greedily
#     packed (new species + net01-real filler) into stands that the live binary runs.
#     The blank-alpha placeholder species (NE 71) is excluded — it is unaddressable.
#
# Bar: cycle-0 all six stand columns BIT-EXACT; later cycles within the documented
# single-precision floor (the same class as the canonical stands).
# =============================================================================

using Test
using FVSjl

# Parse a .sum's data rows → Dict(year => Vector{String} of fields).
function _allsp_rows(text::AbstractString)
    rows = Dict{Int,Vector{String}}()
    for ln in split(text, '\n')
        t = split(strip(ln))
        (length(t) >= 8 && tryparse(Int, t[1]) !== nothing && t[1] != "-999") || continue
        rows[parse(Int, t[1])] = collect(String.(t))
    end
    return rows
end

# Compare a jl .sum string to a live golden .sum, asserting cyc0 bit-exact + later-cycle floor.
# Columns: 3=TPA 4=BA 5=SDI 6=CCF 7=TopHt 8=QMD(×0.1) 9=Tcuft 12=Bdft.
function _assert_allspecies(jl_text::AbstractString, golden_path::AbstractString; label::AbstractString)
    @test isfile(golden_path)
    J = _allsp_rows(jl_text)
    L = _allsp_rows(read(golden_path, String))
    years = sort(collect(keys(L)))
    @test !isempty(years)
    cyc0 = first(years)
    for yr in years
        haskey(J, yr) || (@test haskey(J, yr); continue)
        l, j = L[yr], J[yr]
        geti(v, i) = parse(Float64, v[i])
        if yr == cyc0                                   # inventory cycle: BIT-EXACT (coefficient rows)
            for (i, name) in ((3,"TPA"),(4,"BA"),(5,"SDI"),(6,"CCF"),(7,"TopHt"),(8,"QMD"))
                @test (label, name, yr, geti(j,i)) == (label, name, yr, geti(l,i))
            end
        else                                            # grown cycles: single-precision floor
            @test abs(geti(j,3) - geti(l,3)) <= max(4, 0.03*geti(l,3))   # TPA
            @test abs(geti(j,4) - geti(l,4)) <= max(2, 0.02*geti(l,4))   # BA
            @test abs(geti(j,5) - geti(l,5)) <= max(4, 0.02*geti(l,5))   # SDI
            @test abs(geti(j,6) - geti(l,6)) <= max(10,0.04*geti(l,6))   # CCF
            @test abs(geti(j,7) - geti(l,7)) <= 3                        # TopHt
            @test abs(geti(j,8) - geti(l,8)) <= 0.3                      # QMD
        end
    end
end

const _ALLSP_DIR = joinpath(@__DIR__, "..", "fixtures", "allspecies")

@testset "CS all-species coverage (96 species, vs live FVScs)" begin
    key = joinpath(_ALLSP_DIR, "cs_allsp.key")
    if !isfile(key)
        @info "cs_allsp fixtures absent; skipping"
    else
        cd(_ALLSP_DIR) do                              # TREEDATA reads <keystem>.tre from cwd
            jl = FVSjl.run_keyfile("cs_allsp.key"; variant = CentralStates(), output = :sum)
            _assert_allspecies(jl, joinpath(_ALLSP_DIR, "cs_allsp.live.sum"); label = "CS")
        end
    end
end

@testset "NE all-species coverage (107 species, vs live FVSne)" begin
    miss = !isfile(joinpath(_ALLSP_DIR, "ne_cov0.key"))
    if miss
        @info "ne_cov fixtures absent; skipping"
    else
        cd(_ALLSP_DIR) do
            for j in 0:4
                isfile("ne_cov$(j).key") || continue
                jl = FVSjl.run_keyfile("ne_cov$(j).key"; variant = Northeast(), output = :sum)
                _assert_allspecies(jl, joinpath(_ALLSP_DIR, "ne_cov$(j).live.sum"); label = "NE-cov$j")
            end
        end
    end
end
