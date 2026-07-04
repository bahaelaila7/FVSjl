# =============================================================================
# test_ls_sitesweep.jl — LS site-productivity sweep (deterministic growth), vs live FVSls
#
# The canonical lst01 validates ONE site index (SITECODE 2 / 60). Real analysts run a
# stand across a productivity range, which exercises the LS site-index path — ls_site_index_
# setup!'s 68×68 SICOEF fan-out → the DGF SITEC term → the whole growth trajectory — at
# indices the single-site test never touches. This sweep runs the lst01 inventory at site
# indices 40/50/70/80 under NOTRIPLE (deterministic: no stochastic tripling spread) and
# asserts jl == live FVSls BIT-EXACT across the ENTIRE projection (all cycles, all columns).
#
# This test originally exposed a Δ2-5 TPA terminal residual at non-canonical site indices,
# which was cornered end-to-end (terminal mortality residual → BAMAX BA-cap amplification →
# per-record diameter growth → sp5/sp41 DGF) to a REAL BUG in the LS QMDGE5 cap, now FIXED:
#   FVS caps the STAND-WIDE QMDGE5 in place as it walks species in INDEX order (ls/dgf.f:362-390),
#   so a species sees QMDGE5 capped by ALL LOWER-INDEXED present cap-species — e.g. white pine
#   (sp5, uncapped) sees the 13" cap that jack pine (sp1, cap-13) applied upstream. jl had applied
#   a LOCAL per-tree cap, so sp5 kept its uncapped QMD (14.5 vs live 13.0), biasing its RDBH/RDBHSQ
#   growth terms → a tiny DG error the terminal BA-cap mortality amplified to Δ2-5.
# dgf! now replicates the cumulative species-order cap mutation. The LS mortality (background /
# VARMRT / BAMAX) was proven bit-faithful throughout — the residual was upstream, in growth.
# =============================================================================

using Test, FVSjl

const _LSSITE_DIR = joinpath(@__DIR__, "..", "fixtures", "sitesweep")

function _lssite_rows(txt::AbstractString)
    d = Dict{Int,Vector{Int}}()
    for l in split(txt, "\n")
        f = split(strip(l))
        (length(f) >= 7 && tryparse(Int, f[1]) !== nothing && startswith(f[1], "20")) || continue
        d[parse(Int, f[1])] = [parse(Int, f[i]) for i in 3:7]   # TPA,BA,SDI,CCF,TopHt
    end
    return d
end

@testset "LS site-productivity sweep — deterministic growth BIT-EXACT (vs live FVSls)" begin
    if !isfile(joinpath(_LSSITE_DIR, "ls_nt50.key"))
        @info "sitesweep fixtures absent; skipping"
    else
        for si in (40, 50, 70, 80)
            jl = cd(_LSSITE_DIR) do
                FVSjl.run_keyfile("ls_nt$(si).key"; variant = LakeStates(), output = :sum)
            end
            J = _lssite_rows(jl)
            L = _lssite_rows(read(joinpath(_LSSITE_DIR, "ls_nt$(si).live.sum"), String))
            for y in sort(collect(keys(L)))
                haskey(J, y) || (@test haskey(J, y); continue)
                # BIT-EXACT across all 5 stand columns, every cycle — the QMDGE5 cap fix.
                @test ("SI$si", y, J[y]) == ("SI$si", y, L[y])
            end
        end
    end
end
