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
# Later-cycle tolerance set — set just above the MEASURED live-vs-jl floor. cyc0 is ALWAYS bit-exact
# (the coefficient rows); the grown-cycle drift is the accepted tripling-spread single-precision floor,
# amplified by dense near-SDImax mortality. Defaults fit CS/NE (dense CS mix); LS passes a slightly wider
# set (its 10-yr cycle + denser mix drift the volume/CCF/QMD columns ~1-1.5%). Each `(abs, pct)` bound is
# ~1.3-1.7× the measured floor: tight enough to catch a ≥3% coefficient regression, loose enough not to
# flake on ULP. Columns: 3 TPA · 4 BA · 5 SDI · 6 CCF · 7 TopHt · 8 QMD · 9 Tcuft · 10 Mcuft · 11 Scuft · 12 Bdft.
const _ALLSP_TOL_DEFAULT = (tpa=(3,0.025), ba=(2,0.015), sdi=(3,0.015), ccf=(4,0.020),
                            topht=(2,0.0), qmd=(0.25,0.0), tcuft=(2,0.012), mcuft=(2,0.012),
                            scuft=(2,0.018), bdft=(5,0.018))
const _ALLSP_TOL_LS = (tpa=(3,0.030), ba=(2,0.015), sdi=(3,0.015), ccf=(5,0.025),
                       topht=(2,0.0), qmd=(0.4,0.0), tcuft=(2,0.015), mcuft=(2,0.015),
                       scuft=(3,0.020), bdft=(6,0.020))

function _assert_allspecies(jl_text::AbstractString, golden_path::AbstractString;
                            label::AbstractString, tol = _ALLSP_TOL_DEFAULT)
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
            # stand columns + ALL FOUR VOLUME columns (Tcuft/Mcuft/Scuft/Bdft) — the cyc0 volume is
            # bit-exact for every species here, so this guards the per-species merch/saw/board standards
            # (the SN sweep found green-ash/cypress merch divergences; all four variants are clean at cyc0).
            for (i, name) in ((3,"TPA"),(4,"BA"),(5,"SDI"),(6,"CCF"),(7,"TopHt"),(8,"QMD"),
                              (9,"Tcuft"),(10,"Mcuft"),(11,"Scuft"),(12,"Bdft"))
                @test (label, name, yr, geti(j,i)) == (label, name, yr, geti(l,i))
            end
        else
            chk(i, t) = @test abs(geti(j,i) - geti(l,i)) <= max(t[1], t[2]*geti(l,i))
            chk(3, tol.tpa); chk(4, tol.ba); chk(5, tol.sdi); chk(6, tol.ccf)
            chk(7, tol.topht); chk(8, tol.qmd)
            chk(9, tol.tcuft); chk(10, tol.mcuft); chk(11, tol.scuft); chk(12, tol.bdft)
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

@testset "SN all-species coverage (90 species, vs live FVSsn)" begin
    if !isfile(joinpath(_ALLSP_DIR, "sn_cov0.key"))
        @info "sn_cov fixtures absent; skipping (run test/harness/gen_allspecies_snls.sh)"
    else
        cd(_ALLSP_DIR) do
            for j in 0:9
                isfile("sn_cov$(j).key") || continue
                jl = FVSjl.run_keyfile("sn_cov$(j).key"; variant = Southern(), output = :sum)
                _assert_allspecies(jl, joinpath(_ALLSP_DIR, "sn_cov$(j).live.sum"); label = "SN-cov$j")
            end
        end
    end
end

@testset "LS all-species coverage (67 species, vs live FVSls)" begin
    if !isfile(joinpath(_ALLSP_DIR, "ls_cov0.key"))
        @info "ls_cov fixtures absent; skipping (run test/harness/gen_allspecies_snls.sh)"
    else
        cd(_ALLSP_DIR) do
            for j in 0:9
                isfile("ls_cov$(j).key") || continue
                jl = FVSjl.run_keyfile("ls_cov$(j).key"; variant = LakeStates(), output = :sum)
                _assert_allspecies(jl, joinpath(_ALLSP_DIR, "ls_cov$(j).live.sum");
                                   label = "LS-cov$j", tol = _ALLSP_TOL_LS)
            end
        end
    end
end
