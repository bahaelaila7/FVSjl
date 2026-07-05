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
# Per-variant grown-cycle tolerances, re-measured to the TRUE observed floor across every coverage
# stand / cycle / column (test/harness measurement, not guessed). The four variants split sharply:
#
#   • NE and LS coverage are 100% BIT-EXACT vs the live binary on EVERY column, cycle and species
#     (max|Δ| = 0 everywhere) — so their grown-cycle bound is 0 = `==`, same as cyc0.
#   • SN is bit-exact on BA/TopHt, and — measured per coverage file (2026-07-05) — NEAR-BIT-EXACT on
#     ALL columns for sn_cov0..3 (max|Δ|=0 bar a ≤1 tcuft/mcuft print-knife-edge). ONLY sn_cov4 carries
#     the accumulated DGSCOR + tripling tail (its coverage group holds the WK3-calibrated sp33/65-family
#     species — the same class as the @test_broken). So the wide envelope applies ONLY to cov4; cov0..3
#     use a near-exact bound. A uniform loose SN bound would MASK the bit-exactness of four files (rule #4).
#   • CS is bit-exact on BA/SDI but is a SINGLE unpartitioned 96-species stand, so it accumulates every
#     species' sub-ULP per-cycle DBH-growth + tripling-spread residual into the nonlinear density/volume
#     sums — the ACCEPTED, DOCUMENTED aggregate DGSCOR + tripling class (same family as test_timeint's
#     cuft tail and the sp33/65 WK3 @test_broken). Bounds = the observed DETERMINISTIC envelope (absolute,
#     not a loosened percentage). CS board feet carry the largest (0.95% of a big per-acre Scribner sum)
#     because board is the most nonlinear column; it cannot be isolated further (one stand, no sub-covers).
const _ALLSP_TOL_BITEXACT = (tpa=(0,0.0), ba=(0,0.0), sdi=(0,0.0), ccf=(0,0.0), topht=(0,0.0),
                             qmd=(0.0,0.0), tcuft=(0,0.0), mcuft=(0,0.0), scuft=(0,0.0), bdft=(0,0.0))
const _ALLSP_TOL_NE = _ALLSP_TOL_BITEXACT
const _ALLSP_TOL_LS = _ALLSP_TOL_BITEXACT
# SN is NOT a single envelope: measured per-coverage-file (2026-07-05), sn_cov0..3 are NEAR-BIT-EXACT
# (max|Δ| = 0 on every stand+volume column bar a single ≤1 tcuft/mcuft print-knife-edge), and ONLY
# sn_cov4 carries the accumulated DGSCOR + tripling tail (that coverage group holds the WK3-calibrated
# sp33/65-family species — the same class as the @test_broken). So cov0..3 get the near-exact bound
# (== + print-knife-edge) and the wide envelope applies ONLY to cov4 — a uniform loose bound across all
# five would MASK the bit-exactness of four of them (doctrine #4).
const _ALLSP_TOL_SN_NEAREXACT = (tpa=(0,0.0), ba=(0,0.0), sdi=(0,0.0), ccf=(0,0.0), topht=(0,0.0),
                       qmd=(0.0,0.0), tcuft=(1,0.0), mcuft=(1,0.0), scuft=(0,0.0), bdft=(0,0.0))
const _ALLSP_TOL_SN_DGSCOR = (tpa=(2,0.0), ba=(0,0.0), sdi=(1,0.0), ccf=(1,0.0), topht=(0,0.0),
                       qmd=(0.1,0.0), tcuft=(3,0.0), mcuft=(3,0.0), scuft=(4,0.0), bdft=(54,0.0))
const _ALLSP_TOL_SN = _ALLSP_TOL_SN_DGSCOR   # back-compat alias (fallback)
const _ALLSP_TOL_CS = (tpa=(1,0.0), ba=(0,0.0), sdi=(0,0.0), ccf=(4,0.0), topht=(1,0.0),
                       qmd=(0.1,0.0), tcuft=(21,0.0), mcuft=(20,0.0), scuft=(20,0.0), bdft=(464,0.0))
const _ALLSP_TOL_DEFAULT = _ALLSP_TOL_CS   # CS is the widest-envelope variant (kept as the fallback)

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
            _assert_allspecies(jl, joinpath(_ALLSP_DIR, "cs_allsp.live.sum"); label = "CS", tol = _ALLSP_TOL_CS)
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
                _assert_allspecies(jl, joinpath(_ALLSP_DIR, "ne_cov$(j).live.sum"); label = "NE-cov$j", tol = _ALLSP_TOL_NE)
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
                # cov0..3 are near-bit-exact; only cov4 (the WK3/DGSCOR-tail species group) needs the envelope.
                tol = j >= 4 ? _ALLSP_TOL_SN_DGSCOR : _ALLSP_TOL_SN_NEAREXACT
                _assert_allspecies(jl, joinpath(_ALLSP_DIR, "sn_cov$(j).live.sum"); label = "SN-cov$j", tol = tol)
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
