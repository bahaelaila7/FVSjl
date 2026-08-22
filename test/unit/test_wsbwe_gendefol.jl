# Western Spruce Budworm (WSBWE) GENDEFOL / BUDLITE — bit-exact regression for the
# ported bwelit.f BUDLITE population/defoliation core (`wsbwe_bwelit!`).
#
# GENDEFOL runs a SELF-GENERATING budworm outbreak (vs the manual DEFOL path): the
# WEATHER keyword supplies a weather-stats file, bwewea.f derives per-year params,
# and bwelit.f (+ bwedie.f survival + the bwebkem.f block data) generates the
# defoliation into the FNEW/FOLD foliage cells, which the SHARED BWEAGE/BWEDAM/
# BWEPDM chain (already validated on the manual-DEFOL path) turns into per-tree
# growth-loss. See src/engine/wsbwe.jl and the port audit in
# scratchpad/wsbwe/gendefol/staged/VALIDATION.md.
#
# SYNTHETIC vs FAITHFUL: the weather INPUT (synwx.txt) is synthetic (disclaims
# ecological realism); the MODEL is a faithful port validated bit-exact against the
# live FVSem_wsbwe oracle fed the IDENTICAL synthetic file. This test replays the
# ported `wsbwe_bwelit!` against oracle golden dumps captured on the s248 stand for
# the three signal-bearing outbreak years — 1990 (the dominant pulse), 1991 (the
# tail), 1992 (population crash → the bwelit.f:647 IF(BW<1) dead-larvae early
# return). Every foliage cell (FNEW/FOLD1/FOLD2/FREM) and the carried EGGS must be
# Float32-hex byte-identical to the oracle.
#
# GOLDEN PROVENANCE: test/fixtures/wsbwe/gendefol/bwelit_io_golden.txt (per-cell
# IN/OUT dump from the instrumented FVSem_wsbwe, gdins/bwelit.f) +
# bwe_wea_golden.txt (bwewea.f per-year params). The IDENTICAL synthetic input
# (synwx.txt) is fed to both engines, so the bit-exactness is genuine.

using Test
using FVSjl
const _GDW = FVSjl

_h2f(x::AbstractString) = reinterpret(Float32, parse(UInt32, x; base = 16))
_u32(x::Float32) = reinterpret(UInt32, x)

const _GD_DIR = joinpath(@__DIR__, "..", "fixtures", "wsbwe", "gendefol")

# --- parse the per-year weather golden: DAYS(1..3), DFLUSH, TREEDD per year ---
function _gd_weather()
    wx = Dict{Int,Tuple{NTuple{3,Float32},Float32,Float32}}()
    for ln in readlines(joinpath(_GD_DIR, "bwe_wea_golden.txt"))
        p = split(ln)
        isempty(p) && continue
        yr = parse(Int, p[3])
        v = [_h2f(x) for x in p[5:end]]        # DAYS1 DAYS2 DAYS3 WHOTF DFLUSH TREEDD ...
        wx[yr] = ((v[1], v[2], v[3]), v[5], v[6])
    end
    wx
end

# --- parse the per-cell BWELIT IN/OUT golden into per-year records ---
function _gd_records()
    recs = Any[]
    cur = nothing
    for ln in readlines(joinpath(_GD_DIR, "bwelit_io_golden.txt"))
        p = split(ln)
        isempty(p) && continue
        tag = p[1]
        if tag == "IN"
            cur = Dict{Symbol,Any}()
            cur[:yr] = parse(Int, p[2])
            cur[:EGGS] = _h2f(p[3]); cur[:HOSTST] = _h2f(p[4])
            cur[:FNEW]  = zeros(Float32, 9, 6); cur[:FOLD1] = zeros(Float32, 9, 6)
            cur[:FOLD2] = zeros(Float32, 9, 6); cur[:FREM]  = zeros(Float32, 9, 6)
            cur[:DEFYRS] = zeros(Float32, 6)
        elseif tag == "INF"
            ic = parse(Int, p[2]); ih = parse(Int, p[3])
            cur[:FNEW][ic, ih]  = _h2f(p[4]); cur[:FOLD1][ic, ih] = _h2f(p[5])
            cur[:FOLD2][ic, ih] = _h2f(p[6]); cur[:FREM][ic, ih]  = _h2f(p[7])
            cur[:DEFYRS][ih] = _h2f(p[8])
        elseif tag == "OUT"
            cur[:EGGSout] = _h2f(p[3])
            cur[:oFNEW]  = zeros(Float32, 9, 6); cur[:oFOLD1] = zeros(Float32, 9, 6)
            cur[:oFOLD2] = zeros(Float32, 9, 6); cur[:oFREM]  = zeros(Float32, 9, 6)
        elseif tag == "OUTF"
            ic = parse(Int, p[2]); ih = parse(Int, p[3])
            cur[:oFNEW][ic, ih]  = _h2f(p[4]); cur[:oFOLD1][ic, ih] = _h2f(p[5])
            cur[:oFOLD2][ic, ih] = _h2f(p[6]); cur[:oFREM][ic, ih]  = _h2f(p[7])
            if ic == 9 && ih == 6
                push!(recs, cur); cur = nothing
            end
        end
    end
    recs
end

@testset "WSBWE GENDEFOL/BUDLITE — bwelit.f core bit-exact (dump-replay vs FVSem_wsbwe)" begin
    wx = _gd_weather()
    recs = _gd_records()
    @test !isempty(recs)
    @test Set(r[:yr] for r in recs) == Set([1990, 1991, 1992])

    for c in recs
        (days, dflush, treedd) = wx[c[:yr]]
        FNEW  = copy(c[:FNEW]);  FOLD1 = copy(c[:FOLD1])
        FOLD2 = copy(c[:FOLD2]); FREM  = copy(c[:FREM]); DEF = copy(c[:DEFYRS])
        (eo, _dl, _ran) = _GDW.wsbwe_bwelit!(FNEW, FOLD1, FOLD2, FREM,
                                             c[:EGGS], c[:HOSTST], treedd, dflush, days, DEF)
        # carried eggs bit-exact
        @test _u32(eo) == _u32(c[:EGGSout])
        # every foliage cell bit-exact (Float32-hex)
        mism = 0
        for ic in 1:9, ih in 1:6
            for (m, o) in ((FNEW, c[:oFNEW]), (FOLD1, c[:oFOLD1]),
                           (FOLD2, c[:oFOLD2]), (FREM, c[:oFREM]))
                _u32(m[ic, ih]) == _u32(o[ic, ih]) || (mism += 1)
            end
        end
        @test mism == 0
    end
end
