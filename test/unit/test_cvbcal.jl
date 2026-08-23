# COVER shrub CALIBRATION (covr/cvbcal.f) — bit-exact regression for the ported
# `cvbcal!` correction-factor kernel and its `_cvbcal_apply` application.
#
# CVBCAL computes per-shrub-species height/cover correction factors BHTCF(31)/BPCCF(31)
# from user-observed shrub data (SHRBLAYR keyword → method 1 "by layer"; SHRUBHT/SHRUBPC
# → method 2 "by species"), which cvbrow.f then multiplies into the predicted SH/PBCV
# arrays.  This is a report-only path (DGSD=0): no tree record, no RNG, no growth.
#
# The test replays the ported kernel against golden dumps captured from an instrumented
# FVSem_g16 (single-.o swap of cvbcal.f + cvbrow.f, Float32-hex TRANSFER/Z8.8 dump) on the
# EM S248112 stand for both calibration methods.  For each method the fixture carries the
# oracle's pre-calibration predictions (CVBCAL_X{SH,CV,PB}) + the height-sorted species
# index (CVBCAL_HTIDX), which are fed into `cvbcal!`; the resulting BHTCF/BPCCF and the
# calibrated cycle-0 SH/PB/PBCV (via `_cvbcal_apply`) must be Float32-hex byte-identical to
# the oracle's CVBCAL_B{HTCF,PCCF} and CVBROW_{SH,PB,PBCV}.
#
# Feeding the oracle's own pre-cal arrays isolates the calibration kernel exactly (the
# end-to-end residual is a pre-existing 1-ULP upstream stand-BA difference; see
# scratchpad/cover/cvbcal/staged/VALIDATION.md).  GOLDEN PROVENANCE:
# test/fixtures/cover/cvbcal/oracle_m{1,2}_dump.txt.
#
# Observed inputs mirror the keyfiles em_cal{1,2}.key:
#   method 1 (SHRBLAYR): layers (4.0ft,30%),(2.0ft,20%),(0.5ft,10%) → SUMCVR=60
#   method 2 (SHRUBHT/SHRUBPC): SPBE h1.5/c12, VASC h0.8/c5, PHMA h3.0, ACGL h5.0/c0,
#                               CARX h0.0/c6, XETE c8   (exercises all four apply branches)

using Test
using FVSjl

const _CVBCAL_DIR = joinpath(@__DIR__, "..", "fixtures", "cover", "cvbcal")
_cvbcal_h2f(x::AbstractString) = reinterpret(Float32, parse(UInt32, x; base = 16))

# Load a golden dump into tag → (index → value).  CVBROW_* keeps cyc0/ithn1 only.
function _cvbcal_load(path)
    d = Dict{String,Dict{Int,String}}()
    for ln in eachline(path)
        f = split(ln)
        isempty(f) && continue
        tag = f[1]
        if startswith(tag, "CVBCAL_")
            length(f) >= 3 && (get!(d, tag, Dict{Int,String}())[parse(Int, f[2])] = f[3])
        elseif startswith(tag, "CVBROW_") && length(f) >= 5 && f[2] == "0" && f[3] == "1"
            get!(d, tag, Dict{Int,String}())[parse(Int, f[4])] = f[5]
        end
    end
    return d
end

# Build a minimal ShrubCycle carrying the pre-cal arrays so _cvbcal_apply can run.
_cvbcal_mkcyc(xsh, xcv, xpb, htindx) =
    FVSjl.ShrubCycle(0f0, 0f0, copy(xpb), copy(xsh), zeros(Float32, 31), copy(xcv),
                     zeros(Float32, 31), sum(xcv), copy(htindx))

function _cvbcal_run(cal, d)
    xsh = [_cvbcal_h2f(d["CVBCAL_XSH"][i]) for i in 1:31]
    xcv = [_cvbcal_h2f(d["CVBCAL_XCV"][i]) for i in 1:31]
    xpb = [_cvbcal_h2f(d["CVBCAL_XPB"][i]) for i in 1:31]
    htindx = [parse(Int, d["CVBCAL_HTIDX"][i]) for i in 1:31]
    FVSjl.cvbcal!(cal, copy(xsh), copy(xcv), copy(xpb), htindx, sum(xcv))
    calcyc = FVSjl._cvbcal_apply(cal, _cvbcal_mkcyc(xsh, xcv, xpb, htindx))
    return cal, calcyc, d
end

_hexeq(jl::Float32, orc::AbstractString) =
    reinterpret(UInt32, jl) == parse(UInt32, orc; base = 16)

function _cvbcal_check(cal, calcyc, d)
    ok = 0; tot = 0
    for (jlv, tag) in ((cal.bhtcf, "CVBCAL_BHTCF"), (cal.bpccf, "CVBCAL_BPCCF"),
                       (calcyc.sh, "CVBROW_SH"), (calcyc.pb, "CVBROW_PB"),
                       (calcyc.pbcv, "CVBROW_PBCV"))
        for i in 1:31
            tot += 1
            _hexeq(jlv[i], d[tag][i]) && (ok += 1)
        end
    end
    return ok, tot
end

@testset "COVER CVBCAL — method 1 (SHRBLAYR) bit-exact (dump-replay vs FVSem_g16)" begin
    d = _cvbcal_load(joinpath(_CVBCAL_DIR, "oracle_m1_dump.txt"))
    cal = FVSjl.ShrubCalib()
    cal.lcal1 = true; cal.lcalib = true
    cal.avgbht .= Float32[4.0, 2.0, 0.5]
    cal.avgbpc .= Float32[30.0, 20.0, 10.0]
    cal.nklass = 3; cal.sumcvr = 60.0f0
    cal, calcyc, d = _cvbcal_run(cal, d)
    ok, tot = _cvbcal_check(cal, calcyc, d)
    @test tot == 155
    @test ok == 155            # BHTCF+BPCCF+SH+PB+PBCV all Float32-hex byte-identical
end

@testset "COVER CVBCAL — method 2 (SHRUBHT/SHRUBPC) bit-exact (dump-replay vs FVSem_g16)" begin
    d = _cvbcal_load(joinpath(_CVBCAL_DIR, "oracle_m2_dump.txt"))
    cal = FVSjl.ShrubCalib()
    cal.lcal2 = true; cal.lcalib = true
    # SPBE=5 VASC=6 CARX=7 PHMA=10 XETE=17 ACGL=20
    cal.shrbht[5] = 1.5f0; cal.shrbht[6] = 0.8f0; cal.shrbht[10] = 3.0f0
    cal.shrbht[20] = 5.0f0; cal.shrbht[7] = 0.0f0
    cal.shrbpc[5] = 12.0f0; cal.shrbpc[6] = 5.0f0; cal.shrbpc[17] = 8.0f0
    cal.shrbpc[20] = 0.0f0; cal.shrbpc[7] = 6.0f0
    cal, calcyc, d = _cvbcal_run(cal, d)
    ok, tot = _cvbcal_check(cal, calcyc, d)
    @test tot == 155
    @test ok == 155            # incl PB reset + the SHRBPC<0 / SHRBHT=0 not-observed branches
end
