# Mountain Pine Beetle (LPMPB) impact model — the DEFAULT deterministic Cole
# "rate of loss" mortality path (COLIND / COLDBH / COLMOD / COLMRT), MPBER, the
# variant IDXLP crosswalk, plus the inert-seam guarantee.
#
# DOCTRINE: the deterministic routines are validated BIT-EXACT (Float32) against a
# relinked FVSie_lpmpb oracle (IE object set, base/exmpb.f no-op stub replaced by
# real lpmpb/*.o + wsbwe/txnote.o; scratchpad/lpmpb/build_ie_lpmpb.sh). An
# instrumented colmrt.o (g16 single-.o swap) dumps Float32 bit patterns for the IE
# lodgepole host stand lp_on.key (MPBSTART cycle 1, POPDYN off), the instrumented
# .sum data rows byte-identical to the clean relink. The golden hex constants
# embedded here are that dump (scratchpad/lpmpb/lp_dbg.stderr) verbatim.
#
# Golden IE lodgepole stand — 18 LP records (IND1 order); ICYC=1, IFINT=10,
# MPMXYR=10 ⇒ NUMYRS=10. DBH/PROB below; COLDBH bins them into START[1..10].

using Test
using FVSjl
const _F = FVSjl

_hex(x::Float32) = uppercase(string(reinterpret(UInt32, x); base = 16, pad = 8))
_f(s) = reinterpret(Float32, parse(UInt32, s; base = 16))

# --- golden per-record (DBH-hex, PROB-hex, size class, XT-hex) from FVSie_lpmpb ---
const _REC = [
    ("41400000","40A2F984",6,"40883DB2"), ("41600000","406F7916",7,"406D875E"),
    ("41800000","403758B4",8,"403758AE"), ("41200000","40EAAEFB",5,"3FED570F"),
    ("41100000","4110DDCB",5,"40128198"), ("41300000","40C1F40D",6,"40A22357"),
    ("41500000","408ADDB1",7,"4089BD13"), ("41700000","40509B8A",8,"40509B83"),
    ("41000000","413758B4",4,"3E3CF580"), ("41880000","40226927",9,"401BED40"),
    ("41300000","40C1F40D",6,"40A22357"), ("41500000","408ADDB1",7,"4089BD13"),
    ("41100000","4110DDCB",5,"40128198"), ("41700000","40509B8A",8,"40509B83"),
    ("41400000","40A2F984",6,"40883DB2"), ("41200000","40EAAEFB",5,"3FED570F"),
    ("41600000","406F7916",7,"406D875E"), ("41800000","403758B4",8,"403758AE"),
]
# COLDBH START, COLMOD GREEN(NUMYRS=10), COLMRT PRKILL by class 1..10 (0 elsewhere)
const _G_START  = Dict(4=>"413758B4",5=>"42031AA4",6=>"41B276C8",7=>"41814D1E",8=>"4143FA1F",9=>"40226927")
const _G_GREEN  = Dict(4=>"413464DE",5=>"41C3EA00",6=>"406A3223",7=>"3E065EB0",8=>"36C6E3D7",9=>"3DCF7CEC")
const _G_PRKILL = Dict(4=>"3C83EB13",5=>"3E8172ED",6=>"3F5601C1",7=>"3F7DEBEF",8=>"3F7FFFF8",9=>"3F75C794")

const _LP_DBH = Float32[_f(r[1]) for r in _REC]
const _LP_TPA = Float32[_f(r[2]) for r in _REC]
const _NUMYRS = 10

@testset "LPMPB — Cole rate-of-loss deterministic core (bit-exact vs FVSie_lpmpb)" begin
    # COLIND (colind.f): DBH → size class
    @testset "COLIND" begin
        for r in _REC
            @test _F.mpb_colind(_f(r[1])) == r[3]
        end
        @test _F.mpb_colind(0.5f0) == 1     # < 1" → class 1
        @test _F.mpb_colind(25.0f0) == 10   # > 20" → class 10
    end

    # COLDBH (coldbh.f): START[] = LP trees/ac per class
    start = _F.mpb_coldbh_start(_LP_DBH, _LP_TPA)
    @testset "COLDBH START" begin
        for c in 1:10
            g = get(_G_START, c, "00000000")
            @test _hex(start[c]) == g
        end
    end

    # COLMOD (colmod.f): deterministic epidemic loop → GREEN(NUMYRS)
    green = _F.mpb_colmod(start, _NUMYRS)
    @testset "COLMOD GREEN(NUMYRS)" begin
        for c in 1:10
            g = get(_G_GREEN, c, "00000000")
            @test _hex(green[_NUMYRS, c]) == g
        end
    end

    # COLMRT (colmrt.f): PRKILL + per-record XT
    prkill = _F.mpb_prkill(start, @view green[_NUMYRS, :])
    @testset "COLMRT PRKILL" begin
        for c in 1:10
            g = get(_G_PRKILL, c, "00000000")
            @test _hex(prkill[c]) == g
        end
    end
    @testset "COLMRT per-record XT = PRKILL(colind(DBH))·PROB" begin
        for r in _REC
            xt = prkill[_F.mpb_colind(_f(r[1]))] * _f(r[2])
            @test _hex(xt) == r[4]
        end
    end

    # MPBER (mpber.f): minimum-condition flag on the host stand
    @testset "MPBER minimum condition" begin
        er = _F.mpb_er(_LP_DBH, _LP_TPA, 300.0f0)
        @test er.noer                       # host clears (≥1 TPA LP ≥4.5")
        @test _F.mpb_er(Float32[], Float32[], 300.0f0).noer == false
        # < 1 TPA of LP ≥4.5" → fails
        @test _F.mpb_er(Float32[3.0], Float32[0.5], 100.0f0).noer == false
    end

    # IDXLP crosswalk (mpblkd<v>.f)
    @testset "mpb_idxlp" begin
        @test _F.mpb_idxlp(InlandEmpire()) == 7
        @test _F.mpb_idxlp(CentralRockies()) == 11
        @test _F.mpb_idxlp(EasternMontana()) == 7
        @test _F.mpb_idxlp(Teton()) == 7
        @test _F.mpb_idxlp(Utah()) == 7
    end
end

# --- inert-seam guarantee: a stand with no MPB block is unchanged; a stand that
# parses an MPB block but whose variant/host/schedule doesn't fire projects
# byte-identically (mirror the DFB/DFTM/WPBR inert-seam tests). ---
@testset "LPMPB — inert seam" begin
    # (Full .sum-DELTA / run-based inert test belongs in the integration suite once
    #  the engine seam is wired; see scratchpad/lpmpb/HANDOFF.md. The unit checks
    #  above pin the mortality math bit-exact.)
    @test _F.mpb_defaults!().active == false
    @test _F.mpb_defaults!().lpopdy == false
    @test length(_F.mpb_defaults!().zinmor) == 10
end
