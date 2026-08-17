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

# --- RANSTART: MPOTPR outbreak probability (mpotpr.f) + the MPRANN inclusion draw.
# Goldens are the g16 instrumented FVSie_lpmpb dump for lp_ran.key (MPB/RANSTART 1/DEBUG),
# 5 cycles on the IE lodgepole host. Each cycle: MPOTPR reads stand stats → PROTBK, then
# one MPRANN draw. RELDSP(IDXLP)==RELDEN on this pure-lodgepole stand. ---
@testset "LPMPB — RANSTART MPOTPR probability (bit-exact vs FVSie_lpmpb)" begin
    # (pbalpp, relden, reldsp_lp, a45dbh, cntlp, istdt, icyc, iy, PROTBK-golden)
    mpo = [
        ("3F7FFFFE","428794CB","428794CB","4135B128","42C2FBF0",1,1,1990,"3C10C78E"),
        ("3F800000","429413C2","429413C2","414CC9B9","42AD8EC1",1,2,2000,"3C61867A"),
        ("3F800000","42997F34","42997F34","415FA1B7","429B2614",1,3,2010,"3C888E27"),
        ("3F800000","42993A03","42993A03","41726CCA","4286FC50",1,4,2020,"3C87427C"),
        ("3F800001","4295C866","4295C866","418179C1","426BE526",1,5,2030,"3C6F8B89"),
    ]
    for (pb,rd,rl,a4,cn,istdt,icyc,iy,g) in mpo
        p = _F.mpb_mpotpr(_f(pb), _f(rd), _f(rl), _f(a4), _f(cn), istdt, icyc, iy)
        @test _hex(p) == g
    end
    # minimum-condition gates → PROTBK = 0
    @test _F.mpb_mpotpr(_f("3F800000"), _f("428794CB"), _f("428794CB"), 5.5f0,  100.0f0, 1, 1, 1990) == 0f0  # A45DBH<6
    @test _F.mpb_mpotpr(0.20f0,         _f("428794CB"), _f("428794CB"), 8.0f0,  100.0f0, 1, 1, 1990) == 0f0  # PBALPP<0.25
    @test _F.mpb_mpotpr(_f("3F800000"), 100.0f0,        10.0f0,         8.0f0,  100.0f0, 1, 1, 1990) == 0f0  # RELDSP/RELDEN<0.2
    @test _F.mpb_mpotpr(_f("3F800000"), _f("428794CB"), _f("428794CB"), 8.0f0,   30.0f0, 1, 1, 1990) == 0f0  # CNTLP<40
    @test _F.mpb_mpotpr(_f("3F800000"), _f("428794CB"), _f("428794CB"), 8.0f0,  100.0f0, 5, 1, 1990) == 0f0  # ICYC<ISTDT
end

# MPRANN (mprann.f) — the LPMPB LCG seeded 55329, drawn once per RANSTART-eligible cycle.
@testset "LPMPB — MPRANN draw sequence (seed 55329)" begin
    m = _F.mpb_defaults!()
    for g in ("3EDDB57A","3F5AB21B","3F630A07","3F2734C9","3EF4FB2C")
        @test _hex(_F.mpb_rand!(m)) == g
    end
end

# CURRMORT / INVMORT ICYC=1 GREINF branch (colmod.f:93-99). Golden = instrumented
# FVSie_lpmpb dump for lp_curr.key (CURRMORT classes 3-7 = 2,5,10,3,1). GREINF=0 on the
# loadable stand (no treelist damage codes); CURRMR is the keyword. Classes 4-7 diverge
# from the default epidemic; classes 8-9 (CURRMR=0) stay at the default GREEN/PRKILL.
@testset "LPMPB — CURRMORT GREINF branch (bit-exact vs FVSie_lpmpb)" begin
    st   = _F.mpb_coldbh_start(_LP_DBH, _LP_TPA)
    gdef = _F.mpb_colmod(st, _NUMYRS)                            # default epidemic (no CURRMORT)
    currmr = Float32[0,0,2,5,10,3,1,0,0,0]
    gc = _F.mpb_colmod(st, _NUMYRS; icyc=1, lcurmr=true, currmr=currmr)
    pc = _F.mpb_prkill(st, @view gc[_NUMYRS, :])
    Gc = Dict(4=>"4123BE78",5=>"41250985",6=>"4022B3FC",7=>"3DD0597F",8=>"36C6E3D7",9=>"3DCF7CEC")
    Pc = Dict(4=>"3DDAF637",5=>"3F2F6F66",6=>"3F62D378",7=>"3F7E6380",8=>"3F7FFFF8",9=>"3F75C794")
    for c in 4:9
        @test _hex(gc[_NUMYRS, c]) == Gc[c]
        @test _hex(pc[c]) == Pc[c]
    end
    # GREINF=0 AND CURRMR=0 (INVMORT with no inventory damage) ⇒ byte-identical to default.
    g0 = _F.mpb_colmod(st, _NUMYRS; icyc=1, linvmr=true)         # currmr defaults to zeros
    for c in 1:10
        @test _hex(g0[_NUMYRS, c]) == _hex(gdef[_NUMYRS, c])
    end
    # ICYC≠1 ⇒ default branch even with LCURMR (CURRMORT only auto-fires cycle 1).
    g2 = _F.mpb_colmod(st, _NUMYRS; icyc=2, lcurmr=true, currmr=currmr)
    for c in 1:10
        @test _hex(g2[_NUMYRS, c]) == _hex(gdef[_NUMYRS, c])
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
