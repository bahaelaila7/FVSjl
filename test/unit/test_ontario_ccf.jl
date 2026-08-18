# =============================================================================
# test_ontario_ccf.jl — ON open-grown crown width (cwcalc.f IWHO=1) → CCF (ccfcal.f).
#
# Validated bit-exact vs an instrumented FVSon_g16 (scratchpad/on/instr_cw: cwcalc.f dumps
# ISPC / D / CW / HI-inputs when IWHO=1). All 8 ont01 species reproduce the oracle's open-grown
# crown width to the Float32 bit; the beech (ISPC 28) Hopkins term also nails the oracle's
# TLAT=46.78 / TLONG=92.11 / ELEV=300 m→984.25 ft location (forkod.f US-Superior default +
# metric STDINFO elevation). stand_ccf then overflows the sumout.f I4 field ⇒ the .sum renders
# `****`, matching the oracle byte-for-byte.
# =============================================================================
using FVSjl, Test
const _F = FVSjl
_h2f(h) = reinterpret(Float32, parse(UInt32, h; base=16))
_f2h(x::Float32) = uppercase(string(reinterpret(UInt32, x); base=16, pad=8))

@testset "ON open-grown crown width (cwcalc.f IWHO=1) bit-exact vs FVSon_g16" begin
    # oracle FVSon_cwdump first-cycle dump: ISPC => (D hex, CW hex)
    orc = [(1,"413060AA","41C29F7E"), (5,"413CF9DB","41C47486"), (6,"411D7AE1","4192C3CD"),
           (8,"40FBF7CF","416C9717"), (9,"40BCF9DB","413E9300"), (11,"410A9518","41594FAA"),
           (26,"415C78D5","41F54340"), (28,"416F5E9E","41EB0455")]
    # ON default location (forkod.f US-Superior 915/916): lat 46.78, long 92.11, elev 300 m → 9.8425 (100 ft)
    lat, long, elev = 46.78f0, 92.11f0, 300f0 * 3.28084f0 / 100f0
    for (ispc, dh, cwh) in orc
        cw = _F.on_open_crown_width(ispc, _h2f(dh), lat, long, elev)
        @test _f2h(cw) == cwh
    end
end

@testset "ON stand CCF overflows the sumout I4 field (renders ****)" begin
    for s in _F.each_stand("scratchpad/on/ont01.key"; variant=_F.Ontario(), faithful=true)
        _F.notre!(s); _F.setup_growth!(s)
        # forkod default location wired through kw_stdinfo! (metric elevation + US-Superior lat/long)
        @test s.plot.latitude  == 46.78f0
        @test s.plot.longitude == 92.11f0
        @test _f2h(s.plot.elevation) == _f2h(300f0 * 3.28084f0 / 100f0)
        ccf = _F.stand_ccf(s)
        @test ccf > 9999f0                      # ⇒ I4 overflow ⇒ `****`
        break
    end
end
