# test_bc_port_live.jl — BC (British Columbia, metric) fixes vs the LIVE BC tiered goldens (test/fixtures/tiered/bc: stands.db +
# <stand>_<regime>.key + the .sum / *_Metric DBS tables of the private DB-capable oracle /workspace/.bcwork/dbfix/FVSbc_dbfix,
# ORACLE_SOURCE_AUDIT §7). Each testset names the Fortran it follows and the case it was measured on.
module BcPortLiveTest
using Test
using FVSjl
include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))

function _case(cn, r)
    d = mktempdir()
    txt, db, crashed, err = run_case("BC", cn, r; dir = d)
    return (txt = txt, db = db, crashed = crashed, err = err,
            ms = crashed ? Mismatch[] : compare_case("BC", cn, r, txt, db))
end
_in(ms, files) = [m for m in ms if m.file in files]

# BC FFE (canada/fire/bc fm*.f + the shared fire/base, fire/vbase, fire/ie layers it links, metric/fire I/O): SIMFIRE and
# SALVAGE ran to an "FFE fuel tables not ported" error. Ported: bc/fmvinit.f species/decay/snag parameters, bc/fmcba.f (the
# IE FMCBA on 15-species FULIVE/FULIVI/FUINIE/FUINII), bc/fmbrkt.f, bc/fmcrow.f ISPMAP (+FMCROWE for 11/12/13/15),
# bc/fmsvol.f (Kozak CFVOL + CFTOPK, no cone floor), metric fmin.f SIMFIRE km/h·°C and SALVAGE cm, ICMETRC=1, and the
# metric DBS (FVS_Carbon_Metric, FVS_Mortality_Metric; no BurnReport/PotFire — metric fmfout.f:96 / fmpofl.f:295 windows).
@testset "BC FFE SIMFIRE/SALVAGE (canada/fire/bc, metric/fire)" begin
    for r in ("simfire", "salvage")
        c = _case("YSM029-250", r)
        @test !c.crashed
        @test isempty(_in(c.ms, ("sum", "FVS_Summary_Metric", "FVS_Mortality_Metric", "FVS_Hrv_Carbon")))
        @test !any(m -> m.col == "PRESENCE", c.ms)          # the live table set: no PotFire/BurnReport/InvReference/Error
    end
end

# canada/bc/cfvol.f → min.f → log.f (Kozak 2002 taper), as FMSVOL calls it: LOG's REAL**REAL/ALOG/EXP are glibc
# powf/logf/expf; MIN converts CFVOL's imperial grinit.f thresholds back to metric (SH = 30·CMtoFT·FTtoM, TD =
# 10·CMtoIN·INtoCM — not 0.3 m / 10 cm, so a sub-0.3 m stump takes log.f's partial-first-log branch and DI3 there feeds the
# stump volume); cfvol.f:55 gates merch in INCHES (D ≥ DBHMIN = 12.5·CMtoIN catches the 12.5-cm PL). Replays 499 FMSVOL calls
# dumped from the private oracle (test/fixtures/bc_ffe/fmsvol_calls.txt): VN and VM bit-exact.
@testset "BC Kozak taper volume via FMSVOL (canada/bc/cfvol.f, min.f, log.f; canada/fire/bc/fmsvol.f)" begin
    h2f(x) = reinterpret(Float32, parse(UInt32, x; base = 16))
    n = 0; bad = 0
    for l in eachline(joinpath(@__DIR__, "..", "fixtures", "bc_ffe", "fmsvol_calls.txt"))
        startswith(l, "#") && continue
        f = split(l)
        sp = parse(Int, f[1]); d = h2f(f[2]); h = h2f(f[3]); xht = h2f(f[4])
        ltkil = f[8] == "T"
        vn, vm = FVSjl.bc_fmsvol(nothing, sp, d, h, ltkil ? xht : -1f0)
        n += 1
        (vn == h2f(f[6]) && vm == h2f(f[7])) || (bad += 1)
    end
    @test n == 499
    @test bad == 0
    c = _case("YSM029-250", "salvage")                     # live FVS_Carbon_Metric pools now exact too
    @test isempty(c.ms)
end

end # module
