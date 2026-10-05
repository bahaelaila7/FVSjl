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

end # module
