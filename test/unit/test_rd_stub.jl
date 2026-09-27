# test_rd_stub.jl — RDIN in a build without the Western Root Disease extension is the exrd.f stub.
#
# FVS{sn,cs,ls,ne,ak,ca,oc,op,on}_buildDir link rd/exrd.f instead of rd/rdin.f: RDIN raises ERRGRO(.TRUE.,11)
# ("REQUESTED EXTENSION IS NOT PART OF THIS PROGRAM") and nothing activates (live SN 216786838010854 rootdis: FVS11 +
# four FVS01 rows, .sum = the no-RD run). jl ran the WRD model there with the KT species crosswalk and died with a
# SIGSEGV in rd_iprp! (SN tiered rootdis regime, every stand). Here the RDIN block must leave the run unchanged.
using FVSjl, Test, SQLite, DBInterface

const _RDSTUB_SN = joinpath(@__DIR__, "..", "fixtures", "tiered", "sn")

@testset "RDIN is the exrd.f stub where the WRD extension is not linked" begin
    @test !FVSjl.rd_extension_linked(FVSjl.Southern())
    @test !FVSjl.rd_extension_linked(FVSjl.CentralCalifornia())
    @test FVSjl.rd_extension_linked(FVSjl.InlandEmpire())
    @test FVSjl.rd_extension_linked(FVSjl.WestSierra())

    key = read(joinpath(_RDSTUB_SN, "216786838010854_rootdis.key"), String)
    nord = replace(key, r"RDIN\n.*?\nEnd\n"s => "")
    @test nord != key
    dir = mktempdir()
    cp(joinpath(_RDSTUB_SN, "stands.db"), joinpath(dir, "stands.db"))
    write(joinpath(dir, "rd.key"), key); write(joinpath(dir, "nord.key"), nord)
    rd   = cd(() -> FVSjl.run_keyfile("rd.key";   variant = FVSjl.Southern(), output = :sum), dir)
    base = cd(() -> FVSjl.run_keyfile("nord.key"; variant = FVSjl.Southern(), output = :sum), dir)
    rows(t) = filter(l -> !startswith(l, "-999"), split(strip(t), '\n'))
    @test !isempty(rows(rd)) && rows(rd) == rows(base)
    # errgro.f → dbserror.f: the FVS11 from the stub + FVS01 for RRTYPE/SAREA/RRINIT/End (RECORDS READ 15-18: the
    # DSNout/DSNin filename lines are read without IRECNT+1), each Message the 256-char CMSG — equal to live.
    db = SQLite.DB(joinpath(dir, "out.db"))
    got = sort([(String(r.StandID), String(r.Message)) for r in DBInterface.execute(db, "SELECT * FROM FVS_Error")])
    SQLite.close(db)
    live = sort([(String(m[1]), String(m[2])) for m in (match(r"^([^,]*),(.*)$", l) for l in
                 readlines(joinpath(_RDSTUB_SN, "216786838010854_rootdis.FVS_Error.csv"))[2:end])])
    @test length(live) == 5 && got == live
end
