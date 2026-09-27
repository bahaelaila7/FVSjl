# test_rd_stub.jl — RDIN in a build without the Western Root Disease extension is the exrd.f stub.
#
# FVS{sn,cs,ls,ne,ak,ca,oc,op,on}_buildDir link rd/exrd.f instead of rd/rdin.f: RDIN raises ERRGRO(.TRUE.,11)
# ("REQUESTED EXTENSION IS NOT PART OF THIS PROGRAM") and nothing activates (live SN 216786838010854 rootdis: FVS11 +
# four FVS01 rows, .sum = the no-RD run). jl ran the WRD model there with the KT species crosswalk and died with a
# SIGSEGV in rd_iprp! (SN tiered rootdis regime, every stand). Here the RDIN block must leave the run unchanged.
using FVSjl, Test

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
end
