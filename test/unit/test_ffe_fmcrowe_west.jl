# test_ffe_fmcrowe_west.jl — FFE FMCROWE (fmcrowe.f, the Jenkins/proportion crown for the hardwoods) on a WESTERN variant
# vs the live FVSem_g16 CROWNW (private FMDOUT trace of 684750664126144 SIMFIRE, 2018 record 28).
module FfeFmcroweWestTest
using Test
using FVSjl

# fmcrowe.f:263-266 adds the small-tree bole weight to TTOPW only for VARACD ∈ {SN,LS,NE,CS,ON}, and :396-410 the LILPCE
# correction likewise; SG is the runtime V2T = raw/2000 (fmvinit.f), and the arithmetic is REAL*4 EXP/LOG/**/TAN/ATAN
# with /12/12/12. MEASURED EM aspen 5.9"×31' CR 20 under DBHMIN 7: live CROWNW(0:4) 3.87054/3.594432/8.53118/8.376474/
# 32.608471; jl added the bole (TTOPW 55.7291 vs 53.1106) ⇒ 3.87054/4.193983/9.924952/9.001656/32.608475.
@testset "FMCROWE on a western variant: no eastern bole/LILPCE terms, REAL*4 arithmetic (fmcrowe.f) vs FVSem_g16" begin
    s = FVSjl.StandState(FVSjl.EasternMontana()); FVSjl.init_blockdata!(s, s.variant)
    s.control.merch_init = true; resize!(s.control.sp_dbh_min, 19); fill!(s.control.sp_dbh_min, 7f0)
    xv = FVSjl.crown_biomass(s, 12, 5.9f0, 31f0, 20)
    live = (0x4077B6EE, 0x40660B2C, 0x41087FB7, 0x4106060A, 0x42026F13, 0x00000000)
    @test all(reinterpret(UInt32, xv[k]) == live[k] for k in 1:6)
end

include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))
# End to end: the stand's live carbon rows once the post-fire aspen regrowth dominates (2038, 2058; 2048 still carries a
# 1-ULP-per-record FMCROWW conifer crown residual — Julia exp/log/^ in cr_crownw — 3.6e-6 relative)
@testset "EM 684750664126144 SIMFIRE live carbon with the western FMCROWE vs FVSem_g16" begin
    d = mktempdir()
    txt, db, crashed, _ = run_case("EM", "684750664126144", "simfire"; dir = d)
    ms = compare_case("EM", "684750664126144", "simfire", txt, db)
    @test !crashed
    @test !any(m -> m.file == "FVS_Carbon" && m.col == "Aboveground_Total_Live" && m.year in ("2038", "2058") &&
               abs(parse(Float64, m.got) - parse(Float64, m.gold)) > 1e-6 * abs(parse(Float64, m.gold)), ms)
end
end # module
