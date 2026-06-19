# C2 unit tests — Southern species tables & BLOCK DATA defaults.
# (Full 90-species equality with Oracle A was verified at port time; these pin the
#  shape, anchor values, and the init_blockdata! application as a regression guard.)

using Test
using FVSjl
using FVSjl: SN_ALPHA, SN_FIA, SN_PLANTS, SN_SIGMAR, SN_VALID_HABITAT, SN_NSPECIES,
             init_blockdata!, DEFAULT_TREE_FORMAT

@testset "species table shapes & anchors" begin
    @test SN_NSPECIES == 90
    @test length(SN_ALPHA) == 90
    @test length(SN_FIA) == 90
    @test length(SN_PLANTS) == 90
    @test length(SN_SIGMAR) == 90
    @test length(SN_VALID_HABITAT) == 122
    # anchors (first / mid / last species)
    @test strip(SN_ALPHA[1]) == "FR"   && SN_FIA[1] == "010" && strip(SN_PLANTS[1]) == "ABIES"
    @test strip(SN_ALPHA[13]) == "LP"  && SN_FIA[13] == "131"   # loblolly pine
    @test strip(SN_ALPHA[90]) == "OT"  && SN_FIA[90] == "999"
    @test SN_SIGMAR[1] == 0.4511f0
    @test SN_VALID_HABITAT[95] == 999 && SN_VALID_HABITAT[96] == 0
end

@testset "species resolution (direct + SPCTRN crosswalk)" begin
    s = StandState(Southern()); init_blockdata!(s, s.variant)
    sp, v = s.species, s.variant
    # direct matches against the variant's own codes (format 1=alpha,2=FIA,3=PLANTS)
    @test resolve_species("FR", v, sp)    == (Int32(1), Int32(1))
    @test resolve_species("010", v, sp)   == (Int32(1), Int32(2))
    @test resolve_species("ABIES", v, sp) == (Int32(1), Int32(3))
    @test resolve_species("LP", v, sp)    == (Int32(13), Int32(1))   # loblolly pine
    @test resolve_species("131", v, sp)   == (Int32(13), Int32(2))
    # SPCTRN crosswalk: "BF" (balsam fir, not an SN species) → SN target FR → idx 1
    @test resolve_species("BF", v, sp)[1]    == Int32(1)
    # PLANTS "2TREE" → SN target OT → idx 90 (catch-all other)
    @test resolve_species("2TREE", v, sp)[1] == Int32(90)
    # blank code → OT
    @test resolve_species("   ", v, sp)[1]   == Int32(90)
    # SPCTRN_TABLE shape
    @test length(FVSjl.SPCTRN_TABLE) == 562
    @test FVSjl.spctrn_column(v) == 7
    @test FVSjl.other_species(v) == Int32(90)
end

@testset "init_blockdata! applies SN defaults" begin
    s = StandState(Southern())
    @test s.rng.s0 == 0.0                     # bare state: neutral seed
    init_blockdata!(s, s.variant)
    @test s.rng.s0 == 55329.0                 # blkdat seeds main stream
    @test s.rng.ss == 55329.0f0
    @test s.control.tree_format == DEFAULT_TREE_FORMAT
    @test strip(s.species.alpha[5]) == "SP"   # shortleaf pine
    @test strip(s.species.fia[5]) == "110"
    @test strip(s.species.class_codes[5, 2]) == "SP2"
    @test s.plot.valid_habitat[1] == 10
end
