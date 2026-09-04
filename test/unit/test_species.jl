# C2 unit tests — Southern species tables & BLOCK DATA defaults.
# (Full 90-species equality with Oracle A was verified at port time; these pin the
#  shape, anchor values, and the init_blockdata! application as a regression guard.)

using Test
using FVSjl
using FVSjl: SN_NSPECIES, init_blockdata!, DEFAULT_TREE_FORMAT, coefficients

@testset "species table shapes & anchors (loaded from CSV)" begin
    c = coefficients(Southern())
    @test SN_NSPECIES == 90
    @test length(c.code_alpha) == MAXSP
    @test length(c.code_fia) == MAXSP
    @test length(c.code_plants) == MAXSP
    @test length(c.species[:dg_resid_sd]) == MAXSP   # numeric coeff vectors padded to capacity
    @test length(c.valid_habitat) == 122
    # anchors (first / mid / last species)
    @test strip(c.code_alpha[1]) == "FR"  && c.code_fia[1] == "010" && strip(c.code_plants[1]) == "ABIES"
    @test strip(c.code_alpha[13]) == "LP" && c.code_fia[13] == "131"   # loblolly pine
    @test strip(c.code_alpha[90]) == "OT" && c.code_fia[90] == "999"
    @test c.species[:dg_resid_sd][1] == 0.4511f0
    @test c.valid_habitat[95] == 999 && c.valid_habitat[96] == 0
end

@testset "species resolution (direct + SPCTRN crosswalk)" begin
    s = StandState(Southern()); init_blockdata!(s, s.variant)
    sp, v, co = s.species, s.variant, s.coef
    # direct matches against the variant's own codes (format 1=alpha,2=FIA,3=PLANTS)
    @test resolve_species("FR", v, sp, co)    == (Int32(1), Int32(1))
    @test resolve_species("010", v, sp, co)   == (Int32(1), Int32(2))
    @test resolve_species("ABIES", v, sp, co) == (Int32(1), Int32(3))
    @test resolve_species("LP", v, sp, co)    == (Int32(13), Int32(1))   # loblolly pine
    @test resolve_species("131", v, sp, co)   == (Int32(13), Int32(2))
    # SPCTRN crosswalk: "BF" (balsam fir, not an SN species) → SN target FR → idx 1
    @test resolve_species("BF", v, sp, co)[1]    == Int32(1)
    # PLANTS "2TREE" → SN target OT → idx 90 (catch-all other)
    @test resolve_species("2TREE", v, sp, co)[1] == Int32(90)
    # blank code → OT
    @test resolve_species("   ", v, sp, co)[1]   == Int32(90)
    # crosswalk table shape (loaded from CSV)
    @test length(co.translation) == 562
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

@testset "IE/EM western SPCTRN crosswalk (ie/spctrn.f ASPT: IE=col11, EM=col10)" begin
    # Regression guard for the IE crosswalk-column bug: the IE species_translation.csv had been
    # derived from KT's ASPT column (target_kt), which collapses ~274 unrecognized species to "OT"
    # → other_species(23)=OS. IE's OWN ASPT column (11) maps hardwoods to OH(22) etc. The wrong table
    # starved translated-species stands (e.g. bitter cherry 768) of height/diameter growth (OS grows
    # far slower than OH), a one-directional BA/QMD under-build in the FIA establishment/dense regime.
    for (v, other) in ((InlandEmpire(), Int32(23)), (EasternMontana(), Int32(19)))
        s = StandState(v); init_blockdata!(s, s.variant)
        sp, var, co = s.species, s.variant, s.coef
        # bitter cherry (FIA 768, Prunus emarginata) → OH, NOT the OS/OH catch-all-by-accident
        idx = resolve_species("768", var, sp, co)[1]
        @test strip(sp.alpha[idx]) == "OH"
        # generic hardwood 998 → OH; generic softwood 999 → OS; oak 850 → OH
        @test strip(sp.alpha[resolve_species("998", var, sp, co)[1]]) == "OH"
        @test strip(sp.alpha[resolve_species("999", var, sp, co)[1]]) == "OS"
        @test strip(sp.alpha[resolve_species("850", var, sp, co)[1]]) == "OH"
        # ABAM (FIA 011) → OS (no true fir in IE/EM 23/19-species set)
        @test strip(sp.alpha[resolve_species("011", var, sp, co)[1]]) == "OS"
        @test length(co.translation) == 442
        @test FVSjl.other_species(var) == other
    end
end
