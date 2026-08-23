# =============================================================================
# test_ec_hbdecd_estab.jl — two EC (East Cascades) per-stand fixes surfaced by the
# FVSppe historical-oracle landscape run (reproducing stands: the ecpp.key S248112
# stands, STDINFO forest 606 / plant-association CWS222, + a NOTREES/PLANT stand).
#
# BUG 1 (habitat/plant-association → site resolution): EC STDINFO field 2 may be an
#   ALPHA plant-association code (e.g. "CWS222"), decoded by ec/habtyp.f→hbdecd.f into
#   the KODTYP index that ec_sitset! reads. The keyword reader was parsing the alpha
#   field to habitat_code=0 ⇒ the stand grew on the poor CPS241 default site (SI 0) instead
#   of its real GF / SI≈103 site. ec_hbdecd now decodes it; validated vs FVSec_g16
#   (habitat 129, site species GF(6), SITEAR fan matches the oracle SITECODE dump to NINT).
#
# BUG 2 (ESTAB crash): a PLANT/regen EC stand hit `KeyError: :estab_min_ht` — EC was
#   absent from the establishment XMIN/HHTMAX/essubh dispatch and lacks that shared coef
#   column. Wired EC's ec/blkdat.f XMIN/HHTMAX + ec/essubh.f (SMHTGF MODE=0) base height +
#   ec/estab.f:626 DBH=0.1; the stand now runs (the plantation establishes + grows).
# =============================================================================
using Test
using FVSjl
const F = FVSjl

@testset "EC ec_hbdecd — STDINFO plant-association decode (ec/hbdecd.f)" begin
    # Alpha code → KODTYP index into EC_PCOML (array2 = 0 because the field isn't numeric).
    @test F.ec_hbdecd("CWS222", 0) == 129
    @test F.ec_hbdecd("CWS222  ", 0) == 129            # trailing blanks ignored
    @test F.ec_hbdecd("cws222", 0) == 129              # UPCASE before match
    @test F.ec_hbdecd("CAG112", 0) == 1                # first PCOML entry
    @test F.ec_hbdecd("HQS211", 0) == 155              # last PCOML entry
    # Numeric sequence number is used directly (1..NPA); out of range / unknown → 0 default.
    @test F.ec_hbdecd("", 129) == 129
    @test F.ec_hbdecd("", 156) == 0                    # > NPA (155) → default
    @test F.ec_hbdecd("ZZZ999", 0) == 0                # no alpha match → default (ITYPE=114)
    @test F.ec_hbdecd("0", 0) == 0                     # leading '0' → HBDECD 'DEFAULT'
    @test F.ec_hbdecd("", 0) == 0                      # blank → default
end

@testset "EC site-index resolution for CWS222 (habitat 129) vs FVSec_g16" begin
    s = F.StandState(F.EastCascades())
    p = s.plot
    p.user_forest_code = Int32(606)
    p.habitat_code = Int32(129)          # CWS222 (as ec_hbdecd would resolve it)
    p.site_species = Int32(0)
    fill!(p.sp_site_index, 0f0); fill!(p.sp_sdi_def, 0f0)
    F.ec_site_index_setup!(s)
    si = p.sp_site_index

    @test Int(p.forest_idx) == 1         # forkod: 606 → IFOR 1 (Mt Hood)
    @test Int(p.site_species) == 6       # ecocls CWS222 site species = GF
    # SITEAR fan matches the FVSec_g16 SITECODE dump (rounded to the oracle's NINT print).
    @test round(Int, si[6])  == 103      # GF (site species)
    @test round(Int, si[1])  == 110      # WP
    @test round(Int, si[2])  == 103      # WL
    @test round(Int, si[3])  == 103      # DF
    @test round(Int, si[10]) == 148      # PP
    @test round(Int, si[23]) == 223      # PB
    @test [round(Int, si[i]) for i in 1:12] ==
          [110,103,103,103,150,103,150,148,148,148,103,28]
end

@testset "EC STDINFO alpha field decodes through the keyword reader" begin
    key = tempname() * ".key"
    open(key, "w") do io
        print(io, """STDIDENT
ECHAB    alpha plant-association STDINFO
NOAUTOES
STDINFO        606.0    CWS222      60.0     315.0      30.0      34.0
INVYEAR       1990.0
NUMCYCLE         1.0
NOTREES
PROCESS
STOP
""")
    end
    s = first(F.each_stand(key; variant = F.EastCascades()))
    F.notre!(s); F.setup_growth!(s)
    @test Int(s.plot.habitat_code) == 129        # "CWS222" decoded, not parsed to 0
    @test Int(s.plot.site_species) == 6          # GF
    @test round(Int, s.plot.sp_site_index[6]) == 103
    rm(key; force = true)
end

@testset "EC PLANT establishment stand runs (no :estab_min_ht crash)" begin
    key = tempname() * ".key"
    open(key, "w") do io
        print(io, """STDIDENT
ECPLANT  NOTREES + PLANT regen stand
NOTRIPLE
NOAUTOES
NOTREES
INVYEAR         1990
ESTAB           1992
PLANT           1992         2       400
PLANT           1992        10       400
END
NUMCYCLE           5
PROCESS
STOP
""")
    end
    # BUG 2 was a hard crash here (KeyError :estab_min_ht). It must now run to completion.
    local sumtxt
    @test (sumtxt = F.run_keyfile(key; variant = F.EastCascades())) isa AbstractString
    @test !isempty(sumtxt)
    # The plantation (400 WL + 400 PP) establishes and appears in the summary.
    @test occursin("800", sumtxt)
    rm(key; force = true)
end
