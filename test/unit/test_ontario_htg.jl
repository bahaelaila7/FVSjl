# test_ontario_htg.jl — ON large-tree height growth (canada/on/htgf.f + htont.f) end-to-end vs the
# LIVE oracle. Golden = an instrumented FVSon_g16 (htgf.f, unit 773) on ont01 cyc0→cyc1: per tree
#   I ISPC | DBH DG RMSQD BA QMD10 BA10 HTNOW HT10 DBH10 HTG SITEAR PROB BARK HT  (14 Float32-hex).
# Asserts the SHIPPED height_growth!(::Ontario), fed the oracle's exact stand state, reproduces the
# oracle per-tree HTG bit-exact (Float32-hex) — including HTONT (Penner diameter-height) at both the
# current and grown diameters, the internal BA10/QMD10 rebuild (record order, 0.0054542), the 0.1
# floor and the SCALE·XHMULT·HCOR2 folding. This isolates the height-growth chunk from the upstream
# DG (which carries the DGSD=2.0 OLDRN #206 straddle shared by every western variant).

using Test

@testset "ON — large-tree height growth htgf.f/htont.f (ont01 cyc0→cyc1) bit-exact vs FVSon_g16" begin
    _h2f(h) = reinterpret(Float32, parse(UInt32, h; base=16))
    _u(x::Float32) = reinterpret(UInt32, x)

    # golden fort.773 rows (I ISPC + 14 hex), instrumented FVSon_g16 on ont01
    golden = split.(strip.([
      "      7      1 413060AA 3CBF4E00 4113EE12 467A0001 4114974D 467C3DD2 452D85A6 452F264B 4130C435 41D05280 42700000 453C98F2 3F75FD8B 42795806",
      "      1      5 413CF9DB 3CBEDE00 4113EE12 467A0001 4114974D 467C3DD2 44874DAE 448898F2 413D6198 4125A200 42700000 45244A10 3F6B80AA 42833BCD",
      "      2      6 411D7AE1 3CC6A800 4113EE12 467A0001 4114974D 467C3DD2 44478E98 44499BDC 411DE2E3 41035100 42700000 456C939D 3F747AE1 426C3871",
      "      8      8 40FBF7CF 3E059400 4113EE12 467A0001 4114974D 467C3DD2 4488B252 448B5BC3 410030F5 41AA5C40 42700000 45B8D352 3F7212D7 4251F948",
      "      5      9 40BCF9DB 3CCC7700 4113EE12 467A0001 4114974D 467C3DD2 44447F19 4446E26A 40BDD1F3 4118D440 42700000 46244A10 3F7239AF 421D7AF6",
      "      6     11 410A9518 3CC99200 4113EE12 467A0001 4114974D 467C3DD2 4422D1AE 442481C8 410B043C 40D80D00 42700000 4598BF93 3F6825F8 4244D9B4",
      "      3     26 415C78D5 3CC33C00 4113EE12 467A0001 4114974D 467C3DD2 44820E74 44833476 415CE1D2 41130100 42700000 44F1679A 3F6E0766 42905B62",
      "      4     28 416F5E9E 3CC8A200 4113EE12 467A0001 4114974D 467C3DD2 44C1B01A 44C36D3F 416FC956 415E9280 42700000 44CCCAF3 3F70A3D7 429D7AF6",
    ]))
    rows = map(golden) do f
        (; I=parse(Int,f[1]), ispc=parse(Int,f[2]),
           dbh=_h2f(f[3]), dg=_h2f(f[4]), rmsqd=_h2f(f[5]), ba=_h2f(f[6]),
           qmd10=_h2f(f[7]), ba10=_h2f(f[8]), htnow=_h2f(f[9]), ht10=_h2f(f[10]),
           dbh10=_h2f(f[11]), htg=_h2f(f[12]), sitear=_h2f(f[13]), prob=_h2f(f[14]),
           bark=_h2f(f[15]), ht=_h2f(f[16]))
    end

    # --- on_htont (Penner diameter-height): HTNOW & HT10 bit-exact per tree ---
    sim = rows[1].sitear * FVSjl.ON_FTtoM         # SIM = SITEAR(ISISP)·FTtoM
    for r in rows
        htnow = FVSjl.on_htont(r.ispc, r.dbh,   r.rmsqd, r.ba,   sim)
        ht10  = FVSjl.on_htont(r.ispc, r.dbh10, r.qmd10, r.ba10, sim)
        @test _u(htnow) == _u(r.htnow)
        @test _u(ht10)  == _u(r.ht10)
    end

    # --- SHIPPED height_growth!(::Ontario) on the oracle stand state: per-tree HTG bit-exact ---
    s = FVSjl.StandState(FVSjl.Ontario())
    FVSjl.load_species_coefficients!(s, s.variant)
    t=s.trees; p=s.plot; ctl=s.control
    t.n = maximum(r->r.I, rows)
    p.basal_area = rows[1].ba; p.qmd = rows[1].rmsqd
    p.site_species = Int32(3)                      # ISISP=3 (sitset)
    p.sp_site_index[3] = rows[1].sitear            # SITEAR(3)=60 ft
    for sp in 1:FVSjl.MAXSP; ctl.htg_cor2[sp]=1f0; ctl.sp_size_cap[sp,4]=999f0; end
    for r in rows
        t.species[r.I]=Int32(r.ispc); t.dbh[r.I]=r.dbh; t.height[r.I]=r.ht
        t.diam_growth[r.I]=r.dg; t.tpa[r.I]=r.prob
    end
    FVSjl.height_growth!(s, FVSjl.Ontario(); scale=1f0)
    nok = 0
    for r in rows
        @test _u(t.ht_growth[r.I]) == _u(r.htg)
        nok += _u(t.ht_growth[r.I]) == _u(r.htg)
    end
    @test nok == length(rows)
end
