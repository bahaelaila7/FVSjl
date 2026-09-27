# test_ak_port.jl — AK (SoutheastAlaska) establishment / small-tree growth / volume vs the live FVSak_g16 oracle.
#
# Fixture test/fixtures/southeastalaska/akplant.key = akt01's bare PLANT stand (NOTREES, ESTAB 1992, PLANT 400 WS +
# 400 SS, NOAUTOES, NOTRIPLE, 10 cycles); akplant.live.sum = FVSak_g16 on it.
#  * establishment: ak/estab.f (the tally + the per-regen-plot PLANT trees at the ESSUBH HTCALC height) and
#    ak/esgent.f → ak/regent.f REGENT(.TRUE.) grow the 100 records in their birth cycle, then the cycling
#    REGENT(.FALSE.) small-tree model (HTCALC Hegyi/Payandeh + Chapman-Richards DBH inverse) — every
#    TPA/BA/SDI/CCF/TopHt/QMD row through 2082.
#  * volume start bark: vols.f:150-151 BARK=BRATIO(DBH_start) before D=D+DG/BARK, so NVEL's DBTBH=D·(1−BARK)
#    scales the F32 profile with the start-of-cycle bark — the small-SS TCuFt 0.1 roundings (2012 231 not 240).
#
# Fixture akffe.key/.tre = akt01's "FFE TEST" stand (THINDBH 3" 1993, SNAGINIT, FLAMEADJ + SIMFIRE 2003, SALVAGE,
# DEFULMOD); akffe.live.sum = FVSak_g16. The AK FFE (fire/ak fmvinit/fmcba/fmcfmd/fmbrkt + the FVSpn-identical
# rest) must burn the 2003 fire as live does (fuel model 8, flame 0.81 ft, scorch 1.47 ft) and kill with the AK
# FOFEM bark thickness: 2013 TPA 139 (386 with the fire model off).
#
# akt01.key/.tre (tests/FVSak) + akt01.live.sum: every .sum row of all five stands (unthinned control, THINDBH, the
# shelterwood + ECON, FFE, PLANT) — the full AK growth chain (DGF/DGSCOR calibration, HTGF, REGENT, CROWN, MORTS with
# its PASS density loop, R10 volume) must be exact.
#
# ak_fia.db + akfia.key: a 12-stand AK FIA sub-DB (FVS_STANDINIT_COND/FVS_TREEINIT_COND rows copied from the FIA DB;
# live .sum identical to the full-DB run) with AUTOES on (ingrowth), missing crowns (the LSTART CROWN/DUBSCR dub reads
# the CRATET DENSE's PTBAL point BA: live + inventory-dead at read DBH), HISTORY-8 dead, broken tops. Stands listed in
# `exact` must match every .sum row; the total count guards the rest.
using FVSjl, Test

function _ak_sum_rows(path)
    rows = Dict{Int,Vector{Float64}}()
    for l in eachline(path)
        f = split(l)
        (length(f) > 11 && all(isdigit, f[1]) && length(f[1]) == 4) || continue
        rows[parse(Int, f[1])] = parse.(Float64, f[3:12])      # TPA BA SDI CCF TopHt QMD TCuFt MCuFt SCuFt BdFt
    end
    return rows
end

@testset "AK PLANT stand vs live FVSak (establishment + REGENT + start-bark volume)" begin
    fx = joinpath(@__DIR__, "..", "fixtures", "southeastalaska")
    live = _ak_sum_rows(joinpath(fx, "akplant.live.sum"))
    dir = mktempdir()
    cp(joinpath(fx, "akplant.key"), joinpath(dir, "akplant.key"))
    jl = cd(dir) do
        write("jl.sum", FVSjl.run_keyfile("akplant.key"; variant = FVSjl.SoutheastAlaska(), output = :sum))
        _ak_sum_rows("jl.sum")
    end
    yrs = sort([y for y in keys(live) if y <= 2082])
    @test length(yrs) == 10
    @testset "establishment + small-tree growth (TPA..QMD)" begin
        for y in yrs
            @test get(jl, y, zeros(10))[1:6] == live[y][1:6]
        end
    end
    @testset "volume with the start-of-cycle bark" begin
        for y in yrs
            @test get(jl, y, zeros(10))[7:10] == live[y][7:10]
        end
    end
end

@testset "AK FFE stand vs live FVSak (SIMFIRE 2003 mortality)" begin
    fx = joinpath(@__DIR__, "..", "fixtures", "southeastalaska")
    live = _ak_sum_rows(joinpath(fx, "akffe.live.sum"))
    dir = mktempdir()
    for f in ("akffe.key", "akffe.tre"); cp(joinpath(fx, f), joinpath(dir, f)); end
    jl = cd(dir) do
        write("jl.sum", FVSjl.run_keyfile("akffe.key"; variant = FVSjl.SoutheastAlaska(), output = :sum))
        _ak_sum_rows("jl.sum")
    end
    @test get(jl, 2013, zeros(10))[1] == 139
    for y in (1993, 2003, 2013, 2023)
        @test get(jl, y, zeros(10))[1:6] == live[y][1:6]
    end
end

@testset "AK akt01 vs live FVSak (all stands, every .sum row)" begin
    fx = joinpath(@__DIR__, "..", "fixtures", "southeastalaska")
    function rows_by_stand(path)
        out = Dict{Tuple{Int,Int},Vector{Float64}}(); k = 0
        for l in eachline(path)
            startswith(l, "-999") && (k += 1; continue)
            f = split(l)
            (k > 0 && length(f) > 11 && all(isdigit, f[1]) && length(f[1]) == 4) || continue
            out[(k, parse(Int, f[1]))] = parse.(Float64, f[3:12])
        end
        return out
    end
    live = rows_by_stand(joinpath(fx, "akt01.live.sum"))
    dir = mktempdir()
    for f in ("akt01.key", "akt01.tre"); cp(joinpath(fx, f), joinpath(dir, f)); end
    jl = cd(dir) do
        write("jl.sum", FVSjl.run_keyfile("akt01.key"; variant = FVSjl.SoutheastAlaska(), output = :sum))
        rows_by_stand("jl.sum")
    end
    @test length(live) == 56
    @test count(k -> get(jl, k, nothing) == live[k], collect(keys(live))) == 56
    # every .sum column, incl. removals/after-treatment/accretion/mortality and MAI (evtstv.f TOTREM = Σ the INTEGER
    # removed merch: the shelterwood's 2010 MAI (1537+163)/80 = 21.25 prints "21.2")
    full(path) = Dict((k, l) for (k, l) in enumerate(filter(l -> occursin(r"^(19|20)\d\d ", l), readlines(path))))
    lf = full(joinpath(fx, "akt01.live.sum")); jf = full(joinpath(dir, "jl.sum"))
    @test length(lf) == 56 && all(k -> split(get(jf, k, "")) == split(lf[k]), keys(lf))
end

@testset "AK FIA sample vs live FVSak (12 stands, AUTOES)" begin
    fx = joinpath(@__DIR__, "..", "fixtures", "southeastalaska")
    function rows_by_id(path)
        out = Dict{Tuple{String,Int},Vector{Float64}}(); sid = ""
        for l in eachline(path)
            startswith(l, "-999") && (sid = split(l)[3]; continue)
            f = split(l)
            (sid != "" && length(f) > 11 && all(isdigit, f[1]) && length(f[1]) == 4) || continue
            out[(sid, parse(Int, f[1]))] = parse.(Float64, f[3:12])
        end
        return out
    end
    live = rows_by_id(joinpath(fx, "akfia.live.sum"))
    dir = mktempdir()
    for f in ("akfia.key", "ak_fia.db"); cp(joinpath(fx, f), joinpath(dir, f)); end
    jl = cd(dir) do
        write("jl.sum", FVSjl.run_keyfile("akfia.key"; variant = FVSjl.SoutheastAlaska(), output = :sum))
        rows_by_id("jl.sum")
    end
    exact = ("10705712010497", "10706339010497", "10708179010497", "10708351010497", "1549083042290487",
             "24731081010497", "24739066010497", "644808316126144", "666740939126144", "720755825290487")
    for sid in exact
        ks = [k for k in keys(live) if k[1] == sid]
        @test length(ks) == 7
        @test all(k -> get(jl, k, nothing) == live[k], ks)
    end
    @test count(k -> get(jl, k, nothing) == live[k], collect(keys(live))) >= 81
    # every .sum column of the exact stands — incl. MAI: 720755825290487 has no inventory AGE (0) but trees, so
    # evtstv.f MAIFLG shuts MAI off (0.0) for every row
    function lines_by_id(path)
        out = Dict{Tuple{String,Int},Vector{SubString{String}}}(); sid = ""
        for l in eachline(path)
            startswith(l, "-999") && (sid = split(l)[3]; continue)
            f = split(l)
            (sid != "" && length(f) > 20 && all(isdigit, f[1]) && length(f[1]) == 4) || continue
            out[(sid, parse(Int, f[1]))] = f
        end
        return out
    end
    lf = lines_by_id(joinpath(fx, "akfia.live.sum")); jf = lines_by_id(joinpath(dir, "jl.sum"))
    for sid in exact
        @test all(k -> get(jf, k, nothing) == lf[k], [k for k in keys(lf) if k[1] == sid])
    end
end
