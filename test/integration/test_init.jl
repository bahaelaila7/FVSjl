# C2 — stand initialization: keyword dispatch + tree loading (INTREE) for snt01.

using Test
using FVSjl
using FVSjl: initialize, parse_tree_format, parse_tree_record

include(joinpath(@__DIR__, "..", "oracle", "oracle.jl"))

const KEY = joinpath(Oracle.FVSSN_TESTS, "snt01.key")

@testset "initialize! snt01 stand 1: keyword state" begin
    s, reason = initialize(KEY)
    @test reason == :process                      # stopped at the first PROCESS
    @test s.trees.n == 27                          # 30 records − 2 dead − 1 non-stockable
    @test s.control.ncycle == 10                   # NUMCYCLE 10
    @test s.control.cycle_year[1] == 1990           # INVYEAR 1990
    @test strip(s.plot.stand_id) == "S248112"      # STDIDENT
    @test occursin("UNTHINNED", s.control.title)
    @test s.plot.site_index == 60f0                # SITECODE 63 60.
    @test s.plot.points_inv == 11                  # DESIGN field 4
    @test occursin("T24", s.control.tree_format)   # custom TREEFMT applied
end

@testset "tree attributes: crown conversion + known species" begin
    s, _ = initialize(KEY)
    t = s.trees
    # crown class code → percent: ICR<10 → ICR*10-5
    @test t.crown_pct[1] == 35                     # input 4 → 35
    @test t.crown_pct[2] == 55                     # input 6 → 55
    # SM is a direct alpha match → sugar maple (index 22)
    @test t.species[2] == 22
    @test t.dbh[1] == 11.5f0
    @test t.height[1] == 73.0f0
end

@testset "cycle-0 density stats match Oracle A (.sum stand 1)" begin
    s, _ = initialize(KEY)
    @test s.trees.n == 27                      # 30 records − 2 dead (ITH 6,8) − 1 IMC1=8
    FVSjl.notre!(s)
    g = s.plot.gross_space                      # reciprocal stockable multiplier
    tpa = FVSjl.stand_tpa(s) / g
    ba  = FVSjl.stand_ba(s)  / g
    sdi = FVSjl.stand_sdi(s) / g
    qmd = FVSjl.stand_qmd(s)
    # targets from snt01.sum / cruise summary (per total acre)
    @test isapprox(tpa, 536.05; rtol=0.002)
    @test isapprox(ba,   77.39; rtol=0.002)
    @test round(Int, sdi) == 160
    @test isapprox(qmd, 5.14;  atol=0.05)
    @test round(Int, FVSjl.stand_top_height(s)) == 63   # AVHT40 top height
    @test s.plot.pi == 11f0                      # PI overwritten with IPTINV
end

@testset "species loading matches Oracle A (custom format)" begin
    s, _ = initialize(KEY)
    # re-parse raw species codes with the same custom format the loader used
    fields = parse_tree_format(s.control.tree_format)
    codes = String[]
    for line in eachline(joinpath(Oracle.FVSSN_TESTS, "snt01.tre"))
        isempty(strip(line)) && continue
        occursin("-999", line) && break
        r = parse_tree_record(fields, line)
        r === nothing && break
        r.mort_code == 8 && continue            # non-stockable plot record (skipped by loader)
        (6 <= r.history <= 9) && continue        # dead tree (skipped by loader)
        push!(codes, strip(r.species_code))
    end
    @test length(codes) == s.trees.n
    oracle_idx = Oracle.oracle_a_resolve(codes)
    mine = [Int(s.trees.species[i]) for i in 1:s.trees.n]
    @test mine == oracle_idx
end
