# test_compress.jl — COMPRESS keyword (opt 78, act 250, comprs.f + comcup.f).
#
# COMPRESS reduces the tree list to NCLAS representative records by PC-score clustering, then
# merges each class (PROB-weighted). The 1966 IBM-SSP eigensolver is replaced by
# LinearAlgebra.eigen (project direction), so the exact class PARTITION is not bit-identical to
# Fortran — but the algorithm is faithful and the per-cycle aggregate is conserved. Checks:
#   1. KEYWORD — recognized + scheduled (act 250, params target/PN1/date);
#   2. RECORD COUNT — compress! reduces the live list to exactly NCLAS records;
#   3. TPA CONSERVED — total trees/acre is preserved exactly by the merge;
#   4. FORTRAN AGGREGATE — at the compression cycle, the compressed stand's `.sum` row (TPA / BA /
#      SDI / CCF / TopHt / QMD / cubic volume) is bit-identical to live Fortran (the merge is
#      correct; only the later-cycle trajectory diverges with the substituted eigensolver).

using Test, FVSjl

const _CP_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_cp_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]
_cp_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                                y !== nothing && 1900 < y < 2100)]

@testset "COMPRESS keyword recognition + scheduling" begin
    mkrec(fields, vals, present) =
        FVSjl.KeywordRecord("COMPRESS", "", fields, vals, present, 12, FVSjl.KW_OK, 0)
    s = FVSjl.StandState(FVSjl.Southern())
    FVSjl.kw_compress!(s, mkrec(["1990", "200", "60", fill("", 9)...],
                                Float32[1990, 200, 60, zeros(Float32, 9)...],
                                [true, true, true, falses(9)...]))
    a = last(s.control.schedule)
    @test a.icflag == Int32(250) && a.year == Int32(1990)
    @test a.params[1] == 200f0 && a.params[2] == 60f0
    s2 = FVSjl.StandState(FVSjl.Southern())
    FVSjl.kw_compress!(s2, mkrec(fill("", 12), zeros(Float32, 12), falses(12)))
    d = last(s2.control.schedule)
    @test d.year == Int32(1) && d.params[1] == 1500f0 && d.params[2] == 50f0
    @test !("COMPRESS" in FVSjl.KNOWN_NOOP)
end

@testset "COMPRESS algorithm — record clustering + merge" begin
    key = joinpath(_CP_DIR, "fire_early.key")
    if !isfile(key)
        @test_skip "fire_early scenario not available"
    else
        # RECORD COUNT + TPA CONSERVATION at several targets.
        for nclas in (5, 10, 15, 20)
            s, _ = FVSjl.initialize(key)
            FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
            n0 = s.trees.n
            tpa0 = sum(s.trees.tpa[i] for i in 1:n0)
            @test n0 > nclas                                  # compression actually fires
            @test FVSjl.compress!(s, nclas, 0.5)
            @test s.trees.n == nclas                          # exactly NCLAS records
            tpa1 = sum(s.trees.tpa[i] for i in 1:s.trees.n)
            @test isapprox(tpa0, tpa1; atol = 7f-5)           # TPA conserved to 1 Float32 ULP (merge re-sums tpa in a
                                                              # different order): measured Δ=6.1035e-5 = EXACTLY 2^-14 = eps(Float32(589.65)); rtol 1f-4 was ~1000× padded
        end
        # no-op when the target ≥ the record count
        s2, _ = FVSjl.initialize(key); FVSjl.notre!(s2); FVSjl.setup_growth!(s2)
        @test !FVSjl.compress!(s2, s2.trees.n + 5, 0.5)
        @test s2.trees.n > 0
    end

    # FORTRAN AGGREGATE — compress scenario (NCLAS=15 @ cycle 1): the compression-cycle row's
    # stand aggregates are bit-identical (the merge conserves them) vs live Fortran.
    ckey = joinpath(_CP_DIR, "compress.key")
    sav = joinpath(_CP_DIR, "compress.sum.save")
    if !isfile(ckey) || !isfile(sav)
        @test_skip "compress scenario not available"
    else
        jl = _cp_rows(FVSjl.run_keyfile(ckey; faithful = true))
        ft = _cp_base(sav)
        @test !isempty(jl) && length(jl) == length(ft)
        j = jl[1]; f = ft[1]                                  # compression cycle (1990)
        @test j[1] == f[1]                                    # YEAR
        for col in (3, 4, 5, 6, 7, 8, 9, 10)                  # TPA/BA/SDI/CCF/TopHt/QMD/TCuFt/MCuFt
            @test j[col] == f[col]
        end
    end
end
