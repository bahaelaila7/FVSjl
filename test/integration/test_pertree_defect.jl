# test_pertree_defect.jl — per-tree DEFECT input (basdam.f → FVSsn vols.f) vs live Fortran.
#
# A tree's damage codes carry volume defect (basdam.f): agent 25 = percent defect for both cubic
# and board, 26 = cubic only, 27 = board only; the paired value (0-99) is the percent, packed into
# the tree's DEFECT as CF·1e6 + BF·1e4. The volume model then reads ICDF = DEFECT/1e6 and folds it
# into the cubic-defect reduction (max with the MCDEFECT/CFDEFT curve). `pertree_defect` is the base
# stand with damage agent 26 (cubic defect) + severity 30 applied to every undamaged tree (the
# scenario uses a TREEFMT that widens the severity field to two digits). Two checks:
#   1. FIRES — vs the same trees with the 26 codes cleared, merch cubic (col 10) drops in every row.
#   2. FORTRAN — TPA/BA and the defect-reduced merch cubic match live Fortran (±Float32 noise).

using Test, FVSjl

const _PD_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_pd_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && tryparse(Int, first(split(l))) !== nothing]
_pd_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_pdcol(r, c) = parse(Float64, r[c])

@testset "per-tree DEFECT input vs Fortran" begin
    key = joinpath(_PD_DIR, "pertree_defect.key")
    sav = joinpath(_PD_DIR, "pertree_defect.sum.save")
    tre = joinpath(_PD_DIR, "pertree_defect.tre")
    if !isfile(key) || !isfile(sav)
        @test_skip "pertree_defect scenario not available"
    else
        jl = _pd_rows(FVSjl.run_keyfile(key; faithful = true))

        # 1. the per-tree defect FIRES: clear the agent-26 codes (cols 52-53 → "00") for a twin,
        #    merch cubic (col 10) must rise back up everywhere.
        nokey = tempname() * ".key"; notre = replace(nokey, ".key" => ".tre")
        write(nokey, join(readlines(key), "\n") * "\n")
        open(notre, "w") do io
            for l in eachline(tre)
                L = rpad(l, 67)
                SubString(L, 52, 53) == "26" && (L = L[1:51] * "00" * L[54:end])
                println(io, rstrip(L))
            end
        end
        # the twin must reference its own .tre — rewrite STDIDENT-independent: same key, new dir.
        ndir = mktempdir(); cp(nokey, joinpath(ndir, "t.key"); force = true)
        cp(notre, joinpath(ndir, "t.tre"); force = true)
        nd = _pd_rows(FVSjl.run_keyfile(joinpath(ndir, "t.key"); faithful = true))
        @test length(jl) == length(nd)
        @test all(_pdcol(nd[i], 10) - _pdcol(jl[i], 10) > 100 for i in 1:length(jl))

        # 2. matches live Fortran on the structural + defect-reduced merch-cubic columns.
        ft = _pd_base(sav)
        @test length(jl) == length(ft)
        if length(jl) == length(ft)
            for i in 1:length(jl)
                @test _pdcol(jl[i], 3) == _pdcol(ft[i], 3)     # TPA
                @test _pdcol(jl[i], 4) == _pdcol(ft[i], 4)     # BA
                @test _pdcol(jl[i], 10) == _pdcol(ft[i], 10)   # merch cubic (defect target)
            end
        end
    end
end
