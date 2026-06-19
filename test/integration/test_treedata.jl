# C1 integration test — the .tre parser must byte-match Oracle A on real data.

using Test
using FVSjl
using FVSjl: parse_tree_format, parse_tree_record, DEFAULT_TREE_FORMAT, TreeRecord

include(joinpath(@__DIR__, "..", "oracle", "oracle.jl"))

const TRE = joinpath(Oracle.FVSSN_TESTS, "snt01.tre")

# Convert a FVSjl TreeRecord to the 25-tuple shape FVSjulia.parse_tree_record uses.
as_tuple(r::TreeRecord) = (
    r.plot, r.id, r.tpa, r.history, r.species_code, r.dbh, r.diam_growth,
    r.height, r.top_height, r.ht_growth, r.crown_pct,
    r.damage[1], r.damage[2], r.damage[3], r.damage[4], r.damage[5], r.damage[6],
    r.mort_code, r.cut_code,
    r.pest_vars[1], r.pest_vars[2], r.pest_vars[3], r.pest_vars[4], r.pest_vars[5],
    r.birth_age,
)

@testset "Fortran FORMAT parser" begin
    f = parse_tree_format(DEFAULT_TREE_FORMAT)
    @test length(f) == 25
    # I4 then T1 then I7: id field overlaps the plot field (cols 1-7)
    @test f[1] == FVSjl.FormatField(1, 4, :int, 0)
    @test f[2] == FVSjl.FormatField(1, 7, :int, 0)
    @test f[5].kind == :string          # A3 species code
    @test f[6].kind == :float           # F4.1 DBH, 1 implied decimal
    @test f[6].decimals == 1
end

@testset "tree records match Oracle A (snt01.tre)" begin
    @test isfile(TRE)
    fields = parse_tree_format(DEFAULT_TREE_FORMAT)
    mine = [as_tuple(parse_tree_record(fields, l)) for l in eachline(TRE) if !isempty(strip(l))]
    theirs = Oracle.oracle_a_tree_tuples(TRE)
    @test length(mine) == length(theirs)
    @test length(mine) > 0
    nbad = 0
    for (i, (a, b)) in enumerate(zip(mine, theirs))
        if a != b
            nbad += 1
            nbad <= 5 && @info "record $i differs" mine=a oracle=b
        end
    end
    @test nbad == 0
end
