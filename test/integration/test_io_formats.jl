# C1b — modern I/O formats must losslessly round-trip the record schema, and the
# legacy/modern front-ends must agree (so the engine is format-agnostic).

using Test
using FVSjl
using FVSjl: read_tree_records, read_trees_csv, write_trees_csv, read_tree_file,
             convert_tre_to_csv, TreeRecord, TREE_CSV_HEADER

include(joinpath(@__DIR__, "..", "oracle", "oracle.jl"))

const TRE = joinpath(Oracle.FVSSN_TESTS, "snt01.tre")

# CSV stores stripped species codes; normalize for comparison.
strip_sp(r::TreeRecord) = TreeRecord(
    r.plot, r.id, r.tpa, r.history, rpad(strip(r.species_code), 8), r.dbh, r.diam_growth,
    r.height, r.top_height, r.ht_growth, r.crown_pct, r.damage, r.mort_code, r.cut_code,
    r.pest_vars, r.birth_age)

@testset "CSV trees: lossless round-trip" begin
    recs = read_tree_file(TRE)
    @test length(recs) == 30
    tmp = tempname() * ".csv"
    write_trees_csv(recs, tmp)
    back = read_trees_csv(tmp)
    @test strip_sp.(recs) == strip_sp.(back)
    # header is the documented schema
    @test split(readlines(tmp)[1], ',') == TREE_CSV_HEADER
    rm(tmp; force=true)
end

@testset "format-agnostic loader: .tre and .csv agree" begin
    tmp = tempname() * ".csv"
    convert_tre_to_csv(TRE, tmp)
    from_tre = read_tree_records(TRE)       # legacy path
    from_csv = read_tree_records(tmp)       # modern path
    @test strip_sp.(from_tre) == strip_sp.(from_csv)
    rm(tmp; force=true)
end

# dispatch-relevant signature of a keyword record (what the engine consumes)
kw_sig(r) = (strip(r.name), round.(r.values; digits=4), r.present,
             [strip(f) for f in r.fields], r.status)

@testset "YAML keywords: lossless round-trip + format-agnostic loader" begin
    for keyname in ("sn.key", "snt01.key")
        keypath = joinpath(Oracle.FVSSN_TESTS, keyname)
        isfile(keypath) || continue
        orig = FVSjl.read_keyfile_records(keypath)
        tmp = tempname() * ".yaml"
        FVSjl.convert_key_to_yaml(keypath, tmp)
        back = FVSjl.read_keyword_records(tmp)              # modern path
        @test length(orig) == length(back)
        @test kw_sig.(orig) == kw_sig.(back)
        # legacy loader agrees with direct lexer
        @test kw_sig.(FVSjl.read_keyword_records(keypath)) == kw_sig.(orig)
        rm(tmp; force=true)
    end
end
