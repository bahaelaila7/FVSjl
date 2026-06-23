# C1b — modern I/O formats must losslessly round-trip the record schema, and the
# legacy/modern front-ends must agree (so the engine is format-agnostic).

using Test
using FVSjl
using FVSjl: read_tree_records, read_trees_csv, write_trees_csv, read_tree_file,
             convert_tre_to_csv, convert_csv_to_tre, write_tree_file, TreeRecord, TREE_CSV_HEADER,
             write_keyfile, convert_yaml_to_key, translate_io

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

@testset "reverse translation: modern → legacy round-trips" begin
    # .tre → write back → .tre reproduces the same records (incl. the I4,T1,I7 plot/id overlap)
    recs = read_tree_file(TRE)
    tre2 = tempname() * ".tre"
    write_tree_file(recs, tre2)
    @test strip_sp.(read_tree_file(tre2)) == strip_sp.(recs)
    # full .tre → .csv → .tre chain
    csv = tempname() * ".csv"; tre3 = tempname() * ".tre"
    convert_tre_to_csv(TRE, csv); convert_csv_to_tre(csv, tre3)
    @test strip_sp.(read_tree_file(tre3)) == strip_sp.(recs)
    rm.([tre2, csv, tre3]; force=true)

    # .key → .yaml → .key reproduces the dispatch-relevant keyword signature
    for keyname in ("sn.key", "snt01.key")
        keypath = joinpath(Oracle.FVSSN_TESTS, keyname)
        isfile(keypath) || continue
        orig = FVSjl.read_keyfile_records(keypath)
        yml = tempname() * ".yaml"; key2 = tempname() * ".key"
        FVSjl.convert_key_to_yaml(keypath, yml)
        convert_yaml_to_key(yml, key2)
        @test kw_sig.(FVSjl.read_keyfile_records(key2)) == kw_sig.(orig)
        rm.([yml, key2]; force=true)
    end

    # a keyfile with an inline TREEFMT FORMAT string (free-form supplemental lines) must
    # survive key→yaml→key verbatim — the raw lines are carried, not mangled into keywords.
    fmtkey = joinpath(@__DIR__, "..", "harness", "scenarios", "fire_early.key")
    if isfile(fmtkey)
        orig = FVSjl.read_keyfile_records(fmtkey)
        yml = tempname() * ".yaml"; key2 = tempname() * ".key"
        FVSjl.convert_key_to_yaml(fmtkey, yml); convert_yaml_to_key(yml, key2)
        @test any(occursin("raw:", l) for l in readlines(yml))           # format line stored as raw
        @test kw_sig.(FVSjl.read_keyfile_records(key2)) == kw_sig.(orig)  # incl. the format-string lines
        rm.([yml, key2]; force=true)
    end
end

@testset "translate_io picks direction by extension" begin
    yml = tempname() * ".yaml"; key2 = tempname() * ".key"
    csv = tempname() * ".csv"; tre2 = tempname() * ".tre"
    keypath = joinpath(Oracle.FVSSN_TESTS, "snt01.key")
    if isfile(keypath)
        translate_io(keypath, yml); translate_io(yml, key2)
        @test kw_sig.(FVSjl.read_keyfile_records(key2)) == kw_sig.(FVSjl.read_keyfile_records(keypath))
    end
    translate_io(TRE, csv); translate_io(csv, tre2)
    @test strip_sp.(read_tree_file(tre2)) == strip_sp.(read_tree_file(TRE))
    @test_throws ErrorException translate_io("a.key", "b.csv")   # unsupported cross-kind
    rm.([yml, key2, csv, tre2]; force=true)
end
