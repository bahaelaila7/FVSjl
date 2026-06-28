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

# Semantic equivalence ignores field-TEXT formatting (the semantic form emits "11", a
# .key may carry "11.0") — same keyword, numeric value and presence is the contract.
sig2(r) = (strip(r.name), round.(r.values; digits=4), r.present)

@testset "semantic stand YAML (format: fvs-stand) unravels to canonical keyword order" begin
    yml = tempname() * ".yaml"
    open(yml, "w") do io
        write(io, """
        format: fvs-stand/v1
        stand:
          stdident: "SN SEMANTIC TEST — em-dash id"
          noautoes: true
          design: { plots: 11, nsc: 1 }
          stdinfo: { forest: 80106, habitat: "231Dd", age: 60, aspect: 315, slope: 30, elev: 7 }
          sitecode: { species: 63, index: 60.0 }
          invyr: 1990
          numcycle: 6
          echosum: true
          treatments:
            - thinsdi: { year: 2010, sdi: 200 }
          treelist:
            format: "(I4,T1,I7,F6.0,I1,A3,F4.1)"
        """)
    end
    recs  = FVSjl.read_keyword_records(yml)
    names = [strip(r.name) for r in recs if !isempty(strip(r.name))]
    # canonical order: STDIDENT(+id) → NOAUTOES → DESIGN → STDINFO → SITECODE → INVYEAR →
    # NUMCYCLE → ECHOSUM → THINSDI → TREEFMT(+fmt) → TREEDATA → PROCESS → STOP
    @test "STDIDENT" in names && names[end] == "STOP"
    @test findfirst(==("NOAUTOES"), names) < findfirst(==("DESIGN"), names)
    @test findfirst(==("DESIGN"),   names) < findfirst(==("THINSDI"), names)
    @test findfirst(==("THINSDI"),  names) < findfirst(==("TREEDATA"), names)
    @test findfirst(==("TREEDATA"), names) < findfirst(==("PROCESS"),  names)
    # the em-dash id rides verbatim (multi-byte UTF-8 must not crash the decoder)
    @test any(r -> occursin("em-dash id", r.raw), recs)
    # the THINSDI card maps year→field1, residual SDI→field2 (FVS-term keys)
    thin = recs[findfirst(r -> strip(r.name) == "THINSDI", recs)]
    @test thin.values[1] ≈ 2010 && thin.values[2] ≈ 200
    # field positions match a hand-written equivalent .key (numeric value + presence)
    key = tempname() * ".key"
    open(key, "w") do io
        write(io, """
        STDIDENT
        SN SEMANTIC TEST
        NOAUTOES
        DESIGN                                        11.0       1.0
        STDINFO        80106   231Dd        60.0     315.0      30.0       7.0
        SITECODE          63      60.
        INVYEAR       1990.0
        NUMCYCLE         6.0
        ECHOSUM
        THINSDI       2010.0     200.0
        TREEFMT
        (I4,T1,I7,F6.0,I1,A3,F4.1)
        TREEDATA
        PROCESS
        STOP
        """)
    end
    kr = FVSjl.read_keyfile_records(key)
    # compare the keyword cards (bare-token names) — drop the free-text id/FORMAT lines,
    # which are emitted verbatim and differ by content (the semantic id has an em-dash).
    cards(rs) = [sig2(r) for r in rs if occursin(r"^[A-Z][A-Z0-9]*$", strip(r.name))]
    @test cards(recs) == cards(kr)
    rm.([yml, key]; force=true)
end

@testset "semantic multi-stand + raw_keywords escape hatch" begin
    yml = tempname() * ".yaml"
    open(yml, "w") do io
        write(io, """
        format: fvs-stand/v1
        stands:
          - { stdident: "S1", invyr: 1990, numcycle: 6, treelist: { format: "(I4)" } }
          - stdident: "S2"
            invyr: 1990
            numcycle: 6
            treatments: [ { thinbba: { year: 2010, ba: 80 } } ]
            raw_keywords: [ { FMIN: {} }, { END: {} } ]
            treelist: { format: "(I4)" }
        """)
    end
    recs  = FVSjl.read_keyword_records(yml)
    names = [strip(r.name) for r in recs if !isempty(strip(r.name))]
    @test count(==("STDIDENT"), names) == 2     # two stands
    @test count(==("PROCESS"),  names) == 2     # one PROCESS each
    @test count(==("STOP"),     names) == 1     # single terminating STOP, last
    @test names[end] == "STOP"
    @test "FMIN" in names && "THINBBA" in names # escape-hatch keyword rode along
    # FMIN (raw) emits before the modeled THINBBA in the same stand (documented order)
    @test findfirst(==("FMIN"), names) < findlast(==("THINBBA"), names)
    rm(yml; force=true)
end
