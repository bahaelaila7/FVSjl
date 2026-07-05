# C1 integration test — the keyword lexer must match Oracle A's KEYRDR.

using Test
using FVSjl
using FVSjl: KeywordReader, read_keyword!, KW_OK, KW_EOF, KW_STOP

include(joinpath(@__DIR__, "..", "oracle", "oracle.jl"))

"Read all keyword records from a file with FVSjl (stop at EOF/STOP)."
function fvsjl_keywords(path)
    recs = Tuple{String,Vector{Float32},Vector{Bool}}[]
    open(path) do io
        r = KeywordReader(io)
        while true
            rec = read_keyword!(r)
            rec.status == KW_EOF && break
            rec.status == KW_STOP && break
            push!(recs, (strip(rec.name), copy(rec.values), copy(rec.present)))
        end
    end
    return recs
end

@testset "keyword lexer vs Oracle A" begin
    for keyname in ("sn.key", "snt01.key")
        keypath = joinpath(Oracle.FVSSN_TESTS, keyname)
        isfile(keypath) || continue
        mine   = fvsjl_keywords(keypath)
        theirs = Oracle.oracle_a_keywords(keypath)
        @test length(mine) == length(theirs)
        nbad = 0
        for (i, (a, b)) in enumerate(zip(mine, theirs))
            namematch = a[1] == b[1]
            valmatch  = a[2] == b[2]        # keyword field VALUES parse bit-identically (both read the same .key literals)
            presmatch = a[3] == b[3]
            if !(namematch && valmatch && presmatch)
                nbad += 1
                nbad <= 5 && @info "$keyname record $i differs" mine=a oracle=b
            end
        end
        @test nbad == 0
    end
end
