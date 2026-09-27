# C1 integration test — the keyword lexer must match Oracle A's KEYRDR.

using Test
using FVSjl
using FVSjl: KeywordReader, read_keyword!, KW_OK, KW_EOF, KW_STOP

include(joinpath(@__DIR__, "..", "oracle", "oracle.jl"))

"Read all keyword records from a file with FVSjl (stop at EOF/STOP)."
function fvsjl_keywords(path)
    recs = Tuple{String,Vector{Float32},Vector{Bool},Vector{String}}[]
    open(path) do io
        r = KeywordReader(io)
        while true
            rec = read_keyword!(r)
            rec.status == KW_EOF && break
            rec.status == KW_STOP && break
            push!(recs, (strip(rec.name), copy(rec.values), copy(rec.present), copy(rec.fields)))
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
        nblank = 0
        for (i, (a, b)) in enumerate(zip(mine, theirs))
            namematch = a[1] == b[1]
            # keyword field VALUES parse bit-identically (both read the same .key literals) — except a field with an
            # EMBEDDED blank: keyrdr.f:144 READ(KARD(I),'(G10.0)') is an internal read, BLANK='NULL', so the blanks
            # are ignored ("  7 34" = 734). Oracle A (FVSjulia) reads such a field as 0; there the reference is
            # keyrdr.f itself: the value must equal the field with its blanks removed.
            valmatch  = all(eachindex(a[2])) do k
                a[2][k] == b[2][k] && return true
                f = k <= length(a[4]) ? strip(a[4][k]) : ""
                occursin(' ', f) || return false
                nblank += 1
                a[2][k] == parse(Float32, replace(f, " " => ""))
            end
            presmatch = a[3] == b[3]
            if !(namematch && valmatch && presmatch)
                nbad += 1
                nbad <= 5 && @info "$keyname record $i differs" mine=a oracle=b
            end
        end
        @test nbad == 0
        nblank > 0 && @info "$keyname: $nblank embedded-blank fields decoded per keyrdr.f:144 (Oracle A reads 0)"
    end
end
