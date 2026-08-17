# test_cover.jl — COVER beachhead: CVCW crown-area dump-replay bit-exact vs FVSem_g16.
#
# COVER (covr/vcovr) is a REPORT-ONLY understory/canopy-cover extension (no tree-record
# write, no RNG), active in the shipped western binaries. This beachhead validates the
# CVCW crown-area first computation: CRAREA = (Σ TRECW(i)²·PROB(i))·0.785398.
#
# Golden = cover_cvcw_g16_dump.txt, produced by an instrumented FVSem_g16 (single-.o
# swap of cvcw.f; instrumented COVER report proven byte-identical to clean apart from
# the run timestamp), EM emt01, 2 cycles → 3 CVCW invocations after tripling.

using Test

const _COVER_HERE = @__DIR__

_cover_parse_e22(s) = Float32(parse(Float64, replace(strip(s), "E" => "e")))

# group per-invocation (each CVCWFINAL closes a CVCW call)
function _cover_load_dump(path)
    invs = Vector{Tuple{Int,Int,Vector{Float32},Vector{Float32},Float32}}()
    cw = Float32[]; pr = Float32[]
    for ln in eachline(path)
        f = split(ln)
        if startswith(ln, "CVCWTREE")            # icyc tree isp TRECW P CRAREA
            push!(cw, _cover_parse_e22(f[5])); push!(pr, _cover_parse_e22(f[6]))
        elseif startswith(ln, "CVCWFINAL")        # icyc itrn CRAREA
            push!(invs, (parse(Int, f[2]), parse(Int, f[3]), copy(cw), copy(pr),
                         _cover_parse_e22(f[4])))
            empty!(cw); empty!(pr)
        end
    end
    return invs
end

@testset "COVER — CVCW crown-area dump-replay (EM, bit-exact vs FVSem_g16)" begin
    invs = _cover_load_dump(joinpath(_COVER_HERE, "cover_cvcw_g16_dump.txt"))
    @test !isempty(invs)
    cv = FVSjl.CoverState()
    for (icyc, itrn, cw, pr, orc) in invs
        @test length(cw) == itrn
        jl = FVSjl.cover_cvcw!(cv, cw, pr, itrn)
        @test reinterpret(UInt32, jl) == reinterpret(UInt32, orc)   # hex-exact
    end
end
