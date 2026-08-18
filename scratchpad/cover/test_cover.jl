# test_cover.jl — COVER beachhead: CVCW crown-area dump-replay bit-exact vs FVSem_g16.
#
# Standalone-runnable (no FVSjl module needed) so it can be validated in isolation;
# when wired, register it in test/runtests.jl and drop the `include` shim below.
#
#   JULIA_DEPOT_PATH=/workspace/.julia_depot julia test_cover.jl
#
# Golden = scratchpad/cover/cvcw_g16_dump.txt, produced by an instrumented
# FVSem_g16 (single-.o swap of cvcw.f; instrumented COVER report proven
# byte-identical to clean apart from the run timestamp).

using Test

const _HERE = @__DIR__
isdefined(Main, :CoverState) || include(joinpath(_HERE, "cover.jl"))

parse_e22(s) = Float32(parse(Float64, replace(strip(s), "E" => "e")))

# group per-invocation (each CVCWFINAL closes a CVCW call)
function load_dump(path)
    invs = Vector{Tuple{Int,Int,Vector{Float32},Vector{Float32},Float32}}()
    cw = Float32[]; pr = Float32[]
    for ln in eachline(path)
        f = split(ln)
        if startswith(ln, "CVCWTREE")            # icyc tree isp TRECW P CRAREA
            push!(cw, parse_e22(f[5])); push!(pr, parse_e22(f[6]))
        elseif startswith(ln, "CVCWFINAL")        # icyc itrn CRAREA
            push!(invs, (parse(Int, f[2]), parse(Int, f[3]), copy(cw), copy(pr),
                         parse_e22(f[4])))
            empty!(cw); empty!(pr)
        end
    end
    return invs
end

@testset "COVER CVCW crown-area dump-replay (EM, bit-exact)" begin
    invs = load_dump(joinpath(_HERE, "cvcw_g16_dump.txt"))
    @test !isempty(invs)
    cv = CoverState()
    for (icyc, itrn, cw, pr, orc) in invs
        @test length(cw) == itrn
        jl = cover_cvcw!(cv, cw, pr, itrn)
        @test reinterpret(UInt32, jl) == reinterpret(UInt32, orc)   # hex-exact
    end
end
