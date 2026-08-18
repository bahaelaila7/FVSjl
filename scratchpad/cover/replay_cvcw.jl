# Dump-replay of CVCW (COVER crown-area) transfer function.
# CVCW resets CRAREA=0 each call, loops I=1,ITRN: CRAREA += TRECW^2*P,
# then CRAREA *= 0.785398. Replays each invocation (flush on CVCWFINAL),
# bit-exact in Float32 as Fortran does.

function parse_e22(s::AbstractString)::Float32
    return Float32(parse(Float64, replace(strip(s), "E" => "e")))
end

const CONST = 0.785398f0
buf = Tuple{Float32,Float32}[]
ninv = 0; nok = 0
for ln in eachline(joinpath(@__DIR__, "fort.71"))
    if startswith(ln, "CVCWTREE")
        f = split(ln)
        push!(buf, (parse_e22(f[5]), parse_e22(f[6])))
    elseif startswith(ln, "CVCWFINAL")
        f = split(ln)
        icyc = parse(Int, f[2]); itrn = parse(Int, f[3])
        orc = parse_e22(f[4])
        acc = 0.0f0
        for (trecw, p) in buf
            acc = acc + trecw*trecw*p
        end
        jl = acc * CONST
        ok = reinterpret(UInt32, jl) == reinterpret(UInt32, orc)
        global ninv += 1; global nok += ok ? 1 : 0
        if !ok || ninv <= 6
            println("inv=$ninv icyc=$icyc ntree=$(length(buf))/$itrn jl=$jl g16=$orc ",
                    "hexjl=", string(reinterpret(UInt32, jl), base=16),
                    " hexg16=", string(reinterpret(UInt32, orc), base=16),
                    ok ? "  BIT-EXACT" : "  MISMATCH")
        end
        empty!(buf)
    end
end
println("invocations=$ninv  bit_exact=$nok  ", nok==ninv ? "ALL BIT-EXACT" : "FAILURES=$(ninv-nok)")
