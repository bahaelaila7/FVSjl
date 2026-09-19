# classify.jl — mechanical classification of numeric fast-tier mismatches (report.jl TSV) by precision class, so
# KNOWN_RESIDUALS entries name the MEASURED class, not a guess:
#   jl_not_f32   gold is a Float32 value but jl's stored value is not ⇒ jl writes a Float64 quantity where FVS
#                stores REAL*4 (storage-precision gap in the writer or a Float64 computation)
#   gold_text9   gold has ≤9 significant digits and is not Float32 (FVS list-directed text) while jl is full precision
#   ulp1 / ulp2_16 / ulpbig   both Float32-representable: Float32 ULP distance
#   material     |Δ| > 1 unit and ≥ 1% (or presence/rowcount/structural)
# Usage: julia test/harness/tiered/classify.jl <report.tsv>  → per (variant,regime,file,col) class counts
isf32(x::Float64) = Float64(Float32(x)) == x
ulps(a::Float64, b::Float64) = abs(Int64(reinterpret(Int32, Float32(a))) - Int64(reinterpret(Int32, Float32(b))))
sig9(s) = (t = replace(lowercase(split(s, 'e')[1]), r"[-+.]" => ""); t = lstrip(t, '0'); length(rstrip(t, '0')) <= 9)

function cls(g, j)
    a = tryparse(Float64, g); b = tryparse(Float64, j)
    (a === nothing || b === nothing) && return "structural"
    d = abs(a - b); rel = a == 0 ? Inf : d / abs(a)
    (d > 1 && rel >= 0.01) && return "material"
    if isf32(a) && !isf32(b); return "jl_not_f32"; end
    if !isf32(a) && sig9(g) && !isf32(b); return "gold_text9"; end
    if isf32(a) && isf32(b)
        u = ulps(a, b); return u <= 1 ? "ulp1" : u <= 16 ? "ulp2_16" : "ulpbig"
    end
    rel < 1e-6 ? "tiny" : "small"
end

