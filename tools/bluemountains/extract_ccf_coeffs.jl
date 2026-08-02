# Extract BM CCF coefficients (bm/ccfcal.f RD1/RD2/RD3/RDA/RDB, 18 each) → data/bluemountains/ccf_coeffs_bm.csv.
const LINES = readlines("/workspace/ForestVegetationSimulator/bm/ccfcal.f")
function dv(name)
    i0 = findfirst(l -> occursin(Regex("^\\s+DATA\\s+" * name * "/"), l), LINES)
    body = String[]
    for j in i0:length(LINES)
        s = strip(LINES[j])
        (startswith(s, "C") || startswith(s, "!") || isempty(s)) && continue
        seg = if j == i0
            LINES[j][findfirst('/', LINES[j])+1:end]
        else
            m = match(r"^\s{5}\S(.*)$", LINES[j]); m === nothing ? LINES[j] : m.captures[1]
        end
        ci = findfirst('!', seg); ci !== nothing && (seg = seg[1:ci-1])
        cl = findfirst('/', seg)
        if cl !== nothing; push!(body, seg[1:cl-1]); break; else push!(body, seg); end
    end
    out = Float64[]
    for t in split(join(body, " "), ',')
        t = strip(t); isempty(t) && continue
        if occursin('*', t)
            n, v = split(t, '*'); append!(out, fill(parse(Float64, strip(v)), parse(Int, strip(n))))
        else
            push!(out, parse(Float64, t))
        end
    end
    out
end
open("data/bluemountains/ccf_coeffs_bm.csv", "w") do io
    println(io, "species_index,RD1,RD2,RD3,RDA,RDB")
    a = dv("RD1"); b = dv("RD2"); c = dv("RD3"); d = dv("RDA"); e = dv("RDB")
    for arr in (a,b,c,d,e); @assert length(arr) == 18; end
    for sp in 1:18; println(io, join([sp, a[sp], b[sp], c[sp], d[sp], e[sp]], ",")); end
    println(stderr, "RD1[DF]=", a[3], " RD2[DF]=", b[3], " RDA[DF]=", d[3], " RDB[DF]=", e[3])
end
println("wrote ccf_coeffs_bm.csv (18)")
