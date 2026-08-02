# Extract BM crown Weibull coefficients from bm/crown.f DATA → data/bluemountains/crown_coeffs_bm.csv.
const L = readlines("/workspace/ForestVegetationSimulator/bm/crown.f")
const Q = Char(47); const BANG = Char(33); const STAR = Char(42)
function dv(name)
    i0 = findfirst(l -> occursin(Regex("^\\s+DATA\\s+" * name * "/"), l), L)
    body = String[]
    for j in i0:length(L)
        s = strip(L[j]); (startswith(s, "C") || isempty(s)) && continue
        seg = j == i0 ? String(L[j][findfirst(Q, L[j])+1:end]) :
              (m = match(r"^\s{5}\S(.*)$", L[j]); m === nothing ? String(L[j]) : String(m.captures[1]))
        ci = findfirst(BANG, seg); ci !== nothing && (seg = seg[1:ci-1])
        cl = findfirst(Q, seg)
        if cl !== nothing; push!(body, seg[1:cl-1]); break; else push!(body, seg); end
    end
    out = Float64[]
    for t in split(join(body, " "), ",")
        t = strip(t); isempty(t) && continue
        if occursin(STAR, t)
            nv = split(t, STAR)
            cnt = strip(nv[1]) == "MAXSP" ? 18 : parse(Int, strip(nv[1]))   # BM MAXSP=18
            append!(out, fill(parse(Float64, strip(nv[2])), cnt))
        else
            push!(out, parse(Float64, strip(t)))
        end
    end
    out
end
open("data/bluemountains/crown_coeffs_bm.csv", "w") do io
    println(io, "species_index,WEIBA,WEIBB0,WEIBB1,WEIBC0,WEIBC1,C0,C1,CRNMLT,DLOW,DHI")
    cols = [dv(n) for n in ("WEIBA","WEIBB0","WEIBB1","WEIBC0","WEIBC1","C0","C1","CRNMLT","DLOW","DHI")]
    for (k, arr) in enumerate(cols); @assert length(arr) == 18 "$k len $(length(arr))"; end
    for sp in 1:18; println(io, join(vcat(sp, [c[sp] for c in cols]), ",")); end
end
println("wrote crown_coeffs_bm.csv (18); WEIBB0[DF]=", dv("WEIBB0")[3], " C0[DF]=", dv("C0")[3])
