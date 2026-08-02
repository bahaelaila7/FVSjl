# Extract the UT SBB height coefficients from ut/htgf.f → data/utah/htgf_cof_ut.csv (COF[K,1:9], K=1..55)
# + data/utah/htgf_zbias_ut.csv (AZBIAS/BZBIAS/MPCRSI per species, 24).
# COF(9,55) is EQUIVALENCEd from COF1..COF11, each COF_n(9,5) col-major ⇒ the flat DATA stream of the 11
# blocks in EQUIVALENCE order (COF1→K1-5, COF2→K6-10, … COF11→K51-55) lays out as K-major, 9 L-values each:
# COF(L,K) = flat[(K-1)*9 + L]. Verifiable: re-run vs source.
const SRC = "/workspace/ForestVegetationSimulator/ut/htgf.f"
lines = readlines(SRC)

function extract(name)
    istart = findfirst(l -> occursin(Regex("^\\s*DATA\\s+$(name)\\s*/"), l), lines)
    istart === nothing && error("DATA $name not found")
    buf = IOBuffer(); started = false; i = istart
    while i <= length(lines)
        seg = lines[i]
        if !started
            si = findfirst('/', seg); si === nothing && (i += 1; continue)
            seg = seg[si+1:end]; started = true
        end
        ei = findfirst('/', seg)
        if ei !== nothing; print(buf, seg[1:ei-1]); break; else; print(buf, seg, " "); end
        i += 1
    end
    toks = String[]
    for t0 in split(String(take!(buf)), r"[,\s]+"; keepempty=false)
        t = lstrip(t0, '&'); isempty(t) && continue
        if occursin('*', t); n, v = split(t, '*'); append!(toks, fill(String(v), parse(Int, n)))
        else; push!(toks, String(t)); end
    end
    return toks
end

# COF blocks in EQUIVALENCE order (COF1 at K=1, COF2 at K=6, … per the EQUIVALENCE list).
cofnames = ["COF1","COF2","COF3","COF4","COF5","COF6","COF7","COF8","COF9","COF10","COF11"]
flat = String[]
for nm in cofnames
    t = extract(nm); length(t) == 45 || error("$nm: $(length(t)) tokens, want 45"); append!(flat, t)
end
length(flat) == 495 || error("COF flat: $(length(flat))")
open("data/utah/htgf_cof_ut.csv", "w") do io
    println(io, "K,", join(["COF$l" for l in 1:9], ","))
    for k in 1:55
        row = [flat[(k-1)*9 + l] for l in 1:9]
        println(io, k, ",", join(row, ","))
    end
end

az = extract("AZBIAS"); bz = extract("BZBIAS"); mp = extract("MPCRSI")
open("data/utah/htgf_zbias_ut.csv", "w") do io
    println(io, "species_index,AZBIAS,BZBIAS,MPCRSI")
    for sp in 1:24; println(io, sp, ",", az[sp], ",", bz[sp], ",", mp[sp]); end
end
println("wrote htgf_cof_ut.csv (55 K-rows) + htgf_zbias_ut.csv (24 sp)")
