# Extract the REAL BM species crosswalk from bm/spctrn.f ASPT table → species_translation.csv.
# SPCTRN selects the target by VARIANT: CASE('BM') → SPCOUT = ASPT(I,5). ASPT row = (col1 alpha,
# col2 FIA, col3 PLANTS, col4 AK, col5 BM, col6 CA, ...). The J=1,10 DATA blocks hold cols 1-10, so
# col5 (BM target) is there. The prior CSV was an EM placeholder (col10) → WF/GF mis-mapped to AF.
const SRC = "/workspace/ForestVegetationSimulator/bin/FVSbm_buildDir/spctrn.f"
function extract()
rows = Tuple{String,String,String,String}[]
in_head = false                                # true only inside a `DATA ((ASPT(I,J),J=1,10)...` block
for raw in readlines(SRC)
    if occursin(r"DATA\s*\(\(ASPT\(I,J\),J=1,10\)", raw)
        in_head = true; continue
    elseif occursin(r"DATA\s*\(\(ASPT\(I,J\),J=11,21\)", raw)
        in_head = false; continue
    end
    in_head || continue
    l = strip(raw)
    occursin(r"^&\s*'", l) || continue
    ci = findfirst('!', l); ci !== nothing && (l = l[1:ci-1])       # drop trailing ! comment
    l = strip(replace(l, r"^&" => ""))
    toks = [strip(replace(t, "'" => "")) for t in split(l, ',')]
    length(toks) < 5 && continue
    alpha = toks[1]; fia = toks[2]; plants = toks[3]; bm = toks[5]   # ASPT col5 = BM target
    (isempty(alpha) && isempty(fia) && isempty(plants)) && continue
    push!(rows, (alpha, fia, plants, bm))
end
return rows
end
rows = extract()
open("data/bluemountains/species_translation.csv", "w") do io
    println(io, "code_alpha,code_fia,code_plants,target_bm,target_bm2,target_bm3,target_bm4")
    for (a,f,p,bm) in rows
        println(io, join([a, f, p, bm, bm, bm, bm], ","))
    end
end
println("wrote species_translation.csv ($(length(rows)) rows)")
# sanity: WF/GF/WH should map to GF (grand fir) / WH per BM col5
for (a,f,p,bm) in rows
    (a in ("WF","GF","WH","SF","AF","LP","DF","WL","PP","ES") ) && println("  $a/$f/$p → $bm")
end
