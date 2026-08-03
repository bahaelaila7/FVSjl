# Extract the REAL CI species crosswalk from ci/spctrn.f ASPT table → species_translation.csv.
# SPCTRN selects the target by VARIANT: CASE('CI') → SPCOUT = ASPT(I,7)(1:3) (ci/spctrn.f:1213).
# ASPT row cols: 1 alpha_in, 2 FIA, 3 PLANTS, 4 AK, 5 BM, 6 CA, 7 CI, ... The J=1,10 DATA blocks
# hold cols 1-10, so col7 (CI target) is captured there.
const SRC = "/workspace/ForestVegetationSimulator/bin/FVSci_buildDir/spctrn.f"
function extract()
    rows = Tuple{String,String,String,String}[]
    in_head = false
    for raw in readlines(SRC)
        if occursin(r"DATA\s*\(\(ASPT\(I,J\),J=1,10\)", raw)
            in_head = true; continue
        elseif occursin(r"DATA\s*\(\(ASPT\(I,J\),J=11,21\)", raw)
            in_head = false; continue
        end
        in_head || continue
        l = strip(raw)
        occursin(r"^&\s*'", l) || continue
        ci = findfirst('!', l); ci !== nothing && (l = l[1:ci-1])
        l = strip(replace(l, r"^&" => ""))
        toks = [strip(replace(t, "'" => "")) for t in split(l, ',')]
        length(toks) < 7 && continue
        alpha = toks[1]; fia = toks[2]; plants = toks[3]; cit = toks[7]   # ASPT col7 = CI target
        (isempty(alpha) && isempty(fia) && isempty(plants)) && continue
        push!(rows, (alpha, fia, plants, cit))
    end
    return rows
end
rows = extract()
open("data/centralidaho/species_translation.csv", "w") do io
    println(io, "code_alpha,code_fia,code_plants,target_ci,target_ci2,target_ci3,target_ci4")
    for (a, f, p, c) in rows
        println(io, join([a, f, p, c, c, c, c], ","))
    end
end
println("wrote data/centralidaho/species_translation.csv ($(length(rows)) rows)")
for (a, f, p, c) in rows
    (a in ("WP", "GF", "WH", "RC", "AF", "LP", "DF", "WL", "PP", "ES", "WJ", "MC")) && println("  $a/$f/$p → $c")
end
