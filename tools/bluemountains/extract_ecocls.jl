# Extract the BM ECOCLS eco-class table (bm/ecocls.f DATA, 283 rows) → data/bluemountains/ecocls.csv.
# Each DATA row: 'PA','SCIEN', SDIMX,'SPC', SITE, NUMBR, IFLAG, FVSSEQ  ! comment. IFLAG=1 marks the
# eco-class's SITE SPECIES (the FVS reference). Columns: pa (plant-assoc code), spc (species alpha),
# fvsseq (FVS species #), sdimx (SDI max), site (site index), numbr (#species in eco), iflag.
const SRC = "/workspace/ForestVegetationSimulator/bm/ecocls.f"

rows = Vector{NTuple{7,String}}()
for raw in readlines(SRC)
    # a data row line starts (after &) with a quoted PA code and has the 8-field pattern
    l = strip(raw)
    startswith(l, "&'") || continue
    body = l[2:end]                       # drop &
    ci = findfirst('!', body); ci !== nothing && (body = body[1:ci-1])   # drop trailing ! comment
    body = rstrip(strip(body), ['/', ' ', ','])                          # drop closing / and trailing punct
    # tokenize on commas, honoring that quoted strings contain no commas here
    toks = strip.(split(body, ','))
    length(toks) == 8 || continue
    unq(t) = strip(replace(t, "'" => ""))
    pa = unq(toks[1]); spc = unq(toks[4]); fvsseq = strip(toks[8])
    sdimx = strip(toks[3]); site = strip(toks[5]); numbr = strip(toks[6]); iflag = strip(toks[7])
    push!(rows, (pa, spc, fvsseq, sdimx, site, numbr, iflag))
end

open("data/bluemountains/ecocls.csv", "w") do io
    println(io, "pa,spc,fvsseq,sdimx,site,numbr,iflag")
    for r in rows; println(io, join(r, ",")); end
end
println("wrote data/bluemountains/ecocls.csv ($(length(rows)) rows; want 283)")
