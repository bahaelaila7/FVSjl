# Extract all BM large-tree + small-tree DG coefficient DATA arrays from bm/dgf.f into CSVs.
# Handles Fortran DATA: `&` continuations, interspersed C-comment lines, and N*value repeats.
# Emits:
#   data/bluemountains/dg_coeffs_bm.csv    — 18 species rows, per-species arrays + DGFOR_1..4
#   data/bluemountains/dg_smcoeffs_bm.csv  — 5 small-tree groups, SM* arrays + SMFOR_1..4 + SMHAB_1..5
#   data/bluemountains/dg_maps_bm.csv      — SMMAPS(18); SMMAPH(92,5) + IBSERV(5,3) as JSON-ish rows
const SRC = "/workspace/ForestVegetationSimulator/bm/dgf.f"
const LINES = readlines(SRC)

# Grab the token stream inside `DATA <name>/ ... /`, flattening continuations, skipping C-comments.
function data_values(name::AbstractString)
    # find the DATA line index
    i0 = findfirst(l -> occursin(Regex("^\\s+DATA\\s+" * replace(name, "(" => "\\(", ")" => "\\)", "*" => "\\*")), l), LINES)
    i0 === nothing && error("DATA $name not found")
    # Accumulate the value text line by line from the DATA line to the closing '/'.
    body = String[]
    for j in i0:length(LINES)
        raw = LINES[j]
        stripped = strip(raw)
        (startswith(stripped, "C") || startswith(stripped, "!") || isempty(stripped)) && continue
        seg = raw
        if j == i0
            k = findfirst('/', seg); seg = seg[k+1:end]  # drop `DATA name/`
        else
            # drop leading continuation marker (col-6 '&' or char)
            m = match(r"^\s{5}\S(.*)$", seg)
            seg = m === nothing ? seg : m.captures[1]
        end
        # trailing comment after '!'
        ci = findfirst('!', seg); ci !== nothing && (seg = seg[1:ci-1])
        # closing slash?
        cl = findfirst('/', seg)
        if cl !== nothing
            push!(body, seg[1:cl-1]); break
        else
            push!(body, seg)
        end
    end
    joined = join(body, " ")
    # tokenize, expand N*value
    out = Float64[]
    for tok in split(joined, ',')
        t = strip(tok); isempty(t) && continue
        if occursin('*', t)
            n, v = split(t, '*'); append!(out, fill(parse(Float64, strip(v)), parse(Int, strip(n))))
        else
            push!(out, parse(Float64, t))
        end
    end
    return out
end

names18 = ["DGLD","DGCR","DGCRSQ","DGBA","DGBAL","DGCCFA","DGSIC","DGSITE","DGDBAL",
           "DGLBA","DGPCCF","OBSERV","DGDS","DGCASP","DGSASP","DGSLOP","DGSLSQ","DGEL","DGEL2"]
cols18 = Dict(n => data_values(n) for n in names18)
for (n,v) in cols18; @assert length(v)==18 "$n has $(length(v))"; end
dgfor = data_values("DGFOR");  @assert length(dgfor)==4*18   # (4,18) col-major
# write per-species CSV
open("data/bluemountains/dg_coeffs_bm.csv","w") do io
    println(io, "species_index," * join(names18, ",") * "," * join(["DGFOR_$k" for k in 1:4], ","))
    for sp in 1:18
        vals = [string(cols18[n][sp]) for n in names18]
        for4 = [string(dgfor[(sp-1)*4 + k]) for k in 1:4]   # col-major: DGFOR(k,sp)
        println(io, join(vcat(string(sp), vals, for4), ","))
    end
end

names5 = ["SMLD","SMCR","SMCRSQ","SMDBAL","SMLBA","SMPCCF","SMDS","SMCASP","SMSASP","SMSLOP","SMEL","SMEL2"]
cols5 = Dict(n => data_values(n) for n in names5)
for (n,v) in cols5; @assert length(v)==5 "$n has $(length(v))"; end
smfor = data_values("SMFOR"); @assert length(smfor)==4*5
smhab = data_values("SMHAB"); @assert length(smhab)==5*5
open("data/bluemountains/dg_smcoeffs_bm.csv","w") do io
    println(io, "group," * join(names5, ",") * "," * join(["SMFOR_$k" for k in 1:4],",") * "," * join(["SMHAB_$k" for k in 1:5],","))
    for g in 1:5
        vals = [string(cols5[n][g]) for n in names5]
        f4 = [string(smfor[(g-1)*4+k]) for k in 1:4]         # SMFOR(k,g) col-major
        h5 = [string(smhab[(g-1)*5+k]) for k in 1:5]         # SMHAB(k,g) col-major (INDXH,INDXS)
        println(io, join(vcat(string(g), vals, f4, h5), ","))
    end
end

smmaps = Int.(data_values("SMMAPS")); @assert length(smmaps)==18
smmaph = Int.(data_values("SMMAPH")); @assert length(smmaph)==92*5 "SMMAPH $(length(smmaph))"
ibserv = Int.(data_values("((IBSERV(I,J),I=1,5),J=1,3)")); @assert length(ibserv)==15
open("data/bluemountains/dg_maps_bm.csv","w") do io
    println(io, "# SMMAPS (18): ", join(smmaps, " "))
    println(io, "# IBSERV (5,3) col-major: ", join(ibserv, " "))
    println(io, "smmaph_row,g1,g2,g3,g4,g5")   # SMMAPH(hab, group): 92 rows × 5 groups (col-major)
    for h in 1:92
        row = [smmaph[(g-1)*92 + h] for g in 1:5]
        println(io, join(vcat(string(h), string.(row)), ","))
    end
end
println("wrote dg_coeffs_bm.csv (18), dg_smcoeffs_bm.csv (5), dg_maps_bm.csv (SMMAPS/SMMAPH/IBSERV)")
println("verify DGLD[3](DF)=", cols18["DGLD"][3], " (want 0.57990); SMMAPS=", smmaps)
