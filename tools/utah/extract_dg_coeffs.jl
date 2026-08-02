# Extract the UT diameter-growth coefficient DATA arrays from ut/dgf.f into data/utah/dg_coeffs_ut.csv.
# Parses each `DATA <name>/ ... /` block, expands Fortran `N*val` repeats, and lays the flat value
# stream into per-species rows. 1D arrays (24 values) → one col; 2D arrays (K,24) are stored column-major
# in Fortran ⇒ the flat stream is species-major (K consecutive values per species) ⇒ reshape (K, 24)'.
# Verifiable: re-run vs the source any time. Bit-exact validation happens in diameter_growth.jl (per-tree WK2).

const SRC = "/workspace/ForestVegetationSimulator/ut/dgf.f"
const NSP = 24

# (name, K) — K = values per species (1 for 1D arrays; the leading dim for 2D (K,MAXSP)).
const ARRAYS = [
    ("DGLD",1), ("DGCR",1), ("DGCRSQ",1), ("DGPCCF",1), ("DGBA",1), ("DGBAL",1), ("DGDBAL",1),
    ("ISMAP",1), ("IBSERV",5), ("OBSERV",1), ("DGCCFA",1), ("MAPLOC",6), ("DGFOR",5), ("MAPDSQ",6),
    ("DGDS",4), ("DGCASP",1), ("DGSASP",1), ("DGSLOP",1), ("DGSLSQ",1), ("DGEL",1), ("DGEL2",1),
    ("IDGSIM",7), ("DGSIC",5),
]

lines = readlines(SRC)

"Return the flat token list for `DATA <name>/ ... /` (multi-line, & continuations), repeats expanded."
function extract(name)
    # find the DATA line for this exact name
    istart = findfirst(l -> occursin(Regex("^\\s*DATA\\s+$(name)\\b|^\\s*DATA\\s*\\(\\s*$(name)\\b"), l), lines)
    istart === nothing && error("DATA $name not found")
    # accumulate from the first '/' after the name through the closing '/'
    buf = IOBuffer()
    started = false
    i = istart
    while i <= length(lines)
        l = lines[i]
        # strip a leading comment continuation marker column-6 '&' and trailing comments
        seg = l
        if !started
            si = findfirst('/', seg)
            si === nothing && (i += 1; continue)
            seg = seg[si+1:end]; started = true
        end
        ei = findfirst('/', seg)
        if ei !== nothing
            print(buf, seg[1:ei-1]); break
        else
            print(buf, seg, " ")
        end
        i += 1
    end
    raw = String(take!(buf))
    # tokenize on commas/spaces, strip the '&' continuation markers
    toks = String[]
    for t0 in split(raw, r"[,\s]+"; keepempty=false)
        t = lstrip(t0, '&')                      # column-6 continuation marker, may abut the value (&6*1)
        isempty(t) && continue
        if occursin('*', t)                     # Fortran N*val repeat
            n, v = split(t, '*')
            append!(toks, fill(String(v), parse(Int, n)))
        else
            push!(toks, String(t))
        end
    end
    return toks
end

# Build a wide CSV: one row per species, columns name_k.
open("data/utah/dg_coeffs_ut.csv", "w") do io
    headers = String["species_index"]
    cols = Dict{String,Vector{String}}()
    for (name, K) in ARRAYS
        toks = extract(name)
        length(toks) == K * NSP || error("$name: got $(length(toks)) tokens, want $(K*NSP)")
        for k in 1:K
            h = K == 1 ? name : "$(name)_$k"
            push!(headers, h)
            cols[h] = [toks[(sp-1)*K + k] for sp in 1:NSP]   # species-major flat stream
        end
    end
    println(io, join(headers, ","))
    for sp in 1:NSP
        row = String[string(sp)]
        for h in headers[2:end]; push!(row, cols[h][sp]); end
        println(io, join(row, ","))
    end
end
println("wrote data/utah/dg_coeffs_ut.csv ($NSP rows)")
