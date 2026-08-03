# =============================================================================
# dg_coefficients.jl (bluemountains) — large-tree + small-tree DG coefficient arrays
# (bm/dgf.f DATA + bm/dgdriv.f PSIGSQ + bm/bratio.f BARK1/BARK2). Loaded at include time
# from the CSVs written by tools/bluemountains/extract_dg_coeffs.jl (verified vs source).
# =============================================================================

let
    path = joinpath(BM_DATADIR, "dg_coeffs_bm.csv")
    hdr = split(strip(readlines(path)[1]), ',')
    rows = [split(strip(l), ',') for l in readlines(path)[2:end]]
    ci(name) = findfirst(==(name), hdr)
    col1d(name) = Float32[parse(Float32, rows[sp][ci(name)]) for sp in 1:18]
    col2d(name, L) = Float32[parse(Float32, rows[sp][ci("$(name)_$k")]) for k in 1:L, sp in 1:18]
    for n in ("DGLD","DGCR","DGCRSQ","DGBA","DGBAL","DGCCFA","DGSIC","DGSITE","DGDBAL",
              "DGLBA","DGPCCF","OBSERV","DGDS","DGCASP","DGSASP","DGSLOP","DGSLSQ","DGEL","DGEL2")
        @eval global const $(Symbol("BM_"*n)) = $(col1d(n))
    end
    global const BM_DGFOR = col2d("DGFOR", 4)                 # [4, 18] (location class × sp)
end

let
    path = joinpath(BM_DATADIR, "dg_smcoeffs_bm.csv")
    hdr = split(strip(readlines(path)[1]), ',')
    rows = [split(strip(l), ',') for l in readlines(path)[2:end]]
    ci(name) = findfirst(==(name), hdr)
    col1d(name) = Float32[parse(Float32, rows[g][ci(name)]) for g in 1:5]
    col2d(name, L) = Float32[parse(Float32, rows[g][ci("$(name)_$k")]) for k in 1:L, g in 1:5]
    for n in ("SMLD","SMCR","SMCRSQ","SMDBAL","SMLBA","SMPCCF","SMDS","SMCASP","SMSASP","SMSLOP","SMEL","SMEL2")
        @eval global const $(Symbol("BM_"*n)) = $(col1d(n))
    end
    global const BM_SMFOR = col2d("SMFOR", 4)                 # [4, 5] (location × small-group)
    global const BM_SMHAB = col2d("SMHAB", 5)                 # [5, 5] (INDXH × INDXS)
end

let
    lines = readlines(joinpath(BM_DATADIR, "dg_maps_bm.csv"))
    global const BM_SMMAPS = [parse(Int, x) for x in split(strip(split(lines[1], ": ")[2]))]      # (18)
    global const BM_IBSERV_FLAT = [parse(Int, x) for x in split(strip(split(lines[2], ": ")[2]))] # (5,3) col-major
    # SMMAPH: 92 hab rows × 5 groups
    maph = zeros(Int, 92, 5)
    for l in lines[4:end]
        f = split(strip(l), ','); isempty(f[1]) && continue
        h = parse(Int, f[1]); for g in 1:5; maph[h, g] = parse(Int, f[g+1]); end
    end
    global const BM_SMMAPH = maph                             # [92, 5]
end
BM_IBSERV(isic, j) = BM_IBSERV_FLAT[(j-1)*5 + isic]          # IBSERV(ISIC, J), col-major

# bm/dgdriv.f DATA PSIGSQ (per-species DG serial-correlation prior variance).
const BM_PSIGSQ = Float32[0.0408, 0.0586, 0.1556, 0.0970, 0.0858, 0.0858, 0.0636, 0.0970, 0.0970,
                          0.0636, 0.0408, 0.0408, 0.0898, 0.0898, 0.1433, 0.0898, 0.0636, 0.0898]

# bm/bratio.f — BM bark ratio DIB/D. POWER model (DIB=BARK1·D^BARK2) for most species; per-group forms.
# BARK1/BARK2 live in species_coefficients.csv (:bark1/:bark2). Final clamp [0.80, 0.99].
function bm_bratio(sd, sp::Integer, d::Real)::Float32
    b1 = sd[:bark1][sp]; b2 = sd[:bark2][sp]
    br = if sp in (1,2,3,4,5,7,8,9,10,17)
        # BRATIO = DIB/D where DIB = BARK1·D^BARK2 (bm/bratio.f:70-71). r>1 or ≤0 ⇒ 0.999.
        d > 0f0 ? (r = b1 * Float32(d)^b2 / Float32(d); (r > 1f0 || r <= 0f0) ? 0.999f0 : r) : 0.999f0
    elseif sp in (13,14,18)
        d > 0f0 ? b1 * Float32(d)^b2 / Float32(d) : 0.99f0
    elseif sp == 6
        temd = clamp(Float32(d), 1f0, 19f0); 0.9002f0 - 0.3089f0 * (1f0 / temd)
    elseif sp == 11 || sp == 15
        b1
    elseif sp == 12
        temd = max(Float32(d), 1f0); b1 + b2 * (1f0 / temd)
    elseif sp == 16
        d > 0f0 ? (b1 + b2 * Float32(d)) / Float32(d) : 0.99f0
    else
        0.99f0
    end
    return clamp(br, 0.80f0, 0.99f0)
end
