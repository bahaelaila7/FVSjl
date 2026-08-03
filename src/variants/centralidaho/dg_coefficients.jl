# CI large-tree DDS coefficient arrays (ci/dgf.f), loaded from data/centralidaho/dg_*.csv
# (extracted by tools/centralidaho/extract_dg_coeffs.py). CI DG = standard western Wykoff DDS.
# DEFAULT (main conifer) DDS: CONSPP + DGLD·lnD + DGLBA·ln(BA) + DGDSQ·D² + DGDBAL·PBAL/ln(D+1)
#                            + DGBA·BAL/ln(D+1) + DGPCCF·PCCF. DGDSQ=DGDS (dgf.f:602). No MAPDSQ/MAPCCF.

let
    d = CI_DATADIR
    rows = readlines(joinpath(d, "dg_coeffs_1d.csv"))
    hdr = split(strip(rows[1]), ',')
    ci(n) = findfirst(==(n), hdr)
    col(n) = Float32[parse(Float32, strip(split(strip(rows[1+sp]), ',')[ci(n)])) for sp in 1:19]
    for n in ("DGLD","DGCR","DGCRSQ","DGBAL","DGDBAL","DGBA","DGLBA","DGPCCF","DGEL","DGEL2",
              "DGSLOP","DGSLSQ","DGCASP","DGSASP","DGCCFA","DGDS","OBSERV")
        @eval global const $(Symbol("CI_" * n)) = $(col(n))
    end
    # dg_dgfor.csv: 19 lines (per species), each = 3 location values ⇒ DGFOR[loc, sp]
    readcols(f, nr) = (ls = readlines(joinpath(d, f));
        hcat([Float32[parse(Float32, strip(x)) for x in split(strip(l), ',')] for l in ls]...))
    global const CI_DGFOR  = readcols("dg_dgfor.csv", 3)     # [3, 19]
    global const CI_MAPLOC = round.(Int, readcols("dg_maploc.csv", 6))   # [6, 19]
    global const CI_IBSERV = round.(Int, readcols("dg_ibserv.csv", 6))   # [6, 3]
    # dg_dghab.csv: 12 lines (per hab group), each 19 species ⇒ DGHAB[group, sp]; rows→matrix
    readrows(f) = (ls = readlines(joinpath(d, f));
        permutedims(hcat([Float32[parse(Float32, strip(x)) for x in split(strip(l), ',')] for l in ls]...)))
    global const CI_DGHAB  = readrows("dg_dghab.csv")        # [12, 19]
    global const CI_ICHBCL = round.(Int, readrows("dg_ichbcl.csv"))      # [130, 19]
end

# ci/bratio.f: POWER model DIB = BARK1·D^BARK2, BRATIO = DIB/D (same family as BM). BARK1/BARK2 in
# species_coefficients.csv (:bark1/:bark2). r>1 or ≤0 ⇒ 0.99 fallback; final clamp handled by caller.
@inline function ci_bratio(sd, sp::Int, d::Real)::Float32
    b1 = sd[:bark1][sp]; b2 = sd[:bark2][sp]
    d <= 0f0 && return 0.99f0
    r = b1 * Float32(d)^b2 / Float32(d)
    (r > 1f0 || r <= 0f0) ? 0.99f0 : r
end
