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

# ci/dgdriv.f:133 DATA PSIGSQ — per-species prior variance for the DG self-calibration COR shrinkage
# (dgdriv.f:665 PVAR=PSIGSQ(ISPC)). Missing branch left CI on the SN default (0.0898) ⇒ wrong shrinkage ⇒
# GF(sp4) COR fit to 0 vs live 0.0569 ⇒ DG shortfall ⇒ mortality over-kill. First 10 = the N-Rockies conifers (=KT).
const CI_PSIGSQ = Float32[0.0408, 0.0586, 0.1556, 0.0970, 0.0858, 0.1433, 0.0636, 0.0970, 0.0970, 0.0636,
                          0.0586, 0.0586, 0.1433, 0.07, 0.0898, 0.0586, 0.07, 0.0858, 0.07]

# ci/bratio.f — per-species BRATIO (like BM). BARK1/BARK2 in species_coefficients.csv.
#   CASE(3,5,9,10) DF/WH/AF/PP: DIB=BARK1·D^BARK2; BRATIO=DIB/D (D>0 else 0.97); cap ≤0.97.
#   CASE(14) WJ: 0.9002 − 0.3089/TEMD (TEMD∈[1,19]); clamp [0.80,0.99].
#   CASE(17,19) CW/OH: BARK1 + BARK2/TEMD (TEMD≥1); clamp [0.80,0.99].
#   DEFAULT (rest): BRATIO = BARK1 (constant).
@inline function ci_bratio(sd, sp::Int, d::Real)::Float32
    b1 = sd[:bark1][sp]; b2 = sd[:bark2][sp]; dd = Float32(d)
    if sp == 3 || sp == 5 || sp == 9 || sp == 10
        r = dd > 0f0 ? b1 * dd^b2 / dd : 0.97f0
        return r > 0.97f0 ? 0.97f0 : r
    elseif sp == 14
        temd = clamp(dd, 1f0, 19f0)
        return clamp(0.9002f0 - 0.3089f0 * (1f0 / temd), 0.80f0, 0.99f0)
    elseif sp == 17 || sp == 19
        temd = dd < 1f0 ? 1f0 : dd
        return clamp(b1 + b2 * (1f0 / temd), 0.80f0, 0.99f0)
    else
        return b1                                            # DEFAULT: constant BARK1
    end
end
