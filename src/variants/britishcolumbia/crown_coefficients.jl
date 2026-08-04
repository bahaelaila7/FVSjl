# =============================================================================
# crown_coefficients.jl (britishcolumbia) — BC V3 crown-RATIO model coefficients (canada/bc/crown.f).
#
# CHUNK 5b. Two BEC-zone tables (ICH, IDF), each 15 COLUMNS mapped to species via SEQ (CRKONST%Spp).
# CRCONS resolves per-species site-adjusted CRCON + slopes ONCE per stand (bc_crcons!). The per-cycle
# crown change (crown_ratio_update!) then uses CRNMD: a metric logistic 1/(1+exp(·)) of a linear
# predictor in (H/D, H, D², BAL) with the site-adjusted constant + CRLNCCF·log(RELDEN) density term.
# Column order in the DATA (crown.f:178-272): SEQ, [obs], CON, HTDBH, HT, DBH2, BAL, LNCCF, EL, EL2,
# SLP, SLP2, SASP, CASP, SEE. Stored here column-ordered; bc_crcons! inverts SEQ → per-species (faithful).
# =============================================================================

struct BCCrZone
    zone::String
    seq::Vector{Int}          # column → species (0 = unused column)
    con::Vector{Float32}      # intercept
    htdbh::Vector{Float32}    # H/D slope
    ht::Vector{Float32}       # H slope
    dbh2::Vector{Float32}     # D² slope
    bal::Vector{Float32}      # BAL slope
    lnccf::Vector{Float32}    # ln(CCF/RELDEN) density term
    el::Vector{Float32}       # elevation (site-adj)
    el2::Vector{Float32}      # elevation²
    slp::Vector{Float32}      # slope
    slp2::Vector{Float32}     # slope²
    sasp::Vector{Float32}     # slope·sin(aspect)
    casp::Vector{Float32}     # slope·cos(aspect)
    see::Vector{Float32}      # std err est (CRSD = sqrt(SEE))
end

const _BC_CR_ICH = BCCrZone("ICH",
    [13,12,9,6,11,3,5,2,7,1,10,8,14,15,0],
    Float32[-5.79412,-1.90243,-7.12668,-9.07770,-3.25686,-5.61114,-7.49226,-5.29520,-4.64241,-6.50196,-6.45206,-2.85285,-5.61114,-3.25686,0.0],
    Float32[0.51146,0.0,0.93059,1.66517,0.56474,0.70786,1.13328,1.15971,1.05587,1.22750,0.0,2.75133,0.70786,0.56474,0.0],
    Float32[0.0,0.02921,0.01364,-0.00873,0.01872,0.04186,-0.01291,0.00550,0.03164,0.01112,0.0,0.0,0.04186,0.01872,0.0],
    Float32[0.0,-0.00043,0.00031,0.00017,0.0,-0.00017,0.00030,0.0,-0.00021,0.0,-0.00043,0.0,-0.00017,0.0,0.0],
    Float32[0.0,0.01200,0.01943,0.00904,0.0,0.01275,0.00939,0.0,0.0,0.0,0.0,0.0,0.01275,0.0,0.0],
    Float32[1.05450,0.41115,0.76411,1.35115,0.68275,0.80422,1.09699,0.88346,0.66560,0.95538,0.0,0.0,0.80422,0.68275,0.0],
    Float32[0.0,0.0,0.18717,-0.04259,-0.21167,-0.00873,0.0,0.0,0.12483,0.0,2.06051,0.0,-0.00873,-0.21167,0.0],
    Float32[0.0,0.0,-0.00911,0.0,0.01092,0.0,-0.00482,0.0,-0.00819,0.0,-0.13668,0.0,0.0,0.01092,0.0],
    Float32[0.0,0.0,0.0,0.0,0.0,0.0,1.37891,0.0,-0.49553,0.0,0.0,0.0,0.0,0.0,0.0],
    Float32[0.0,0.0,0.0,0.0,0.0,0.0,-1.65551,-0.52292,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
    Float32[0.0,0.0,0.0,0.19629,0.0,-0.09833,0.23237,0.0,0.14207,0.0,0.0,0.23259,-0.09833,0.0,0.0],
    Float32[-1.22540,0.0,0.42349,0.35682,0.11238,0.28971,0.40582,-0.26042,0.11301,0.0,0.0,0.0,0.28971,0.11238,0.0],
    Float32[0.14,0.12,0.18,0.17,0.13,0.06,0.17,0.11,0.11,0.16,0.09,0.19,0.06,0.13,0.0],
)

const _BC_CR_IDF = BCCrZone("IDF",
    [13,9,6,11,3,2,7,10,8,12,1,14,15,0,0],
    Float32[-1.11300,-0.18008,-5.14375,-4.27986,-6.13692,-2.91609,-4.55575,-1.81808,-0.88873,-5.79412,-6.50196,-6.13692,-4.27986,0.0,0.0],
    Float32[-0.60066,0.0,0.96161,0.45446,1.76502,0.63874,0.72156,1.50771,1.08889,0.51146,1.22750,1.76502,0.45446,0.0,0.0],
    Float32[0.09484,-0.00938,-0.02965,0.04155,0.02994,0.02136,0.06559,0.01341,0.04229,0.0,0.01112,0.02994,0.04155,0.0,0.0],
    Float32[-0.00168,0.0,0.00082,-0.00068,-0.00010,-0.00039,-0.00065,0.0,0.0,0.0,0.0,-0.00010,-0.00068,0.0,0.0],
    Float32[0.02445,0.0,0.02964,-0.01616,-0.00590,0.0,0.00617,0.0,0.02837,0.0,0.0,-0.00590,-0.01616,0.0,0.0],
    Float32[0.21742,0.0,0.0,0.93227,0.83551,0.68627,0.49905,0.49928,0.0,1.05450,0.95538,0.83551,0.93227,0.0,0.0],
    Float32[0.0,0.03312,0.64904,-0.12703,0.0,-0.16511,0.26356,-0.33412,-0.30236,0.0,0.0,0.0,-0.12703,0.0,0.0],
    Float32[0.0,-0.00420,-0.04465,0.00557,0.0,0.00636,-0.01682,0.01516,0.01057,0.0,0.0,0.0,0.00557,0.0,0.0],
    Float32[0.0,0.0,0.83931,-1.86690,0.0,0.0,-0.44586,-1.77894,-2.37636,0.0,0.0,0.0,-1.86690,0.0,0.0],
    Float32[0.0,0.0,0.0,2.34112,0.11122,0.0,0.0,2.98975,1.63206,0.0,0.0,0.11122,2.34112,0.0,0.0],
    Float32[0.0,0.0,-0.63232,0.0,-0.09148,0.0,0.0,0.0,0.0,0.0,0.0,-0.09148,0.0,0.0,0.0],
    Float32[0.0,0.0,0.76681,0.50400,-0.10358,0.69217,0.20520,0.0,1.06805,-1.22540,0.0,-0.10358,0.50400,0.0,0.0],
    Float32[0.12,0.19,0.13,0.15,0.16,0.12,0.15,0.14,0.16,0.14,0.16,0.16,0.15,0.0,0.0],
)

"""Crown-model BEC zone (crown.f:657-673): SBSdw2/SBPS→IDF, SBS→ICH, else match ICH/IDF, default ICH."""
function bc_cr_zone(zone::AbstractString, series::AbstractString)
    (occursin("SBSdw2", series) || occursin("SBPS", zone)) && return _BC_CR_IDF
    occursin("SBS", zone) && return _BC_CR_ICH
    occursin("IDF", zone) && return _BC_CR_IDF
    _BC_CR_ICH
end

"""
    bc_crcons!(s) -> (crcon, crhtdbh, crht, crdbh2, crbal, crlnccf, crsd) each Vector{Float32}(15)

CRCONS (crown.f:634-720): resolve the site-adjusted per-species crown constants ONCE per stand.
CRCON = CON + EL·(ELEV·100·FTtoM) + EL2·(ELEV·100·FTtoM)² + SLP·SLOPE + SLP2·SLOPE² +
        SASP·sin(ASP)·SLOPE + CASP·cos(ASP)·SLOPE.  Slopes copied through the SEQ→species map.
"""
function bc_crcons!(s::StandState)
    p = s.plot
    zone, series = bc_stand_zone(s)
    z = bc_cr_zone(zone, series)
    elm = Float32(p.elevation) * 100f0 * BC_FTtoM            # ELEV(100ft) → m
    slp = Float32(p.slope); asp = Float32(p.aspect)
    sa = sin(asp) * slp; ca = cos(asp) * slp
    nsp = 15
    crcon = zeros(Float32, nsp); crhtdbh = zeros(Float32, nsp); crht = zeros(Float32, nsp)
    crdbh2 = zeros(Float32, nsp); crbal = zeros(Float32, nsp); crlnccf = zeros(Float32, nsp)
    crsd = zeros(Float32, nsp)
    @inbounds for i in 1:15
        j = z.seq[i]; (j < 1 || j > nsp) && continue         # unused column
        crcon[j] = z.con[i] + z.el[i]*elm + z.el2[i]*elm*elm +
                   z.slp[i]*slp + z.slp2[i]*slp*slp + z.sasp[i]*sa + z.casp[i]*ca
        crhtdbh[j] = z.htdbh[i]; crht[j] = z.ht[i]; crdbh2[j] = z.dbh2[i]
        crbal[j] = z.bal[i]; crlnccf[j] = z.lnccf[i]; crsd[j] = sqrt(z.see[i])
    end
    return (crcon, crhtdbh, crht, crdbh2, crbal, crlnccf, crsd)
end

"""
    bc_crnmd(con, htdbh, ht, dbh2, bal, d_in, h_ft, bal_ftac) -> crown ratio ∈ (0,1)

CRNMD (crown.f:730, per-cycle path, D≥3in ⇒ no <2cm dub): metric logistic
YCR = con + htdbh·(Hm/Dm) + ht·Hm + dbh2·Dm² + bal·BALm; CR = 1/(1+exp(clamp(YCR,-9.21,9.21))).
"""
@inline function bc_crnmd(con::Float32, htdbh::Float32, ht::Float32, dbh2::Float32, bal::Float32,
                          d_in::Float32, h_ft::Float32, bal_ftac::Float32)
    dm = d_in * BC_INtoCM; hm = h_ft * BC_FTtoM; balm = bal_ftac * BC_FT2pACRtoM2pHA
    ycr = con + htdbh*(hm/dm) + ht*hm + dbh2*dm*dm + bal*balm
    ycr = min(9.21f0, max(-9.21f0, ycr))
    return 1f0 / (1f0 + exp(ycr))
end
