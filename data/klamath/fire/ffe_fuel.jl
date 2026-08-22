# NC (Klamath Mountains) FFE initial surface-fuel loading — nc/fmcba.f DATA FULIVE/FULIVI/FUINIE/FUINII.
# WESTERN structure like CR, but NC is a CALIFORNIA/WESTSIDE variant: the live + initial-dead fuel pools
# are interpolated from the TOP TWO cover-type species (COVCA(1..2), weighted by their share of the top-2
# basal area COVCAWT) — nc/fmcba.f:236-431 — NOT from the single dominant COVTYP the interior western
# variants (CR/IE/EM/CI/TT/UT/BM) use. Tables indexed by NC species 1..12; columns (herb, shrub) for the
# live tables, and the 11 FFE size classes (<.25, .25-1, 1-3, 3-6, 6-12, 12-20, 20-35, 35-50, >50, Litter,
# Duff) for the dead tables. Established (60% cover) vs initiating (10% cover); interpolated by PERCOV.
# Values ported verbatim from nc/fmcba.f (from J. Brown / Ottmar 2000). Reuses `_cr_algslp2` (2-pt ALGSLP).
#
# NC species (JSP): 1 OS  2 SP  3 DF  4 WF  5 MA  6 IC  7 BO  8 TO  9 RF  10 PP  11 OH  12 RW

const _NC_FULIVE = Float32[    # established stands (60% cover)  (herb, shrub)  — nc/fmcba.f:92-103
    0.20 0.20;  #  1 other conifer (use DF)
    0.20 0.10;  #  2 sugar pine (use lodgepole pine NI)
    0.20 0.20;  #  3 Douglas-fir
    0.15 0.10;  #  4 white fir
    0.20 0.20;  #  5 madrone (use DF)
    0.20 0.20;  #  6 incense cedar (use DF)
    0.23 0.22;  #  7 California black oak (Gambel oak — Ottmar)
    0.25 0.25;  #  8 tanoak (use QA — Ottmar 2000b)
    0.15 0.10;  #  9 red fir (use white fir)
    0.20 0.25;  # 10 ponderosa pine
    0.25 0.25;  # 11 other hardwood (tanoak — QA Ottmar)
    0.20 0.20;  # 12 coast redwood (use DF)
]

const _NC_FULIVI = Float32[    # initiating stands (10% cover)  (herb, shrub)  — nc/fmcba.f:108-119
    0.40 2.00;  #  1 other conifer (use DF)
    0.40 1.00;  #  2 sugar pine (use lodgepole pine NI)
    0.40 2.00;  #  3 Douglas-fir
    0.30 2.00;  #  4 white fir
    0.40 2.00;  #  5 madrone (use DF)
    0.40 2.00;  #  6 incense cedar (use DF)
    0.55 0.35;  #  7 California black oak (Gambel — Ottmar)
    0.18 2.00;  #  8 tanoak (use QA — Ottmar 2000b, modified)
    0.30 2.00;  #  9 red fir (use white fir)
    0.25 1.00;  # 10 ponderosa pine
    0.18 2.00;  # 11 other hardwood (tanoak — QA Ottmar, modified)
    0.40 2.00;  # 12 coast redwood
]

const _NC_FUINIE = Float32[    # established (60% cover), 11 size classes  — nc/fmcba.f:125-136
    0.9 0.9 1.6 3.5 3.5 0.0 0.0 0.0 0.0 0.6 10.0;  #  1 other conifer (use DF)
    0.9 0.9 1.2 7.0 8.0 0.0 0.0 0.0 0.0 0.6 15.0;  #  2 sugar pine (use lodgepole pine)
    0.9 0.9 1.6 3.5 3.5 0.0 0.0 0.0 0.0 0.6 10.0;  #  3 Douglas-fir
    0.7 0.7 3.0 7.0 7.0 0.0 0.0 0.0 0.0 0.6 25.0;  #  4 white fir
    0.9 0.9 1.6 3.5 3.5 0.0 0.0 0.0 0.0 0.6 10.0;  #  5 madrone (use DF)
    0.9 0.9 1.6 3.5 3.5 0.0 0.0 0.0 0.0 0.6 10.0;  #  6 incense cedar (use DF)
    0.3 0.7 1.4 0.2 0.1 0.0 0.0 0.0 0.0 3.9  0.0;  #  7 California black oak (Gambel — Ottmar)
    0.2 0.6 2.4 3.6 5.6 0.0 0.0 0.0 0.0 1.4 16.8;  #  8 tanoak (aspen — Ottmar)
    0.7 0.7 3.0 7.0 7.0 0.0 0.0 0.0 0.0 0.6 25.0;  #  9 red fir (use white fir)
    0.9 0.9 1.2 7.0 8.0 0.0 0.0 0.0 0.0 0.6 15.0;  # 10 ponderosa pine (use sugar pine)
    0.2 0.6 2.4 3.6 5.6 0.0 0.0 0.0 0.0 1.4 16.8;  # 11 other hardwood (aspen — Ottmar)
    0.9 0.9 1.6 3.5 3.5 0.0 0.0 0.0 0.0 0.6 10.0;  # 12 coast redwood (use DF)
]

const _NC_FUINII = Float32[    # initiating (10% cover), 11 size classes  — nc/fmcba.f:142-153
    0.5 0.5 1.0 1.4 1.4 0.0 0.0 0.0 0.0 0.3  5.0;  #  1 other conifer (use DF)
    0.6 0.7 0.8 2.8 3.2 0.0 0.0 0.0 0.0 0.3  7.0;  #  2 sugar pine (use lodgepole pine)
    0.5 0.5 1.0 1.4 1.4 0.0 0.0 0.0 0.0 0.3  5.0;  #  3 Douglas-fir
    0.5 0.5 2.0 2.8 2.8 0.0 0.0 0.0 0.0 0.3 12.0;  #  4 white fir
    0.1 0.4 5.0 2.2 2.3 0.0 0.0 0.0 0.0 0.8  5.6;  #  5 madrone (use DF)
    0.5 0.5 1.0 1.4 1.4 0.0 0.0 0.0 0.0 0.3  5.0;  #  6 incense cedar (use DF)
    0.1 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.9  0.0;  #  7 California black oak (Gambel — Ottmar)
    0.1 0.4 5.0 2.2 2.3 0.0 0.0 0.0 0.0 0.8  5.6;  #  8 tanoak (aspen — Ottmar)
    0.5 0.5 2.0 2.8 2.8 0.0 0.0 0.0 0.0 0.3 12.0;  #  9 red fir (use white fir)
    0.6 0.7 0.8 2.8 3.2 0.0 0.0 0.0 0.0 0.3  7.0;  # 10 ponderosa pine (use sugar pine)
    0.1 0.4 5.0 2.2 2.3 0.0 0.0 0.0 0.0 0.8  5.6;  # 11 other hardwood (aspen — Ottmar)
    0.5 0.5 1.0 1.4 1.4 0.0 0.0 0.0 0.0 0.3  5.0;  # 12 coast redwood (use DF)
]

# NC live herb/shrub fuel (nc/fmcba.f:369-378): FLIVE(I) = Σ_{j=1,2} ALGSLP(PERCOV, [10,60],
# [FULIVI(I,COVCA(j))·COVCAWT(j), FULIVE(I,COVCA(j))·COVCAWT(j)]) — the TOP-2 cover-type interpolation.
@inline function nc_live_fuel_loading(covca::NTuple{2,Int}, covcawt::NTuple{2,Float32}, percov::Float32)::NTuple{2,Float32}
    herb = 0f0; shrub = 0f0
    @inbounds for j in 1:2
        c = covca[j]
        (1 <= c <= 12) || continue
        wt = covcawt[j]
        herb  += _cr_algslp2(percov, 10f0, 60f0, _NC_FULIVI[c, 1] * wt, _NC_FULIVE[c, 1] * wt)
        shrub += _cr_algslp2(percov, 10f0, 60f0, _NC_FULIVI[c, 2] * wt, _NC_FULIVE[c, 2] * wt)
    end
    return (herb, shrub)
end

# NC initial dead surface fuel (nc/fmcba.f:421-431): STFUEL(ISZ,2) = Σ_{j=1,2} ALGSLP(PERCOV, [10,60],
# [FUINII(ISZ,COVCA(j))·COVCAWT(j), FUINIE(ISZ,COVCA(j))·COVCAWT(j)]) — the TOP-2 cover-type interpolation.
function nc_dead_fuel_loading(covca::NTuple{2,Int}, covcawt::NTuple{2,Float32}, percov::Float32)::Vector{Float32}
    out = zeros(Float32, 11)
    @inbounds for isz in 1:11
        for j in 1:2
            c = covca[j]
            (1 <= c <= 12) || continue
            wt = covcawt[j]
            out[isz] += _cr_algslp2(percov, 10f0, 60f0, _NC_FUINII[c, isz] * wt, _NC_FUINIE[c, isz] * wt)
        end
    end
    return out
end

# =============================================================================
# nc_cwcalc — NC forest-grown crown width (base/cwidth.f → cwcalc.f, IWHO=0). Needed by FMCBA's
# PERCOV (percent canopy cover) — the generic `crown_width` returns the 0.5 default for every NC species,
# collapsing PERCOV≈0 and mis-selecting the initiating-stand fuel loads. NCMAP maps NC species 1..12 →
# a 5-char CWEQN (FIA code + eqn#). nct01 (forest 505 = Region 5) skips the R6 forest-specific BF section
# ⇒ BF=1 (same as CR); the R6-forest (600s) BF factor is a follow-up. Equation forms reuse `_cr_r6m2`
# (Crookston R6 model 2) and the Donnelly a·D^b power form. Values ported verbatim from nc/cwcalc.f.
const _NC_CWMAP = ("12205", "11705", "20205", "01505", "36102", "08105",
                   "81802", "63102", "02006", "12205", "81802", "21104")

# Donnelly (R6) power form a·D^b with OMIND=1 small-tree scaling and a cap (cwcalc.f a·D^b cases).
@inline function _nc_donnelly(a::Float32, b::Float32, d::Float32, cap::Float32)::Float32
    cw = d >= 1f0 ? a * fpow(d, b) : (a * fpow(1f0, b)) * d
    cw > cap && (cw = cap)
    return cw
end

"""
    nc_cwcalc(sp, d, h, cr, barea, el, hi) -> Float32

NC forest-grown crown width (ft) for NC species `sp` (1..12), DBH `d` (in), height `h` (ft),
crown ratio `cr` (%), stand basal area `barea`, elevation `el` (100s ft), Hopkins index `hi`.
Region-5 (forest 505) ⇒ BF=1. Errors loudly on NC species whose CWEQN is not yet ported.
"""
# NC/Klamath forest 505 is REGION-5, where cwcalc.f branches to R5CRWD (ws/r5crwd.f) — a function of sp/D/H only,
# NOT the R6M2 Crookston models nc_cwcalc uses. R5CRWD is FIA-keyed and shared with WS, so the FVS_TreeList CRWDTH
# reuses ws_r5crwd via this NC-species→WS-species map (by FIA code): OS(299)→WS42, SP(117)→1, DF(202)→2, WF(015)→3,
# MA(361)→38, IC(081)→5, BO(818)→31, TO(631)→34, RF(020)→7, PP(122)→8, OH(998)→43, RW(211)→23. Validated 29/29
# vs FVSnc_clean nct01. (nc_cwcalc's R6M2 is the FFE/PERCOV path — a separate, currently-unported NC FFE concern.)
const _NC_TO_WS_R5 = Int[42, 1, 2, 3, 38, 5, 31, 34, 7, 8, 43, 23]
@inline nc_r5crwd(sp::Int, d::Float32, h::Float32)::Float32 =
    (1 <= sp <= 12) ? ws_r5crwd(_NC_TO_WS_R5[sp], d, h) : 0f0

function nc_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32,
                   el::Float32, hi::Float32)::Float32
    (1 <= sp <= 12) || return 0f0
    eqn = _NC_CWMAP[sp]
    cl  = cr * h * 0.01f0
    ba1 = barea + 1f0
    if     eqn == "12205"; return _cr_r6m2(4.7762f0,0.74126f0,-0.28734f0,0.17137f0,-0.00602f0,-0.00209f0, d,h,cl,ba1,el,13f0,75f0,50f0)  # OS/PP
    elseif eqn == "11705"; return _cr_r6m2(3.5930f0,0.63503f0,-0.22766f0,0.17827f0, 0.04267f0,-0.00290f0, d,h,cl,ba1,el, 5f0,75f0,56f0)  # sugar pine
    elseif eqn == "20205"; return _cr_r6m2(6.0227f0,0.54361f0,-0.20669f0,0.20395f0,-0.00644f0,-0.00378f0, d,h,cl,ba1,el, 1f0,75f0,80f0)  # Douglas-fir
    elseif eqn == "01505"; return _cr_r6m2(5.0312f0,0.53680f0,-0.18957f0,0.16199f0, 0.04385f0,-0.00651f0, d,h,cl,ba1,el, 2f0,75f0,35f0)  # white fir
    elseif eqn == "02006"; return _nc_donnelly(3.1146f0,0.5780f0, d, 65f0)                                                               # red fir
    else
        error("nc_cwcalc: crown-width equation $(eqn) (NC species $(sp)) not yet ported — nct01 exercises " *
              "only SP/DF/WF/RF/PP/OS; the MA/IC/BO/TO/OH/RW equations are a follow-up crown-width chunk.")
    end
end
