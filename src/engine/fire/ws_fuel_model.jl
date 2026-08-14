# =============================================================================
# fire/ws_fuel_model.jl — WS (WestSierra) California-CWHR surface fuel-model selection
#
# Ported from ws/fmcfmd.f (FMCFMD) + ws/cwhr.f (== nc/cwhr.f, verified identical). WS uses the same
# California-Wildlife-Habitat-Relationships classifier as NC/CA: classify the stand into a CWHR SIZE(1..6) ×
# DENSITY(S/P/M/D) structural stage, map (forest-type IFT × stage JSS) → a base fuel model via the WS
# CWHRFMD 12×18 table, blend density sub-models (LDYNFM), add the natural-fuel candidates 10/12/13 (+ the
# post-activity 11), and resolve the final weighted set with `_fmdyn` over WS's 15-model XPTS. WS differs from
# NC in: the CWHR constants (CWXPTS/CCBP/DBHBP), the 12-row × 18-col CWHRFMD matrix (which uses the extended
# fuel models 25 & 26 at IPTR positions 14/15), and the 43-species IFT forest-type classification.
# Root cause of the wst01+SIMFIRE over-kill (#220/#221): WS previously fell through to the SN default selection
# → high-intensity models 6+10 (flame 4.20/scorch 17.59) vs the oracle's 8+11 (flame 1.65/scorch 4.21).
# Reuses `_ca_cwhr` (shared CWHR core), `_fmdyn`, `_nc_findjss` (WS FINDJSS == NC's, verified). Values verbatim
# from ws/fmcfmd.f. The California FMSHRUB shrub-delay (fmcfmd.f:413) is not modelled here (as in NC).
# =============================================================================

# ws/fmcfmd.f DATA XPTS — SMALL/LARGE fuel intercepts per POSITION (1..15). Positions 14/15 = models 25/26.
const _WS_FMD_XPTS = Float32[
    5.0 15.0;  5.0 15.0;  5.0 15.0;  5.0 15.0;  5.0 15.0;  5.0 15.0;  5.0 15.0;
    5.0 15.0;  5.0 15.0; 15.0 30.0; 15.0 30.0; 30.0 60.0; 45.0 100.0;
    5.0 15.0;  5.0 15.0]                        # pos 14 = FMD 25, pos 15 = FMD 26

# ws/fmcfmd.f DATA CWXPTS / CCBP / DBHBP — the CWHR canopy-cover breakpoints (differ from NC).
const _WS_CWXPTS = Float32[26.8 187.5; 42.2 435.8; 70.0 350.0; 180.0 180.0]
const _WS_CCBP   = (25f0, 40f0, 60f0)
const _WS_DBHBP  = (1f0, 6f0, 11f0, 24f0)
const _WS_CWHR_IPTR = [10, 20, 30, 40]

# ws/fmcfmd.f DATA IPTR — POSITION → fuel-model number. Models 25/26 sit at positions 14/15.
const _WS_IPTR = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 25, 26]

# FINDMOD (ws/fmcfmd.f:428): fuel-model NUMBER → its position in IPTR (1..15), else 8. Handles 25/26 → 14/15.
@inline function _ws_findmod(jmod::Integer)::Int
    @inbounds for k in 1:15
        _WS_IPTR[k] == jmod && return k
    end
    return 8
end

# ws/fmcfmd.f DATA CWHRFMD(12,18) — rows = forest type IFT (1..12), cols = structural stage JSS (1..18).
const _WS_CWHRFMD = Int[
    9 2 2 9 9 2 2 2 9 2 2 8 8 2 2 8 8 10;    #  1 pine E
    9 5 5 9 9 26 26 25 9 26 26 8 8 26 26 8 8 10;  #  2 pine W (P)
    8 8 8 8 8 11 11 8 8 8 8 8 8 8 8 8 8 10;  #  3 red fir
    8 8 8 8 8 11 11 11 8 8 8 8 8 8 8 8 8 10; #  4 white fir E
    8 5 5 8 8 11 11 8 8 8 8 8 8 8 8 8 8 10;  #  5 white fir W
    8 5 5 8 8 5 5 8 8 11 11 9 8 11 11 9 8 10;  #  6 Douglas-fir
    8 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 10;    #  7 giant sequoia
    9 9 9 9 9 2 2 2 9 2 2 2 9 2 2 2 9 10;    #  8 Jeffrey pine
    8 5 5 9 9 11 11 11 9 9 9 9 9 9 9 9 9 10; #  9 hardwoods
    8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 10;    # 10 lodgepole pine
    9 5 5 9 9 26 26 25 9 26 26 8 8 26 26 8 8 10;  # 11 mixed pine
    8 9 9 8 8 26 26 11 8 5 5 8 8 5 5 8 8 10] # 12 mixed conifer

# ws_cwhr — reuse the shared CA-CWHR core with WS's constants (base cwhr.f is identical to NC's).
ws_cwhr(s::StandState, cws::Vector{Float32}) =
    _ca_cwhr(s, cws, _WS_DBHBP, _WS_CCBP, _WS_CWXPTS, _WS_CWHR_IPTR)

# ws/fmcfmd.f:272-313 — forest-type row IFT (1..12) from per-species BA% + the site-species site index.
# XMCON = mixed-conifer BA (sp {1,5,10:27,42}); XHARD = hardwood BA (sp {28:41,43}).
@inline function _ws_ift(bapct::Vector{Float32}, sitear_isisp::Float32)::Int
    xmcon = 0f0; xhard = 0f0
    @inbounds for i in 1:43
        if i == 1 || i == 5 || (10 <= i <= 27) || i == 42
            xmcon += bapct[i]
        elseif (28 <= i <= 41) || i == 43
            xhard += bapct[i]
        end
    end
    if bapct[8] >= 80f0                                  # PP
        return sitear_isisp < 55f0 ? 1 : 2
    end
    if (bapct[7] + bapct[3]) >= 80f0                     # RF + WF
        if bapct[7] >= bapct[3]
            return 3                                     # red fir
        else
            return sitear_isisp < 55f0 ? 4 : 5           # white fir E / W
        end
    end
    bapct[2] >= 80f0 && return 6                          # Douglas-fir
    bapct[4] >= 80f0 && return 7                          # giant sequoia
    bapct[6] >= 80f0 && return 8                          # Jeffrey pine
    xhard    >= 80f0 && return 9                          # hardwoods
    bapct[9] >= 80f0 && return 10                         # lodgepole pine
    xmcon    >= 80f0 && return 12                         # mixed conifer
    return bapct[8] >= bapct[3] ? 11 : 12                 # mixed pine / mixed conifer
end

"""
    ws_select_fuel_models(s, mois, sm, lg) -> Vector{(model, weight)}

WS California-CWHR fuel-model selection (ws/fmcfmd.f). `sm`/`lg` are the SMALL/LARGE down-wood loads
(the FMDYN point). Uses ws_r5crwd for the per-tree CRWDTH that drives the CWHR canopy cover.
"""
function ws_select_fuel_models(s::StandState, mois::AbstractMatrix{Float32}, sm::Float32, lg::Float32)
    t = s.trees
    eqwt = zeros(Float32, 15)              # WS ICLSS = 15 positions (models 1-13, 25, 26)

    # per-tree crown widths (CRWDTH) via R5CRWD (fn of sp/D/H only)
    cws = zeros(Float32, t.n)
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        cws[i] = ws_r5crwd(Int(t.species[i]), t.dbh[i], t.height[i])
    end

    # per-species BA% + the site species' site index (SITEAR(ISISP))
    bapct = zeros(Float32, 43); stndba = 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        stndba += t.tpa[i] * t.dbh[i]^2 * 0.0054542f0
    end
    if stndba > 0.001f0
        @inbounds for i in 1:t.n
            t.tpa[i] > 0f0 || continue
            bapct[Int(t.species[i])] += 100f0 * (t.tpa[i] * t.dbh[i]^2 * 0.0054542f0) / stndba
        end
    end
    isisp = Int(s.plot.site_species)
    sitear = (1 <= isisp <= length(s.plot.sp_site_index)) ? s.plot.sp_site_index[isisp] : 0f0

    # CWHR structural stage + forest type
    sz, dn, cwhr_mod, cwhr_wt = ws_cwhr(s, cws)
    jss = _nc_findjss(sz, dn)                             # WS FINDJSS == NC's (verified)
    ift = _ws_ift(bapct, sitear)

    # base model from (IFT, JSS)
    if 1 <= ift <= 12 && 1 <= jss <= 18
        eqwt[_ws_findmod(_WS_CWHRFMD[ift, jss])] = 1f0
    else
        eqwt[8] = 1f0
    end

    # dynamic density blending (LDYNFM; only SS 2-5)
    if sz in ('2', '3', '4', '5')
        imd = zeros(Int, 2)
        for i in 1:2
            if cwhr_wt[i] > 0f0
                dc = cwhr_mod[i] == 10 ? 'S' : cwhr_mod[i] == 20 ? 'P' :
                     cwhr_mod[i] == 30 ? 'M' : cwhr_mod[i] == 40 ? 'D' : ' '
                jc = _nc_findjss(sz, dc)
                imd[i] = (1 <= ift <= 12 && 1 <= jc <= 18) ? _ws_findmod(_WS_CWHRFMD[ift, jc]) : 0
            end
        end
        if imd[1] == 0 && imd[2] == 0
            eqwt[8] = 1f0
        elseif imd[1] != 0 && imd[2] != 0 && imd[1] != imd[2]
            eqwt[imd[1]] = cwhr_wt[1]; eqwt[imd[2]] = cwhr_wt[2]
        end
    end

    # FM10/FM11 sharing (ws/fmcfmd.f:383-404) — the 5-yr post-activity FM11 (AFWT/SLCHNG/LATFUEL) is not yet
    # wired (as in NC); with no recent activity AFWT=0 ⇒ FM10=1. FM11 present via the CWHR table bypasses this.
    if eqwt[11] == 0f0
        eqwt[10] = 1f0
    else
        eqwt[10] = 0f0
    end
    eqwt[12] = 1f0
    eqwt[13] = 1f0

    return _fmdyn(sm, lg, eqwt, _WS_FMD_XPTS; iptr = _WS_IPTR)
end
