# =============================================================================
# fire/ca_fuel_model.jl — CA (CentralCalifornia) California-CWHR surface fuel-model selection
#
# Ported from ca/fmcfmd.f (FMCFMD) + ca/cwhr.f. CA uses the same California-CWHR classifier as NC/WS. CA's CWHR
# constants (CWXPTS/CCBP/DBHBP) + IPTR (1..13) + FINDJSS/FINDMOD + the FM10/11/12/13 natural-fuel tail are all
# IDENTICAL to NC's (verified) — so this reuses the NC framework (_ca_cwhr with NC consts, _nc_findjss/_nc_findmod,
# _NC_FMD_XPTS). CA-specific: the CWHRFMD 11×18 matrix (ca/fmcfmd.f; FMD_R5 == FMD_R6, one table) and the
# 50-species IFT forest-type classification. Crown widths for the CWHR canopy cover come from ca_cwcalc (the R6
# Crookston CAMAP; cat01 forest 610). Root cause of the CA fire under-kill: without this, CA fell to the default
# SN fuel-model selection → weak fire (oracle's near-total 530→2 needs CA's high-intensity models). Verbatim.
# =============================================================================

# ca/fmcfmd.f DATA FMD_R5/FMD_R6 — 11 forest-type rows × 18 structural-stage cols (the two tables are identical).
const _CA_CWHRFMD = Int[
    5 6 6 6 6 2 2 9 9 2 2 2 9 2 2 9 9 10;      #  1 pine
    5 5 5 8 8 11 11 8 8 8 8 8 8 8 8 8 8 10;    #  2 red fir
    5 5 5 8 8 11 11 11 8 8 8 8 8 8 8 8 8 10;   #  3 white fir E
    5 5 5 8 8 11 11 8 8 8 8 8 8 8 8 8 8 10;    #  4 white fir W
    5 5 5 6 6 6 6 8 8 11 11 9 8 11 11 9 8 10;  #  5 Douglas-fir
    5 2 2 6 6 2 2 2 9 2 2 2 9 2 2 2 9 10;      #  6 Jeffrey pine
    5 5 5 6 6 11 11 11 9 9 9 9 9 9 9 9 9 10;   #  7 hardwoods
    8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 10;      #  8 lodgepole pine
    5 5 5 6 6 6 6 6 9 9 9 8 8 8 8 8 8 10;      #  9 mixed pine
    5 5 5 6 6 6 6 6 8 6 6 8 8 6 6 8 8 10;      # 10 mixed conifer
    5 5 5 6 6 6 6 6 8 6 6 8 8 6 6 8 8 10]      # 11 other softwood

# ca_cwhr — reuse the shared CA-CWHR core with NC's constants (CA's CWXPTS/CCBP/DBHBP == NC's, verified).
ca_cwhr(s::StandState, cws::Vector{Float32}) =
    _ca_cwhr(s, cws, _NC_DBHBP, _NC_CCBP, _NC_CWXPTS, _NC_CWHR_IPTR)

# ca/fmcfmd.f:236-258 — forest-type row IFT (1..11) from per-species BA% + the site-species site index.
# XMCON = mixed-conifer BA (sp {1:3,8:11,13:14,16:17,19:25}); XHARD = hardwood BA (sp {26:49}).
@inline function _ca_ift(bapct::Vector{Float32}, sitear_isisp::Float32)::Int
    xmcon = 0f0; xhard = 0f0
    @inbounds for i in 1:50
        if (1 <= i <= 3) || (8 <= i <= 11) || (13 <= i <= 14) || (16 <= i <= 17) || (19 <= i <= 25)
            xmcon += bapct[i]
        elseif 26 <= i <= 49
            xhard += bapct[i]
        end
    end
    if bapct[18] >= 80f0                                   # pine (PP)
        return 1
    end
    if (bapct[4] + bapct[5] + bapct[6]) >= 80f0            # WF + RF + SF
        if (bapct[5] + bapct[6]) >= bapct[4]
            return 2                                        # red fir
        else
            return sitear_isisp < 55f0 ? 3 : 4              # white fir E / W
        end
    end
    bapct[7]  >= 80f0 && return 5                           # Douglas-fir
    bapct[15] >= 80f0 && return 6                           # Jeffrey pine
    xhard     >= 80f0 && return 7                           # hardwoods
    bapct[12] >= 80f0 && return 8                           # lodgepole pine
    xmcon     >= 80f0 && return 10                          # mixed conifer
    return bapct[18] > bapct[4] ? 9 : 11                    # mixed pine / other softwood
end

"""
    ca_select_fuel_models(s, mois, sm, lg) -> Vector{(model, weight)}

CA California-CWHR fuel-model selection (ca/fmcfmd.f). `sm`/`lg` are the SMALL/LARGE down-wood loads (the FMDYN
point). Uses ca_cwcalc (R6 Crookston CAMAP) for the per-tree CRWDTH that drives the CWHR canopy cover.
"""
function ca_select_fuel_models(s::StandState, mois::AbstractMatrix{Float32}, sm::Float32, lg::Float32)
    t = s.trees
    eqwt = zeros(Float32, 13)                              # CA ICLSS = 13 models (ca/fmcfmd.f)

    # per-tree crown widths (CRWDTH) — same load-time BAREA=1 cyc1 clamp as NC/WS (Crookston uses BAREA)
    ba = (s.control.cycle <= Int32(1)) ? 1f0 : s.plot.basal_area
    el = s.plot.elevation; hi = _cr_hopkins(s.plot.latitude, s.plot.longitude, s.plot.elevation)
    cws = zeros(Float32, t.n)
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        cws[i] = ca_cwcalc(Int(t.species[i]), t.dbh[i], t.height[i], Float32(t.crown_pct[i]), ba, el, hi)
    end

    # per-species BA% + the site species' site index (SITEAR(ISISP))
    bapct = zeros(Float32, 50); stndba = 0f0
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

    sz, dn, cwhr_mod, cwhr_wt = ca_cwhr(s, cws)
    jss = _nc_findjss(sz, dn)                              # CA FINDJSS == NC's (verified)
    ift = _ca_ift(bapct, sitear)

    # base model from (IFT, JSS)
    if 1 <= ift <= 11 && 1 <= jss <= 18
        eqwt[_nc_findmod(_CA_CWHRFMD[ift, jss])] = 1f0
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
                imd[i] = (1 <= ift <= 11 && 1 <= jc <= 18) ? _nc_findmod(_CA_CWHRFMD[ift, jc]) : 0
            end
        end
        if imd[1] == 0 && imd[2] == 0
            eqwt[8] = 1f0
        elseif imd[1] != 0 && imd[2] != 0 && imd[1] != imd[2]
            eqwt[imd[1]] = cwhr_wt[1]; eqwt[imd[2]] = cwhr_wt[2]
        end
    end

    # FM10/FM11 sharing + FM12/FM13 natural fuels (ca/fmcfmd.f:362-404, identical to NC/WS)
    if eqwt[11] == 0f0
        eqwt[10] = 1f0
    else
        eqwt[10] = 0f0
    end
    eqwt[12] = 1f0
    eqwt[13] = 1f0

    return _fmdyn(sm, lg, eqwt, _NC_FMD_XPTS)
end
