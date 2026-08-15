# =============================================================================
# so_fuel_model.jl — SO (SouthCentralOregon) FFE dynamic fuel-model selection (so/fmcfmd.f).
#
# so/fmcfmd.f is REGION-DEPENDENT: KODFOR 500-599 / 701 selects the R5-California SOSPDM dominant-species
# path (an RNG-jittered dominant finder → models 2/5/8/9/10); every other (R6-Oregon) forest uses the
# 8-plant-group path ported here. The reference stand is forest 601 (DESCHUTES, R6), so it rides the
# Oregon branch. The R5-California SOSPDM branch is a documented follow-on (not exercised by the R6
# reference stand, and its RNG dominant-finder is out of scope for the surface-fire smoke).
#
# Oregon branch (so/fmcfmd.f:400-621): IPAG = IPASO(ITYPE) collapses the 92 SORNEC plant associations into
# 8 super-groups (1 dry PP, 2 wet PP, 3 dry MC, 4 wet MC, 5 hemlock, 6 dry LP, 7 wet LP, 8 juniper); each
# builds an EQWT[1..14] over the 14 fuel-model candidates from PERCOV / disturbance-age / ACTCBH / the
# FMSSTAGE structure class (IFMST), plus the always-on natural-fuel models {10,12,13} and the post-activity
# {11,14}. FMDYN then resolves the weighted candidates over the (SMALL,LARGE) down-wood point → FMD.
# =============================================================================

# so/fmcfmd.f DATA IPASO — 92 SORNEC plant associations → 8 Oregon super-groups.
const _SO_IPASO = Int8[
    3, 3, 3, 4, 4, 4, 4, 4, 6, 6,
    7, 6, 6, 6, 6, 6, 6, 6, 6, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 6,
    7, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 5, 1, 2, 2, 2, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 2, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 7, 4, 3, 3, 3,
    4, 3, 3, 3, 3, 3, 3, 3, 3, 7,
    7, 7]

# so/fmcfmd.f XPTS — the (SMALL@LARGE=0, LARGE@SMALL=0) intercepts for each of the 14 fuel models. ICLSS=14
# (SO adds model 14 = activity twin of 12, so it extends the 13-model WC/EC table by one row).
const _SO_FMD_XPTS = Float32[
     5.0  15.0;   5.0  15.0;   5.0  15.0;   5.0  15.0;   5.0  15.0;   5.0  15.0;   5.0  15.0;
     5.0  15.0;   5.0  15.0;  15.0  30.0;  15.0  30.0;  30.0  60.0;  45.0 100.0;  30.0  60.0]

"""
    so_select_fuel_models(s, mois, sm, lg) -> Vector{Tuple{Int,Float32}}

SO dynamic fuel-model selection (so/fmcfmd.f Oregon branch). Builds the 14-candidate EQWT from the plant
super-group IPAG, PERCOV, and (for the LP/juniper groups) the FMSSTAGE structure class, then FMDYN-resolves
over the (SMALL,LARGE) down-wood point. The reference stand (forest 601, ITYPE 49 = CPS111 = dry PP) rides
IPAG=1.
"""
function so_select_fuel_models(s::StandState, mois::AbstractMatrix{Float32}, sm::Float32, lg::Float32)
    eqwt = zeros(Float32, 14)
    percov = s.fire.percov
    itype = Int(s.plot.habitat_input); itype <= 0 && (itype = 49)   # so/habtyp.f DEFAULT CPS111 = ITYPE 49
    ipag = (1 <= itype <= 92) ? Int(_SO_IPASO[itype]) : 1
    alg(x, x1, x2) = _cr_algslp2(x, Float32(x1), Float32(x2), 0f0, 1f0)   # ALGSLP(x,[x1,x2],[0,1],2)

    # FMSSTAGE structure class (IFMST) — only the LP (6,7) and juniper (8) groups use it.
    ifmst = (ipag == 6 || ipag == 7 || ipag == 8) ? Int(structure_class(s).class) : 0
    # ACTCBH (active crown base height) feeds the LP groups (6,7). It is computed at fire time (fmcfir), not
    # persisted — the R6 PP reference stand (IPAG=1) does not use it, so 0 here; the LP/juniper ACTCBH path is
    # a documented follow-on (same deferral class as the R5-California SOSPDM branch).
    actcbh = 0f0
    # DSTLG (years since last disturbance) — undisturbed stands (no harvest/burn) ⇒ large; only the dry-PP
    # WT1(1)>0 (PERCOV<45) branch reads it, which the reference (PERCOV 66.6) does not take.
    dstlg = 999f0

    if ipag == 1                          # DRY PONDEROSA PINE
        w2 = alg(percov, 35, 45); w1 = 1f0 - w2
        w2 > 0f0 && (eqwt[9] += w2)
        if w1 > 0f0
            d2 = alg(dstlg, 5, 7); d1 = 1f0 - d2
            d1 > 0f0 && (eqwt[2] += w1 * d1)
            d2 > 0f0 && (eqwt[6] += w1 * d2)
        end
    elseif ipag == 2                      # WET PONDEROSA PINE
        w2 = alg(percov, 60, 80); w1 = 1f0 - w2
        w1 > 0f0 && (eqwt[5] += w1)
        w2 > 0f0 && (eqwt[9] += w2)
    elseif ipag == 3 || ipag == 4         # DRY / WET MIXED CONIFER (identical 40/60 → 5/8)
        w2 = alg(percov, 40, 60); w1 = 1f0 - w2
        w1 > 0f0 && (eqwt[5] += w1)
        w2 > 0f0 && (eqwt[8] += w2)
    elseif ipag == 5                      # HEMLOCK
        eqwt[8] = 1f0
    elseif ipag == 6                      # DRY LODGEPOLE PINE
        w2 = alg(actcbh, 2, 4); w1 = 1f0 - w2
        if w1 > 0f0
            ifmst <= 1 ? (eqwt[5] += w1) : (eqwt[10] += w1)
        end
        w2 > 0f0 && (eqwt[8] += w2)
    elseif ipag == 7                      # WET LODGEPOLE PINE
        if ifmst <= 1
            w2 = alg(actcbh, 2, 4); w1 = 1f0 - w2
            w1 > 0f0 && (eqwt[5] += w1)
            w2 > 0f0 && (eqwt[8] += w2)
        elseif ifmst > 1 && ifmst < 5
            eqwt[8] = 1f0
        else
            w2 = alg(percov, 40, 60); w1 = 1f0 - w2
            w1 > 0f0 && (eqwt[3] += w1)
            w2 > 0f0 && (eqwt[8] += w2)
        end
    else                                  # ipag == 8 JUNIPER SHRUBLAND
        w2 = alg(percov, 60, 80); w1 = 1f0 - w2
        if w1 > 0f0
            ifmst <= 1 ? (eqwt[1] += w1) : (eqwt[6] += w1)
        end
        w2 > 0f0 && (eqwt[8] += w2)
    end

    # Post-activity models 11/14 (AFWT/LATFUEL, so/fmcfmd.f:633-640) — deferred (no activity fuel ⇒ AFWT=0,
    # LATFUEL false); models 10/12/13 are the always-on natural-fuel candidates.
    eqwt[10] = 1f0
    eqwt[12] = 1f0
    eqwt[13] = 1f0

    return _fmdyn(sm, lg, eqwt, _SO_FMD_XPTS)
end
