# =============================================================================
# fire/wc_fuel_model.jl — WC (WestCascades) FIRE-VPN cover-type-metagroup fuel-model selection.
#
# Ported from wc/fmcfmd.f (VARACD 'WC' branch). This is the Pacific-Northwest "FIRE-VPN" structure, NOT the
# California-CWHR classifier used by NC/WS/CA. It pools the 39 species into 6 cover metagroups
# (SF/DF/MH/RA/LP/WO), takes the top 2 by basal area (RDPSRT), computes QMD80 (the QMD of the lower-80%-BA
# trees, borrowed from WS-FFE CWHR — wc/fmcfmd.f:236) when SF or DF dominates, then applies a set of
# ALGSLP weighting rules keyed by QMD80, PERCOV, and the habitat's forb/grass/shrub (PNFGS→MAPFGS) and
# moist (PNWET→MAPDRY) flags to build the candidate fuel-model weights. Models 10/12/13 are always natural-
# fuel candidates; model 11 is the 5-yr post-activity fuel-jump (deferred — no activity ⇒ AFWT=0). Resolved
# by `_fmdyn` over WC's 13-model XPTS. Values verbatim from wc/fmcfmd.f + wc/fmcba.f (MAPFGS/MAPDRY).
# =============================================================================

# wc/fmcfmd.f DATA XPTS — SMALL/LARGE fuel intercepts per model (1..13). Differs from the default (model 10
# = (15,30) here vs (10,30)); ICLSS=13 (no model 14).
const _WC_FMD_XPTS = Float32[
    5.0 15.0;  5.0 15.0;  5.0 15.0;  5.0 15.0;  5.0 15.0;  5.0 15.0;  5.0 15.0;
    5.0 15.0;  5.0 15.0; 15.0 30.0; 15.0 30.0; 30.0 60.0; 45.0 100.0]

# wc/fmcba.f DATA MAPFGS(139) — each R6 habitat code (ITYPE) → forb(1)/grass(2)/shrub(3), else 0 (→grass).
const _WC_MAPFGS = Int[
    0,0,0,0,0,0,0,3,3,3, 2,2,3,3,3,1,1,1,1,1, 1,1,1,1,3,3,3,3,3,3,
    3,3,3,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,3,3,3,
    3,3,3,3,3,3,3,3,3,3, 3,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,3,3,3,3,3,
    3,3,3,3,3,3,3,3,3,1, 3,3,3,3,3,3,3,3,3,3,
    3,3,3,3,3,3,3,3,3,3, 3,3,1,0,3,3,3,3,3,3, 3,3,3,3,3,3,3,3,3,3,
    3,3,3,0,3,3,3,3,3]

# wc/fmcba.f DATA MAPDRY(139) — each R6 habitat code → dry(0)/mesic(1)/moist(2). LWET = (IWET == 2).
const _WC_MAPDRY = Int[
    0,1,1,1,1,2,0,0,1,0, 0,0,0,0,0,0,2,2,2,1, 1,0,2,2,1,1,1,1,1,0,
    1,0,1,2,2,0,1,2,0,0, 2,1,1,1,1,1,1,1,2,2,
    2,2,2,0,0,2,1,0,1,0, 0,2,2,2,2,1,1,2,2,0, 2,2,1,2,2,0,2,1,1,1,
    1,1,1,1,1,1,0,0,0,0, 1,0,1,0,0,0,0,0,2,1,
    2,2,2,2,2,2,0,2,1,1, 1,1,2,2,0,1,0,1,2,2, 2,1,0,0,2,1,1,2,0,2,
    2,1,1,2,2,1,0,0,0]

# pn/fmcba.f DATA MAPFGS(75) / MAPDRY(75) — PN's own R6-habitat forb/grass/shrub + moist maps (DIFFER from
# WC's 139-code arrays; pn/fmcba.f MXR6CODE=75). fmcfmd.f itself is byte-identical to WC (same XPTS, same
# _wc_covgrp cover-metagroup pooling, same ALGSLP rules) — only these fmcba.f habitat arrays are PN-specific.
const _PN_MAPFGS = Int[
    3,3,3,3,3,3,3,1,1,1, 1,1,0,3,3,3,3,3,3,3, 3,3,3,3,3,1,1,1,1,1,
    3,1,0,1,3,3,3,3,3,3, 3,3,3,3,3,3,3,3,3,1,
    3,3,3,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,1,1,1, 3,3,3,3,3]
const _PN_MAPDRY = Int[
    0,0,0,0,0,0,0,2,0,0, 1,2,1,1,0,1,1,0,1,2, 1,1,2,0,1,2,2,1,2,2,
    2,0,1,2,1,1,0,1,0,0, 0,2,1,1,1,1,2,0,0,2,
    0,1,0,0,0,2,1,2,1,2, 2,1,1,0,2,2,2,2,2,2, 2,2,2,2,2]

# wc/fmcfmd.f:174-190 — pool species into the 6 cover metagroups (VARACD 'WC' branch).
@inline function _wc_covgrp(sp::Int)::Int
    (sp == 1 || sp == 19 || sp == 10 || sp == 18) && return 1   # SFCT: SF,WH,ES,RC
    (sp == 20 || sp == 4 || sp == 31)             && return 3   # MHCT: MH,AF,WB
    sp == 22                                       && return 4   # RACT: RA
    sp == 11                                       && return 5   # LPCT: LP
    (sp == 28 || sp == 25)                         && return 6   # WOCT: WO,GC
    return 2                                                     # DFCT: DF,GF,WP + everything else
end

# OP (Olympic NWO) species→cover-group map — op/fmcfmd.f the ELSE (VARACD≠'WC') branch. Differs from WC in
# 2 species (SS sp6→SFCT vs WC's ES sp10; sp24→WOCT). The selection RULES are BYTE-IDENTICAL to wc/fmcfmd.f.
@inline function _op_covgrp(sp::Int)::Int
    (sp == 1 || sp == 19 || sp == 6 || sp == 18)  && return 1   # SFCT: SF,WH,SS,RC
    (sp == 20 || sp == 4 || sp == 31)             && return 3   # MHCT: MH,AF,WB
    sp == 22                                       && return 4   # RACT: RA
    sp == 11                                       && return 5   # LPCT: LP
    (sp == 28 || sp == 24 || sp == 25)             && return 6   # WOCT: WO,TO,GC
    return 2                                                     # DFCT + everything else
end

"""
    wc_select_fuel_models(s, mois, sm, lg) -> Vector{(model, weight)}

WC FIRE-VPN cover-type-metagroup fuel-model selection (wc/fmcfmd.f). `sm`/`lg` are the SMALL/LARGE down-wood
loads (the FMDYN point); PERCOV comes from the FireState (fmcba). Source-faithful transcription.
"""
function wc_select_fuel_models(s::StandState, mois::AbstractMatrix{Float32}, sm::Float32, lg::Float32)
    t = s.trees
    eqwt = zeros(Float32, 13)                                  # WC ICLSS = 13
    percov = s.fire.percov
    itype = Int(s.plot.habitat_input)

    # PNFGS / PNWET (wc/fmcba.f entry points): habitat forb/grass/shrub + moist flags. PN carries its own
    # 75-code habitat arrays (pn/fmcba.f); WC uses the 139-code arrays. fmcfmd.f logic is shared (identical).
    _mapfgs = (s.variant isa PacificNorthwest || s.variant isa Olympic) ? _PN_MAPFGS : _WC_MAPFGS  # OP calls PNFGS
    _mapdry = (s.variant isa PacificNorthwest || s.variant isa Olympic) ? _PN_MAPDRY : _WC_MAPDRY
    icov = (1 <= itype <= length(_mapfgs)) ? _mapfgs[itype] : 0
    lforb  = icov == 1
    lgrass = icov == 0 || icov == 2
    lshrub = icov == 3
    iwet = (1 <= itype <= length(_mapdry)) ? _mapdry[itype] : 0
    lwet = iwet == 2

    # per-species FFE basal area (FMTBA) and cover-group BA (CTBA[1..6]).
    nsp = length(coef_col(s.coef, :dbh_min))
    fmtba = zeros(Float32, nsp)
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        sp = Int(t.species[i])
        fmtba[sp] += t.tpa[i] * t.dbh[i] * t.dbh[i] * 0.0054542f0
    end
    ctba = zeros(Float32, 6); stndba = 0f0
    @inbounds for sp in 1:nsp
        fmtba[sp] > 0f0 || continue
        stndba += fmtba[sp]
        ctba[(s.variant isa Olympic ? _op_covgrp(sp) : _wc_covgrp(sp))] += fmtba[sp]
    end

    # top-2 cover groups by BA (RDPSRT descending), rescaled to sum 1 (wc/fmcfmd.f:213-233).
    ict = [0, 0]; ictwt = Float32[0f0, 0f0]
    if t.n > 0 && stndba > 0.001f0
        idx = collect(1:6)
        rdpsrt!(6, ctba, idx, true)
        ict[1] = idx[1]; ict[2] = idx[2]
        x1 = ctba[idx[1]] + ctba[idx[2]]
        ictwt[1] = ctba[idx[1]] / x1
        ictwt[2] = ctba[idx[2]] / x1
    end

    # QMD80: QMD of the lower-80%-BA trees, computed only if SF or DF is a top-2 cover (wc/fmcfmd.f:238-270).
    lqmd = (ict[1] == 1 || ict[1] == 2 || ict[2] == 1 || ict[2] == 2)
    qmd80 = 0f0
    if lqmd
        targba = stndba * 0.80f0
        sumba = 0f0; ltrees = 0f0; qsum = 0f0
        ixs = collect(1:t.n)
        rdpsrt!(t.n, t.dbh, ixs, true)                        # descending by DBH
        @inbounds for j in t.n:-1:1                           # walk from smallest DBH up
            i = ixs[j]
            t.tpa[i] > 0f0 || continue
            treeba = t.tpa[i] * t.dbh[i] * t.dbh[i] * 0.0054542f0
            if (sumba + treeba) < targba
                sumba += treeba
                qsum  += t.tpa[i] * t.dbh[i] * t.dbh[i]
                ltrees += t.tpa[i]
            else
                sng1 = treeba > 0f0 ? (targba - sumba) / treeba : 0f0
                ltrees += t.tpa[i] * sng1
                qsum   += t.dbh[i] * t.dbh[i] * t.tpa[i] * sng1
                break
            end
        end
        ltrees > 0f0 && (qmd80 = sqrt(qsum / ltrees))
    end

    # ALGSLP weighting rules per top-2 cover group (wc/fmcfmd.f:282-489).
    alg(v, x1, x2) = _cr_algslp2(Float32(v), Float32(x1), Float32(x2), 0f0, 1f0)
    @inbounds for ii in 1:2
        j = ict[ii]; w1 = ictwt[ii]
        (j >= 1 && w1 > 0f0) || continue
        if j == 1 || j == 2                                   # SFCT / DFCT
            wt2b = alg(qmd80, 3.0, 5.0); wt2a = 1f0 - wt2b
            wt2a > 0f0 && (eqwt[5] += w1 * wt2a)
            if wt2b > 0f0
                if j == 1                                     # SFCT
                    if percov < 60f0
                        wt3b = alg(percov, 40.0, 60.0); wt3a = 1f0 - wt3b
                        wt3a > 0f0 && (eqwt[5] += w1 * wt2b * wt3a)
                        wt3b > 0f0 && (lshrub ? (eqwt[5] += w1 * wt2b * wt3b) : (eqwt[8] += w1 * wt2b * wt3b))
                    else
                        wt3c = alg(percov, 70.0, 90.0); wt3b = 1f0 - wt3c
                        wt3b > 0f0 && (eqwt[8] += w1 * wt2b * wt3b)
                        wt3c > 0f0 && (eqwt[8] += w1 * wt2b * wt3c)
                    end
                else                                          # DFCT
                    wt3b = alg(percov, 70.0, 90.0); wt3a = 1f0 - wt3b
                    if wt3a > 0f0
                        if lgrass;     eqwt[2] += w1 * wt2b * wt3a
                        elseif lshrub; eqwt[5] += w1 * wt2b * wt3a
                        elseif lforb;  eqwt[8] += w1 * wt2b * wt3a
                        end
                    end
                    wt3b > 0f0 && (eqwt[8] += w1 * wt2b * wt3b)
                end
            end
        elseif j == 3                                         # MHCT
            wt2b = alg(percov, 70.0, 90.0); wt2a = 1f0 - wt2b
            wt2a > 0f0 && (lshrub ? (eqwt[5] += w1 * wt2a) : (eqwt[8] += w1 * wt2a))
            wt2b > 0f0 && (eqwt[8] += w1 * wt2b)
        elseif j == 4                                         # RACT
            wt2b = alg(percov, 40.0, 60.0); wt2a = 1f0 - wt2b
            wt2a > 0f0 && (eqwt[5] += w1 * wt2a)
            wt2b > 0f0 && (eqwt[9] += w1 * wt2b)
        elseif j == 5                                         # LPCT
            if lwet
                wt2b = alg(percov, 40.0, 60.0); wt2a = 1f0 - wt2b
                wt2a > 0f0 && (eqwt[5] += w1 * wt2a)
                wt2b > 0f0 && (eqwt[8] += w1 * wt2b)
            else
                eqwt[5] += w1
            end
        elseif j == 6                                         # WOCT
            wt2b = alg(percov, 40.0, 60.0); wt2a = 1f0 - wt2b
            if wt2a > 0f0
                lshrub ? (eqwt[2] += w1 * wt2a) : (eqwt[1] += w1 * wt2a)
            end
            if wt2b > 0f0
                if lshrub;     eqwt[5] += w1 * wt2b
                elseif lgrass; eqwt[2] += w1 * wt2b
                elseif lforb;  eqwt[8] += w1 * wt2b
                end
            end
        end
    end

    # post-activity model 11 (AFWT/LATFUEL) — deferred (no activity ⇒ AFWT=0); models 10/12/13 natural fuels.
    eqwt[10] = 1f0
    eqwt[12] = 1f0
    eqwt[13] = 1f0

    return _fmdyn(sm, lg, eqwt, _WC_FMD_XPTS)
end
