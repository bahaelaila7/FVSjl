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


# _so_sospdm(s) — so/fmcfmd.f SOSPDM: the jittered dominant-species weights WDOM(1:MAXSP). 50 passes; each draws RANN
# once per tree record (I = 1..ITRN) for DT = D·(1 + 0.2·(2·XRAN − 1)), sums DT²·P by species, and credits the species
# with the most (LP 7 when there is no BA; ties keep the lower index). The main RANN stream is saved and restored.
function _so_sospdm(s::StandState)
    t = s.trees; nsp = 33
    wdom = zeros(Float32, nsp); xba = zeros(Float32, nsp)
    saved = rannget(s.rng)
    for _ in 1:50
        fill!(xba, 0f0)
        @inbounds for i in 1:t.n
            xran = rann!(s.rng)
            dt = t.dbh[i] * (1f0 + 0.2f0 * ((xran * 2f0) - 1f0))
            k = Int(t.species[i])
            (1 <= k <= nsp) && (xba[k] = xba[k] + (dt * dt * t.tpa[i]))
        end
        basp = 7
        @inbounds for k in 1:nsp
            xba[k] > xba[basp] && (basp = k)
        end
        wdom[basp] += 1f0
    end
    rannput!(s.rng, saved)
    @inbounds for i in 1:nsp; wdom[i] = wdom[i] / 50f0; end
    return wdom
end

# so/fmcfmd.f:248-411 — the CALIFORNIA (Region-5, KODFOR 500-599 or 701) rules: per SOSPDM dominant-species weight WD,
# the species group's percent-cover / average-height / snag-density / relative-density fuel-model weights.
function _so_r5_eqwt!(eqwt::Vector{Float32}, s::StandState, percov::Float32)
    t = s.trees; fs = s.fire
    alg(x, x1, x2) = _cr_algslp2(Float32(x), Float32(x1), Float32(x2), 0f0, 1f0)
    wd = _so_sospdm(s)
    # RMSQD, AVH, RELDEN / RELDSP(7,9) of the current tree list (dense.f / avht40.f quantities)
    tp = 0f0; sd2 = 0f0; relden = 0f0; rel7 = 0f0; rel9 = 0f0
    ifor = Int(s.plot.forest_idx)
    @inbounds for i in 1:t.n
        p = t.tpa[i]; d = t.dbh[i]; tp += p; sd2 += d * (d * p)
        c = so_tree_ccf(Int(t.species[i]), d, t.height[i]; ifor = ifor) * p
        relden += c
        t.species[i] == 7 && (rel7 += c); t.species[i] == 9 && (rel9 += c)
    end
    rmsqd = tp > 0f0 ? sqrt(sd2 / tp) : 0f0
    avh = s.plot.avg_height
    @inbounds for i in 1:33
        w = wd[i]; w > 0f0 || continue
        if i == 1 || i == 2 || i == 6 || i == 11
            w2 = alg(percov, 25, 35); w1 = 1f0 - w2
            w1 > 0f0 && (eqwt[2] = eqwt[2] + w1 * w)
            w2 > 0f0 && (eqwt[9] = eqwt[9] + w2 * w)
        elseif i == 10 || i == 4 || i == 3 || i == 5 || i == 8
            xsng = 0f0
            sn = fs.snags
            for j in eachindex(sn.dbh)
                sn.dbh[j] > rmsqd && (xsng = xsng + sn.den_hard[j] + sn.den_soft[j])
            end
            a2 = alg(avh, 30, 40); a1 = 1f0 - a2
            b2 = alg(percov, 55, 65); b1 = 1f0 - b2
            c2 = alg(xsng, 3.5, 4.5); c1 = 1f0 - c2
            if a2 > 0f0 || b2 > 0f0
                if c1 > 0f0
                    eqwt[9] = eqwt[9] + w * a2 * b1 * c1
                    eqwt[9] = eqwt[9] + w * a2 * b2 * c1
                    eqwt[9] = eqwt[9] + w * a1 * b2 * c1
                end
                if c2 > 0f0
                    eqwt[10] = eqwt[10] + w * a2 * b1 * c2
                    eqwt[10] = eqwt[10] + w * a2 * b2 * c2
                    eqwt[10] = eqwt[10] + w * a1 * b2 * c2
                end
            end
            (a1 > 0f0 && b1 > 0f0) && (eqwt[5] = eqwt[5] + w * a1 * b1)
        elseif i == 7 || i == 9
            xa = relden > 1f-6 ? max(rel7, rel9) / relden : 1f0
            w12 = alg(xa, 0.85, 0.95); w11 = 1f0 - w12
            w12 > 0f0 && (eqwt[8] = eqwt[8] + w * w12)
            if w11 > 0f0
                xm = 0f0; xmlp = 0f0
                for j in 1:t.n
                    x = t.mort_pa[j] * t.dbh[j]^2
                    # fmcfmd.f tests ISP(I) — the species of RECORD I (I is the species loop index), not of record J
                    (i <= t.n && t.species[i] == 7) && (xmlp += x)
                    xm += x
                end
                xr = xm > 1f-6 ? xmlp / xm : 0f0
                v2 = alg(xr, 0.45, 0.55); v1 = 1f0 - v2
                v1 > 0f0 && (eqwt[8] = eqwt[8] + w * w11 * v1)
                v2 > 0f0 && (eqwt[10] = eqwt[10] + w * w11 * v2)
            end
        else
            eqwt[8] = eqwt[8] + w                              # FM 8: R5 default when no listed species dominates
        end
    end
    return eqwt
end

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
    kf = Int(s.plot.user_forest_code)
    if (500 <= kf < 600) || kf == 701                           # CALIFORNIA (so/fmcfmd.f:248)
        _so_r5_eqwt!(eqwt, s, percov)
        eqwt[10] = 1f0; eqwt[12] = 1f0; eqwt[13] = 1f0
        return _fmdyn(sm, lg, eqwt, _SO_FMD_XPTS)
    end
    itype = Int(s.plot.habitat_input); itype <= 0 && (itype = 49)   # so/habtyp.f DEFAULT CPS111 = ITYPE 49
    ipag = (1 <= itype <= 92) ? Int(_SO_IPASO[itype]) : 1
    alg(x, x1, x2) = _cr_algslp2(x, Float32(x1), Float32(x2), 0f0, 1f0)   # ALGSLP(x,[x1,x2],[0,1],2)

    # FMSSTAGE structure class (IFMST) — only the LP (6,7) and juniper (8) groups use it.
    ifmst = (ipag == 6 || ipag == 7 || ipag == 8) ? Int(structure_class(s).class) : 0
    # ACTCBH (active crown base height) feeds the LP groups (6,7). It is computed at fire time (fmcfir), not
    # persisted — the R6 PP reference stand (IPAG=1) does not use it, so 0 here; the LP/juniper ACTCBH path is
    # a documented follow-on (same deferral class as the R5-California SOSPDM branch).
    actcbh = 0f0
    # DSTLG = MIN(IYR−HARVYR, IYR−BURNYR) (so/fmcfmd.f:415): years since the last harvest (fmscut.f HARVYR, 0 if none)
    # or the last carried burn with SCH > PBSCOR (fmburn.f BURNYR, −1 if none). MEASURED FVSso_g16 15184869010497
    # simfire 2020 (the fire year, PERCOV 40.5 < 45): live Fuel_Mod3 = 2 (DSTLG 0 ⇒ WT2=0 ⇒ FM2), jl FM6 from DSTLG=999.
    dstlg = let yr = Int(current_cycle_year(s)), burnyr = -1, fs = s.fire
        for br in fs.burn_reports
            (get(br, :carried, true)::Bool && (br.scorch::Float32) > fs.params.pb_scor) && (burnyr = max(burnyr, Int(br.year::Int)))
        end
        Float32(min(yr - Int(fs.harvyr), yr - burnyr))
    end

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
