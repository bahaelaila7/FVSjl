# =============================================================================
# standstats.jl — per-acre expansion (NOTRE) and stand summary statistics
#
# Ported from: base/notre.f (expansion) + base/stats.f / sumout (summary columns).
#
# NOTRE turns each record's tally (PROB) into trees-per-acre using the sampling
# design (variable-radius BAF plots vs fixed plots). The stand statistics
# (TPA, basal area, QMD, ...) are then simple weighted reductions over the trees —
# pure loops, autovectorizable, no allocation.
# =============================================================================

const BA_PER_TREE = 0.005454154f0     # ft² of basal area per inch² of DBH

"""
    finalize_design!(state)

INITRE end-of-keywords finalization (initre.f:360-414). The crucial, easy-to-miss
step: FVS overwrites the variable named `PI` with the plot count `IPTINV` — that is
the divisor NOTRE actually uses (NOT π). Also finalizes the stockable proportion
GROSPC into the reciprocal multiplier NOTRE applies.
"""
function finalize_design!(s::StandState)
    p = s.plot
    p.points_inv <= 0 && (p.points_inv = Int32(1))
    p.pi = Float32(p.points_inv)                      # PI := IPTINV (NOTRE divisor!)
    p.sample_weight < 0f0 && (p.sample_weight = Float32(p.points_inv))
    if p.gross_space < 0f0                             # default (set to -1 at stand start)
        g = (p.pi - Float32(p.nonstockable)) / p.pi
        g > 1f0 && (g = 1f0)
        (p.pi - Float32(p.nonstockable)) <= 0f0 && (g = 1f0)
        p.gross_space = g
    end
    p.gross_space = 1f0 / p.gross_space               # → reciprocal multiplier
    return s
end

"""
    notre!(state)

Compute trees-per-acre (`trees.tpa`) for every live record from the sampling
design. Mirrors notre.f: variable-plot trees (DBH ≥ BRK) scale by BAF/DBH²,
fixed-plot trees (DBH < BRK) by FPA, all times the stockable proportion GROSPC.
"""
function notre!(s::StandState)
    p, t = s.plot, s.trees
    fp  = p.fixed_plot_inv / p.pi
    p.total_fixed_plot > 0f0 && (fp = 1f0 / p.total_fixed_plot)
    vp  = p.baf * 183.3465f0 / p.pi
    fp2 = 0f0
    if p.baf <= 0f0
        vp = 0f0; fp2 = -p.baf / p.pi
    end
    brk = p.min_dbh_var_plot
    # expand live records and the dead partition (n+1:n+ndead) alike — dead trees
    # carry their expanded TPA into the backdated calibration BA.
    @inbounds for i in 1:(t.n + t.ndead)
        pr = t.tpa[i]; d = t.dbh[i]
        pr <= 0f0 && (pr = 1f0)
        pr = d < brk ? pr * fp : pr * vp / (d * d) + pr * fp2
        pr <= 0f0 && (pr = 9.0f-25)
        t.tpa[i] = pr * p.gross_space
    end
    # base species-sort key = original record index (FVS chain order pre-tripling);
    # tripling derives child keys so the per-tree RNG draw order matches the oracle.
    @inbounds for i in 1:(t.n + t.ndead)
        t.sort_key[i] = Float64(i)
    end
    return s
end

"""
    stand_tpa(state)   -> total trees per acre        (TPROB)
    stand_ba(state)    -> basal area ft²/acre         (BA)
    stand_qmd(state)   -> quadratic mean diameter, in (RMSQD)

Weighted reductions over the live trees. Pure (no mutation).
"""
function stand_tpa(s::StandState)
    t = s.trees; tot = 0f0
    @inbounds for i in 1:t.n; tot += t.tpa[i]; end
    return tot
end

function stand_ba(s::StandState)
    t = s.trees; ba = 0f0
    @inbounds for i in 1:t.n; ba += t.tpa[i] * BA_PER_TREE * t.dbh[i]^2; end
    return ba
end

function stand_qmd(s::StandState)
    t = s.trees; sd2 = 0f0; tpa = 0f0
    @inbounds for i in 1:t.n
        sd2 += t.tpa[i] * t.dbh[i]^2
        tpa += t.tpa[i]
    end
    return tpa > 0f0 ? sqrt(sd2 / tpa) : 0f0
end

"""
    stand_sdimax(s) -> Float32

BA-weighted stand maximum SDI (`SDICAL`, base/sdical.f, pre-CLMAXDEN). General across variants —
the per-species SDImax (`plot.sp_sdi_def`) is variant coefficient data; the averaging is the same
base algorithm. Used by the mortality SDImax cap and the structure-stage PCTSMX demotion (BTSDIX).
"""
function stand_sdimax(s::StandState)
    t = s.trees; p = s.plot
    num = 0f0; totba = 0f0
    @inbounds for i in 1:t.n
        tb = 0.0054542f0 * t.dbh[i]^2 * t.tpa[i]
        num   += p.sp_sdi_def[t.species[i]] * tb
        totba += tb
    end
    return totba <= 0f0 ? 1f0 : num / totba
end

"""
    stand_top_height(state)

Average height of the largest-diameter 40 trees/acre (AVHT40, the summary "top
height"). Trees are taken in descending-DBH order; the last one is prorated to
hit exactly 40 TPA. (Uses a sort — fine for once-per-cycle stats, not the hotpath.)
"""
function stand_top_height(s::StandState)
    t = s.trees
    t.n == 0 && return 0f0
    # avht40.f sorts IND with FVS's RDPSRT (Scowen quickersort, descending DBH) — NOT a stable sort. The
    # tie-break among equal-DBH trees decides WHICH tree lands at the 40-TPA boundary (and so its height
    # enters AVH), so a stable `sortperm!` (ascending-index ties) diverges from live on tie-heavy stands.
    # (The DG `point_basal_area!` also sorts by DBH but its BAL is an order-independent sum ⇒ tie order is
    # inert there; only AVH exposes it.) Use the ported `_rdpsrt!` to match FVS's IND tie-break exactly.
    idx = view(s.scratch.stat_idx, 1:t.n)
    dbhv = view(t.dbh, 1:t.n)
    # FVS cratet.f computes the AVHT40 top-height IND by a DOUBLE sort: IND1 = a fresh RDPSRT(.TRUE.), then
    # RDPSRT(.FALSE.) re-sorts that pre-ordered IND. RDPSRT is unstable, so the `.FALSE.` pass SWAPS equal-DBH
    # ties (the later-read record lands at the 40-tpa boundary), which a single sort does NOT — this is the
    # cycle-0 top-height divergence on tie-heavy stands (equal-DBH, different-height trees at the boundary).
    _rdpsrt!(dbhv, idx)                 # LSEQ=.TRUE. → IND1
    _rdpsrt!(dbhv, idx; lseq = false)   # LSEQ=.FALSE. → re-sort preserving IND1 (swaps ties, matches FVS)
    avh = 0f0; ssumn = 0f0
    for k in 1:t.n
        ii = Int(idx[k])
        p = t.tpa[ii]
        ssumn + p > 40f0 && (p = 40f0 - ssumn)
        ssumn += p
        avh += t.height[ii] * p
        ssumn >= 40f0 && break
    end
    return ssumn > 0f0 ? avh / ssumn : 0f0
end

"""
    point_basal_area!(state)

Fill `density.point_ba[ip]` (PTBAA, per-point basal area) AND `density.point_bal[i]`
(PTBALT, the BA in trees LARGER than tree i on the same point). For each subplot,
trees are taken in descending-DBH order and BA accumulates; PTBALT[i] is the sum
before tree i. The diameter-growth competition term uses PTBALT (= pbal). Per-tree
BA = tpa·0.005454154·DBH²·PI/GROSPC. (Sorts per point — once-per-cycle, not hotpath.)
"""
function point_basal_area!(s::StandState)
    p, t = s.plot, s.trees
    pb = s.density.point_ba; pbal = s.density.point_bal
    fill!(pb, 0f0)
    scale = p.pi / p.gross_space
    npts = 0
    @inbounds for i in 1:t.n
        npts = max(npts, Int(t.plot_id[i]))
        pbal[i] = 0f0
    end
    order = sortperm!(view(s.scratch.stat_idx, 1:t.n), view(t.dbh, 1:t.n); rev = true)  # descending DBH (stable)
    @inbounds for i in order
        ip = Int(t.plot_id[i])
        pbal[i] = pb[ip]                                # BA already accumulated = larger trees
        pb[ip] += t.tpa[i] * BA_PER_TREE * t.dbh[i]^2 * scale
    end
    return s
end

"""
    point_density!(state)

Fill `density.point_ccf[ip]` (PCCF) and `density.point_tpa[ip]` (PTPA) — the per-point crown
competition factor and trees-per-acre (dense.f:210-211). Each tree contributes its open-grown crown
area (the same CCFT as `stand_ccf`) and TPA to its OWN subplot, scaled by PI/GROSPC (`p.pi/p.gross_space`,
the same scale `point_basal_area!` uses) so the point value is the gross per-acre density on that point.
Consumed by the multi-point regen crown ratio (regent.f:178 `CR=0.89722−0.0000461·PCCF`) and the
TCONDMLT point weights (cuts.f:1074 `+PBAWT·PTBAA+PCCFWT·PCCF+PTPAWT·PTPA`). (Once per cycle, not hotpath.)
"""
function point_density!(s::StandState)
    p, t = s.plot, s.trees
    pccf = s.density.point_ccf; ptpa = s.density.point_tpa
    fill!(pccf, 0f0); fill!(ptpa, 0f0)
    pi_f = p.pi; gross = p.gross_space
    @inbounds for i in 1:t.n
        ip = Int(t.plot_id[i])
        (1 <= ip <= length(pccf)) || continue
        sp2 = s.species.code2[t.species[i]]
        cw  = crown_width(s.coef, sp2, t.dbh[i], t.height[i], 90, 1,
                          p.latitude, p.longitude, p.elevation)
        ccft = t.dbh[i] > 0.1f0 ? 0.001803f0 * cw * cw * t.tpa[i] : 0.001f0 * t.tpa[i]
        # dense.f:210-211 accumulates each term as `CCFT*PI/GROSPC` — i.e. (ccft·pi)/gross evaluated
        # left-to-right, NOT ccft·(pi/gross) with a precomputed reciprocal-scale. The two differ by ~1
        # Float32 ULP per term; on the dense estab_pccf points that sub-ULP tips a regen-crown INT(CR·100+0.5)
        # boundary. Match FVS's exact op order to deconfound (doctrine #8).
        pccf[ip] += ccft * pi_f / gross
        ptpa[ip] += t.tpa[i] * pi_f / gross
    end
    return s
end

"""
    stand_pct!(state)

Fill `trees.crown_ratio[i]` with PCT, the stand basal-area percentile (PCTILE,
pctile.f via dense.f): trees in descending-DBH order, `PCT[i] = (BA of tree i and
all smaller) / total · 100`. So `1 − PCT/100` is the fraction of stand BA in larger
trees, which the diameter-growth competition term uses. (Despite the field name,
this is FVS's PCT array, not the crown ratio — the crown ratio is `crown_pct`/ICR.)
"""
function stand_pct!(s::StandState)
    t = s.trees; n = t.n
    n == 0 && return s
    order = sortperm!(view(s.scratch.stat_idx, 1:n), view(t.dbh, 1:n); rev = true)  # largest DBH first
    pct = t.crown_ratio
    cum = 0f0
    @inbounds for k in n:-1:1                            # accumulate from smallest up
        ii = order[k]
        cum += t.dbh[ii]^2 * t.tpa[ii]
        pct[ii] = cum
    end
    if cum > 0f0
        @inbounds for ii in 1:n
            pct[ii] = pct[ii] / cum * 100f0
        end
    end
    return s
end

"""
    stand_ccf(state)

Crown competition factor (RELDEN): Σ over trees of the open-grown crown area
(CCFCAL/ccfcal.f): `0.001803·crownwidth²·tpa` (or `0.001·tpa` for DBH ≤ 0.1).
"""
function stand_ccf(s::StandState)
    p, t = s.plot, s.trees
    ccf = 0f0
    @inbounds for i in 1:t.n
        sp = t.species[i]
        sp2 = s.species.code2[sp]
        cw = crown_width(s.coef, sp2, t.dbh[i], t.height[i], 90, 1,
                         p.latitude, p.longitude, p.elevation)
        ccf += t.dbh[i] > 0.1f0 ? 0.001803f0 * cw * cw * t.tpa[i] : 0.001f0 * t.tpa[i]
    end
    return ccf
end

"""
    stand_sdi(s)

Reported `.sum` stand density index (SDICLS, sdical.f), following the SDICALC method flag
`zeide_sdi` (LZEIDE, SN default Zeide) — the SAME flag the SDImax mortality uses, so the two
stay consistent. **Zeide:** Σ TPA·(D/10)^1.605 over `D ≥ DBHZEIDE` (sdical.f:326). **Reineke:**
the `SDI = SPROB·A + B·SDSQ` Taylor form over `D ≥ DBHSTAGE` (sdical.f:281-327). Defaults
(Zeide, threshold 0) reproduce the prior behavior.
"""
function stand_sdi(s::StandState)
    t = s.trees
    if s.control.zeide_sdi
        thr = s.control.dbh_zeide; sdi = 0f0
        @inbounds for i in 1:t.n
            # sdical.f:326 `(DBH/10.)**1.605` — FVS `**` is gfortran powf, NOT Julia's openlibm `^` (differ ~0.07%);
            # route through the companion (doctrine #8) so the reported/MYSDI Zeide SDI matches FVS bit-exactly.
            t.dbh[i] >= thr && (sdi += t.tpa[i] * fpow(t.dbh[i] / 10f0, 1.605f0))
        end
        return sdi
    end
    return stand_sdi_reineke(s)
end

"Reineke/STAGE stand SDI (SDIC = SPROB*A + B*SDSQ, sdical.f:47-61/105) — the form FVS's CROWN uses."
function stand_sdi_reineke(s::StandState)
    t = s.trees
    thr = s.control.dbh_stage; sprob = 0f0; sdsq = 0f0
    @inbounds for i in 1:t.n
        if t.dbh[i] >= thr
            sprob += t.tpa[i]; sdsq += t.dbh[i]^2 * t.tpa[i]
        end
    end
    sprob <= 0f0 && return 0f0
    mdsq = sdsq / sprob
    # sdical.f:281-282 `(10.0**(-1.605))*…*((SDSQ/SPROB)**(1.605/2.))` — all FVS `**` = gfortran powf, route via
    # the companion (doctrine #8) not Julia's openlibm `^`. This feeds CROWN's SDI ⇒ keep SN/NE/CS/LS bit-exact.
    a = fpow(10f0, -1.605f0) * (1f0 - 1.605f0 / 2f0) * fpow(mdsq, 1.605f0 / 2f0)
    b = fpow(10f0, -1.605f0) * (1.605f0 / 2f0) * fpow(mdsq, 1.605f0 / 2f0 - 1f0)
    return sprob * a + b * sdsq
end
