# =============================================================================
# diameter_growth.jl (centralstates) — CS large-tree diameter growth (cs/dgf.f)
#
# The one genuinely-new CS model: an SN-family ln(DDS) regression (NOT NE's BAL-
# potential iteration). DDS (Δ squared inside-bark DBH) is predicted from DBH, its
# inverse and square, the crown ratio and its square, stand BA in trees ≥5", the
# BAL (BA in larger trees, from the BA percentile), site index, and two QMD-relative
# diameter terms (D/QMD≥5 and D²/QMD≥5). Coefficients: data/centralstates/dg_coeffs.csv
# (cs/dgf.f DATA). The bark inside-out conversion is done INSIDE the routine (the CS
# regression predicts OUTSIDE-bark Δd²), then WK2 holds ln(inside-bark DDS) — the same
# convention SN/NE downstream (calibrate / DGDRIV growth) consume.
#
#   DDS_ob = DGCON+COR + INTERC + VDBHC/D + DBHC·D + DBH2C·D² + RDBHC·(D/QMD5)
#            + RDBHSQC·(D²/QMD5) + CRWNC·CR + CRSQC·CR² + SBAC·BAGE5 + BALC·BAL
#            + SITEC·SITEAR                                              (cs/dgf.f:510)
#   DIAGRO = √(D² + e^DDS_ob) − D ;  DIAGRI = DIAGRO·BRATIO              (:547)
#   WK2    = ln((D·BR+DIAGRI)² − (D·BR)²)                               (:550)
# =============================================================================

# Per-species QMD≥5 caps (cs/dgf.f:460-478) and crown-ratio caps (:482-494), by species index.
# The QMD cap MUTATES the stand-wide QMDGE5 IN PLACE as the Fortran walks species in index order, so
# a cap-species lowers it for every LATER species — dgf! replicates that cumulative walk (see below),
# `_cs_qmd_cap` is the single-species cap kernel it folds. (Crown caps stay per-tree-local, as in the
# Fortran.) The old per-tree-local QMD cap was bit-identical only with ≤1 capping species (cst01); the
# LS site sweep proved it wrong for multi-cap-species stands, so both CS and LS now do the cumulative walk.
@inline function _cs_qmd_cap(spc::Int, q::Float32)::Float32
    spc == 50 && return min(q, 12f0)
    (spc == 3 || (10 <= spc <= 13)) && return min(q, 13f0)
    ((14 <= spc <= 17) || spc == 28 || (53 <= spc <= 56)) && return min(q, 25f0)
    spc == 24 && return min(q, 40f0)
    spc == 51 && return min(q, 11f0)
    ((44 <= spc <= 46) || spc == 59) && return min(q, 20f0)
    (48 <= spc <= 49) && return min(q, 30f0)
    spc == 91 && return min(q, 17f0)
    return q
end

@inline function _cs_cr_cap(spc::Int, cr::Float32)::Float32
    spc == 7 && return min(cr, 50f0)
    (8 <= spc <= 10) && return min(cr, 75f0)
    (spc == 28 || spc == 41) && return min(cr, 60f0)
    (44 <= spc <= 46) && return min(cr, 80f0)
    (spc == 32 || (78 <= spc <= 84)) && return min(cr, 85f0)
    return cr
end

"""
    dgf!(s, ::CentralStates)

CS variant `dgf!` hook — fill `scratch.wk[2, i]` with ln(inside-bark DDS) for every live
tree (the cs/dgf.f regression + its OB→IB bark conversion). COR enters through `conspp`
(DGCON+COR), exactly like SN's `dgf!`. Consumed by the shared DGDRIV calibrate / growth pass.
"""
function dgf!(s::StandState, ::CentralStates)
    p, t, c, sd = s.plot, s.trees, s.calib, s.coef.species
    wk2 = view(s.scratch.wk, 2, :)
    interc = sd[:dg_interc]; vdbh = sd[:dg_vdbh]; dbhc = sd[:dg_dbh]; dbh2 = sd[:dg_dbh2]
    rdbh = sd[:dg_rdbh]; rdbhsq = sd[:dg_rdbhsq]; crwn = sd[:dg_crwn]; crsq = sd[:dg_crsq]
    sbac = sd[:dg_sba]; balc = sd[:dg_bal]; sitec = sd[:dg_site]
    ba_a = c.bark_a; ba_b = c.bark_b

    # BAGE5 (stand BA in trees ≥5") and QMDGE5 (their QMD) — cs/dgf.f:415-436.
    sdqge5 = 0f0; tt = 0f0; bage5 = 0f0
    @inbounds for i in 1:t.n
        d = t.dbh[i]
        d < 5f0 && continue
        pr = t.tpa[i]
        sdqge5 += pr * d * d
        tt += pr
        bage5 += d * d * pr * 0.005454154f0
    end
    qmdge5 = sdqge5 > 0f0 ? sqrt(sdqge5 / tt) : 0f0
    bage5 <= 0f0 && (bage5 = 10f0)
    ba_v = p.basal_area

    # Per-species effective QMDGE5 — FVS mutates the stand-wide QMDGE5 IN PLACE as it walks species in
    # INDEX order (cs/dgf.f:460-478, DO 200 ISPC=1,MAXSP), so a species sees the value capped by ALL
    # LOWER-INDEXED present cap-species. Formerly a per-tree LOCAL cap (bit-identical only with ≤1 cap
    # species, true for cst01) — the LS sweep proved it wrong for multi-cap-species stands (an uncapped
    # species missed a lower-indexed cap-species' mutation). Replicate the cumulative walk. No-op for a
    # single capping species; a no-op entirely for cst01.
    nsp = nspecies(s.variant)
    present = falses(nsp)
    @inbounds for i in 1:t.n
        (t.dbh[i] > 0f0 && t.tpa[i] > 0f0) && (present[Int(t.species[i])] = true)
    end
    qmd_eff = fill(qmdge5, nsp)
    qcur = qmdge5
    @inbounds for spc in 1:nsp
        present[spc] && (qcur = _cs_qmd_cap(spc, qcur))
        qmd_eff[spc] = qcur
    end

    @inbounds for i in 1:t.n
        d = t.dbh[i]
        d <= 0f0 && continue
        sp = Int(t.species[i])
        conspp = c.dg_const[sp] + c.dg_cor[sp]                 # DGCON + COR
        cr = Float32(t.crown_pct[i]);  cr <= 0f0 && (cr = 10f0)
        cr = _cs_cr_cap(sp, cr)
        q  = qmd_eff[sp]                                       # cumulative species-order cap (see above)
        reldbh = q > 0f0 ? d / q : 0f0
        reldbhsq = q > 0f0 ? d * d / q : 0f0
        bal = (1f0 - t.crown_ratio[i] / 100f0) * ba_v          # PCT = BA percentile

        dds = conspp + interc[sp] +
              vdbh[sp]   * (1f0 / d) +
              dbhc[sp]   * d +
              dbh2[sp]   * d * d +
              rdbh[sp]   * reldbh +
              rdbhsq[sp] * reldbhsq +
              crwn[sp]   * cr +
              crsq[sp]   * cr * cr +
              sbac[sp]   * bage5 +
              balc[sp]   * bal +
              sitec[sp]  * p.sp_site_index[sp]

        dds < -9.21f0 && (dds = -9.21f0)
        # OB→IB bark conversion (cs/dgf.f:547-550).
        diagro = sqrt(d * d + exp(dds)) - d
        bark = bark_ratio(ba_a, ba_b, sp, d)
        diagri = diagro * bark
        db = d * bark
        dds = log((db + diagri)^2 - db * db)
        dds < -9.21f0 && (dds = -9.21f0)
        wk2[i] = dds
    end
    return s
end

"""
    cs_dgcons!(s)

CS per-stand DG setup (cs/dgf.f ENTRY DGCONS): DGCON=0 (the intercept lives in the INTERC
coefficient; site dependence is the SITEC·SITEAR term, not a precomputed constant), ATTEN ←
OBSERV (per-species observation count, the DGDRIV calibration prior weight), SMCON=0, plus the
per-stand bark copy (`calib.bark_a/b` ← intercept 0 / slope BKRAT) the CFTOPK + DGDRIV read.
"""
function cs_dgcons!(s::StandState)
    c = s.calib; sd = s.coef.species; ctl = s.control
    ba = sd[:bark_intercept]; bb = sd[:bark_slope]; obs = sd[:dg_observ]
    @inbounds for sp in 1:MAXSP
        c.bark_a[sp] = ba[sp]; c.bark_b[sp] = bb[sp]
        c.dg_const[sp] = 0f0
        c.atten[sp] = obs[sp]              # ATTEN = OBSERV(ISPC) — cs/dgf.f:590
        # READCORD/REUSCORD (LDCOR2): add ln(COR2) to DGCON (cs/dgf.f:597-598), same as the SN dgcons! path.
        ctl.dg_cor2_on && ctl.dg_cor2[sp] > 0f0 && (c.dg_const[sp] += log(ctl.dg_cor2[sp]))
    end
    return s
end

"""
    _cs_init_crowns!(s)

CS CRATET: dub the INITIAL crown ratio for inventory trees with no input crown (`crown_pct==0`),
using the eastern (TWIGS) crown model, so the calibration `dgf!` sees a real crown (its `CRWNC·CR +
CRSQC·CR²` terms) rather than the 0→10 fallback. FVS runs this dub AFTER DGDRIV's DENSE backdate, so
the stand BA the crown model reads is the BACKDATED-dbh total per-acre BA (crown.f COMMON `BA`) —
computed here on the backdated diameters, then the current dbh is restored (the crown formula's own D
stays the current dbh, matching crown.f). No-op if every tree already has a crown. (SN's analogue is
`init_crown_ratios!`; NE's DGF uses BAL not crown, so NE needs no crown at calibration.)
"""
function _cs_init_crowns!(s::StandState)
    t = s.trees; n = t.n
    n == 0 && return s
    any(@views t.crown_pct[1:n] .== 0) || return s
    saved_dbh = Float32[t.dbh[i] for i in 1:n]
    _backdate_dbh!(s)
    bd_ba = 0f0
    @inbounds for i in 1:n; d = t.dbh[i]; bd_ba += d * d * t.tpa[i] * 0.005454154f0; end
    @inbounds for i in 1:n; t.dbh[i] = saved_dbh[i]; end
    crown_ratio_update!(s, s.variant; ba_override = bd_ba, lstart = true)
    return s
end
