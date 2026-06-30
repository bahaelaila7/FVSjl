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
# NOTE: in the Fortran the QMD cap mutates the stand-wide QMDGE5 IN PLACE inside the species loop,
# so a capping species would lower it for every later species (an order-dependent latent quirk).
# We apply the cap to a per-tree LOCAL copy instead — bit-identical whenever at most one capping
# species is present (true for cst01, whose species 8/19/43/47/60 hit no QMD cap). Revisit if a
# multi-capping-species stand ever diverges. (Crown caps are already per-tree-local in the Fortran.)
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

    @inbounds for i in 1:t.n
        d = t.dbh[i]
        d <= 0f0 && continue
        sp = Int(t.species[i])
        conspp = c.dg_const[sp] + c.dg_cor[sp]                 # DGCON + COR
        cr = Float32(t.crown_pct[i]);  cr <= 0f0 && (cr = 10f0)
        cr = _cs_cr_cap(sp, cr)
        q  = _cs_qmd_cap(sp, qmdge5)
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
    c = s.calib; sd = s.coef.species
    ba = sd[:bark_intercept]; bb = sd[:bark_slope]; obs = sd[:dg_observ]
    @inbounds for sp in 1:MAXSP
        c.bark_a[sp] = ba[sp]; c.bark_b[sp] = bb[sp]
        c.dg_const[sp] = 0f0
        c.atten[sp] = obs[sp]              # ATTEN = OBSERV(ISPC) — cs/dgf.f:590
    end
    return s
end
