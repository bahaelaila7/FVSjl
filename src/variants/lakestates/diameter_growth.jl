# =============================================================================
# diameter_growth.jl (lakestates) — LS large-tree diameter growth (ls/dgf.f)
#
# LS's large-tree DG is the SAME SN-family ln(DDS) regression as CS (NOT NE's BAL-
# potential iteration): DDS (Δ squared inside-bark DBH) predicted from DBH, its
# inverse and square, crown ratio + its square, stand BA in trees ≥5" (BAGE5), BAL
# (BA in larger trees, from the BA percentile), site index, and two QMD-relative
# diameter terms (D/QMD≥5, D²/QMD≥5). Coefficients: data/lakestates/dg_coeffs.csv
# (ls/dgf.f DATA — the 68-species INTERC/VDBHC/DBHC/DBH2C/RDBHC/RDBHSQC/CRWNC/CRSQC/
# SBAC/BALC/SITEC + OBSERV). The bark inside-out conversion is done INSIDE the routine;
# WK2 holds ln(inside-bark DDS) for the shared DGDRIV calibrate / growth pass.
#
#   DDS_ob = DGCON+COR + INTERC + VDBHC/D + DBHC·D + DBH2C·D² + RDBHC·(D/QMD5)
#            + RDBHSQC·(D²/QMD5) + CRWNC·CR + CRSQC·CR² + SBAC·BAGE5 + BALC·BAL
#            + SITEC·SITEAR                                            (ls/dgf.f:260-271)
#   DIAGRO = √(D² + e^DDS_ob) − D ;  DIAGRI = DIAGRO·BRATIO            (:288)
#   WK2    = ln((D·BR+DIAGRI)² − (D·BR)²)                              (:291)
#
# The model body is identical to CS `dgf!`; only the per-species QMD/CR caps differ
# (ls/dgf.f:236-251 SELECT CASE(ISPC)).
# =============================================================================

# Per-species QMD≥5 caps (ls/dgf.f:238-244) — like CS, the Fortran mutates the stand-wide QMDGE5 in
# place inside the species loop (an order-dependent latent quirk). We apply the cap to a per-tree LOCAL
# copy — bit-identical whenever at most one capping species is present (lst01 QMDGE5 stays below 13, so
# no cap fires). Revisit if a multi-capping-species stand ever diverges.
@inline function _ls_qmd_cap(spc::Int, q::Float32)::Float32
    (spc == 1 || spc == 2 || spc == 10 || (37 <= spc <= 39) || spc == 59) && return min(q, 13f0)  # JP,SC,TA,BH,PH,SH,BG
    spc == 11 && return min(q, 15f0)                                                                # WC
    (spc == 17 || (30 <= spc <= 33)) && return min(q, 25f0)                                         # EC,WO,SW,BR,CK
    return q
end

# Per-species crown-ratio caps (ls/dgf.f:245-250).
@inline function _ls_cr_cap(spc::Int, cr::Float32)::Float32
    spc == 17 && return min(cr, 60f0)                    # EC
    (spc == 60 || spc == 65) && return min(cr, 85f0)     # SY, BL
    return cr
end

"""
    dgf!(s, ::LakeStates)

LS variant `dgf!` hook — fill `scratch.wk[2, i]` with ln(inside-bark DDS) for every live tree
(the ls/dgf.f regression + its OB→IB bark conversion). COR enters through `conspp` (DGCON+COR),
exactly like SN/CS. Consumed by the shared DGDRIV calibrate / growth pass.
"""
function dgf!(s::StandState, ::LakeStates)
    p, t, c, sd = s.plot, s.trees, s.calib, s.coef.species
    wk2 = view(s.scratch.wk, 2, :)
    interc = sd[:dg_interc]; vdbh = sd[:dg_vdbh]; dbhc = sd[:dg_dbh]; dbh2 = sd[:dg_dbh2]
    rdbh = sd[:dg_rdbh]; rdbhsq = sd[:dg_rdbhsq]; crwn = sd[:dg_crwn]; crsq = sd[:dg_crsq]
    sbac = sd[:dg_sba]; balc = sd[:dg_bal]; sitec = sd[:dg_site]
    ba_a = c.bark_a; ba_b = c.bark_b

    # BAGE5 (stand BA in trees ≥5") and QMDGE5 (their QMD) — ls/dgf.f:204-222.
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

    # Per-species effective QMDGE5. FVS caps the stand-wide QMDGE5 IN PLACE as it walks species in
    # INDEX order (ls/dgf.f:362-390, DO 200 ISPC=1,MAXSP): each PRESENT cap-species mutates the shared
    # QMDGE5 (a running min), so a species sees the value capped by ALL LOWER-INDEXED present cap-species,
    # not just its own cap. jl previously applied a LOCAL per-tree cap (_ls_qmd_cap(sp, qmdge5)), so an
    # UNCAPPED species (e.g. white pine sp5) missed the cap that a lower-indexed cap-species (jack pine
    # sp1, cap 13) had already applied to the shared value — biasing its RDBH/RDBHSQ terms (live-verified:
    # sp5 QMD5 13.0 not jl's uncapped 14.51). Replicate the cumulative walk here.
    nsp = nspecies(s.variant)
    present = falses(nsp)
    @inbounds for i in 1:t.n
        (t.dbh[i] > 0f0 && t.tpa[i] > 0f0) && (present[Int(t.species[i])] = true)
    end
    qmd_eff = fill(qmdge5, nsp)
    qcur = qmdge5
    @inbounds for spc in 1:nsp
        present[spc] && (qcur = _ls_qmd_cap(spc, qcur))
        qmd_eff[spc] = qcur
    end

    @inbounds for i in 1:t.n
        d = t.dbh[i]
        d <= 0f0 && continue
        sp = Int(t.species[i])
        conspp = c.dg_const[sp] + c.dg_cor[sp]                 # DGCON + COR
        cr = Float32(t.crown_pct[i]);  cr <= 0f0 && (cr = 10f0)
        cr = _ls_cr_cap(sp, cr)
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
        # OB→IB bark conversion (ls/dgf.f:288-291).
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
    ls_dgcons!(s)

LS per-stand DG setup (ls/dgf.f ENTRY DGCONS): DGCON=0 (the intercept lives in INTERC; site
dependence is the SITEC·SITEAR term), ATTEN ← OBSERV (per-species observation count = the DGDRIV
calibration prior weight), plus the per-stand bark copy (`calib.bark_a/b` ← intercept 0 / slope
BKRAT) that CFTOPK + DGDRIV read. Identical to `cs_dgcons!`.
"""
function ls_dgcons!(s::StandState)
    c = s.calib; sd = s.coef.species; ctl = s.control
    ba = sd[:bark_intercept]; bb = sd[:bark_slope]; obs = sd[:dg_observ]
    @inbounds for sp in 1:MAXSP
        c.bark_a[sp] = ba[sp]; c.bark_b[sp] = bb[sp]
        c.dg_const[sp] = 0f0
        c.atten[sp] = obs[sp]              # ATTEN = OBSERV(ISPC)
        ctl.dg_cor2_on && ctl.dg_cor2[sp] > 0f0 && (c.dg_const[sp] += log(ctl.dg_cor2[sp]))
    end
    return s
end

"""
    _ls_init_crowns!(s)

LS CRATET: dub the INITIAL crown ratio for inventory trees with no input crown (`crown_pct==0`)
using the eastern (TWIGS) crown model, so the calibration `dgf!` sees a real crown (its `CRWNC·CR +
CRSQC·CR²` terms). FVS runs this after DGDRIV's DENSE backdate, so the crown model reads the
BACKDATED-dbh per-acre BA. Identical to `_cs_init_crowns!`.
"""
function _ls_init_crowns!(s::StandState)
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
