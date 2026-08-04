# =============================================================================
# crown.jl (britishcolumbia) — BC per-tree CCF + stand RELDEN (canada/bc/ccfcal.f). CHUNK 5 (CCF part).
#
# Standard western CCF (same form as KT): D≥10" → RD1(IC)+D·RD2(IC)+D²·RD3(IC); D<10" → RDA(IC)·D^RDB(IC);
# clamp ≥0.001; ×P (TPA). NI-derived RD arrays (11 IC classes) for sp1-10,14 (IC=sp for 1-10, IC=3=FD for 14).
# ⚠ PN-derived (sp11-13,15) uses pnRD/INDCCF — TODO (not exercised by all_BC = sp14). Crown-width (B1-6) TODO.
# This CCF is the density SPINE: feeds RELDEN (dgf DGCCF term) + regent per-subcycle density (RDJ).
# =============================================================================

# NI-derived CCF coefficients, 11 IC classes (ccfcal.f:87-100).
const BC_CCF_RD1 = Float32[.03,.02,.11,.04,.03,.03,.01925,.03,.03,.03,.03]
const BC_CCF_RD2 = Float32[.0167,.0148,.0333,.0270,.0215,.0238,.01676,.0173,.0216,.0180,.0215]
const BC_CCF_RD3 = Float32[.00230,.00338,.00259,.00405,.00363,.00490,.00365,.00259,.00405,.00281,.00363]
const BC_CCF_RDA = Float32[0.009884,0.007244,0.017299,0.015248,0.011109,0.008915,0.009187,0.007875,0.011402,0.007813,0.011109]
const BC_CCF_RDB = Float32[1.6667,1.8182,1.5571,1.7333,1.7250,1.7800,1.7600,1.7360,1.7560,1.7680,1.7250]

"""IC species→coefficient-class map (ccfcal.f:162-171): sp1-10 → sp; sp14(OC) → 3(FD). ⚠ sp11-13,15 (PN) TODO."""
@inline function _bc_ccf_ic(sp::Integer)
    (1 <= sp <= 10) ? sp : (sp == 14 ? 3 : 0)
end

"""
    bc_tree_ccf(sp, d, p) -> per-tree CCF contribution (ccfcal.f MODE 1, NI species)

D in inches, P = trees/acre. Returns 0 for PN species (11-13,15) until pnRD ported.
"""
@inline function bc_tree_ccf(sp::Integer, d::Real, p::Real)
    ic = _bc_ccf_ic(sp)
    ic < 1 && return 0f0                                     # PN species (TODO)
    D = Float32(d)
    ccf = D >= 10f0 ? BC_CCF_RD1[ic] + D*BC_CCF_RD2[ic] + D*D*BC_CCF_RD3[ic] :
                      BC_CCF_RDA[ic] * (D ^ BC_CCF_RDB[ic])
    ccf < 0.001f0 && (ccf = 0.001f0)
    return ccf * Float32(p)
end

"""Stand CCF = Σ per-tree CCF over live trees (RELDEN); mirrors FVS ccfcal accumulation."""
function bc_stand_ccf(s::StandState)
    t = s.trees; acc = 0f0
    @inbounds for i in 1:t.n
        d = t.dbh[i]; (d <= 0f0 || t.tpa[i] <= 0f0) && continue
        acc += bc_tree_ccf(Int(t.species[i]), d, t.tpa[i])
    end
    return acc
end

# --- crown-RATIO update (canada/bc/crown.f, V3/CRNMD) — chunk 5b. Coeffs in crown_coefficients.jl. ---
"""
    crown_ratio_update!(s, ::BritishColumbia; fint, lstart, crown_sdi, ...)

BC V3 crown ratio (crown.f). Per-cycle CRNMD logistic PCR/DCR change: `XCRCON = CRCON(sp) +
CRLNCCF·log(RELDEN)`; `EXPPCR = CRNMD(...)`; backdated `EXPDCR` (D-DG/BARK, H-HTG, OLDPCT, OLDBA, RDM1);
`CHG = EXPPCR−EXPDCR` bounded ±1%/yr; `ICRI = ICR + CHG·100`; clamp [5,95]. Trees with backdated D<3in
keep their (regent-set) crown (GOTO 60). No CRMAX cap (BC source has it commented out). ⚠ V2 (non-ICH/IDF/
SBS/SBPS) + LSTART DUBSCR (missing-crown dubbing, RNG-cornered) deferred — all_BC carries crowns & is V3.
"""
function crown_ratio_update!(s::StandState, ::BritishColumbia; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0.0f0, kwargs...)
    p, t, c = s.plot, s.trees, s.calib
    t.n == 0 && return s
    zone, series = bc_stand_zone(s)
    bc_lv2atv(zone) && return bc_v2_crown_ratio_update!(s; fint = fint, lstart = lstart)
    crcon, crhtdbh, crht, crdbh2, crbal, crlnccf, crsd = bc_crcons!(s)
    ba = p.basal_area; relden = p.relative_density
    reldm1 = p.relative_density_prev; oba = p.old_ba; rdm1 = reldm1
    if reldm1 < 100f0; oba = ba; rdm1 = relden; end
    lnrd = log(max(0.01f0, relden))
    ba_a = c.bark_a; ba_b = c.bark_b
    @inbounds for i in 1:t.n
        t.tpa[i] <= 0f0 && continue
        icr = Int(t.crown_pct[i])
        (lstart && icr > 0) && continue
        icr < 0 && (t.crown_pct[i] = Int32(-icr); continue)
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        d <= 0f0 && continue
        d * BC_INtoCM < 2f0 && continue                        # <2cm → dub path; per-cycle keeps crown
        bark = bark_ratio(ba_a, ba_b, sp, d)
        (!lstart && (d - t.diam_growth[i]/bark) < 3f0) && continue   # backdated D<3in: keep regent crown (GOTO 60)
        pp = t.crown_ratio[i]; pp < 0.01f0 && (pp = 0.01f0)
        balf = (1f0 - pp/100f0) * ba
        xcrcon = crcon[sp] + crlnccf[sp] * lnrd
        exppcr = bc_crnmd(xcrcon, crhtdbh[sp], crht[sp], crdbh2[sp], crbal[sp], d, h, balf)
        expdcr = 0f0
        if !lstart
            db = d - t.diam_growth[i]/bark; db <= 0f0 && (db = d)
            hb = h - t.ht_growth[i]; hb <= 0f0 && (hb = h)
            pb = t.old_crown_pct[i]; pb <= 0f0 && (pb = t.crown_ratio[i])   # OLDPCT fallback (crown.f:474-476)
            pb < 0.01f0 && (pb = 0.01f0)
            balb = (1f0 - pb/100f0) * oba
            # V3 backdated CR reuses XCRCON (crown.f:486) — DCRCON is the V2-only branch; density term not re-backdated.
            expdcr = bc_crnmd(xcrcon, crhtdbh[sp], crht[sp], crdbh2[sp], crbal[sp], db, hb, balb)
        end
        chg = exppcr - expdcr
        if !lstart || icr > 0
            pdifpy = chg / Float32(icr) / fint * 100f0
            pdifpy > 0.01f0  && (chg = Float32(icr) * 0.01f0 * fint / 100f0)
            pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint / 100f0)
        end
        icri = trunc(Int, Float32(icr) + chg*100f0 + 0.50005f0)     # CRNMLT=1, DLOW/DHI defaults
        icri > 95 && (icri = 95)
        icri < 5 && (icri = 5)                                       # CRNMLT==1 lower bound
        t.crown_pct[i] = Int32(icri)
    end
    return s
end

# --- V2 (LV2ATV) crown ratio (crown.f LV2ATV branches) — exp-form logistic, PARM(sp,14) + CRHAB. ---
# XCRCON = CRCON + PARM1·BA + PARM2·BA² + PARM3·ln(BA) + PARM4·RELDEN + PARM5·RELDEN² + PARM6·ln(RELDEN);
# PCR = XCRCON + PARM7·D + PARM8·D² + PARM9·ln(D) + PARM10·H + PARM11·H² + PARM12·ln(H) + PARM13·P + PARM14·ln(P);
# EXPPCR = exp(PCR); backdated DCR uses OBA/RDM1/X1/X2 + backdated D/H/OLDPCT. CRCON = CRHAB(MAPHAB(ITYPE,sp),sp)
# (habitat intercept; ITYPE=4 grinit default). Small-tree bypass D<3in. Post-EXPPCR (chg bound ±1%/yr, ICRI
# [5,95]) is shared with V3. Coeffs = crown.f DATA PARM(MAXSP,14) + CRHAB(14,MAXSP), extracted col-major.
const BC_CR_PARM = ([  # [sp][1..14]
    Float32[0.0,0.0,-0.34566,0.0,0.0,0.0,0.03882,-0.0007,0.0,0.0,0.0,-0.21217,0.00301,0.0],
    Float32[-0.00204,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.30066,0.0,0.0,-0.59302,0.0,0.19558],
    Float32[0.0,0.0,0.0,0.0,0.0,-0.15334,0.0,0.0,0.3384,0.0,0.0,-0.59685,0.0,0.16488],
    Float32[-0.00183,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.24293,0.0,0.0,-0.25601,0.0,0.0726],
    Float32[0.0,-1.902f-6,0.0,0.0,0.0,0.0,0.03027,-0.00055,0.0,0.0,0.0,-0.25776,0.0,0.06887],
    Float32[0.0,0.0,0.17479,-0.00183,0.0,0.0,-0.0056,0.0,0.0,0.0,0.0,0.0,0.0,0.1105],
    Float32[0.0,0.0,0.0,0.0,0.0,-0.18555,0.0,0.0,0.53172,-0.02989,0.00011,0.0,0.0042,0.0],
    Float32[-0.00203,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.29699,0.0,0.0,-0.38334,0.0,0.09918],
    Float32[-0.0019,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.23372,0.0,0.0,-0.28433,0.001903,0.0],
    Float32[-0.002165,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.26558,0.0,0.0,-0.31555,0.0,0.16072],
    Float32[-0.00264,0.0,0.0,0.0,5.116f-6,0.0,0.0,0.0,0.0,0.0,0.0,-0.25138,0.0,0.0514],
    Float32[-0.00264,0.0,0.0,0.0,5.116f-6,0.0,0.0,0.0,0.0,0.0,0.0,-0.25138,0.0,0.0514],
    Float32[-0.00264,0.0,0.0,0.0,5.116f-6,0.0,0.0,0.0,0.0,0.0,0.0,-0.25138,0.0,0.0514],
    Float32[0.0,0.0,0.0,0.0,0.0,-0.15334,0.0,0.0,0.3384,0.0,0.0,-0.59685,0.0,0.16488],
    Float32[-0.00264,0.0,0.0,0.0,5.116f-6,0.0,0.0,0.0,0.0,0.0,0.0,-0.25138,0.0,0.0514],
]...,)
# crown.f has its OWN MAPHAB(30,MAXSP) (DISTINCT from dgf.f's) — ICRHAB=MAPHAB(ITYPE,sp). ITYPE=4 (grinit
# default, confirmed via CRHAB14 instrument: sp14→ICRHAB=2→CRHAB(2,14)=0.7271). Column ITYPE=4, sp1..15.
# ⚠ Full crown MAPHAB[30,15] deferred until a non-default-ITYPE V2 stand appears.
const BC_CR_ICRHAB = Int[2,2,4,2,1,1,2,2,2,1,1,1,1,2,1]   # crown MAPHAB(ITYPE=4, sp)
const BC_CR_CRHAB = ([  # [sp][1..14 habitat]
    Float32[0.8884,0.7309,0.9347,0.9888,0.9945,1.1126,1.0263,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
    Float32[0.06533,0.03441,0.2307,0.1661,-0.1253,-0.05018,0.11005,0.08113,0.1782,0.03919,0.2107,0.0,0.0,0.0],
    Float32[0.8643,0.7271,0.984,0.8127,0.8874,0.7055,0.7708,0.7849,0.8038,0.8742,0.8232,0.8415,0.9759,0.0],
    Float32[-0.2304,-0.5421,-0.4343,-0.3759,-0.4129,-0.4879,-0.2674,-0.1941,0.0,0.0,0.0,0.0,0.0,0.0],
    Float32[-0.2413,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
    Float32[-1.6053,-1.7128,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
    Float32[-0.3785,-0.4142,-0.3985,-0.2987,-0.381,-0.4087,-0.3577,-0.2994,-0.2486,-0.2863,-0.1968,-0.4931,-0.2676,-0.5625],
    Float32[0.05351,-0.05031,0.1075,-0.1872,0.01729,0.03667,0.01885,0.09102,0.1371,0.08368,0.123,-0.02365,0.0,0.0],
    Float32[0.09453,-0.0774,0.07113,0.2039,0.06176,0.1513,0.09086,0.158,0.09229,0.01551,0.0,0.0,0.0,0.0],
    Float32[-0.9436,-0.8654,-0.8849,-0.9067,-0.8783,-1.0103,-1.0268,-1.005,-1.0301,0.0,0.0,0.0,0.0,0.0],
    Float32[0.4649,0.3211,0.197,0.2295,0.3383,0.345,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
    Float32[0.4649,0.3211,0.197,0.2295,0.3383,0.345,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
    Float32[0.4649,0.3211,0.197,0.2295,0.3383,0.345,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
    Float32[0.8643,0.7271,0.984,0.8127,0.8874,0.7055,0.7708,0.7849,0.8038,0.8742,0.8232,0.8415,0.9759,0.0],
    Float32[0.4649,0.3211,0.197,0.2295,0.3383,0.345,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
]...,)

"""BC V2 crown ratio (crown.f LV2ATV) — PARM exp-form PCR/DCR. CRCON=CRHAB(MAPHAB(4,sp),sp); D<3in bypass."""
function bc_v2_crown_ratio_update!(s::StandState; fint::Float32 = 10.0f0, lstart::Bool = false)
    p, t = s.plot, s.trees
    t.n == 0 && return s
    ba = p.basal_area; relden = p.relative_density
    reldm1 = p.relative_density_prev; oba = p.old_ba; rdm1 = reldm1
    if reldm1 < 100f0; oba = ba; rdm1 = relden; end          # crown.f:344-348 (no density impact if RDM1<100)
    lnba = log(ba); lnrd = log(relden)
    x1 = oba > 0f0 ? log(oba) : 0f0; x2 = rdm1 > 0f0 ? log(rdm1) : 0f0
    ba_a = s.calib.bark_a; ba_b = s.calib.bark_b
    it = BC_V2_ITYPE
    @inbounds for i in 1:t.n
        t.tpa[i] <= 0f0 && continue
        icr = Int(t.crown_pct[i])
        (lstart && icr > 0) && continue
        icr < 0 && (t.crown_pct[i] = Int32(-icr); continue)
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        d <= 0f0 && continue
        d < 3f0 && continue                                   # V2 small-tree bypass (crown.f:432)
        bark = bark_ratio(ba_a, ba_b, sp, d)
        (!lstart && (d - t.diam_growth[i]/bark) < 3f0) && continue   # backdated D<3in → keep (crown.f:440)
        prm = BC_CR_PARM[sp]
        crcon = BC_CR_CRHAB[sp][BC_CR_ICRHAB[sp]]            # crown's own MAPHAB (ITYPE=4), NOT dgf's
        xcrcon = crcon + prm[1]*ba + prm[2]*ba*ba + prm[3]*lnba +
                 prm[4]*relden + prm[5]*relden*relden + prm[6]*lnrd
        pp = t.crown_ratio[i]; pp < 0.01f0 && (pp = 0.01f0)
        pcr = xcrcon + prm[7]*d + prm[8]*d*d + prm[9]*log(d) + prm[10]*h + prm[11]*h*h +
              prm[12]*log(h) + prm[13]*pp + prm[14]*log(pp)
        exppcr = exp(pcr)
        expdcr = 0f0
        if !lstart
            dcrcon = crcon + prm[1]*oba + prm[2]*oba*oba + prm[3]*x1 +
                     prm[4]*rdm1 + prm[5]*rdm1*rdm1 + prm[6]*x2
            db = d - t.diam_growth[i]/bark; db <= 0f0 && (db = d)
            hb = h - t.ht_growth[i]; hb <= 0f0 && (hb = h)
            pb = t.old_crown_pct[i]; pb <= 0f0 && (pb = t.crown_ratio[i]); pb < 0.01f0 && (pb = 0.01f0)
            dcr = dcrcon + prm[7]*db + prm[8]*db*db + prm[9]*log(db) + prm[10]*hb + prm[11]*hb*hb +
                  prm[12]*log(hb) + prm[13]*pb + prm[14]*log(pb)
            expdcr = exp(dcr)
        end
        chg = exppcr - expdcr
        if !lstart || icr > 0
            pdifpy = chg / Float32(icr) / fint * 100f0
            pdifpy > 0.01f0  && (chg = Float32(icr) * 0.01f0 * fint / 100f0)
            pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint / 100f0)
        end
        icri = trunc(Int, Float32(icr) + chg*100f0 + 0.50005f0)
        icri > 95 && (icri = 95)
        icri < 5 && (icri = 5)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end
