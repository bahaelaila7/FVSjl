# =============================================================================
# crown.jl (inlandempire) — IE per-tree CCF (ie/ccfcal.f MODE=1) + crown width (MODE=2).
# stand CCF = Σ CCFT is the RELDEN the DG (dgf!) and height (htgf) chunks read (western polynomial,
# NOT the eastern crown-width→area path). Per-species dispatch (3 classes) + a large-tree/small-tree split.
#   CCFT poly (D large):   RD1(sp) + D*RD2(sp) + D^2*RD3(sp)
#   CCFT power (D small):  RDA(sp) * D^RDB(sp)
# =============================================================================

const IE_RD1 = Float32[.03,.02,.11,.04,.03,.03,.01925,.03,.03,.03,.03, .02,.01925,.03, .01925,.01925,.01925, .03,.03,.03,.03,.03,.03]
const IE_RD2 = Float32[.0167,.0148,.0333,.0270,.0215,.0238,.01676,.0173,.0216,.0180,.0215, .0148,.01676,.0216, .01676,.01676,.01676, .0238,.0215,.0238,.0238,.0215,.0215]
const IE_RD3 = Float32[.00230,.00338,.00259,.00405,.00363,.00490,.00365,.00259,.00405,.00281,.00363, .00338,.00365,.00405, .00365,.00365,.00365, .00490,.00363,.00490,.00490,.00363,.00363]
const IE_RDA = Float32[0.009884,0.007244,0.017299,0.015248,0.011109,0.008915,0.009187,0.007875,0.011402,0.007813,0.011109, .007244,.009187,.011402, .009187,.009187,.009187, .008915,.011109,.008915,.008915,.011109,.011109]
const IE_RDB = Float32[1.6667,1.8182,1.5571,1.7333,1.7250,1.7800,1.7600,1.7360,1.7560,1.7680,1.7250, 1.8182,1.7600,1.7560, 1.7600,1.7600,1.7600, 1.7800,1.7250,1.7800,1.7800,1.7250,1.7250]

# crown-width B1..B6 (ie/ccfcal.f MODE=2, exp form for sp 1:10,12:22; sp 11/23 special).
const IE_CWB1 = Float32[1.04050,1.02478,1.01685,1.03030,1.02460,1.03597,1.03992,1.02687,1.02886,1.02687,0.0, 1.02478,1.03992,1.02886, 1.03992,1.03992,1.03992, 1.03597,1.02460, 1.03597,1.03597, 1.02460,0.0]
const IE_CWB2 = Float32[1.27990,0.99889,1.48372,1.14079,1.35223,1.46111,1.58777,1.28027,1.01255,1.49085,0.0, 0.99889,1.58777,1.01255, 1.58777,1.58777,1.58777, 1.46111,1.35223, 1.46111,1.46111, 1.35223,0.0]
const IE_CWB3 = Float32[0.11941,0.19422,0.27378,0.20904,0.24844,0.26289,0.30812,0.22490,0.30374,0.18620,0.0, 0.19422,0.30812,0.30374, 0.30812,0.30812,0.30812, 0.26289,0.24844, 0.26289,0.26289, 0.24844,0.0]
const IE_CWB4 = Float32[0.42745,0.59423,0.49646,0.38787,0.41212,0.18779,0.64934,0.47075,0.37093,0.68272,0.0, 0.59423,0.64934,0.37093, 0.64934,0.64934,0.64934, 0.18779,0.41212, 0.18779,0.18779, 0.41212,0.0]
const IE_CWB5 = Float32[0.0,-0.09078,-0.18669,0.0,-0.10436,0.0,-0.38964,-0.15911,-0.13731,-0.28242,0.0, -0.09078,-0.38964,-0.13731, -0.38964,-0.38964,-0.38964, 0.0,-0.10436, 0.0,0.0, -0.10436,0.0]
const IE_CWB6 = Float32[-0.07182,-0.02341,-0.01509,0.0,0.03539,0.0,0.0,0.0,0.0,0.0,0.0, -0.02341,0.0,0.0, 0.0,0.0,0.0, 0.03539,0.0, 0.0,0.0, 0.03539,0.0]

"""Per-tree CCF contribution (ie/ccfcal.f MODE=1), excluding the ×TPA factor."""
@inline function ie_tree_ccf(sp::Integer, d::Real)::Float32
    d <= 0f0 && return 0f0
    D = Float32(d)
    if sp <= 12 || sp == 14 || sp == 23
        return D >= 10f0 ? IE_RD1[sp] + D*IE_RD2[sp] + D*D*IE_RD3[sp] : IE_RDA[sp] * D^IE_RDB[sp]
    elseif sp == 19 || sp == 22
        D >= 10f0 && return IE_RD1[sp] + D*IE_RD2[sp] + D*D*IE_RD3[sp]
        return D > 0.1f0 ? IE_RDA[sp] * D^IE_RDB[sp] : 0.001f0
    else  # sp 13,15,16,17,18,20,21
        D >= 1f0 && return IE_RD1[sp] + D*IE_RD2[sp] + D*D*IE_RD3[sp]
        return D > 0.1f0 ? IE_RDA[sp] * D^IE_RDB[sp] : 0.001f0
    end
end

"""IE crown width (ie/ccfcal.f MODE=2). JCR = crown ratio (%); barea = BA (or OLDBA under thin/fire).
IFOR==5 (Colville) uses R6CRWD (deferred). Returns width clamped [0.1, 99.9]."""
@inline function ie_crown_width(sp::Integer, d::Real, h::Real, jcr::Integer, barea::Real)::Float32
    (jcr <= 0 || h <= 0 || d <= 0 || barea <= 0) && return 0.1f0
    D = Float32(d); H = Float32(h)
    cl = Float32(jcr) * H * 0.01f0
    cl <= 0f0 && return 0.1f0
    cw = if sp == 11 || sp == 23
        if H <= 5f0
            0.8f0 * H * max(0.5f0, jcr * 0.01f0)
        elseif H >= 15f0
            6.90396f0 * D^0.55645f0 * H^(-0.28509f0) * cl^0.20430f0
        else
            c1 = 0.8f0 * H * max(0.5f0, jcr * 0.01f0)
            c2 = 6.90396f0 * D^0.55645f0 * H^(-0.28509f0) * cl^0.20430f0
            w = (H - 5f0) * 0.1f0
            c1 * (1f0 - w) + c2 * w
        end
    else
        ba = barea < 1f0 ? 1f0 : Float32(barea)
        IE_CWB1[sp] * exp(IE_CWB2[sp] + IE_CWB3[sp]*log(cl) + IE_CWB4[sp]*log(D) +
                          IE_CWB5[sp]*log(H) + IE_CWB6[sp]*log(ba))
    end
    cw > 99.9f0 && (cw = 99.9f0)
    cw < 0.1f0 && (cw = 0.1f0)
    return cw
end

# --- crown ratio (ie/crown.f) — chunk 5b. Tables in crown_coefficients.jl (dumped from live). -------------
"""IE per-species crown constant: CRCON(sp) = CRHAB(MAPHAB(ITYPE,sp), sp) (ie/crown.f:757)."""
@inline function ie_crcon(itype::Int, sp::Int)::Float32
    IE_CRHAB[clamp(Int(IE_CRMAPHAB[itype, sp]), 1, 14), sp]
end

"""
    crown_ratio_update!(s, ::InlandEmpire; fint, lstart, crown_sdi, ...)

IE crown (ie/crown.f). Per-species variant branch: NIVAR (sp≤12,14,23) logistic PCR/DCR change; CRVAR (19,22)
+ LPIJU (15,16) linear crown-length; UTTVAR (13,17,18,20,21) Weibull. Bounds NIVAR [5,95], others [10,95].
Cycle update (D-backdated≥3 for NIVAR). LSTART dub via DUBSCR/BACHLO (RNG-cornered) for missing crowns.
"""
function crown_ratio_update!(s::StandState, ::InlandEmpire; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0.0f0, kwargs...)
    p, t, c = s.plot, s.trees, s.calib
    t.n == 0 && return s
    itype = Int(p.habitat_input); (itype < 1 || itype > 30) && (itype = 1)
    ba = p.basal_area; relden = p.relative_density
    lnba = ba > 0f0 ? log(ba) : 0f0; lnrd = relden > 0f0 ? log(relden) : 0f0
    reldm1 = p.relative_density_prev; oba = p.old_ba; rdm1 = reldm1
    if reldm1 < 100f0; oba = ba; rdm1 = relden; end
    x1 = (!lstart && oba > 0f0) ? log(oba) : 0f0
    x2 = (!lstart && rdm1 > 0f0) ? log(rdm1) : 0f0
    dgsd = s.control.dg_sd
    ba_a = c.bark_a; ba_b = c.bark_b
    # ISORT: descending-DBH rank (ie/crown.f:152, for UTTVAR Weibull X). IND is DBH-descending order.
    nlim = t.n + (lstart ? Int(t.ndead) : 0)
    @inbounds for i in 1:nlim
        t.tpa[i] <= 0f0 && continue
        icr = Int(t.crown_pct[i])
        (lstart && icr > 0) && continue
        icr < 0 && (t.crown_pct[i] = Int32(-icr); continue)
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        d <= 0f0 && continue
        bark = bark_ratio(ba_a, ba_b, sp, d)
        crcon = ie_crcon(itype, sp)
        nivar = sp <= 12 || sp == 14 || sp == 23
        crvar = sp == 19 || sp == 22
        lpiju = sp == 15 || sp == 16
        uttvar = !(nivar || crvar || lpiju)
        local icri::Int
        b7=IE_CRPARM[sp,7]; b8=IE_CRPARM[sp,8]; b9=IE_CRPARM[sp,9]; b10=IE_CRPARM[sp,10]
        b11=IE_CRPARM[sp,11]; b12=IE_CRPARM[sp,12]; b13=IE_CRPARM[sp,13]; b14=IE_CRPARM[sp,14]
        if crvar || lpiju
            hf = h + t.ht_growth[i]
            cl = crvar ? (5.17281f0 + 0.32552f0*hf - 0.01675f0*ba) : (-0.59373f0 + 0.67703f0*hf)
            cl < 1f0 && (cl = 1f0); cl > hf && (cl = hf)
            crnew = (cl / hf) * 100f0
            icri = _ie_crown_label53(crnew, icr, lstart, fint, sp, d)
        elseif nivar && (!lstart && (d - t.diam_growth[i]/bark) < 3f0)
            continue                                    # cycling: backdated D<3 keeps its crown (GOTO 60)
        elseif nivar
            xcrcon = crcon + IE_CRPARM[sp,1]*ba + IE_CRPARM[sp,2]*ba*ba + IE_CRPARM[sp,3]*lnba +
                     IE_CRPARM[sp,4]*relden + IE_CRPARM[sp,5]*relden*relden + IE_CRPARM[sp,6]*lnrd
            pp = t.crown_ratio[i]; pp < 0.01f0 && (pp = 0.01f0)      # P = PCT
            pcr = xcrcon + b7*d + b8*d*d + b9*log(d) + b10*h + b11*h*h + b12*log(h) + b13*pp + b14*log(pp)
            exppcr = exp(pcr); expdcr = 0f0
            if !lstart
                dcrcon = crcon + IE_CRPARM[sp,1]*oba + IE_CRPARM[sp,2]*oba*oba + IE_CRPARM[sp,3]*x1 +
                         IE_CRPARM[sp,4]*rdm1 + IE_CRPARM[sp,5]*rdm1*rdm1 + IE_CRPARM[sp,6]*x2
                db = d - t.diam_growth[i]/bark; db <= 0f0 && (db = d)
                hb = h - t.ht_growth[i]; hb <= 0f0 && (hb = h)
                # OLDPCT: previous cycle's PCT (gradd.f:267 snapshot). crown.f:480 falls back to current PCT
                # when OLDPCT<=0 (or OLDPCT>PCT under thinning — the thin branch is omitted; no removals here).
                pb = t.old_crown_pct[i]; pb <= 0f0 && (pb = t.crown_ratio[i])
                pb < 0.01f0 && (pb = 0.01f0)
                dcr = dcrcon + b7*db + b8*db*db + b9*log(db) + b10*hb + b11*hb*hb + b12*log(hb) + b13*pb + b14*log(pb)
                expdcr = exp(dcr)
            end
            chg = exppcr - expdcr
            if !lstart || icr > 0
                pdifpy = chg / Float32(icr) / fint * 100f0
                pdifpy > 0.01f0  && (chg = Float32(icr) * 0.01f0 * fint / 100f0)
                pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint / 100f0)
            end
            icri = trunc(Int, Float32(icr) + chg*100f0 + 0.50005f0)   # DLOW=0/DHI=99/CRNMLT=1 defaults
            (lstart && dgsd >= 1f0) && (icri = trunc(Int, bachlo(s.rng, Float32(icri), IE_CRSD)))
            # CRMAX cap (crown.f:556-568), skipped at LSTART or ICR==0
            if !(lstart || icr == 0)
                crln = h * Float32(icr) / 100f0
                crmax = (crln + t.ht_growth[i]) / (h + t.ht_growth[i]) * 100f0
                icri < 10 && (icri = trunc(Int, crmax + 0.5f0))
                Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
            end
        else
            # UTTVAR Weibull (sp13,17,18,20,21): needs ISORT/RANN/DUBSCR — deferred (not in iet01).
            lstart || continue
            continue
        end
        # final bounds (crown.f:382-390)
        icri > 95 && (icri = 95)
        if nivar
            icri < 5 && (icri = 5)                       # CRNMLT==1 default
        else
            icri < 10 && (icri = 10); icri < 1 && (icri = 1)
        end
        t.crown_pct[i] = Int32(icri)
    end
    return s
end

# ie/crown.f label 53: CRVAR/LPIJU change bound + ICRI (CHG=CRNEW-ICR, PDIFPY ±1%/yr, CRMAX cap).
@inline function _ie_crown_label53(crnew::Float32, icr::Int, lstart::Bool, fint::Float32, sp::Int, d::Float32)
    chg = crnew - Float32(icr)
    if !lstart || icr > 0
        pdifpy = chg / Float32(icr) / fint
        pdifpy > 0.01f0  && (chg = Float32(icr) * 0.01f0 * fint)
        pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
        crnew = Float32(icr) + chg                       # CRNMLT=1 default
    end
    return trunc(Int, crnew + 0.5f0)
end
