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
        return D >= 10f0 ? IE_RD1[sp] + D*IE_RD2[sp] + D*D*IE_RD3[sp] : IE_RDA[sp] * fpow(D, IE_RDB[sp])
    elseif sp == 19 || sp == 22
        D >= 10f0 && return IE_RD1[sp] + D*IE_RD2[sp] + D*D*IE_RD3[sp]
        return D > 0.1f0 ? IE_RDA[sp] * fpow(D, IE_RDB[sp]) : 0.001f0
    else  # sp 13,15,16,17,18,20,21
        D >= 1f0 && return IE_RD1[sp] + D*IE_RD2[sp] + D*D*IE_RD3[sp]
        return D > 0.1f0 ? IE_RDA[sp] * fpow(D, IE_RDB[sp]) : 0.001f0
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
            6.90396f0 * fpow(D, 0.55645f0) * fpow(H, -0.28509f0) * fpow(cl, 0.20430f0)
        else
            c1 = 0.8f0 * H * max(0.5f0, jcr * 0.01f0)
            c2 = 6.90396f0 * fpow(D, 0.55645f0) * fpow(H, -0.28509f0) * fpow(cl, 0.20430f0)
            w = (H - 5f0) * 0.1f0
            c1 * (1f0 - w) + c2 * w
        end
    else
        ba = barea < 1f0 ? 1f0 : Float32(barea)
        IE_CWB1[sp] * fexp(IE_CWB2[sp] + IE_CWB3[sp]*flog(cl) + IE_CWB4[sp]*flog(D) +
                          IE_CWB5[sp]*flog(H) + IE_CWB6[sp]*flog(ba))
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

# ie/crown.f:177-180 RMAIAS/RMAILM — stand MAI constants for the UTTVAR (aspen-group) DUBSCR TMAI term.
# RMAIAS = min(ADJMAI(746, SITEAR(18), 10), 128); RMAILM = ADJMAI(101, SITEAR(13), 10). (crown.f's RMAILM>128
# guard is the well-known `RMAI=128` typo — it never caps RMAILM — so RMAILM is left uncapped, faithfully.)
_ie_rmai_as(s::StandState)::Float32 = (r = _adjmai(746, s.plot.sp_site_index[18], 10f0); r > 128f0 ? 128f0 : r)
_ie_rmai_lm(s::StandState)::Float32 = _adjmai(101, s.plot.sp_site_index[13], 10f0)

"""
    crown_ratio_update!(s, ::InlandEmpire; fint, lstart, crown_sdi, ...)

IE crown (ie/crown.f). Per-species variant branch: NIVAR (sp≤12,14,23) logistic PCR/DCR change; CRVAR (19,22)
+ LPIJU (15,16) linear crown-length; UTTVAR (13,17,18,20,21) Weibull. Bounds NIVAR [5,95], others [10,95].
Cycle update (D-backdated≥3 for NIVAR). LSTART dub via DUBSCR/BACHLO (RNG-cornered) for missing crowns.
"""
function crown_ratio_update!(s::StandState, ::InlandEmpire; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0.0f0, kwargs...)
    p, t, c = s.plot, s.trees, s.calib
    if t.n == 0
        # crown.f:191 IF((ITRN.LE.0).AND.(IREC2.LT.MAXTP1)) GO TO 74 — with NO live trees FVS still dubs the missing
        # crowns of the inventory's standing-dead records (each NIVAR/UTTVAR one main-stream DUBSCR draw). Returning
        # here skipped those draws: an all-dead FIA stand (1627671438290487) entered cycle 1 with jl's main RANN 69
        # draws behind live, desyncing every birth-cycle REGENT ZZRAN.
        lstart && t.ndead > 0 && _ie_dead_crown_dub!(s, p.basal_area, s.control.dg_sd)
        return s
    end
    itype = Int(p.habitat_input); (itype < 1 || itype > 30) && (itype = 1)
    ba = p.basal_area; relden = p.relative_density
    lnba = ba > 0f0 ? flog(ba) : 0f0; lnrd = relden > 0f0 ? flog(relden) : 0f0
    reldm1 = p.relative_density_prev; oba = p.old_ba; rdm1 = reldm1
    if reldm1 < 100f0; oba = ba; rdm1 = relden; end
    x1 = (!lstart && oba > 0f0) ? flog(oba) : 0f0
    x2 = (!lstart && rdm1 > 0f0) ? flog(rdm1) : 0f0
    dgsd = s.control.dg_sd
    ba_a = c.bark_a; ba_b = c.bark_b
    # ISORT: descending-DBH rank (ie/crown.f:152, for UTTVAR Weibull X). IND is DBH-descending order.
    # #158-class species-major RNG order: FVS ie/crown.f processes the LIVE trees SPECIES-MAJOR
    # (`DO 70 ISPC=1,MAXSP; DO 60 I3=I1,I2; I=IND1(I3)`), so the per-tree NIVAR/DUBSCR BACHLO crown draw
    # (line ~143, rejection-sampled with a species-specific SD) is consumed in species order. jl dubbed in raw
    # tree-index order ⇒ on a multi-species seedling cohort the per-tree crown draws were mis-assigned ⇒ wrong
    # crown ⇒ small-tree height over-/under-grows ⇒ compounding dense-stand divergence. Only the lstart DUB
    # path draws RNG, so iterate species-major there; cycling (no draw) keeps natural order. Within a species,
    # tree-index order = FVS IND1 (stable species bucket). The deterministic per-tree crown update is
    # order-independent, so this is a no-op except on the RNG stream. IE-only (dispatches on ::InlandEmpire).
    # NOTE: the LIVE loop is live-only (1:t.n); the cycle-0 DEAD records are dubbed in a SEPARATE reverse-index
    # pass below (crown.f:634 DO 79), NOT folded into this species-major order — see the dead-dub loop.
    # ISORT for the UTTVAR (aspen-group) Weibull X: whole-stand DBH rank (ie/crown.f:261-263 ISORT(IND)=ITRN-JJ+1,
    # so largest DBH → ITRN, smallest → 1). Ranks on t.dbh directly — post-growth at cycling (simulate.jl applies DG
    # before CROWN), inventory DBH at LSTART — matching FVS's DBH(I). Once/cycle; local buffers, not the hot path.
    isort = crown_isort(s; lstart = lstart)   # shared (crown_init.jl): adds cratet.f's IND1-seeded LSTART sort
    # RMAIAS/RMAILM (ie/crown.f:177-180) — stand MAI constants; only the UTTVAR D<1 LSTART DUBSCR uses them.
    rmai_as = lstart ? _ie_rmai_as(s) : 0f0
    rmai_lm = lstart ? _ie_rmai_lm(s) : 0f0
    order = species_major_order(s)   # ie/crown.f DO 70 ISPC … IND1 — every call (RANN d≤0 draws; post-TRIPLE lineage order)
    @inbounds for i in order
        t.tpa[i] <= 0f0 && continue
        icr = Int(t.crown_pct[i])
        (lstart && icr > 0) && continue
        icr < 0 && (t.crown_pct[i] = Int32(-icr); continue)
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        nivar = sp <= 12 || sp == 14 || sp == 23
        crvar = sp == 19 || sp == 22
        lpiju = sp == 15 || sp == 16
        uttvar = !(nivar || crvar || lpiju)
        d <= 0f0 && !uttvar && continue        # non-UTTVAR need D>0; UTTVAR D≤0 cycling draws RANN (crown.f:431)
        bark = d > 0f0 ? bark_ratio(ba_a, ba_b, sp, d) : 0f0
        crcon = ie_crcon(itype, sp)
        local icri::Int
        b7=IE_CRPARM[sp,7]; b8=IE_CRPARM[sp,8]; b9=IE_CRPARM[sp,9]; b10=IE_CRPARM[sp,10]
        b11=IE_CRPARM[sp,11]; b12=IE_CRPARM[sp,12]; b13=IE_CRPARM[sp,13]; b14=IE_CRPARM[sp,14]
        if crvar || lpiju
            hf = h + t.ht_growth[i]
            cl = crvar ? (5.17281f0 + 0.32552f0*hf - 0.01675f0*ba) : (-0.59373f0 + 0.67703f0*hf)
            cl < 1f0 && (cl = 1f0); cl > hf && (cl = hf)
            crnew = (cl / hf) * 100f0
            icri = _ie_crown_label53(crnew, icr, lstart, fint, sp, d, h, t.ht_growth[i])
        elseif nivar && (!lstart && (d - t.diam_growth[i]/bark) < 3f0)
            continue                                    # cycling: backdated D<3 keeps its crown (GOTO 60)
        elseif nivar && lstart && d < 3.0f0
            # ie/crown.f:415 — NIVAR D<3 at LSTART branches to stmt 58 → CALL DUBSCR (crown.f:601), NOT the
            # PCR-change path. DUBSCR is a distinct logistic model with its OWN stochastic draw that carries a
            # |FCR|≤CRSD outer rejection loop (dubscr.f:101-104), so it consumes a DIFFERENT number of BACHLO
            # uniforms than the PCR path's single unbounded draw. Routing D<3 seedlings through the PCR path
            # (the old code) mis-consumed the main RNG stream by one bachlo per rejection ⇒ every downstream
            # draw (the regent NIVAR ZZRAN + the frozen DG serial-correlation OLDRN) desynced vs the oracle ⇒
            # dense small-tree stands scattered ±. ie_dubscr restores both the crown VALUE and the draw count.
            icri = ie_dubscr(s.rng, sp, d, h, ba, dgsd)
        elseif nivar
            xcrcon = crcon + IE_CRPARM[sp,1]*ba + IE_CRPARM[sp,2]*ba*ba + IE_CRPARM[sp,3]*lnba +
                     IE_CRPARM[sp,4]*relden + IE_CRPARM[sp,5]*relden*relden + IE_CRPARM[sp,6]*lnrd
            pp = t.crown_ratio[i]; pp < 0.01f0 && (pp = 0.01f0)      # P = PCT
            pcr = xcrcon + b7*d + b8*d*d + b9*flog(d) + b10*h + b11*h*h + b12*flog(h) + b13*pp + b14*flog(pp)
            exppcr = fexp(pcr); expdcr = 0f0
            if !lstart
                dcrcon = crcon + IE_CRPARM[sp,1]*oba + IE_CRPARM[sp,2]*oba*oba + IE_CRPARM[sp,3]*x1 +
                         IE_CRPARM[sp,4]*rdm1 + IE_CRPARM[sp,5]*rdm1*rdm1 + IE_CRPARM[sp,6]*x2
                db = d - t.diam_growth[i]/bark; db <= 0f0 && (db = d)
                hb = h - t.ht_growth[i]; hb <= 0f0 && (hb = h)
                # OLDPCT: previous cycle's PCT (gradd.f:267 snapshot). crown.f:480 falls back to current PCT
                # when OLDPCT<=0 (or OLDPCT>PCT under thinning — the thin branch is omitted; no removals here).
                pb = t.old_crown_pct[i]; pb <= 0f0 && (pb = t.crown_ratio[i])
                pb < 0.01f0 && (pb = 0.01f0)
                dcr = dcrcon + b7*db + b8*db*db + b9*flog(db) + b10*hb + b11*hb*hb + b12*flog(hb) + b13*pb + b14*flog(pb)
                expdcr = fexp(dcr)
            end
            chg = exppcr - expdcr
            if !lstart || icr > 0
                pdifpy = chg / Float32(icr) / fint * 100f0
                pdifpy > 0.01f0  && (chg = Float32(icr) * 0.01f0 * fint / 100f0)
                pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint / 100f0)
            end
            icri = trunc(Int, Float32(icr) + chg*100f0 + 0.50005f0)   # DLOW=0/DHI=99/CRNMLT=1 defaults
            (lstart && dgsd >= 1f0) && (icri = trunc(Int, bachlo(s.rng, Float32(icri), IE_CRSD)))
            # NIVAR ends with `GO TO 55` (crown.f "END OF NI BLOCK") — it SKIPS the CRMAX cap (which lives on the
            # label-53 CRVAR/LPIJU path, now in _ie_crown_label53). jl formerly applied CRMAX here (NIVAR) and
            # omitted it for CRVAR/LPIJU — inverted. No CRMAX on the NIVAR path.
        else
            # UTTVAR (sp 13,17,18,20,21 — aspen group AS/MM/PB + LM/PY): ie/crown.f:412-447 rank-Weibull crown,
            # falling through to label 53 (same CHG-bound/CRMAX as CRVAR/LPIJU). Previously a NO-OP (crown frozen at
            # the inventory value) — the documented ~7% aspen BA deficit: FVS re-predicts the crown EVERY cycle for
            # ALL DBH (UTTVAR never routes to the NIVAR D<3 freeze), converging toward the SDI-driven mean ACRNEW
            # (which recedes as RELSDI climbs), bounded ±1%/yr.
            if lstart && d < 1f0
                # crown.f:413 — D<1 at LSTART → GOTO 58 → DUBSCR (UTTVAR BCR5/6/8/9/10 with TPCCF/AVH/TMAI). Draws
                # one main-stream BACHLO (|FCR|≤CRSD rejection) when DGSD≥1 — restores the LSTART draw COUNT vs FVS.
                pt = Int(t.plot_id[i])
                tpccf = (1 <= pt <= length(s.density.point_ccf)) ? s.density.point_ccf[pt] : 0f0
                tmai = (sp == 13 || sp == 17) ? rmai_lm : rmai_as
                icri = ie_dubscr(s.rng, sp, d, h, ba, dgsd; tpccf = tpccf, avh = p.avg_height, tmai = tmai)
            else
                relsdi = p.sp_sdi_def[sp] > 0f0 ? crown_sdi / p.sp_sdi_def[sp] : 1f0   # SDIAC/SDIDEF (≤1.5)
                relsdi > 1.5f0 && (relsdi = 1.5f0)
                acrnew = IE_CRC0[sp] + IE_CRC1[sp] * relsdi * 100f0
                A = IE_WEIBA[sp]
                B = IE_WEIBB0[sp] + IE_WEIBB1[sp]*acrnew; B < 1f0 && (B = 1f0)   # crown.f:358
                C = IE_WEIBC0[sp] + IE_WEIBC1[sp]*acrnew; C < 2f0 && (C = 2f0)   # crown.f:359
                scale = 1f0 - 0.00167f0*(relden - 100f0)                          # crown.f:425-427
                scale > 1f0 && (scale = 1f0); scale < 0.30f0 && (scale = 0.30f0)
                x = d > 0f0 ? (Float32(isort[i]) / Float32(t.n)) * scale : rann!(s.rng) * scale   # crown.f:428-433
                x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
                crnew = (A + B * fpow(-flog(1f0 - x), 1f0/C)) * 10f0                   # crown.f:436,443
                icri = _ie_crown_label53(crnew, icr, lstart, fint, sp, d, h, t.ht_growth[i])
            end
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
    # ---- CYCLE-0 DEAD-TREE CROWN DUB (crown.f:633-694 `DO 79 I=IREC2,MAXTRE`) ----
    # FVS dubs MISSING crowns on the standing-dead (HISTORY 6/7/8) records in a SEPARATE pass AFTER the live
    # loop, iterating the dead partition in STORAGE-INDEX order (IREC2→MAXTRE). FVS stores dead records from
    # MAXTRE downward in inventory-READ order, so IREC2 = the LAST dead read and MAXTRE = the FIRST — i.e. the
    # DO 79 order is the REVERSE of jl's dead storage (t.n+1 = first dead read … t.n+ndead = last). Verified
    # bit-for-bit on stand 530569588126144: the oracle's DO 79 (D,H,species) sequence == jl dead indices
    # (t.n+ndead):-1:(t.n+1). Each NIVAR/UTTVAR dead tree draws one DUBSCR error (a main-stream BACHLO with the
    # |FCR|≤CRSD rejection loop) REGARDLESS OF DBH (unlike the live loop's D<3 gate); CRVAR/LPIJU use the
    # deterministic CL dub (no draw). The OLD code folded the dead into the live species-major order and routed
    # them through the D≥3 PCR path (single unbounded BACHLO) or, for UTTVAR sp17, skipped them entirely — both
    # desynced the main RNG stream by ~33 draws vs the oracle before the first REGENT ZZRAN, scattering every
    # small-tree seedling stand's ZZRAN height-growth draw (#137/#206 class). Dead crown VALUES are not reported
    # in the .sum, so the DUBSCR mean (NIVAR BCR0-3 form; UTTVAR TPCCF/AVH/TMAI terms omitted) is immaterial —
    # only the SD (=CRSD[sp], faithful) and hence the draw COUNT/ORDER matter for the RNG sync. crown.f:684-687
    # bounds the dead ICRI to [10,95] for ALL species (not the live loop's NIVAR<5). IE-only dispatch.
    lstart && _ie_dead_crown_dub!(s, ba, dgsd)
    return s
end

# crown.f:633-694 DO 79 (see the comment block above its call site): cycle-0 dub of MISSING crowns on the
# standing-dead records, in FVS's IREC2→MAXTRE order (= reverse of jl's dead storage).
function _ie_dead_crown_dub!(s::StandState, ba::Float32, dgsd::Float32)
    t = s.trees
    if t.ndead > 0
        @inbounds for i in (t.n + Int(t.ndead)):-1:(t.n + 1)
            t.tpa[i] <= 0f0 && continue
            icr = Int(t.crown_pct[i]); icr > 0 && continue
            sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
            d <= 0f0 && continue
            local icri::Int
            if sp == 19 || sp == 22                                   # CRVAR: deterministic CL, INT(CR*100) (no +0.5)
                cl = 5.17281f0 + 0.32552f0*h - 0.01675f0*ba
                cl < 1f0 && (cl = 1f0); cl > h && (cl = h)
                icri = trunc(Int, (cl / h) * 100f0)
            elseif sp == 15 || sp == 16                               # LPIJU: deterministic CL
                cl = -0.59373f0 + 0.67703f0*h
                cl < 1f0 && (cl = 1f0); cl > h && (cl = h)
                icri = trunc(Int, (cl / h) * 100f0)
            else                                                       # NIVAR + UTTVAR → DUBSCR (draws)
                icri = ie_dubscr(s.rng, sp, d, h, ba, dgsd)
            end
            icri > 95 && (icri = 95); icri < 10 && (icri = 10)         # crown.f:684-686 dead bounds (all species)
            t.crown_pct[i] = Int32(icri)
        end
    end
    return s
end

# ie/dubscr.f DUBSCR — small-tree (D<3 NIVAR) crown-ratio dub. Logistic model with a stochastic error that
# is REJECTION-BOUNDED to |FCR|≤CRSD(sp) (dubscr.f:100-104) — distinct from the D≥3 PCR path's single
# unbounded BACHLO(icri,6.35). For NIVAR the BCR5/6/8/9/10 (TPCCF/AVH/TMAI) terms are all zero, so only
# BCR0-3 enter. Returns ICRI = INT(CR*100+0.5), CR ∈ [0.05,0.95].
@inline function ie_dubscr(rng, sp::Int, d::Float32, h::Float32, ba::Float32, dgsd::Float32;
                           tpccf::Float32 = 0f0, avh::Float32 = 0f0, tmai::Float32 = 0f0)::Int
    # dubscr.f:105 — full logistic incl. the UTTVAR TPCCF/AVH/TMAI terms (BCR5/6/8/9/10 are zero for
    # NIVAR/CRVAR/LPIJU, so the NIVAR/dead-dub callers that omit these kwargs stay bit-exact).
    cr = IE_DUB_BCR0[sp] + IE_DUB_BCR1[sp]*d + IE_DUB_BCR2[sp]*h + IE_DUB_BCR3[sp]*ba +
         IE_DUB_BCR5[sp]*tpccf + (h > 0f0 ? IE_DUB_BCR6[sp]*(avh/h) : 0f0) + IE_DUB_BCR8[sp]*avh +
         IE_DUB_BCR9[sp]*(ba*tpccf) + IE_DUB_BCR10[sp]*tmai
    sd = IE_DUB_CRSD[sp]
    fcr = 0.0f0
    if dgsd >= 1.0f0
        while true
            fcr = bachlo(rng, 0.0f0, sd)
            abs(fcr) <= sd && break                     # dubscr.f:104 IF(ABS(FCR).GT.SD) GO TO 10 (redraw)
        end
    end
    abs(cr + fcr) >= 86.0f0 && (cr = 86.0f0)            # dubscr.f:105 overflow guard
    crf = 1.0f0 / (1.0f0 + fexp(cr + fcr))
    crf < 0.05f0 && (crf = 0.05f0); crf > 0.95f0 && (crf = 0.95f0)
    return trunc(Int, crf * 100.0f0 + 0.5f0)
end

# ie/crown.f label 53: CRVAR/LPIJU change bound + ICRI (CHG=CRNEW-ICR, PDIFPY ±1%/yr, CRMAX cap).
@inline function _ie_crown_label53(crnew::Float32, icr::Int, lstart::Bool, fint::Float32, sp::Int, d::Float32,
                                   h::Float32, htg::Float32)
    chg = crnew - Float32(icr)
    if !lstart || icr > 0
        pdifpy = chg / Float32(icr) / fint
        pdifpy > 0.01f0  && (chg = Float32(icr) * 0.01f0 * fint)
        pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
        crnew = Float32(icr) + chg                       # CRNMLT=1 default
    end
    icri = trunc(Int, crnew + 0.5f0)
    # CRMAX cap (crown.f:556-568) is on the LABEL-53 (CRVAR/LPIJU) path — the NIVAR block ends with GO TO 55 and
    # SKIPS it (crown.f "END OF NI BLOCK"). Skipped at LSTART or ICR==0; CRNMLT=1 default so the ICRI<10 guard fires.
    if !(lstart || icr == 0)
        crln = h * Float32(icr) / 100f0
        crmax = (crln + htg) / (h + htg) * 100f0
        icri < 10 && (icri = trunc(Int, crmax + 0.5f0))
        Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
    end
    return icri
end

# ie/cratet.f:513 `OLDPCT(I)=PCT(I)` — seed the first grow cycle's crown DCR percentile (OLDPCT) from the
# BACKDATED percentile, not the plain inventory PCT. At the initial CRATET the PCT saved into OLDPCT is the one
# computed by the BACKDATING density pass (cratet.f:217-219 `LBKDEN=IDG<2; CALL DENSE`, dense.f), NOT the later
# non-backdated PCT (dense.f:273) that overwrites PCT for the current cycle. The backdating pass (dense.f):
#   • backdates each LIVE tree's diameter to start-of-growth: WK3 = SQRT(D²·R), R = 1-(2·D·G-G²)/D² = (1-G/D)²
#     with G = DG/bark (dense.f:118-128) ⇒ WK3 = D-G; a missing increment (DG≤0 → cratet.f:211 DG=-1) takes
#     R = BAGR, the stand-mean BA-growth ratio over the valid-DG live trees (dense.f:89-117) ⇒ WK3 = D·√BAGR.
#   • builds WK5 = WK3²·PROB and runs PCTILE (pctile.f) over ITRN records, which at this point STILL INCLUDES the
#     standing-dead (cratet.f deletes them only AFTER this DENSE, lines 224-255). HISTORY 6/7 dead enter with
#     WK3 = DBH and PROB inflated ×FINT/FINTM (notre.f:122-124); HISTORY 8/9 dead enter with WK3 = 0 (dense.f:86).
#     PCTILE's INDEX is the CURRENT-DBH descending sort (cratet.f:189).
# Seeding OLDPCT from the plain current PCT instead over-states each tree's start-of-cycle percentile rank ⇒ DCR's
# b13·PB + b14·log(PB) term too high ⇒ EDCR too high / CHG too low ⇒ crown ICR biased low on a subset ⇒ pre-fire
# TREES gap. Verified vs FVSie_g16 (ie_simfire): tree I=3 OLDPCT 34.2→25.2 (oracle backdated), all 27 live match.
"""
    ie_dub_aspen_birthage!(s)

Dub the effective tree AGE (ABIRTH) for inventory aspen-group trees at init, per FVS cratet.f:544-563
(`CALL FINDAG` for ISPC ∈ {18,20,21} = AS/MM/PB when ABIRTH≤0) with ie/findag.f's single Sheppard aspen
age-from-height equation `SITAGE = (H·2.54·12/26.9825)^(1/1.1752)` (H in feet). Only Climate-FVS reads
birth_age in IE (apply_climate_dds! + inlandempire/regent.jl clim_treemult → BIRTHYR=THISYR−ABIRTH → the
Leites XDF/XPP/XWL transfer distance), so this is byte-identical for climate-off IE. Without it, mature
inventory aspen carried birth_age=0 ⇒ BIRTHYR=now ⇒ XRELGR≡1 ⇒ under-applied the climate growth multiplier
(the oracle carries ABIRTH=SITAGE 30-85 yr, XRELGR≈1.03). Aged +FINT/cycle thereafter (gradd.f:205).
"""
function ie_dub_aspen_birthage!(s::StandState)
    t = s.trees
    @inbounds for i in 1:t.n
        sp = Int(t.species[i])
        (sp == 18 || sp == 20 || sp == 21) || continue
        (t.birth_age[i] <= 0f0 && t.height[i] > 0f0) || continue
        t.birth_age[i] = fpow(t.height[i] * 2.54f0 * 12.0f0 / 26.9825f0, 1.0f0 / 1.1752f0)   # ie/findag.f
    end
    return s
end

function ie_seed_backdated_oldpct!(s::StandState)
    t = s.trees; nlive = t.n
    nlive == 0 && return s
    ntot = nlive + Int(t.ndead)
    ba_a = s.calib.bark_a; ba_b = s.calib.bark_b
    # BAGR = mean over the valid-DG live trees of R=1-(2·D·G-G²)/D², G=DG/bark, skipping DG≤0 and G>D (dense.f:89-117).
    bagr = 0f0; sn = 0
    @inbounds for i in 1:nlive
        dg = t.diam_growth[i]; dg <= 0f0 && continue                 # DG≤0 → -1 (missing), excluded (dense.f:100)
        d = t.dbh[i]; d <= 0f0 && continue
        g = dg / bark_ratio(ba_a, ba_b, Int(t.species[i]), d)
        g > d && continue                                            # growth ≥ current diameter → excluded (dense.f:105)
        bagr += 1f0 - (2f0*d*g - g*g)/(d*d); sn += 1
    end
    have_bagr = sn > 0; have_bagr && (bagr /= Float32(sn))
    # WK5 = WK3²·PROB (CHAR for PCTILE). Live: WK3 = backdated diameter (dense.f:118-128). Dead: HISTORY 6/7 →
    # WK3 = DBH with PROB×FINT/FINTM (notre.f); HISTORY 8/9 → WK3 = 0 (dense.f:86); DBH is kept for the sort key.
    fintr = s.control.growth_fintm > 0f0 ? s.control.growth_fint / s.control.growth_fintm : 1f0
    char = Vector{Float32}(undef, ntot)
    key  = Vector{Float32}(undef, ntot)
    @inbounds for i in 1:nlive
        d = t.dbh[i]; key[i] = d
        if !have_bagr || d <= 0f0
            wk3 = d                                                  # SN≤0 ⇒ no backdating (dense.f:112 GO TO 6)
        else
            dg = t.diam_growth[i] > 0f0 ? t.diam_growth[i] : -1f0    # cratet.f:211 missing DG → -1
            g = dg / bark_ratio(ba_a, ba_b, Int(t.species[i]), d)
            g > d && (g = d)                                         # dense.f:124
            r = 1f0 - (2f0*d*g - g*g)/(d*d)
            (g < 0f0 || r <= 0f0) && (r = bagr)                      # dense.f:127
            wk3 = sqrt(d*d*r)
        end
        char[i] = wk3*wk3 * t.tpa[i]
    end
    @inbounds for j in 1:Int(t.ndead)
        i = nlive + j; d = t.dbh[i]; key[i] = d
        h = Int(t.history[i])
        # HISTORY 8/9 (older dead, IMC=9) → WK3=0; HISTORY 6/7 (recent dead, IMC=7) → WK3=DBH, PROB inflated
        char[i] = h >= 8 ? 0f0 : (d*d * t.tpa[i] * fintr)
    end
    # INDEX = current-DBH descending (cratet.f:189 RDPSRT), then PCTILE cumulative-from-smallest (pctile.f).
    idx = Vector{Int32}(undef, ntot)
    _rdpsrt!(key, idx)
    percnt = Vector{Float32}(undef, ntot)
    if ntot == 1
        percnt[1] = 100f0
    else
        cum = 0f0
        @inbounds for k in ntot:-1:1                                 # accumulate from the smallest up (pctile.f:44-54)
            ii = Int(idx[k]); cum += char[ii]; percnt[ii] = cum
        end
        tot = cum
        if tot > 0f0
            @inbounds for ii in 1:ntot; percnt[ii] = percnt[ii] / tot * 100f0; end
        end
        percnt[Int(idx[1])] = 100f0                                  # largest = 100 (pctile.f:73)
    end
    @inbounds for i in 1:nlive; t.old_crown_pct[i] = percnt[i]; end
    return s
end
