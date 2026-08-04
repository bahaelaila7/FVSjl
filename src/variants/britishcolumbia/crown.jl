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
    bc_lv2atv(zone) && error("BC crown_ratio_update!: V2 regime not yet ported (zone=$zone).")
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
