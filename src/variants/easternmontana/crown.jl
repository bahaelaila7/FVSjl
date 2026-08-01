# =============================================================================
# crown.jl (easternmontana) — per-tree CCF (em/ccfcal.f MODE=1).
#
# EM CCF is the direct per-species polynomial (Paine-Hann / NI INT-133 form), summed
# over the tree list = stand CCF (like KT/IE ccfcal, not the CR crown-width→area path).
#   large:  CCFT = RD1 + D·RD2 + D²·RD3
#   small:  CCFT = RDA · D^RDB
# The large/small DBH threshold is species-group dependent (em/ccfcal.f SELECT CASE):
#   sp 5 (LL):                 D≥10 → poly, else RDA·D^RDB           (no 0.001 floor)
#   sp 1-4,6-10,12,17,18:      D≥1  → poly; D>0.1 → RDA·D^RDB; else 0.001
#   sp 11,13-16,19:            D≥10 → poly; D>0.1 → RDA·D^RDB; else 0.001
# Coefficients transcribed verbatim from em/ccfcal.f DATA (species 1=WB..19=OH).
# =============================================================================

const EM_RD1 = Float32[0.0186, 0.0392, 0.0388, 0.01925, 0.03, 0.01925, 0.01925, 0.03, 0.0172,
                       0.0219, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.0204, 0.03]
const EM_RD2 = Float32[0.0146, 0.0180, 0.0269, 0.01676, 0.0216, 0.01676, 0.01676, 0.0173, 0.00876,
                       0.0169, 0.0215, 0.0238, 0.0215, 0.0215, 0.0215, 0.0215, 0.0238, 0.0246, 0.0215]
const EM_RD3 = Float32[0.00288, 0.00207, 0.00466, 0.00365, 0.00405, 0.00365, 0.00365, 0.00259, 0.00112,
                       0.00325, 0.00363, 0.00490, 0.00363, 0.00363, 0.00363, 0.00363, 0.00490, 0.0074, 0.00363]
const EM_RDA = Float32[0.009884, 0.007244, 0.017299, 0.009187, 0.011402, 0.009187, 0.009187, 0.007875,
                       0.011402, 0.007813, 0.011109, 0.008915, 0.011109, 0.011109, 0.011109, 0.011109,
                       0.008915, 0.011109, 0.011109]
const EM_RDB = Float32[1.6667, 1.8182, 1.5571, 1.7600, 1.7560, 1.7600, 1.7600, 1.7360, 1.7560,
                       1.7780, 1.7250, 1.7800, 1.7250, 1.7250, 1.7250, 1.7250, 1.7800, 1.7250, 1.7250]

# em/ccfcal.f MODE=1: per-tree CCF (before the ×P expansion the caller applies).
@inline function em_tree_ccf(sp::Integer, d::Real)::Float32
    d <= 0f0 && return 0f0
    poly()  = EM_RD1[sp] + d * EM_RD2[sp] + d * d * EM_RD3[sp]
    small() = EM_RDA[sp] * d ^ EM_RDB[sp]
    if sp == 5                                   # LL: D≥10 poly, else small (no floor)
        return d >= 10f0 ? poly() : small()
    elseif sp == 11 || (13 <= sp <= 16) || sp == 19   # GA/CW/BA/PW/NC/OH: D≥10 threshold
        return d >= 10f0 ? poly() : (d > 0.1f0 ? small() : 0.001f0)
    else                                          # 1-4,6-10,12,17,18: D≥1 threshold
        return d >= 1f0 ? poly() : (d > 0.1f0 ? small() : 0.001f0)
    end
end

# =============================================================================
# EM crown-ratio update (em/crown.f) — the IE/NI Weibull/DCR change-in-crown model.
# Single PARM(14) vector (IE/NI sp9 form, applied to all species) + habitat intercept
# CRHAB[MAPHAB[ITYPE]] + CRSD=6.35. Mirrors KT's crown_ratio_update! DCR structure
# (crown change = exp(PCR) − exp(DCR-backdated), bounded ±1%/yr; d<3 keeps crown while cycling).
# =============================================================================
# em/crown.f PARM(14): 1-6 density(BA,BA²,lnBA,RELDEN,RELDEN²,lnRELDEN), 7-14(D,D²,lnD,H,H²,lnH,P,lnP)
const EM_CRPARM = Float32[-0.00190, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.23372, 0.0, 0.0, -0.28433, 0.001903, 0.0]
const EM_CRHAB  = Float32[0.09453, -0.07740, 0.07113, 0.2039, 0.06176, 0.1513, 0.09086, 0.1580, 0.09229, 0.01551, 0.0, 0.0, 0.0, 0.0]
const EM_CR_MAPHAB = Int[2,2,2,2,2,2,2,2,2,2,2,2,2,3,4,4,4,4,5,6,6,7,6,1,8,1,9,10,10,6]  # ITYPE→CRHAB idx (13*2,3,4*4,5,6,6,7,6,1,8,1,9,2*10,6)
const _EM_CRSD = 6.35f0

function crown_ratio_update!(s::StandState, ::EasternMontana; fint::Float32 = 10.0f0, lstart::Bool = false, kwargs...)
    p, t = s.plot, s.trees
    t.n == 0 && return s
    itype = Int(p.habitat_input); it = (1 <= itype <= 30) ? itype : 1
    crcon = EM_CRHAB[clamp(EM_CR_MAPHAB[it], 1, 14)]
    ba = p.basal_area; relden = p.relative_density
    lnba = ba > 0f0 ? log(ba) : 0f0; lnrd = relden > 0f0 ? log(relden) : 0f0
    reldm1 = p.relative_density_prev; oba = p.old_ba
    if reldm1 < 100f0; oba = ba; reldm1 = relden; end
    x1 = (!lstart && oba > 0f0) ? log(oba) : 0f0
    x2 = (!lstart && reldm1 > 0f0) ? log(reldm1) : 0f0
    dgsd = s.control.dg_sd
    ba_a = s.calib.bark_a; ba_b = s.calib.bark_b
    P = EM_CRPARM
    nlim = t.n + (lstart ? Int(t.ndead) : 0)
    @inbounds for i in 1:nlim
        t.tpa[i] <= 0f0 && continue
        icr = Int(t.crown_pct[i])
        (lstart && icr > 0) && continue
        icr < 0 && (t.crown_pct[i] = Int32(-icr); continue)
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        bark = bark_ratio(ba_a, ba_b, sp, d)
        local icri::Int
        if d >= 3.0f0
            xcrcon = crcon + P[1]*ba + P[2]*ba*ba + P[3]*lnba + P[4]*relden + P[5]*relden*relden + P[6]*lnrd
            pp = t.crown_ratio[i]; pp < 0.01f0 && (pp = 0.01f0)
            pcr = xcrcon + P[7]*d + P[8]*d*d + P[9]*log(d) + P[10]*h + P[11]*h*h + P[12]*log(h) + P[13]*pp + P[14]*log(pp)
            exppcr = exp(pcr)
            if lstart
                icri = trunc(Int, icr + exppcr*100f0 + 0.50005f0)
                dgsd >= 1.0f0 && (icri = trunc(Int, bachlo(s.rng, Float32(icri), _EM_CRSD)))
            else
                dcrcon = crcon + P[1]*oba + P[2]*oba*oba + P[3]*x1 + P[4]*reldm1 + P[5]*reldm1*reldm1 + P[6]*x2
                db = d - t.diam_growth[i]/bark; db <= 0f0 && (db = d)
                hb = h - t.ht_growth[i]; hb <= 0f0 && (hb = h)
                pb = t.crown_ratio[i]; pb < 0.01f0 && (pb = 0.01f0)
                dcr = dcrcon + P[7]*db + P[8]*db*db + P[9]*log(db) + P[10]*hb + P[11]*hb*hb + P[12]*log(hb) + P[13]*pb + P[14]*log(pb)
                chg = exppcr - exp(dcr)
                if icr > 0
                    pdifpy = chg / Float32(icr) / fint * 100f0
                    pdifpy > 0.01f0  && (chg = Float32(icr) * 0.01f0 * fint / 100f0)
                    pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint / 100f0)
                end
                icri = trunc(Int, Float32(icr) + chg*100f0 + 0.50005f0)
            end
        else
            lstart || continue                       # cycling: D<3 keeps its crown
            icri = icr > 0 ? icr : 40                 # minimal lstart dub (emt01 inventory crowns present)
        end
        icri > 95 && (icri = 95); icri < 5 && (icri = 5)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end
