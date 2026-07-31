# =============================================================================
# crown.jl (kootenai) — KT per-tree CCF (kt/ccfcal.f MODE=1). The stand CCF = Σ CCFT is the RELDEN
# that the DG (dgf!) and height (htgf) chunks read. UNLIKE the eastern/CR stand_ccf (crown-width → area),
# KT's CCF is a DIRECT per-species polynomial (Bush/Crookston, kt/ccfcal.f):
#   D >= 10:  CCFT = RD1(sp) + D*RD2(sp) + D^2*RD3(sp)
#   D <  10:  CCFT = RDA(sp) * D^RDB(sp)
#   CCFT *= P (tree TPA);  stand CCF = Σ CCFT.
# Crown-WIDTH (ccfcal MODE=2, B1..B6) is the crown-ratio/crown model — separate; lands with crown_ratio!.
# =============================================================================

const KT_RD1 = Float32[0.03,0.02,0.11,0.04,0.03,0.03,0.01925,0.03,0.03,0.03,0.03]
const KT_RD2 = Float32[0.0167,0.0148,0.0333,0.0270,0.0215,0.0238,0.01676,0.0173,0.0216,0.0180,0.0215]
const KT_RD3 = Float32[0.00230,0.00338,0.00259,0.00405,0.00363,0.00490,0.00365,0.00259,0.00405,0.00281,0.00363]
const KT_RDA = Float32[0.009884,0.007244,0.017299,0.015248,0.011109,0.008915,0.009187,0.007875,0.011402,0.007813,0.011109]
const KT_RDB = Float32[1.6667,1.8182,1.5571,1.7333,1.7250,1.7800,1.7600,1.7360,1.7560,1.7680,1.7250]

"""Per-tree CCF contribution (kt/ccfcal.f MODE=1), excluding the ×TPA factor."""
@inline function kt_tree_ccf(sp::Integer, d::Real)::Float32
    d <= 0f0 && return 0f0
    if d >= 10f0
        return KT_RD1[sp] + d * KT_RD2[sp] + d * d * KT_RD3[sp]
    else
        return KT_RDA[sp] * d ^ KT_RDB[sp]
    end
end
# crown.f PARM(sp,1:14): 1-6 density terms (BA,BA²,lnBA,RELDEN,RELDEN²,lnRELDEN), 7-14 = B7..B14 (D,D²,lnD,H,H²,lnH,P,lnP). [sp][col]
const KT_CRPARM = Float32[
    0f0 0f0 -0.34566f0 0f0 0f0 0f0 0.03882f0 -0.0007f0 0f0 0f0 0f0 -0.21217f0 0.00301f0 0f0;
    -0.00204f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0 0.30066f0 0f0 0f0 -0.59302f0 0f0 0.19558f0;
    0f0 0f0 0f0 0f0 0f0 -0.15334f0 0f0 0f0 0.3384f0 0f0 0f0 -0.59685f0 0f0 0.16488f0;
    -0.00183f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0 0.24293f0 0f0 0f0 -0.25601f0 0f0 0.0726f0;
    0f0 -1.902f-6 0f0 0f0 0f0 0f0 0.03027f0 -0.00055f0 0f0 0f0 0f0 -0.25776f0 0f0 0.06887f0;
    0f0 0f0 0.17479f0 -0.00183f0 0f0 0f0 -0.0056f0 0f0 0f0 0f0 0f0 0f0 0f0 0.1105f0;
    0f0 0f0 0f0 0f0 0f0 -0.18555f0 0f0 0f0 0.53172f0 -0.02989f0 0.00011f0 0f0 0.0042f0 0f0;
    -0.00203f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0 0.29699f0 0f0 0f0 -0.38334f0 0f0 0.09918f0;
    -0.0019f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0 0.23372f0 0f0 0f0 -0.28433f0 0.001903f0 0f0;
    -0.002165f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0 0.26558f0 0f0 0f0 -0.31555f0 0f0 0.16072f0;
    -0.00264f0 0f0 0f0 0f0 5.116f-6 0f0 0f0 0f0 0f0 0f0 0f0 -0.25138f0 0f0 0.0514f0]

# crown.f CRHAB(icrhab,sp) habitat intercepts → CRCON. [sp][icrhab 1..14]
const KT_CRHAB = Float32[
    0.8884f0 0.7309f0 0.9347f0 0.9888f0 0.9945f0 1.1126f0 1.0263f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0;
    0.06533f0 0.03441f0 0.2307f0 0.1661f0 -0.1253f0 -0.05018f0 0.11005f0 0.08113f0 0.1782f0 0.03919f0 0.2107f0 0f0 0f0 0f0;
    0.8643f0 0.7271f0 0.984f0 0.8127f0 0.8874f0 0.7055f0 0.7708f0 0.7849f0 0.8038f0 0.8742f0 0.8232f0 0.8415f0 0.9759f0 0f0;
    -0.2304f0 -0.5421f0 -0.4343f0 -0.3759f0 -0.4129f0 -0.4879f0 -0.2674f0 -0.1941f0 0f0 0f0 0f0 0f0 0f0 0f0;
    -0.2413f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0;
    -1.6053f0 -1.7128f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0;
    -0.3785f0 -0.4142f0 -0.3985f0 -0.2987f0 -0.381f0 -0.4087f0 -0.3577f0 -0.2994f0 -0.2486f0 -0.2863f0 -0.1968f0 -0.4931f0 -0.2676f0 -0.5625f0;
    0.05351f0 -0.05031f0 0.1075f0 -0.1872f0 0.01729f0 0.03667f0 0.01885f0 0.09102f0 0.1371f0 0.08368f0 0.123f0 -0.02365f0 0f0 0f0;
    0.09453f0 -0.0774f0 0.07113f0 0.2039f0 0.06176f0 0.1513f0 0.09086f0 0.158f0 0.09229f0 0.01551f0 0f0 0f0 0f0 0f0;
    -0.9436f0 -0.8654f0 -0.8849f0 -0.9067f0 -0.8783f0 -1.0103f0 -1.0268f0 -1.005f0 -1.0301f0 0f0 0f0 0f0 0f0 0f0;
    0.4649f0 0.3211f0 0.197f0 0.2295f0 0.3383f0 0.345f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0 0f0]

# crown.f MAPHAB(itype,sp) → icrhab (row into CRHAB). [itype 1..30][sp]
const KT_CR_MAPHAB = Int[
    2 2 2 2 1 1 2 2 2 2 1;
    2 2 2 2 1 1 2 2 2 2 1;
    2 2 2 2 1 1 2 2 2 4 1;
    2 2 4 2 1 1 2 2 2 1 1;
    2 2 4 2 1 1 2 2 2 1 1;
    2 2 4 2 1 1 2 2 2 1 1;
    2 2 6 2 1 1 4 2 2 5 1;
    2 3 7 2 1 1 5 3 2 6 1;
    2 2 4 2 1 1 5 2 2 1 1;
    2 4 8 1 1 1 2 1 2 1 1;
    2 4 8 1 1 1 2 1 2 1 1;
    2 5 5 2 1 1 6 2 2 8 1;
    3 6 9 3 1 1 7 4 2 7 2;
    4 7 10 4 1 1 8 5 3 9 2;
    4 7 10 4 1 1 8 5 4 9 2;
    4 7 10 4 1 1 8 5 4 9 2;
    5 8 11 5 1 2 9 6 4 3 3;
    5 8 11 5 1 2 9 6 4 3 3;
    5 4 8 6 1 2 10 7 5 3 4;
    6 1 1 1 1 1 11 8 6 1 1;
    6 10 12 7 1 1 11 8 6 1 1;
    1 9 12 7 1 1 12 9 7 1 1;
    6 10 13 7 1 1 11 8 6 1 5;
    1 1 1 1 1 1 1 10 1 1 1;
    7 11 3 8 1 1 13 11 8 1 6;
    1 1 1 1 1 1 1 1 1 1 1;
    6 1 3 7 1 1 14 1 9 1 1;
    6 1 1 1 1 1 3 12 10 1 1;
    6 1 1 1 1 1 3 12 10 1 1;
    6 2 1 1 1 1 11 8 6 1 1]

const KT_DUB_BCR0 = Float32[-0.44316f0, -0.83965f0, -0.89122f0, -0.62646f0, -0.49548f0, 0.11847f0, -0.32466f0, -0.92007f0, -0.89014f0, -0.17561f0, -0.49548f0]
const KT_DUB_BCR1 = Float32[-0.48446f0, -0.16106f0, -0.18082f0, -0.06141f0, 0.00012f0, -0.39305f0, -0.20108f0, -0.22454f0, -0.18026f0, -0.33847f0, 0.00012f0]
const KT_DUB_BCR2 = Float32[0.05825f0, 0.04161f0, 0.05186f0, 0.0236f0, 0.00362f0, 0.02783f0, 0.04219f0, 0.03248f0, 0.02233f0, 0.05699f0, 0.00362f0]
const KT_DUB_BCR3 = Float32[0.00513f0, 0.00602f0, 0.00454f0, 0.00505f0, 0.00456f0, 0.00626f0, 0.00436f0, 0.0062f0, 0.00614f0, 0.00692f0, 0.00456f0]
const KT_DUB_CRSD = Float32[0.9476f0, 0.7396f0, 0.8706f0, 0.9203f0, 0.945f0, 0.8012f0, 0.7707f0, 0.9721f0, 0.8871f0, 0.8866f0, 0.945f0]
const KT_CRSD_LARGE = 6.35f0   # crown.f CRSD (LSTART dub stochastic spread)

# kt/dubscr.f: small-tree (D<3) crown-ratio dub — logistic. Stochastic FCR when DGSD≥1 (RNG-cornered).
function kt_dubscr(sp::Int, d::Float32, h::Float32, ba::Float32, dgsd::Float32, rng)::Int
    cr = KT_DUB_BCR0[sp] + KT_DUB_BCR1[sp]*d + KT_DUB_BCR2[sp]*h + KT_DUB_BCR3[sp]*ba
    sdv = KT_DUB_CRSD[sp]
    fcr = 0.0f0
    if dgsd >= 1.0f0
        while true
            fcr = bachlo(rng, 0.0f0, sdv)
            abs(fcr) <= sdv && break
        end
    end
    abs(cr + fcr) >= 86.0f0 && (cr = 86.0f0)
    crf = 1.0f0 / (1.0f0 + exp(cr + fcr))
    crf < 0.05f0 && (crf = 0.05f0); crf > 0.95f0 && (crf = 0.95f0)
    return trunc(Int, crf * 100.0f0 + 0.5f0)
end

"""
    crown_ratio_update!(s, ::Kootenai; fint, lstart, ...)

KT crown (kt/crown.f): DUB missing crowns at LSTART, UPDATE crowns each cycle. Large trees (D≥3) use the
PCR/DCR change model; D<3 use DUBSCR (LSTART only; cycling keeps them). Bounds [5,95]. The LSTART dub carries a
stochastic BACHLO (DGSD≥1) — RNG-cornered. Backdate P uses the current PCT (OLDPCT not yet threaded).
"""
function crown_ratio_update!(s::StandState, ::Kootenai; fint::Float32 = 10.0f0, lstart::Bool = false, kwargs...)
    p, t, c = s.plot, s.trees, s.calib
    t.n == 0 && return s
    itype = Int(p.habitat_input); it = (1 <= itype <= 30) ? itype : 1
    ba = p.basal_area; relden = p.relative_density
    lnba = ba > 0f0 ? log(ba) : 0f0; lnrd = relden > 0f0 ? log(relden) : 0f0
    reldm1 = p.relative_density_prev; oba = p.old_ba
    if reldm1 < 100f0; oba = ba; reldm1 = relden; end
    x1 = (!lstart && oba > 0f0) ? log(oba) : 0f0
    x2 = (!lstart && reldm1 > 0f0) ? log(reldm1) : 0f0
    dgsd = s.control.dg_sd
    ba_a = c.bark_a; ba_b = c.bark_b
    nlim = t.n + (lstart ? Int(t.ndead) : 0)
    @inbounds for i in 1:nlim
        t.tpa[i] <= 0f0 && continue
        icr = Int(t.crown_pct[i])
        (lstart && icr > 0) && continue
        icr < 0 && (t.crown_pct[i] = Int32(-icr); continue)
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        bark = bark_ratio(ba_a, ba_b, sp, d)
        crcon = KT_CRHAB[sp, clamp(Int(KT_CR_MAPHAB[it, sp]), 1, 14)]
        local icri::Int
        if d >= 3.0f0
            xcrcon = crcon + KT_CRPARM[sp,1]*ba + KT_CRPARM[sp,2]*ba*ba + KT_CRPARM[sp,3]*lnba +
                     KT_CRPARM[sp,4]*relden + KT_CRPARM[sp,5]*relden*relden + KT_CRPARM[sp,6]*lnrd
            b7=KT_CRPARM[sp,7]; b8=KT_CRPARM[sp,8]; b9=KT_CRPARM[sp,9]; b10=KT_CRPARM[sp,10]
            b11=KT_CRPARM[sp,11]; b12=KT_CRPARM[sp,12]; b13=KT_CRPARM[sp,13]; b14=KT_CRPARM[sp,14]
            pp = t.crown_ratio[i]; pp < 0.01f0 && (pp = 0.01f0)
            pcr = xcrcon + b7*d + b8*d*d + b9*log(d) + b10*h + b11*h*h + b12*log(h) + b13*pp + b14*log(pp)
            exppcr = exp(pcr)
            local chg::Float32
            if lstart
                chg = exppcr
                icri = trunc(Int, icr + chg*100f0 + 0.50005f0)
                if dgsd >= 1.0f0
                    icri = trunc(Int, bachlo(s.rng, Float32(icri), KT_CRSD_LARGE))
                end
            else
                dcrcon = crcon + KT_CRPARM[sp,1]*oba + KT_CRPARM[sp,2]*oba*oba + KT_CRPARM[sp,3]*x1 +
                         KT_CRPARM[sp,4]*reldm1 + KT_CRPARM[sp,5]*reldm1*reldm1 + KT_CRPARM[sp,6]*x2
                db = d - t.diam_growth[i]/bark; db <= 0f0 && (db = d)
                hb = h - t.ht_growth[i]; hb <= 0f0 && (hb = h)
                pb = t.crown_ratio[i]; pb < 0.01f0 && (pb = 0.01f0)     # OLDPCT ≈ current PCT
                dcr = dcrcon + b7*db + b8*db*db + b9*log(db) + b10*hb + b11*hb*hb + b12*log(hb) + b13*pb + b14*log(pb)
                expdcr = exp(dcr)
                chg = exppcr - expdcr
                if icr > 0                                             # bound ±1%/yr (crown.f:329-332)
                    pdifpy = chg / Float32(icr) / fint * 100f0
                    pdifpy > 0.01f0  && (chg = Float32(icr) * 0.01f0 * fint / 100f0)
                    pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint / 100f0)
                end
                icri = trunc(Int, Float32(icr) + chg*100f0 + 0.50005f0)
            end
        else
            lstart || continue                                        # cycling: D<3 keeps its crown
            icri = kt_dubscr(sp, d, h, ba, dgsd, s.rng)
        end
        icri > 95 && (icri = 95); icri < 5 && (icri = 5)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end
