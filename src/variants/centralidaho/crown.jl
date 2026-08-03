# =============================================================================
# crown.jl (centralidaho) — CI per-tree CCF (ci/ccfcal.f MODE=1). Chunk 5 (CCF portion).
# The stand CCF = Σ CCFT·TPA is the RELDEN the DG (dgf!) CONSPP term reads. Same western
# polynomial as KT/IE (Bush/Crookston): D≥10 → RD1+D·RD2+D²·RD3; D<10 → RDA·D^RDB. CI's RD*
# coefficients (ci/ccfcal.f DATA) match KT for the shared N-Rockies conifers.
# (Crown-WIDTH MODE=2 B1..B6 + crown-ratio model land separately with crown_ratio_update!.)
# =============================================================================

const CI_RD1 = Float32[0.03,0.02,0.11,0.04,0.03,0.03,0.01925,0.03,0.03,0.03,0.01925,0.01925,0.03,0.01925,0.0204,0.01925,0.03,0.03,0.03]
const CI_RD2 = Float32[0.0167,0.0148,0.0333,0.027,0.0215,0.0238,0.01676,0.0173,0.0216,0.018,0.01676,0.01676,0.0238,0.01676,0.0246,0.01676,0.0215,0.0215,0.0215]
const CI_RD3 = Float32[0.0023,0.00338,0.00259,0.00405,0.00363,0.0049,0.00365,0.00259,0.00405,0.00281,0.00365,0.00365,0.0049,0.00365,0.0074,0.00365,0.00363,0.00363,0.00363]
const CI_RDA = Float32[0.009884,0.007244,0.017299,0.015248,0.011109,0.008915,0.009187,0.007875,0.011402,0.007813,0.009187,0.009187,0.008915,0.009187,0.0,0.009187,0.011109,0.011109,0.011109]
const CI_RDB = Float32[1.6667,1.8182,1.5571,1.7333,1.725,1.78,1.76,1.736,1.756,1.768,1.76,1.76,1.78,1.76,0.0,1.76,1.725,1.725,1.725]

"""Per-tree CCF contribution (ci/ccfcal.f MODE=1), excluding the ×TPA factor."""
@inline function ci_tree_ccf(sp::Integer, d::Real)::Float32
    dd = Float32(d)
    if 11 <= sp <= 16                                   # ccfcal.f CASE(11..16): poly at D≥1.0 (NOT ≥10)
        if dd >= 1.0f0
            return CI_RD1[sp] + dd * CI_RD2[sp] + dd * dd * CI_RD3[sp]
        elseif dd > 0.1f0
            return sp == 15 ? dd * (CI_RD1[sp] + CI_RD2[sp] + CI_RD3[sp]) : CI_RDA[sp] * dd ^ CI_RDB[sp]
        else
            return sp == 15 ? dd * (CI_RD1[sp] + CI_RD2[sp] + CI_RD3[sp]) : 0.001f0
        end
    elseif sp == 17 || sp == 19                         # ccfcal.f CASE(17,19): poly at D≥10, 0.001 floor at D≤0.1
        if dd >= 10f0
            return CI_RD1[sp] + dd * CI_RD2[sp] + dd * dd * CI_RD3[sp]
        elseif dd > 0.1f0
            return CI_RDA[sp] * dd ^ CI_RDB[sp]
        else
            return 0.001f0
        end
    else                                                # DEFAULT (1-10,18): poly at D≥10, else RDA·D^RDB
        dd <= 0f0 && return 0f0
        return dd >= 10f0 ? (CI_RD1[sp] + dd * CI_RD2[sp] + dd * dd * CI_RD3[sp]) : CI_RDA[sp] * dd ^ CI_RDB[sp]
    end
end

# ci/crown.f Weibull crown ratio (chunk 5b). ACRNEW=C0+C1·RELSDI·100; A=WEIBA; B=WEIBB0+WEIBB1·ACRNEW;
# C=WEIBC0 (WEIBC1=0); crnew=(A+B·(−ln(1−x))^(1/C))·10, x=ISORT/N·SCALE. CRNMLT=1/DLOW=0/DHI=99 (no-op mods).
const CI_WEIBA = Float32[2.0,0.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,0.0,1.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0]
const CI_WEIBB0 = Float32[-2.12713,0.07609,-1.19297,-1.19297,0.06593,-1.38636,0.07609,-0.91567,-0.91567,0.24916,-0.82631,-0.82631,-0.08414,0.0,-0.2383,-0.82631,0.0,0.07609,0.0]
const CI_WEIBB1 = Float32[1.10526,1.10184,1.12928,1.12928,1.09624,1.16801,1.10184,1.06469,1.06469,1.04831,1.06217,1.06217,1.14765,0.0,1.18016,1.06217,0.0,1.10184,0.0]
const CI_WEIBC0 = Float32[2.77,3.01,3.42,3.42,3.71,3.02,3.01,3.5,3.5,4.36,3.31429,3.31429,2.775,0.0,3.04,3.31429,0.0,3.01,0.0]
const CI_CRC0 = Float32[7.16846,5.50719,5.52653,5.52653,6.61291,6.17373,5.50719,6.774,6.12779,6.41166,6.19911,6.19911,4.01678,0.0,4.62512,6.19911,0.0,7.238,0.0]
const CI_CRC1 = Float32[-0.02375,-0.01833,0.0,0.0,-0.02182,-0.01795,-0.01833,0.0,-0.01269,-0.02041,-0.02216,-0.02216,-0.01516,0.0,-0.01604,-0.02216,0.0,0.0,0.0]

# --- ci/dubscr.f : small-tree (DBH<1") + regen crown-ratio dub. Logistic in DBH/H/BA (+ TPCCF/AVH/MAI for
# the 4 MAI-species 11,12,13,16). CR = 1/(1+exp(arg)); sp15 uses the linear MC map. Bounded [0.05,0.95].
# The BACHLO random error (active because DGSD=1.7≥1, ci/grinit.f:176) is DEFERRED — deterministic mean only
# (stochastic ZZRAN/BACHLO class, accepted-cornered). RMAI (BCR10, sp 11/12/13/16 only) also deferred → tmai=0.
const CI_BCR0  = Float32[-0.44316,-0.83965,-0.89122,-0.62646,-0.49548,0.11847,-0.32466,-0.92007,-0.89014,-0.17561,-1.66949,-1.66949,-0.426688,-2.19723,5.0,-1.66949,0.0,-0.49548,0.0]
const CI_BCR1  = Float32[-0.48446,-0.16106,-0.18082,-0.06141,0.00012,-0.39305,-0.20108,-0.22454,-0.18026,-0.33847,-0.209765,-0.209765,-0.093105,0.0,0.0,-0.209765,0.0,0.00012,0.0]
const CI_BCR2  = Float32[0.05825,0.04161,0.05186,0.02360,0.00362,0.02783,0.04219,0.03248,0.02233,0.05699,0.0,0.0,0.022409,0.0,0.0,0.0,0.0,0.00362,0.0]
const CI_BCR3  = Float32[0.00513,0.00602,0.00454,0.00505,0.00456,0.00626,0.00436,0.00620,0.00614,0.00692,0.003359,0.003359,0.002633,0.0,0.0,0.003359,0.0,0.00456,0.0]
const CI_BCR5  = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.011032,0.011032,0.0,0.0,0.0,0.011032,0.0,0.0,0.0]
const CI_BCR6  = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.045532,0.0,0.0,0.0,0.0,0.0,0.0]
const CI_BCR8  = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.017727,0.017727,0.0,0.0,0.0,0.017727,0.0,0.0,0.0]
const CI_BCR9  = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.000053,-0.000053,0.000022,0.0,0.0,-0.000053,0.0,0.0,0.0]
const CI_BCR10 = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.014098,0.014098,-0.013115,0.0,0.0,0.014098,0.0,0.0,0.0]
const CI_CRSD  = Float32[0.9476,0.7396,0.8706,0.9203,0.9450,0.8012,0.7707,0.9721,0.8871,0.8866,0.5,0.5,0.9310,0.2,0.5,0.5,0.0,0.9450,0.0]

"""
ci/dubscr.f DUBSCR — dubbed crown ratio (fraction 0.05–0.95) for a small (<1") or regen tree.
A BACHLO random error FCR (rejected if |FCR|>SD) perturbs the LOGIT argument before the transform
(DGSD=1.7≥1 ⇒ the draw always fires). Consumes the main RANN stream, matching live CROWN order.
"""
@inline function ci_dubscr(rng, sp::Integer, d::Real, h::Real, ba::Real, tpccf::Real, avh::Real, tmai::Real)::Float32
    hf = Float32(h); hf <= 0f0 && (hf = 0.1f0)
    cr = CI_BCR0[sp] + CI_BCR1[sp]*Float32(d) + CI_BCR2[sp]*hf + CI_BCR3[sp]*Float32(ba) +
         CI_BCR5[sp]*Float32(tpccf) + CI_BCR6[sp]*(Float32(avh)/hf) + CI_BCR8[sp]*Float32(avh) +
         CI_BCR9[sp]*(Float32(ba)*Float32(tpccf)) + CI_BCR10[sp]*Float32(tmai)
    sd = CI_CRSD[sp]
    fcr = 0f0
    while true                                     # dubscr.f label 10: FCR=BACHLO(0,SD); reject |FCR|>SD
        fcr = bachlo(rng, 0f0, sd)                 # DGSD=1.7≥1 ⇒ always draws (sd=0 ⇒ bachlo returns 0)
        abs(fcr) > sd && continue
        break
    end
    if sp == 15
        cr = ((cr - 1f0)*10f0 + 1f0) / 100f0
    else
        abs(cr + fcr) >= 86f0 && (cr = 86f0)       # overflow guard (faithful: sets +86 regardless of sign)
        cr = 1f0 / (1f0 + exp(cr + fcr))
    end
    cr < 0.05f0 && (cr = 0.05f0); cr > 0.95f0 && (cr = 0.95f0)
    return cr
end

function crown_ratio_update!(s::StandState, ::CentralIdaho; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    relden = p.relative_density; sdiac = crown_sdi
    sd = s.coef.species
    p_pccf = s.density.point_ccf
    key = Vector{Float32}(undef, n); idx = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        bk = ci_bratio(sd, Int(t.species[i]), t.dbh[i])
        key[i] = t.dbh[i] + t.diam_growth[i] / bk; idx[i] = Int32(i)
    end
    _rdpsrt!(key, idx; lseq = false)
    isort = Vector{Int32}(undef, n)
    @inbounds for jj in 1:n; isort[idx[jj]] = Int32(n - jj + 1); end
    @inbounds for i in 1:n
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (sp < 1 || sp > 19) && continue
        (lstart && t.crown_pct[i] > 0) && continue
        if sp == 17 || sp == 19                         # ci/crown.f CASE(17,19): CW/OH crown model at ALL sizes
            hf = h + t.ht_growth[i]; hf <= 0f0 && (hf = 0.1f0)   # HF=H+HTG (HTG=0 at the lstart dub)
            cl = 5.17281f0 + 0.32552f0*hf - 0.01675f0*p.basal_area
            cl < 1f0 && (cl = 1f0); cl > hf && (cl = hf)
            icri = trunc(Int, (cl/hf)*100f0 + 0.5f0)
            if !lstart && t.crown_pct[i] > 0            # cycling: limit change to 1%/yr (label 53)
                icr0 = Float32(t.crown_pct[i]); chg = Float32(icri) - icr0; pdifpy = chg/icr0/fint
                pdifpy > 0.01f0 && (chg = icr0*0.01f0*fint); pdifpy < -0.01f0 && (chg = icr0*(-0.01f0)*fint)
                icri = trunc(Int, icr0 + chg + 0.5f0)
            end
            icri > 95 && (icri = 95); icri < 10 && (icri = 10)
            t.crown_pct[i] = Int32(icri)
            continue
        end
        # ci/crown.f:290 — DBH<1" at LSTART jumps to label 58: dub via DUBSCR (small-tree logistic),
        # NOT the Weibull path, and NOT skipped. Floor to 10 (CRNMLT=1). (17/19 use the Weibull-loop
        # CW/OH model at all sizes; here they fall through to Weibull like the rest — sp15 handled in dubscr.)
        if d < 1f0 && lstart
            pt = Int(t.plot_id[i])
            tpccf = (1 <= pt <= length(p_pccf)) ? p_pccf[pt] : 0f0
            cr = ci_dubscr(s.rng, sp, d, t.height[i], p.basal_area, tpccf, p.avg_height, 0f0)
            icri = trunc(Int, cr*100f0 + 0.5f0)
            icri < 10 && (icri = 10); icri > 95 && (icri = 95); icri < 1 && (icri = 1)
            t.crown_pct[i] = Int32(icri)
            continue
        end
        icr = Int(t.crown_pct[i])
        relsdi = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0
        relsdi > 1.5f0 && (relsdi = 1.5f0)
        acrnew = CI_CRC0[sp] + CI_CRC1[sp] * relsdi * 100f0
        A = CI_WEIBA[sp]
        B = CI_WEIBB0[sp] + CI_WEIBB1[sp] * acrnew; B < 1f0 && (B = 1f0)
        C = CI_WEIBC0[sp]; C < 2f0 && (C = 2f0)
        scale = 1f0 - 0.00167f0 * (relden - 100f0)
        scale > 1f0 && (scale = 1f0); scale < 0.30f0 && (scale = 0.30f0)
        x = d > 0f0 ? (Float32(isort[i]) / Float32(n)) * scale : 0.5f0 * scale
        x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
        crnew = (A + B * (-log(1f0 - x))^(1f0 / C)) * 10f0
        if !(lstart || icr == 0)
            chg = crnew - Float32(icr); pdifpy = chg / Float32(icr) / fint
            pdifpy > 0.01f0 && (chg = Float32(icr) * 0.01f0 * fint)
            pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
            crnew = Float32(icr) + chg                      # CRNMLT=1, DLOW=0/DHI=99 ⇒ no-op branch
        end
        icri = trunc(Int, crnew + 0.5f0)
        if !(lstart || icr == 0)
            crln = h * Float32(icr) / 100f0; htg = t.ht_growth[i]
            crmax = (h + htg) > 0f0 ? (crln + htg) / (h + htg) * 100f0 : 100f0
            Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
        end
        icri > 95 && (icri = 95); icri < 10 && (icri = 10)   # ci/crown.f:405-406 label-59 floor (CRNMLT=1)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end
