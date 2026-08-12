# =============================================================================
# crown.jl (teton) — TT per-tree CCF (tt/ccfcal.f MODE=1). Chunk 3 prerequisite
# (RELDEN = stand CCF + PCCF feed the DG), also reused by chunk 5 (crown).
#
# tt/ccfcal.f: BREAK=10 for PP(10)/NC(15)/OH(18), else 1. Per-tree CCF (before ×P):
#   D ≥ BREAK           : RD1 + D·RD2 + D²·RD3                          (large-tree poly)
#   0.1 < D < BREAK     : sp13/16 → D·(RD1+RD2+RD3) ; else RDA·D^RDB   (small-tree)
#   D ≤ 0.1             : sp13/16 → D·(RD1+RD2+RD3) ; sp10 → RDA·D^RDB ; else 0.001
# Returns the per-tree CCF (the caller multiplies by tpa for the stand/point CCF sum),
# matching em_tree_ccf. Verified vs raw tt/ccfcal.f DATA (WB RD1=.01925, DF RD1=.11).
# =============================================================================

const TT_RD1 = Float32[0.01925, 0.01925, 0.11, 0.01925, 0.03, 0.03, 0.01925, 0.03, 0.03, 0.03, 0.01925, 0.01925, 0.0204, 0.03, 0.03, 0.0204, 0.01925, 0.03]
const TT_RD2 = Float32[0.0168, 0.0168, 0.0333, 0.01676, 0.0173, 0.0238, 0.0168, 0.0173, 0.0216, 0.018, 0.01676, 0.01676, 0.0246, 0.0238, 0.0215, 0.0246, 0.0168, 0.0215]
const TT_RD3 = Float32[0.00365, 0.00365, 0.00259, 0.00365, 0.00259, 0.0049, 0.00365, 0.00259, 0.00405, 0.00281, 0.00365, 0.00365, 0.0074, 0.0049, 0.00363, 0.0074, 0.00365, 0.00363]
const TT_RDA = Float32[0.009187, 0.009187, 0.017299, 0.009187, 0.007875, 0.008915, 0.009187, 0.007875, 0.011402, 0.007813, 0.009187, 0.009187, 0.0, 0.008915, 0.011109, 0.0, 0.009187, 0.011109]
const TT_RDB = Float32[1.76, 1.76, 1.5571, 1.76, 1.736, 1.78, 1.76, 1.736, 1.756, 1.768, 1.76, 1.76, 0.0, 1.78, 1.725, 0.0, 1.76, 1.725]

@inline function tt_tree_ccf(sp::Integer, d::Real)::Float32
    d <= 0f0 && return 0f0
    poly()  = TT_RD1[sp] + d * TT_RD2[sp] + d * d * TT_RD3[sp]
    small() = TT_RDA[sp] * d ^ TT_RDB[sp]
    brk = (sp == 10 || sp == 15 || sp == 18) ? 10f0 : 1f0
    if d >= brk
        return poly()
    elseif d > 0.1f0
        return (sp == 13 || sp == 16) ? d * (TT_RD1[sp] + TT_RD2[sp] + TT_RD3[sp]) : small()
    else
        (sp == 13 || sp == 16) && return d * (TT_RD1[sp] + TT_RD2[sp] + TT_RD3[sp])
        sp == 10 && return small()
        return 0.001f0
    end
end

# ============================ crown ratio (tt/crown.f) — chunk 5 ============================
# Rank-based Weibull crown-ratio model. ACRNEW=C0+C1·RELSDI·100 (mean CR%); Weibull A/B/C from ACRNEW;
# per tree X=(ISORT/ITRN)·SCALE (DBH rank percentile), CRNEW=(A+B·(−ln(1−X))^(1/C))·10; CHG bounded ±1%/yr;
# CRMAX cap from HT+HTG. ISORT = whole-stand DBH rank (1=smallest…ITRN=largest) via _rdpsrt!.
const TT_WEIBA  = Float32[1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
const TT_WEIBB0 = Float32[-0.82631, -0.82631, -0.24217, 0.0, -0.90648, -0.08414, 0.17162, -0.90648, -0.89553, 0.24916, 0.0, 0.0, -0.2383, -0.08414, 0.0, -0.2383, -0.26595, 0.0]
const TT_WEIBB1 = Float32[1.06217, 1.06217, 0.96529, 0.0, 1.08122, 1.14765, 1.07338, 1.08122, 1.07728, 1.04831, 0.0, 0.0, 1.18016, 1.14765, 0.0, 1.18016, 0.98326, 0.0]
const TT_WEIBC0 = Float32[3.31429, 3.31429, -7.94832, 0.0, 3.48889, 2.775, 3.15, 3.48889, 1.74621, 4.36, 0.0, 0.0, 3.04, 2.775, 0.0, 3.04, -7.00555, 0.0]
const TT_WEIBC1 = Float32[0.0, 0.0, 1.93832, 0.0, 0.0, 0.0, 0.0, 0.0, 0.29052, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.60411, 0.0]
const TT_CRC0   = Float32[6.19911, 6.19911, 7.46296, 0.0, 6.81087, 4.01678, 6.00567, 6.81087, 7.65751, 6.41166, 0.0, 0.0, 4.62512, 4.01678, 0.0, 4.62512, 7.9281, 0.0]
const TT_CRC1   = Float32[-0.02216, -0.02216, -0.02944, 0.0, -0.01037, -0.01516, -0.0352, -0.01037, -0.03513, -0.02041, 0.0, 0.0, -0.01604, -0.01516, 0.0, -0.01604, -0.06298, 0.0]

@inline _tt_crown_diagr(sp::Int) = sp == 4 || sp == 11 || sp == 12 || sp == 15 || sp == 18

# ===================== small-tree crown dub (tt/dubscr.f) — DBH<1" at LSTART ======================
# crown.f:237 `IF(D.LT.1.0.AND.LSTART) GO TO 58` routes sub-1" inventory seedlings to label 58:
# sp15/18 (NC/OH) use a crown-length form; else CALL DUBSCR (logistic on BA/PCCF/AVH/RMAI). jl formerly
# applied the rank-Weibull to these seedlings (wrong). Coefficients verbatim from tt/dubscr.f DATA (no
# BCR4/BCR7 term). CR variation across trees is driven by the per-POINT PCCF (TPCCF=PCCF(point)).
const TT_BCR0  = Float32[-1.66949,-1.66949,-0.426688,-2.19723,-0.426688,-0.426688,-1.66949,-0.426688,-0.426688,-0.17561,-2.19723,-2.19723,5.0,-0.426688,-0.426688,5.0,-2.19723,-0.426688]
const TT_BCR1  = Float32[-0.209765,-0.209765,-0.093105,0.0,-0.093105,-0.093105,-0.209765,-0.093105,-0.093105,-0.33847,0.0,0.0,0.0,-0.093105,-0.093105,0.0,0.0,-0.093105]
const TT_BCR2  = Float32[0.0,0.0,0.022409,0.0,0.022409,0.022409,0.0,0.022409,0.022409,0.05699,0.0,0.0,0.0,0.022409,0.022409,0.0,0.0,0.022409]
const TT_BCR3  = Float32[0.003359,0.003359,0.002633,0.0,0.002633,0.002633,0.003359,0.002633,0.002633,0.00692,0.0,0.0,0.0,0.002633,0.002633,0.0,0.0,0.002633]
const TT_BCR5  = Float32[0.011032,0.011032,0.0,0.0,0.0,0.0,0.011032,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const TT_BCR6  = Float32[0.0,0.0,-0.045532,0.0,-0.045532,-0.045532,0.0,-0.045532,-0.045532,0.0,0.0,0.0,0.0,-0.045532,-0.045532,0.0,0.0,-0.045532]
const TT_BCR8  = Float32[0.017727,0.017727,0.0,0.0,0.0,0.0,0.017727,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const TT_BCR9  = Float32[-0.000053,-0.000053,0.000022,0.0,0.000022,0.000022,-0.000053,0.000022,0.000022,0.0,0.0,0.0,0.0,0.000022,0.000022,0.0,0.0,0.000022]
const TT_BCR10 = Float32[0.014098,0.014098,-0.013115,0.0,-0.013115,-0.013115,0.014098,-0.013115,-0.013115,0.0,0.0,0.0,0.0,-0.013115,-0.013115,0.0,0.0,-0.013115]
const TT_CRSD  = Float32[0.5,0.5,0.6957,0.2,0.6957,0.931,0.6124,0.6957,0.6957,0.8866,0.2,0.2,0.5,0.931,0.6957,0.5,0.2,0.6957]

@inline function _tt_dubscr(rng, sp::Integer, d::Real, h::Real, ba::Real, tpccf::Real, avh::Real, tmai::Real)::Float32
    hf = Float32(h); hf <= 0f0 && (hf = 0.1f0)
    cr = TT_BCR0[sp] + TT_BCR1[sp]*Float32(d) + TT_BCR2[sp]*hf + TT_BCR3[sp]*Float32(ba) +
         TT_BCR5[sp]*Float32(tpccf) + TT_BCR6[sp]*(Float32(avh)/hf) + TT_BCR8[sp]*Float32(avh) +
         TT_BCR9[sp]*(Float32(ba)*Float32(tpccf)) + TT_BCR10[sp]*Float32(tmai)
    sd = TT_CRSD[sp]
    fcr = 0f0
    while true                                          # dubscr.f label 10: FCR=BACHLO(0,SD); reject |FCR|>SD
        fcr = bachlo(rng, 0f0, sd)                      # DGSD=2.0≥1 ⇒ always draws
        abs(fcr) > sd && continue
        break
    end
    if sp == 13 || sp == 16                             # dubscr.f CASE(13,16): BI/MC linear form
        cr = ((cr - 1f0)*10f0 + 1f0)/100f0
    else
        abs(cr + fcr) >= 86f0 && (cr = 86f0)            # overflow guard (faithful: +86 regardless of sign)
        cr = 1f0/(1f0 + exp(cr + fcr))
    end
    cr < 0.05f0 && (cr = 0.05f0); cr > 0.95f0 && (cr = 0.95f0)
    return cr
end

# tt/maical.f RMAI = min(ADJMAI(ISPNUM(ISISP), SITEAR(ISISP), 10.0), 128). Stand constant; the
# BCR10*RMAI term matters for TT (WB/LM/LP BCR10≠0) — CI could pass 0 (its species have BCR10=0).
const TT_MAI_ISPNUM = Int32[101,101,202,101,93,746,108,93,19,122,101,101,746,746,101,746,101,101]
const _ADJMAI_ISP = Int32[19,41,42,71,81,94,95,242,263,264,98,298,
                          101,103,108,109,120,124, 119, 116,117,122, 201,202, 73, 92,93,
                          11,15,17,19,20,21,22, 211,212]
const _ADJMAI_IMAP = Int32[1,1,1,1,1,1,1,1,1,1,1,1, 2,2,2,2,2,2, 3, 4,4,4, 5,5, 6, 7,7,
                           9,9,9,9,9,9,9, 10,10]

function _adjmai(inspec::Integer, si::Float32, points::Float32)::Float32   # base/adjmai.f
    inum = 0
    if inspec >= 300
        inum = 8                                         # FIA ≥300 ⇒ hardwood group 8
    else
        for i in 1:32                                    # adjmai.f DO 2 I=1,32 (positions 33-36 unsearched)
            if _ADJMAI_ISP[i] == inspec; inum = _ADJMAI_IMAP[i]; break; end
        end
        inum == 0 && return 0f0                          # illegal site species → no MAI
    end
    a = 0f0
    if inum == 1
        si < 33f0 && return 0f0; a = -63.689706f0 + 1.9402941f0*si
    elseif inum == 2
        si < 11f0 && return 0f0; a = -12.0388f0 + 1.18672f0*si
    elseif inum == 3
        a = 5.972615f0 + 1.857675f0*si
    elseif inum == 4
        a = 2.305357f0 + 0.033890056f0*si + 0.0090108543f0*si^2
    elseif inum == 5
        si < 29f0 && return 0f0; a = -10.303313f0 + 0.032929911f0*si + 0.012207163f0*si^2 - 0.00003543129f0*si^3
    elseif inum == 6
        si < 11f0 && return 0f0; a = -6.0892857f0 + 0.45178571f0*si + 0.014464286f0*si^2
    elseif inum == 7
        si < 10f0 && return 0f0; a = -18.4f0 + 1.92f0*si
    elseif inum == 8
        si < 32f0 && return 0f0; a = -53.892857f0 + 1.7178571f0*si
    elseif inum == 9
        a = -4.89001f0 + 311.29546f0*((exp((si/170f0 - 1f0)^3/0.343f0) - 0.055f0)/0.95f0)
    elseif inum == 10
        si < 62f0 && return 0f0; a = 157.94643f0 - 1.78125f0*si + 0.014330357f0*si^2
    end
    a = a * points/10f0
    a < 0f0 && (a = 0f0)
    return a
end

function _tt_rmai(s::StandState)::Float32                # tt/maical.f
    p = s.plot
    isisp = Int(p.site_species); isisp == 0 && (isisp = 3)
    sssi = p.sp_site_index[isisp]; sssi == 0f0 && (sssi = 140f0)
    rmai = _adjmai(TT_MAI_ISPNUM[isisp], sssi, 10f0)
    rmai > 128f0 && (rmai = 128f0)
    return rmai
end

function crown_ratio_update!(s::StandState, ::Teton; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    relden = p.relative_density
    sdiac = crown_sdi
    # ISORT: whole-stand DBH rank via RDPSRT (1=smallest … n=largest). RDPSRT sorts DESCENDING (IND(1)=largest),
    # ISORT(IND(JJ)) = n−JJ+1  ⇒ largest→n, smallest→1. (once/cycle — local buffers, not the hot path.)
    key = Vector{Float32}(undef, n); idx = Vector{Int32}(undef, n)
    # rank on GROWN DBH (tt/crown.f runs after DG is applied: DBH(I) is post-growth). dbh += DG/bark.
    @inbounds for i in 1:n
        bk = bark_ratio(s.calib.bark_a, s.calib.bark_b, Int(t.species[i]), t.dbh[i])
        key[i] = t.dbh[i] + t.diam_growth[i] / bk
        idx[i] = Int32(i)
    end
    _rdpsrt!(key, idx; lseq = false)                    # descending: idx[1] = largest DBH
    isort = Vector{Int32}(undef, n)
    @inbounds for jj in 1:n; isort[idx[jj]] = Int32(n - jj + 1); end
    p_pccf = s.density.point_ccf
    rmai = lstart ? _tt_rmai(s) : 0f0                    # RMAI stand constant, only used by the lstart dub
    @inbounds for i in 1:n
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (lstart && t.crown_pct[i] > 0) && continue      # inventory crown present → keep (CROWN:222)
        # crown.f:237 — DBH<1" at LSTART routes to label 58 (small-tree dub), NOT the Weibull path.
        if d < 1f0 && lstart
            if sp == 15 || sp == 18                     # NC/OH crown-length form (label 58 CASE(15,18))
                hf = h + t.ht_growth[i]; hf <= 0f0 && (hf = 0.1f0)
                cl = 5.17281f0 + 0.32552f0*hf - 0.01675f0*p.basal_area
                cl < 1f0 && (cl = 1f0); cl > hf && (cl = hf)
                icri = trunc(Int, (cl/hf)*100f0 + 0.5f0)
            else                                        # label 58 CASE DEFAULT → DUBSCR
                pt = Int(t.plot_id[i])
                tpccf = (1 <= pt <= length(p_pccf)) ? p_pccf[pt] : 0f0
                cr = _tt_dubscr(s.rng, sp, d, h, p.basal_area, tpccf, p.avg_height, rmai)
                icri = trunc(Int, cr*100f0 + 0.5f0)
            end
            icri > 95 && (icri = 95); icri < 10 && (icri = 10); icri < 1 && (icri = 1)  # label 59 (CRNMLT=1)
            t.crown_pct[i] = Int32(icri)
            continue
        end
        _tt_crown_diagr(sp) && continue                 # PM/UJ/RM/NC/OH crown-length form (not ttt01) — skip
        icr = Int(t.crown_pct[i])
        # per-species Weibull params from mean crown ratio
        relsdi = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0
        relsdi > 1.5f0 && (relsdi = 1.5f0)
        (sp == 17 && relsdi > 1f0) && (relsdi = 1f0)
        acrnew = TT_CRC0[sp] + TT_CRC1[sp] * relsdi * 100f0
        A = TT_WEIBA[sp]
        B = TT_WEIBB0[sp] + TT_WEIBB1[sp] * acrnew
        B < (sp == 10 ? 3f0 : 1f0) && (B = sp == 10 ? 3f0 : 1f0)
        C = TT_WEIBC0[sp] + TT_WEIBC1[sp] * acrnew
        C < 2f0 && (C = 2f0)
        scale = 1f0 - 0.00167f0 * (relden - 100f0)
        scale > 1f0 && (scale = 1f0); scale < 0.30f0 && (scale = 0.30f0)
        x = d > 0f0 ? (Float32(isort[i]) / Float32(n)) * scale : 0.5f0 * scale   # d≤0 uses RANN (not in ttt01)
        x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
        crnew = (A + B * (-log(1f0 - x))^(1f0 / C)) * 10f0
        # change bounded ±1%/yr (skip when lstart or icr==0 → CRNEW stands)
        if !(lstart || icr == 0)
            chg = crnew - Float32(icr)
            pdifpy = chg / Float32(icr) / fint
            pdifpy > 0.01f0 && (chg = Float32(icr) * 0.01f0 * fint)
            pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
            crnew = Float32(icr) + chg
        end
        icri = trunc(Int, crnew + 0.5f0)
        # CRMAX cap (cycling only)
        if !(lstart || icr == 0)
            htg = t.ht_growth[i]
            crln = h * Float32(icr) / 100f0
            crmax = (crln + htg) / (h + htg) * 100f0
            (icri < 10 && Float32(icri) <= crmax) && (icri = trunc(Int, crmax + 0.5f0))  # CRNMLT=1
            Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
        end
        icri < 0 && (icri = 0); icri > 100 && (icri = 100)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end
