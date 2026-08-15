# =============================================================================
# crown.jl (westcascades) — WC crown ratio (wc/crown.f + wc/dubscr.f). Chunk 5.
#
# Rank-based Weibull crown-ratio (Wykoff/Prognosis), coefficients indexed by the 16 crown
# groups (crown_imap, wc/crown.f IMAP[39]). Per species group:
#   ACRNEW = C0 + C1·RELSDI·100        (mean CR%, RELSDI = SDIAC/SDIDEF, capped 1.5)
#   A = WEIBA ; B = max(WEIBB0+WEIBB1·ACRNEW, 3) ; C = max(WEIBC0+WEIBC1·ACRNEW, 2)
# Per tree: X = (ISORT/ITRN)·SCALE, SCALE = clamp(1−0.00167·(RELDEN−100), 0.30, 1.0);
#   CRNEW = (A + B·(−ln(1−X))^(1/C))·10 ; ±1%/yr change limit vs old ICR; CRMAX cap; [10,95].
# Species 17 (RW, redwood) = logistic on HDR/PRD/(D/QMD): X=1/(1+exp(...)); CRNEW=X·100. No RW in
# wct01, ported source-faithful (June-2021 RW edit). d<1" at LSTART → wc/dubscr.f (6 BCR groups + RW).
# =============================================================================

# CRCONS 16-group coefficients (wc/crown.f DATA WEIBA/WEIBB0/WEIBB1/WEIBC0/WEIBC1/C0/C1).
# Group order: 1=PSF 2=WF 3=RF 4=NF 5=WWP 6=PP 7=DF 8=WRC 9=WH 10=MH 11=RW 12=MAP 13=ALD 14=OH 15=IC 16=LP
const WC_WEIBA  = Float32[0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0]
const WC_WEIBB0 = Float32[-0.173100, 0.130939, -0.981113, -0.135807, 0.019948, -0.036696, -0.082379, 0.179839, 0.490848, 0.162672, 0.196054, -0.818809, -1.112738, -0.238295, -0.811424, -0.131210]
const WC_WEIBB1 = Float32[1.080573, 1.093406, 1.092273, 1.147712, 1.108738, 1.132792, 1.137459, 1.084924, 1.014138, 1.073404, 1.073909, 1.054176, 1.123138, 1.180163, 1.056190, 1.159760]
const WC_WEIBC0 = Float32[1.062168, 1.355139, 1.326047, 3.017494, 2.621230, 2.876094, 2.914892, 0.122967, 3.164558, 3.288501, 0.345647, -2.366108, 2.533158, 3.044134, -3.831124, 2.598238]
const WC_WEIBC1 = Float32[0.445799, 0.350472, 0.318386, 0.000000, 0.186734, 0.000000, 0.000000, 0.567784, 0.000000, 0.000000, 0.620145, 1.202413, 0.000000, 0.000000, 1.401938, 0.000000]
const WC_CRC0   = Float32[5.614200, 5.212394, 4.860467, 5.568864, 4.279655, 5.073273, 5.067560, 5.570928, 5.488532, 6.484942, 5.417431, 4.420000, 4.120478, 4.625125, 5.200550, 4.890318]
const WC_CRC1   = Float32[-0.016547, -0.011623, -0.006173, -0.021293, -0.002484, -0.020988, -0.010484, -0.012043, -0.007173, -0.023248, -0.011608, -0.010660, -0.006357, -0.016042, -0.014890, -0.018837]

# ------------------- wc/dubscr.f — small-tree (d<1") / dead-tree crown dub -------------------
# IMAP2[39]→6 BCR groups. Group 1=SF..ES, 2=DF, 3=YC/IC/RC/WH/MH, 4=LP..PY, 5=hardwoods(const ICR=5),
# 6=J(const ICR=9). CR code 0-9 from BCR0+BCR1·H+BCR2·BA; RW (sp17) logistic on HDR/PRD/(D/QMD).
const WC_DUB_IMAP = Int32[1,1,1,1,1,1,1,3,3,1,4,4,4,4,4,2,2,3,3,3,5,5,5,5,5,5,5,5,6,4,4,4,4,5,5,5,5,5,5]
const WC_DUB_BCR0 = Float32[8.042774, 8.477025, 7.558538, 6.489813, 5.000000, 9.000000]
const WC_DUB_BCR1 = Float32[0.007198, -0.018033, -0.015637, -0.029815, 0.000000, 0.000000]
const WC_DUB_BCR2 = Float32[-0.016163, -0.018140, -0.009064, -0.009276, 0.000000, 0.000000]
const WC_DUB_CRSD = Float32[1.3167, 1.3756, 1.9658, 2.0426, 0.5, 0.5]

@inline function _wc_dubscr(rng, sp::Integer, d::Real, h::Real, ba::Real, prd::Real, qmdplt::Real)::Float32
    g = Int(WC_DUB_IMAP[sp])
    if sp == 17                                    # redwood logistic (wc/dubscr.f RW branch)
        hdr = Float32(h) * 12f0 / Float32(d)
        cr = -1.021064f0 + 0.309296f0 * log(max(hdr, 1f-6)) + 0.869720f0 * Float32(prd) -
             0.116274f0 * (Float32(d) / Float32(qmdplt))
        sd = 0.15f0
        fcr = 0f0
        while true
            fcr = bachlo(rng, 0f0, sd)
            abs(fcr) > sd && continue
            break
        end
        cr = 1f0 / (1f0 + exp(cr + fcr))
    else
        cr = WC_DUB_BCR0[g] + WC_DUB_BCR1[g] * Float32(h) + WC_DUB_BCR2[g] * Float32(ba)
        sd = WC_DUB_CRSD[g]
        fcr = 0f0
        while true
            fcr = bachlo(rng, 0f0, sd)
            abs(fcr) > sd && continue
            break
        end
        cr = ((cr + fcr) - 1f0) * 10f0 / 100f0 + 1f0 / 100f0   # ((CR-1)*10+1)/100
    end
    cr > 0.95f0 && (cr = 0.95f0); cr < 0.05f0 && (cr = 0.05f0)
    return cr
end

function crown_ratio_update!(s::StandState, ::WestCascades; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    sd = s.coef.species
    cimap = sd[:crown_imap]
    relden = p.relative_density
    sdiac = crown_sdi
    ba = p.basal_area
    qmd = stand_qmd(s)
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    # ISORT: whole-stand DBH rank on the GROWN diameter (crown.f runs after DG). RDPSRT sorts
    # descending (idx[1]=largest); ISORT(idx[jj]) = n−jj+1 ⇒ largest→n, smallest→1.
    key = Vector{Float32}(undef, n); idx = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        bk = wc_bratio(sd, Int(t.species[i]), t.dbh[i])
        key[i] = t.dbh[i] + t.diam_growth[i] / bk; idx[i] = Int32(i)
    end
    _rdpsrt!(key, idx; lseq = false)
    isort = Vector{Int32}(undef, n)
    @inbounds for jj in 1:n; isort[idx[jj]] = Int32(n - jj + 1); end
    @inbounds for i in 1:n
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (sp < 1 || sp > 39) && continue
        (lstart && t.crown_pct[i] > 0) && continue        # crown.f:229 keep inventory crown
        icr = Int(t.crown_pct[i])
        # RW plot-level QMD / relative density (only used by the RW branch + dubscr; cheap to always set)
        prd = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0   # PRD≈RELSDI proxy (point SDI not ported)
        qmdplt = qmd > 1f0 ? qmd : 1f0
        # crown.f:278 — d<1" at LSTART routes to label 58 (DUBSCR), not the Weibull path.
        if d < 1f0 && lstart
            icr != 0 && continue                          # crown.f:382 IF(ICR.NE.0) GO TO 60
            cr = _wc_dubscr(s.rng, sp, d, h, ba, prd, qmdplt)
            icri = trunc(Int, cr * 100f0 + 0.5f0)
            icri > 95 && (icri = 95); icri < 10 && (icri = 10); icri < 1 && (icri = 1)
            t.crown_pct[i] = Int32(icri)
            continue
        end
        grp = Int(cimap[sp])
        relsdi = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0
        relsdi > 1.5f0 && (relsdi = 1.5f0)
        acrnew = WC_CRC0[grp] + WC_CRC1[grp] * relsdi * 100f0
        local crnew::Float32
        if sp == 17                                       # redwood logistic (wc/crown.f CASE(17))
            hdr = d > 0f0 ? h * 12f0 / d : 1f0
            xl = -1.021064f0 + 0.309296f0 * log(max(hdr, 1f-6)) + 0.869720f0 * prd - 0.116274f0 * (d / qmdplt)
            x = 1f0 / (1f0 + exp(xl))
            x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
            crnew = x * 10f0                               # CASE(17): CRNEW = X*10
        else
            A = WC_WEIBA[grp]
            B = WC_WEIBB0[grp] + WC_WEIBB1[grp] * acrnew
            B < 3f0 && (B = 3f0)
            C = WC_WEIBC0[grp] + WC_WEIBC1[grp] * acrnew
            C < 2f0 && (C = 2f0)
            scale = 1f0 - 0.00167f0 * (relden - 100f0)
            scale > 1f0 && (scale = 1f0); scale < 0.30f0 && (scale = 0.30f0)
            x = d > 0f0 ? (Float32(isort[i]) / Float32(n)) * scale : 0.5f0 * scale   # d≤0 uses RANN (not in wct01)
            x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
            crnew = A + B * (-log(1f0 - x))^(1f0 / C)
        end
        crnew *= 10f0
        # ±1%/yr change limit (skip when lstart or icr==0)
        if !(lstart || icr == 0)
            chg = crnew - Float32(icr)
            pdifpy = chg / Float32(icr) / fint
            pdifpy > 0.01f0 && (chg = Float32(icr) * 0.01f0 * fint)
            pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
            crnew = Float32(icr) + chg                     # CRNMLT=1 ⇒ no band multiplier
        end
        icri = trunc(Int, crnew + 0.5f0)
        # CRMAX cap (cycling only)
        if !(lstart || icr == 0)
            htg = t.ht_growth[i]
            crln = h * Float32(icr) / 100f0
            crmax = (crln + htg) / (h + htg) * 100f0
            (icri < 10) && (icri = trunc(Int, crmax + 0.5f0))          # CRNMLT=1
            Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
        end
        # topkill (LSTART & ITRUNC≠0) — no truncation state on wct01 trees (inert), source-faithful.
        if lstart && t.trunc[i] != 0
            hn = Float32(t.norm_ht[i]) / 100f0
            hd = hn - Float32(t.trunc[i]) / 100f0
            cl = (Float32(icri) / 100f0) * hn - hd
            icri = trunc(Int, (cl * 100f0 / hn) + 0.5f0)
        end
        icri > 95 && (icri = 95); icri < 10 && (icri = 10); icri < 1 && (icri = 1)  # label 59 (CRNMLT=1)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end

# wc/ccfcal.f MODE=1 — per-tree crown competition factor (CCFT, before ×P). 16 CCF groups (INDCCF).
# stand_ccf sums CCFT·P = RELDEN; point_density sums per point (PCCF). D<1 uses the linear-extrapolated form.
const WC_CCF_INDCCF = Int[4,4,4,13,9,9,13,6,6,14,15,10,1,1,10,3,3,6,12,12,8,2,2,8,16,8,8,8,6,6,15,15,8,8,8,8,8,8,8]
const WC_CCF_RD1 = Float32[0.0392,0.03561,0.0388,0.0690,0.0392,0.0194,0.0212,0.0204,0.0172,0.0219,0.0388,0.03758,0.02453,0.03,0.01925,0.0160]
const WC_CCF_RD2 = Float32[0.0180,0.02731,0.0269,0.0225,0.0180,0.0142,0.0167,0.0246,0.00876,0.01676,0.0269,0.0233,0.0115,0.0173,0.0168,0.0167]
const WC_CCF_RD3 = Float32[0.00207,0.00524,0.00466,0.00183,0.00207,0.00261,0.00330,0.0074,0.00112,0.00325,0.00466,0.00361,0.00134,0.00259,0.00365,0.00434]

@inline function wc_tree_ccf(sp::Integer, d::Real)::Float32
    (sp < 1 || sp > 39) && return 0f0
    ic = WC_CCF_INDCCF[sp]; D = Float32(d)
    return D < 1.0f0 ? D * (WC_CCF_RD1[ic] + WC_CCF_RD2[ic] + WC_CCF_RD3[ic]) :
                       WC_CCF_RD1[ic] + WC_CCF_RD2[ic] * D + WC_CCF_RD3[ic] * D * D
end

# ---------------------------------------------------------------------------
# FFE crown-biomass group per species (wc/fmcrow.f:103 DATA ISPMAP). wc/fmcrow.f:157-162 routes
# CASE(24,26,27,34,35,36,37,39)=PB/AS/CW/DG/HT/CH/WI/OT → FMCROWE (eastern Jenkins TOTABV), all others
# → FMCROWW (western crown-width, shared cr_crownw). WC FFE chunk F1.
# ---------------------------------------------------------------------------
const WC_ISPMAP = Int[
   4,  4,  4,  1,  4,  0,  4,  8, 20, 18,
  11, 15, 15, 15, 13,  3, 19,  7,  6, 24,
   5, 23, 10, 43, 17, 41, 17, 17, 16,  1,
  14, 11,  7, 56, 57, 61, 64,  0, 41]
@inline wc_uses_fmcrowe(sp::Integer) = (sp == 24 || sp == 26 || sp == 27 || sp == 34 ||
                                        sp == 35 || sp == 36 || sp == 37 || sp == 39)
