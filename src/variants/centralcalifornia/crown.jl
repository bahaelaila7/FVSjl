# =============================================================================
# crown.jl (centralcalifornia) — CA crown ratio (ca/crown.f + ca/dubscr.f). Chunk 5.
#
# Rank-based Weibull crown-ratio (Wykoff/Prognosis), coefficients indexed by the 17 crown groups
# (ca/crown.f IMAP[50]). Per species group:
#   ACRNEW = C0 + C1·RELSDI·100        (mean CR%, RELSDI = SDIAC/SDIDEF, capped 1.5)
#   A = WEIBA ; B = max(WEIBB0+WEIBB1·ACRNEW, 3) ; C = max(WEIBC0+WEIBC1·ACRNEW, 2)
# Per tree: X = (ISORT/ITRN)·SCALE, SCALE = clamp(1.5−RELSDI, 0.30, 1.0)  [CA-specific — WC uses a
#   RELDEN form]; CRNEW = (A + B·(−ln(1−X))^(1/C))·10 ; ±1%/yr change limit; CRMAX cap; [10,95].
# Species 23/50 (GS giant-sequoia, RW coast-redwood) = logistic on HDR/PRD/(D/QMDPLT): X=1/(1+exp(...));
#   CRNEW=X·10. PRD (point relative density) + point-QMDPLT feed ONLY the GS/RW branch (ca/crown.f SDICAL/
#   XMAXPT); no GS/RW in cat01 ⇒ the stand-SDI/stand-QMD proxies are inert there (point SDICAL deferred,
#   like WC). d<1" at LSTART → ca/dubscr.f (6 BCR groups + GS/RW logistic, group 7).
# =============================================================================

# ca/crown.f CRCONS 17-group coefficients (DATA WEIBA/WEIBB0/WEIBB1/WEIBC0/WEIBC1/C0/C1).
const CA_WEIBA  = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0]
const CA_WEIBB0 = Float32[0.52909,0.25115,0.52909,0.48464,0.08402,0.29964,0.06607,0.25667,0.16601,0.03685,0.25667,0.49085,0.16267,-0.81881,-1.11274,-0.23830,-0.13121]
const CA_WEIBB1 = Float32[1.00677,1.05987,1.00677,1.01272,1.10297,1.05398,1.10705,1.06474,1.08150,1.09499,1.06474,1.01414,1.07340,1.05418,1.12314,1.18016,1.15976]
const CA_WEIBC0 = Float32[-3.48211,0.33383,-3.48211,-2.78353,0.91078,-1.09270,2.04714,0.11729,0.91420,4.01340,0.11729,3.16456,3.28850,-2.36611,2.53316,3.04413,2.59824]
const CA_WEIBC1 = Float32[1.38780,0.63833,1.38780,1.27283,0.45819,0.80687,0.15070,0.61681,0.45768,0.04946,0.61681,0.00000,0.00000,1.20241,0.00000,0.00000,0.00000]
const CA_CRC0   = Float32[7.48846,6.92893,7.48846,7.44422,3.64292,5.12357,6.82187,5.95912,6.14578,6.04928,5.95912,5.48853,6.48494,4.42000,4.12048,4.62512,4.89032]
const CA_CRC1   = Float32[-0.02899,-0.04053,-0.02899,-0.04779,-0.00317,-0.01042,-0.02247,-0.01812,-0.02781,-0.01091,-0.01812,-0.00717,-0.02325,-0.01066,-0.00636,-0.01604,-0.01884]

# ca/crown.f DATA IMAP — 50 species → 17 crown groups.
const CA_CROWN_IMAP = Int[6,6,6,4,9,9,3,12,12,13,13,17,13,13,10,2,2,10,10,10,1,1,1,1,3,7,7,7,7,7,7,7,7,14,16,15,5,16,16,16,16,8,16,16,16,16,16,16,16,1]

# ------------------- ca/dubscr.f — small-tree (d<1") / dead-tree crown dub -------------------
# IMAP[50]→7 groups. 1=firs/spruces, 2=DF, 3=cedars/hemlocks, 4=pines/PY, 5=hardwoods(const ICR=5),
# 6=juniper(const ICR=9), 7=GS/RW logistic. CR code 0-9 from BCR0+BCR1·H+BCR2·BA; GS(23)/RW(50) logistic.
const CA_DUB_IMAP = Int32[3,3,3,1,1,1,2,3,3,4,4,4,4,4,4,4,4,4,4,4,6,1,7,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,7]
const CA_DUB_BCR0 = Float32[8.042774,8.477025,7.558538,6.489813,5.000000,9.000000,0.000000]
const CA_DUB_BCR1 = Float32[0.007198,-0.018033,-0.015637,-0.029815,0.000000,0.000000,0.000000]
const CA_DUB_BCR2 = Float32[-0.016163,-0.018140,-0.009064,-0.009276,0.000000,0.000000,0.000000]
const CA_DUB_CRSD = Float32[1.3167,1.3756,1.9658,2.0426,0.5,0.5,0.15]

@inline function _ca_dubscr(rng, sp::Integer, d::Real, h::Real, ba::Real, prd::Real, qmdplt::Real)::Float32
    g = Int(CA_DUB_IMAP[sp])
    if sp == 23 || sp == 50                        # GS/RW logistic (ca/dubscr.f)
        hdr = Float32(h) * 12f0 / Float32(d)
        cr = -1.021064f0 + 0.309296f0 * log(max(hdr, 1f-6)) + 0.869720f0 * Float32(prd) -
             0.116274f0 * (Float32(d) / Float32(qmdplt))
        sd = CA_DUB_CRSD[g]
        fcr = 0f0
        while true
            fcr = bachlo(rng, 0f0, sd)
            abs(fcr) > sd && continue
            break
        end
        cr = 1f0 / (1f0 + exp(cr + fcr))
    else
        cr = CA_DUB_BCR0[g] + CA_DUB_BCR1[g] * Float32(h) + CA_DUB_BCR2[g] * Float32(ba)
        sd = CA_DUB_CRSD[g]
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

# ca/ccfcal.f MODE=1 — per-tree CCF (÷TPA) = 0.001803·CRWD5², CRWD5 = ca/r5crwd.f R5CRWD open-grown crown
# width. CA's R5CRWD (MAPCA species map + WB1/WB2/WB3/DX1/DX2/IEQN/SPLINE/SM) is BYTE-IDENTICAL to OC's copy
# (both the shared vws/r5crwd.f) — verified MAPCA 27,7,16,4,… == OC_R5CRWD_MAPCA ⇒ reuse (doctrine #5).
@inline ca_tree_ccf(sp::Int, d::Float32, h::Float32) = oc_tree_ccf(sp, d, h)

function crown_ratio_update!(s::StandState, ::CentralCalifornia; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    sd = s.coef.species
    relden = p.relative_density
    sdiac = crown_sdi
    ba = p.basal_area
    qmd = stand_qmd(s)
    # ISORT: whole-stand DBH rank on the GROWN diameter (crown.f runs after DG). RDPSRT sorts
    # descending (idx[1]=largest); ISORT(idx[jj]) = n−jj+1 ⇒ largest→n, smallest→1.
    key = Vector{Float32}(undef, n); idx = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        isp = Int(t.species[i])
        bk = wc_bratio(sd[:bark1][isp], sd[:bark2][isp], Int(sd[:bark_imap][isp]), t.dbh[i])
        key[i] = t.dbh[i] + t.diam_growth[i] / bk; idx[i] = Int32(i)
    end
    _rdpsrt!(key, idx; lseq = false)
    isort = Vector{Int32}(undef, n)
    @inbounds for jj in 1:n; isort[idx[jj]] = Int32(n - jj + 1); end
    @inbounds for i in 1:n
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (sp < 1 || sp > 50) && continue
        (lstart && t.crown_pct[i] > 0) && continue        # crown.f:227 keep inventory crown
        icr = Int(t.crown_pct[i])
        # GS/RW plot-level QMD / relative density (only used by the GS/RW branch + dubscr). Point SDICAL
        # deferred (no GS/RW in cat01) ⇒ stand-SDI/stand-QMD proxies, source-faithful for the exercised path.
        prd = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0
        qmdplt = qmd > 1f0 ? qmd : 1f0
        # crown.f:271 — d<1" at LSTART routes to statement 58 (DUBSCR), not the Weibull path.
        if d < 1f0 && lstart
            icr != 0 && continue                          # crown.f:379 IF(ICR.NE.0) GO TO 60
            cr = _ca_dubscr(s.rng, sp, d, h, ba, prd, qmdplt)
            icri = trunc(Int, cr * 100f0 + 0.5f0)
            icri > 95 && (icri = 95); icri < 10 && (icri = 10); icri < 1 && (icri = 1)
            t.crown_pct[i] = Int32(icri)
            continue
        end
        grp = CA_CROWN_IMAP[sp]
        relsdi = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0
        relsdi > 1.5f0 && (relsdi = 1.5f0)
        acrnew = CA_CRC0[grp] + CA_CRC1[grp] * relsdi * 100f0
        local crnew::Float32
        if sp == 23 || sp == 50                           # GS/RW logistic (ca/crown.f CASE(23,50))
            hdr = d > 0f0 ? h * 12f0 / d : 1f0
            xl = -1.021064f0 + 0.309296f0 * log(max(hdr, 1f-6)) + 0.869720f0 * prd - 0.116274f0 * (d / qmdplt)
            x = 1f0 / (1f0 + exp(xl))
            x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
            crnew = x * 10f0                               # CASE(23,50): CRNEW = X*10
        else
            A = CA_WEIBA[grp]
            B = CA_WEIBB0[grp] + CA_WEIBB1[grp] * acrnew
            B < 3f0 && (B = 3f0)
            C = CA_WEIBC0[grp] + CA_WEIBC1[grp] * acrnew
            C < 2f0 && (C = 2f0)
            scale = 1.5f0 - relsdi                         # ca/crown.f:292 (CA-specific)
            scale > 1f0 && (scale = 1f0); scale < 0.30f0 && (scale = 0.30f0)
            x = d > 0f0 ? (Float32(isort[i]) / Float32(n)) * scale : 0.5f0 * scale
            x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
            crnew = A + B * (-log(1f0 - x))^(1f0 / C)
        end
        crnew *= 10f0
        # ±1%/yr change limit (skip when lstart or icr==0); CRNMLT=1 ⇒ no band multiplier.
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
            (icri < 10) && (icri = trunc(Int, crmax + 0.5f0))          # CRNMLT=1
            Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
        end
        # topkill (LSTART & ITRUNC≠0) — no truncation state on cat01 trees (inert), source-faithful.
        if lstart && t.trunc[i] != 0
            hn = Float32(t.norm_ht[i]) / 100f0
            hd = hn - Float32(t.trunc[i]) / 100f0
            cl = (Float32(icri) / 100f0) * hn - hd
            icri = trunc(Int, (cl * 100f0 / hn) + 0.5f0)
        end
        icri > 95 && (icri = 95); icri < 10 && (icri = 10); icri < 1 && (icri = 1)  # statement 59 (CRNMLT=1)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end

# ---------------------------------------------------------------------------
# CA FFE crown-biomass species map (ca/fmcrow.f:109 DATA ISPMAP) — the Jenkins/FMCROWE
# crown-biomass group per species. ca/fmcrow.f:167-169 routes CASE(35,39,40,41,43,44,45,46)
# → FMCROWE (eastern Jenkins TOTABV), all others → FMCROWW (western, shared cr_crownw).
# ca/fmcroww.f is byte-identical (md5 f764dce1) to CR/WC/WS's, so CA reuses cr_crownw.
# Fixes #229 (the CA crown-fire under-kill): without this, CentralCalifornia fell to the
# Jenkins FMCROWE path (via ls_spi) ⇒ crown biomass 6-20× too low ⇒ canopy_bulk_density
# actcbh=12/cbd=0.049 vs oracle 4/0.129 ⇒ surface fire instead of the oracle's passive crown.
# (Consumed by crown_biomass.jl once CentralCalifornia is wired into the cr_crownw dispatch.)
# ---------------------------------------------------------------------------
const CA_ISPMAP = Int[
   7, 20,  7,  4,  4,  4,  3,  6, 24, 14,
  11, 11, 11, 11, 15, 15, 15, 13, 13, 11,
  16, 18, 19,  7, 11, 17, 17, 21, 17, 21,
  21, 21, 17,  5, 44, 23, 10, 17, 56, 29,
  46, 17, 60, 41, 17, 64, 17, 17, 21, 19]
@inline ca_uses_fmcrowe(sp::Integer) =
    (sp == 35 || sp == 39 || sp == 40 || sp == 41 || sp == 43 || sp == 44 || sp == 45 || sp == 46)
