# =============================================================================
# crown.jl (olympic) — OP crown ratio (op/crown.f + op/dubscr.f), COOPERATING with ORGANON.
#
# Rank-based Weibull crown-ratio (Wykoff/Prognosis) for the FVS-native (IORG=0) trees; the valid-
# ORGANON (IORG=1) trees take CRNEW = ANINT(CR2·100) from the stashed ORGANON crown (op/crown.f:279-286).
# Coefficients indexed by the 17 crown groups (op/crown.f IMAP[39]); crown-group order:
#   1=SF 2=WF 3=RF 4=NF 5=WP 6=PP 7=DF 8=RC 9=WH 10=MH 11=RW 12=BM 13=RA 14=OH 15=IC 16=LP 17=SS.
# Per FVS-native tree:  X = (ISORT/ITRN)·SCALE, SCALE = clamp(1 − 0.00167·(RELDEN−100), 0.30, 1.0)
#   [op/crown.f:306 — the WC RELDEN form]; CRNEW = (A + B·(−ln(1−X))^(1/C))·10 ; ±1%/yr limit; CRMAX
#   cap; [1,95]. RW (sp17) uses the logistic HDR/PRD/(D/QMDPLT) form (no RW in S248112). d<1" at LSTART
#   → op/dubscr.f (WC 6-BCR groups + RW logistic).
# =============================================================================

# op/crown.f DATA IMAP — 39 species → 17 crown groups.
const OP_CROWN_IMAP = Int[
    1,2,2,3,3,17,4,15,11,11,16,6,5,5,6,7,11,8,9,10,12,13,14,14,14,14,14,14,14,11,11,11,11,14,14,14,14,14,14]

# op/crown.f CRCONS 17-group coefficients (DATA WEIBA/WEIBB0/WEIBB1/WEIBC0/WEIBC1/C0/C1).
const OP_WEIBA  = Float32[0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0]
const OP_WEIBB0 = Float32[-0.171680,0.130939,-0.981113,-0.135807,0.019948,-0.036696,-0.012061,-0.062693,0.073435,0.162672,0.196054,-0.818809,0.035786,-0.238295,-0.811424,-0.131210,-0.107413]
const OP_WEIBB1 = Float32[1.161549,1.093406,1.092273,1.147712,1.108738,1.132792,1.119712,1.139657,1.107183,1.073404,1.073909,1.054176,1.121389,1.180163,1.056190,1.159760,1.140775]
const OP_WEIBC0 = Float32[2.8263,1.355139,1.326047,3.017494,2.621230,2.876094,3.2126,1.7664,2.6237,3.288501,0.345647,-2.366108,2.0408,3.044134,-3.831124,2.598238,3.0712]
const OP_WEIBC1 = Float32[0.0,0.350472,0.318386,0.0,0.186734,0.0,0.0,0.0,0.0,0.0,0.620145,1.202413,0.0,0.0,1.401938,0.0,0.0]
const OP_CRC0   = Float32[5.073342,5.212394,4.860467,5.568864,4.279655,5.073273,5.666442,4.481330,5.671345,6.484942,5.417431,4.420000,4.656659,4.625125,5.200550,4.890318,5.812879]
const OP_CRC1   = Float32[-0.011430,-0.011623,-0.006173,-0.021293,-0.002484,-0.020988,-0.025199,-0.018092,-0.023463,-0.023248,-0.011608,-0.010660,-0.022612,-0.016042,-0.014890,-0.018837,-0.028504]

# op/dubscr.f (WC copy) — 39 species → 6 BCR groups; RW (sp17) uses the logistic branch instead.
const OP_DUB_IMAP = Int32[1,1,1,1,1,1,1,3,3,1,4,4,4,4,4,2,2,3,3,3,5,5,5,5,5,5,5,5,6,4,4,4,4,5,5,5,5,5,5]
const OP_DUB_BCR0 = Float32[8.042774,8.477025,7.558538,6.489813,5.0,9.0]
const OP_DUB_BCR1 = Float32[0.007198,-0.018033,-0.015637,-0.029815,0.0,0.0]
const OP_DUB_BCR2 = Float32[-0.016163,-0.018140,-0.009064,-0.009276,0.0,0.0]
const OP_DUB_CRSD = Float32[1.3167,1.3756,1.9658,2.0426,0.5,0.5]

"op/dubscr.f — d<1\" (or dead) missing-crown dub. Draws BACHLO(0,SD,RANN) rejecting |FCR|>SD."
@inline function _op_dubscr(rng, sp::Integer, d::Real, h::Real, ba::Real, prd::Real, qmdplt::Real)::Float32
    if sp == 17                                            # RW — logistic (op/dubscr.f:21-24)
        hdr = Float32(h) * 12f0 / Float32(d)
        cr = -1.021064f0 + 0.309296f0 * log(max(hdr, 1f-6)) + 0.869720f0 * Float32(prd) -
             0.116274f0 * (Float32(d) / Float32(qmdplt))
        sd = 0.15f0
        fcr = 0f0
        while true
            fcr = bachlo(rng, 0f0, sd); abs(fcr) > sd && continue; break
        end
        cr = 1f0 / (1f0 + exp(cr + fcr))
    else
        g = Int(OP_DUB_IMAP[sp])
        cr = OP_DUB_BCR0[g] + OP_DUB_BCR1[g] * Float32(h) + OP_DUB_BCR2[g] * Float32(ba)
        sd = OP_DUB_CRSD[g]
        fcr = 0f0
        while true
            fcr = bachlo(rng, 0f0, sd); abs(fcr) > sd && continue; break
        end
        cr = ((cr + fcr) - 1f0) * 10f0 / 100f0 + 1f0 / 100f0    # ((CR-1)*10+1)/100
    end
    cr > 0.95f0 && (cr = 0.95f0); cr < 0.05f0 && (cr = 0.05f0)
    return cr
end

function crown_ratio_update!(s::StandState, ::Olympic; fint::Float32 = 5.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t, c = s.plot, s.trees, s.calib
    n = t.n; n == 0 && return s
    relden = p.relative_density
    sdiac = crown_sdi
    ba = p.basal_area
    qmd = stand_qmd(s)
    org_ran = length(c.op_iorg) == n
    # ISORT: whole-stand DBH rank on the GROWN diameter (op/crown.f runs after DG; DBH+DG/BRATIO).
    # RDPSRT sorts descending (idx[1]=largest); ISORT(idx[jj]) = n−jj+1 ⇒ largest→n, smallest→1.
    key = Vector{Float32}(undef, n); idx = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        isp = Int(t.species[i])
        bk = op_bratio(isp, t.dbh[i])
        key[i] = t.dbh[i] + t.diam_growth[i] / bk; idx[i] = Int32(i)
    end
    _rdpsrt!(key, idx; lseq = false)
    isort = Vector{Int32}(undef, n)
    @inbounds for jj in 1:n; isort[idx[jj]] = Int32(n - jj + 1); end
    @inbounds for i in 1:n
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (sp < 1 || sp > 39) && continue
        (lstart && t.crown_pct[i] > 0) && continue          # op/crown.f:226 keep inventory crown
        icr = Int(t.crown_pct[i])
        iorg = org_ran && c.op_iorg[i] == 1
        prd = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0     # inert unless RW present
        qmdplt = qmd > 1f0 ? qmd : 1f0
        # op/crown.f:275 — d<1" at LSTART routes to statement 58 (DUBSCR).
        if d < 1f0 && lstart
            icr != 0 && continue                            # op/crown.f:407 IF(ICR.NE.0) GO TO 60
            cr = _op_dubscr(s.rng, sp, d, h, ba, prd, qmdplt)
            icri = trunc(Int, cr * 100f0 + 0.5f0)
            icri > 95 && (icri = 95); icri < 10 && (icri = 10); icri < 1 && (icri = 1)
            t.crown_pct[i] = Int32(icri)
            continue
        end
        grp = OP_CROWN_IMAP[sp]
        relsdi = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0
        relsdi > 1.5f0 && (relsdi = 1.5f0)
        acrnew = OP_CRC0[grp] + OP_CRC1[grp] * relsdi * 100f0
        local crnew::Float32
        local chg::Float32 = 0f0
        if iorg                                             # ORGANON crown (op/crown.f:279-286)
            crnew = round(c.op_cr2[i] * 100f0, RoundNearestTiesAway)   # ANINT(CR2*100)
            chg = crnew - Float32(icr)
        elseif sp == 17                                     # RW logistic (op/crown.f CASE(17))
            hdr = d > 0f0 ? h * 12f0 / d : 1f0
            xl = -1.021064f0 + 0.309296f0 * log(max(hdr, 1f-6)) + 0.869720f0 * prd - 0.116274f0 * (d / qmdplt)
            x = 1f0 / (1f0 + exp(xl))
            x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
            crnew = x * 10f0                                # CASE(17): CRNEW = X*10
            crnew *= 10f0
            if !(lstart || icr == 0)
                chg = crnew - Float32(icr)
                pdifpy = chg / Float32(icr) / fint
                pdifpy > 0.01f0 && (chg = Float32(icr) * 0.01f0 * fint)
                pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
            end
        else                                                # DEFAULT Weibull (op/crown.f:305-314)
            A = OP_WEIBA[grp]
            B = OP_WEIBB0[grp] + OP_WEIBB1[grp] * acrnew; B < 3f0 && (B = 3f0)
            C = OP_WEIBC0[grp] + OP_WEIBC1[grp] * acrnew; C < 2f0 && (C = 2f0)
            scale = 1f0 - 0.00167f0 * (relden - 100f0)      # op/crown.f:306 (WC RELDEN form)
            scale > 1f0 && (scale = 1f0); scale < 0.30f0 && (scale = 0.30f0)
            x = d > 0f0 ? (Float32(isort[i]) / Float32(n)) * scale : 0.5f0 * scale   # RANN path (d=0) inert here
            x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
            crnew = A + B * (-log(1f0 - x))^(1f0 / C)
            crnew *= 10f0
            if !(lstart || icr == 0)
                chg = crnew - Float32(icr)
                pdifpy = chg / Float32(icr) / fint
                pdifpy > 0.01f0 && (chg = Float32(icr) * 0.01f0 * fint)
                pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
            end
        end
        # op/crown.f:347-352 statement 41 — CRNEW = ICR + CHG (CRNMLT=1, DLOW/DHI cover all d).
        if !(lstart || icr == 0) || iorg
            crnew = Float32(icr) + chg
        end
        icri = trunc(Int, crnew + 0.5f0)                    # 9052
        # CRMAX cap (op/crown.f:365-378, cycling with ICR>0).
        if !(lstart || icr == 0)
            htg = t.ht_growth[i]
            crln = h * Float32(icr) / 100f0
            crmax = (crln + htg) / (h + htg) * 100f0
            # op/crown.f:383-385 — the ICRI<10 bump is SKIPPED for ORGANON trees.
            (!iorg && icri < 10) && (icri = trunc(Int, crmax + 0.5f0))
            Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
        end
        # topkill (LSTART & ITRUNC≠0) — inert on S248112 (no truncation).
        if lstart && t.trunc[i] != 0
            hn = Float32(t.norm_ht[i]) / 100f0
            hd = hn - Float32(t.trunc[i]) / 100f0
            cl = (Float32(icri) / 100f0) * hn - hd
            icri = trunc(Int, (cl * 100f0 / hn) + 0.5f0)
        end
        icri > 95 && (icri = 95)                            # statement 59
        # op/crown.f:426-428 — the ICRI<10→10 bump is SKIPPED for ORGANON trees.
        (!iorg && icri < 10) && (icri = 10)
        icri < 1 && (icri = 1)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end
