# =============================================================================
# regent.jl (britishcolumbia) — BC small-tree growth engine hook (canada/bc/regent.f, V3).
#
# CHUNK 6. `small_tree_growth!(s, stash, ::BritishColumbia; fint)` overrides DG/HTG for small trees,
# mirroring KT (src/variants/kootenai/regent.jl:289-411) with BC swaps:
#   • height increment = bc_v3_sthg (VALIDATED regression) accumulated over subcycles;
#   • ZZRAN (V3): HTGR = max(0, HTGR1 + ZZRAN·ST_COEF%SD)  [cornered class, seed 55329];
#   • DBH-from-height = bc_st_dbh (POWER form HHT1·(H-4.5)^HHT2 + DADJ);
#   • per-subcycle density via bc_tree_ccf (chunk-5 CCF spine).
# ⚠ V2 regime + SBS logistic branch TODO. Bounds XMINV3/XMAXV3: sp14=2/4 confirmed, others fill(2/4) TODO.
# =============================================================================

"""Small-tree site setup (regent MORCON-analog): resolve ST_COEF ip + RHCON per species (V3)."""
function bc_regcons!(s::StandState)
    ctl = s.control
    zone, series = bc_stand_zone(s)
    nsp = nspecies(BritishColumbia())
    ip = zeros(Int, nsp); rhcon = zeros(Float32, nsp)
    @inbounds for sp in 1:nsp
        ip[sp] = bc_resolve_stcoef(sp, series, zone)
        # V3 RHCON (regent.f:2112): 0 base, ln(RCOR2) when small-tree HT calibrated (default RCOR2=1 ⇒ 0)
        rhcon[sp] = (ctl.regh_cor2_on && ctl.regh_cor2[sp] > 0f0) ? log(ctl.regh_cor2[sp]) : 0f0
    end
    return ip, rhcon
end

"""BC `small_tree_growth!` — V3 small-tree height + DBH increment (regent.f main body)."""
function small_tree_growth!(s::StandState, stash, ::BritishColumbia; fint::Float32 = 10.0f0)
    p, t = s.plot, s.trees
    t.n == 0 && return s
    n = t.n
    zone, series = bc_stand_zone(s)
    ip, rhcon = bc_regcons!(s)
    ba = p.basal_area; relden = p.relative_density; avh = p.avg_height
    aspect = p.aspect; slope = p.slope; dgsd = s.control.dg_sd
    regyr = 5.0f0
    # V2 (LV2ATV) small-tree: CON = RHCON(habitat) + HCOR(htg_cor_small); exp-form HTGRL; XMAXV2 bounds;
    # multiplicative ZZRAN (HSIGMA). No ST_COEF/ip needed. (regent.f LV2ATV branches.)
    v2 = bc_lv2atv(zone)
    nsp = nspecies(BritishColumbia())
    con_v2 = zeros(Float32, nsp)
    if v2
        regch_v2 = bc_v2_regch(aspect, slope)
        @inbounds for sp in 1:nsp
            con_v2[sp] = bc_v2_rhcon(sp, regch_v2) + s.calib.htg_cor_small[sp]
        end
    end
    # subcycle count + lengths (regent.f:186-203, mirror KT)
    ntyr = Int(round(fint)); iyr = Int(regyr)
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1); nper < 1 && (nper = 1)
    kper = zeros(Int, nper); itot = ntyr; nn = nper
    @inbounds for i in 1:nper
        if nn == 1; kper[i] = itot; break; end
        kper[i] = itot ÷ nn; itot -= kper[i]; nn -= 1
    end
    # per-subcycle stand density from LARGE trees (regent.f:236-256, via bc_tree_ccf)
    banext = fill(ba, nper); rdnext = fill(relden, nper)
    if nper > 1
        @inbounds for i in 1:n
            d1 = t.dbh[i]; d1 < 3.0f0 && continue
            sp = Int(t.species[i]); pr = t.tpa[i]
            bark = bc_bratio(sp)
            d2 = d1 + t.diam_growth[i] / bark
            b1 = 0.005454154f0*d1*d1; b2 = 0.005454154f0*d2*d2
            # CCFCAL (regent.f:1173-1174) computes the per-tree CCF WITH the tree's TPA (C=ccf·P), so CI=P·Δccf/10
            # and CI/P·PN recovers Δccf/10·PN. Using P=1 here dropped the ×TPA ⇒ RDNEXT barely moved (70.4 vs the
            # oracle's −163.79) ⇒ the small-tree HTGRL CCF term collapsed ⇒ ~4× height under-growth. Use pr.
            c1 = bc_tree_ccf(sp, d1, pr); c2 = bc_tree_ccf(sp, d2, pr)
            bi = (b2 - b1) / 10f0; ci = (c2 - c1) / 10f0
            k = 0
            for j in 2:nper
                k += kper[j-1]; pn = pr * 0.985f0^k
                rdnext[j] += k * ci / pr * pn
                banext[j] += k * bi * pn
            end
        end
    end
    delmax = (avh/36f0)*(0.01232f0*relden - 1.75f0); delmax > 0f0 && (delmax = 0f0)
    wk3 = Float32[t.height[i] for i in 1:n]
    # subcycle height accumulation (regent.f:1273-1458). Each subcycle ALSO adds the SMALL trees' own CCF/BA
    # contribution to the NEXT subcycle's density (regent.f:1450-1458) — omitting it collapsed RDNEXT (e.g.
    # 74 vs the oracle's −163.79 for garbage-height trees, whose power-form D1 from H1 dwarfs the actual DBH),
    # which killed the small-tree HTGRL CCF term ⇒ ~4× height under-growth. KY = cumulative years thru subcycle j.
    ky = 0
    @inbounds for j in 1:nper
        baj = banext[j]; rdj = rdnext[j]; kpj = Float32(kper[j]); ky += kper[j]
        decay = 0.985f0 ^ ky
        for i in 1:n
            sp = Int(t.species[i]); d = t.dbh[i]
            xmx_sp = v2 ? BC_RG_V2_XMAX[sp] : BC_RG_XMAX[sp]
            (d >= xmx_sp || t.tpa[i] <= 0f0 || (!v2 && ip[sp] < 1)) && continue
            h1 = wk3[i]; pct = t.crown_ratio[i]; pr = t.tpa[i]
            bal = (1f0 - pct/100f0) * baj
            incr = v2 ? exp(max(-40f0, con_v2[sp] + BC_RG_V2_RHLH[sp]*log(h1) +
                                BC_RG_V2_RHCCF[sp]*rdj + BC_RG_V2_RHBAL[sp]*bal)) :
                        bc_v3_sthg(sp, ip[sp], h1, bal, rdj, rhcon[sp], aspect, slope)
            xrhgro = active_multiplier(s.control, :regh, sp, current_cycle_year(s))
            h2 = h1 + incr * (kpj/regyr) * xrhgro
            wk3[i] = h2
            # small-tree density contribution to RDNEXT(j+1)/BANEXT(j+1) (regent.f:1437-1458). Skip last subcycle
            # and D≥3in (large trees handled by the pre-loop projection); needs H2>4.5 for the power-form DBH.
            if j < nper && d < 3f0 && h2 > 4.5f0
                relh = abs(avh - 4.5f0) < 0.01f0 ? 0f0 : clamp((h1 - 4.5f0)/(avh - 4.5f0), 0f0, 1f0)
                dadj = delmax*relh*relh - 2f0*delmax*relh + 0.65f0
                d1pf = h1 > 4.5f0 ? BC_RG_HHT1[sp]*(h1 - 4.5f0)^BC_RG_HHT2[sp] + dadj : BC_RG_DIAM[sp] + dadj
                d2pf = BC_RG_HHT1[sp]*(h2 - 4.5f0)^BC_RG_HHT2[sp] + dadj
                xrdgro = active_multiplier(s.control, :regd, sp, current_cycle_year(s))
                dgj = (d2pf - d1pf) * xrdgro; dgj < 0f0 && (dgj = 0f0)
                d2 = d + dgj
                c1 = bc_tree_ccf(sp, d1pf, pr); c2 = bc_tree_ccf(sp, d2, pr)
                b1 = 0.005454154f0 * d * d
                rdnext[j+1] += Float32(ky) * (c2 - c1) / 10f0 * decay
                banext[j+1] += (0.005454154f0*d2*d2 - b1) * pr * decay
            end
        end
    end
    # final: HTGR1 + ZZRAN + blend + DBH-dub (regent.f:473-560), species order for RNG determinism
    order = sortperm(view(t.species, 1:n); alg = Base.Sort.MergeSort)
    @inbounds for oi in 1:n
        i = order[oi]
        sp = Int(t.species[i]); d = t.dbh[i]
        xmn = v2 ? BC_RG_V2_XMIN[sp] : BC_RG_XMIN[sp]; xmx = v2 ? BC_RG_V2_XMAX[sp] : BC_RG_XMAX[sp]
        (d >= xmx || t.tpa[i] <= 0f0 || (!v2 && ip[sp] < 1)) && continue
        h = t.height[i]
        htgr1 = wk3[i] - h; htgr1 < 0f0 && (htgr1 = 0f0)
        zzran = 0f0
        if dgsd >= 1f0
            while true
                zzran = bachlo(s.rng, 0f0, 1f0)
                (zzran <= 1f0 && zzran >= -1.5f0) && break
            end
        end
        # V2: MULTIPLICATIVE error HTGR=HTGR1·exp(ZZRAN·HSIGMA) (regent.f:1544); V3: additive + ST_COEF.SD
        htgr = v2 ? htgr1 * exp(zzran * BC_RG_V2_HSIGMA) :
                    max(0f0, htgr1 + zzran * BC_STCOEF[ip[sp]].SD)
        xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)
        htg = htgr*(1f0 - xwt) + xwt*t.ht_growth[i]
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        t.ht_growth[i] = htg
        if d < 3.0f0                                              # DBH dub (regent.f:1355-1445)
            relh = (h - 4.5f0)/(avh - 4.5f0); relh = clamp(relh, 0f0, 1f0)
            dadj = delmax*relh*relh - 2f0*delmax*relh + 0.65f0
            d1 = bc_st_dbh(sp, h, dadj); dk = bc_st_dbh(sp, h + htg, dadj)
            xrdgro = active_multiplier(s.control, :regd, sp, current_cycle_year(s))
            # dg = outside-bark increment (dk−d1). The oracle's inside-bark DDS round-trip (regent.f:1618-1625)
            # cancels for SCALE=YR/FINT=1 when applied to the outside-bark DBH, so keep the direct increment.
            dg = (dk - d1) * xrdgro; dg < 0f0 && (dg = 0f0)
            t.diam_growth[i] = dg
        end
    end
    return s
end

# --- V2 (LV2ATV) small-tree coefficients (regent.f DATA) — verified via instrument (RHCON(14)=1.090813). ---
# RHCON(sp) = REGCH + RHSC(sp) + RHHAB(IRHHAB,sp); REGCH = RHGL(IGL=2) + (RSAB0+RSAB1·cosA+RSAB2·sinA)·slope.
# IRHHAB = regent.f's OWN MAPHAB(ITYPE=4,sp) [4th distinct subsystem MAPHAB]. See v2_smalltree_data.txt.
const BC_RG_V2_RHLH  = Float32[0.4214,0.2716,0.3907,0.3487,0.3417,0.2354,0.5843,0.2827,0.374,0.4485,0.2354,0.2354,0.2354,0.3907,0.2354]
const BC_RG_V2_RHCCF = Float32[-0.00591,-0.00654,-0.00591,-0.00391,-0.00391,-0.00391,-0.00654,-0.00391,-0.00391,-0.00654,-0.00391,-0.00391,-0.00391,-0.00591,-0.00391]
const BC_RG_V2_RHBAL = Float32[-0.37199,-0.41532,-0.40043,-0.25355,-0.34693,-0.12013,-0.24172,-0.253,-0.22957,-0.47299,-0.25349,-0.25349,-0.25349,-0.40043,-0.25349]
const BC_RG_V2_RHSC  = Float32[1.47,1.6204,1.4932,0.9981,1.0202,0.8953,1.2336,1.0964,1.0667,1.7311,0.8953,0.8953,0.8953,1.4932,0.8953]
const BC_RG_V2_IRHHAB = Int[3,3,4,3,1,1,5,1,4,3,3,3,3,4,3]     # regent MAPHAB(ITYPE=4, sp)
const BC_RG_V2_RHHAB = ([  # [sp][1..6]
    Float32[-0.2146,-0.0941,-0.3141,0,0,0], Float32[-0.2146,-0.0941,-0.3296,0,0,0],
    Float32[-0.2146,-0.0941,-0.5401,-0.3948,0,0], Float32[-0.2146,-0.0941,-0.2776,0,0,0],
    Float32[-0.2146,-0.0941,0,0,0,0], Float32[-0.2146,-0.0941,0,0,0,0],
    Float32[-0.2146,-0.0941,-0.2484,-0.5134,-0.3495,0], Float32[-0.2146,-0.0941,-0.3431,0,0,0],
    Float32[-0.2146,-0.0941,-0.4916,-0.3582,0,0], Float32[-0.2146,-0.0941,-0.4345,0,0,0],
    Float32[-0.2146,-0.0941,-0.3738,0,0,0], Float32[-0.2146,-0.0941,-0.3738,0,0,0],
    Float32[-0.2146,-0.0941,-0.3738,0,0,0], Float32[-0.2146,-0.0941,-0.5401,-0.3948,0,0],
    Float32[-0.2146,-0.0941,-0.3738,0,0,0],
]...,)
const BC_RG_V2_XMAX = Float32[10,10,10,10,10,10,5,10,10,10,10,10,10,10,10]
const BC_RG_V2_XMIN = Float32[2,2,2,2,2,2,1,2,2,2,2,2,2,2,2]
const BC_RG_V2_RHGL = Float32[-0.2785,-0.0480,0.0]   # RHGL(IGL); IGL=2 (grinit.f:204)
const BC_RG_V2_RSAB = Float32[-0.10987,0.22157,-0.12432]   # RSAB0/1/2 (aspect/slope)
const BC_RG_V2_HSIGMA = 0.59f0

"""V2 stand-level REGCH = RHGL(2) + (RSAB0 + RSAB1·cosA + RSAB2·sinA)·slope (regent.f:2020)."""
bc_v2_regch(aspect::Real, slope::Real) = BC_RG_V2_RHGL[2] +
    (BC_RG_V2_RSAB[1] + BC_RG_V2_RSAB[2]*cos(Float32(aspect)) + BC_RG_V2_RSAB[3]*sin(Float32(aspect))) * Float32(slope)

"""V2 per-species RHCON = REGCH + RHSC(sp) + RHHAB(IRHHAB,sp) (regent.f:2024, NI case, no RCOR2)."""
bc_v2_rhcon(sp::Integer, regch::Real) =
    Float32(regch) + BC_RG_V2_RHSC[sp] + BC_RG_V2_RHHAB[sp][BC_RG_V2_IRHHAB[sp]]
