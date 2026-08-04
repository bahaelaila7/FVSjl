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
            c1 = bc_tree_ccf(sp, d1, 1f0); c2 = bc_tree_ccf(sp, d2, 1f0)
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
    # subcycle height accumulation (regent.f:321-460)
    @inbounds for j in 1:nper
        baj = banext[j]; rdj = rdnext[j]; kpj = Float32(kper[j])
        for i in 1:n
            sp = Int(t.species[i]); d = t.dbh[i]
            (d >= BC_RG_XMAX[sp] || t.tpa[i] <= 0f0 || ip[sp] < 1) && continue
            h1 = wk3[i]; pct = t.crown_ratio[i]
            bal = (1f0 - pct/100f0) * baj
            htgrl = bc_v3_sthg(sp, ip[sp], h1, bal, rdj, rhcon[sp], aspect, slope)
            xrhgro = active_multiplier(s.control, :regh, sp, current_cycle_year(s))
            wk3[i] = h1 + htgrl * (kpj/regyr) * xrhgro
        end
    end
    # final: HTGR1 + ZZRAN + blend + DBH-dub (regent.f:473-560), species order for RNG determinism
    order = sortperm(view(t.species, 1:n); alg = Base.Sort.MergeSort)
    @inbounds for oi in 1:n
        i = order[oi]
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= BC_RG_XMAX[sp] || t.tpa[i] <= 0f0 || ip[sp] < 1) && continue
        h = t.height[i]; xmn = BC_RG_XMIN[sp]; xmx = BC_RG_XMAX[sp]
        htgr1 = wk3[i] - h; htgr1 < 0f0 && (htgr1 = 0f0)
        zzran = 0f0
        if dgsd >= 1f0
            while true
                zzran = bachlo(s.rng, 0f0, 1f0)
                (zzran <= 1f0 && zzran >= -1.5f0) && break
            end
        end
        htgr = max(0f0, htgr1 + zzran * BC_STCOEF[ip[sp]].SD)     # V3 ZZRAN (regent.f:1547)
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
            dg = (dk - d1) * xrdgro; dg < 0f0 && (dg = 0f0)
            t.diam_growth[i] = dg
        end
    end
    return s
end
