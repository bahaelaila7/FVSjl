# =============================================================================
# regent.jl (pacificnorthwest) — PN small-tree growth. Chunk 6.
#
# PN reuses WC's regent.f + vwc/smhgdg.f (byte-identical) ⇒ the WC_RG_* size/diam DATA + the SMHGDG
# DMAX/BETA/ALPHA (WC_SMH_*) + wc_dgbnd are shared. TWO PN differences: (1) `pn_smhgdg` omits the DF
# Curtis→King SI line (VARACD≠'WC' ⇒ raw SITEAR); (2) PN's forest-dependent htdbh (htdbh_coeffs_pn.csv,
# 6 tables) for missing-height dubbing + redwood DBH. Otherwise the driver is identical to WC's.
# =============================================================================

# pn/htdbh.f — 6-forest × 39-species Curtis-Arney P2/P3/P4 (data/pacificnorthwest/htdbh_coeffs_pn.csv).
let
    P2 = zeros(Float32, 6, 39); P3 = zeros(Float32, 6, 39); P4 = zeros(Float32, 6, 39)
    for l in readlines(joinpath(PN_DATADIR, "htdbh_coeffs_pn.csv"))[2:end]
        f = split(strip(l), ','); (isempty(f) || isempty(f[1])) && continue
        fi = parse(Int, f[1]); sp = parse(Int, f[2])
        P2[fi, sp] = parse(Float32, f[3]); P3[fi, sp] = parse(Float32, f[4]); P4[fi, sp] = parse(Float32, f[5])
    end
    global const PN_HTDBH_P2 = P2; global const PN_HTDBH_P3 = P3; global const PN_HTDBH_P4 = P4
end

# PN forkod IFOR (1..6) is used directly as the htdbh forest-table index (no BLM remap).
@inline _pn_htdbh_ifor(ifor::Int)::Int = clamp(ifor, 1, 6)

@inline function pn_htdbh_height(ifor::Int, sp::Int, d::Float32)::Float32
    p2 = PN_HTDBH_P2[ifor, sp]; p3 = PN_HTDBH_P3[ifor, sp]; p4 = PN_HTDBH_P4[ifor, sp]
    if d >= 3.0f0
        return 4.5f0 + p2 * exp(-1f0 * p3 * d^p4)
    else
        return ((4.5f0 + p2 * exp(-1f0 * p3 * 3.0f0^p4) - 4.51f0) * (d - 0.3f0) / 2.7f0) + 4.51f0
    end
end
@inline function pn_htdbh_dbh(ifor::Int, sp::Int, h::Float32)::Float32
    p2 = PN_HTDBH_P2[ifor, sp]; p3 = PN_HTDBH_P3[ifor, sp]; p4 = PN_HTDBH_P4[ifor, sp]
    hat3 = 4.5f0 + p2 * exp(-1f0 * p3 * 3.0f0^p4)
    if h >= hat3
        return exp(log((log(h - 4.5f0) - log(p2)) / (-1f0 * p3)) * (1f0 / p4))
    else
        return (((h - 4.51f0) * 2.7f0) / (4.5f0 + p2 * exp(-1f0 * p3 * 3.0f0^p4) - 4.51f0)) + 0.3f0
    end
end

# vwc/smhgdg.f — identical to wc_smhgdg EXCEPT NO DF Curtis→King line (raw SITEAR for PN).
@inline function pn_smhgdg(sp::Int, h::Float32, d::Float32, cr::Float32, ptbal::Float32,
                           ptba::Float32, si::Float32, avht::Float32)
    if sp == 17
        htmax = 2.242202f0 * si
        (htmax - h <= 1.0f0) && return (0.0f0, 0.0f0)
        age1 = (1.0f0 / -0.010742f0) * log(1.0f0 - (h / 2.242202f0 / si)^(1.0f0 / 0.919076f0))
        age2 = age1 + 5.0f0
        h1 = 2.242202f0 * si * (1.0f0 - exp(-0.010742f0 * age1))^0.919076f0
        h2 = 2.242202f0 * si * (1.0f0 - exp(-0.010742f0 * age2))^0.919076f0
        return (h2 - h1, 0.0f0)
    end
    relht = 0.0f0
    avht > 0.0f0 && (relht = h / avht)
    relht > 1.5f0 && (relht = 1.5f0)
    # (PN: NO sp==16 Curtis→King — raw SITEAR)
    ptbal2 = log(ptbal + 2.71f0); ptba2 = log(ptba + 2.71f0); relht2 = sqrt(relht)
    boost = ptba < 100.0f0 ? 1.0f0 / (1.0f0 + exp(-3.1f0 + 0.18f0 * ptba)) : 0.0f0
    hbh = h >= 4.5f0 ? 4.5f0 : h
    a = WC_SMH_ALPHA[sp]; beta = WC_SMH_BETA[sp]
    dgs = WC_SMH_DMAX[sp] / (1.0f0 + exp(a[1] + a[2]*ptba + a[3]*ptba2 + a[4]*ptbal + a[5]*ptbal2 +
              a[6]*boost + a[7]*cr + a[8]*relht + a[9]*relht2 + a[10]*si))
    hg5 = dgs / beta
    local dg5::Float32
    if h < 4.5f0
        dg5 = (h + hg5) >= 4.5f0 ? beta * (h + hg5 - 4.5f0) : 0.0f0
    else
        dg5 = dgs
    end
    return (hg5, dg5)
end

function small_tree_growth!(s::StandState, stash, ::PacificNorthwest; fint::Float32 = 10.0f0)
    p, t, c = s.plot, s.trees, s.calib
    n = t.n; n == 0 && return s
    sd = s.coef.species; dens = s.density
    avh = p.avg_height; dgsd = s.control.dg_sd
    ifor = _pn_htdbh_ifor(Int(p.forest_idx))
    scale = fint / _WC_RG_REGYR; scale2 = _WC_RG_REGYR / fint
    avht = avh
    @inbounds for i in 1:n
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= WC_RG_XMAX[sp] || t.tpa[i] <= 0.0f0) && continue
        h = t.height[i]
        cr = Float32(t.crown_pct[i]) * 0.01f0
        ip = Int(t.plot_id[i])
        ptbal = dens.point_bal[i]
        ptba = (1 <= ip <= length(dens.point_ba)) ? dens.point_ba[ip] : 0.0f0
        si = p.sp_site_index[sp]
        con = exp(c.htg_cor_small[sp])
        wk4 = t.htimlt[i]
        hg1, dg1 = pn_smhgdg(sp, h, d, cr, ptbal, ptba, si, avht)
        hk = h + hg1; dk = d + dg1
        hg2, dg2 = pn_smhgdg(sp, hk, dk, cr, ptbal, ptba, si, avht)
        htgr = hg1 + hg2; dgr = dg1 + dg2
        zzran = 0.0f0
        if dgsd >= 1.0f0
            while true
                zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                (zzran <= 0.5f0 && zzran >= -2.0f0) && break
            end
        end
        htgr = (htgr + zzran * 0.1f0) * scale * con * wk4
        htgr < 0.1f0 && (htgr = 0.1f0)
        xmn = WC_RG_XMIN[sp]; xmx = WC_RG_XMAX[sp]
        xwt = d <= xmn ? 0.0f0 : (d - xmn) / (xmx - xmn)
        if sp == 17
            lthg = t.ht_growth[i]
            htgr = (htgr + lthg) / 2.0f0
            t.ht_growth[i] = htgr * (1.0f0 - xwt) + xwt * lthg
        end
        htg = htgr * (1.0f0 - xwt) + xwt * t.ht_growth[i]
        htg < 0.1f0 && (htg = 0.1f0)
        cap = s.control.sp_size_cap[sp, 4]
        if (h + htg) > cap
            htg = cap - h; htg < 0.1f0 && (htg = 0.1f0)
        end
        t.ht_growth[i] = htg
        if d >= WC_RG_DGMIN[sp]
            _wc_rg_stash!(stash, t, i); continue
        end
        hkk = h + htg
        if hkk < 4.5f0
            t.diam_growth[i] = 0.0f0
            t.dbh[i] = d + 0.001f0 * hkk
            _wc_rg_stash!(stash, t, i); continue
        end
        bark = wc_bratio(sd, sp, d)
        local dgk::Float32
        if sp == 17
            dk2 = pn_htdbh_dbh(ifor, sp, hkk)
            dkk = h <= 4.5f0 ? d : pn_htdbh_dbh(ifor, sp, h)
            xdwt = d <= xmn ? 0.0f0 : (d - xmn) / (7.0f0 - xmn)
            dgsm = (dk2 - dkk) * bark; dgsm < 0.0f0 && (dgsm = 0.0f0)
            dds = dgsm * (2.0f0 * bark * d + dgsm) * scale2
            dgsm = sqrt((d * bark)^2 + dds) - bark * d
            dgk = dgsm * (1.0f0 - xdwt) + t.diam_growth[i] * xdwt
        else
            dgk = dgr * scale * wk4
            if d < 0.0f0 || dgk < 0.0f0
                dgk = htg * 0.2f0 * bark
            else
                dgk = dgk * bark
            end
            dgk < 0.0f0 && (dgk = 0.1f0)
            dgmx = WC_RG_DGMAX[sp] * scale
            dgk > dgmx && (dgk = dgmx)
            dds = dgk * (2.0f0 * bark * d + dgk) * scale2
            dgk = sqrt((d * bark)^2 + dds) - bark * d
        end
        (t.dbh[i] + dgk) < WC_RG_DIAM[sp] && (dgk = WC_RG_DIAM[sp] - t.dbh[i])
        dgk = wc_dgbnd(sp, t.dbh[i], dgk, s.control.sp_size_cap[sp, 1], s.control.sp_size_cap[sp, 3])
        t.diam_growth[i] = dgk
        _wc_rg_stash!(stash, t, i)
    end
    return s
end

regenerate!(s::StandState, ::PacificNorthwest; kwargs...) = s
