# =============================================================================
# regent.jl (eastcascades) — EC small-tree growth (ec/regent.f + ec/smhtgf.f) + htdbh. Chunk 6.
#
# EC's small-tree model is BESPOKE (not the WC smhgdg reuse): ec/smhtgf.f is HEIGHT-growth-only with
# inline per-species Chapman-Richards/curve literals (WP/RC Chapman-Richards; MH/OS ×3.281 metric), and
# small-tree DIAMETER growth lives in ec/regent.f (per-species H→DBH inversion bounded by DGMAX·SCALE,
# SLO/SHI SI clamp, no DGMIN). Missing-height dubbing + small-tree DBH inversion use the forest-dependent
# Curtis-Arney htdbh (ec/htdbh.f: MTHOOD/OKANOG/GIFFPC/WENATC by IFOR/IGL).
# =============================================================================

# ec/htdbh.f — 4 physical forest tables (MTHOOD/OKANOG/GIFFPC/WENATC), stored in htdbh_coeffs_ec.csv by
# resolved CSV forest_idx {1=MTHOOD, 2/4/7=OKANOG, 3=WENATC(IGL1), 5/6=GIFFPC(IGL2,3)} × 32 species.
let
    P2 = zeros(Float32, 7, 32); P3 = zeros(Float32, 7, 32); P4 = zeros(Float32, 7, 32)
    for l in readlines(joinpath(EC_DATADIR, "htdbh_coeffs_ec.csv"))[2:end]
        f = split(strip(l), ','); (isempty(f) || isempty(f[1])) && continue
        fi = parse(Int, f[1]); sp = parse(Int, f[2])
        P2[fi, sp] = parse(Float32, f[3]); P3[fi, sp] = parse(Float32, f[4]); P4[fi, sp] = parse(Float32, f[5])
    end
    global const EC_HTDBH_P2 = P2; global const EC_HTDBH_P3 = P3; global const EC_HTDBH_P4 = P4
end

# ec/htdbh.f SELECT CASE(IFOR)/CASE(IGL) → CSV forest_idx. forkod passes the CORRECTED IFOR (1..4) + IGL.
@inline function _ec_htdbh_fidx(ifor::Int, igl::Int)::Int
    ifor == 1 && return 1                 # MTHOOD
    (ifor == 2 || ifor == 4) && return 2  # OKANOG
    (ifor == 3 && (igl == 2 || igl == 3)) && return 5  # GIFFPC
    return 3                              # WENATC (IFOR 3, IGL 1) / default
end

@inline function ec_htdbh_height(fidx::Int, sp::Int, d::Float32)::Float32
    p2 = EC_HTDBH_P2[fidx, sp]; p3 = EC_HTDBH_P3[fidx, sp]; p4 = EC_HTDBH_P4[fidx, sp]
    if d >= 3.0f0
        return 4.5f0 + p2 * exp(-1f0 * p3 * d^p4)
    else
        return ((4.5f0 + p2 * exp(-1f0 * p3 * 3.0f0^p4) - 4.51f0) * (d - 0.3f0) / 2.7f0) + 4.51f0
    end
end
@inline function ec_htdbh_dbh(fidx::Int, sp::Int, h::Float32)::Float32
    p2 = EC_HTDBH_P2[fidx, sp]; p3 = EC_HTDBH_P3[fidx, sp]; p4 = EC_HTDBH_P4[fidx, sp]
    hat3 = 4.5f0 + p2 * exp(-1f0 * p3 * 3.0f0^p4)
    if h >= hat3
        return exp(log((log(h - 4.5f0) - log(p2)) / (-1f0 * p3)) * (1f0 / p4))
    else
        return (((h - 4.51f0) * 2.7f0) / (hat3 - 4.51f0)) + 0.3f0
    end
end

@inline ec_htdbh_ifor(p)::Int = _ec_htdbh_fidx(Int(p.forest_idx), Int(p.geo_location))

regenerate!(s::StandState, ::EastCascades; kwargs...) = s
# mortality!(::EastCascades) = the shared Stage/Reineke driver mortality!(::AbstractVariant)
# (southern/mortality.jl); EC's only variant hook is _varmrt_efftr!(::EastCascades) in mortality.jl.


# ── ec/regent.f + ec/smhtgf.f small-tree growth (chunk 6).
const EC_RG_DGMAX = Float32[2.8,2.8,2.4,3.6,2.5,2.5,3.5,3.6,3.6,2.8,5.0,2.8,5.0,5.0,5.0,2.5,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,2.8,5.0]
const EC_RG_XMAX  = Float32[4.0,4.0,4.0,4.0,10.0,4.0,5.0,4.0,6.0,6.0,4.0,6.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,4.0,6.0,4.0]
const EC_RG_XMIN  = Float32[2.0,2.0,2.0,2.0,2.0,2.0,1.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0]
const EC_RG_SLO   = Float32[20,50,50,50,15,50,30,40,50,70,0,15,0,0,0,50,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,0]
const EC_RG_SHI   = Float32[80,110,110,110,30,110,70,120,150,140,999,30,999,999,999,110,999,999,999,999,999,999,999,999,999,999,999,999,999,999,30,999]
const EC_RG_DIAM  = Float32[0.4,0.3,0.3,0.3,0.2,0.3,0.4,0.3,0.3,0.5,0.2,0.2,0.2,0.4,0.3,0.3,0.3,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2]
const EC_RG_AB    = Float64[1.11436,-0.011493,0.43012e-4,-0.72221e-7,0.5607e-10,-0.1641e-13]
const EC_RG_REGYR = 10.0f0

# ec/smhtgf.f — small-tree potential HEIGHT increment over DTIME years (MODE 1 = established, effective-age).
@inline function ec_smhtgf(ispc::Int, h::Float32, dtime::Float32, si::Float32)::Float32
    if ispc == 1                                       # WP — Chapman-Richards
        c1=0.375045f0; c2=0.92503f0; c3=-0.020796f0; c4=2.48811f0
        effage = log((1f0 - (c1/si*h)^(1f0/c4))/c2)/c3; agepdt = effage + dtime
        return (si/c1)*(1f0-c2*exp(c3*agepdt))^c4 - (si/c1)*(1f0-c2*exp(c3*effage))^c4
    elseif ispc == 5                                   # RC — Chapman-Richards
        c1=0.752842f0; c2=1.0f0; c3=-0.0174f0; c4=1.4711f0
        e = (1f0 - (c1/si*h)^(1f0/c4))/c2
        effage = e > 0f0 ? log(e)/c3 : 100f0; agepdt = effage + dtime
        return (si/c1)*(1f0-c2*exp(c3*agepdt))^c4 - (si/c1)*(1f0-c2*exp(c3*effage))^c4
    elseif ispc == 2 || ispc == 17
        return ((-3.97245f0 + 0.50995f0*si)/(28.11668f0-0.05661f0*si))*dtime
    elseif ispc == 3
        return ((2.0f0 + 0.420f0*si)/(28.5f0 - 0.05f0*si))*dtime
    elseif ispc == 4
        return ((-0.6667f0 + 0.4333f0*si)/(28.5f0 - 0.05f0*si))*dtime
    elseif ispc == 6 || ispc == 16
        return ((-1.0470f0 + 0.4220f0*si)/(28.7739f0 - 0.0597f0*si))*dtime
    elseif ispc == 7
        return (0.3277f0 + 0.01296f0*si)*dtime
    elseif ispc == 8
        return ((-8.0f0 + 0.35f0*si)/(53.72545f0 - 0.274509f0*si))*dtime
    elseif ispc == 9
        return ((6.0f0 + 0.14f0*si)/(33.882f0 - 0.06588f0*si))*dtime
    elseif ispc == 10
        return ((-1.0f0 + 0.32857f0*si)/(28.0f0 - 0.042857f0*si))*dtime
    elseif ispc == 11
        return ((-5.74874f0 + 0.54576f0*si)/(26.15767f0 - 0.03596f0*si))*dtime
    elseif ispc == 12 || ispc == 31
        return ((0.965758f0 + 0.082969f0*si)/(55.249612f0 - 1.288852f0*si))*dtime*3.280833f0
    elseif ispc == 15
        return ((11.26677f0 + 0.12027f0*si)/(27.93806f0 - 0.02873f0*si))*dtime
    elseif ispc == 22
        return (-0.007025f0 + 0.056794f0*si)*dtime
    elseif ispc == 28
        return (((-37.60812f0*log(1f0-(si/114.24569f0)^0.44444f0))*0.01f0)-0.1f0)*dtime
    else                                               # 13,14,18:21,23:27,29,30,32
        return ((1.47043f0 + 0.23317f0*si)/(31.56252f0 - 0.05586f0*si))*dtime
    end
end

# ec/essubh.f — planted/subsequent-tree base height: SMHTGF called with MODE=0, DTIME=AGE. For the two
# Chapman-Richards species (WP=1, RC=5) MODE=0 forces EFFAGE=0 (HHT = curve(AGE) − curve(0)); every other
# species is linear (coef·AGE, H unused) so ec_smhtgf already gives the MODE=0 value. SI = the species' SITEAR.
@inline function ec_essubh_hht(sp::Int, si::Float32, age::Float32)::Float32
    age <= 0f0 && return 0f0                            # SMHTGF: DTIME≤0 → no growth (ec/smhtgf.f:97)
    if sp == 1
        c1=0.375045f0; c2=0.92503f0; c3=-0.020796f0; c4=2.48811f0
        return (si/c1)*(1f0-c2*exp(c3*age))^c4 - (si/c1)*(1f0-c2)^c4     # EFFAGE=0
    elseif sp == 5
        c1=0.752842f0; c2=1.0f0; c3=-0.0174f0; c4=1.4711f0
        return (si/c1)*(1f0-c2*exp(c3*age))^c4 - (si/c1)*(1f0-c2)^c4     # EFFAGE=0
    else
        return ec_smhtgf(sp, 0f0, age, si)             # linear species: H unused ⇒ = coef·AGE
    end
end

function small_tree_growth!(s::StandState, stash, ::EastCascades; fint::Float32 = 10.0f0)
    p, t, c = s.plot, s.trees, s.calib
    n = t.n; n == 0 && return s
    sd = s.coef.species
    avh = p.avg_height; dgsd = s.control.dg_sd
    relden = p.relative_density
    ifor = ec_htdbh_ifor(p)
    isisp = Int(p.site_species); (isisp < 1 || isisp > 32) && (isisp = 10)
    scale = fint / EC_RG_REGYR; scale2 = EC_RG_REGYR / fint
    # density modifier PCTRED (ec/regent.f): X = AVH·(RELDEN/100), capped 300; 5th-order polynomial.
    xden = avh * (relden / 100f0); xden > 300f0 && (xden = 300f0)
    pctred = Float32(EC_RG_AB[1] + xden*(EC_RG_AB[2] + xden*(EC_RG_AB[3] + xden*(EC_RG_AB[4] + xden*(EC_RG_AB[5] + xden*EC_RG_AB[6])))))
    pctred > 1f0 && (pctred = 1f0); pctred < 0.01f0 && (pctred = 0.01f0)
    @inbounds for i in 1:n
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= EC_RG_XMAX[sp] || t.tpa[i] <= 0f0) && continue
        h = t.height[i]
        # site index for smhtgf = the tree's own SITEAR, clamped to the site species' [SLO+0.5, SHI]
        si = p.sp_site_index[sp]
        si > EC_RG_SHI[isisp] && (si = EC_RG_SHI[isisp])
        si <= EC_RG_SLO[isisp] && (si = EC_RG_SLO[isisp] + 0.5f0)
        con = exp(c.htg_cor_small[sp])                  # RHCON=1 default ⇒ CON = exp(HCOR)
        x = Float32(t.crown_pct[i]) / 100f0
        vigor = 150f0 * x^3 * exp(-6f0*x) + 0.3f0; vigor > 1f0 && (vigor = 1f0)
        pothtg = ec_smhtgf(sp, h, 10f0, si)
        htgr = pothtg * pctred * vigor * con
        zzran = 0f0
        if dgsd >= 1f0
            while true
                zzran = bachlo(s.rng, 0f0, 1f0)
                (zzran <= 0.5f0 && zzran >= -2f0) && break
            end
        end
        wcform = sp in EC_HT_WCFORM
        if wcform
            htgr = (htgr + zzran*0.1f0) * scale; htgr < 0.1f0 && (htgr = 0.1f0)
        else
            htgr = (htgr + zzran*0.1f0) * scale; htgr < 0f0 && (htgr = 0f0)
        end
        xmn = EC_RG_XMIN[sp]; xmx = EC_RG_XMAX[sp]
        xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)
        htg = htgr*(1f0 - xwt) + xwt*t.ht_growth[i]
        htg < 0.1f0 && (htg = 0.1f0)
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = cap - h; htg < 0.1f0 && (htg = 0.1f0))
        t.ht_growth[i] = htg
        if d >= 3f0                                     # DBH from large-tree dgf; only HTG overridden
            _ec_rg_stash!(stash, t, i); continue
        end
        bark = wc_bratio(sd, sp, d)
        hk = h + htg
        if hk <= 4.5f0
            t.diam_growth[i] = 0f0; t.dbh[i] = d + 0.001f0*hk
            _ec_rg_stash!(stash, t, i); continue
        end
        dk  = ec_htdbh_dbh(ifor, sp, hk)
        dkk = h <= 4.5f0 ? d : ec_htdbh_dbh(ifor, sp, h)
        local dg::Float32
        if dk < 0f0 || dkk < 0f0
            dg = htg * 0.2f0 * bark                      # xrdgro = 1
            dk = d + dg
        else
            dg = (dk - dkk) * bark
        end
        dg < 0f0 && (dg = 0f0)
        dgmx = EC_RG_DGMAX[sp] * scale
        dg > dgmx && (dg = dgmx)
        dds = dg*(2f0*bark*d + dg)*scale2
        dg = sqrt((d*bark)^2 + dds) - bark*d
        (t.dbh[i] + dg) < EC_RG_DIAM[sp] && (dg = EC_RG_DIAM[sp] - t.dbh[i])
        dg = wc_dgbnd(sp, t.dbh[i], dg, s.control.sp_size_cap[sp, 1], s.control.sp_size_cap[sp, 3])
        t.diam_growth[i] = dg
        _ec_rg_stash!(stash, t, i)
    end
    return s
end

@inline function _ec_rg_stash!(stash, t, i::Int)
    if stash !== nothing && !isempty(stash.dgU) && i <= length(stash.dgU)
        stash.dgU[i] = t.diam_growth[i]; stash.dgL[i] = t.diam_growth[i]
        stash.htgU[i] = t.ht_growth[i]; stash.htgL[i] = t.ht_growth[i]
        !isempty(stash.is_small) && (stash.is_small[i] = true)
    end
end
