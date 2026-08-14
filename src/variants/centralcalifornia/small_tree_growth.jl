# =============================================================================
# small_tree_growth.jl (centralcalifornia) — CA small-tree growth (ca/regent.f + smhtgf.f). Chunk 6.
#
# Height-first small-tree model: SMHTGF gives a 5-yr height increment (5 MAPSP groups — pines, firs+DF,
# black oak, tanoak, coast redwood); ×SCALE(=FINT/REGYR) → FINT-yr HTGR, + ZZRAN·0.1 stochastic (DGSD≥1),
# blended with the large-tree HTG by XWT=(D−XMN)/(XMX−XMN). DBH is derived from the HEIGHT crossing 4.5':
# DK=htdbh⁻¹(H+HTG), DKK=htdbh⁻¹(H); DG = ((DK−DKK)·bark → DDS) blended with the large-tree DG by
# XDWT=(D−1.5)/1.5. Only D<DGMIN gets small-tree DG; D≥XMAX skips regent (keeps large-tree HTG/DG).
# =============================================================================

const CA_RG_XMAX  = Float32[fill(4.0f0,22); 10.0f0; fill(4.0f0,26); 10.0f0]   # ca/regent.f DATA XMAX (22*4,10,26*4,10)
const CA_RG_XMIN  = fill(2.0f0, 50)                                            # DATA XMIN/MAXSP*2.0/
const CA_RG_DGMIN = Float32[fill(3.0f0,22); 7.0f0; fill(3.0f0,26); 7.0f0]     # DATA DGMIN (22*3,7,26*3,7)
const CA_RG_REGYR = 5.0f0
const CA_RG_DIAM = Float32[
    0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.5, 0.5,0.4,0.5,0.5,0.5,0.5,0.3,0.5,0.5,0.5,
    0.3,0.3,0.3,0.3,0.3,0.2,0.2,0.2,0.2,0.2, 0.2,0.2,0.2,0.2,0.3,0.1,0.1,0.2,0.1,0.3,
    0.4,0.2,0.2,0.1,0.1,0.1,0.2,0.2,0.2,0.3]
# ca/smhtgf.f DATA SPADJF + MAPSP (50 sp → 5 height groups: 1=pines 2=firs+DF 3=blackoak 4=tanoak 5=redwood).
const CA_SMH_SPADJF = Float32[
    1.0,1.0,0.9,1.1,1.1,1.1,1.1,0.8,0.9,0.9, 1.0,1.0,1.0,1.0,1.0,1.1,1.1,1.0,1.1,0.9,
    1.0,0.9,1.0,0.8,1.0,1.1,0.9,1.1,1.1,1.0, 1.1,1.0,1.1,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
    1.1,1.0,1.1,1.2,1.2,1.1,0.8,1.0,1.0,1.0]
const CA_SMH_MAPSP = Int[
    2,2,2,2,2,2,2,2,2,1, 1,1,1,1,1,1,1,1,1,1, 1,2,5,2,1,3,3,3,3,3,
    3,3,3,4,3,3,4,3,4,3, 3,4,3,3,3,3,2,4,3,5]

# ca/smhtgf.f — 5-yr small-tree height increment. `cr` = ICR/10 (0–9.5, as regent passes it).
function ca_smhtgf(sp::Int, d::Float32, h::Float32, cr::Float32, ba::Float32, bal::Float32,
                   si::Float32, relht::Float32)::Float32
    factor = CA_SMH_SPADJF[sp] * (0.80f0 + 0.004f0*(si - 50f0))
    tembal = bal < 5f0 ? 5f0 : bal
    grp = CA_SMH_MAPSP[sp]
    local htgr::Float32
    if grp == 1                                     # pines
        htgr = exp(0.7452f0 - 0.003271f0*bal - 0.1632f0*cr + 0.0217f0*cr*cr + 0.00536f0*si) * factor * 1.75f0
    elseif grp == 2                                 # firs incl. Douglas-fir
        domhtgr = 5f0*(2.2227f0 + 0.4314f0*si)/(29.0f0 - 0.05f0*si)
        crr = cr/10f0
        crmod = 1f0 - exp(-4.26558f0*crr)
        rhmod = exp(2.54119f0*(relht^0.250537f0 - 1f0))
        htgr = domhtgr * 1.016605f0*crmod*rhmod
    elseif grp == 3                                 # black oak
        htgr = exp(3.817f0 - 0.7829f0*log(tembal)) * factor
    elseif grp == 4                                 # tanoak
        htgr = exp(3.385f0 - 0.5898f0*log(tembal)) * factor
    else                                            # coast redwood
        htmax = 2.242202f0*si
        if htmax - h <= 1f0
            htgr = 0f0
        else
            age1 = 1f0/-0.010742f0 * log(1f0 - (h/2.242202f0/si)^(1f0/0.919076f0))
            age2 = age1 + 5f0
            h1 = 2.242202f0*si*(1f0-exp(-0.010742f0*age1))^0.919076f0
            h2 = 2.242202f0*si*(1f0-exp(-0.010742f0*age2))^0.919076f0
            htgr = h2 - h1
        end
    end
    htgr < 0.1f0 && (htgr = 0.1f0)
    return htgr
end

# ca/dgbnd.f — size-cap DG bound only (no exp-DGMAX, unlike WC).
@inline function ca_dgbnd(sp::Int, dbh::Float32, ddg::Float32, sizcap1::Float32, sizcap3::Float32)::Float32
    if (dbh + ddg) > sizcap1 && sizcap3 < 1.5f0
        ddg = sizcap1 - dbh
        ddg < 0.01f0 && (ddg = 0.01f0)
    end
    return ddg
end

@inline function _ca_rg_stash!(stash, t, i::Int)
    if stash !== nothing && !isempty(stash.dgU) && i <= length(stash.dgU)
        stash.dgU[i] = t.diam_growth[i]; stash.dgL[i] = t.diam_growth[i]
        stash.htgU[i] = t.ht_growth[i]; stash.htgL[i] = t.ht_growth[i]
        !isempty(stash.is_small) && (stash.is_small[i] = true)
    end
end

function small_tree_growth!(s::StandState, stash, ::CentralCalifornia; fint::Float32 = 10.0f0)
    p, t, c = s.plot, s.trees, s.calib
    n = t.n; n == 0 && return s
    sd = s.coef.species
    avh = p.avg_height; ba = p.basal_area; dgsd = s.control.dg_sd
    scale = fint / CA_RG_REGYR                       # SCALE = FNT/REGYR (non-estab FNT=FINT)
    scale2 = s.control.year / fint                   # SCALE2 = YR/FNT
    @inbounds for i in 1:n
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= CA_RG_XMAX[sp] || t.tpa[i] <= 0f0) && continue
        h = t.height[i]
        cr = Float32(t.crown_pct[i]) / 10f0          # ca/regent.f CR = ICR/10
        bal = ((100f0 - t.crown_ratio[i]) / 100f0) * ba   # PCT (stand BA percentile) lives in crown_ratio
        relht = avh <= 0f0 ? 1f0 : h/avh
        relht > 1.05f0 && (relht = 1.05f0)
        si = p.sp_site_index[sp]
        bark = wc_bratio(sd[:bark1][sp], sd[:bark2][sp], Int(sd[:bark_imap][sp]), d)
        con = exp(c.htg_cor_small[sp])               # RHCON(=1)·exp(HCOR)
        htgrr = ca_smhtgf(sp, d, h, cr, ba, bal, si, relht)
        htgr = htgrr * con
        zzran = 0f0
        if dgsd >= 1f0
            while true
                zzran = bachlo(s.rng, 0f0, 1f0)
                (zzran <= 0.5f0 && zzran >= -2.0f0) && break
            end
        end
        htgr = (htgr + zzran*0.1f0) * scale          # XRHGRO=1
        xmn = CA_RG_XMIN[sp]; xmx = CA_RG_XMAX[sp]
        xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)
        if sp == 23 || sp == 50                      # GS/RW blend with the large-tree HTG
            lthg = t.ht_growth[i]
            htgr = (htgr + lthg)/2f0
            t.ht_growth[i] = htgr*(1f0-xwt) + xwt*lthg
        else
            t.ht_growth[i] = htgr*(1f0-xwt) + xwt*t.ht_growth[i]
        end
        t.ht_growth[i] < 0.1f0 && (t.ht_growth[i] = 0.1f0)
        cap = s.control.sp_size_cap[sp, 4]
        (h + t.ht_growth[i] > cap) && (t.ht_growth[i] = max(cap - h, 0.1f0))
        # --- small-tree DBH (only D<DGMIN; else keep large-tree DG) ---
        if d >= CA_RG_DGMIN[sp]
            _ca_rg_stash!(stash, t, i); continue
        end
        htg = t.ht_growth[i]; hk = h + htg
        if hk <= 4.5f0
            t.diam_growth[i] = 0f0
            t.dbh[i] = d + 0.001f0*hk
            _ca_rg_stash!(stash, t, i); continue
        end
        dk = ca_htdbh_dbh(sp, hk)
        dkk = h <= 4.5f0 ? d : ca_htdbh_dbh(sp, h)
        if sp == 23 || sp == 50
            xdwt = d <= xmn ? 0f0 : (d - xmn)/(CA_RG_DGMIN[sp] - xmn)
        else
            xdwt = d <= 1.5f0 ? 0f0 : d >= 3f0 ? 1f0 : (d - 1.5f0)/1.5f0
        end
        dgsm = (dk - dkk)*bark; dgsm < 0f0 && (dgsm = 0f0)   # XRDGRO=1
        dds = dgsm*(2f0*bark*d + dgsm)*scale2
        dgsm = sqrt((d*bark)^2 + dds) - bark*d; dgsm < 0f0 && (dgsm = 0f0)
        dglt = t.diam_growth[i]
        dgk = dgsm*(1f0-xdwt) + dglt*xdwt
        (t.dbh[i] + dgk) < CA_RG_DIAM[sp] && (dgk = CA_RG_DIAM[sp] - t.dbh[i])
        dgk = ca_dgbnd(sp, t.dbh[i], dgk, s.control.sp_size_cap[sp,1], s.control.sp_size_cap[sp,3])
        t.diam_growth[i] = dgk
        _ca_rg_stash!(stash, t, i)
    end
    return s
end
