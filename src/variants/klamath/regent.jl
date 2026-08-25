# =============================================================================
# regent.jl (klamath) — NC small-tree growth (nc/regent.f + htgr5.f + htdbh.f). Chunk 6.
#   nc_htgr5(sp,ssite,baa,relht,cr,h) — small-tree height increment (3 methods, htgr5.f)
#   nc_htdbh(sp,d,h,mode)             — Curtis-Arney HT-DBH (SISKIY, htdbh.f; mode1 H→D)
#   small_tree_growth!(::Klamath)     — HTGR5→XWT-blend w/ large-tree HTG + small-tree DG (D<DGMIN), blended
# For NC LHTDRG=.FALSE. ⇒ regent recomputes the small-tree DK/DKK via HTDBH (SISKIY), NOT the HT1/HT2 Wykoff.
# =============================================================================

const NC_REGYR = 5.0f0
const NC_ST_XMAX = Float32[5,5,5,5,5,5,5,5,5,5,5,10]
const NC_ST_XMIN = Float32[2,2,2,2,2,2,2,2,2,2,2,2]
const NC_ST_DGMIN = Float32[3,3,3,3,3,3,3,3,3,3,3,7]
const NC_ST_DIAM = Float32[0.3,0.4,0.3,0.3,0.2,0.2,0.2,0.2,0.3,0.5,0.2,0.3]
# htgr5.f
const NC_HG_IMETH = Int[1,1,1,1,2,1,2,2,1,1,2,3]
const NC_HG_HCON = Float32[-2.193,-2.193,-2.193,-2.193,3.560,-2.193,3.817,3.385,-2.193,-2.193,3.385,-2.193]
const NC_HG_HBA  = Float32[-0.00828,-0.00828,-0.00828,-0.00828,-0.54648,-0.00828,-0.78296,-0.58984,-0.00828,-0.00828,-0.54984,-0.00828]
# htdbh.f SISKIY(sp,1:3) = P2,P3,P4
const NC_HD_P2 = Float32[523.0987,819.8690,523.0987,604.8450,160.6821,1530.3300,48.6795,679.1972,202.8860,1348.0419,679.1972,595.1068]
const NC_HD_P3 = Float32[5.7243,6.4531,5.7243,5.9835,4.1677,7.0811,8.9420,5.5698,8.7469,7.0463,5.5698,5.8103]
const NC_HD_P4 = Float32[-0.4109,-0.3434,-0.4109,-0.3789,-0.4954,-0.2544,-1.4832,-0.3074,-0.8317,-0.3076,-0.3074,-0.3821]

"nc/htgr5.f — small-tree height increment (3 methods by species). CR is the crown ratio on the ICR/10 (0–10)
scale, NOT the 0–1 proportion (regent.f:183 `CR=ICR(I)/10.0`). Final floor `IF(HTGR.LE.0.0)HTGR=0.01`."
@inline function nc_htgr5(sp::Int, ssite::Float32, baa::Float32, relht::Float32, cr::Float32, h::Float32)
    im = NC_HG_IMETH[sp]
    htgr = if im == 2
        b = baa <= 5.0f0 ? 5.0f0 : baa
        exp(NC_HG_HCON[sp] + NC_HG_HBA[sp] * log(b))
    elseif im == 1
        NC_HG_HCON[sp] + relht * 4.292f0 + 0.0566f0 * cr * cr +
            0.1699f0 * h + NC_HG_HBA[sp] * baa + 0.00768f0 * ssite
    else                                              # sp12 redwood — site-age curve
        htmax = 2.242202f0 * ssite
        if htmax - h <= 1.0f0
            0.0f0
        else
            age1 = (1.0f0 / -0.010742f0) * log(1.0f0 - (h / 2.242202f0 / ssite)^(1.0f0 / 0.919076f0))
            age2 = age1 + 5.0f0
            h1 = 2.242202f0 * ssite * (1.0f0 - exp(-0.010742f0 * age1))^0.919076f0
            h2 = 2.242202f0 * ssite * (1.0f0 - exp(-0.010742f0 * age2))^0.919076f0
            h2 - h1
        end
    end
    htgr <= 0.0f0 && (htgr = 0.01f0)                  # htgr5.f: IF(HTGR.LE.0.0)HTGR=0.01
    return htgr
end

"nc/htdbh.f — Curtis-Arney HT-DBH (SISKIY, all forests). mode 1: H→D."
@inline function nc_htdbh_d(sp::Int, h::Float32)
    p2 = NC_HD_P2[sp]; p3 = NC_HD_P3[sp]; p4 = NC_HD_P4[sp]
    hlim = 4.5f0 + p2 * exp(-p3 * 3.0f0^p4)            # H at D=3 (curve→linear break)
    if h > hlim
        return exp(log((log(h - 4.5f0) - log(p2)) / (-p3)) / p4)
    else
        return ((h - 4.51f0) * 2.7f0) / (hlim - 4.51f0) + 0.3f0
    end
end

function small_tree_growth!(s::StandState, stash, ::Klamath; fint::Float32 = 10.0f0)
    p, t = s.plot, s.trees
    dens = s.density; sd = s.coef.species
    ba = p.basal_area; avh = p.avg_height
    scale = fint / NC_REGYR
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    @inbounds for i in 1:t.n
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]
        d >= NC_ST_XMAX[sp] && continue
        h = t.height[i]
        ssite = p.sp_site_index[sp]
        pt_i = Int(t.plot_id[i])
        tpccf = (1 <= pt_i <= length(dens.point_ccf)) ? dens.point_ccf[pt_i] : 0f0
        relht = (h > 0f0 && avh > 0f0) ? h / avh : 1f0
        tpccf <= 75.0f0 && (relht = 1.0f0 - ((relht - 1.0f0) / 75.0f0) * tpccf)
        relht > 1.5f0 && (relht = 1.5f0)
        cr = Float32(t.crown_pct[i]) * 0.1f0                   # regent.f:183 CR=ICR(I)/10.0 (0–10 scale, NOT /100)
        htgr = nc_htgr5(sp, ssite, ba, relht, cr, h) * scale   # CON=1, XRHMLT=1
        # height: XWT blend with the large-tree HTG (already in t.ht_growth[i])
        xmn = NC_ST_XMIN[sp]; xmx = NC_ST_XMAX[sp]
        xwt = d <= xmn ? 0f0 : (d - xmn) / (xmx - xmn)
        lthg = t.ht_growth[i]
        htg = sp == 12 ? ((htgr + lthg) / 2f0) * (1f0 - xwt) + xwt * lthg :
                         htgr * (1f0 - xwt) + xwt * lthg
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        t.ht_growth[i] = htg
        # small-tree DG (D < DGMIN): HT-DBH (SISKIY) DK/DKK
        if d < NC_ST_DGMIN[sp]
            hk = h + htg
            if hk <= 4.5f0
                t.diam_growth[i] = 0f0
            else
                dk = nc_htdbh_d(sp, hk)
                dkk = h <= 4.5f0 ? d : nc_htdbh_d(sp, h)
                bark = bark_ratio(bark_a, bark_b, sp, d)
                dgsm = (dk < 0f0 || dkk < 0f0) ? htg * 0.2f0 * bark : (dk - dkk) * bark
                dgsm < 0f0 && (dgsm = 0f0)
                # blend small-tree DG with large-tree DG via XWT
                dglt = t.diam_growth[i]
                dg = dgsm * (1f0 - xwt) + xwt * dglt
                (d + dg) < NC_ST_DIAM[sp] && (dg = NC_ST_DIAM[sp] - d)
                t.diam_growth[i] = dg
            end
        end
    end
    return s
end
