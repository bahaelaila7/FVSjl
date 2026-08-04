# =============================================================================
# height_growth.jl (britishcolumbia) — BC large-tree height-growth engine hook (canada/bc/htgf.f).
#
# CHUNK 4. `height_growth!(s, ::BritishColumbia; scale)` fills trees.ht_growth per live tree,
# mirroring dgf!/IE height_growth!. V3 regime (LV2ATV=.FALSE., zones ICH/IDF/SBS/SBPS):
#   HTG(i) = V3HTG(sp, HT, DBH, DG) · Y,  Y = SCALE · XHMULT(sp) · MISHGF   (htgf.f:1612-1635)
# V3HTG (bc_v3_htg) already folds MtoFT·MLT (MLT=HCOR2). MISHGF (mistletoe) = 1 for BC.
# HTG uses chunk-3 DG output (t.diam_growth). ⚠ V2 regime (other zones, htgf.f:1618-1622 exp-form) TODO.
# =============================================================================

"""height_growth!(s, ::BritishColumbia; scale) — per-tree periodic height increment → trees.ht_growth (V3)."""
function height_growth!(s::StandState, ::BritishColumbia; scale::Float32 = 1.0f0)
    p, t, ctl = s.plot, s.trees, s.control
    zone, series = bc_stand_zone(s)
    bc_lv2atv(zone) && return bc_v2_height_growth!(s; scale = scale)   # V2 (ESSF/MS/PP): exp-form HTG
    cur_year = current_cycle_year(s)
    nsp = nspecies(BritishColumbia())
    ip = zeros(Int, nsp)
    @inbounds for sp in 1:nsp
        ip[sp] = bc_resolve_lthg(sp, series, zone)
    end
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i])
        (sp < 1 || sp > nsp || ip[sp] < 1) && continue          # LTH%FIT false ⇒ HTG 0
        d = t.dbh[i]; hti = t.height[i]
        (d <= 0f0 || hti <= 0f0) && continue
        dg = t.diam_growth[i]
        mlt = ctl.htg_cor2[sp] > 0f0 ? ctl.htg_cor2[sp] : 1f0    # LTH%MLT = HCOR2 (htgf.f:1810)
        v3  = bc_v3_htg(sp, ip[sp], hti, d, dg; mlt = mlt)
        xht = active_multiplier(ctl, :htg, sp, cur_year)         # XHMULT (HTGMULT keyword; MISHGF=1 for BC)
        t.ht_growth[i] = v3 * scale * xht
    end
    return s
end

# --- V2 (LV2ATV) height growth (htgf.f:1618-1622) — exp-form: HTG = (exp(CON + HDGCOF·ln(DG)) + BIAS)·Y. ---
# CON = HTCON(sp) + H2COF·HTI² + HGLD(sp)·ln(D) + HGLH·ln(HTI); HTCON(sp) = HGHC(IHT) + HGSC(sp) [+ln(HCOR2)].
# IHT = htgf.f's OWN 1-D MAPHAB(ITYPE) [ITYPE=4 → IHT=2]; H2COF/HDGCOF/HGHC indexed by IHT (habitat class 1-8).
# Y = SCALE·XHMULT·MISHGF (SCALE=FINT/YR, YR=10; MISHGF=1 for BC). Imperial (HTI/HTG in ft). Coeffs = htgf.f DATA.
const BC_HTG_MAPHAB = Int[1,1,2,2,2,2,2,2,2,3,3,4,5,6,7,7,7,7,4,4,1,4,4,8,8,8,1,1,1,1]  # MAPHAB(ITYPE)→IHT
const BC_HGHC  = Float32[2.03035,1.72222,1.19728,1.81759,2.14781,1.76998,2.21104,1.74090]  # HGHC(IHT)
const BC_HGLDD = Float32[0.62144,1.02372,0.85493,0.75756,0.46238,0.49643,0.37042,0.34003]  # HDGCOF(IHT)
const BC_HGH2  = Float32[-13.358f-5,-3.809f-5,-3.715f-5,-2.607f-5,-5.200f-5,-1.605f-5,-3.631f-5,-4.460f-5]  # H2COF(IHT)
const BC_HGSC  = Float32[-.5342,.1433,.1641,-.6458,-.6959,-.9941,-.6004,.2089,-.5478,.7316,-.9941,-.9941,-.9941,.1641,-.9941]
const BC_HGLD  = Float32[-.04935,-.3899,-.4574,-.09775,-.1555,-.1219,-.2454,-.5720,-.1997,-.5657,-.1219,-.1219,-.1219,-.4574,-.1219]
const BC_HTG_BIAS = 0.4809f0
const BC_HTG_HGLH = 0.23315f0

"""BC V2 height growth (htgf.f LV2ATV) — exp(CON + HDGCOF·ln(DG)) + BIAS, scaled by Y. Imperial ft."""
function bc_v2_height_growth!(s::StandState; scale::Float32 = 1.0f0)
    t, ctl = s.trees, s.control
    cur_year = current_cycle_year(s)
    iht = BC_HTG_MAPHAB[BC_V2_ITYPE]                          # ITYPE=4 → IHT=2 (htgf's own MAPHAB)
    hghch = BC_HGHC[iht]; h2cof = BC_HGH2[iht]; hdgcof = BC_HGLDD[iht]
    nsp = nspecies(BritishColumbia())
    htcon = zeros(Float32, nsp)
    @inbounds for sp in 1:nsp
        htcon[sp] = hghch + BC_HGSC[sp]
        ctl.htg_cor2[sp] > 0f0 && (htcon[sp] += log(ctl.htg_cor2[sp]))   # LHCOR2 calib (htgf.f:1800)
    end
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; hti = t.height[i]
        (d <= 0f0 || hti <= 0f0) && continue
        dg = t.diam_growth[i]
        con = htcon[sp] + h2cof*hti*hti + BC_HGLD[sp]*log(d) + BC_HTG_HGLH*log(hti)
        hexp = dg > 0f0 ? exp(con + hdgcof*log(dg)) : 0f0     # DG>0 (ln); DG≤0 ⇒ exp term 0 (ln→-Inf)
        xht = active_multiplier(ctl, :htg, sp, cur_year)      # XHMULT (MISHGF=1 for BC)
        t.ht_growth[i] = (hexp + BC_HTG_BIAS) * scale * xht    # Y = SCALE·XHT·MISHGF
    end
    return s
end
