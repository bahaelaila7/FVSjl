# =============================================================================
# regent.jl (westcascades) — WC small-tree height + diameter growth. Chunk 6.
#
# Ports wc/regent.f (REGENT driver) + vwc/smhgdg.f (SMHGDG Gould-Harrington 2011 small-tree
# height/diameter model) + wc/htdbh.f (6-forest Curtis-Arney HT-DBH) + wc/dgbnd.f (DG cap).
#
# small_tree_growth!(s, stash, ::WestCascades) overrides HTG/DG for trees below XMAX, blended with
# the large-tree prediction by XWT=(D−XMIN)/(XMAX−XMIN). Unlike BM's SMHTGF potential-height model,
# WC calls SMHGDG TWICE per tree (each returns a 5-yr HG5/DG5 increment ⇒ the 10-yr REGYR increment),
# then HTGR=(HTGR+ZZRAN·0.1)·XRHGRO·SCALE·CON·WK4. CON=RHCON(=1)·exp(HCOR). Small-tree DG (D<DGMIN=3,
# RW 7): DG=DGR·SCALE·WK4·BARK, DDS round-trip (identity at FINT=10), DIAM floor, then DGBND.
#
# MEASURED vs FVSwc_g16 DEBUG REGENT (LSTART calibration dump, wct01): SMHGDG bit-exact — e.g. tree
# IT=13 (WF, D=0.1,H=3.0,PTBA=128.31,PTBAL=128.26,CR=0.65,SI=43.60,AVHT=0.5·AVH) ⇒ DGS=0.294894,
# HG5=1.35583, matching the live LEAVING-SMHGDG dump. HCOR≡0 (all 39 scale factors 1.00 on wct01).
# =============================================================================

const _WC_RG_REGYR = 10.0f0

# wc/regent.f DATA (species 1..39: SF WF GF AF RF __ NF YC IC ES LP JP SP WP PP DF RW RC WH MH BM RA
#   WA PB GC AS CW WO WJ LL WB KP PY DG HT CH WI __ OT).
const WC_RG_DGMAX = fill(5.0f0, 39)                                   # DATA DGMAX/39*5.0/
const WC_RG_XMAX  = Float32[fill(4f0,10)...,3f0,fill(4f0,5)...,10f0,fill(4f0,22)...]  # 10*4,3,5*4,10,22*4
const WC_RG_XMIN  = Float32[fill(2f0,10)...,1f0,fill(2f0,28)...]      # 10*2,1,28*2
const WC_RG_DGMIN = Float32[fill(3f0,16)...,7f0,fill(3f0,22)...]     # 16*3,7,22*3
const WC_RG_DIAM  = Float32[0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.2,0.2,0.3,  # 7*0.3,2*0.2,0.3,
                            0.4,0.4,0.4,0.4,0.5,0.3,                   #   4*0.4,0.5,0.3,
                            0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,  # 13*0.2,
                            0.3,0.4,0.4,0.2,0.2,0.2,0.2,0.2,0.2,0.2]   # 0.3,2*0.4,7*0.2

# vwc/smhgdg.f DMAX/BETA (39 species) + ALPHA (10 coefficients × 39 species).
const WC_SMH_DMAX = Float32[1.704,1.496,1.639,1.196,1.515,3.396,2.939,1.54,1.683,1.885,1.653,1.798,2.474,2.474,1.798,5.373,2.849,2.79,3.419,1.383,3.094,3.094,2.011,2.166,3.094,2.475,3.713,0.986,1.219,0.623,0.807,0.586,0.86,1.003,1.89,2.166,2.166,0,5.373]
const WC_SMH_BETA = Float32[0.2474,0.2175,0.1797,0.2056,0.2168,0.2168,0.2822,0.2168,0.2815,0.1704,0.1682,0.2168,0.2168,0.2168,0.2369,0.1635,0.1727,0.1829,0.1727,0.3029,0.2168,0.2168,0.2168,0.2168,0.2168,0.2168,0.2168,0.2168,0.2168,0.2168,0.2168,0.1682,0.2168,0.2168,0.2168,0.2168,0.2168,0,0.1635]
# ALPHA[sp] = (a0..a9); indexed a[k+1] for ALPHA(k,ISPC), k=0..9.
const WC_SMH_ALPHA = NTuple{10,Float32}[
  (2.9445,0,0,0.0068,0,0,-0.1895,0,-1.4049,-0.0168),        # 1  SF
  (1.7536,0,0.2928,0.0009,0,-0.0446,-2.0349,0,-1.3839,-0.0033), # 2 WF
  (2.3571,0.0052,0,0.0006,0,-0.4269,-1.2219,0,0,-0.017),    # 3  GF
  (2.5839,0,0.041,0.002,0,-0.0152,-2.206,0,-0.5915,-0.0009),# 4  AF
  (2.4743,0,0,0.0032,0,-0.8934,-2.2709,0,-1.069,0),         # 5  RF
  (3.8205,0,0.0523,0.0051,0,-0.4102,-1.6968,0,-1.4001,-0.0109), # 6 __
  (0.3376,0,0,0.0101,0,0,0,0,0,-0.0043),                    # 7  NF
  (-2.0216,0.0063,0,0,0.7175,0,0,0,0,0),                    # 8  YC
  (0.5996,0,0,0.008,0,0,0,-1.0479,0,0),                     # 9  IC
  (0.0452,0.008,0,0.0071,0,0,0,0,0,0),                      # 10 ES
  (1.74,0,0.3718,0.0027,0,-0.1712,-2.1359,0,-0.7266,-0.0074), # 11 LP
  (1.8451,0,0,0.0167,0,-1.4737,0,0,-0.4103,-0.0112),        # 12 JP
  (3.8085,0,0,0.0023,0,-0.4265,-2.0913,0,-1.3932,-0.0093),  # 13 SP
  (3.8085,0,0,0.0023,0,-0.4265,-2.0913,0,-1.3932,-0.0093),  # 14 WP
  (1.8451,0,0,0.0167,0,-1.4737,0,0,-0.4103,-0.0112),        # 15 PP
  (2.4473,0,0,0.0098,0,-0.429,-0.171,0,-0.1879,-0.011),     # 16 DF
  (2.9527,0,0,0.0066,0,0,-0.4734,0,-0.7394,-0.0207),        # 17 RW
  (1.6815,0,0,0.0068,0,0,0,0,-0.6049,-0.0121),              # 18 RC
  (2.9527,0,0,0.0066,0,0,-0.4734,0,-0.7394,-0.0207),        # 19 WH
  (2.6762,0.0024,0,0.0006,0,-0.4309,-1.6205,0,-0.593,-0.0051), # 20 MH
  (-1.2421,0.0124,0,0,0.4161,0,0,0,0,0),                    # 21 BM
  (1.4593,0,0,0.0085,0,-0.6,0,0,-1.228,0),                  # 22 RA
  (-1.19,0.0158,0,0,0.66,0,0,0,0,0),                        # 23 WA
  (-1.2421,0.0124,0,0,0.7813,0,0,0,0,0),                    # 24 PB
  (-1.2421,0.0124,0,0,0.6382,0,0,0,0,0),                    # 25 GC
  (-1.2421,0.0124,0,0,0.6013,0,0,0,0,0),                    # 26 AS
  (-1.2421,0.0124,0,0,0.6013,0,0,0,0,0),                    # 27 CW
  (-2.191,0,0,0,0.7191,-3.1321,0,0,0,0),                    # 28 WO
  (0.3755,0.012,0,0,0,0,0,0,0,0),                           # 29 WJ
  (1.0527,0,0.358,0.0019,0,0,-0.6008,0,-0.7451,-0.0101),    # 30 LL
  (2.4949,0,0,0.0049,0,-0.2085,-1.7001,0,-0.7952,-0.0177),  # 31 WB
  (-0.8085,0,0.5001,0,0,0,0,0,0,-0.0081),                   # 32 KP
  (1.5156,0,0,0.0012,0,0,-0.5478,0,-0.6123,0),              # 33 PY
  (-3.8345,0,0,0,1.0701,0,0,0,0,0),                         # 34 DG
  (3.5521,0,0,0.0002,0,0,-0.5932,0,-0.5029,-0.0038),        # 35 HT
  (-1.2421,0.0124,0,0,0.7312,0,0,0,0,0),                    # 36 CH
  (-1.2421,0.0124,0,0,0.6598,0,0,0,0,0),                    # 37 WI
  (0,0,0,0,0,0,0,0,0,0),                                    # 38 __
  (2.4473,0,0,0.0098,0,-0.3575,-0.171,0,-0.1879,-0.011),    # 39 OT
]

# wc/htdbh.f — 6-forest × 39-species Curtis-Arney HT-DBH P2/P3/P4 (data/westcascades/htdbh_coeffs_wc.csv).
# forest_ifor 1..6 = GIFFPC MTBSNQ MTHOOD ROGRIV UMPQUA WILLAM. JFOR→ifor: 1→1,2→2,3/7→3,4/10→4,5/9→5,6/8→6.
let
    P2 = zeros(Float32, 6, 39); P3 = zeros(Float32, 6, 39); P4 = zeros(Float32, 6, 39)
    for l in readlines(joinpath(WC_DATADIR, "htdbh_coeffs_wc.csv"))[2:end]
        f = split(strip(l), ','); (isempty(f) || isempty(f[1])) && continue
        fi = parse(Int, f[1]); sp = parse(Int, f[2])
        P2[fi, sp] = parse(Float32, f[3]); P3[fi, sp] = parse(Float32, f[4]); P4[fi, sp] = parse(Float32, f[5])
    end
    global const WC_HTDBH_P2 = P2; global const WC_HTDBH_P3 = P3; global const WC_HTDBH_P4 = P4
end

# JFOR (forkod.f index 1..11: 603,605,606,610,615,618,708,709,710,711,613) → htdbh table 1..6.
# 1(GP)→1, 2(MBS)→2, 3(Hood)/7(Salem)→3, 4(Rogue)/10(Medford)→4, 5(Umpqua)/9(Roseburg)→5, 6(Willam)/8(Eugene)→6.
@inline function _wc_htdbh_ifor(jfor::Int)::Int
    jfor == 1 ? 1 : jfor == 2 ? 2 : (jfor == 3 || jfor == 7) ? 3 :
    (jfor == 4 || jfor == 10) ? 4 : (jfor == 5 || jfor == 9) ? 5 : 6
end

"wc/htdbh.f MODE=0 (D→H): dub a missing height for species `sp` at DBH `d` on forest-table `ifor`."
@inline function wc_htdbh_height(ifor::Int, sp::Int, d::Float32)::Float32
    p2 = WC_HTDBH_P2[ifor, sp]; p3 = WC_HTDBH_P3[ifor, sp]; p4 = WC_HTDBH_P4[ifor, sp]
    if d >= 3.0f0
        return 4.5f0 + p2 * exp(-1f0 * p3 * d^p4)
    else
        return ((4.5f0 + p2 * exp(-1f0 * p3 * 3.0f0^p4) - 4.51f0) * (d - 0.3f0) / 2.7f0) + 4.51f0
    end
end

"wc/htdbh.f MODE=1 (H→D): estimate DBH from height `h` (used by REGENT for redwood, sp 17)."
@inline function wc_htdbh_dbh(ifor::Int, sp::Int, h::Float32)::Float32
    p2 = WC_HTDBH_P2[ifor, sp]; p3 = WC_HTDBH_P3[ifor, sp]; p4 = WC_HTDBH_P4[ifor, sp]
    hat3 = 4.5f0 + p2 * exp(-1f0 * p3 * 3.0f0^p4)
    if h >= hat3
        return exp(log((log(h - 4.5f0) - log(p2)) / (-1f0 * p3)) * (1f0 / p4))
    else
        return (((h - 4.51f0) * 2.7f0) / (4.5f0 + p2 * exp(-1f0 * p3 * 3.0f0^p4) - 4.51f0)) + 0.3f0
    end
end

"wc/dgbnd.f — cap the diameter increment (Oliver-Cochran DF envelope; redwood passes through)."
@inline function wc_dgbnd(sp::Int, dbh::Float32, ddg::Float32, sizcap1::Float32, sizcap3::Float32)::Float32
    if sp != 17
        temdbh = dbh > 150f0 ? 150f0 : dbh
        dgmax = 7.92f0 * exp(-0.03f0 * temdbh)
        ddg > dgmax && (ddg = dgmax)
        ddg < 0f0 && (ddg = 0f0)
    end
    if (dbh + ddg) > sizcap1 && sizcap3 < 1.5f0
        ddg = sizcap1 - dbh
        ddg < 0.01f0 && (ddg = 0.01f0)
    end
    return ddg
end

# vwc/smhgdg.f — one 5-yr SMHGDG call for species `sp` (MODE=1, IT>0), returning (HG5, DG5). Inputs:
#   h,d = current height/diameter; cr = crown fraction; ptbal = PTBALT; ptba = PTBAA; si = SITEAR;
#   avht = (5/FINT)·AVH+((FINT−5)/FINT)·ATAVH (the caller resolves this). Redwood (17) is height-only here.
@inline function wc_smhgdg(sp::Int, h::Float32, d::Float32, cr::Float32, ptbal::Float32,
                           ptba::Float32, si::Float32, avht::Float32)
    if sp == 17                                             # redwood — Chapman-Richards site curve, DG5=0
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
    sp == 16 && (si = 5.21486f0 + 0.66486f0 * si)          # WC Douglas-fir: Curtis→King SI
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

@inline function _wc_rg_stash!(stash, t, i::Int)
    if stash !== nothing && !isempty(stash.dgU) && i <= length(stash.dgU)
        stash.dgU[i] = t.diam_growth[i]; stash.dgL[i] = t.diam_growth[i]
        stash.htgU[i] = t.ht_growth[i]; stash.htgL[i] = t.ht_growth[i]
        !isempty(stash.is_small) && (stash.is_small[i] = true)
    end
end

function small_tree_growth!(s::StandState, stash, ::WestCascades; fint::Float32 = 10.0f0)
    p, t, c = s.plot, s.trees, s.calib
    n = t.n; n == 0 && return s
    sd = s.coef.species; dens = s.density
    avh = p.avg_height; dgsd = s.control.dg_sd
    ifor = _wc_htdbh_ifor(Int(p.forest_idx))
    scale = fint / _WC_RG_REGYR                             # SCALE = FINT/REGYR
    scale2 = _WC_RG_REGYR / fint                            # SCALE2 = YR/FNT  (YR=10)
    # ATAVH = AVH during the cycling growth call (grincr.f:318 sets ATAVH=AVH before TREGRO/REGENT) ⇒
    # SMHGDG's AVHT = (5/FINT)·AVH+((FINT−5)/FINT)·ATAVH = AVH. (The LSTART calibration uses ATAVH=0.)
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
        con = exp(c.htg_cor_small[sp])                      # RHCON(=1)·exp(HCOR)
        wk4 = t.htimlt[i]                                   # per-tree height multiplier (1 for inventory)
        # --- SMHGDG twice (each 5-yr) ⇒ 10-yr HTGR/DGR (regent.f:234-255) ---
        hg1, dg1 = wc_smhgdg(sp, h, d, cr, ptbal, ptba, si, avht)
        hk = h + hg1; dk = d + dg1
        hg2, dg2 = wc_smhgdg(sp, hk, dk, cr, ptbal, ptba, si, avht)
        htgr = hg1 + hg2; dgr = dg1 + dg2
        # --- ZZRAN reject-loop (regent.f:257-260); DGSD=1.7 ≥ 1 ⇒ draw ---
        zzran = 0.0f0
        if dgsd >= 1.0f0
            while true
                zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                (zzran <= 0.5f0 && zzran >= -2.0f0) && break
            end
        end
        htgr = (htgr + zzran * 0.1f0) * scale * con * wk4    # XRHGRO=1
        htgr < 0.1f0 && (htgr = 0.1f0)
        # --- XWT blend with the large-tree HTG ---
        xmn = WC_RG_XMIN[sp]; xmx = WC_RG_XMAX[sp]
        xwt = d <= xmn ? 0.0f0 : (d - xmn) / (xmx - xmn)
        if sp == 17
            lthg = t.ht_growth[i]
            htgr = (htgr + lthg) / 2.0f0
            t.ht_growth[i] = htgr * (1.0f0 - xwt) + xwt * lthg   # regent.f:282 (pre-final RW blend)
        end
        htg = htgr * (1.0f0 - xwt) + xwt * t.ht_growth[i]
        htg < 0.1f0 && (htg = 0.1f0)
        cap = s.control.sp_size_cap[sp, 4]
        if (h + htg) > cap
            htg = cap - h; htg < 0.1f0 && (htg = 0.1f0)
        end
        t.ht_growth[i] = htg
        # --- small-tree DG (regent.f:304-437): only D<DGMIN (large-tree DG kept otherwise) ---
        if d >= WC_RG_DGMIN[sp]
            _wc_rg_stash!(stash, t, i); continue
        end
        hkk = h + htg
        if hkk < 4.5f0
            t.diam_growth[i] = 0.0f0
            t.dbh[i] = d + 0.001f0 * hkk                     # DBH(K)=D+0.001·HK (regent.f:317)
            _wc_rg_stash!(stash, t, i); continue
        end
        bark = wc_bratio(sd, sp, d)
        local dgk::Float32
        if sp == 17                                          # redwood — Curtis-Arney DBH from HT
            dk2 = wc_htdbh_dbh(ifor, sp, hkk)
            dkk = h <= 4.5f0 ? d : wc_htdbh_dbh(ifor, sp, h)
            xdwt = d <= xmn ? 0.0f0 : (d - xmn) / (7.0f0 - xmn)
            dgsm = (dk2 - dkk) * bark; dgsm < 0.0f0 && (dgsm = 0.0f0)
            dds = dgsm * (2.0f0 * bark * d + dgsm) * scale2
            dgsm = sqrt((d * bark)^2 + dds) - bark * d
            dgk = dgsm * (1.0f0 - xdwt) + t.diam_growth[i] * xdwt
        else
            dgk = dgr * scale * wk4                           # DG(K)=DGR·SCALE·WK4
            if d < 0.0f0 || dgk < 0.0f0
                dgk = htg * 0.2f0 * bark                      # XRDGRO=1
            else
                dgk = dgk * bark                             # XRDGRO=1
            end
            dgk < 0.0f0 && (dgk = 0.1f0)
            dgmx = WC_RG_DGMAX[sp] * scale
            dgk > dgmx && (dgk = dgmx)
            dds = dgk * (2.0f0 * bark * d + dgk) * scale2     # 10-yr DDS (identity at FINT=10)
            dgk = sqrt((d * bark)^2 + dds) - bark * d
        end
        (t.dbh[i] + dgk) < WC_RG_DIAM[sp] && (dgk = WC_RG_DIAM[sp] - t.dbh[i])
        dgk = wc_dgbnd(sp, t.dbh[i], dgk, s.control.sp_size_cap[sp, 1], s.control.sp_size_cap[sp, 3])
        t.diam_growth[i] = dgk
        _wc_rg_stash!(stash, t, i)
    end
    return s
end

# regenerate!(::WestCascades) — WC has no cyc0 regeneration hook on the pure-growth path (ESTAB is
# a separate keyword-gated driver). Provide a no-op so the abstract dispatch resolves.
regenerate!(s::StandState, ::WestCascades; kwargs...) = s
