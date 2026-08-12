# =============================================================================
# regent.jl (bluemountains) — BM small-tree height + diameter growth (bm/regent.f + bm/smhtgf.f). Chunk 6.
#
# small_tree_growth!(s, stash, ::BlueMountains) overrides HTG/DG for trees below XMAX, blended with the
# large-tree prediction by XWT=(D−XMIN)/(XMAX−XMIN). Height: HTGR = POTHTG·PCTRED·VIGOR·CON, where POTHTG
# = bm_smhtgf(sp,SI,H,DTIME=10) (per-species small-tree height curves); LM(12)=SI/5; AS(15) Sheppard.
# + ZZRAN·0.1 (reject until ∈[−2,0.5]), ·XRHGRO(=1)·SCALE(=FINT/REGYR). CON=RHCON(=1)·exp(HCOR).
# Small-tree DG (D<BKPT=3, WJ 99): HK=H+HTG; DK/DKK ht-dbh curve → DGMX clamp → DDS → DG, DIAM floor.
# =============================================================================

# bm/regent.f DATA (species order WP WL DF GF MH WJ LP ES AF PP WB LM PY YC AS CW OS OH).
const BM_RG_DGMAX = Float32[2.8,2.8,2.4,3.6,2.5,2.0,3.5,3.6,3.6,2.8,2.8,2.8,5.0,5.0,2.5,5.0,2.8,5.0]
const BM_RG_XMAX  = Float32[3,2,4,4,2,99,4,4,4,5,3,4,4,4,4,4,5,4]
const BM_RG_XMIN  = Float32[2,1,2,2,1,90,2,2,2,1,1.5,2,2,2,2,2,1,2]
const BM_RG_DIAM  = Float32[0.4,0.3,0.3,0.3,0.2,0.3,0.4,0.3,0.3,0.5,0.4,0.4,0.2,0.2,0.2,0.2,0.5,0.2]
const BM_RG_AB    = Float32[1.11436, -0.011493, 0.43012f-4, -0.72221f-7, 0.5607f-10, -0.1641f-13]
const _BM_RG_REGYR = 10.0f0

# bm/smhtgf.f — small-tree potential height growth (POTHTG) over DTIME years for species sp on site si.
function bm_smhtgf(sp::Int, si::Float32, h::Float32, dtime::Float32)::Float32
    if sp == 1                                            # WP — Chapman-Richards
        c1 = 0.375045f0; c2 = 0.92503f0; c3 = -0.020796f0; c4 = 2.48811f0
        arg = (1.0f0 - (c1 / si * h)^(1f0 / c4)) / c2
        effage = arg > 0f0 ? log(arg) / c3 : 0f0
        agepdt = effage + dtime
        return (si / c1) * (1f0 - c2 * exp(c3 * agepdt))^c4 - (si / c1) * (1f0 - c2 * exp(c3 * effage))^c4
    elseif sp == 2
        return ((-3.9725f0 + 0.50995f0 * si) / (28.1168f0 - 0.05661f0 * si)) * dtime
    elseif sp == 3
        return ((2.0f0 + 0.420f0 * si) / (28.5f0 - 0.05f0 * si)) * dtime
    elseif sp == 4
        return ((4.2435f0 + 0.1510f0 * si) / (19.0184f0 - 0.0570f0 * si)) * dtime
    elseif sp == 5
        return ((0.965758f0 + 0.082969f0 * si) / (55.249612f0 - 1.288852f0 * si)) * dtime
    elseif sp == 6
        s = clamp(si, 5.5f0, 75.0f0)
        return (s / 5.0f0) * (s * 1.5f0 - h) / (s * 1.5f0)
    elseif sp == 7
        return (0.02008805f0 * si) * dtime
    elseif sp == 8
        return ((0.09211f0 + 0.208517f0 * si) / (43.358f0 - 0.168166f0 * si)) * dtime
    elseif sp == 9
        return ((6.0f0 + 0.14f0 * si) / (33.882f0 - 0.06588f0 * si)) * dtime
    elseif sp == 10 || sp == 17
        return ((-1.0f0 + 0.32857f0 * si) / (28.0f0 - 0.042857f0 * si)) * dtime
    elseif sp == 11
        return ((0.02008805f0 * si) * dtime) * 1.6f0
    elseif sp == 12
        return 0.5f0
    elseif sp == 15
        return 5.0f0
    else                                                  # PY/YC/CW/OH (13,14,16,18)
        return ((1.47043f0 + 0.23317f0 * si) / (31.56252f0 - 0.05586f0 * si)) * dtime
    end
end

# bm/essubh.f — subsequent/planted-tree base height. HHT = SMHTGF(sp, HHT, H, MODE=0, DTIME=AGE): the
# small-tree height-at-total-age curve. NO EMSQR/DILATE/ELEV (bm/essubh.f explicitly discards them, lines
# 47-49). MODE=0 sets EFFAGE=0 (start-from-zero), so only the non-linear WP(sp1) differs from the regent
# (MODE=1) call; all linear species are coef·AGE with H unused. Deterministic ⇒ bit-exact-portable.
function bm_essubh_hht(sp::Int, si::Float32, age::Float32)::Float32
    age <= 0f0 && return 0f0                                # SMHTGF: DTIME≤0 → HHT=0 (bm/smhtgf.f:76)
    if sp == 1                                              # WP — MODE=0: EFFAGE=0 (not H-derived)
        c1 = 0.375045f0; c2 = 0.92503f0; c3 = -0.020796f0; c4 = 2.48811f0
        return (si / c1) * (1f0 - c2 * exp(c3 * age))^c4 - (si / c1) * (1f0 - c2)^c4
    end
    return bm_smhtgf(sp, si, 0f0, age)                      # linear/fixed species: H unused ⇒ = coef·AGE
end

@inline function _bm_rg_stash!(stash, t, i::Int)
    if stash !== nothing && !isempty(stash.dgU) && i <= length(stash.dgU)
        stash.dgU[i] = t.diam_growth[i]; stash.dgL[i] = t.diam_growth[i]
        stash.htgU[i] = t.ht_growth[i]; stash.htgL[i] = t.ht_growth[i]
        !isempty(stash.is_small) && (stash.is_small[i] = true)
    end
end

function small_tree_growth!(s::StandState, stash, ::BlueMountains; fint::Float32 = 10.0f0)
    p, t, c = s.plot, s.trees, s.calib
    n = t.n; n == 0 && return s
    sd = s.coef.species; slo = sd[:site_lo]; shi = sd[:site_hi]
    relden = p.relative_density; avh = p.avg_height
    dgsd = s.control.dg_sd
    scale = fint / _BM_RG_REGYR                           # SCALE = FINT/REGYR
    # PCTRED (density modifier), stand-level (bm/regent.f:198-201)
    xd = avh * (relden / 100.0f0); xd > 300.0f0 && (xd = 300.0f0)
    ab = BM_RG_AB
    pctred = ab[1] + xd*(ab[2] + xd*(ab[3] + xd*(ab[4] + xd*(ab[5] + xd*ab[6]))))
    pctred > 1.0f0 && (pctred = 1.0f0); pctred < 0.01f0 && (pctred = 0.01f0)
    @inbounds for i in 1:n
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= BM_RG_XMAX[sp] || t.tpa[i] <= 0.0f0) && continue
        h = t.height[i]
        sitear = p.sp_site_index[sp]
        si = sitear; si > shi[sp] && (si = shi[sp]); si <= slo[sp] && (si = slo[sp] + 0.5f0)
        relsi = (si - slo[sp]) / (shi[sp] - slo[sp]); rsimod = 0.5f0 * (1.0f0 + relsi)
        con = exp(c.htg_cor_small[sp])                    # RHCON(=1)·exp(HCOR)
        xcr = Float32(t.crown_pct[i]) / 100.0f0
        vigor = 150.0f0 * xcr^3 * exp(-6.0f0 * xcr) + 0.3f0; vigor > 1.0f0 && (vigor = 1.0f0)
        sp == 6 && (vigor = 1.0f0 - (1.0f0 - vigor) / 3.0f0)   # WJ pinyon (bm/regent.f:277)
        # POTHTG + HTGR
        local htgr::Float32
        if sp == 12
            htgr = (si / 5.0f0) * pctred * vigor * con
        elseif sp == 15                                   # aspen Sheppard (bm/regent.f:294-305)
            age = (h * 2.54f0 * 12.0f0 / 26.9825f0)^(1.0f0 / 1.1752f0)
            hite1 = 26.9825f0 * age^1.1752f0; hite2 = 26.9825f0 * (age + 10.0f0)^1.1752f0
            htgr = (hite2 - hite1) / (2.54f0 * 12.0f0) * rsimod * con * 2.40f0 * 0.75f0
        else
            pothtg = bm_smhtgf(sp, si, h, _BM_RG_REGYR)    # DTIME=TEMT=10
            htgr = pothtg * pctred * vigor * con
        end
        # ZZRAN reject-loop (bm/regent.f:308-310)
        zzran = 0.0f0
        if dgsd >= 1.0f0
            while true
                zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                (zzran <= 0.5f0 && zzran >= -2.0f0) && break
            end
        end
        htgr = (htgr + zzran * 0.1f0) * scale             # XRHGRO=1
        htgr < 0.1f0 && (htgr = 0.1f0)
        # XWT blend with the large-tree HTG
        xmn = BM_RG_XMIN[sp]; xmx = BM_RG_XMAX[sp]
        xwt = d <= xmn ? 0.0f0 : (d - xmn) / (xmx - xmn)
        htg = htgr * (1.0f0 - xwt) + xwt * t.ht_growth[i]; htg < 0.1f0 && (htg = 0.1f0)
        t.ht_growth[i] = htg
        # ---- small-tree DG (bm/regent.f:388-460). BKPT=3 (WJ 99). HK=H+HTG; DK/DKK ht-dbh → DGMX → DDS.
        bkpt = sp == 6 ? 99.0f0 : 3.0f0
        d >= bkpt && (_bm_rg_stash!(stash, t, i); continue)
        hk = h + htg
        bark = bm_bratio(sd, sp, d)
        if hk <= 4.5f0
            t.diam_growth[i] = 0.0f0
        else
            local dk::Float32, dkk::Float32
            if sp == 7                                    # LP — fixed ht-dbh
                dk = -9.8752f0 / (log(hk - 4.5f0) - 4.8656f0) - 1.0f0
                dkk = h <= 4.5f0 ? d : -9.8752f0 / (log(h - 4.5f0) - 4.8656f0) - 1.0f0
            elseif sp == 6                                # WJ — linear site
                dk = (hk - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
                dkk = h < 4.5f0 ? d : (h - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0)
            elseif (!BM_LHTDRG[sp] || c.ht_dbh_iabflg[sp] == 1) && _bm_has_htdbh(Int(p.forest_idx), sp)
                # regent.f:522 — .NOT.LHTDRG OR (LHTDRG & IABFLG==1) ⇒ HTDBH (Curtis-Arney, forest-dependent).
                # BM conifers DF/GF/ES/WL have LHTDRG=false, so they use HTDBH — NOT the Wykoff AX/BX. But ONLY when
                # the species actually has Curtis-Arney coeffs (P2>0): AS(15)/LM(12)/WB(11) have P2=0 (bm/regent.f
                # routes them to the AX/BX / SMDGF branches, never HTDBH) — without this guard bm_htdbh hits
                # log(P2=0)=-Inf ⇒ DomainError crash on real-FIA aspen/limber-pine stands.
                ifor = Int(p.forest_idx)
                dk = bm_htdbh(ifor, sp, hk)
                dkk = h <= 4.5f0 ? d : bm_htdbh(ifor, sp, h)
            else                                          # AX/BX (bm/regent.f CASE 1:5,8:9,12,15 = incl AS/LM)
                bx = sd[:ht2][sp]; ax = c.ht_dbh_aa[sp]
                dk = bx / (log(hk - 4.5f0) - ax) - 1.0f0
                dkk = h <= 4.5f0 ? d : bx / (log(h - 4.5f0) - ax) - 1.0f0
            end
            dgk = (dk - dkk) * bark                        # XRDGRO=1
            dgk < 0.0f0 && (dgk = 0.0f0)
            dgmx = BM_RG_DGMAX[sp] * scale
            sp == 11 && (dgmx = fint * 0.2f0)              # WB (bm/regent.f:239)
            dgk > dgmx && (dgk = dgmx)
            scale2 = _BM_RG_REGYR / fint                   # YR/FINT
            dds = dgk * (2.0f0 * bark * d + dgk) * scale2
            dgk = sqrt((d * bark)^2 + dds) - bark * d
            (d + dgk) < BM_RG_DIAM[sp] && (dgk = BM_RG_DIAM[sp] - d)
            t.diam_growth[i] = dgk
        end
        _bm_rg_stash!(stash, t, i)
    end
    return s
end

# bm/esgent.f (CALL REGENT(.TRUE.,ITRNIN)) — grow the JUST-ESTABLISHED regen IN its birth cycle. BM was OMITTED
# from the esgent dispatch (simulate.jl had CR/TT/EM/UT/CI), so planted/established BM seedlings never got their
# first-cycle height growth (BM BARE-PLANT: persistent TopHt lag ~5 ft). Same class as EM #137 / UT #184 / CI #185.
# Mirrors small_tree_growth!'s POTHTG/PCTRED/VIGOR/CON height + DK/DKK DBH over the birth subperiod (subyr=FINT−
# GENTIM=5), applying HT/DBH directly (esgent.f HT(I)=HT(I)+HTG(I)·WK4). Gated to the new records nstart+1:n.
function bm_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0)
    p, t, c = s.plot, s.trees, s.calib
    nstart >= t.n && return s
    sd = s.coef.species; slo = sd[:site_lo]; shi = sd[:site_hi]
    relden = p.relative_density; avh = p.avg_height; dgsd = s.control.dg_sd
    gentim = max(fint - 5.0f0, 0.0f0)
    bscale = (fint - gentim) / _BM_RG_REGYR              # birth-cycle fraction (WK4; =0.5 for fint=10)
    xd = avh * (relden / 100.0f0); xd > 300.0f0 && (xd = 300.0f0)
    ab = BM_RG_AB
    pctred = ab[1] + xd*(ab[2] + xd*(ab[3] + xd*(ab[4] + xd*(ab[5] + xd*ab[6]))))
    pctred > 1.0f0 && (pctred = 1.0f0); pctred < 0.01f0 && (pctred = 0.01f0)
    @inbounds for i in (nstart+1):t.n
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= BM_RG_XMAX[sp] || t.tpa[i] <= 0.0f0) && continue
        h = t.height[i]
        sitear = p.sp_site_index[sp]
        si = sitear; si > shi[sp] && (si = shi[sp]); si <= slo[sp] && (si = slo[sp] + 0.5f0)
        relsi = (si - slo[sp]) / (shi[sp] - slo[sp]); rsimod = 0.5f0 * (1.0f0 + relsi)
        con = exp(c.htg_cor_small[sp])
        xcr = Float32(t.crown_pct[i]) / 100.0f0
        vigor = 150.0f0 * xcr^3 * exp(-6.0f0 * xcr) + 0.3f0; vigor > 1.0f0 && (vigor = 1.0f0)
        sp == 6 && (vigor = 1.0f0 - (1.0f0 - vigor) / 3.0f0)
        local htgr::Float32
        if sp == 12
            htgr = (si / 5.0f0) * pctred * vigor * con
        elseif sp == 15
            age = (h * 2.54f0 * 12.0f0 / 26.9825f0)^(1.0f0 / 1.1752f0)
            hite1 = 26.9825f0 * age^1.1752f0; hite2 = 26.9825f0 * (age + 10.0f0)^1.1752f0
            htgr = (hite2 - hite1) / (2.54f0 * 12.0f0) * rsimod * con * 2.40f0 * 0.75f0
        else
            pothtg = bm_smhtgf(sp, si, h, _BM_RG_REGYR)
            htgr = pothtg * pctred * vigor * con
        end
        zzran = 0.0f0
        if dgsd >= 1.0f0
            while true
                zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                (zzran <= 0.5f0 && zzran >= -2.0f0) && break
            end
        end
        htgr = (htgr + zzran * 0.1f0) * bscale           # birth-cycle subperiod (was scale=fint/REGYR)
        htgr < 0.1f0 && (htgr = 0.1f0)
        xmn = BM_RG_XMIN[sp]; xmx = BM_RG_XMAX[sp]
        xwt = d <= xmn ? 0.0f0 : (d - xmn) / (xmx - xmn)
        htg = htgr * (1.0f0 - xwt); htg < 0.1f0 && (htg = 0.1f0)   # new tree: large-tree HTG(K)=0
        hk = h + htg
        t.height[i] = hk; t.ht_growth[i] = htg
        bkpt = sp == 6 ? 99.0f0 : 3.0f0
        (d >= bkpt || hk <= 4.5f0) && continue
        bark = bm_bratio(sd, sp, d)
        local dk::Float32, dkk::Float32
        if sp == 7
            dk = -9.8752f0 / (log(hk - 4.5f0) - 4.8656f0) - 1.0f0
            dkk = h <= 4.5f0 ? d : -9.8752f0 / (log(h - 4.5f0) - 4.8656f0) - 1.0f0
        elseif sp == 6
            dk = (hk - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
            dkk = h < 4.5f0 ? d : (h - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0)
        elseif (!BM_LHTDRG[sp] || c.ht_dbh_iabflg[sp] == 1) && _bm_has_htdbh(Int(p.forest_idx), sp)
            ifor = Int(p.forest_idx)
            dk = bm_htdbh(ifor, sp, hk)
            dkk = h <= 4.5f0 ? d : bm_htdbh(ifor, sp, h)
        else
            bx = sd[:ht2][sp]; ax = c.ht_dbh_aa[sp]
            dk = bx / (log(hk - 4.5f0) - ax) - 1.0f0
            dkk = h <= 4.5f0 ? d : bx / (log(h - 4.5f0) - ax) - 1.0f0
        end
        dgk = (dk - dkk) * bark; dgk < 0.0f0 && (dgk = 0.0f0)
        dgmx = BM_RG_DGMAX[sp] * bscale
        sp == 11 && (dgmx = fint * 0.2f0 * bscale)
        dgk > dgmx && (dgk = dgmx)
        dds = dgk * (2.0f0 * bark * d + dgk)
        dgk = sqrt((d * bark)^2 + dds) - bark * d
        (d + dgk) < BM_RG_DIAM[sp] && (dgk = BM_RG_DIAM[sp] - d)
        dgk > 0.0f0 && (t.dbh[i] = d + dgk; t.diam_growth[i] = dgk)
    end
    return s
end

# --- bm/htdbh.f Curtis-Arney ht-dbh (forest-dependent P2/P3/P4), used for LHTDRG=false species ---
# BM LHTDRG (grinit.f:129,151-154): FALSE for all except WJ(6)/WB(11)/LM(12)/AS(15). regent.f:522 uses
# HTDBH when .NOT.LHTDRG OR (LHTDRG & IABFLG==1); else the Wykoff AX/BX. So BM conifers (DF/GF/ES/WL, LHTDRG=false)
# use HTDBH — jl's regent must too (was using AX/BX for all ⇒ over-grew seedlings).
const BM_LHTDRG = Bool[false,false,false,false,false,true,false,false,false,false,true,true,false,false,true,false,false,false]
let
    path = joinpath(BM_DATADIR, "htdbh_coeffs_bm.csv")
    rows = [split(strip(l), ',') for l in readlines(path)[2:end]]
    P2 = zeros(Float32,4,18); P3 = zeros(Float32,4,18); P4 = zeros(Float32,4,18)
    for r in rows
        fi=round(Int,parse(Float32,r[1])); sp=round(Int,parse(Float32,r[2]))
        P2[fi,sp]=parse(Float32,r[3]); P3[fi,sp]=parse(Float32,r[4]); P4[fi,sp]=parse(Float32,r[5])
    end
    global const BM_HTDBH_P2 = P2; global const BM_HTDBH_P3 = P3; global const BM_HTDBH_P4 = P4
end

# True iff species `sp` has Curtis-Arney HTDBH coefficients (P2>0) for forest `ifor` (clamped 1..4 like bm_htdbh).
# AS(15)/LM(12)/WB(11) have all-zero rows ⇒ false ⇒ the caller must NOT use bm_htdbh (log(P2=0)=-Inf crash).
@inline function _bm_has_htdbh(ifor::Int, sp::Int)::Bool
    (ifor < 1 || ifor > 4) && (ifor = 3)
    return BM_HTDBH_P2[ifor, sp] > 0f0
end

# bm/htdbh.f MODE=1 (HT→DBH), Curtis-Arney with a linear small-tree segment below HAT3 (= _ut_htdbh_dbh form).
@inline function bm_htdbh(ifor::Int, sp::Int, h::Float32)::Float32
    (ifor < 1 || ifor > 4) && (ifor = 3)
    p2 = BM_HTDBH_P2[ifor,sp]; p3 = BM_HTDBH_P3[ifor,sp]; p4 = BM_HTDBH_P4[ifor,sp]
    hat3 = 4.5f0 + p2 * exp(-1f0 * p3 * 3.0f0^p4)
    if h >= hat3
        return exp(log((log(h - 4.5f0) - log(p2)) / (-1f0 * p3)) * (1f0 / p4))
    else
        return ((h - 4.51f0) * 2.7f0) / (4.5f0 + p2 * exp(-1f0 * p3 * 3.0f0^p4) - 4.51f0) + 0.3f0
    end
end
