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
    # TRIPLING (bm/regent.f:288 label 2 … :634-637 `L=L+1; K=ITRN+2*I-2+L; GO TO 2`): on a tripling cycle each
    # small-tree record is grown THREE times — central (L=0), upper copy (L=1, K=ITRN+2I-1 → stash U) and lower
    # copy (L=2, K=ITRN+2I → stash L) — and every pass re-enters label 2, so each draws its OWN ZZRAN (label 3,
    # :357-359) and gets its OWN HTG(K) (+ its own DG(K) when D<BKPT). Drawing once and copying the central values
    # left jl 2 ZZRAN draws per small tree short on every tripling cycle ⇒ the whole main rann! stream desynced for
    # every downstream consumer (DGSCOR, MISTOE spread …). D≥BKPT jumps to 23 (no DG(K)): the copies keep their own
    # large-tree tripled DG (dgU/dgL from diameter_growth!) — only their HTG is REGENT's.
    ntrip = (stash !== nothing && !isempty(stash.htgU)) ? 3 : 1
    # Visit order = FVS DO 30 ISPC / DO 25 I3=ISCT(ISPC,1..2) over IND1 (species-major, lineage order), which fixes
    # which record consumes which ZZRAN draw.
    species_sort!(s)
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    @inbounds for sp_o in 1:MAXSP, i3 in (isct[sp_o, 1] == 0 ? (1:0) : (Int(isct[sp_o, 1]):Int(isct[sp_o, 2])))
        i = Int(ind1[i3])
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= BM_RG_XMAX[sp] || t.tpa[i] <= 0.0f0) && continue
        h = t.height[i]
        sitear = p.sp_site_index[sp]
        # regent.f:230-232 clamps SI to the species SITERANGE [SLO,SHI], but that clamped SI drives ONLY the
        # LM(12) POTHTG=SI/5 and the aspen(15) RELSI. For every SMHTGF species (CASE 1:11,13:14,16:18) SMHTGF
        # IGNORES regent's SI and re-reads the RAW SITEAR(I) itself (smhtgf.f:78, no clamp except its own WJ
        # branch) — so POTHTG uses the UNCLAMPED site index. jl was passing the clamped `si` to bm_smhtgf,
        # over-growing species whose SITEAR falls outside [SLO,SHI] (e.g. LP SITEAR 23.9 clamped up to SLO+0.5=
        # 30.5 ⇒ POTHTG 6.13 vs live 4.80 ⇒ seedling DBH/HT one-directionally high). Pass the raw SITEAR.
        si = sitear; si > shi[sp] && (si = shi[sp]); si <= slo[sp] && (si = slo[sp] + 0.5f0)
        relsi = (si - slo[sp]) / (shi[sp] - slo[sp]); rsimod = 0.5f0 * (1.0f0 + relsi)
        con = exp(c.htg_cor_small[sp])                    # RHCON(=1)·exp(HCOR)
        xcr = Float32(t.crown_pct[i]) / 100.0f0
        vigor = 150.0f0 * xcr^3 * exp(-6.0f0 * xcr) + 0.3f0; vigor > 1.0f0 && (vigor = 1.0f0)
        sp == 6 && (vigor = 1.0f0 - (1.0f0 - vigor) / 3.0f0)   # WJ pinyon (bm/regent.f:277)
        # POTHTG + HTGR (identical for every tripled pass: H/D/VIGOR are the central record's)
        local htgr0::Float32
        if sp == 12
            htgr0 = (si / 5.0f0) * pctred * vigor * con    # LM: regent.f:311 uses the CLAMPED SI
        elseif sp == 15                                   # aspen Sheppard (bm/regent.f:294-305)
            age = (h * 2.54f0 * 12.0f0 / 26.9825f0)^(1.0f0 / 1.1752f0)
            hite1 = 26.9825f0 * age^1.1752f0; hite2 = 26.9825f0 * (age + 10.0f0)^1.1752f0
            htgr0 = (hite2 - hite1) / (2.54f0 * 12.0f0) * rsimod * con * 2.40f0 * 0.75f0
        else
            pothtg = bm_smhtgf(sp, sitear, h, _BM_RG_REGYR)   # SMHTGF reads raw SITEAR (unclamped); DTIME=TEMT=10
            htgr0 = pothtg * pctred * vigor * con
        end
        xmn = BM_RG_XMIN[sp]; xmx = BM_RG_XMAX[sp]
        xwt = d <= xmn ? 0.0f0 : (d - xmn) / (xmx - xmn)
        large_htg = t.ht_growth[i]                        # HTG(K): htgf gives the copies the central TEMHTG
        cap = s.control.sp_size_cap[sp, 4]
        bkpt = sp == 6 ? 99.0f0 : 3.0f0
        small = d < bkpt
        bark = small ? bm_bratio(sd, sp, d) : 0.0f0
        for l in 0:(ntrip - 1)
            # ZZRAN reject-loop (bm/regent.f:357-359) — one per tripled pass
            zzran = 0.0f0
            if dgsd >= 1.0f0
                while true
                    zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                    (zzran <= 0.5f0 && zzran >= -2.0f0) && break
                end
            end
            htgr = (htgr0 + zzran * 0.1f0) * scale        # XRHGRO=1
            htgr < 0.1f0 && (htgr = 0.1f0)
            htg = htgr * (1.0f0 - xwt) + xwt * large_htg; htg < 0.1f0 && (htg = 0.1f0)
            if h + htg > cap                              # regent.f:378-381 SIZCAP(ISPC,4)
                htg = cap - h; htg < 0.1f0 && (htg = 0.1f0)
            end
            # ---- small-tree DG (bm/regent.f:388-460). BKPT=3 (WJ 99). HK=H+HTG; DK/DKK ht-dbh → DGMX → DDS.
            dgk = 0.0f0
            if small
                hk = h + htg
                if hk > 4.5f0
                    local dk::Float32, dkk::Float32
                    # LP(7)/PP(10) have a fixed/linear ht-dbh in the regent.f CASE block (CASE 7 / CASE 10,17), but
                    # LHTDRG(7)=LHTDRG(10)=.FALSE. ⇒ regent.f:519 UNCONDITIONALLY overrides DK/DKK with HTDBH (Curtis-
                    # Arney), so those fixed formulas are DEAD CODE. LP therefore uses HTDBH like DF/GF/ES — jl's fixed
                    # -9.8752 LP branch (dk≈1.14) over-grew LP seedling DBH ~2.7× vs HTDBH (dk≈0.5). Drop it; LP falls
                    # into the HTDBH branch below (all 4 BM forests have LP P2>0).
                    if sp == 6                                    # WJ — linear site
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
                        # bm/regent.f:411-414: AX = IABFLG(ISPC)==1 ? HT1(ISPC) : AA(ISPC). AA is the CRATET-fitted
                        # Wykoff intercept (0 until calibrated); IABFLG=1 (default) means Wykoff calibration did NOT
                        # occur, so the RAW HT1 intercept must be used. The AX/BX branch is reached by AS(15)/LM(12)
                        # (P2=0 ⇒ no HTDBH), which carry IABFLG=1 ⇒ AA=0 ⇒ ax=0 gave DK<0 ⇒ DG floored to 0.1
                        # (seedling DBH frozen, e.g. BM 22960605010497 aspen 0.1→2.2 became 0.1→0.2). Match FVS +
                        # the cr/regent.f:158 & sprout.f:644 pattern.
                        bx = sd[:ht2][sp]; ax = c.ht_dbh_iabflg[sp] == 0 ? c.ht_dbh_aa[sp] : sd[:ht1][sp]
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
                end
            end
            if l == 0
                t.ht_growth[i] = htg
                small && (t.diam_growth[i] = dgk)
            elseif l == 1
                stash.htgU[i] = htg; !isempty(stash.is_small) && (stash.is_small[i] = true)
                small && (stash.dgU[i] = dgk)
            else
                stash.htgL[i] = htg
                small && (stash.dgL[i] = dgk)
            end
        end
    end
    return s
end

# bm/esgent.f (CALL REGENT(.TRUE.,ITRNIN)) — grow the JUST-ESTABLISHED regen IN its birth cycle. BM was OMITTED
# from the esgent dispatch (simulate.jl had CR/TT/EM/UT/CI), so planted/established BM seedlings never got their
# first-cycle height growth (BM BARE-PLANT: persistent TopHt lag ~5 ft). Same class as EM #137 / UT #184 / CI #185.
# Mirrors small_tree_growth!'s POTHTG/PCTRED/VIGOR/CON height + DK/DKK DBH over the birth subperiod (subyr=FINT−
# GENTIM=5), applying HT/DBH directly (esgent.f HT(I)=HT(I)+HTG(I)·WK4). Gated to the new records nstart+1:n.
function bm_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0,
                    atavh::Float32 = -1.0f0, atrelden::Float32 = -1.0f0)
    p, t, c = s.plot, s.trees, s.calib
    nstart >= t.n && return s
    sd = s.coef.species; slo = sd[:site_lo]; shi = sd[:site_hi]
    relden = p.relative_density; avh = p.avg_height; dgsd = s.control.dg_sd
    gentim = max(fint - 5.0f0, 0.0f0)
    bscale = (fint - gentim) / _BM_RG_REGYR              # birth-cycle fraction (WK4; =0.5 for fint=10)
    # REGENT(LESTB) PCTRED reads a MID-PERIOD blend of the CURRENT (post-growth) and the START-of-cycle
    # (post-thin, pre-growth) CCF/top-height (bm/regent.f:193-202): CCF=(5/FINT)*RELDEN+((FINT-5)/FINT)*ATCCF,
    # AVHT=(5/FINT)*AVH+((FINT-5)/FINT)*ATAVH, X=AVHT*(CCF/100). ATCCF/ATAVH = grincr.f:318-320 post-thin
    # (== simulate.jl es_at_relden/es_at_avh). Using the current-only stand inflated X (the grown top height)
    # and drove PCTRED to its floor => birth-cycle HTG collapse => one-directional PLANT under-production.
    # #194-class start-of-cycle fix, mirror of IE/UT. Missing (-1) sentinel => legacy current-only (defensive).
    w0 = fint > 0f0 ? 5.0f0 / fint : 1.0f0; w1 = 1.0f0 - w0
    ccf = atrelden >= 0f0 ? w0*relden + w1*atrelden : relden
    avht = atavh >= 0f0 ? w0*avh + w1*atavh : avh
    xd = avht * (ccf / 100.0f0); xd > 300.0f0 && (xd = 300.0f0)
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
            pothtg = bm_smhtgf(sp, sitear, h, _BM_RG_REGYR)  # SMHTGF reads raw SITEAR (unclamped)
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
        d >= bkpt && continue                         # bm/regent.f:390 D>=BKPT ⇒ GO TO 23 (large-tree)
        # REGENT(LESTB) birth diameter is the ABSOLUTE dubbed DK, NOT the DDS growth-increment reconstruction
        # (bm/regent.f:542-556 `IF(LESTB) … DBH(K)=DK ; IF(DBH<DIAM .OR. HK<4.5)DBH=DIAM ; DBH=DBH+0.001*HK ;
        # DG(K)=DBH(K)`). All BM planted species (13/14/16/18 have LHTDRG=false) take the plain DBH(K)=DK arm.
        # The former code ran the LESTB=F increment path (DG=(DK−DKK)*BARK, DDS-rescaled, capped at
        # DGMX=DGMAX·SCALE): the half-cycle DGMX cap clipped the synchronized PLANT cohort under DK. Mirror of
        # the IE/UT fix. HK≤4.5 ⇒ DBH=D+0.001*HK, DG=0 (bm/regent.f:392-394).
        if hk > 4.5f0
            local dk::Float32
            if sp == 6                                    # LP(7) fixed formula is dead code (HTDBH override); see small_tree_growth!
                dk = (hk - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
            elseif (!BM_LHTDRG[sp] || c.ht_dbh_iabflg[sp] == 1) && _bm_has_htdbh(Int(p.forest_idx), sp)
                dk = bm_htdbh(Int(p.forest_idx), sp, hk)
            else
                bx = sd[:ht2][sp]; ax = c.ht_dbh_iabflg[sp] == 0 ? c.ht_dbh_aa[sp] : sd[:ht1][sp]   # AX = IABFLG==1?HT1:AA (bm/regent.f:411)
                dk = bx / (log(hk - 4.5f0) - ax) - 1.0f0
            end
            dbh = dk
            dbh < BM_RG_DIAM[sp] && (dbh = BM_RG_DIAM[sp])   # regent.f:554 IF(DBH<DIAM)DBH=DIAM
            dbh = dbh + hk * 0.001f0                          # regent.f:555 DBH=DBH+0.001*HK
            t.dbh[i] = dbh; t.diam_growth[i] = dbh            # regent.f:556 DG(K)=DBH(K)
        else
            t.dbh[i] = d + 0.001f0 * hk; t.diam_growth[i] = 0.0f0   # regent.f:392-394 HK≤4.5
        end
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

# bm/htdbh.f MODE=0 (DBH→HT), Curtis-Arney with the linear D<3 segment. IFOR 1/2/3 → Malheur/Ochoco/
# Umatilla, anything else → Wallowa-Whitman (the Fortran ELSE branch).
@inline function bm_htdbh_height(ifor::Int, sp::Int, d::Float32)::Float32
    (ifor < 1 || ifor > 3) && (ifor = 4)
    p2 = BM_HTDBH_P2[ifor,sp]; p3 = BM_HTDBH_P3[ifor,sp]; p4 = BM_HTDBH_P4[ifor,sp]
    d >= 3f0 && return 4.5f0 + p2 * exp(-1f0 * p3 * d^p4)
    return ((4.5f0 + p2 * exp(-1f0 * p3 * (3f0^p4)) - 4.51f0) * (d - 0.3f0) / 2.7f0) + 4.51f0
end

# bm/blkdat.f:168-177 Wykoff HT-DBH HT1/HT2 (the cratet dub/AA-fit coefficients). The species CSV ht1/ht2
# columns are CR-template placeholders and differ (e.g. WJ 4.192/−5.1651 vs blkdat 3.2/−5.0).
const BM_BLK_HT1 = Float32[5.035, 5.043, 4.929, 4.874, 4.874, 3.2, 4.954, 5.035, 4.875, 4.993,
                           4.192, 4.192, 5.188, 5.143, 4.4421, 5.152, 4.993, 5.152]
const BM_BLK_HT2 = Float32[-10.674, -9.123, -10.744, -10.405, -10.405, -5.0, -9.177, -10.674, -9.568, -12.43,
                           -5.1651, -5.1651, -13.801, -13.497, -6.5405, -13.576, -12.43, -13.576]

# bm/cratet.f:385-412 (live DO 130) / :486-515 (dead DO 145) missing-height / top-kill dub for D>0.1:
# Wykoff H=exp(AX+HT2/(D+1))+4.5 (AX=AA if IABFLG==0 else HT1), WC small-tree forms for sp 13/14 & 16/18 at
# D<5, the PP/OS D<3 linear form — then, for every species except WJ/WB/LM/AS (6/11/12/15), the inventory
# HTDBH curve OVERRIDES H whenever .NOT.LHTDRG or IABFLG==1. The caller applies the 4.5 floor.
function bm_cratet_dub(ifor::Int, sp::Int, d::Float32, icr::Integer,
                       lhtdrg::Bool, iabflg::Integer, aa::Float32)::Float32
    ax = iabflg == 0 ? aa : BM_BLK_HT1[sp]
    h = if d < 5f0 && (sp == 13 || sp == 14)
        exp(1.5907f0 + 0.3040f0 * d)
    elseif d < 5f0 && (sp == 16 || sp == 18)
        0.0994f0 + 4.9767f0 * d
    else
        exp(ax + BM_BLK_HT2[sp] / (d + 1f0)) + 4.5f0
    end
    if (sp == 10 || sp == 17) && d < 3f0
        jcr = icr <= 0 ? 4 : clamp((icr - 1) ÷ 10 + 1, 1, 7)
        h = 8.31485f0 + 3.03659f0 * d - 0.59200f0 * jcr
    end
    (sp == 6 || sp == 11 || sp == 12 || sp == 15) && return h
    (!lhtdrg || iabflg == 1) && (h = bm_htdbh_height(ifor, sp, d))
    return h
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
