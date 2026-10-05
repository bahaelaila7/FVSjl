# =============================================================================
# small_tree_growth.jl (southcentraloregon) — SO small-tree growth (so/regent.f + so/smhtgf.f). Chunk 6.
#
# STEP 1 (this file): so_smhtgf — the per-species small-tree 5/10-yr HEIGHT increment (so/smhtgf.f SELECT
# CASE(ISPC)). Height-first small-tree model: SMHTGF gives HHT, regent then converts to DBH for trees
# crossing 4.5' and blends with the large-tree HTG via XWT. Most species are a linear-SI increment
# ((a+b·SI)/(c+d·SI))·DTIME·factor; WP(1)/RC(18) are the Chapman-Richards age-inversion form; WJ(11) a
# SI-ratio form; RF(9) the CA Ritchie&Hann DOMHTGR·SMHMOD; WO(27) the CA blackoak BAL form. MODE 0=from
# ESSUBH (estab), 1=from REGENT (invert age from H for the nonlinear species). MEASURED vs FVSso_g16 smhtgf.
#
# STEP 2 (TODO — the regent DBH driver): small_tree_growth!(::SouthCentralOregon) mirroring the CA driver,
# using so/regent.f's OWN height-DBH regression (HTT1/HTT2 per-crown-group, HT1/HT2, AA, IABFLG from
# so/blkdat.f — NOT the chunk-4a Curtis so_htdbh) + XMAX/XMIN/DIAM/REGYR. Needed to unblock grow_cycle!.
# =============================================================================

# so/smhtgf.f species case sets.
const SO_SMHT_WFGRP = (4, 12)      # WF/GF (·1.2)
const SO_SMHT_PYGRP = (20, 21, 23, 25, 26, 28, 29, 30, 31, 33)  # WC PY-form

# so/smhtgf.f — 5/10-yr small-tree height increment (ft). `cr1` = ICR (integer percent, 0..100).
function so_smhtgf(sp::Int, d::Float32, h::Float32, cr1::Float32, dtime::Float32, mode::Int;
                   si::Float32, ba::Float32, pct::Float32, avh::Float32)::Float32
    S = si
    if sp == 1 || sp == 18                                  # WP / RC — Chapman-Richards age inversion
        c1, c2, c3, c4 = sp == 1 ? (0.375045f0, 0.92503f0, -0.020796f0, 2.48811f0) :
                                   (0.752842f0, 1.0f0, -0.0174f0, 1.4711f0)
        effage = mode == 1 ? flog((1f0 - fpow(c1 / S * h, 1f0 / c4)) / c2) / c3 : 0f0
        agepdt = effage + dtime
        hht1 = (S / c1) * fpow(1f0 - c2 * fexp(c3 * effage), c4)
        hht2 = (S / c1) * fpow(1f0 - c2 * fexp(c3 * agepdt), c4)
        return hht2 - hht1
    elseif sp == 2 || sp == 10                              # SP / PP
        return ((-1.0f0 + 0.32857f0 * S) / (28.0f0 - 0.042857f0 * S)) * dtime
    elseif sp == 3 || sp == 32                              # DF / OS
        return ((2.0f0 + 0.420f0 * S) / (28.5f0 - 0.05f0 * S)) * dtime * 1.1f0
    elseif sp in SO_SMHT_WFGRP                              # WF / GF
        return ((4.2435f0 + 0.1510f0 * S) / (19.0184f0 - 0.0570f0 * S)) * dtime * 1.2f0
    elseif sp == 5                                          # MH (metric → ft)
        return ((0.965758f0 + 0.082969f0 * S) / (55.249612f0 - 1.288852f0 * S)) * dtime * 3.280833f0 * 1.60f0
    elseif sp == 6                                          # IC
        return ((4.2435f0 + 0.1510f0 * S) / (19.0184f0 - 0.0570f0 * S)) * dtime * 1.3f0
    elseif sp == 7 || sp == 16                              # LP / WB
        h0 = 0.02008805f0 * S * dtime
        return sp == 16 ? 1.6f0 * h0 : h0
    elseif sp == 8                                          # ES
        return ((0.09211f0 + 0.208517f0 * S) / (43.358f0 - 0.168166f0 * S)) * dtime * 1.35f0
    elseif sp == 9                                          # RF — CA Ritchie&Hann
        relht = avh <= 0f0 ? 1f0 : h / avh
        relht > 1.05f0 && (relht = 1.05f0)
        domhtgr = 5f0 * (2.2227f0 + 0.4314f0 * S) / (29.0f0 - 0.05f0 * S)
        cr = cr1 / 100f0
        crmod = 1.0f0 - fexp(-4.26558f0 * cr)
        rhmod = fexp(2.54119f0 * (fpow(relht, 0.250537f0) - 1.0f0))
        smhmod = 1.016605f0 * crmod * rhmod                 # smhtgf.f: SMHMOD formed first, HHT=DOMHTGR*SMHMOD
        return domhtgr * smhmod
    elseif sp == 11                                        # WJ — SI-ratio (SI clamp [SLO+.5, SHI])
        sj = S; sj > 75f0 && (sj = 75f0); sj <= 5f0 && (sj = 5.5f0)
        return (sj / 5.0f0) * (sj * 1.5f0 - h) / (sj * 1.5f0)
    elseif sp == 13                                        # AF
        return ((6.0f0 + 0.14f0 * S) / (33.882f0 - 0.06588f0 * S)) * dtime
    elseif sp == 14                                        # SF (EC form)
        return ((-0.6667f0 + 0.4333f0 * S) / (28.5f0 - 0.05f0 * S)) * dtime
    elseif sp == 15                                        # NF (WC form)
        return ((11.26677f0 + 0.12027f0 * S) / (27.93806f0 - 0.02873f0 * S)) * dtime
    elseif sp == 17                                        # WL (EC form)
        return ((-3.9725f0 + 0.50995f0 * S) / (28.1168f0 - 0.05661f0 * S)) * dtime
    elseif sp == 19                                        # WH (WC form)
        return ((-5.74874f0 + 0.54576f0 * S) / (26.15767f0 - 0.03596f0 * S)) * dtime
    elseif sp in SO_SMHT_PYGRP                             # PY + WC catch-all
        return ((1.47043f0 + 0.23317f0 * S) / (31.56252f0 - 0.05586f0 * S)) * dtime
    elseif sp == 22                                        # RA (WC form)
        return (-0.007025f0 + 0.056794f0 * S) * dtime * 1.20f0
    elseif sp == 27                                        # WO — CA blackoak (BAL form)
        bal = ((100.0f0 - pct) / 100.0f0) * ba
        factor = 0.80f0 + 0.004f0 * (S - 50.0f0)
        bal < 5.0f0 && (bal = 5.0f0)
        return fexp(3.817f0 - 0.7829f0 * flog(bal)) * factor
    end
    return 0f0
end

# =============================================================================
# STEP 2 — small_tree_growth!(::SouthCentralOregon): the so/regent.f DBH driver.
#
# Height-first: so_smhtgf (STEP 1) gives a 10-yr HHT increment; regent modifies it by a stand-density
# modifier (PCTRED, from AB poly on AVHT·CCF) and a per-tree crown VIGOR, adds ZZRAN·0.1 (DGSD≥1 fires
# for SO), scales to FINT (SCALE=FNT/REGYR), and blends with the large-tree HTG by XWT=(D−XMIN)/(XMAX−XMIN).
# DBH is derived from the HEIGHT crossing 4.5' for D<3" (D<99 for WJ11).
#
# ★ KEY (measured in so/regent.f, NOT inferred): SO sets LHTDRG(ISPC)=.FALSE. for ALL species (grinit.f:102).
#   regent.f:557-569 then OVERRIDES the DK/DKK computed by the HTT1/HTT2/HT1/HT2/RDCON regression forms with
#   `CALL HTDBH(IFOR,ISPC,DK,HK,1)` — i.e. the chunk-4a Curtis so_htdbh_dbh — for EVERY species except WJ(11)/
#   WB(16)/AS(24) (which GOTO 300 and keep their special forms). So the ~825-value HTT/HT/RDCON transcription
#   is dead code under SO's default config; the DK conversion collapses to so_htdbh_dbh + three specials.
#   (WJ = SITEAR linear; WB = EM-SMDGF form; AS = HT1(24)/HT2(24) ln-form.)
# =============================================================================

const SO_RG_XMAX  = Float32[3,5,4,4,2,4,4,4,4,5, 99,4,4,4,4,3,4,10,4,4, 4,4,4,4,4,4,4,4,4,4, 4,4,4]
const SO_RG_XMIN  = Float32[2,1,2,2,1,2,2,2,2,1, 90,2,2,2,2,1.5,2,2,2,2, 2,2,2,2,2,2,2,2,2,2, 2,2,2]
const SO_RG_DIAM  = Float32[0.4,0.4,0.3,0.3,0.2,0.2,0.4,0.3,0.2,0.5, 0.3,0.3,0.3,0.3,0.3,0.4,0.3,0.2,0.2,0.2,
                            0.2,0.3,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2, 0.2,0.3,0.2]
const SO_RG_DGMAX = Float32[2.8,2.8,2.4,3.6,2.5,2.5,3.5,3.6,3.6,2.8, 2.0,3.6,3.6,3.6,5.0,2.8,2.8,2.5,5.0,5.0,
                            5.0,5.0,5.0,2.5,5.0,5.0,2.8,5.0,5.0,5.0, 5.0,2.4,5.0]
const SO_RG_AB    = Float32[1.11436, -0.011493, 0.43012f-4, -0.72221f-7, 0.5607f-10, -0.1641f-13]  # density-modifier poly
const SO_RG_HT1_24 = 4.4421f0     # blkdat.f HT1(24)=AS  (AX for the ln-form)
const SO_RG_HT2_24 = -6.5405f0    # blkdat.f HT2(24)=AS  (BX for the ln-form)
# so/regent.f:520 CASE(15,19:23,25,26,28:31,33) — WC-variant hardwoods (double-DDS quirk, DGBND applied).
const SO_RG_HARDWOOD = Set([15, 19, 20, 21, 22, 23, 25, 26, 28, 29, 30, 31, 33])

# so/regent.f label 5 → 23 for one slot K: returns (DBH(K) set directly (<0 ⇒ none), DG(K)). D≥BKPT ⇒ GO TO 23
# (DG(K) = the slot's DGDRIV DG `dglt`, untouched); HK≤4.5 ⇒ DG=0, DBH(K)=D+0.001·HK (then DGBND); SH/WO (9,27)
# jump to 23 after their XDWT blend (no DGBND). DBH(K)=D for every slot (dgdriv.f tripling copies DBH(I)).
function _so_regent_dg(s::StandState, sp::Int, ifor::Int, d::Float32, h::Float32, htg::Float32, icr::Float32,
                       si_raw::Float32, bkpt::Float32, scale::Float32, scale2::Float32, dgmx::Float32,
                       dglt::Float32, xrdgro::Float32, i::Int)
    t = s.trees; sd = s.coef.species
    d >= bkpt && return (-1f0, dglt)
    hk = h + htg
    if hk <= 4.5f0
        return (d + 0.001f0*hk, dg_bound(nothing, nothing, sp, d + 0.001f0*hk, 0f0, s.control.sp_size_cap))
    end
    bark = so_bratio(sd, sp, d)
    # --- DK/DKK (regent.f:424-569): htdbh override for all but WJ/WB/AS ---
    local dk::Float32, dkk::Float32
    if sp == 11                                        # WJ — SITEAR linear
        dk = (hk - 4.5f0)*10f0/(si_raw - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
        dkk = (h - 4.5f0)*10f0/(si_raw - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0)
        h < 4.5f0 && (dkk = d)
    elseif sp == 16                                    # WB — EM SMDGF form (PPCCF=1)
        pt = Int(t.plot_id[i])
        tpccf = (pt >= 1 && pt <= length(s.density.point_ccf)) ? s.density.point_ccf[pt] : 0f0
        tpccf > 300f0 && (tpccf = 300f0); tpccf < 25f0 && (tpccf = 25f0)
        hl = h - 4.5f0
        dkk = 0.000231f0*hl*icr - 0.00005f0*hl*tpccf + 0.001711f0*icr + 0.17023f0*hl + 0.3f0
        hlk = hk - 4.5f0
        dk = 0.000231f0*hlk*icr - 0.00005f0*hlk*tpccf + 0.001711f0*icr + 0.17023f0*hlk + 0.3f0
    elseif sp == 24                                    # AS — HT1(24)/HT2(24) ln-form (GOTO 300, no override)
        dk = (SO_RG_HT2_24/(flog(hk - 4.5f0) - SO_RG_HT1_24)) - 1f0
        dkk = h <= 4.5f0 ? d : (SO_RG_HT2_24/(flog(h - 4.5f0) - SO_RG_HT1_24)) - 1f0
    else                                               # all others: htdbh override (LHTDRG=false)
        dk = so_htdbh_dbh(ifor, sp, hk)
        dkk = h <= 4.5f0 ? d : so_htdbh_dbh(ifor, sp, h)
    end
    h < 4.5f0 && (dkk = d)                             # regent.f:605
    # --- DG (LESTB=false, regent.f:597-692) ---
    local dg::Float32
    if sp in SO_RG_HARDWOOD                            # WC hardwoods — double-DDS quirk, then DGBND
        dg = (dk < 0f0 || dkk < 0f0) ? htg*0.2f0*bark*xrdgro : (dk - dkk)*bark*xrdgro
        dg < 0f0 && (dg = 0.1f0); dg > dgmx && (dg = dgmx)
        dds = dg*(2f0*bark*d + dg)*scale2
        dg = sqrt(fpow(d*bark, 2f0) + dds) - bark*d
        dg < 0f0 && (dg = 0f0); dg > dgmx && (dg = dgmx)  # common block (681-688): re-converts
        dds = dg*(2f0*bark*d + dg)*scale2
        dg = sqrt(fpow(d*bark, 2f0) + dds) - bark*d
    elseif sp == 9 || sp == 27                         # SH/WO — XDWT blend; GO TO 23 skips DGBND
        xdwt = d <= 1.5f0 ? 0f0 : d >= 3f0 ? 1f0 : (d - 1.5f0)/1.5f0
        dgsm = (dk - dkk)*bark*xrdgro; dgsm < 0f0 && (dgsm = 0f0)
        dds = dgsm*(2f0*bark*d + dgsm)*scale2
        dgsm = sqrt(fpow(d*bark, 2f0) + dds) - bark*d; dgsm < 0f0 && (dgsm = 0f0)
        dg = dgsm*(1f0 - xdwt) + dglt*xdwt
        (d + dg) < SO_RG_DIAM[sp] && (dg = SO_RG_DIAM[sp] - d)
        return (-1f0, dg)
    else                                               # DEFAULT
        dg = (dk < 0f0 || dkk < 0f0) ? htg*0.2f0*bark*xrdgro : (dk - dkk)*bark*xrdgro
        dg < 0f0 && (dg = 0f0); dg > dgmx && (dg = dgmx)
        (sp == 16 && (d + dg) < SO_RG_DIAM[sp]) && (dg = SO_RG_DIAM[sp] - d)
        dds = dg*(2f0*bark*d + dg)*scale2
        dg = sqrt(fpow(d*bark, 2f0) + dds) - bark*d
    end
    (d + dg) < SO_RG_DIAM[sp] && (dg = SO_RG_DIAM[sp] - d)
    return (-1f0, dg_bound(nothing, nothing, sp, d, dg, s.control.sp_size_cap))   # so/dgbnd.f = SIZCAP cap
end

function small_tree_growth!(s::StandState, stash, ::SouthCentralOregon; fint::Float32 = 10.0f0)
    p, t, c = s.plot, s.trees, s.calib
    n = t.n; n == 0 && return s
    sd = s.coef.species
    # AVH = the COMMON AVHT40 of the last DENSE (cycle start: CRATET's IND at cycle 0, gradd.f:186's after)
    avh = p.avg_height; ba = p.basal_area; dgsd = s.control.dg_sd
    relden = p.relative_density; ifor = Int(p.forest_idx); yr = htg_period(s.variant)
    fnt = fint                                             # LESTB=false (cycling): FNT=FINT
    # PCTRED density modifier (regent.f:220-225) — computed once from stand AVHT·CCF.
    xden = avh * (relden / 100f0); xden > 300f0 && (xden = 300f0)
    pctred = SO_RG_AB[1] + xden*(SO_RG_AB[2] + xden*(SO_RG_AB[3] + xden*(SO_RG_AB[4] +
             xden*(SO_RG_AB[5] + xden*SO_RG_AB[6]))))
    pctred > 1f0 && (pctred = 1f0); pctred < 0.01f0 && (pctred = 0.01f0)
    trip = stash !== nothing && !isempty(stash.dgU)
    yr_now = current_cycle_year(s)                         # so/regent.f MULTS(3/6, IY(ICYC))
    # so/regent.f:232-269 is SPECIES-MAJOR (DO 30 ISPC … I=IND1(I3)); the per-tree ZZRAN draw must follow it.
    @inbounds for i in species_major_order(s)
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= SO_RG_XMAX[sp] || t.tpa[i] <= 0f0) && continue
        h = t.height[i]
        icr = Float32(t.crown_pct[i])                      # ICR (integer crown percent 0..100)
        slo = SO_SITELO[sp]; shi = SO_SITEHI[sp]
        si_raw = p.sp_site_index[sp]                        # SITEAR (raw; so_smhtgf/WJ use this)
        si_c = si_raw; si_c > shi && (si_c = shi); si_c <= slo && (si_c = slo + 0.5f0)  # regent SI clamp (AS RELSI)
        con = fexp(c.htg_cor_small[sp])                     # RHCON(=1)·fexp(HCOR); HCOR=0 (no small-tree calib)
        xrhgro = active_multiplier(s.control, :regh, sp, yr_now)   # XRHGRO=XRHMLT(ISPC) (REGHMULT)
        xrdgro = active_multiplier(s.control, :regd, sp, yr_now)   # XRDGRO=XRDMLT(ISPC) (REGDMULT)
        regyr = (sp == 9 || sp == 27) ? 5f0 : 10f0
        scale = fnt / regyr; scale2 = yr / fnt
        dgmx = sp == 16 ? fint*0.2f0 : SO_RG_DGMAX[sp]*scale
        # --- VIGOR modifier from crown ratio (regent.f:303-310) ---
        xv = icr / 100f0
        vigor = 150f0 * fpow(xv, 3f0) * fexp(-6f0*xv) + 0.3f0
        vigor > 1f0 && (vigor = 1f0)
        sp == 11 && (vigor = 1f0 - (1f0 - vigor)/3f0)      # WJ: cut knock-down by 2/3
        # --- POTHTG → HTGR ---
        local htgr::Float32
        if sp == 24                                        # AS: Sheppard aspen (findag effective age)
            sitage, _, _, _ = so_findag(24, ifor, h, si_raw)
            relsi = (si_c - slo) / (shi - slo); rsimod = 0.5f0 * (1f0 + relsi)
            hite1 = 26.9825f0 * fpow(sitage, 1.1752f0)
            hite2 = 26.9825f0 * fpow(sitage + 10f0, 1.1752f0)
            htgr = (hite2 - hite1) / (2.54f0*12f0) * rsimod * con
            htgr *= 2.40f0; htgr *= 0.75f0                 # Dixon 8-27-92 −25%
        else
            pothtg = so_smhtgf(sp, d, h, icr, 10f0, 1; si = si_raw, ba = ba,
                               pct = t.crown_ratio[i], avh = avh)
            htgr = (sp == 9 || sp == 27) ? pothtg*con : pothtg*pctred*vigor*con  # SH/WO skip PCTRED·VIGOR
        end
        xmn = SO_RG_XMIN[sp]; xmx = SO_RG_XMAX[sp]
        xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)
        cap = s.control.sp_size_cap[sp, 4]
        bkpt = sp == 11 ? 99f0 : 3f0                       # BKPT: 3", 99 for WJ(11)
        large_htg = t.ht_growth[i]; dg_main = t.diam_growth[i]
        htgr0 = htgr
        # so/regent.f:708-711 — with LTRIP each tripled copy K=ITRN+2I−2+L (L=1,2) reruns label 2 → 23: a FRESH
        # ZZRAN, the blend with the copy's own HTG(K), SIZCAP and (D<BKPT) its own DG(K)/DBH(K). H and D stay
        # record I's; DBH(K)=DBH(I) (dgdriv.f). POTHTG (label 2) is deterministic in H/D/ICR ⇒ computed once.
        for l in 0:(trip ? 2 : 0)
            lthg = l == 0 ? large_htg : (stash.htg_copy[i] ? (l == 1 ? stash.htgU[i] : stash.htgL[i]) : large_htg)
            dglt = l == 0 ? dg_main : (l == 1 ? stash.dgU[i] : stash.dgL[i])
            # --- random ht component (label 4; SO DGSD=2.0 ≥ 1 ⇒ fires every tree) ---
            zzran = 0f0
            if dgsd >= 1f0
                while true
                    zzran = bachlo(s.rng, 0f0, 1f0)
                    (zzran <= 0.5f0 && zzran >= -2.0f0) && break
                end
            end
            htgr = (htgr0 + zzran*0.1f0) * xrhgro * scale
            # --- blend small & large tree HTG ---
            htg = htgr*(1f0 - xwt) + xwt*lthg
            htg < 0.1f0 && (htg = 0.1f0)
            (h + htg > cap) && (htg = max(cap - h, 0.1f0))
            # WB DK/DKK read ICR(K) (so/regent.f:494-512). A tripled copy K = ITRN+2I−2+L is the copy's FUTURE slot, which
            # TRIPLE fills only later ⇒ FVS reads that slot's current contents: 0 if never used, else what TREDEL left
            # (t.stale_icr) — the same quirk as bm/regent.f (MEASURED FVSso_g16 15184869010497 2020 WB K=416: ICR 0).
            icrk = l == 0 ? icr : Float32(t.stale_icr[n + 2i - 2 + l])
            dbhk, dg = _so_regent_dg(s, sp, ifor, d, h, htg, icrk, si_raw, bkpt, scale, scale2, dgmx, dglt, xrdgro, i)
            if l == 0
                t.ht_growth[i] = htg; t.diam_growth[i] = dg
                dbhk >= 0f0 && (t.dbh[i] = dbhk)
            else
                l == 1 ? (stash.htgU[i] = htg) : (stash.htgL[i] = htg)
                l == 1 ? (stash.dgU[i] = dg) : (stash.dgL[i] = dg)
                l == 1 ? (stash.dbhU[i] = dbhk >= 0f0 ? dbhk : d) : (stash.dbhL[i] = dbhk >= 0f0 ? dbhk : d)
                stash.is_small[i] = true
            end
        end
    end
    return s
end

"""so/regent.f LSTART small-tree HEIGHT calibration (label 40; called from cratet.f:681). PCTRED from
X=AVH·RELDEN/100 (AVH = cratet.f:639 AVHT40, RELDEN = the cratet.f:171 backdating DENSE's). Per species (IND1
order): current DBH<5, backdated H=HT−HTG (IHTG<2) ≥0.01, measured HTG≥0.001 ⇒ EDH = POTHTG·PCTRED·VIGOR·RHCON
(SH/WO 9,27: POTHTG·RHCON; AS 24: the Sheppard curve from H on the clamped-SI RELSI ×2.4×0.75), TERM = HTG·SCALE3
(SCALE3 = REGYR/FINTH, REGYR 5 for 9/27 else 10); CORNEW = mean TERM / mean EDH when N≥NCALHT(5), ≤0 ⇒ 1E−4,
trapped to [0.0821, 12.1825] else 1. HCOR = ln(CORNEW) → htg_cor_init (dgdriv.f:178/199 attenuation applies it)."""
function so_regent_hcor_init!(s::StandState, isct::AbstractMatrix, ind1::AbstractVector,
                              saved_dbh::AbstractVector, avh::Float32)
    p, t, c = s.plot, s.trees, s.calib
    t.n == 0 && return s
    s.control.growth_ifinth == 0 && return s
    ba = c.cratet_ba; relden = c.cratet_relden
    finth = s.control.growth_finth > 0f0 ? s.control.growth_finth : 5f0
    x = avh * (relden / 100f0); x > 300f0 && (x = 300f0)
    pctred = SO_RG_AB[1] + x*(SO_RG_AB[2] + x*(SO_RG_AB[3] + x*(SO_RG_AB[4] + x*(SO_RG_AB[5] + x*SO_RG_AB[6]))))
    pctred > 1f0 && (pctred = 1f0); pctred < 0.01f0 && (pctred = 0.01f0)
    @inbounds for sp in 1:size(isct, 1)
        i1 = Int(isct[sp, 1]); i1 == 0 && continue
        i2 = Int(isct[sp, 2])
        slo = SO_SITELO[sp]; shi = SO_SITEHI[sp]
        si_raw = p.sp_site_index[sp]
        si = si_raw; si > shi && (si = shi); si <= slo && (si = slo + 0.5f0)
        regyr = (sp == 9 || sp == 27) ? 5f0 : 10f0
        scale3 = regyr / finth
        rhcon = (s.control.regh_cor2_on && s.control.regh_cor2[sp] > 0f0) ? s.control.regh_cor2[sp] : 1f0
        snp = 0f0; snx = 0f0; sny = 0f0; nh = 0
        for k in i1:i2
            i = Int(ind1[k])
            h = t.height[i]
            s.control.growth_ihtg < 2 && (h = h - t.ht_growth[i])
            (saved_dbh[i] >= 5f0 || h < 0.01f0) && continue
            icr = Float32(t.crown_pct[i])
            xv = icr / 100f0
            vigor = 150f0 * fpow(xv, 3f0) * fexp(-6f0*xv) + 0.3f0
            vigor > 1f0 && (vigor = 1f0)
            sp == 11 && (vigor = 1f0 - (1f0 - vigor)/3f0)
            local edh::Float32
            if sp == 24                                       # AS — Sheppard (regent.f CASE(24))
                relsi = (si - slo) / (shi - slo); rsimod = 0.5f0 * (1f0 + relsi)
                ag1 = fpow(h*12f0*2.54f0/26.9825f0, 0.8509f0)
                h2 = (26.9825f0 * fpow(ag1 + 10f0, 1.1752f0)) / (2.54f0*12f0)
                edh = (h2 - h) * rsimod * rhcon
                edh *= 2.4f0; edh *= 0.75f0
                edh < 0f0 && (edh = 0f0)
            else
                pct = i <= length(c.cratet_pct) ? c.cratet_pct[i] : t.crown_ratio[i]
                pothtg = so_smhtgf(sp, saved_dbh[i], h, icr, 10f0, 1; si = si_raw, ba = ba, pct = pct, avh = avh)
                edh = (sp == 9 || sp == 27) ? pothtg * rhcon : pothtg * pctred * vigor * rhcon
            end
            hg = t.ht_growth[i]; hg < 0.001f0 && continue
            pr = t.tpa[i]
            snp += pr; snx += edh * pr; sny += hg * scale3 * pr; nh += 1
        end
        nh < 5 && continue
        snx /= snp; sny /= snp
        cornew = sny / snx
        cornew <= 0f0 && (cornew = 1f-4)
        (cornew < 0.0821f0 || cornew > 12.1825f0) && (cornew = 1f0)
        c.htg_cor_init[sp] = flog(cornew)
    end
    return s
end

"""
    so_esgent!(s, nstart; fint, atavh, atrelden, relden_pre, avh_pre, ba_pre, pccf_pre)

strp/esgent.f → so/regent.f REGENT(.TRUE.,ITRNIN) for the records ESTAB created this cycle (nstart+1:n), then
esgent.f's WK4 step (`esgent_finish!`). SO was missing from the esgent dispatch, so planted/established SO seedlings
kept their ESSUBH base height and 0.1" DBH through their birth cycle.

REGENT(LESTB): FNT=FINT−5 (LSKIPH when FINT≤5, regent.f:199-206); PCTRED from the mid-period blend
CCF=(5/FINT)·RELDEN+((FINT−5)/FINT)·ATCCF, AVHT likewise with ATAVH (:211-225). Species-major over the new records
(DO 30 ISPC / DO 25 I3, I<ITRNIN skipped): the open-grown crown draw (:282-289, PCCF of the gradd.f:192 DENSE),
then (unless LSKIPH) VIGOR, POTHTG=SMHTGF(MODE 1, H, ICR, 10) (AS: the Sheppard curve on SITAGE=ABIRTH), HTGR=
POTHTG·PCTRED·VIGOR·CON (SH/WO: POTHTG·CON), the ZZRAN draw, HTGR=(HTGR+0.1·ZZRAN)·XRHGRO·SCALE, XWT=0, HTG≥0.1 and
the SIZCAP cap. DBH (:416-597): HK≤4.5 ⇒ DBH=D+0.001·HK, DG=0; else DBH=DK (floored at DIAM), +0.001·HK, DG=DBH;
then DGBND. RELDEN/BA/AVH/PCCF are the post-growth PRE-regen values (ESTAB runs before gradd.f:244's DENSE); the
new record's PCT is 0 (estab.f), so BAL=BA inside SMHTGF.
"""
function so_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0, atavh::Float32 = 0f0, atrelden::Float32 = 0f0,
                    relden_pre::Float32 = 0f0, avh_pre::Float32 = 0f0, ba_pre::Float32 = 0f0,
                    pccf_pre::Vector{Float32} = Float32[])
    p, t, c = s.plot, s.trees, s.calib
    nstart >= t.n && return s
    sd = s.coef.species
    lskiph = fint <= 5f0
    fnt = lskiph ? fint : fint - 5f0
    ccf = relden_pre; avht = avh_pre
    if !lskiph && fnt > 0f0
        ccf = (5f0 / fint) * relden_pre + ((fint - 5f0) / fint) * atrelden
        avht = (5f0 / fint) * avh_pre + ((fint - 5f0) / fint) * atavh
    end
    xden = avht * (ccf / 100f0); xden > 300f0 && (xden = 300f0)
    pctred = SO_RG_AB[1] + xden*(SO_RG_AB[2] + xden*(SO_RG_AB[3] + xden*(SO_RG_AB[4] +
             xden*(SO_RG_AB[5] + xden*SO_RG_AB[6]))))
    pctred > 1f0 && (pctred = 1f0); pctred < 0.01f0 && (pctred = 0.01f0)
    ifor = Int(p.forest_idx); dgsd = s.control.dg_sd
    yr_now = current_cycle_year(s)
    newidx = sort(collect((nstart + 1):t.n); by = i -> (Int(t.species[i]), i))   # esgent.f:49 SPESRT → IND1
    @inbounds for i in newidx
        sp = Int(t.species[i]); d = t.dbh[i]
        d >= SO_RG_XMAX[sp] && continue
        ip = Int(t.plot_id[i])
        pccf = (1 <= ip <= length(pccf_pre)) ? pccf_pre[ip] : s.density.point_ccf[ip]
        ran = 0f0
        while true
            ran = bachlo(s.rng, 0f0, 1f0); (-1f0 <= ran <= 1f0) && break
        end
        cr0 = 0.89722f0 - 0.0000461f0 * pccf
        cr0 = cr0 + 0.07985f0 * ran
        cr0 > 0.90f0 && (cr0 = 0.90f0); cr0 < 0.20f0 && (cr0 = 0.20f0)
        icr0 = unsafe_trunc(Int32, cr0 * 100f0 + 0.5f0)
        t.crown_pct[i] = icr0; t.crown_ratio[i] = Float32(icr0)
        icr = Float32(icr0)
        h = t.height[i]
        slo = SO_SITELO[sp]; shi = SO_SITEHI[sp]
        si_raw = p.sp_site_index[sp]
        si_c = si_raw; si_c > shi && (si_c = shi); si_c <= slo && (si_c = slo + 0.5f0)
        regyr = (sp == 9 || sp == 27) ? 5f0 : 10f0
        scale = fnt / regyr
        local htg::Float32
        if lskiph
            htg = 0f0
        else
            con = fexp(c.htg_cor_small[sp])
            xrhgro = active_multiplier(s.control, :regh, sp, yr_now)
            xv = icr / 100f0
            vigor = 150f0 * fpow(xv, 3f0) * fexp(-6f0*xv) + 0.3f0
            vigor > 1f0 && (vigor = 1f0)
            sp == 11 && (vigor = 1f0 - (1f0 - vigor)/3f0)
            local htgr::Float32
            if sp == 24                                           # AS: SITAGE = ABIRTH under LESTB (:320-321)
                sitage = t.birth_age[i]
                relsi = (si_c - slo) / (shi - slo); rsimod = 0.5f0 * (1f0 + relsi)
                hite1 = 26.9825f0 * fpow(sitage, 1.1752f0)
                hite2 = 26.9825f0 * fpow(sitage + 10f0, 1.1752f0)
                htgr = (hite2 - hite1) / (2.54f0*12f0) * rsimod * con
                htgr *= 2.40f0; htgr *= 0.75f0
            else
                pothtg = so_smhtgf(sp, d, h, icr, 10f0, 1; si = si_raw, ba = ba_pre, pct = 0f0, avh = avh_pre)
                htgr = (sp == 9 || sp == 27) ? pothtg*con : pothtg*pctred*vigor*con
            end
            zzran = 0f0
            if dgsd >= 1f0
                while true
                    zzran = bachlo(s.rng, 0f0, 1f0)
                    (zzran <= 0.5f0 && zzran >= -2.0f0) && break
                end
            end
            htgr = (htgr + zzran*0.1f0) * xrhgro * scale
            htg = htgr; htg < 0.1f0 && (htg = 0.1f0)             # XWT=0 under LESTB (:394)
            cap = s.control.sp_size_cap[sp, 4]
            (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        end
        # label 5: DBH (BKPT 3", WJ 99")
        bkpt = sp == 11 ? 99f0 : 3f0
        if d < bkpt
            hk = h + htg
            local dbhk::Float32, dgk::Float32
            if hk <= 4.5f0
                dgk = 0f0; dbhk = d + 0.001f0*hk
            else
                local dk::Float32
                if sp == 11
                    dk = (hk - 4.5f0)*10f0/(si_raw - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
                elseif sp == 16
                    tpccf = pccf; tpccf > 300f0 && (tpccf = 300f0); tpccf < 25f0 && (tpccf = 25f0)
                    hlk = hk - 4.5f0
                    dk = 0.000231f0*hlk*icr - 0.00005f0*hlk*tpccf + 0.001711f0*icr + 0.17023f0*hlk + 0.3f0
                elseif sp == 24
                    dk = (SO_RG_HT2_24/(flog(hk - 4.5f0) - SO_RG_HT1_24)) - 1f0
                else
                    dk = so_htdbh_dbh(ifor, sp, hk)
                end
                dbhk = dk                                         # LESTB: DBH(K)=DK (:574-596)
                (dbhk < SO_RG_DIAM[sp] || hk < 4.5f0) && (dbhk = SO_RG_DIAM[sp])
                dbhk = dbhk + 0.001f0*hk
                dgk = dbhk
                (dbhk + dgk) < SO_RG_DIAM[sp] && (dgk = SO_RG_DIAM[sp] - dbhk)
            end
            dgk = dg_bound(nothing, nothing, sp, dbhk, dgk, s.control.sp_size_cap)   # DGBND(ISPC,DBH(K),DG(K))
            t.dbh[i] = dbhk; t.diam_growth[i] = dgk
        end
        esgent_finish!(t, i, htg, _SO_ES_HHTMAX[sp])
    end
    return s
end
