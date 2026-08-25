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
        effage = mode == 1 ? log((1f0 - (c1 / S * h)^(1f0 / c4)) / c2) / c3 : 0f0
        agepdt = effage + dtime
        hht1 = (S / c1) * (1f0 - c2 * exp(c3 * effage))^c4
        hht2 = (S / c1) * (1f0 - c2 * exp(c3 * agepdt))^c4
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
        crmod = 1.0f0 - exp(-4.26558f0 * cr)
        rhmod = exp(2.54119f0 * (relht^0.250537f0 - 1.0f0))
        return domhtgr * 1.016605f0 * crmod * rhmod
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
        return exp(3.817f0 - 0.7829f0 * log(bal)) * factor
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

function small_tree_growth!(s::StandState, stash, ::SouthCentralOregon; fint::Float32 = 10.0f0)
    p, t, c = s.plot, s.trees, s.calib
    n = t.n; n == 0 && return s
    sd = s.coef.species
    avh = stand_top_height(s); ba = p.basal_area; dgsd = s.control.dg_sd
    relden = p.relative_density; ifor = Int(p.forest_idx); yr = s.control.year
    fnt = fint                                             # LESTB=false (cycling): FNT=FINT
    # PCTRED density modifier (regent.f:220-225) — computed once from stand AVHT·CCF.
    xden = avh * (relden / 100f0); xden > 300f0 && (xden = 300f0)
    pctred = SO_RG_AB[1] + xden*(SO_RG_AB[2] + xden*(SO_RG_AB[3] + xden*(SO_RG_AB[4] +
             xden*(SO_RG_AB[5] + xden*SO_RG_AB[6]))))
    pctred > 1f0 && (pctred = 1f0); pctred < 0.01f0 && (pctred = 0.01f0)
    @inbounds for i in 1:n
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= SO_RG_XMAX[sp] || t.tpa[i] <= 0f0) && continue
        h = t.height[i]
        icr = Float32(t.crown_pct[i])                      # ICR (integer crown percent 0..100)
        slo = SO_SITELO[sp]; shi = SO_SITEHI[sp]
        si_raw = p.sp_site_index[sp]                        # SITEAR (raw; so_smhtgf/WJ use this)
        si_c = si_raw; si_c > shi && (si_c = shi); si_c <= slo && (si_c = slo + 0.5f0)  # regent SI clamp (AS RELSI)
        con = exp(c.htg_cor_small[sp])                     # RHCON(=1)·exp(HCOR); HCOR=0 (no small-tree calib)
        regyr = (sp == 9 || sp == 27) ? 5f0 : 10f0
        scale = fnt / regyr; scale2 = yr / fnt
        dgmx = sp == 16 ? fint*0.2f0 : SO_RG_DGMAX[sp]*scale
        # --- VIGOR modifier from crown ratio (regent.f:303-310) ---
        xv = icr / 100f0
        vigor = 150f0 * fpow(xv, 3f0) * exp(-6f0*xv) + 0.3f0
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
        # --- random ht component (SO DGSD=2.0 ≥ 1 ⇒ fires every tree, #206 cornered) ---
        zzran = 0f0
        if dgsd >= 1f0
            while true
                zzran = bachlo(s.rng, 0f0, 1f0)
                (zzran <= 0.5f0 && zzran >= -2.0f0) && break
            end
        end
        htgr = (htgr + zzran*0.1f0) * scale                # XRHGRO=1
        # --- blend small & large tree HTG ---
        xmn = SO_RG_XMIN[sp]; xmx = SO_RG_XMAX[sp]
        xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)
        t.ht_growth[i] = htgr*(1f0 - xwt) + xwt*t.ht_growth[i]
        t.ht_growth[i] < 0.1f0 && (t.ht_growth[i] = 0.1f0)
        cap = s.control.sp_size_cap[sp, 4]
        (h + t.ht_growth[i] > cap) && (t.ht_growth[i] = max(cap - h, 0.1f0))
        # --- small-tree DBH: only D < BKPT (3", 99 for WJ11); else keep large-tree DG ---
        bkpt = sp == 11 ? 99f0 : 3f0
        if d >= bkpt
            _ca_rg_stash!(stash, t, i); continue
        end
        htg = t.ht_growth[i]; hk = h + htg
        if hk <= 4.5f0
            t.diam_growth[i] = 0f0
            t.dbh[i] = d + 0.001f0*hk
            _ca_rg_stash!(stash, t, i); continue
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
            dk = (SO_RG_HT2_24/(log(hk - 4.5f0) - SO_RG_HT1_24)) - 1f0
            dkk = h <= 4.5f0 ? d : (SO_RG_HT2_24/(log(h - 4.5f0) - SO_RG_HT1_24)) - 1f0
        else                                               # all others: htdbh override (LHTDRG=false)
            dk = so_htdbh_dbh(ifor, sp, hk)
            dkk = h <= 4.5f0 ? d : so_htdbh_dbh(ifor, sp, h)
        end
        h < 4.5f0 && (dkk = d)                             # regent.f:605
        # --- DG (LESTB=false, regent.f:597-692) ---
        local dg::Float32
        if sp in SO_RG_HARDWOOD                            # WC hardwoods — double-DDS quirk, then DGBND
            dg = (dk < 0f0 || dkk < 0f0) ? htg*0.2f0*bark : (dk - dkk)*bark
            dg < 0f0 && (dg = 0.1f0); dg > dgmx && (dg = dgmx)
            dds = dg*(2f0*bark*d + dg)*scale2
            dg = sqrt((d*bark)^2 + dds) - bark*d
            dg < 0f0 && (dg = 0f0); dg > dgmx && (dg = dgmx)  # common block (681-688): re-converts
            dds = dg*(2f0*bark*d + dg)*scale2
            dg = sqrt((d*bark)^2 + dds) - bark*d
        elseif sp == 9 || sp == 27                         # SH/WO — XDWT blend; GO TO 23 skips DGBND
            xdwt = d <= 1.5f0 ? 0f0 : d >= 3f0 ? 1f0 : (d - 1.5f0)/1.5f0
            dgsm = (dk - dkk)*bark; dgsm < 0f0 && (dgsm = 0f0)
            dds = dgsm*(2f0*bark*d + dgsm)*scale2
            dgsm = sqrt((d*bark)^2 + dds) - bark*d; dgsm < 0f0 && (dgsm = 0f0)
            dg = dgsm*(1f0 - xdwt) + t.diam_growth[i]*xdwt
            (t.dbh[i] + dg) < SO_RG_DIAM[sp] && (dg = SO_RG_DIAM[sp] - t.dbh[i])
            t.diam_growth[i] = dg
            _ca_rg_stash!(stash, t, i); continue
        else                                               # DEFAULT
            dg = (dk < 0f0 || dkk < 0f0) ? htg*0.2f0*bark : (dk - dkk)*bark
            dg < 0f0 && (dg = 0f0); dg > dgmx && (dg = dgmx)
            (sp == 16 && (t.dbh[i] + dg) < SO_RG_DIAM[sp]) && (dg = SO_RG_DIAM[sp] - t.dbh[i])
            dds = dg*(2f0*bark*d + dg)*scale2
            dg = sqrt((d*bark)^2 + dds) - bark*d
        end
        (t.dbh[i] + dg) < SO_RG_DIAM[sp] && (dg = SO_RG_DIAM[sp] - t.dbh[i])
        dg = dg_bound(nothing, nothing, sp, t.dbh[i], dg, s.control.sp_size_cap)  # so/dgbnd.f = SIZCAP cap
        t.diam_growth[i] = dg
        _ca_rg_stash!(stash, t, i)
    end
    return s
end
