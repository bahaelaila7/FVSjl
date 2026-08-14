# =============================================================================
# small_tree_growth.jl (westsierra) — WS small-tree HEIGHT+DBH (ws/regent.f + ws/smhtgf.f). Chunk 6.
#
# Mirrors SO small_tree_growth.jl (regent DBH driver + smhtgf height increment), with WS specifics:
#   • ws_smhtgf: species-GROUPED HTGR (raw ISPC groups: pines/firs/oak/tanoak/LP-WB/KP/RW-GS). CR = ICR/10.
#   • ★ WS-native CASE DEFAULT species use HTGRR = ws_smhtgf DIRECTLY (·CON) — NO PCTRED·VIGOR (density enters
#     via BAL inside smhtgf). ONLY MC(41)/GB(21) apply POTHTG·PCTRED·VIGOR. (SO applies PCTRED·VIGOR to DEFAULT;
#     WS does NOT — a load-bearing measure catch, ws/regent.f:311-347.)
#   • DK/DKK by MSP=SMTMAP group: MSP1 pines AX=−0.6197/BX=0.2626 linear (DK=AX+BX·HK); MSP2 firs −0.6096/0.2433
#     linear; MSP3,4 oak/tanoak AX=4.80420/BX=−9.92422 ln-form (DK=BX/(ln(HK−4.5)−AX)−1). CA-surrogate (4,9:10,
#     12,14:17,19:20,23,25:27) via ws_htdbh_dbh (chunk 4b, fires even LHTDRG=false); MC(41)/GB(21) own forms.
#   • DG: XDWT=(D−1.5)/1.5 clamp[0,1] blend of small (DGSM from DK−DKK) and large (DGLT) DG; then DGBND (SIZCAP).
# MEASURED vs FVSws_g16 regent per-tree HTGR + DG on wst01 (SP/DF/WF/RF + WB/TO regen).
# =============================================================================

# ws/regent.f DATA (index = ISPC).
const WS_RG_XMAX = Float32[
  3.5,3.5,3.5,10.0,3.5, 3.5,3.5,3.5,4.0,4.0, 3.5,4.0,3.5,4.0,4.0, 4.0,4.0,3.5,4.0,4.0,
  199.,3.5,10.0,3.5,4.0, 4.0,4.0,3.5,3.5,3.5, 3.5,3.5,3.5,3.5,3.5, 3.5,3.5,3.5,3.5,3.5, 4.0,3.5,3.5]
const WS_RG_XMIN = Float32[
  2,2,2,2,2, 2,2,2,2,2, 2,2,2,2,2, 2,2,2,2,2, 99,2,2,2,2, 2,2,1,1,1, 1,1,1,2,2, 2,2,2,2,1, 2,2,1]
const WS_RG_DIAM = Float32[
  0.3,0.3,0.3,0.2,0.2, 0.3,0.3,0.5,0.4,0.4, 0.3,0.5,0.3,0.5,0.5, 0.5,0.5,0.5,0.5,0.5,
  0.4,0.3,0.2,0.3,0.5, 0.5,0.5,0.4,0.4,0.4, 0.4,0.4,0.4,0.2,0.2, 0.2,0.2,0.2,0.2,0.4, 0.2,0.4,0.4]
# ws/regent.f SMTMAP (species → smhtgf/DK group: 1=pines,2=firs,3=black-oak,4=tanoak; 0=CA/RW-GS-other).
const WS_RG_SMTMAP = Int[
  1,2,2,0,1, 1,2,1,0,0, 1,0,2,0,0, 0,0,1,0,0, 0,2,0,1,0, 0,0,3,3,3, 3,3,3,4,4, 4,4,4,4,3, 0,1,3]
const WS_RG_AB = Float32[1.11436, -0.011493, 0.43012f-4, -0.72221f-7, 0.5607f-10, -0.1641f-13]  # PCTRED poly (=SO)

# ws/smhtgf.f species groups (raw ISPC).
const WS_SM_PINE    = Set([1,5,6,8,11,18,24,42])
const WS_SM_PINE175 = Set([1,5,11,24])
const WS_SM_FIR     = Set([2,3,7,13,22])
const WS_SM_OAK     = Set([28,29,30,31,32,33,40,43])
const WS_SM_TANOAK  = Set([34,35,36,37,38,39])
const WS_SM_CASURR  = Set([9,10,12,14,15,16,17,19,20,25,26,27])   # LP/WB + KP-group (same smhtgf exp·factor·1.75)
# regent CA-surrogate species that use ws_htdbh_dbh for DK/DKK.
const WS_RG_CASURR  = Set([4,9,10,12,14,15,16,17,19,20,23,25,26,27])

# ws/smhtgf.f — small-tree HEIGHT increment. nspc = raw ISPC; cr = ICR/10 (regent passes ICR/10). floor 0.1.
@inline function ws_smhtgf(nspc::Int, d::Float32, cr::Float32, ba::Float32, bal::Float32, si::Float32, h::Float32)::Float32
    tembal = bal < 5f0 ? 5f0 : bal
    factor = 0.80f0 + 0.004f0 * (si - 50f0)
    htgr = 0f0
    if nspc in WS_SM_PINE
        htgr = exp(0.7452f0 - 0.003271f0*bal - 0.1632f0*cr + 0.0217f0*cr*cr + 0.00536f0*si)
        htgr *= (nspc in WS_SM_PINE175 ? 1.75f0 : 1.50f0) * factor
    elseif nspc in WS_SM_FIR
        htgr = exp(-0.2495f0 - 0.00111f0*bal + 0.0100f0*cr*cr)
        htgr = (nspc == 2 || nspc == 22 ? (htgr + 1f0)*2.5f0 : (htgr + 0.75f0)*2.0f0) * factor
    elseif nspc in WS_SM_OAK
        htgr = exp(3.817f0 - 0.7829f0*log(tembal))
    elseif nspc in WS_SM_TANOAK
        htgr = exp(3.385f0 - 0.5898f0*log(tembal))
    elseif nspc in WS_SM_CASURR
        htgr = exp(0.7452f0 - 0.003271f0*bal - 0.1632f0*cr + 0.0217f0*cr*cr + 0.00536f0*si) * factor * 1.75f0
    elseif nspc == 4 || nspc == 23                              # RW/GS Chapman-Richards
        htmax = 2.242202f0 * si
        if htmax - h <= 1f0
            htgr = 0f0
        else
            age1 = (1f0 / -0.010742f0) * log(1f0 - fpow(h/2.242202f0/si, 1f0/0.919076f0))
            age2 = age1 + 5f0
            h1 = 2.242202f0*si*fpow(1f0 - exp(-0.010742f0*age1), 0.919076f0)
            h2 = 2.242202f0*si*fpow(1f0 - exp(-0.010742f0*age2), 0.919076f0)
            htgr = h2 - h1
        end
    end
    htgr < 0.1f0 && (htgr = 0.1f0)
    return htgr
end

# DK/DKK helper — breast-height crossing DBH for a target height (regent.f:419-511, MSP-grouped for WS-native).
@inline function _ws_regent_dk(ifor::Int, sp::Int, msp::Int, hk::Float32, sitear::Float32)::Float32
    if sp in WS_RG_CASURR                                       # CA-surrogate → ws_htdbh (chunk 4b)
        return ws_htdbh_dbh(ifor, sp, hk)
    elseif sp == 41                                             # MC linear
        return 3.1020f0 + 0.0210f0 * hk
    elseif sp == 21                                             # GB linear
        dk = (hk - 4.5f0) * 10f0 / (sitear - 4.5f0)
        return dk < 0.1f0 ? 0.1f0 : dk
    else                                                        # WS-native MSP groups
        if msp == 1                                             # pines
            return -0.6197f0 + 0.2626f0 * hk
        elseif msp == 2                                         # firs
            return -0.6096f0 + 0.2433f0 * hk
        else                                                    # msp 3/4 oak/tanoak ln-form
            return -9.92422f0 / (log(hk - 4.5f0) - 4.80420f0) - 1f0
        end
    end
end

function small_tree_growth!(s::StandState, stash, ::WestSierra; fint::Float32 = 10.0f0)
    p, t, c = s.plot, s.trees, s.calib
    n = t.n; n == 0 && return s
    sd = s.coef.species
    avh = stand_top_height(s); ba = p.basal_area; dgsd = s.control.dg_sd
    relden = p.relative_density; ifor = Int(p.forest_idx); yr = s.control.year
    fnt = fint
    # PCTRED density modifier (ws/regent.f:243-248) — used only by MC(41)/GB(21).
    xden = avh * (relden / 100f0); xden > 300f0 && (xden = 300f0)
    pctred = WS_RG_AB[1] + xden*(WS_RG_AB[2] + xden*(WS_RG_AB[3] + xden*(WS_RG_AB[4] +
             xden*(WS_RG_AB[5] + xden*WS_RG_AB[6]))))
    pctred > 1f0 && (pctred = 1f0); pctred < 0.01f0 && (pctred = 0.01f0)
    @inbounds for i in 1:n
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= WS_RG_XMAX[sp] || t.tpa[i] <= 0f0) && continue
        h = t.height[i]
        icr = Float32(t.crown_pct[i])
        crf = icr / 10f0                                        # CR = ICR/10 passed to smhtgf
        si = p.sp_site_index[sp]
        con = exp(c.htg_cor_small[sp])                          # RHCON·exp(HCOR); HCOR=0 ⇒ 1
        regyr = 10f0
        scale = fnt / regyr; scale2 = yr / fnt
        msp = WS_RG_SMTMAP[sp]
        # --- HTGRR: MC/GB use POTHTG·PCTRED·VIGOR; WS-native use ws_smhtgf DIRECTLY (no PCTRED·VIGOR) ---
        local htgrr::Float32
        if sp == 41                                            # MC
            xv = icr / 100f0; vigor = 150f0*fpow(xv,3f0)*exp(-6f0*xv) + 0.3f0; vigor > 1f0 && (vigor = 1f0)
            pothtg = ((1.47043f0 + 0.23317f0*si) / (31.56252f0 - 0.05586f0*si)) * 10f0
            htgrr = pothtg * pctred * vigor
        elseif sp == 21                                        # GB
            xv = icr / 100f0; vigor = 150f0*fpow(xv,3f0)*exp(-6f0*xv) + 0.3f0; vigor > 1f0 && (vigor = 1f0)
            vigor = 1f0 - (1f0 - vigor)/3f0
            pothtg = ((si/5f0) * (si*1.5f0 - h)/(si*1.5f0)) * 0.83f0
            htgrr = pothtg * pctred * vigor
        else                                                   # WS-native CASE DEFAULT — smhtgf direct
            bal = ((100f0 - t.crown_ratio[i]) / 100f0) * ba
            htgrr = ws_smhtgf(sp, d, crf, ba, bal, si, h)
        end
        htgr = htgrr * con
        # --- random ht component (WS DGSD=2.0 ≥ 1 ⇒ fires; redraw if >0.5 or <−2) ---
        zzran = 0f0
        if dgsd >= 1f0
            while true
                zzran = bachlo(s.rng, 0f0, 1f0)
                (zzran <= 0.5f0 && zzran >= -2.0f0) && break
            end
        end
        htgr = sp == 21 ? (htgr + zzran*0.2f0) * scale : (htgr + zzran*0.1f0) * scale   # XRHGRO=1 (GB ·WK4 omitted)
        # --- blend small & large tree HTG ---
        xmn = WS_RG_XMIN[sp]; xmx = WS_RG_XMAX[sp]
        xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)
        t.ht_growth[i] = htgr*(1f0 - xwt) + xwt*t.ht_growth[i]
        t.ht_growth[i] < 0.1f0 && (t.ht_growth[i] = 0.1f0)
        cap = s.control.sp_size_cap[sp, 4]
        (cap > 0f0 && h + t.ht_growth[i] > cap) && (t.ht_growth[i] = max(cap - h, 0.1f0))
        # --- small-tree DBH: D < XMAX only (already gated); DK/DKK breast-height crossing ---
        htg = t.ht_growth[i]; hk = h + htg
        bark = ws_bratio(sd, sp, d)
        dk = _ws_regent_dk(ifor, sp, msp, hk, si)
        dkk = h <= 4.5f0 ? d : _ws_regent_dk(ifor, sp, msp, h, si)
        # --- DG conversion (LESTB=false): XDWT blend of small (DK−DKK) and large (DGLT) DG ---
        xdwt = d <= 1.5f0 ? 0f0 : d >= 3f0 ? 1f0 : (d - 1.5f0)/1.5f0
        dgsm = (dk - dkk) * bark
        dgsm < 0f0 && (dgsm = 0f0)
        dds = dgsm*(2f0*bark*d + dgsm)*scale2
        dgsm = sqrt((d*bark)^2 + dds) - bark*d
        dgsm < 0f0 && (dgsm = 0f0)
        dglt = t.diam_growth[i]
        dg = dgsm*(1f0 - xdwt) + dglt*xdwt
        (t.dbh[i] + dg) < WS_RG_DIAM[sp] && (dg = WS_RG_DIAM[sp] - t.dbh[i])
        dg = dg_bound(nothing, nothing, sp, t.dbh[i], dg, s.control.sp_size_cap)   # ws/dgbnd.f = SIZCAP
        t.diam_growth[i] = dg
        _ca_rg_stash!(stash, t, i)
    end
    return s
end
