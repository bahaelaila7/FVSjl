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
"""ws/regent.f:417-541 — (DK, DKK): the DBH the tree would have at HK (end) and at H (start).
First the species' own form: CA-surrogate/RW/GS Wykoff inverse HT2/(ln(H−4.5)−AX)−1 (AX = calibrated AA;
HT1 only matters when IABFLG=1, where the HTDBH override below replaces it anyway); MC(41) linear
3.1020+0.0210·H with DK≥DKK+0.01; GB(21) the CR/UT site line (H−4.5)·10/(SITEAR−4.5) floored at 0.1; WS-native
MSP pine/fir linear or the oak ln-form. DKK=D below breast height. Then for the CA/SO species (4,9:10,12,
14:17,19:20,23,25:27,41) the inventory HTDBH replaces DK/DKK when .NOT.LHTDRG or (LHTDRG and IABFLG=1)."""
function _ws_regent_dk_dkk(s::StandState, ifor::Int, sp::Int, msp::Int, d::Float32, h::Float32, hk::Float32,
                           sitear::Float32)
    local dk::Float32, dkk::Float32
    if sp in WS_RG_CASURR || sp == 4 || sp == 23
        bx = coef_col(s.coef, :wykoff_ht2)[sp]; ax = s.calib.ht_dbh_aa[sp]
        dk = (bx / (log(hk - 4.5f0) - ax)) - 1f0
        dkk = h <= 4.5f0 ? d : (bx / (log(h - 4.5f0) - ax)) - 1f0
    elseif sp == 41
        dkk = 3.1020f0 + 0.0210f0 * h; dkk < 0f0 && (dkk = d)
        dk = 3.1020f0 + 0.0210f0 * hk; dk < dkk && (dk = dkk + 0.01f0)
    elseif sp == 21
        dk = (hk - 4.5f0) * 10f0 / (sitear - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
        dkk = (h - 4.5f0) * 10f0 / (sitear - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0)
        h < 4.5f0 && (dkk = d)
    else
        ihdw = msp == 3 || msp == 4
        ax, bx = msp == 1 ? (-0.6197f0, 0.2626f0) : msp == 2 ? (-0.6096f0, 0.2433f0) : (4.80420f0, -9.92422f0)
        dk = ihdw ? bx / (log(hk - 4.5f0) - ax) - 1f0 : ax + bx * hk
        dkk = h <= 4.5f0 ? d : (ihdw ? bx / (log(h - 4.5f0) - ax) - 1f0 : ax + bx * h)
    end
    if (sp in WS_RG_CASURR || sp == 4 || sp == 23 || sp == 41) &&
       (!s.control.ht_drag_sp[sp] || s.calib.ht_dbh_iabflg[sp] == 1)
        dk = ws_htdbh_dbh(ifor, sp, hk)
        dkk = h <= 4.5f0 ? d : ws_htdbh_dbh(ifor, sp, h)
    end
    return (dk, dkk)
end

function small_tree_growth!(s::StandState, stash, ::WestSierra; fint::Float32 = 10.0f0)
    p, t, c = s.plot, s.trees, s.calib
    n = t.n; n == 0 && return s
    cw = clim_wk4(s, Float32(current_cycle_year(s)) + fint / 2f0)   # CLGMULT WK4 (ws/regent.f:359); nothing ⇒ 1
    sd = s.coef.species
    avh = stand_top_height(s); ba = p.basal_area; dgsd = s.control.dg_sd
    relden = p.relative_density; ifor = Int(p.forest_idx); yr = s.control.year
    fnt = fint
    # PCTRED density modifier (ws/regent.f:243-248) — used only by MC(41)/GB(21).
    xden = avh * (relden / 100f0); xden > 300f0 && (xden = 300f0)
    pctred = WS_RG_AB[1] + xden*(WS_RG_AB[2] + xden*(WS_RG_AB[3] + xden*(WS_RG_AB[4] +
             xden*(WS_RG_AB[5] + xden*WS_RG_AB[6]))))
    pctred > 1f0 && (pctred = 1f0); pctred < 0.01f0 && (pctred = 0.01f0)
    # ws/regent.f:255-277 is SPECIES-MAJOR (DO 30 ISPC … I=IND1(I3)); the per-tree ZZRAN draw must follow it.
    @inbounds for i in species_major_order(s)
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= WS_RG_XMAX[sp] || t.tpa[i] <= 0f0) && continue
        h = t.height[i]
        icr = Float32(t.crown_pct[i])
        crf = icr / 10f0                                        # CR = ICR/10 passed to smhtgf
        si = p.sp_site_index[sp]
        con = exp(c.htg_cor_small[sp])                          # RHCON·exp(HCOR); HCOR=0 ⇒ 1
        # ws/regent.f:265-270 — GB(21)/MC(41) equations are 10-yr (UT/SO); every WS-native species' SMHTGF
        # increment is 5-yr (CA/WS). A flat 10 halved every native small tree's height growth on 10-yr cycles
        # (measured vs FVSws_g16, DGSTDEV 0: OS/WF/DF HtG 7.46/7.30/10.01 vs jl 3.73/3.65/5.01).
        regyr = (sp == 21 || sp == 41) ? 10f0 : 5f0
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
        # ws/regent.f:357-362: GB (21) ·WK4(I) = the CLGMULT climate multiplier (1 without CLIMATE). XRHGRO=1.
        htgr = sp == 21 ? (htgr + zzran*0.2f0) * scale * (cw === nothing ? 1f0 : cw[i]) : (htgr + zzran*0.1f0) * scale
        # --- blend small & large tree HTG ---
        xmn = WS_RG_XMIN[sp]; xmx = WS_RG_XMAX[sp]
        xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)
        t.ht_growth[i] = htgr*(1f0 - xwt) + xwt*t.ht_growth[i]
        t.ht_growth[i] < 0.1f0 && (t.ht_growth[i] = 0.1f0)
        cap = s.control.sp_size_cap[sp, 4]
        (cap > 0f0 && h + t.ht_growth[i] > cap) && (t.ht_growth[i] = max(cap - h, 0.1f0))
        # --- small-tree DBH (ws/regent.f:405-665) ---
        # BKPT: GB grows its DBH in REGENT at every size (99"), RW/GS below 7", everyone else below 3". At or
        # above BKPT FVS jumps to label 23: the large-tree DG stays untouched (no DIAM floor, no DGBND).
        bkpt = sp == 21 ? 99f0 : (sp == 4 || sp == 23) ? 7f0 : 3f0
        if d >= bkpt
            _ca_rg_stash!(stash, t, i); continue
        end
        htg = t.ht_growth[i]; hk = h + htg
        bark = ws_bratio(sd, sp, d)
        if hk <= 4.5f0
            t.diam_growth[i] = 0f0; t.dbh[i] = d + 0.001f0 * hk           # regent.f:413-415
        else
            dk, dkk = _ws_regent_dk_dkk(s, ifor, sp, msp, d, h, hk, si)
            local dg::Float32
            if sp == 21                                               # UT (orig. CR): DG by subtraction, ×SCALE2
                dgmx = 2f0 * scale
                if dk < 0f0 || dkk < 0f0
                    dg = htg * 0.2f0 * bark
                else
                    dg = (dk - dkk) * bark                            # ×XRDGRO (=1)
                end
                dg < 0f0 && (dg = 0f0); dg > dgmx && (dg = dgmx)
                dds = dg * (2f0 * bark * d + dg) * scale2
                dg = sqrt((d * bark)^2 + dds) - bark * d
            elseif sp == 41                                           # SO (orig. WC): no DDS rescale
                h < 4.5f0 && (dkk = d)
                dgmx = 5f0 * scale
                dg = (dk < 0f0 || dkk < 0f0) ? htg * 0.2f0 * bark : (dk - dkk) * bark
                (s.control.ht_drag_sp[41] && s.calib.ht_dbh_iabflg[41] == 0) && (dg = 0.1f0 * htg)
                dg < 0f0 && (dg = 0.1f0); dg > dgmx && (dg = dgmx)
            else                                                      # CA/WS: blend small & large-tree DG
                xdwt = if sp == 4 || sp == 23
                    xmn = WS_RG_XMIN[sp]; d <= xmn ? 0f0 : (d - xmn) / (bkpt - xmn)
                else
                    d <= 1.5f0 ? 0f0 : d >= 3f0 ? 1f0 : (d - 1.5f0) / 1.5f0
                end
                dgsm = (dk - dkk) * bark
                dgsm < 0f0 && (dgsm = 0f0)
                dds = dgsm * (2f0 * bark * d + dgsm) * scale2
                dgsm = sqrt((d * bark)^2 + dds) - bark * d
                dgsm < 0f0 && (dgsm = 0f0)
                dg = dgsm * (1f0 - xdwt) + t.diam_growth[i] * xdwt
            end
            (t.dbh[i] + dg) < WS_RG_DIAM[sp] && (dg = WS_RG_DIAM[sp] - t.dbh[i])
            t.diam_growth[i] = dg
        end
        t.diam_growth[i] = dg_bound(nothing, nothing, sp, t.dbh[i], t.diam_growth[i], s.control.sp_size_cap)   # DGBND
        _ca_rg_stash!(stash, t, i)
    end
    return s
end
