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

# ws/regent.f label 4 → 23 for one slot K: returns (DBH(K) set directly (<0 ⇒ none), DG(K)). D≥BKPT ⇒ GO TO 23
# (DG(K) = the slot's DGDRIV DG `dglt`); HK≤4.5 ⇒ DG=0, DBH(K)=D+0.001·HK; else the species DG, the DIAM floor
# and DGBND. DBH(K)=D for every slot (dgdriv.f tripling copies DBH(I)).
function _ws_regent_dg(s::StandState, ifor::Int, sp::Int, msp::Int, d::Float32, h::Float32, htg::Float32,
                       si::Float32, bkpt::Float32, scale::Float32, scale2::Float32, dglt::Float32, xrdgro::Float32)
    d >= bkpt && return (-1f0, dglt)
    hk = h + htg
    dbhk = d
    local dg::Float32
    if hk <= 4.5f0
        dg = 0f0; dbhk = d + 0.001f0 * hk                              # regent.f:413-415
    else
        bark = ws_bratio(s.coef.species, sp, d)
        dk, dkk = _ws_regent_dk_dkk(s, ifor, sp, msp, d, h, hk, si)
        if sp == 21                                               # UT (orig. CR): DG by subtraction, ×SCALE2
            dgmx = 2f0 * scale
            dg = (dk < 0f0 || dkk < 0f0) ? htg * 0.2f0 * bark * xrdgro : (dk - dkk) * bark * xrdgro
            dg < 0f0 && (dg = 0f0); dg > dgmx && (dg = dgmx)
            dds = dg * (2f0 * bark * d + dg) * scale2
            dg = sqrt((d * bark)^2 + dds) - bark * d
        elseif sp == 41                                           # SO (orig. WC): no DDS rescale
            h < 4.5f0 && (dkk = d)
            dgmx = 5f0 * scale
            dg = (dk < 0f0 || dkk < 0f0) ? htg * 0.2f0 * bark * xrdgro : (dk - dkk) * bark * xrdgro
            (s.control.ht_drag_sp[41] && s.calib.ht_dbh_iabflg[41] == 0) && (dg = 0.1f0 * htg * xrdgro)
            dg < 0f0 && (dg = 0.1f0); dg > dgmx && (dg = dgmx)
        else                                                      # CA/WS: blend small & large-tree DG
            xdwt = if sp == 4 || sp == 23
                xmn = WS_RG_XMIN[sp]; d <= xmn ? 0f0 : (d - xmn) / (bkpt - xmn)
            else
                d <= 1.5f0 ? 0f0 : d >= 3f0 ? 1f0 : (d - 1.5f0) / 1.5f0
            end
            dgsm = (dk - dkk) * bark * xrdgro
            dgsm < 0f0 && (dgsm = 0f0)
            dds = dgsm * (2f0 * bark * d + dgsm) * scale2
            dgsm = sqrt((d * bark)^2 + dds) - bark * d
            dgsm < 0f0 && (dgsm = 0f0)
            dg = dgsm * (1f0 - xdwt) + dglt * xdwt
        end
        (d + dg) < WS_RG_DIAM[sp] && (dg = WS_RG_DIAM[sp] - d)
    end
    dg = dg_bound(nothing, nothing, sp, dbhk, dg, s.control.sp_size_cap)   # DGBND(ISPC,DBH(K),DG(K))
    return (hk <= 4.5f0 ? dbhk : -1f0, dg)
end

function small_tree_growth!(s::StandState, stash, ::WestSierra; fint::Float32 = 10.0f0)
    p, t, c = s.plot, s.trees, s.calib
    n = t.n; n == 0 && return s
    cw = clim_wk4(s, Float32(current_cycle_year(s)) + fint / 2f0)   # CLGMULT WK4 (ws/regent.f:359); nothing ⇒ 1
    sd = s.coef.species
    # AVH = the COMMON AVHT40 of the last DENSE (cycle start: CRATET's IND at cycle 0, gradd.f:186's after)
    avh = p.avg_height; ba = p.basal_area; dgsd = s.control.dg_sd
    relden = p.relative_density; ifor = Int(p.forest_idx); yr = s.control.year
    fnt = fint
    # PCTRED density modifier (ws/regent.f:243-248) — used only by MC(41)/GB(21).
    xden = avh * (relden / 100f0); xden > 300f0 && (xden = 300f0)
    pctred = WS_RG_AB[1] + xden*(WS_RG_AB[2] + xden*(WS_RG_AB[3] + xden*(WS_RG_AB[4] +
             xden*(WS_RG_AB[5] + xden*WS_RG_AB[6]))))
    pctred > 1f0 && (pctred = 1f0); pctred < 0.01f0 && (pctred = 0.01f0)
    trip = stash !== nothing && !isempty(stash.dgU)
    yr_now = current_cycle_year(s)                              # ws/regent.f MULTS(3/6, IY(ICYC))
    # ws/regent.f:255-277 is SPECIES-MAJOR (DO 30 ISPC … I=IND1(I3)); the per-tree ZZRAN draw must follow it.
    @inbounds for i in species_major_order(s)
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= WS_RG_XMAX[sp] || t.tpa[i] <= 0f0) && continue
        h = t.height[i]
        icr = Float32(t.crown_pct[i])
        crf = icr / 10f0                                        # CR = ICR/10 passed to smhtgf
        si = p.sp_site_index[sp]
        con = exp(c.htg_cor_small[sp])                          # RHCON·exp(HCOR); HCOR=0 ⇒ 1
        xrhgro = active_multiplier(s.control, :regh, sp, yr_now)   # XRHGRO=XRHMLT(ISPC) (REGHMULT)
        xrdgro = active_multiplier(s.control, :regd, sp, yr_now)   # XRDGRO=XRDMLT(ISPC) (REGDMULT)
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
        xmn = WS_RG_XMIN[sp]; xmx = WS_RG_XMAX[sp]
        xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)
        cap = s.control.sp_size_cap[sp, 4]
        # BKPT: GB grows its DBH in REGENT at every size (99"), RW/GS below 7", everyone else below 3". At or
        # above BKPT FVS jumps to label 23: the large-tree DG stays untouched (no DIAM floor, no DGBND).
        bkpt = sp == 21 ? 99f0 : (sp == 4 || sp == 23) ? 7f0 : 3f0
        wk4 = sp == 21 ? (cw === nothing ? 1f0 : cw[i]) : 1f0
        large_htg = t.ht_growth[i]; dg_main = t.diam_growth[i]
        # ws/regent.f:666-669 — with LTRIP each tripled copy K=ITRN+2I−2+L (L=1,2) reruns label 2 → 23: a FRESH
        # ZZRAN, the blend with the copy's own HTG(K), SIZCAP and (D<BKPT) its own DG(K)/DBH(K). H and D stay
        # record I's; DBH(K)=DBH(I) (dgdriv.f). HTGRR (above label 2) is computed once per record.
        for l in 0:(trip ? 2 : 0)
            lthg = l == 0 ? large_htg : (stash.htg_copy[i] ? (l == 1 ? stash.htgU[i] : stash.htgL[i]) : large_htg)
            dglt = l == 0 ? dg_main : (l == 1 ? stash.dgU[i] : stash.dgL[i])
            htgr = htgrr * con                                  # label 2
            # --- random ht component (label 3; WS DGSD=2.0 ≥ 1 ⇒ fires; redraw if >0.5 or <−2) ---
            zzran = 0f0
            if dgsd >= 1f0
                while true
                    zzran = bachlo(s.rng, 0f0, 1f0)
                    (zzran <= 0.5f0 && zzran >= -2.0f0) && break
                end
            end
            # ws/regent.f:357-362: GB (21) ·WK4(I) = the CLGMULT climate multiplier (1 without CLIMATE).
            htgr = sp == 21 ? (htgr + zzran*0.2f0) * xrhgro * scale * wk4 : (htgr + zzran*0.1f0) * xrhgro * scale
            # --- blend small & large tree HTG (ws/regent.f:377-389; RW/GS average with the large-tree HTG first) ---
            local htg::Float32
            if sp == 4 || sp == 23
                htgr = (htgr + lthg) / 2f0
                htg = htgr*(1f0 - xwt) + xwt*lthg
            else
                htg = htgr*(1f0 - xwt) + xwt*lthg
            end
            htg < 0.1f0 && (htg = 0.1f0)
            (cap > 0f0 && h + htg > cap) && (htg = max(cap - h, 0.1f0))
            dbhk, dg = _ws_regent_dg(s, ifor, sp, msp, d, h, htg, si, bkpt, scale, scale2, dglt, xrdgro)
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

"""ws/regent.f LSTART small-tree HEIGHT calibration (label 40; called from cratet.f:770). PCTRED from
X=AVH·RELDEN/100 (AVH = cratet.f:728 AVHT40, RELDEN = the cratet.f:247 backdating DENSE's). Per species (IND1
order): current DBH<5, backdated H=HT−HTG (IHTG<2) ≥0.01, measured HTG≥0.001 ⇒ EDH·RHCON with EDH = MC(41)
POTHTG·PCTRED·VIGOR, GB(21) the same with the 2/3 VIGOR knock-down, else SMHTGF(D, ICR/10, BA, BAL, SITEAR,
H_back); TERM = HTG·SCALE3 (SCALE3 = REGYR/FINTH, REGYR 10 for 21/41 else 5); CORNEW = mean TERM / mean EDH when
N≥NCALHT(5), ≤0 ⇒ 1E−4, trapped to [0.0821, 12.1825] else 1. HCOR = ln(CORNEW) → htg_cor_init."""
function ws_regent_hcor_init!(s::StandState, isct::AbstractMatrix, ind1::AbstractVector,
                              saved_dbh::AbstractVector, avh::Float32)
    p, t, c = s.plot, s.trees, s.calib
    t.n == 0 && return s
    s.control.growth_ifinth == 0 && return s
    ba = c.cratet_ba; relden = c.cratet_relden
    finth = s.control.growth_finth > 0f0 ? s.control.growth_finth : 5f0
    x = avh * (relden / 100f0); x > 300f0 && (x = 300f0)
    pctred = WS_RG_AB[1] + x*(WS_RG_AB[2] + x*(WS_RG_AB[3] + x*(WS_RG_AB[4] + x*(WS_RG_AB[5] + x*WS_RG_AB[6]))))
    pctred > 1f0 && (pctred = 1f0); pctred < 0.01f0 && (pctred = 0.01f0)
    @inbounds for sp in 1:size(isct, 1)
        i1 = Int(isct[sp, 1]); i1 == 0 && continue
        i2 = Int(isct[sp, 2])
        si = p.sp_site_index[sp]
        regyr = (sp == 21 || sp == 41) ? 10f0 : 5f0
        scale3 = regyr / finth
        rhcon = (s.control.regh_cor2_on && s.control.regh_cor2[sp] > 0f0) ? s.control.regh_cor2[sp] : 1f0
        snp = 0f0; snx = 0f0; sny = 0f0; nh = 0
        for k in i1:i2
            i = Int(ind1[k])
            d = saved_dbh[i]; h = t.height[i]
            pct = i <= length(c.cratet_pct) ? c.cratet_pct[i] : t.crown_ratio[i]
            bal = ((100f0 - pct) / 100f0) * ba
            icr = Float32(t.crown_pct[i]); cr = icr / 10f0
            s.control.growth_ihtg < 2 && (h = h - t.ht_growth[i])
            (d >= 5f0 || h < 0.01f0) && continue
            xv = icr / 100f0
            vigor = 150f0 * fpow(xv, 3f0) * exp(-6f0*xv) + 0.3f0
            vigor > 1f0 && (vigor = 1f0)
            local edh::Float32
            if sp == 41
                edh = ((1.47043f0 + 0.23317f0*si) / (31.56252f0 - 0.05586f0*si)) * 10f0 * pctred * vigor
            elseif sp == 21
                vigor = 1f0 - (1f0 - vigor)/3f0
                edh = ((si/5f0) * (si*1.5f0 - h)/(si*1.5f0)) * 0.83f0 * pctred * vigor
            else
                edh = ws_smhtgf(sp, d, cr, ba, bal, si, h)
            end
            edh *= rhcon
            hg = t.ht_growth[i]; hg < 0.001f0 && continue
            pr = t.tpa[i]
            snp += pr; snx += edh * pr; sny += hg * scale3 * pr; nh += 1
        end
        nh < 5 && continue
        snx /= snp; sny /= snp
        cornew = sny / snx
        cornew <= 0f0 && (cornew = 1f-4)
        (cornew < 0.0821f0 || cornew > 12.1825f0) && (cornew = 1f0)
        c.htg_cor_init[sp] = log(cornew)
    end
    return s
end
