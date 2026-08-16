# =============================================================================
# diameter_growth.jl (olympic) — OP FVS-NATIVE large-tree DDS (op/dgf.f) for the
# NON-ORGANON species. Chunk 1a.
#
# In OP only GF(3)/DF(16) are ORGANON-grown; the valid-ORGANON set is
# {3,16,18,19,21,22,23,28,33,34,37} (op/dgdriv.f). EVERY OTHER species (LP/PP/SP/ES/
# WF/…) grows the FVS-native op Wykoff ln(DDS) here — a standard western Wykoff DDS,
# WC-derived (op/dgf.f cites the WC/CA data). This ports op/dgf.f verbatim:
#   op_dgcons(ispc,jspc,ifor,elev,slope,asp,sitear) — ENTRY DGCONS (site-dependent DGCON)
#   op_dgf_dds(...)  — the per-tree LN(DDS) DEFAULT branch (+ RA/RW special forms)
#   op_dgcons!(s) / dgf!(s,::Olympic) — StandState wrappers (ORGANON trees are grown by
#     the NWO engine — this path only touches the FVS-native IORG=0 trees).
#
# op/dgf.f maps 39 FVS species → 20 coefficient groups via MAPSPC. DDS DEFAULT (op/dgf.f:466):
#   DDS = CONSPP + DGLD·ln(D) + CR·(DGCR + CR·DGCRSQ) + DGDSQ·D² + DGDBAL·BAL/ln(D+1)
#       + DGPCCF·PCCF + DGHAH·RELHT + DGLBA·ln(BA) + DGBAL·BAL + DGBA·BA
#   TDDS = EXP(DDS); DDS = LN(TDDS/2) (10-yr→5-yr real-scale); clamp DDS ≥ −9.21.
# CONSPP = DGCON(ISPC) + COR(ISPC) (COR=0 in the calibration pass; RW(17) defers COR).
# DGCONS (op/dgf.f:501): DGCON(ISPC) = DGFOR(ISPFOR,JSPC) + DGEL·TEMEL + DGEL2·TEMEL²
#   + DGSITE·ln(XSITE) + SASP, with SASP = (DGSASP·sinASP + DGCASP·cosASP + DGSLOP)·SLOPE
#   + DGSLSQ·SLOPE²; XSITE=SITEAR(ISPC) (WH JSPC10 ×3.281, WO JSPC19 log-transform);
#   TEMEL=ELEV (JSPC14 capped 30); ISPFOR=MAPLOC(IFOR,JSPC), ISPDSQ=MAPDSQ(IFOR,JSPC),
#   DGDSQ(JSPC)=DGDS(ISPDSQ,JSPC).
#
# VALIDATED per-tree vs live FVSop_clean (stand S248112, DEBUG DGF dump, cyc0): the 19
# FVS-native trees (WF/ES/LP/SP/PP) reproduce the oracle LN(DDS) to its F7.4 print — see
# test/unit/test_op_native_growth.jl. op/grinit.f sets DGSD=0 and ICL4=0, so OP has NO OLDRN
# serial-correlation and NO record tripling — it is fully DETERMINISTIC, and every cycle's DDS
# is a genuine bit-exact target (not a straddle). Multi-cycle validation is blocked oracle-side:
# FVSop_clean SIGSEGVs in the FVS-native fvsvol path on cycle ≥1.
# =============================================================================

# =============================================================================
# op/htdbh.f — the CRATET missing-height dub (Curtis-Arney HT-DBH). op/htdbh.f is BYTE-IDENTICAL to
# pn/htdbh.f (verified), so the 6-forest × 39-species P2/P3/P4 tables are REUSED from PN
# (PN_HTDBH_P2/P3/P4, regent.jl). LHTDRG is .FALSE. for all OP species (op/grinit.f:106), so
# missing-height dubbing ALWAYS routes through HTDBH MODE=0 (op/cratet.f:692-694). Faithful to the
# htdbh.f MODE=0 branch: D≥3 Curtis-Arney (D<3 linear splice), DF(ISPC 16) on the Siuslaw-family
# forests {2,4,6} splines at 5.0" (op/htdbh.f:228,252), and the D≥100 Sitka-spruce(ISPC 6) OLY/QUIN
# linear override (op/htdbh.f:233). IFOR = op forkod forest index (1..6) = p.forest_idx directly.
# =============================================================================
@inline function op_htdbh_height(ifor::Int, ispc::Int, d::Float32)::Float32
    (ifor < 1 || ifor > 6) && (ifor = 1)
    p2 = PN_HTDBH_P2[ifor, ispc]; p3 = PN_HTDBH_P3[ifor, ispc]; p4 = PN_HTDBH_P4[ifor, ispc]
    if (ifor == 2 || ifor == 4 || ifor == 6) && ispc == 16          # DF splines at 5.0" (op/htdbh.f:252)
        if d >= 5.0f0
            return 4.5f0 + p2 * exp(-1f0 * p3 * d^p4)
        else
            return ((4.5f0 + p2*exp(-1f0*p3*(5.0f0^p4)) - 4.51f0)*(d - 0.3f0)/4.7f0) + 4.51f0
        end
    end
    if d >= 3.0f0
        h = 4.5f0 + p2 * exp(-1f0 * p3 * d^p4)
        (d >= 100f0 && ispc == 6 && (ifor == 1 || ifor == 3)) && (h = 0.25f0*d + 248f0)  # Sitka OLY
        return h
    else
        return ((4.5f0 + p2*exp(-1f0*p3*(3.0f0^p4)) - 4.51f0)*(d - 0.3f0)/2.7f0) + 4.51f0
    end
end

# op/dgf.f MAPSPC(39): FVS species 1..39 → coefficient group JSPC (1..20).
const OP_DG_MAPSPC = Int[
    1, 2, 2, 3, 4,18, 4,15,11,11,16, 6, 5, 5, 6, 7,20, 8, 9,10,
   12,13,14,14,14,14,14,19,14,11,11,11,11,14,14,14,14,14,14]

# op/dgf.f per-group coefficient vectors (JSPC 1..20).
const OP_DGLD = Float32[
    0.919402,0.905119,0.993986,0.904253,0.844690,0.738750,0.802905,0.744005,0.641956,0.857131,
    0.879338,1.024186,0.511442,0.889596,0.816880,0.478504,0.949631,1.049845,1.66609,0.0]
const OP_DGCR = Float32[
    1.290568,1.754811,1.522401,4.123101,1.597250,3.454857,1.936912,0.771395,1.471926,1.505513,
    1.970052,0.459387,0.623093,1.732535,2.471226,1.905011,1.826879,1.632468,0.0,0.0]
const OP_DGCRSQ = Float32[
    0.125823,0.0,0.0,-2.689340,0.0,-1.773805,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const OP_DGSITE = Float32[
    0.541881,0.318254,0.349888,0.684939,0.404010,1.011504,0.495162,0.708166,0.634098,0.208040,
    0.252853,1.965888,0.237269,0.227307,0.244694,0.391327,0.375175,0.0,0.14995,0.0]
const OP_DGDBAL = Float32[
    -0.002133,-0.005355,-0.002979,-0.006368,-0.003726,-0.013091,-0.001827,-0.016240,-0.012589,-0.004101,
    -0.004215,-0.010222,-0.027074,-0.001265,-0.005950,-0.004706,-0.005350,-0.000086,0.0,0.0]
const OP_DGLBA = Float32[
    -0.136818,0.0,0.0,0.0,0.0,-0.131185,-0.129474,-0.130036,-0.085525,0.0, 0.0,0.0,-0.481983,0.0,0.0,0.0,0.0,-0.198636,0.0,0.0]
const OP_DGBA = Float32[
    0.0,0.0,-0.000137,0.0,0.0,0.0,0.0,0.0,0.0,0.0, -0.000173,0.0,0.0,-0.000981,-0.000147,-0.000114,0.000040,0.0,-0.00204,0.0]
const OP_DGBAL = Float32[
    0.0,0.0,0.0,0.0,0.0,0.0,-0.001639,0.003883,0.002385,0.0, 0.0,0.0,0.008903,0.0,0.0,0.0,0.0,-0.002319,-0.00326,0.0]
const OP_DGPCCF = Float32[
    0.0,0.0,0.0,-0.000471,-0.000257,-0.000593,0.0,0.0,0.0,-0.000201, 0.0,-0.000757,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const OP_DGHAH = Float32[
    0.0,-0.000661,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
# DGCONS site terms
const OP_DGCASP = Float32[
    -0.217205,0.0,-0.782418,-0.374512,0.0,0.0,0.014165,-0.106936,-0.056608,-0.104495,
     0.0,0.0,0.022254,0.085958,-0.023186,0.207853,-0.935870,-0.221095,0.0,0.0]
const OP_DGSASP = Float32[
    0.096326,0.0,0.022160,-0.207659,0.0,0.0,0.003263,-0.106020,0.061254,-0.126130,
    0.0,0.0,-0.085538,-0.863980,0.679903,0.378860,0.202507,0.100081,0.0,0.0]
const OP_DGSLOP = Float32[
    -0.265612,0.0,0.319956,0.400223,0.0,0.0,-0.340401,-0.303490,0.736143,0.411602,
    0.0,0.0,0.0,0.0,0.0,-0.066440,0.0,-0.169141,0.0,0.0]
const OP_DGSLSQ = Float32[
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.082191,0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const OP_DGEL = Float32[
    -0.023858,-0.003051,-0.003773,-0.069045,-0.023376,-0.003784,-0.009845,-0.009564,-0.018444,-0.003809,
     0.0,-0.012111,0.0,-0.075986,0.0,-0.005414,0.323546,0.007009,0.0,0.0]
const OP_DGEL2 = Float32[
    0.0,0.0,0.0,0.000608,0.0,0.0000666,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.001193,0.0,0.0,-0.003130,0.0,0.0,0.0]

# op/dgf.f DGFOR(3,20): location-class constants [JSPC][locclass 1..3] (Fortran col-major → rows here).
const OP_DGFOR = Float32[
    -0.627531  0.0       0.0;
    -0.643920  0.0       0.0;
    -1.888949 -1.276180  0.0;
    -1.401865 -1.127977  0.0;
    -0.589570 -0.909553  0.0;
    -2.922255  0.0       0.0;
    -0.739354 -0.199200  0.0;
    -0.688250 -0.405590  0.0;
    -0.594460 -0.522658  0.0;
    -1.052161 -0.793945  0.0;
    -1.310067 -1.432659  0.0;
    -7.753469 -8.279266  0.0;
     4.253807  3.913250  3.507520;
    -0.107648 -0.098335  0.0;
    -1.277664 -1.178041  0.0;
    -0.524624 -0.803095  0.0;
    -9.211184 -9.800653  0.0;
     2.075598  2.100904  0.0;
    -1.33299   0.0       0.0;
     0.0       0.0       0.0]
# op/dgf.f DGDS(3,20): DBH² coefficients [JSPC][dsqclass 1..3].
const OP_DGDS = Float32[
    -0.0002641  0.0        0.0;
    -0.0003137  0.0        0.0;
    -0.0002621  0.0        0.0;
    -0.0003996  0.0        0.0;
    -0.0000596  0.0        0.0;
    -0.0004708  0.0        0.0;
    -0.0000896 -0.0000641  0.0;
    -0.0000572 -0.0000862  0.0;
    -0.0001736 -0.0001040  0.0;
    -0.0002214  0.0        0.0;
    -0.0001323  0.0        0.0;
    -0.0001737  0.0        0.0;
    -0.0005099  0.0        0.0;
     0.0        0.0        0.0;
    -0.0002536  0.0        0.0;
     0.0        0.0        0.0;
    -0.0003552  0.0        0.0;
    -0.0002123 -0.0001361  0.0;
    -0.00154    0.0        0.0;
     0.0        0.0        0.0]
# op/dgf.f MAPLOC(6,20) and MAPDSQ(6,20): forest (1..6) → location/dbh² class per JSPC.
const OP_DG_MAPLOC = Int[
    1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;  1 2 1 2 2 2;  1 2 1 2 2 2;
    1 1 1 1 1 1;  1 2 2 2 2 2;  1 2 1 2 2 2;  1 2 1 2 2 2;  1 1 1 1 1 1;
    1 2 1 2 2 2;  1 2 1 2 2 2;  1 2 3 3 3 3;  1 2 1 2 2 2;  1 2 1 2 2 2;
    1 2 1 2 2 2;  1 1 1 1 1 1;  1 2 1 2 2 2;  1 1 1 1 1 1;  1 1 1 1 1 1]
const OP_DG_MAPDSQ = Int[
    1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;
    1 1 1 1 1 1;  1 2 1 2 2 2;  1 2 1 2 2 2;  1 2 1 2 2 2;  1 1 1 1 1 1;
    1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 2 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;
    1 1 1 1 1 1;  1 1 1 1 1 1;  1 2 1 2 2 2;  1 1 1 1 1 1;  1 1 1 1 1 1]

"""
    op_dg_dgdsq(ispc, ifor) -> Float32

op/dgf.f DGCONS — the DBH² coefficient DGDSQ(JSPC)=DGDS(MAPDSQ(IFOR,JSPC),JSPC) for FVS species
`ispc`. RW(17) forces DGDSQ=0.
"""
@inline function op_dg_dgdsq(ispc::Int, ifor::Int)
    ispc == 17 && return 0.0f0
    jspc = OP_DG_MAPSPC[ispc]
    (ifor < 1 || ifor > 6) && (ifor = 1)
    ispdsq = OP_DG_MAPDSQ[jspc, ifor]
    return OP_DGDS[jspc, ispdsq]
end

"""
    op_dgcon(ispc, ifor, elev, slope, aspect, sitear) -> Float32

op/dgf.f ENTRY DGCONS — the site-dependent diameter-growth constant DGCON(ISPC) for FVS species
`ispc` (1..39), stand `elev`/`slope`/`aspect` and this species' site index `sitear=SITEAR(ISPC)`.
RW(17) uses the redwood log-site form. `ifor` is the FVS forest index (1..6).
"""
@inline function op_dgcon(ispc::Int, ifor::Int, elev::Float32, slope::Float32,
                          aspect::Float32, sitear::Float32)
    if ispc == 17                                            # REDWOOD
        return -3.502444f0 + 0.415435f0*log(sitear)
    end
    jspc = OP_DG_MAPSPC[ispc]
    (ifor < 1 || ifor > 6) && (ifor = 1)
    ispfor = OP_DG_MAPLOC[jspc, ifor]
    sasp = (OP_DGSASP[jspc]*sin(aspect) + OP_DGCASP[jspc]*cos(aspect) + OP_DGSLOP[jspc])*slope +
           OP_DGSLSQ[jspc]*slope*slope
    xsite = sitear
    jspc == 10 && (xsite = xsite*3.281f0)
    jspc == 19 && (xsite = -37.60812f0*log(1f0 - (xsite/114.24569f0)^0.4444f0))
    temel = elev
    (jspc == 14 && temel > 30f0) && (temel = 30f0)
    return OP_DGFOR[jspc, ispfor] + OP_DGEL[jspc]*temel + OP_DGEL2[jspc]*temel*temel +
           OP_DGSITE[jspc]*log(xsite) + sasp
end

"""
    op_dgf_dds(ispc, conspp, dgdsq, d, cr, bal, pccf, relht, ba) -> Float32

op/dgf.f main body — the per-tree LN(DDS) for a FVS-native (IORG=0) tree, DEFAULT branch.
`conspp` = DGCON(ISPC)+COR(ISPC), `dgdsq` = op_dg_dgdsq(ispc,ifor), `cr` = crown ratio 0..1,
`relht` = min(HT/AVH,1.5). Returns LN(DDS) after the 10→5-yr TDDS/2 rescale and the −9.21 clamp.
RA(JSPC=13) and RW(ISPC=17) use their own forms (op_dgf_dds_ra / handled in dgf!).
"""
@inline function op_dgf_dds(ispc::Int, conspp::Float32, dgdsq::Float32, d::Float32,
                            cr::Float32, bal::Float32, pccf::Float32, relht::Float32, ba::Float32)
    jspc = OP_DG_MAPSPC[ispc]
    dds = conspp + OP_DGLD[jspc]*log(d) + cr*(OP_DGCR[jspc] + cr*OP_DGCRSQ[jspc]) +
          dgdsq*d*d + OP_DGDBAL[jspc]*bal/log(d + 1f0)
    dds = dds + OP_DGPCCF[jspc]*pccf + OP_DGHAH[jspc]*relht + OP_DGLBA[jspc]*log(ba) +
          OP_DGBAL[jspc]*bal + OP_DGBA[jspc]*ba
    tdds = exp(dds)
    dds = log(tdds/2f0)
    dds < -9.21f0 && (dds = -9.21f0)
    return dds
end

"""
    op_dgcons!(s)

op/dgf.f ENTRY DGCONS over the stand — fill `s.calib.dg_const[ispc]` = DGCON for every FVS-native
species (op_dgcon). ORGANON species (GF/DF/…) are grown by the NWO engine; their dg_const is inert.
"""
function op_dgcons!(s::StandState)
    c = s.calib; p = s.plot
    ifor = Int(p.forest_idx); (ifor < 1 || ifor > 6) && (ifor = 1)
    elev = p.elevation; slope = p.slope; asp = p.aspect
    @inbounds for ispc in 1:39
        si = ispc <= length(p.sp_site_index) ? p.sp_site_index[ispc] : 0f0
        si <= 0f0 && (si = 1f0)
        c.dg_const[ispc] = op_dgcon(ispc, ifor, elev, slope, asp, si)
    end
    return s
end

"""
    dgf!(s, ::Olympic)

op/dgf.f per-tree WK2 = LN(DDS) for the FVS-native (IORG=0) trees. ORGANON trees (IORG=1) get their
DDS from the NWO engine (organon_nwo.jl) and are skipped here. RA/RW special forms handled inline.
Requires `op_dgcons!(s)` to have filled `s.calib.dg_const`.
"""
function dgf!(s::StandState, ::Olympic)
    p, t, c = s.plot, s.trees, s.calib
    dens = s.density
    wk2 = view(s.scratch.wk, 2, :)
    ba = p.basal_area; avh = p.avg_height
    ifor = Int(p.forest_idx); (ifor < 1 || ifor > 6) && (ifor = 1)
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        ispc = Int(t.species[i])
        # op/dgdriv.f IORG gate — ORGANON-grown trees get their DDS from the NWO engine, skip here.
        (ispc in OP_ORGANON_VALID && t.height[i] > 4.5f0 && d >= 0.1f0) && continue
        pt_i = Int(t.plot_id[i])
        pccf = (1 <= pt_i <= length(dens.point_ccf)) ? dens.point_ccf[pt_i] : 0f0
        cr = Float32(t.crown_pct[i]) * 0.01f0
        pctfrac = 1f0 - t.crown_ratio[i] / 100f0              # 1 − PCT/100
        bal = pctfrac * ba
        relht = avh > 0f0 ? min(t.height[i]/avh, 1.5f0) : 0f0
        cor = c.dg_cor[ispc]
        conspp = ispc == 17 ? c.dg_const[ispc] : c.dg_const[ispc] + cor
        dgdsq = op_dg_dgdsq(ispc, ifor)
        wk2[i] = op_dgf_dds(ispc, conspp, dgdsq, d, cr, bal, pccf, relht, ba)
    end
    return s
end
