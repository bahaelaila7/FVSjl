# =============================================================================
# organon_dgf.jl — OC (Oregon Coast) FVS-native Wykoff diameter growth for the NON-ORGANON
# (IORG=0 surrogate / no-big-6) trees (chunk C9).
#
# Ported from oc/dgf.f (CA-family Wykoff DDS). Every FVS species maps (MAPSPC) to one of 13
# diameter-growth equation groups; DGCONS builds the per-species site constant DGCON once, DGF
# evaluates the per-tree ln(DDS). This is the path oc/dgdriv.f runs for EVERY tree; the ORGANON
# EXECUTE then OVERWRITES WK2 for the valid ORGANON trees (IORG=1). So the IORG=0 trees keep this
# FVS-native DDS. The GS/RW (sp 23/50) alternate functional form is not exercised by ocmin and is a
# follow-up; ocmin's non-ORGANON trees (LP, BR, sub-4.5-ft DF/GF) all use the standard branch.
#
# forkod: KODFOR (STDINFO field-1) → IFOR via the OC JFOR table (oc/forkod.f:61); ocmin 711 = JFOR[9]
# ⇒ IFOR=9 (Medford). Load-bearing for DGCON's DGFOR location class (only DF's group-4 MAPLOC varies).
#
# MEASURED bit-exact vs FVSoc_clean DEBUG-DGF (ocmin): DGCON DF=1.14167/GF=0.07280/LP=0.44505/
# BR=−0.02756; per-tree growth-cycle ln(DDS) for the 10 IORG=0 trees. See docs/OC_VARIANT_PORT_AUDIT.md (C9).
# =============================================================================

# oc/forkod.f:61 JFOR (KODFOR → IFOR index). 518→IFOR 11 remaps to 5 (Trinity→Shasta-Trinity).
const OC_JFOR = Int[505,506,508,511,514,610,611,710,711,712,518]
"oc/forkod.f — decode KODFOR (location code) → IFOR (1..10). 0/unknown ⇒ 9 (Medford default region)."
function oc_forkod(kodfor::Integer)
    for (i, jf) in enumerate(OC_JFOR)
        kodfor == jf && (return i == 11 ? 5 : i)          # 518→11→5 mapping (forkod.f:212-218)
    end
    return 9
end

# oc/dgf.f MAPSPC(50): FVS species index (1..50) → one of 13 DG equation groups.
const OC_DGF_MAPSPC = Int[
    1,1,1,2,3,3,4,7,7,6,5,6,5,5,9,7,8,9,9,5,
    5,2,12,5,9,10,10,10,10,10,10,10,10,10,10,10,11,11,11,10,
    10,13,10,10,10,10,5,10,10,12]

# oc/dgf.f 13-group coefficient arrays.
const OC_DGF_DGLD   = Float32[0.950418,1.182104,1.186676,0.716226,1.077154,1.218279,0.886150,0.825682,0.738750,1.310111,0.955569,0.0,0.99531]
const OC_DGF_DGCR   = Float32[1.815305,2.856578,2.763519,3.272451,-0.276387,3.167164,1.478650,1.675208,3.454857,0.271183,0.0,0.0,2.08524]
const OC_DGF_DGCRSQ = Float32[0.0,-1.093354,-0.871061,-1.642904,1.063732,-1.568333,0.0,0.0,-1.773805,0.0,0.0,0.0,-0.98396]
const OC_DGF_DGSITE = Float32[0.820451,0.365679,0.492695,0.759305,0.0,0.566946,0.963375,0.724300,1.011504,0.213526,1.334008,0.0,0.00659]
const OC_DGF_DGDBAL = Float32[-0.005433,-0.005992,-0.003728,-0.008787,0.0,0.0,-0.006263,-0.002133,-0.013091,0.0,-0.005893,0.0,-0.00147]
const OC_DGF_DGLBA  = Float32[-0.000016,-0.058039,-0.122905,-0.028564,0.0,-0.267873,-0.129146,-0.203636,-0.131185,0.0,-0.408462,0.0,0.0]
const OC_DGF_DGBAL  = Float32[0.0,0.0,0.0,0.0,-0.000893,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const OC_DGF_DGPCCF = Float32[-0.000779,-0.001014,0.0,-0.000224,0.0,-0.000338,0.0,0.0,-0.000593,-0.000473,0.0,0.0,-0.00018]
const OC_DGF_DGHAH  = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.50155]
const OC_DGF_DGDS   = Float32[-0.0002385,-0.0006362,-0.0004572,-0.0002723,0.0,-0.0014178,-0.0002528,-0.0000731,-0.0004708,-0.0003048,0.0,0.0,-0.000373]
const OC_DGF_DGCASP = Float32[0.0,-0.315227,-0.444594,-0.151727,0.649870,0.0,-0.280294,-0.179510,0.0,0.0,0.0,0.0,-0.19935]
const OC_DGF_DGSASP = Float32[0.0,0.097350,0.139180,0.018681,0.951834,0.0,-0.014463,-0.562259,0.0,0.0,0.0,0.0,-0.03587]
const OC_DGF_DGSLOP = Float32[0.0,-0.206267,0.0,-0.339369,0.0,0.0,-0.581722,-0.544867,0.0,0.0,0.0,0.0,0.73530]
const OC_DGF_DGSLSQ = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.99561]
const OC_DGF_DGEL   = Float32[0.0,0.0301,0.0248,-0.0141,0.0,0.0,0.0,0.0,-0.003784,0.0049,0.0,0.0,0.0]
const OC_DGF_DGELSQ = Float32[0.0,-0.00030732,-0.00033429,0.00024083,0.0,0.0,0.0,0.0,0.00006660,-0.00008781,0.0,0.0,0.0]
const OC_DGF_OBSERV = Float32[613.,3759.,2062.,5400.,84.,372.,561.,253.,2482.,306.,336.,8928.,6504.]
# DGFOR(5,13): [location class 1..5, group 1..13].
const OC_DGF_DGFOR = Float32[
    -3.428338 -2.108357 -2.073942 -1.877695  0.564402 -2.058828 -2.397678 -1.626879 -2.922255 -1.958189 -3.344700  0.0     -0.94563;
    -3.966547  0.0       -1.943608 -2.099646  0.0       -1.596998  0.0        0.0       0.0       0.0       0.0       0.0      0.0;
     0.0       0.0        0.0      -2.211587  0.0        0.0        0.0        0.0       0.0       0.0       0.0       0.0      0.0;
     0.0       0.0        0.0      -1.955301  0.0        0.0        0.0        0.0       0.0       0.0       0.0       0.0      0.0;
     0.0       0.0        0.0      -2.078432  0.0        0.0        0.0        0.0       0.0       0.0       0.0       0.0      0.0]
# MAPLOC(10,13): [IFOR 1..10, group 1..13].
const OC_DGF_MAPLOC = Int[
    1 1 1 1 1 1 1 1 1 1 1 1 1;
    1 1 1 2 1 1 1 1 1 1 1 1 1;
    1 1 1 2 1 1 1 1 1 1 1 1 1;
    1 1 2 1 1 2 1 1 1 1 1 1 1;
    2 1 1 3 1 1 1 1 1 1 1 1 1;
    1 1 1 4 1 1 1 1 1 1 1 1 1;
    1 1 1 5 1 1 1 1 1 1 1 1 1;
    1 1 1 4 1 1 1 1 1 1 1 1 1;
    1 1 1 4 1 1 1 1 1 1 1 1 1;
    1 1 1 5 1 1 1 1 1 1 1 1 1]

"""
    oc_dgcons!(s)

oc/dgf.f DGCONS — the per-species site-dependent DG constant `DGCON`, `DGDSQ`, and `ATTEN`,
resolved once. `DGCON = DGFOR[ispfor,jspc] + DGEL·ELEV + DGELSQ·ELEV² + DGSITE·ln(SITEAR) +
(DGSASP·sin+DGCASP·cos+DGSLOP)·SLOPE + DGSLSQ·SLOPE²`; `ispfor = MAPLOC[IFOR,jspc]`. GS/RW use a
distinct constant (not exercised by ocmin — deferred).
"""
function oc_dgcons!(s::StandState)
    p = s.plot; c = s.calib
    ifor = Int(p.forest_idx); (ifor < 1 || ifor > 10) && (ifor = 9)
    elev = p.elevation; slope = p.slope; aspect = p.aspect
    sasp_sin = sin(aspect); sasp_cos = cos(aspect)
    @inbounds for sp in 1:50
        jspc = OC_DGF_MAPSPC[sp]
        if sp == 23 || sp == 50   # GS/RW — alt form (deferred; ocmin has none)
            c.dg_const[sp] = -3.502444f0 + 0.415435f0*flog(p.sp_site_index[sp])
            c.dg_dsq[sp] = 0f0
        else
            ispfor = OC_DGF_MAPLOC[ifor, jspc]; (ispfor < 1 || ispfor > 5) && (ispfor = 1)
            sasp = (OC_DGF_DGSASP[jspc]*sasp_sin + OC_DGF_DGCASP[jspc]*sasp_cos + OC_DGF_DGSLOP[jspc])*slope +
                   OC_DGF_DGSLSQ[jspc]*slope*slope
            c.dg_const[sp] = OC_DGF_DGFOR[ispfor, jspc] + OC_DGF_DGEL[jspc]*elev +
                             OC_DGF_DGELSQ[jspc]*elev*elev + OC_DGF_DGSITE[jspc]*flog(p.sp_site_index[sp]) + sasp
            c.dg_dsq[sp] = OC_DGF_DGDS[jspc]
        end
        c.atten[sp] = OC_DGF_OBSERV[jspc]
        c.bark_a[sp] = 0f0; c.bark_b[sp] = 0f0   # OC uses oc_bratio directly (not the linear a+b·d)
    end
    return s
end

"""
    dgf!(s, ::OregonCoast)

oc/dgf.f main body — per-tree WK2 = ln(DDS) for EVERY tree (the ORGANON trees' WK2 is overwritten
later). Standard Wykoff branch: `DDS = CONSPP + DGLD·lnD + CR·(DGCR+CR·DGCRSQ) + DGDSQ·D² +
DGDBAL·BAL/ln(D+1) + DGPCCF·PCCF + DGHAH·relHt + DGLBA·lnBA + DGBAL·BAL`, CONSPP=DGCON+COR,
BAL=(1−PCT/100)·BA, relHt=min(HT/AVH,1.5). Then the 10→5-yr halving `ln(exp(DDS)/2)`, floored −9.21.
Tanoak (sp 42) is a 5-yr eqn → ×2 first. (GS/RW deferred.)
"""
function dgf!(s::StandState, ::OregonCoast)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    wk2 = view(s.scratch.wk, 2, :)
    ba = p.basal_area; avh = p.avg_height
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        sp = Int(t.species[i]); jspc = OC_DGF_MAPSPC[sp]
        conspp = c.dg_const[sp] + c.dg_cor[sp]
        cr = Float32(t.crown_pct[i])*0.01f0
        bal = (1f0 - t.crown_ratio[i]/100f0)*ba
        pt = Int(t.plot_id[i])
        pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0f0
        relht = avh > 0f0 ? t.height[i]/avh : 0f0
        relht > 1.5f0 && (relht = 1.5f0)
        dds = conspp + OC_DGF_DGLD[jspc]*flog(d) +
              cr*(OC_DGF_DGCR[jspc] + cr*OC_DGF_DGCRSQ[jspc]) +
              c.dg_dsq[sp]*d*d + OC_DGF_DGDBAL[jspc]*bal/flog(d + 1f0)
        dds = dds + OC_DGF_DGPCCF[jspc]*pccf + OC_DGF_DGHAH[jspc]*relht +
              OC_DGF_DGLBA[jspc]*flog(ba) + OC_DGF_DGBAL[jspc]*bal
        sp == 42 && (dds = flog(fexp(dds)*2f0))                  # tanoak 5→10yr
        dds = flog(fexp(dds)/2f0)                                # 10→5yr
        dds < -9.21f0 && (dds = -9.21f0)
        wk2[i] = dds
    end
    return s
end

# --- R5CRWD (bin/FVSoc_buildDir/r5crwd.f) — CA-family crown width, for the point CCF (PCCF) ------
# MAPCA(50): FVS species index → R5 crown-width equation index (1..35).
const OC_R5CRWD_MAPCA = Int[
    27, 7,16, 4, 5, 5, 1,14, 6,24,
    27,23,32,27,22, 2, 3,21,27,33,
    25,13,27,15,27,34, 9,26,28,19,
     8,35,28,28,28,17,11,20,28,18,
    28,10,28,18,18,28,28,28,28, 1]
const OC_R5CRWD_WB1 = Float32[6.81,-1.476,-0.997,5.82,6.71,4.72,7.11,10.0,5.0,10.0,1.0,6.19,6.50,4.57,4.2,4.00,8.00,0.50,3.08,2.98,2.24,1.52,1.91,2.37,4.31,4.49,6.0,2.0,27.030,12.733,9.0684,3.9347,3.8273,5.3732,4.5628]
const OC_R5CRWD_WB2 = Float32[0.732,1.01,0.92,0.591,0.421,0.608,0.470,1.20,1.69,1.05,1.43,1.01,1.80,1.41,1.42,1.60,1.53,1.62,1.92,1.55,0.763,0.891,0.784,0.736,0.628,0.688,0.6,1.5,0.8612,2.249,0.4702,0.7086,0.7624,0.7707,0.6925]
const OC_R5CRWD_WB3 = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.014,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const OC_R5CRWD_DX1 = Float32[3.62,3.5,3.5,3.26,3.5,3.5,3.5,2.5,2.5,2.23,3.11,3.5,3.5,3.5,3.5,3.5,2.5,2.5,2.5,2.15,3.77,3.5,3.5,3.5,3.5,2.5,3.5,2.5,3.5,2.5,3.5,3.5,3.5,2.5,2.5]
const OC_R5CRWD_DX2 = Float32[1.370,0.338,0.329,1.103,1.063,0.852,1.192,2.700,2.190,1.630,1.008,1.548,2.400,1.624,1.560,1.700,2.630,1.220,2.036,1.646,0.7756,0.5754,0.6492,0.8496,1.6684,2.2175,1.1,1.4,5.5672,4.2956,3.1654,1.7618,1.9108,3.2150,2.2816]
const OC_R5CRWD_IEQN = Int[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3,2,2,2,2,2,2,1,1,1,1,2,2,2,2,2]
const OC_R5CRWD_SPLINE = Float32[5.0,7.4,7.6,5.0,5.0,5.0,5.0,5.0,5.0,13.4,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0]
const OC_R5CRWD_SM = Float32[0.7778,0.7778,0.7778,0.7778,0.7778,0.7778,0.7778,0.5556,0.5556,0.5556,0.5556,0.7778,0.7778,0.7778,0.7778,0.7778,0.5556,0.5556,0.5556,0.5556,0.7778,0.7778,0.7778,0.7778,0.7778,0.5556,0.7778,0.5556,0.7778,0.5556,0.7778,0.7778,0.7778,0.5556,0.5556]

"bin/FVSoc_buildDir/r5crwd.f R5CRWD — CA-family crown width (ft) for FVS species `sp`, DBH `d`, HT `h`."
@inline function oc_r5crwd(sp::Int, d::Float32, h::Float32)
    idx = OC_R5CRWD_MAPCA[sp]
    ity = OC_R5CRWD_IEQN[idx]; spdbh = OC_R5CRWD_SPLINE[idx]
    if d >= spdbh
        if ity == 1
            return OC_R5CRWD_WB1[idx] + OC_R5CRWD_WB2[idx]*d
        elseif ity == 2
            return OC_R5CRWD_WB1[idx]*fpow(d, OC_R5CRWD_WB2[idx])
        else
            return OC_R5CRWD_WB1[idx] + OC_R5CRWD_WB2[idx]*d + OC_R5CRWD_WB3[idx]*d*d
        end
    elseif h >= 4.5f0
        return OC_R5CRWD_DX1[idx] + OC_R5CRWD_DX2[idx]*d
    else
        return OC_R5CRWD_SM[idx]*h
    end
end

"oc/ccfcal.f MODE=1 — per-tree CCF (÷TPA): `0.001803·CRWD5²` (CRWD5 = R5 crown width)."
@inline oc_tree_ccf(sp::Int, d::Float32, h::Float32) = 0.001803f0*(oc_r5crwd(sp, d, h)^2)
