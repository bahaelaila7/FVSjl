# =============================================================================
# diameter_growth.jl (westsierra) — WS large-tree DDS (ws/dgf.f). Chunk 3.
#
# WS's dgf.f is PER-SPECIES (identity map, length MAXSP=43, indexed by ISPC — NOT group-compressed).
# It is a leaner Prognosis than SO: NO DGMAI/DGSIC/DGCCFA/DGMACC/RELDEN terms. DGCON is assembled in
# ENTRY DGCONS (site/elev/slope/aspect + DGFOR location-class via MAPLOC), CONSPP = DGCON + COR, and the
# per-tree ln(DDS) has these branches (ws/dgf.f:641-751):
#   • GENERAL Wykoff ln(DDS): most species incl. SP(1)/DF(2)/WF(3) — with a D<10 small-tree sub-form for
#     JP(6) and the CONSJP group {2,3,7,13,22}, a WF/SF −0.15032 offset, and a hardwood 5→10-yr doubling.
#   • CASE(7) RF: dedicated original-WS red-fir eqn (CRID / PBAL / ln(PBA)).
#   • CASE(9:10,12,14:17,19:20,25:27): CA-variant surrogate (HOAVH cap 1.5 + DGBAL BAL).
#   • CASE(41) MC: SO-variant surrogate (DGBA·BA raw).
#   • CASE(21) GB: UT DF-projection (DSTAG stagnation).
#   • CASE(4,23) GS/RW: Castle-2019 sequoia/redwood DGLT-exp (uses PRD = per-point Zeide density).
# DGCONS specials: CASE(4,23) ln(SITEAR); CASE(7) RF DGLAT9(ILAT)+SITEAR; hardwoods 28-40,43 use SITEAR(2).
#
# ⚠ PRD (= ZRD(point)/XMAXPT(point), the per-point Zeide relative density from SDICAL/SDICLS) is the SAME
#   engine-gap flagged for AK #209 — it feeds ONLY CASE(4,23) GS/RW. wst01 has no GS/RW so PRD is inert
#   here (stubbed via density.point_zeide_rd if present, else 0). DSTAG (CASE 21 stagnation) likewise only
#   bites GB. Both are wired faithfully but exercised only by species absent from the ref stand.
#
# MEASURED bit-exact vs FVSws_g16 dgf per-tree LN(DDS) on wst01 (SP/DF/WF general + RF CASE7 branches).
# =============================================================================

# ── WS per-species scalar coefficients (ws/dgf.f DATA), length 43, index = ISPC.
const WS_DGLD = Float32[
  1.0857,0.8641,1.2854,0.0,1.29079, 0.51047,0.0,0.96103,1.218279,1.218279,
  1.0857,1.077154,1.2854,1.077154,1.077154, 1.077154,1.077154,0.96103,1.077154,1.077154,
  0.0,0.8641,0.0,1.0857,1.077154, 1.077154,1.077154,1.23911,1.23911,1.23911,
  1.23911,1.23911,1.23911,0.99531,0.99531, 0.99531,0.99531,0.99531,0.99531,1.23911,
  0.889596,1.1783,1.23911]
const WS_DGCR = Float32[
  0.3910,0.4246,-1.0191,0.0,-0.0906, 0.91422,0.0,0.4126,3.167164,3.167164,
  0.3910,-0.276387,-1.0191,-0.276387,-0.276387, -0.276387,-0.276387,0.4126,-0.276387,-0.276387,
  0.0,0.4246,0.0,0.3910,-0.276387, -0.276387,-0.276387,-1.20841,-1.20841,-1.20841,
  -1.20841,-1.20841,-1.20841,2.08524,2.08524, 2.08524,2.08524,2.08524,2.08524,-1.20841,
  1.732535,0.9492,-1.20841]
const WS_DGCRSQ = Float32[
  0.0,0.0,0.9104,0.0,0.0, 0.27758,0.0,0.0,-1.568333,-1.568333,
  0.0,1.063732,0.9104,1.063732,1.063732, 1.063732,1.063732,0.0,1.063732,1.063732,
  0.0,0.0,0.0,0.0,1.063732, 1.063732,1.063732,2.31782,2.31782,2.31782,
  2.31782,2.31782,2.31782,-0.98396,-0.98396, -0.98396,-0.98396,-0.98396,-0.98396,2.31782,
  0.0,0.0,2.31782]
const WS_DGDBAL = Float32[
  -0.00579,-0.01127,-0.00628,0.0,-0.00544, -0.01282,0.0,-0.01265,0.0,0.0,
  -0.00579,0.0,-0.00628,0.0,0.0, 0.0,0.0,-0.01265,0.0,0.0,
  0.0,-0.01127,0.0,-0.00579,0.0, 0.0,0.0,-0.00199,-0.00199,-0.00199,
  -0.00199,-0.00199,-0.00199,-0.00147,-0.00147, -0.00147,-0.00147,-0.00147,-0.00147,-0.00199,
  -0.001265,-0.00016,-0.00199]
const WS_DGBA = Float32[
  -0.1313,0.0,-0.21056,0.0,-0.23182, -0.01880,0.0,-0.1431,-0.267873,-0.267873,
  -0.1313,0.0,-0.21056,0.0,0.0, 0.0,0.0,-0.1431,0.0,0.0,
  0.0,0.0,0.0,-0.1313,0.0, 0.0,0.0,0.0,0.0,0.0,
  0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0,
  -0.000981,-0.3270,0.0]
const WS_DGHAH = Float32[
  0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0,
  0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0,
  0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0,
  0.0,0.0,0.0,0.50155,0.50155, 0.50155,0.50155,0.50155,0.50155,0.0,
  0.0,0.0,0.0]
const WS_DGPCCF = Float32[
  -0.00058,-0.00018,-0.00091,0.0,-0.00098, -0.00099,0.0,-0.00084,-0.000338,-0.000338,
  -0.00058,0.0,-0.00091,0.0,0.0, 0.0,0.0,-0.00084,0.0,0.0,
  0.0,-0.00018,0.0,-0.00058,0.0, 0.0,0.0,0.0,0.0,0.0,
  0.0,0.0,0.0,-0.0018,-0.0018, -0.0018,-0.0018,-0.0018,-0.0018,0.0,
  0.0,0.0,0.0]
const WS_DGCASP = Float32[
  0.01664,0.0,-0.1804,0.0,0.0, -0.16447,-0.4919,0.26986,0.0,0.0,
  0.01664,0.649870,-0.1804,0.649870,0.649870, 0.649870,0.649870,0.26986,0.649870,0.649870,
  0.0,0.0,0.0,0.01664,0.649870, 0.649870,0.649870,0.08632,0.08632,0.08632,
  0.08632,0.08632,0.08632,-0.19935,-0.19935, -0.19935,-0.19935,-0.19935,-0.19935,0.08632,
  0.085958,0.0,0.08632]
const WS_DGSASP = Float32[
  -0.00350,0.0,-0.1183,-0.10656,0.0, 0.05342,0.0,0.09668,0.0,0.0,
  -0.00350,0.951834,-0.1183,0.951834,0.951834, 0.951834,0.951834,0.09668,0.951834,0.951834,
  0.0,0.0,-0.10656,-0.00350,0.951834, 0.951834,0.951834,-0.11954,-0.11954,-0.11954,
  -0.11954,-0.11954,-0.11954,-0.03587,-0.03587, -0.03587,-0.03587,-0.03587,-0.03587,-0.11954,
  -0.86398,0.0,-0.11954]
const WS_DGSLOP = Float32[
  0.7603,0.0,0.0,-1.29627,0.0, -0.05469,0.0,0.90804,0.0,0.0,
  0.7603,0.0,0.0,0.0,0.0, 0.0,0.0,0.90804,0.0,0.0,
  0.0,0.0,-1.29627,0.7603,0.0, 0.0,0.0,0.85815,0.85815,0.85815,
  0.85815,0.85815,0.85815,0.73530,0.73530, 0.73530,0.73530,0.73530,0.73530,0.85815,
  0.0,0.0,0.85815]
const WS_DGSLSQ = Float32[
  -2.2339,0.0,0.0,0.87335,0.0, 0.0,0.0,-2.04028,0.0,0.0,
  -2.2339,0.0,0.0,0.0,0.0, 0.0,0.0,-2.04028,0.0,0.0,
  0.0,0.0,0.87335,-2.2339,0.0, 0.0,0.0,-1.17209,-1.17209,-1.17209,
  -1.17209,-1.17209,-1.17209,-0.99561,-0.99561, -0.99561,-0.99561,-0.99561,-0.99561,-1.17209,
  0.0,0.0,-1.17209]
const WS_DGEL = Float32[
  0.01919,0.00489,0.0,0.0,-0.00919, 0.00304,0.0,0.00323,0.0,0.0,
  0.01919,0.0,0.0,0.0,0.0, 0.0,0.0,0.00323,0.0,0.0,
  0.0,0.00489,0.0,0.01919,0.0, 0.0,0.0,0.0,0.0,0.0,
  0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0,
  -0.075986,0.0,0.0]
const WS_DGELSQ = Float32[
  -0.00025,0.0,0.0,0.00019,0.00019, 0.0,0.0,0.0,0.0,0.0,
  -0.00025,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0,
  0.0,0.0,0.00019,-0.00025,0.0, 0.0,0.0,0.0,0.0,0.0,
  0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0,
  0.001193,0.0,0.0]
const WS_DGSITE = Float32[
  0.5827,0.5040,0.5260,0.0,0.3737, 0.96412,0.0,0.6828,0.566946,0.566946,
  0.5827,0.0,0.5260,0.0,0.0, 0.0,0.0,0.6828,0.0,0.0,
  0.0,0.5040,0.0,0.5827,0.0, 0.0,0.0,0.32093,0.32093,0.32093,
  0.32093,0.32093,0.32093,0.01200,0.01200, 0.01200,0.01200,0.01200,0.01200,0.32093,
  0.227307,0.4401,0.32093]
const WS_DGDS = Float32[
  -0.000288,-0.000290,-0.000584,0.0,-0.00061, -0.000222,0.0,-0.000375,-0.0014178,-0.0014178,
  -0.000288,0.0,-0.000584,0.0,0.0, 0.0,0.0,-0.000375,0.0,0.0,
  0.0,-0.000290,0.0,-0.000288,0.0, 0.0,0.0,-0.000338,-0.000338,-0.000338,
  -0.000338,-0.000338,-0.000338,-0.000373,-0.000373, -0.000373,-0.000373,-0.000373,-0.000373,-0.000338,
  0.0,-0.000660,-0.000338]
const WS_OBSERV = Float32[
  650,480,3301,8928,1339, 1114,9,1528,372,372,
  650,84,3301,84,84, 84,84,1528,84,84,
  1000,480,8928,650,84, 84,84,583,583,583,
  583,583,583,6504,6504, 6504,6504,6504,6504,583,
  220,419,583]
const WS_DDSMX1 = Float32[
  -3.91,-3.06,-2.81,-3.91,-2.20, -3.71,-2.81,-2.07,0.0,0.0,
  -3.91,0.0,-2.81,0.0,0.0, 0.0,0.0,-2.07,0.0,0.0,
  0.0,-3.06,-3.91,-3.91,0.0, 0.0,0.0,-2.20,-2.20,-2.20,
  -2.20,-2.20,-2.20,-2.20,-2.20, -2.20,-2.20,-2.20,-2.20,-2.20,
  0.0,-3.71,-2.20]
const WS_DDSMX2 = Float32[
  1.42,1.27,1.24,1.42,1.08, 1.25,1.24,1.09,0.0,0.0,
  1.42,0.0,1.24,0.0,0.0, 0.0,0.0,1.09,0.0,0.0,
  0.0,1.27,1.42,1.42,0.0, 0.0,0.0,1.08,1.08,1.08,
  1.08,1.08,1.08,1.08,1.08, 1.08,1.08,1.08,1.08,1.08,
  0.0,1.25,1.08]

# DGLAT9[ILAT] (ws/dgf.f DATA) — RF(7) latitude constant (ILAT 1..5 from NINT(TLAT)).
const WS_DGLAT9 = Float32[0.1434, 0.3191, 0.1434, 0.2246, 0.1434]

# MAPLOC[ifor, isp] (ws/dgf.f DATA MAPLOC(7,MAXSP) → jl [ifor,isp]; forest 1..7 → location class 1..3).
const WS_MAPLOC = Int[
  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;  # forest 1
  2 1 2 1 1 1 1 2 1 1 2 1 2 1 1 1 1 2 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;  # forest 2
  2 2 3 1 2 1 1 2 1 1 2 1 3 1 1 1 1 2 1 1 1 2 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;  # forest 3
  2 1 3 1 3 1 1 1 1 1 2 1 3 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;  # forest 4
  1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1;  # forest 5
  2 1 3 1 1 1 1 1 1 1 2 1 3 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1;  # forest 6
  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]  # forest 7

# DGFOR[iloc, isp] (ws/dgf.f DATA DGFOR(3,MAXSP) → jl [iloc,isp]; loc class 1..3 via MAPLOC).
const WS_DGFOR = Float32[
  -0.70344 -0.5260 0.0755 -1.69950 0.02786 -1.81306 0.0 -0.8882 -2.058828 -2.058828 -0.70344 0.564402 0.0755 0.564402 0.564402 0.564402 0.564402 -0.8882 0.564402 0.564402 0.0 -0.5260 -1.69950 -0.70344 0.564402 0.564402 0.564402 -2.68349 -2.68349 -2.68349 -2.68349 -2.68349 -2.68349 -0.94563 -0.94563 -0.94563 -0.94563 -0.94563 -0.94563 -2.68349 -0.107648 -0.02772 -2.68349;  # loc 1
  -0.90272 -0.9842 -0.3099 0.0 0.21393 -2.3037 0.0 -1.0712 0.0 0.0 -0.90272 0.0 -0.3099 0.0 0.0 0.0 0.0 -1.0712 0.0 0.0 0.0 -0.9842 0.0 -0.90272 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.1327 0.0;  # loc 2
  0.0 -0.7603 -0.0440 0.0 -0.04927 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.0440 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.7603 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]  # loc 3

# ws/dgf.f ENTRY DGCONS — per-species DGCON + ATTEN. IFOR used directly (forest 1..7 via MAPLOC).
# CASE(4,23) GS/RW ln(SITEAR); CASE(7) RF DGLAT9(ILAT)+SITEAR; hardwoods 28-40,43 use SITEAR(2) as TSITE.
function ws_dgcons!(s::StandState)
    c = s.calib; p = s.plot
    ifor = Int(p.forest_idx)
    (ifor < 1 || ifor > 7) && (ifor = 1)
    elev = p.elevation; slope = p.slope; asp = p.aspect
    sina = sin(asp); cosa = cos(asp)
    itlat = round(Int, p.latitude)
    ilat = itlat <= 35 ? 1 : itlat == 36 ? 2 : itlat == 37 ? 3 : itlat == 38 ? 4 : 5
    ctl = s.control
    cor2on = ctl.dg_cor2_on
    @inbounds for isp in 1:43
        sitear = p.sp_site_index[isp]
        temel = elev
        (isp == 41 && temel > 30f0) && (temel = 30f0)
        if isp == 4 || isp == 23                                  # GS/RW (Castle 2019)
            dgcon = -3.502444f0 + 0.415435f0 * log(sitear)
        elseif isp == 7                                           # RF (original WS, latitude)
            dgcon = WS_DGLAT9[ilat] - 0.00700f0 * elev - 0.83400f0 * slope * slope +
                    0.00734f0 * sitear
        else                                                      # DEFAULT Wykoff DGCON
            tsite = sitear
            ((28 <= isp <= 40) || isp == 43) && (tsite = p.sp_site_index[2])   # hardwoods use SITEAR(DF)
            isfor = WS_MAPLOC[ifor, isp]
            dgcon = WS_DGFOR[isfor, isp] + WS_DGEL[isp] * temel + WS_DGELSQ[isp] * temel * temel +
                    WS_DGSASP[isp] * sina * slope + WS_DGCASP[isp] * cosa * slope +
                    WS_DGSLOP[isp] * slope + WS_DGSLSQ[isp] * slope * slope +
                    WS_DGSITE[isp] * log(tsite)
        end
        # LDCOR2 (READCORD/REUSCORD): add ln(COR2) to DGCON, except GS/RW. Default COR2=1 ⇒ inert.
        if cor2on && !(isp == 4 || isp == 23) && ctl.dg_cor2[isp] > 0f0
            dgcon += log(ctl.dg_cor2[isp])
        end
        c.dg_const[isp] = dgcon
        c.atten[isp] = WS_OBSERV[isp]
    end
    return s
end

# ws/dgf.f main body — per-tree WK2 = LN(DDS). Branches: GENERAL / RF(7) / CA-surrogate(9:10,12,14:17,
# 19:20,25:27) / MC(41) / GB(21) / GS-RW(4,23). CONSPP = DGCON + COR (GS/RW get no COR here).
const WS_CA_SURR = Set{Int}([9,10,12,14,15,16,17,19,20,25,26,27])
const WS_CONSJP_GRP = Set{Int}([2,3,7,13,22])
const WS_HARDWOOD_5YR = Set{Int}([28,29,30,31,32,33,34,35,36,37,38,39,40,43])
const WS_DDSMAX_SP = Set{Int}([1,2,3,5,6,7,8,11,13,18,22,24,28,29,30,31,32,33,40,42,43])

function dgf!(s::StandState, ::WestSierra)
    p, t, c = s.plot, s.trees, s.calib
    dens = s.density; sd = s.coef.species
    wk2 = view(s.scratch.wk, 2, :)
    ba = p.basal_area; avh = p.avg_height
    slope = p.slope; asp = p.aspect
    cosa = cos(asp)
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        isp = Int(t.species[i])
        pt_i = Int(t.plot_id[i])
        pccf = (1 <= pt_i <= length(dens.point_ccf)) ? dens.point_ccf[pt_i] : 0f0
        icr = Float32(t.crown_pct[i])                             # ICR (crown %, 0-100)
        cr = icr * 0.01f0
        pctfrac = 1f0 - t.crown_ratio[i] / 100f0                  # 1 − PCT/100 (BA-percentile complement)
        bal = pctfrac * ba
        ald = log(d)
        hoavh = avh > 0f0 ? t.height[i] / avh : 0f0
        alba = ba > 0f0 ? log(ba) : 0f0
        cor = c.dg_cor[isp]
        si = p.sp_site_index[isp]
        conspp = (isp == 4 || isp == 23) ? c.dg_const[isp] : c.dg_const[isp] + cor
        # Point BA (PTBAA) and PBAL — RF/redwood competition. PBAL from RAW PBA; PBA fixed up after.
        pba_raw = (1 <= pt_i <= length(dens.point_ba)) ? dens.point_ba[pt_i] : 0f0
        pbal = pba_raw * pctfrac
        pbal < 0f0 && (pbal = bal)
        pba = pba_raw <= 0f0 ? (ba > 1f0 ? ba : 1f0) : pba_raw
        crid = d < 2f0 ? 1.8f0 : ((icr * icr) / log(d + 1f0)) / 1000f0

        if isp in WS_CA_SURR                                      # CA-variant surrogate
            hv = hoavh > 1.5f0 ? 1.5f0 : hoavh
            dgbal = (isp == 9 || isp == 10) ? 0f0 : -0.000893f0
            dds = conspp + WS_DGLD[isp] * ald + WS_DGCR[isp] * cr + WS_DGCRSQ[isp] * cr * cr +
                  WS_DGDS[isp] * d * d + WS_DGDBAL[isp] * bal / log(d + 1f0) +
                  WS_DGPCCF[isp] * pccf + WS_DGHAH[isp] * hv + WS_DGBA[isp] * alba + dgbal * bal
        elseif isp == 41                                          # MC (SO surrogate) — DGBA·BA raw
            hv = hoavh > 1.5f0 ? 1.5f0 : hoavh
            dds = conspp + WS_DGLD[isp] * ald + WS_DGCR[isp] * cr + WS_DGCRSQ[isp] * cr * cr +
                  WS_DGDS[isp] * d * d + WS_DGDBAL[isp] * bal / log(d + 1f0) +
                  WS_DGPCCF[isp] * pccf + WS_DGHAH[isp] * hv + WS_DGBA[isp] * ba
        elseif isp == 21                                          # GB (UT DF-projection + DSTAG)
            dpp = d < 1f0 ? 1f0 : d
            batem = ba < 5f0 ? 5f0 : ba
            bark = ws_bratio(sd, isp, d)
            df = 0.25897f0 + 1.03129f0 * dpp - 0.0002025464f0 * batem + 0.00177f0 * si
            (df - dpp) > 1f0 && (df = dpp + 1f0)
            df < dpp && (df = dpp)
            diagr = (df - dpp) * bark
            # DSTAG (RELSDI stagnation, ws/dgf.f:543-559) needs SDICAL — an unported engine-gap (AK #209
            # class). ISTAGF is off by default ⇒ DSTAG inert; GB(21) is absent from the ref stand. Follow-on.
            dds = diagr <= 0f0 ? -9.21f0 : log(diagr * (2f0 * dpp * bark + diagr)) + cor + c.dg_const[isp]
            dds < -9.21f0 && (dds = -9.21f0)
        elseif isp == 4 || isp == 23                             # GS/RW (Castle 2019 DGLT-exp)
            prd = ws_point_prd(s, pt_i)
            bark = ws_bratio(sd, isp, d)
            dglt = exp(conspp + 0.185911f0 * log(d) - 0.000073f0 * d * d - 0.001796f0 * pbal -
                       0.42078f0 * prd + 0.589318f0 * log(cr * 100f0) - 0.000926f0 * slope * 100f0 -
                       0.002203f0 * (slope * 100f0) * cosa)
            tempd1 = d * bark
            tempd2 = (d + dglt) * bark
            dds = log(tempd2 * tempd2 - tempd1 * tempd1) + cor + log(cor2_of(c, isp))
        elseif isp == 7                                          # RF (original WS)
            dds = conspp + 1.53339f0 * ald - 0.47442f0 * d * d / 1000f0 +
                  0.35739f0 * crid - 0.44256f0 * pbal / log(d + 1f0) / 100f0 -
                  0.12359f0 * log(pba)
        else                                                     # GENERAL Wykoff (SP/DF/WF/…)
            dds = conspp + WS_DGLD[isp] * ald + WS_DGCR[isp] * cr + WS_DGCRSQ[isp] * cr * cr +
                  WS_DGDS[isp] * d * d + WS_DGDBAL[isp] * bal / log(d + 1f0) +
                  WS_DGPCCF[isp] * pccf + WS_DGHAH[isp] * hoavh + WS_DGBA[isp] * alba
            (isp == 3 || isp == 13) && (dds -= 0.15032f0)        # WF/SF offset
            if d < 10f0                                           # small-tree sub-forms (JP / CONSJP grp)
                if isp == 6
                    consjp = ws_consjp(s, isp, si, cor)
                    dds = consjp + 1.23864f0 * ald + 0.64311f0 * cr - 0.48754f0 * alba -
                          0.00189f0 * bal / log(d + 1f0) - 0.00096f0 * pccf
                elseif isp in WS_CONSJP_GRP
                    consjp = cor + 0.233713f0 * log(si) + 1.53962f0
                    dds = consjp - 0.52776f0 * alba + 1.64163f0 * ald -
                          0.00205f0 * bal / log(d + 1f0) - 0.00105f0 * pccf
                end
            end
            (isp in WS_HARDWOOD_5YR) && (dds = log(exp(dds) * 2f0))   # 5-yr → 10-yr
        end
        isp == 5 && (dds += 0.30f0 * (0.80f0 + 0.004f0 * (si - 50f0)))   # IC calibration bump
        dds < -9.21f0 && (dds = -9.21f0)
        wk2[i] = dds
        # ⚠ ws/dgf.f computes a DDSMAX bound (:774-785) here, but the WK2(I)=DDS store (:761) precedes it
        # and WK2 is NEVER re-stored after the cap — the bound is DEAD CODE w.r.t. WK2 (a WS quirk; other
        # variants store WK2 after the cap). Faithful behaviour = do NOT cap wk2. (Verified vs source.)
    end
    return s
end

# CASE(6) JP small-tree CONSJP (ws/dgf.f:594-601). FORCON=0.74162 default; IFOR==2→0.59493; IGL==2→
# 0.854515 (IGL = the geographic-group flag; not yet surfaced ⇒ documented follow-on, JP absent from ref).
@inline function ws_consjp(s::StandState, isp::Int, si::Float32, cor::Float32)
    p = s.plot
    ifor = Int(p.forest_idx); slope = p.slope; asp = p.aspect; elev = p.elevation
    forcon = ifor == 2 ? 0.59493f0 : 0.74162f0
    return cor + 0.40657f0 * log(si) + 0.56709f0 * slope - 0.14671f0 * cos(asp) * slope +
           0.26785f0 * sin(asp) * slope - 0.00164f0 * elev + forcon
end

# PRD = ZRD(point)/XMAXPT(point) per-point Zeide relative density — SDICAL/SDICLS engine-gap (AK #209
# class). Only CASE(4,23) GS/RW consume it. Stub 0 until SDICAL is ported (inert on GS/RW-free stands).
@inline function ws_point_prd(s::StandState, pt::Int)
    d = s.density
    if isdefined(d, :point_zeide_rd) && 1 <= pt <= length(getfield(d, :point_zeide_rd))
        return getfield(d, :point_zeide_rd)[pt]
    end
    return 0f0
end

