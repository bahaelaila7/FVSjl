# =============================================================================
# diameter_growth.jl (centralcalifornia) — CA large-tree DDS (ca/dgf.f). Chunk 3.
#
# CA's dgf.f is GROUP-COMPRESSED (like WC/PN, unlike EC): 50 species map via MAPSPC onto 13 diameter-
# growth equations; every dgf coefficient array is length 13 and indexed by the growth group JSPC. DGDSQ
# is the one per-species (ISPC) term, loaded as DGDS[JSPC] in DGCONS. Three DDS branches:
#   • STANDARD (all species except the two below): the Wykoff ln(DDS) model
#       DDS = CONSPP + DGLD·ln(D) + CR·(DGCR + CR·DGCRSQ) + DGDSQ·D² + DGDBAL·BAL/ln(D+1)
#           + DGPCCF·PCCF + DGHAH·RELHT + DGLBA·ln(BA) + DGBAL·BAL          (CONSPP = DGCON + COR)
#   • TANOAK (ISPC 42, group 13): 5-yr equation → 10-yr basis  (DDS += ln 2)
#   • REDWOOD / GIANT SEQUOIA (ISPC 50 / 23, group 12): a DIFFERENT functional form on outside-bark
#     diameter increment (DGLT), converted to ΔDIB² via BRATIO; uses the point Zeide relative density PRD
#     (SDICAL/SDICLS) — deferred to the Zeide-SDI chunk, so the RW/GS branch is scaffolded but NOT yet
#     cyc-exact (cat01 has no RW/GS ⇒ inert for the beachhead).
# MEASURED bit-exact vs FVSca_dump (instrumented ca/dgf.f) per-tree LN(DDS) on cat01 (standard + TANOAK).
# =============================================================================

# ── ca/dgf.f DATA MAPSPC(50): species ISPC → growth group JSPC (1..13).
const CA_MAPSPC = Int[
    1, 1, 1, 2, 3, 3, 4, 7, 7, 6,
    5, 6, 5, 5, 9, 7, 8, 9, 9, 5,
    5, 2,12, 5, 9,10,10,10,10,10,
   10,10,10,10,10,10,11,11,11,10,
   10,13,10,10,10,10, 5,10,10,12]

# ── ca/dgf.f per-group scalar coefficients (length 13, index = JSPC).
const CA_DGLD   = Float32[0.950418,1.182104,1.186676,0.716226,1.077154,1.218279,0.886150,0.825682,0.738750,1.310111,0.955569,0.0,0.99531]
const CA_DGCR   = Float32[1.815305,2.856578,2.763519,3.272451,-0.276387,3.167164,1.478650,1.675208,3.454857,0.271183,0.0,0.0,2.08524]
const CA_DGCRSQ = Float32[0.0,-1.093354,-0.871061,-1.642904,1.063732,-1.568333,0.0,0.0,-1.773805,0.0,0.0,0.0,-0.98396]
const CA_DGSITE = Float32[0.820451,0.365679,0.492695,0.759305,0.0,0.566946,0.963375,0.724300,1.011504,0.213526,1.334008,0.0,0.00659]
const CA_DGDBAL = Float32[-0.005433,-0.005992,-0.003728,-0.008787,0.0,0.0,-0.006263,-0.002133,-0.013091,0.0,-0.005893,0.0,-0.00147]
const CA_DGLBA  = Float32[-0.000016,-0.058039,-0.122905,-0.028564,0.0,-0.267873,-0.129146,-0.203636,-0.131185,0.0,-0.408462,0.0,0.0]
const CA_DGBAL  = Float32[0.0,0.0,0.0,0.0,-0.000893,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const CA_DGPCCF = Float32[-0.000779,-0.001014,0.0,-0.000224,0.0,-0.000338,0.0,0.0,-0.000593,-0.000473,0.0,0.0,-0.00018]
const CA_DGHAH  = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.50155]
const CA_DGDS   = Float32[-0.0002385,-0.0006362,-0.0004572,-0.0002723,0.0,-0.0014178,-0.0002528,-0.0000731,-0.0004708,-0.0003048,0.0,0.0,-0.000373]
const CA_OBSERV = Float32[613.,3759.,2062.,5400.,84.,372.,561.,253.,2482.,306.,336.,8928.,6504.]
# DGCONS site/topography terms (length 13, index = JSPC).
const CA_DGCASP = Float32[0.0,-0.315227,-0.444594,-0.151727,0.649870,0.0,-0.280294,-0.179510,0.0,0.0,0.0,0.0,-0.19935]
const CA_DGSASP = Float32[0.0,0.097350,0.139180,0.018681,0.951834,0.0,-0.014463,-0.562259,0.0,0.0,0.0,0.0,-0.03587]
const CA_DGSLOP = Float32[0.0,-0.206267,0.0,-0.339369,0.0,0.0,-0.581722,-0.544867,0.0,0.0,0.0,0.0,0.73530]
const CA_DGSLSQ = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.99561]
const CA_DGEL   = Float32[0.0,0.0301,0.0248,-0.0141,0.0,0.0,0.0,0.0,-0.003784,0.0049,0.0,0.0,0.0]
const CA_DGELSQ = Float32[0.0,-0.00030732,-0.00033429,0.00024083,0.0,0.0,0.0,0.0,0.00006660,-0.00008781,0.0,0.0,0.0]

# DGFOR[jspc, isfor] (ca/dgf.f DATA DGFOR(5,13) → jl [jspc, isfor]; ISPFOR∈1..5 via MAPLOC).
const CA_DGFOR = Float32[
   -3.428338 -3.966547  0.0        0.0        0.0;       # 1  IC/PC/RC
   -2.108357  0.0        0.0        0.0        0.0;       # 2  WF/BR
   -2.073942 -1.943608   0.0        0.0        0.0;       # 3  RF/SH
   -1.877695 -2.099646  -2.211587  -1.955301  -2.078432; # 4  DF
    0.564402  0.0        0.0        0.0        0.0;       # 5  KP…
   -2.058828 -1.596998   0.0        0.0        0.0;       # 6  LP/WB
   -2.397678  0.0        0.0        0.0        0.0;       # 7  SP/WH/MH
   -1.626879  0.0        0.0        0.0        0.0;       # 8  WP
   -2.922255  0.0        0.0        0.0        0.0;       # 9  PP…
   -1.958189  0.0        0.0        0.0        0.0;       # 10 OAKS…
   -3.344700  0.0        0.0        0.0        0.0;       # 11 MA/GC/DG
    0.0       0.0        0.0        0.0        0.0;       # 12 GS/RW
   -0.94563   0.0        0.0        0.0        0.0]       # 13 TO

# MAPLOC[jspc, ifor] (ca/dgf.f DATA MAPLOC(10,13) → jl [jspc, ifor]; loc class 1..5, IFOR 1..10).
const CA_MAPLOC = Int[
    1 1 1 1 2 1 1 1 1 1;   # 1
    1 1 1 1 1 1 1 1 1 1;   # 2
    1 1 1 2 1 1 1 1 1 1;   # 3
    1 2 2 1 3 4 5 4 4 5;   # 4  DF: forest→loc class
    1 1 1 1 1 1 1 1 1 1;   # 5
    1 1 1 2 1 1 1 1 1 1;   # 6
    1 1 1 1 1 1 1 1 1 1;   # 7
    1 1 1 1 1 1 1 1 1 1;   # 8
    1 1 1 1 1 1 1 1 1 1;   # 9
    1 1 1 1 1 1 1 1 1 1;   # 10
    1 1 1 1 1 1 1 1 1 1;   # 11
    1 1 1 1 1 1 1 1 1 1;   # 12
    1 1 1 1 1 1 1 1 1 1]   # 13

const CA_RWGS = (23, 50)   # ISPC of Giant Sequoia / Redwood (group 12, special functional form)

"CA DGCONS: per-species (50) DGCON (ca/dgf.f ENTRY DGCONS), group-compressed. Stored in c.dg_const[isp]."
function ca_dgcons!(s::StandState)
    c = s.calib; p = s.plot
    ifor = Int(p.forest_idx)
    elev = p.elevation; slope = p.slope; asp = p.aspect
    sina = sin(asp); cosa = cos(asp)
    @inbounds for isp in 1:50
        jspc = CA_MAPSPC[isp]
        xsite = p.sp_site_index[isp]
        if isp == 23 || isp == 50                                  # GS / RW: ln(SITEAR) form, DGDSQ=0
            c.dg_const[isp] = -3.502444f0 + 0.415435f0 * log(xsite)
        else
            isfor = CA_MAPLOC[jspc, ifor]
            sasp = (CA_DGSASP[jspc] * sina + CA_DGCASP[jspc] * cosa + CA_DGSLOP[jspc]) * slope +
                   CA_DGSLSQ[jspc] * slope * slope
            c.dg_const[isp] = CA_DGFOR[jspc, isfor] + CA_DGEL[jspc] * elev + CA_DGELSQ[jspc] * elev * elev +
                              CA_DGSITE[jspc] * log(xsite) + sasp
        end
    end
    return s
end

# ca/dgf.f main body — per-tree WK2 = LN(DDS). Branches: STANDARD, TANOAK(42), RW/GS(23/50).
function dgf!(s::StandState, ::CentralCalifornia)
    p, t, c = s.plot, s.trees, s.calib
    dens = s.density; sd = s.coef.species
    wk2 = view(s.scratch.wk, 2, :)
    ba = p.basal_area; avh = p.avg_height
    slope = p.slope; asp = p.aspect; cosa = cos(asp)
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        isp = Int(t.species[i]); jspc = CA_MAPSPC[isp]
        pt_i = Int(t.plot_id[i])
        pccf = (1 <= pt_i <= length(dens.point_ccf)) ? dens.point_ccf[pt_i] : 0f0
        cr = Float32(t.crown_pct[i]) * 0.01f0
        pctfrac = 1f0 - t.crown_ratio[i] / 100f0
        bal = pctfrac * ba
        cor = c.dg_cor[isp]
        relht = avh > 0f0 ? min(t.height[i] / avh, 1.5f0) : 0f0
        if isp == 23 || isp == 50                                   # GS / RW (group 12): OB increment form
            conspp = c.dg_const[isp]                                # NB: no COR here (added to DDS below)
            pbal = (1 <= pt_i <= length(dens.point_ba)) ? dens.point_ba[pt_i] * pctfrac : bal
            pbal < 0f0 && (pbal = bal)
            prd = 0f0    # TODO(zeide-chunk): PRD = ZRD(pt)/XMAXPT(pt), point Zeide relative density (SDICLS)
            dglt = exp(conspp + 0.185911f0 * log(d) - 0.000073f0 * d * d - 0.001796f0 * pbal -
                       0.42078f0 * prd + 0.589318f0 * log(cr * 100f0) - 0.000926f0 * slope * 100f0 -
                       0.002203f0 * (slope * 100f0) * cosa)
            brat = wc_bratio(sd[:bark1][isp], sd[:bark2][isp], Int(sd[:bark_imap][isp]), d)
            tempd1 = d * brat; dup = d + dglt; tempd2 = dup * brat
            dds = log(tempd2 * tempd2 - tempd1 * tempd1) + cor + log(cor2_of(c, isp))
        else                                                        # STANDARD Wykoff ln(DDS)
            conspp = c.dg_const[isp] + cor
            dds = conspp + CA_DGLD[jspc] * log(d) +
                  cr * (CA_DGCR[jspc] + cr * CA_DGCRSQ[jspc]) +
                  CA_DGDS[jspc] * d * d + CA_DGDBAL[jspc] * bal / log(d + 1f0)
            dds += CA_DGPCCF[jspc] * pccf + CA_DGHAH[jspc] * relht +
                   CA_DGLBA[jspc] * log(ba) + CA_DGBAL[jspc] * bal
        end
        isp == 42 && (dds = log(exp(dds) * 2f0))                    # TANOAK: 5-yr → 10-yr basis
        dds < -9.21f0 && (dds = -9.21f0)
        wk2[i] = dds
    end
    return s
end
