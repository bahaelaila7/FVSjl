# =============================================================================
# diameter_growth.jl (pacificnorthwest) — PN large-tree DDS (pn/dgf.f). Chunk 3.
#
# Standard westside Wykoff LN(DDS), same structure as WC, but PN is a "REFIT 060396" with 20 species
# GROUPS (WC had 19): SS is a NEW group 18 ⇒ WO→group 19, RW→group 20. Differences from wc_dgcons!:
#   • the WO King's-SI XSITE transform is JSPC==19 (was 18); • NO JFOR (BLM→NF) remap — pn/dgf.f:510
#     uses IFOR directly (pn_jfor = identity); • PN_DGFOR is 3-loc-class (20×3), PN_DGDS col-2 nonzero
#     for g7/g8/g9/g18. DEFAULT DDS + RA(group 13) + RW(ISPC 17) branches are byte-identical to WC (bark
#     via the shared wc_bratio, which reads the per-species CSV columns ⇒ variant-generic). ES(jspc 10)
#     ×3.281 m→ft, group-14 elev cap, RW ISPC-17 ln-site const — all identical to WC.
# MEASURED (FVSpn_g16 dgf dump): validated bit-exact vs the live cyc1 growth-pass LN(DDS) (see harness).
# =============================================================================

# ── PN 20-group scalar coefficients (pn/dgf.f DATA), index = JSPC = PN_MAPSPC[ISPC].
const PN_DGLD   = Float32[0.919402,0.905119,0.993986,0.904253,0.844690,0.738750,0.802905,0.744005,0.641956,0.857131,0.879338,1.024186,0.511442,0.889596,0.816880,0.478504,0.949631,1.049845,1.66609,0.0]
const PN_DGCR   = Float32[1.290568,1.754811,1.522401,4.123101,1.597250,3.454857,1.936912,0.771395,1.471926,1.505513,1.970052,0.459387,0.623093,1.732535,2.471226,1.905011,1.826879,1.632468,0.0,0.0]
const PN_DGCRSQ = Float32[0.125823,0.0,0.0,-2.689340,0.0,-1.773805,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const PN_DGSITE = Float32[0.541881,0.318254,0.349888,0.684939,0.404010,1.011504,0.495162,0.708166,0.634098,0.208040,0.252853,1.965888,0.237269,0.227307,0.244694,0.391327,0.375175,0.0,0.14995,0.0]
const PN_DGDBAL = Float32[-0.002133,-0.005355,-0.002979,-0.006368,-0.003726,-0.013091,-0.001827,-0.016240,-0.012589,-0.004101,-0.004215,-0.010222,-0.027074,-0.001265,-0.005950,-0.004706,-0.005350,-0.000086,0.0,0.0]
const PN_DGLBA  = Float32[-0.136818,0.0,0.0,0.0,0.0,-0.131185,-0.129474,-0.130036,-0.085525,0.0,0.0,0.0,-0.481983,0.0,0.0,0.0,0.0,-0.198636,0.0,0.0]
const PN_DGBA   = Float32[0.0,0.0,-0.000137,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.000173,0.0,0.0,-0.000981,-0.000147,-0.000114,0.000040,0.0,-0.00204,0.0]
const PN_DGBAL  = Float32[0.0,0.0,0.0,0.0,0.0,0.0,-0.001639,0.003883,0.002385,0.0,0.0,0.0,0.008903,0.0,0.0,0.0,0.0,-0.002319,-0.00326,0.0]
const PN_DGPCCF = Float32[0.0,0.0,0.0,-0.000471,-0.000257,-0.000593,0.0,0.0,0.0,-0.000201,0.0,-0.000757,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const PN_DGHAH  = Float32[0.0,-0.000661,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const PN_DGEL   = Float32[-0.023858,-0.003051,-0.003773,-0.069045,-0.023376,-0.003784,-0.009845,-0.009564,-0.018444,-0.003809,0.0,-0.012111,0.0,-0.075986,0.0,-0.005414,0.323546,0.007009,0.0,0.0]
const PN_DGEL2  = Float32[0.0,0.0,0.0,0.000608,0.0,0.0000666,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.001193,0.0,0.0,-0.003130,0.0,0.0,0.0]
const PN_DGSASP = Float32[0.096326,0.0,0.022160,-0.207659,0.0,0.0,0.003263,-0.106020,0.061254,-0.126130,0.0,0.0,-0.085538,-0.863980,0.679903,0.378860,0.202507,0.100081,0.0,0.0]
const PN_DGCASP = Float32[-0.217205,0.0,-0.782418,-0.374512,0.0,0.0,0.014165,-0.106936,-0.056608,-0.104495,0.0,0.0,0.022254,0.085958,-0.023186,0.207853,-0.935870,-0.221095,0.0,0.0]
const PN_DGSLOP = Float32[-0.265612,0.0,0.319956,0.400223,0.0,0.0,-0.340401,-0.303490,0.736143,0.411602,0.0,0.0,0.0,0.0,0.0,-0.066440,0.0,-0.169141,0.0,0.0]
const PN_DGSLSQ = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.082191,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]

const PN_MAPSPC = Int[1,2,2,3,4,18,4,15,11,11,16,6,5,5,6,7,20,8,9,10,12,13,14,14,14,14,14,19,14,11,11,11,11,14,14,14,14,14,14]

# DGFOR[jspc, isfor] (pn/dgf.f DATA DGFOR(3,20) → jl [jspc,isfor], 3 location classes).
const PN_DGFOR = Float32[
    -0.627531  0.0       0.0;       -0.643920  0.0       0.0;       -1.888949 -1.276180  0.0;
    -1.401865 -1.127977  0.0;       -0.589570 -0.909553  0.0;       -2.922255  0.0       0.0;
    -0.739354 -0.199200  0.0;       -0.688250 -0.405590  0.0;       -0.594460 -0.522658  0.0;
    -1.052161 -0.793945  0.0;       -1.310067 -1.432659  0.0;       -7.753469 -8.279266  0.0;
     4.253807  3.913250  3.507520;  -0.107648 -0.098335  0.0;       -1.277664 -1.178041  0.0;
    -0.524624 -0.803095  0.0;       -9.211184 -9.800653  0.0;        2.075598  2.100904  0.0;
    -1.33299   0.0       0.0;        0.0       0.0       0.0]

# MAPLOC[jspc, ifor] (pn/dgf.f DATA MAPLOC(6,20) → jl [jspc,ifor]; values 1..3 = loc class; IFOR used raw).
const PN_MAPLOC = Int[
    1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;  1 2 1 2 2 2;  1 2 1 2 2 2;
    1 1 1 1 1 1;  1 2 2 2 2 2;  1 2 1 2 2 2;  1 2 1 2 2 2;  1 1 1 1 1 1;
    1 2 1 2 2 2;  1 2 1 2 2 2;  1 2 3 3 3 3;  1 2 1 2 2 2;  1 2 1 2 2 2;
    1 2 1 2 2 2;  1 1 1 1 1 1;  1 2 1 2 2 2;  1 1 1 1 1 1;  1 1 1 1 1 1]
const PN_MAPDSQ = Int[
    1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;
    1 1 1 1 1 1;  1 2 1 2 2 2;  1 2 1 2 2 2;  1 2 1 2 2 2;  1 1 1 1 1 1;
    1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 2 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;
    1 1 1 1 1 1;  1 1 1 1 1 1;  1 2 1 2 2 2;  1 1 1 1 1 1;  1 1 1 1 1 1]
# DGDS[jspc, idsq] (pn/dgf.f DATA DGDS(3,20) → jl [jspc,idsq]; col 2 nonzero for g7/g8/g9/g18).
const PN_DGDS = Float32[
    -0.0002641  0.0        0.0;   -0.0003137  0.0        0.0;   -0.0002621  0.0        0.0;
    -0.0003996  0.0        0.0;   -0.0000596  0.0        0.0;   -0.0004708  0.0        0.0;
    -0.0000896 -0.0000641  0.0;   -0.0000572 -0.0000862  0.0;   -0.0001736 -0.0001040  0.0;
    -0.0002214  0.0        0.0;   -0.0001323  0.0        0.0;   -0.0001737  0.0        0.0;
    -0.0005099  0.0        0.0;    0.0        0.0        0.0;   -0.0002536  0.0        0.0;
     0.0        0.0        0.0;   -0.0003552  0.0        0.0;   -0.0002123 -0.0001361  0.0;
    -0.00154    0.0        0.0;    0.0        0.0        0.0]

# pn/dgf.f ENTRY DGCONS — per-species DGCON. PN uses IFOR directly (no JFOR remap); WO King's-SI = jspc 19.
function pn_dgcons!(s::StandState)
    c = s.calib; p = s.plot
    ifor = Int(p.forest_idx)
    elev = p.elevation; slope = p.slope; asp = p.aspect
    sina = sin(asp); cosa = cos(asp)
    @inbounds for isp in 1:39
        jspc = PN_MAPSPC[isp]
        si = p.sp_site_index[isp]
        if isp == 17                                   # REDWOOD — ln-site const, no DGDSQ (== WC)
            c.dg_const[isp] = -3.502444f0 + 0.415435f0 * log(max(si, 1f0))
            continue
        end
        isfor = PN_MAPLOC[jspc, ifor]
        sasp = (PN_DGSASP[jspc] * sina + PN_DGCASP[jspc] * cosa + PN_DGSLOP[jspc]) * slope +
               PN_DGSLSQ[jspc] * slope * slope
        xsite = si
        jspc == 10 && (xsite = xsite * 3.281f0)                                    # ES: m→ft
        jspc == 19 && (xsite = -37.60812f0 * log(1f0 - (xsite / 114.24569f0)^0.4444f0))  # WO King's DF SI (PN jspc 19)
        temel = elev
        (jspc == 14 && temel > 30f0) && (temel = 30f0)                             # group-14 elev cap
        c.dg_const[isp] = PN_DGFOR[jspc, isfor] + PN_DGEL[jspc] * temel +
                          PN_DGEL2[jspc] * temel * temel + PN_DGSITE[jspc] * log(max(xsite, 1f0)) + sasp
    end
    return s
end

# pn/dgf.f main body — per-tree WK2 = LN(DDS). DEFAULT + RA(group 13) + RW(ISPC 17), identical to WC.
function dgf!(s::StandState, ::PacificNorthwest)
    p, t, c = s.plot, s.trees, s.calib
    dens = s.density; sd = s.coef.species
    wk2 = view(s.scratch.wk, 2, :)
    ba = p.basal_area; avh = p.avg_height
    slope = p.slope; asp = p.aspect
    ifor = Int(p.forest_idx)
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        isp = Int(t.species[i]); jspc = PN_MAPSPC[isp]
        pt_i = Int(t.plot_id[i])
        pccf = (1 <= pt_i <= length(dens.point_ccf)) ? dens.point_ccf[pt_i] : 0f0
        ptba = (1 <= pt_i <= length(dens.point_ba)) ? dens.point_ba[pt_i] : ba
        cr = Float32(t.crown_pct[i]) * 0.01f0
        pctfrac = 1f0 - t.crown_ratio[i] / 100f0
        bal = pctfrac * ba
        cor = c.dg_cor[isp]
        if jspc == 13                                  # RED ALDER (pn/dgf.f DIAGR eq == WC)
            et = Int(sd[:bark_imap][22]); brat = wc_bratio(sd[:bark1][22], sd[:bark2][22], et, d)
            const0 = 3.250531f0 - 0.003029f0 * ba
            diagr = d <= 18f0 ? const0 - 0.166496f0 * d + 0.004618f0 * d * d :
                                const0 - (const0 / 10f0) * (d - 18f0)
            diagr < 0.1f0 && (diagr = 0.1f0)
            dds = log(diagr * (2f0 * d * brat + diagr)) + log(cor2_of(c, isp)) + cor
        elseif isp == 17                               # REDWOOD (pn/dgf.f DGLT exp eq == WC)
            conspp = c.dg_const[isp]
            pbal = ptba * pctfrac; pbal < 0f0 && (pbal = bal)
            prd = 0f0                                   # PRD point-Zeide (TODO precompute; 0 baseline, same as WC)
            dglt = exp(conspp + 0.185911f0 * log(d) - 0.000073f0 * d * d - 0.001796f0 * pbal -
                       0.42078f0 * prd + 0.589318f0 * log(cr * 100f0) - 0.000926f0 * slope * 100f0 -
                       0.002203f0 * (slope * 100f0) * cos(asp))
            et = Int(sd[:bark_imap][isp]); brat = wc_bratio(sd[:bark1][isp], sd[:bark2][isp], et, d)
            t1 = d * brat; t2 = (d + dglt) * brat
            dds = log(t2 * t2 - t1 * t1) + cor + log(cor2_of(c, isp))
        else                                           # DEFAULT
            conspp = c.dg_const[isp] + cor
            relht = avh > 0f0 ? min(t.height[i] / avh, 1.5f0) : 0f0
            dgdsq = PN_DGDS[jspc, PN_MAPDSQ[jspc, ifor]]
            dds = conspp + PN_DGLD[jspc] * log(d) + cr * (PN_DGCR[jspc] + cr * PN_DGCRSQ[jspc]) +
                  dgdsq * d * d + PN_DGDBAL[jspc] * bal / log(d + 1f0) +
                  PN_DGPCCF[jspc] * pccf + PN_DGHAH[jspc] * relht + PN_DGLBA[jspc] * log(ba) +
                  PN_DGBAL[jspc] * bal + PN_DGBA[jspc] * ba
        end
        dds < -9.21f0 && (dds = -9.21f0)
        wk2[i] = dds
    end
    return s
end
