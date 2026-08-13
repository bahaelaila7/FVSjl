# =============================================================================
# diameter_growth.jl (westcascades) — WC large-tree DDS (wc/dgf.f). Chunk 3.
#   wc_dgcons!(s)              — per-species (39) DGCON (wc/dgf.f ENTRY DGCONS): DGFOR(ISPFOR,JSPC)
#                               + DGEL·EL + DGEL2·EL² + DGSITE·ln(XSITE) + SASP ; RW (ISPC17) = ln-site.
#   dgf!(s, ::WestCascades)    — per-tree WK2 = LN(DDS) (wc/dgf.f main body).
# WC uses 19 SPECIES GROUPS mapped from 39 species via MAPSPC. DEFAULT branch (JSPC≠13, ISPC≠17):
#   DDS = CONSPP + DGLD·lnD + CR·(DGCR + CR·DGCRSQ) + DGDSQ·D² + DGDBAL·BAL/ln(D+1)
#       + DGPCCF·PCCF + DGHAH·RELHT + DGLBA·lnBA + DGBAL·BAL + DGBA·BA ; CONSPP = DGCON + COR.
#   RA (group 13) = red-alder DIAGR eq ; RW (ISPC 17) = DGLT exp eq (both use BRATIO bark).
#
# VALIDATION (2026-08-13, FVSwc_clean wct01 DEBUG DGF, block-3 growth pass):
#   • DEFAULT DDS formula: 27/27 trees BIT-EXACT vs live 9001 LN(DDS) (worst 0.00007 = CONSPP
#     F11.4 print rounding), across groups 2,5,6,7,11,16 (species WF/GF, SP, PP, DF, ES, LP).
#   • DGCON reconstruction: bit-exact vs live 9030 dump (e.g. DF group7 DGCON=1.01298;
#     WO group18 King's-SI XSITE transform 48.79→43.466, DGCON=-0.76738).
#   • RA (group 13) + RW (ISPC 17) branches: ported source-faithful, NOT live-validated (no
#     RA/RW in wct01/nct01) — pending a red-alder / redwood stand.
# =============================================================================

# 19-group scalar coefficients (wc/dgf.f DATA), index = JSPC = WC_MAPSPC[ISPC].
const WC_DGLD   = Float32[0.527758,0.905119,0.993986,0.904253,0.844690,0.738750,0.534138,0.843013,0.722462,0.857131,0.879338,1.024186,0.0,0.889596,0.816880,0.478504,0.949631,1.66609,0.0]
const WC_DGCR   = Float32[2.982807,1.754811,1.522401,4.123101,1.597250,3.454857,1.636854,2.878032,2.160348,1.505513,1.970052,0.459387,0.0,1.732535,2.471226,1.905011,1.826879,0.0,0.0]
const WC_DGCRSQ = Float32[-1.331331,0.0,0.0,-2.689340,0.0,-1.773805,-0.045578,-1.631418,-0.834196,0,0,0,0,0,0,0,0,0,0]
const WC_DGSITE = Float32[0.534255,0.318254,0.349888,0.684939,0.404010,1.011504,1.020863,0.139734,0.380416,0.208040,0.252853,1.965888,0.0,0.227307,0.244694,0.391327,0.375175,0.14995,0.0]
const WC_DGDBAL = Float32[-0.011247,-0.005355,-0.002979,-0.006368,-0.003726,-0.013091,-0.009363,-0.003923,-0.004065,-0.004101,-0.004215,-0.010222,0.0,-0.001265,-0.005950,-0.004706,-0.005350,0.0,0.0]
const WC_DGLBA  = Float32[-0.030730,0,0,0,0,-0.131185,0,0,0,0,0,0,0,0,0,0,0,0,0]
const WC_DGBA   = Float32[0,0,-0.000137,0,0,0,-0.000215,0,0,0,-0.000173,0,0,-0.000981,-0.000147,-0.000114,0.000040,-0.00204,0.0]
const WC_DGBAL  = Float32[0.002839,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.00326,0.0]
const WC_DGPCCF = Float32[0,0,0,-0.000471,-0.000257,-0.000593,0,-0.000552,0,-0.000201,0,-0.000757,0,0,0,0,0,0,0]
const WC_DGHAH  = Float32[0,-0.000661,0,0,0,0,0,0,-0.000358,0,0,0,0,0,0,0,0,0,0]
const WC_DGEL   = Float32[-0.048852,-0.003051,-0.003773,-0.069045,-0.023376,-0.003784,-0.037591,-0.050081,-0.040067,-0.003809,0.0,-0.012111,0.0,-0.075986,0.0,-0.005414,0.323546,0.0,0.0]
const WC_DGEL2  = Float32[0.000478,0,0,0.000608,0,0.00006660,0.000549,0.000660,0.000395,0,0,0,0,0.001193,0,0,-0.003130,0,0]
const WC_DGSASP = Float32[0,0,0.022160,-0.207659,0,0,-0.038992,0,0,-0.126130,0,0,0,-0.863980,0.679903,0.378860,0.202507,0,0]
const WC_DGCASP = Float32[0,0,-0.782418,-0.374512,0,0,-0.080943,0,0,-0.104495,0,0,0,0.085958,-0.023186,0.207853,-0.935870,0,0]
const WC_DGSLOP = Float32[0.245548,0,0.319956,0.400223,0,0,0.077787,0,0.421486,0.411602,0,0,0,0,0,-0.066440,0,0,0]
const WC_DGSLSQ = Float32[0,0,0,0,0,0,-0.215778,0,-0.693610,0,0,0,0,0,0,0,0,0,0]
# MAPSPC: 39 species → 19 groups (wc/dgf.f DATA MAPSPC).
const WC_MAPSPC = Int[1,2,2,3,17,17,4,15,11,11,16,6,5,5,6,7,19,8,9,10,12,13,14,14,14,14,14,18,14,11,11,11,11,14,14,14,14,14,14]
# DGFOR[jspc, isfor] (6 loc classes × 19 groups) — wc/dgf.f DATA DGFOR (column-major → rows here).
const WC_DGFOR = Float32[
    -0.619069 -0.479015 -0.291244 0.0 -0.420228 -0.746419;   # g1
    -0.643920 0.0 0.0 0.0 0.0 0.0;                            # g2
    -1.888949 -1.276180 0.0 0.0 0.0 0.0;                      # g3
    -1.401865 -1.127977 0.0 0.0 0.0 0.0;                      # g4
    -0.589570 -0.909553 0.0 0.0 0.0 0.0;                      # g5
    -2.922255 0.0 0.0 0.0 0.0 0.0;                            # g6
    -2.750874 -2.787499 -2.672664 -2.533437 -2.693964 -2.718852; # g7
     0.412763  0.645645 0.0 0.0 0.0 0.0;                      # g8
    -0.298310 -0.147675 -0.006413 0.0 0.0 0.0;                # g9
    -1.052161 -0.793945 0.0 0.0 0.0 0.0;                      # g10
    -1.310067 -1.432659 0.0 0.0 0.0 0.0;                      # g11
    -7.753469 -8.279266 0.0 0.0 0.0 0.0;                      # g12
     0.0 0.0 0.0 0.0 0.0 0.0;                                 # g13
    -0.107648 -0.098335 0.0 0.0 0.0 0.0;                      # g14
    -1.277664 -1.178041 0.0 0.0 0.0 0.0;                      # g15
    -0.524624 -0.803095 0.0 0.0 0.0 0.0;                      # g16
    -9.211184 -9.800653 0.0 0.0 0.0 0.0;                      # g17
    -1.33299 0.0 0.0 0.0 0.0 0.0;                             # g18
     0.0 0.0 0.0 0.0 0.0 0.0]                                 # g19
# MAPLOC[jspc, jfor] (forest → location class) — wc/dgf.f DATA MAPLOC.
const WC_MAPLOC = Int[
    1 2 3 4 5 6;  1 1 1 1 1 1;  1 1 1 2 1 1;  1 1 2 2 2 2;  1 1 1 1 1 2;
    1 1 1 1 1 1;  1 2 3 4 5 6;  1 1 2 1 1 1;  1 1 2 3 2 1;  1 1 2 2 1 1;
    1 1 1 1 1 2;  1 1 1 1 1 2;  1 1 1 1 1 1;  1 1 1 1 1 2;  1 1 1 1 1 2;
    1 1 1 2 2 2;  1 1 1 1 2 1;  1 1 1 1 1 1;  1 1 1 1 1 1]
# MAPDSQ[jspc, jfor] (forest → DBH² coeff index) — all 1 except g1 forest4=2.
const WC_MAPDSQ = Int[
    1 1 1 2 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;
    1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;
    1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;
    1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1;  1 1 1 1 1 1]
# DGDS[jspc, idsq] (2 DBH² coeffs × 19 groups) — col 2 all 0 (wc/dgf.f, blowup term disabled).
const WC_DGDS = Float32[
    -0.0001983 0.0; -0.0003137 0.0; -0.0002621 0.0; -0.0003996 0.0; -0.0000596 0.0;
    -0.0004708 0.0; -0.0001039 0.0; -0.0000644 0.0; -0.0001546 0.0; -0.0002214 0.0;
    -0.0001323 0.0; -0.0001737 0.0;  0.0 0.0;  0.0 0.0; -0.0002536 0.0;
     0.0 0.0; -0.0003552 0.0; -0.00154 0.0;  0.0 0.0]

"WC IFOR→JFOR (national-forest/BLM-unit alignment, wc/dgf.f DGCONS)."
@inline function wc_jfor(ifor::Int)
    ifor <= 6 && return ifor
    ifor == 7 && return 3
    ifor == 8 && return 6
    ifor == 9 && return 5
    ifor == 10 && return 4
    return 3
end

"WC bark BRATIO (wc/bratio.f eqtypes) — a+b·D (type1/2) or POWER (type3), returned as DIB/D."
@inline function wc_bratio(a::Float32, b::Float32, eqtype::Int, d::Float32)
    d <= 0f0 && return 0.99f0
    eqtype == 1 && return (d - (a + b * d)) / d
    eqtype == 2 && return (a + b * d) / d
    return (a * d^b) / d
end

"WC DGCONS: per-species (39) DGCON (wc/dgf.f ENTRY DGCONS). Stored in c.dg_const[sp]."
function wc_dgcons!(s::StandState)
    c = s.calib; p = s.plot
    jfor = wc_jfor(Int(p.forest_idx))
    elev = p.elevation; slope = p.slope; asp = p.aspect
    sina = sin(asp); cosa = cos(asp)
    @inbounds for isp in 1:39
        jspc = WC_MAPSPC[isp]
        si = p.sp_site_index[isp]
        if isp == 17                                  # REDWOOD — ln-site const, no DGDSQ
            c.dg_const[isp] = -3.502444f0 + 0.415435f0 * log(max(si, 1f0))
            continue
        end
        isfor = WC_MAPLOC[jspc, jfor]
        sasp = (WC_DGSASP[jspc] * sina + WC_DGCASP[jspc] * cosa + WC_DGSLOP[jspc]) * slope +
               WC_DGSLSQ[jspc] * slope * slope
        xsite = si
        jspc == 10 && (xsite = xsite * 3.281f0)                                    # ES: m→ft site
        jspc == 18 && (xsite = -37.60812f0 * log(1f0 - (xsite / 114.24569f0)^0.4444f0))  # WO: King's DF SI
        temel = elev
        (jspc == 14 && temel > 30f0) && (temel = 30f0)                             # group-14 elev cap
        c.dg_const[isp] = WC_DGFOR[jspc, isfor] + WC_DGEL[jspc] * temel +
                          WC_DGEL2[jspc] * temel * temel + WC_DGSITE[jspc] * log(max(xsite, 1f0)) + sasp
    end
    return s
end

"WC dgf! — per-tree WK2 = LN(DDS) (wc/dgf.f main body; DEFAULT + RA group-13 + RW ISPC-17)."
function dgf!(s::StandState, ::WestCascades)
    p, t, c = s.plot, s.trees, s.calib
    dens = s.density; sd = s.coef.species
    wk2 = view(s.scratch.wk, 2, :)
    ba = p.basal_area; avh = p.avg_height
    slope = p.slope; asp = p.aspect
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        isp = Int(t.species[i]); jspc = WC_MAPSPC[isp]
        pt_i = Int(t.plot_id[i])
        pccf = (1 <= pt_i <= length(dens.point_ccf)) ? dens.point_ccf[pt_i] : 0f0
        ptba = (1 <= pt_i <= length(dens.point_ba)) ? dens.point_ba[pt_i] : ba
        cr = Float32(t.crown_pct[i]) * 0.01f0
        pctfrac = 1f0 - t.crown_ratio[i] / 100f0       # 1 − PCT/100 (crown_ratio field holds PCT)
        bal = pctfrac * ba
        cor = c.dg_cor[isp]
        if jspc == 13                                  # RED ALDER (wc/dgf.f DIAGR eq)
            et = Int(sd[:bark_imap][22]); brat = wc_bratio(sd[:bark1][22], sd[:bark2][22], et, d)
            const0 = 3.250531f0 - 0.003029f0 * ba
            diagr = d <= 18f0 ? const0 - 0.166496f0 * d + 0.004618f0 * d * d :
                                const0 - (const0 / 10f0) * (d - 18f0)
            diagr < 0.1f0 && (diagr = 0.1f0)
            dds = log(diagr * (2f0 * d * brat + diagr)) + log(cor2_of(c, isp)) + cor
        elseif isp == 17                               # REDWOOD (wc/dgf.f DGLT exp eq)
            conspp = c.dg_const[isp]                    # RW: COR applied AFTER, not in CONSPP
            pbal = ptba * pctfrac; pbal < 0f0 && (pbal = bal)
            prd = 0f0                                   # PRD point-Zeide (TODO precompute; 0 baseline)
            dglt = exp(conspp + 0.185911f0 * log(d) - 0.000073f0 * d * d - 0.001796f0 * pbal -
                       0.42078f0 * prd + 0.589318f0 * log(cr * 100f0) - 0.000926f0 * slope * 100f0 -
                       0.002203f0 * (slope * 100f0) * cos(asp))
            et = Int(sd[:bark_imap][isp]); brat = wc_bratio(sd[:bark1][isp], sd[:bark2][isp], et, d)
            t1 = d * brat; t2 = (d + dglt) * brat
            dds = log(t2 * t2 - t1 * t1) + cor + log(cor2_of(c, isp))
        else                                           # DEFAULT — validated bit-exact vs live
            conspp = c.dg_const[isp] + cor
            relht = avh > 0f0 ? min(t.height[i] / avh, 1.5f0) : 0f0
            dgdsq = WC_DGDS[jspc, WC_MAPDSQ[jspc, wc_jfor(Int(p.forest_idx))]]
            dds = conspp + WC_DGLD[jspc] * log(d) + cr * (WC_DGCR[jspc] + cr * WC_DGCRSQ[jspc]) +
                  dgdsq * d * d + WC_DGDBAL[jspc] * bal / log(d + 1f0) +
                  WC_DGPCCF[jspc] * pccf + WC_DGHAH[jspc] * relht + WC_DGLBA[jspc] * log(ba) +
                  WC_DGBAL[jspc] * bal + WC_DGBA[jspc] * ba
        end
        dds < -9.21f0 && (dds = -9.21f0)               # shared final clamp (wc/dgf.f:445)
        wk2[i] = dds
    end
    return s
end

# COR2 (READCORD/REUSCORD multiplicative calib; 1.0 baseline unless a WC keyfile sets it).
@inline cor2_of(c, isp::Int) = (isdefined(c, :dg_cor2) && length(c.dg_cor2) >= isp) ? c.dg_cor2[isp] : 1f0
