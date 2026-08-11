# =============================================================================
# diameter_growth.jl (klamath) — NC large-tree DDS (nc/dgf.f). Chunk 3.
#   nc_bratio(a,b,eqtype,d)   — NC bark (BRKRAT eqtype 1/2/3; type3 = POWER, sp12 redwood)
#   nc_dgcons!(s)             — per-species DGCON (nc/dgf.f ENTRY DGCONS; 3 branches)
#   dgf!(s, ::Klamath)        — per-tree WK2 = LN(DDS) (nc/dgf.f main body; 3 DDS branches, 5-yr TDDS/2)
# DEFAULT (sp1,3,4,5,7,8,10,11): DDS = CONSPP + DGLD·lnD + CR·(DGCR+CR·DGCRSQ) + DGDSQ·D² + DGDBAL·BAL/ln(D+1)
#   + DGPCCF·PCCF + DGHAH·HOAVH ; sp4(WF) −0.15032. CONSPP = DGCON + COR + DGCCFA·ALRD + DGBA·ALBA.
# sp2,6,9: DDS = CONSPP + DGLD2·lnD + DGDSQ2·D²/1000 + DGCR2·CRID + DGDBA2·PBAL/ln(D+1)/100 + DGBA2·ALPBA.
# sp12 RW: DGLT=EXP(CONSPP+…); DDS=LN((D+DGLT)²·BRAT²−(D·BRAT)²)+COR+LN(COR2).  All → TDDS=EXP(DDS); DDS=LN(TDDS/2).
# =============================================================================

const NC_DGLD   = Float32[0.88425,0.0,0.86990,1.01718,1.14082,0.0,1.23911,0.99531,0.0,0.96865,0.99531,0.0]
const NC_DGCR   = Float32[2.83271,0.0,2.96040,3.01884,2.82796,0.0,-1.20841,2.08524,0.0,1.5466,2.08524,0.0]
const NC_DGCRSQ = Float32[-0.84141,0.0,-1.08219,-1.12464,-2.14739,0.0,2.31782,-0.98396,0.0,0.07152,-0.98396,0.0]
const NC_DGDBAL = Float32[-0.00358,0.0,-0.00443,-0.00257,-0.00126,0.0,-0.00199,-0.00147,0.0,-0.00408,-0.00147,0.0]
const NC_DGBA   = Float32[0.0,0.0,-0.01744,-0.16596,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const NC_DGHAH  = Float32[0.0,0.0,0.0,0.0,0.56348,0.0,0.0,0.50155,0.0,0.0,0.50155,0.0]
const NC_DGPCCF = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.00018,0.0,-0.00002,-0.0018,0.0]
const NC_DGCCFA = Float32[-0.06784,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
# set-2 (sp2,6,9)
const NC_DGLD2  = Float32[1.52284,1.26883,1.92426,1.53339,1.26883,1.41389,1.41389,1.52284,1.53339,1.63568,1.41389,0.0]
const NC_DGCR2  = Float32[0.51033,0.27986,0.40047,0.35739,0.27986,0.32660,0.32660,0.51033,0.35739,0.27516,0.32660,0.0]
const NC_DGDSQ2 = Float32[-0.26590,-0.35325,-0.44612,-0.47442,-0.35325,-0.48938,-0.48938,-0.26590,-0.47442,-0.54162,-0.48938,0.0]
const NC_DGBA2  = Float32[-0.35579,0.0,-0.24596,-0.12359,0.0,-0.25287,-0.25287,-0.35579,-0.12359,-0.20902,-0.25287,0.0]
const NC_DGDBA2 = Float32[0.0,-0.79922,0.0,-0.44256,-0.79922,-0.16000,-0.16000,0.0,-0.44256,-0.32497,-0.16000,0.0]
# forest-dependent DGCON: DGFOR[locclass,sp] via MAPLOC[ifor,sp]; DGDSQ = DGDS[idx,sp] via MAPDSQ[ifor,sp]
const NC_DGFOR = Float32[
    -2.00201 -2.19449 -1.84083 0.0 0.0 0.0;      # sp1 (6 loc classes)
     0.0      0.0      0.0     0.0 0.0 0.0;       # sp2
    -2.54402 -2.41928 -2.75656 0.0 0.0 0.0;       # sp3
    -1.88042 -2.06853 -1.69815 0.0 0.0 0.0;       # sp4
    -1.69950  0.0      0.0     0.0 0.0 0.0;       # sp5
     0.0      0.0      0.0     0.0 0.0 0.0;       # sp6
    -2.68349  0.0      0.0     0.0 0.0 0.0;       # sp7
    -0.94563  0.0      0.0     0.0 0.0 0.0;       # sp8
     0.0      0.0      0.0     0.0 0.0 0.0;       # sp9
    -4.6744   0.0      0.0     0.0 0.0 0.0;       # sp10
    -0.94563  0.0      0.0     0.0 0.0 0.0;       # sp11
     0.0      0.0      0.0     0.0 0.0 0.0]       # sp12  (indexed [sp, locclass])
const NC_DGDS = Float32[
    -0.000328 -0.000248 0.0 0.0;                  # sp1 (4 idx)
     0.0       0.0      0.0 0.0;                   # sp2
    -0.000313  0.0      0.0 0.0;                   # sp3
    -0.000356 -0.000268 0.0 0.0;                   # sp4
    -0.000875  0.0      0.0 0.0;                   # sp5
     0.0       0.0      0.0 0.0;                   # sp6
    -0.000338  0.0      0.0 0.0;                   # sp7
    -0.000373  0.0      0.0 0.0;                   # sp8
     0.0       0.0      0.0 0.0;                   # sp9
    -0.000728  0.0      0.0 0.0;                   # sp10
    -0.000373  0.0      0.0 0.0;                   # sp11
     0.0       0.0      0.0 0.0]                    # sp12  (indexed [sp, idx])
# MAPLOC[ifor,sp]/MAPDSQ[ifor,sp]: 7 forests × 12 sp. Only sp1/sp4 vary by forest; rest all 1.
const NC_MAPLOC = Int[  # [sp, forest] (forest 1..7)
    1 1 1 2 3 3 2;      # sp1
    1 1 1 1 1 1 1;      # sp2
    1 1 1 2 3 3 2;      # sp3
    1 1 1 2 3 3 2;      # sp4
    1 1 1 1 1 1 1;      # sp5
    1 1 1 1 1 1 1;      # sp6
    1 1 1 1 1 1 1;      # sp7
    1 1 1 1 1 1 1;      # sp8
    1 1 1 1 1 1 1;      # sp9
    1 1 1 1 1 1 1;      # sp10
    1 1 1 1 1 1 1;      # sp11
    1 1 1 1 1 1 1]      # sp12
const NC_MAPDSQ = Int[  # [sp, forest]
    1 1 1 1 2 2 1;      # sp1
    1 1 1 1 1 1 1;      # sp2
    1 1 1 1 1 1 1;      # sp3
    1 1 1 1 2 2 1;      # sp4
    1 1 1 1 1 1 1;1 1 1 1 1 1 1;1 1 1 1 1 1 1;1 1 1 1 1 1 1;1 1 1 1 1 1 1;1 1 1 1 1 1 1;1 1 1 1 1 1 1;1 1 1 1 1 1 1]
# sp2/6/9 site DGCON: DGLAT2(ILAT=5,sp) + DGEL2·ELEV + DGSLP2·SLOPE + DGSLQ2·SLOPE² + DGSITE·SITEAR
const NC_DGLAT2_5 = Float32[0.1630,-0.4297,-0.1043,0.1434,-0.4297,0.0540,0.0540,0.1630,0.1434,-0.3995,0.0540,0.0]
const NC_DGSLP2 = Float32[0,0,0,0,0,0,0,0,0,0.80370,0,0]
const NC_DGSLQ2 = Float32[0,0,0,-0.83400,0,0,0,0,-0.83400,0,0,0]
const NC_DGEL2  = Float32[0,0,0,0,0,0,0,0,-0.00700,0,0,0]
const NC_DGSITE = Float32[0.47932,0.01401,0.56356,0.47360,0.20189,0.01200,0.32093,0.00659,0.00734,1.10842,0.00659,0.0]
# DEFAULT-branch DGCON slope/aspect terms (nc/dgf.f DATA DGSASP/DGCASP/DGSLOP/DGSLSQ) — dgf.f:478-485.
const NC_DGSASP = Float32[-0.02884,0.0,-0.040708,-0.01560,-0.10656,0.0,-0.11954,-0.03587,0.0,0.0,-0.03587,0.0]
const NC_DGCASP = Float32[-0.14319,0.0,-0.16836,-0.15630,-0.19174,0.0,0.08632,-0.19935,0.0,0.0,-0.19935,0.0]
const NC_DGSLOP = Float32[0.63500,0.0,0.46468,0.58937,-1.29627,0.0,0.85815,0.73530,0.0,0.0,0.73530,0.0]
const NC_DGSLSQ = Float32[-1.09400,0.0,-0.87145,-1.05045,0.87335,0.0,-1.17209,-0.99561,0.0,0.0,-0.99561,0.0]

# NC bark BRATIO: eqtype 1 DBT=a+b·D→(D−DBT)/D ; 2 DIB=a+b·D→DIB/D ; 3 DIB=a·D^b→a·D^(b−1) [POWER].
@inline function nc_bratio(a::Float32, b::Float32, eqtype::Int, d::Float32)
    d <= 0f0 && return 0.99f0
    if eqtype == 1
        return (d - (a + b * d)) / d
    elseif eqtype == 2
        return (a + b * d) / d
    else
        return (a * d^b) / d
    end
end

"NC DGCONS: per-species DGCON (nc/dgf.f ENTRY DGCONS), 3 branches. Stored in c.dg_const[sp]."
function nc_dgcons!(s::StandState)
    c = s.calib; p = s.plot
    ifor = Int(p.forest_idx); (ifor < 1 || ifor > 7) && (ifor = 1)
    elev = p.elevation; slope = p.slope
    @inbounds for sp in 1:12
        si = p.sp_site_index[sp]
        if sp == 12                                   # redwood
            dgcon = -3.502444f0 + 0.415435f0 * log(max(si, 1f0))
        elseif sp == 2 || sp == 6 || sp == 9          # SP/IC/RF: site form
            dgcon = NC_DGLAT2_5[sp] + NC_DGEL2[sp] * elev + NC_DGSLP2[sp] * slope +
                    NC_DGSLQ2[sp] * slope * slope + NC_DGSITE[sp] * si
        else                                          # default: DGFOR + elev²/slope/aspect + DGSITE·ln(SITEAR(3))
            asp = p.aspect; si3 = max(p.sp_site_index[3], 1f0)   # nc/dgf.f:478-485; SITEAR(3)=DF site, ALL sp
            dgcon = NC_DGFOR[sp, NC_MAPLOC[sp, ifor]] +
                    NC_DGEL2[sp] * elev * elev +
                    (NC_DGSASP[sp] * sin(asp) + NC_DGCASP[sp] * cos(asp) + NC_DGSLOP[sp]) * slope +
                    NC_DGSLSQ[sp] * slope * slope +
                    NC_DGSITE[sp] * log(si3)
        end
        c.dg_const[sp] = dgcon
        c.bark_a[sp] = 0f0; c.bark_b[sp] = 0.9f0       # NC uses nc_bratio directly in dgf! (not the linear cache)
    end
    return s
end

"NC dgf! — per-tree WK2 = LN(DDS) (nc/dgf.f main body; 3 branches, 5-yr TDDS/2)."
function dgf!(s::StandState, ::Klamath)
    p, t, c = s.plot, s.trees, s.calib
    dens = s.density; sd = s.coef.species
    wk2 = view(s.scratch.wk, 2, :)
    ba = p.basal_area; avh = p.avg_height
    relden = p.relative_density
    alrd = relden > 0f0 ? log(relden) : 0f0
    alba = ba > 0f0 ? log(ba) : 0f0
    slope = p.slope; asp = p.aspect
    ifor = Int(p.forest_idx); (ifor < 1 || ifor > 7) && (ifor = 1)
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        sp = Int(t.species[i])
        pt_i = Int(t.plot_id[i])
        pccf = (1 <= pt_i <= length(dens.point_ccf)) ? dens.point_ccf[pt_i] : 0f0
        ptba = (1 <= pt_i <= length(dens.point_ba)) ? dens.point_ba[pt_i] : ba
        icr = Float32(t.crown_pct[i]); cr = icr * 0.01f0
        pctfrac = 1f0 - t.crown_ratio[i] / 100f0       # 1 − PCT/100 (crown_ratio field holds PCT)
        ald = log(d)
        bal = pctfrac * ba
        pbal = pctfrac * ptba; pbal <= 0f0 && (pbal = bal)
        crid = d < 2f0 ? 1.8f0 : ((icr * icr) / log(d + 1f0)) / 1000f0
        alpba = ptba > 0f0 ? log(ptba) : 0f0
        cor = c.dg_cor[sp]                             # COR additive calib (0 baseline)
        if sp == 12                                    # REDWOOD
            si = p.sp_site_index[sp]
            prd = 0f0                                   # PRD point-Zeide (TODO precompute; 0 baseline)
            conspp = c.dg_const[sp]
            dglt = exp(conspp + 0.185911f0 * log(d) - 0.000073f0 * d * d - 0.001796f0 * pbal -
                       0.42078f0 * prd + 0.589318f0 * log(cr * 100f0) - 0.000926f0 * slope * 100f0 -
                       0.002203f0 * (slope * 100f0) * cos(asp))
            brat = nc_bratio(sd[:bark1][sp], sd[:bark2][sp], Int(sd[:bark_imap][sp]), d)
            t1 = d * brat; t2 = (d + dglt) * brat
            dds = log(t2 * t2 - t1 * t1) + cor          # + LN(COR2)=0 baseline
            dds = log(exp(dds) / 2f0)                    # redwood: TDDS/2 (nc/dgf.f:393-394)
        elseif sp == 2 || sp == 6 || sp == 9           # SP/IC/RF (set-2)
            conspp = c.dg_const[sp] + cor + NC_DGCCFA[sp] * alrd + NC_DGBA[sp] * alba
            dds = conspp + NC_DGLD2[sp] * ald + NC_DGDSQ2[sp] * d * d / 1000f0 +
                  NC_DGCR2[sp] * crid + NC_DGDBA2[sp] * pbal / log(d + 1f0) / 100f0 +
                  NC_DGBA2[sp] * alpba
            dds < -8.52f0 && (dds = -8.52f0)
            dds = log(exp(dds) / 2f0)                    # sp2/6/9: TDDS/2 (nc/dgf.f:410-411)
        else                                           # DEFAULT — NO TDDS/2 (nc/dgf.f has none here)
            dgdsq = NC_DGDS[sp, NC_MAPDSQ[sp, ifor]]
            conspp = c.dg_const[sp] + cor + NC_DGCCFA[sp] * alrd + NC_DGBA[sp] * alba
            hoavh = avh > 0f0 ? min(t.height[i] / avh, 1.5f0) : 1f0
            dds = conspp + NC_DGLD[sp] * ald + cr * (NC_DGCR[sp] + cr * NC_DGCRSQ[sp]) +
                  dgdsq * d * d + NC_DGDBAL[sp] * bal / log(d + 1f0) +
                  NC_DGPCCF[sp] * pccf + NC_DGHAH[sp] * hoavh
            sp == 4 && (dds -= 0.15032f0)
        end
        dds < -9.21f0 && (dds = -9.21f0)                # shared final clamp (nc/dgf.f:429)
        wk2[i] = dds
    end
    return s
end
