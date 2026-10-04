# =============================================================================
# height_growth.jl (westsierra) — WS large-tree HEIGHT GROWTH (ws/htgf.f). Chunk 4a.
#
# WS's htgf.f has SEVERAL branches, but the WS-NATIVE main species (the outer CASE DEFAULT — SP/DF/WF/PP/IC/
# JP/RF/WP/BD/RW-no.../oaks…) use a SIMPLE ADDITIVE LINEAR regression (NOT the SO/CA potential-height/HTCALC
# model). HTCON(ISPC) is the site intercept assembled in ENTRY HTCONS; the large-tree HTG is explicitly NOT
# log-linear (ws/htgf.f:1121) so NO exp() and NO HCOR calibration into HTCON — HCOR2 rides in multiplicatively
# as XHT2 at the end.
#
#   HTG = HTCON + HGDG2·DG + HGRDG2·√DG + HGBA2·BA + HGBAI2·((D+DG)²−D²) + HGLBA2·ln(BA)
#         + HGBLT2·BAL + HGBAD2·BAL/D + HGCR2·ICR + HGDSQ·D²                      (BAL=((100−PCT)/100)·BA)
#   × species multiplier {42→1.5, 2/22→1.25, 3/13→1.5, 5→1.2, 6→1.5, 7→1.2, 8/18→1.15, 34:39→0.75}
#   floor 0.5 ; DG<1 & D>30 → ·DG ; JP/OS BAI cap ; HTMAX cap ; × SCALE·XHT·XHT2 ; × MISHGF(=1) ; SIZCAP.
#   HTCON = HGLAT2(ILAT,ISPC) + HGEL2·ELEV + HGELQ2·ELEV² + HGSL2·SLOPE + HGSI·SITEAR  (HTCONS).
#
# ⚠ CHUNK 4a scope: this implements the WS-native CASE DEFAULT linear branch + HTCONS (wst01 = SP1/DF2/WF3/RF7,
#   all CASE DEFAULT). The GB(21) even/uneven Alexander blend, MC(41) Curtis, the CA-surrogate FINDAG/HTCALC
#   Ritchie-Hann branch (9:10,12,14:17,19:20,25:27) and RW/GS(4,23) Castle-2019 LTHTG are STUBBED (need
#   ws_findag; absent from wst01) = chunk 4b. MEASURED bit-exact vs FVSws_g16 htgf per-tree HTG on wst01.
# =============================================================================

# ── WS per-species HTG coefficients (ws/htgf.f DATA), length 43, index = ISPC.
const WS_HGDG2 = Float32[
  0,0,0,0,0, 2.9620,0,3.7974,0,0, 0,0,0,0,0, 0,0,3.7974,0,0,
  0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,2.9620,0]
const WS_HGRDG2 = Float32[
  10.1890,15.3820,8.7587,10.1890,6.5757, 0,8.7587,0,0,0,
  10.1890,0,8.7587,0,0, 0,0,0,0,0, 0,15.3820,10.1890,10.1890,0,
  0,0,6.5757,6.5757,6.5757, 6.5757,6.5757,6.5757,6.5757,6.5757,
  6.5757,6.5757,6.5757,6.5757,6.5757, 0,0,6.5757]
const WS_HGBAI2 = Float32[
  0,-0.0694,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
  0,-0.0694,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0]
const WS_HGBA2 = Float32[
  0,-0.0285,0,0,0.0080, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
  0,-0.0285,0,0,0, 0,0,0.0080,0.0080,0.0080, 0.0080,0.0080,0.0080,0.0080,0.0080,
  0.0080,0.0080,0.0080,0.0080,0.0080, 0,0,0.0080]
const WS_HGLBA2 = Float32[
  1.2845,6.1190,1.2324,1.2845,0, 0,1.2324,1.4843,0,0,
  1.2845,0,1.2324,0,0, 0,0,1.4843,0,0, 0,6.1190,1.2845,1.2845,0,
  0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0]
const WS_HGBLT2 = Float32[
  0,0,0,0,-0.0079, 0,0,-0.0142,0,0, 0,0,0,0,0, 0,0,-0.0142,0,0,
  0,0,0,0,0, 0,0,-0.0079,-0.0079,-0.0079, -0.0079,-0.0079,-0.0079,-0.0079,-0.0079,
  -0.0079,-0.0079,-0.0079,-0.0079,-0.0079, 0,0,-0.0079]
const WS_HGBAD2 = Float32[
  0,0,-0.0399,0,0, 0,-0.0399,0,0,0, 0,0,-0.0399,0,0, 0,0,0,0,0,
  0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0]
const WS_HGCR2 = Float32[
  0,0,0.0363,0,0, 0,0.0363,0,0,0, 0,0,0.0363,0,0, 0,0,0,0,0,
  0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0]
const WS_HGDSQ = Float32[
  -0.0098,0,-0.0098,-0.0098,0, -0.0082,-0.0098,-0.0086,0,0,
  -0.0098,0,-0.0098,0,0, 0,0,-0.0086,0,0, 0,0,-0.0098,-0.0098,0,
  0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,-0.0082,0]
const WS_MXHTG1 = Float32[
  5.4918,5.4191,5.4283,5.4918,5.2011, 5.4916,5.4283,5.4916,0,0,
  5.4918,0,5.4283,0,0, 0,0,5.4916,0,0, 0,5.4191,5.4918,5.4918,0,
  0,0,4.6570,4.6570,4.6570, 4.6570,4.6570,4.6570,5.2011,5.2011,
  5.2011,5.2011,5.2011,5.2011,4.6570, 0,5.4916,4.6570]
const WS_MXHTG2 = Float32[
  -12.6438,-8.8274,-9.1641,-12.6438,-7.7610, -9.5992,-9.1641,-9.5992,0,0,
  -12.6438,0,-9.1641,0,0, 0,0,-9.5992,0,0, 0,-8.8274,-12.6438,-12.6438,0,
  0,0,-21.9333,-21.9333,-21.9333, -21.9333,-21.9333,-21.9333,-7.7610,-7.7610,
  -7.7610,-7.7610,-7.7610,-7.7610,-21.9333, 0,-9.5992,-21.9333]

# HTCONS site-intercept coefficients (ws/htgf.f ENTRY HTCONS DATA).
const WS_HGEL2 = Float32[
  0,0,-0.0453,0,0, 0,-0.0453,0,0,0, 0,0,-0.0453,0,0, 0,0,0,0,0,
  0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0]
const WS_HGELQ2 = Float32[
  0,0,0,0,-0.0007, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0,
  0,0,0,0,0, 0,0,-0.0007,-0.0007,-0.0007, -0.0007,-0.0007,-0.0007,-0.0007,-0.0007,
  -0.0007,-0.0007,-0.0007,-0.0007,-0.0007, 0,0,-0.0007]
const WS_HGSL2 = Float32[
  0,0,3.2180,0,0, 0,3.2180,0,0,0, 0,0,3.2180,0,0, 0,0,0,0,0,
  0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0]
const WS_HGSI = Float32[
  0.0881,0.0974,0.0713,0.0881,0.0643, 0.0687,0.0713,0.0598,0,0,
  0.0881,0,0.0713,0,0, 0,0,0.0598,0,0, 0,0.0974,0.0881,0.0881,0,
  0,0,0.0643,0.0643,0.0643, 0.0643,0.0643,0.0643,0.0643,0.0643,
  0.0643,0.0643,0.0643,0.0643,0.0643, 0,0.0687,0.0643]

# HGLAT2[isp, ilat] (ws/htgf.f DATA HGLAT2(5,MAXSP) → jl [isp,ilat]; species intercept by latitude class 1..5).
const WS_HGLAT2 = Float32[
  -14.053 -14.053 -14.053 -14.053 -14.053;   #  1 SP
  -39.027 -39.027 -39.027 -39.027 -39.027;   #  2 DF
  -12.246 -10.685 -10.685 -10.685 -10.685;   #  3 WF
  -14.053 -14.053 -14.053 -14.053 -14.053;   #  4 GS
   -3.831   1.711  -3.831  -2.898  -3.831;   #  5 IC
   -0.268  -0.268  -0.268   3.097  -0.268;   #  6 JP
  -12.246 -10.685 -10.685 -10.685 -10.685;   #  7 RF
   -5.554  -5.554  -5.554  -4.282  -5.554;   #  8 PP
    0.0     0.0     0.0     0.0     0.0;      #  9 LP
    0.0     0.0     0.0     0.0     0.0;      # 10 WB
  -14.053 -14.053 -14.053 -14.053 -14.053;   # 11 WP
    0.0     0.0     0.0     0.0     0.0;      # 12 PM
  -12.246 -10.685 -10.685 -10.685 -10.685;   # 13 SF
    0.0     0.0     0.0     0.0     0.0;      # 14 KP
    0.0     0.0     0.0     0.0     0.0;      # 15 FP
    0.0     0.0     0.0     0.0     0.0;      # 16 CP
    0.0     0.0     0.0     0.0     0.0;      # 17 LM
   -5.554  -5.554  -5.554  -4.282  -5.554;   # 18 MP
    0.0     0.0     0.0     0.0     0.0;      # 19 GP
    0.0     0.0     0.0     0.0     0.0;      # 20 WE
    0.0     0.0     0.0     0.0     0.0;      # 21 GB
  -39.027 -39.027 -39.027 -39.027 -39.027;   # 22 BD
  -14.053 -14.053 -14.053 -14.053 -14.053;   # 23 RW
  -14.053 -14.053 -14.053 -14.053 -14.053;   # 24 MH
    0.0     0.0     0.0     0.0     0.0;      # 25 WJ
    0.0     0.0     0.0     0.0     0.0;      # 26 UJ
    0.0     0.0     0.0     0.0     0.0;      # 27 CJ
   -3.831   1.711  -3.831  -2.898  -3.831;   # 28 LO
   -3.831   1.711  -3.831  -2.898  -3.831;   # 29 CY
   -3.831   1.711  -3.831  -2.898  -3.831;   # 30 BL
   -3.831   1.711  -3.831  -2.898  -3.831;   # 31 BO
   -3.831   1.711  -3.831  -2.898  -3.831;   # 32 VO
   -3.831   1.711  -3.831  -2.898  -3.831;   # 33 IO
   -3.831   1.711  -3.831  -2.898  -3.831;   # 34 TO
   -3.831   1.711  -3.831  -2.898  -3.831;   # 35 GC
   -3.831   1.711  -3.831  -2.898  -3.831;   # 36 AS
   -3.831   1.711  -3.831  -2.898  -3.831;   # 37 CL
   -3.831   1.711  -3.831  -2.898  -3.831;   # 38 MA
   -3.831   1.711  -3.831  -2.898  -3.831;   # 39 DG
   -3.831   1.711  -3.831  -2.898  -3.831;   # 40 BM
    0.0     0.0     0.0     0.0     0.0;      # 41 MC
   -0.268  -0.268  -0.268   3.097  -0.268;   # 42 OS
   -3.831   1.711  -3.831  -2.898  -3.831]   # 43 OH

# Species multiplier applied to the raw linear HTG (ws/htgf.f:836-853).
@inline function _ws_htg_mult(isp::Int)
    isp == 42 && return 1.5f0
    (isp == 2 || isp == 22) && return 1.25f0
    (isp == 3 || isp == 13) && return 1.5f0
    isp == 5 && return 1.2f0
    isp == 6 && return 1.5f0
    isp == 7 && return 1.2f0
    (isp == 8 || isp == 18) && return 1.15f0
    (34 <= isp <= 39) && return 0.75f0
    return 1.0f0
end

# WS-native CASE DEFAULT species = everything NOT in the surrogate branches.
const WS_HTG_SURR = Set{Int}([21, 41, 4, 9, 10, 12, 14, 15, 16, 17, 19, 20, 23, 25, 26, 27])
# ENTRY HTCONS zeroes HTCON only for GB(21), MC(41) and the CA-surrogate set; RW/GS (4,23) fall into its CASE
# DEFAULT (HTCON computed), although htgf routes them to the CA branch and never reads it.
const WS_HTCON_ZERO = Set{Int}([21, 41, 9, 10, 12, 14, 15, 16, 17, 19, 20, 25, 26, 27])
# htgf CASE(4,9:10,12,14:17,19:20,23,25:27) — species using the CA variant's equations.
const WS_HTG_CA = Set{Int}([4, 9, 10, 12, 14, 15, 16, 17, 19, 20, 23, 25, 26, 27])

# ws/findag.f DATA AGEMAX(43) / HTMAX(43).
const WS_FINDAG_AGEMAX = Float32[
    0, 0, 0, 0, 0, 0, 0, 0, 400, 400,   0, 400, 0, 400, 400, 400, 400, 0, 400, 400,
    210, 0, 0, 0, 400, 400, 400, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   400, 0, 0]
const WS_FINDAG_HTMAX = Float32[i == 41 ? 20 : 0 for i in 1:43]

"""ws/findag.f — (SITAGE, SITHT, AGMAX, HTMAX1) for tree height `h` on species `ispc`'s site curve.
MC(41) at/above HTMAX1=20 ft jumps straight to SITAGE=AGMAX+(H−HTMAX1)/0.10. GB(21) solves the Alexander
curve in 5-yr steps from AP=10 (BAUTBA = BAU/BA; BAU is all-zero in WS — ws/grinit.f, no BADIST).
Everyone else steps AG by 2 on ws_htcalc until within TOLER=2 ft (or past H, or the curve flattens)."""
function ws_findag(ispc::Int, h::Float32, sindx::Float32, ifor::Int, d1::Float32, bautba::Float32)
    agmax = WS_FINDAG_AGEMAX[ispc]; htmax1 = WS_FINDAG_HTMAX[ispc]
    sitage = 0f0; sitht = 0f0
    if ispc == 41 && h >= htmax1
        return (agmax + (h - htmax1) / 0.10f0, h, agmax, htmax1)
    end
    if ispc == 21
        bautba < 0f0 && (bautba = 0f0)
        site = sindx; tol = 2f0; ap = 10f0
        local tage::Float32
        while true
            agetem = ap < 30f0 ? 30f0 : ap
            hh = (2.75780f0 * fpow(site, 0.83312f0)) *
                 fpow(1f0 - fexp(-0.015701f0 * agetem), 22.71944f0 * fpow(site, -0.63557f0)) + 4.5f0
            ratio = 1f0 - bautba; ratio < 0.728f0 && (ratio = 0.728f0)
            hh *= ratio
            if abs(hh - h) < tol || hh > h
                tsite = site < 20f0 ? 20f0 : site
                tage = agetem + 4.5f0 / (-0.22f0 + 0.0155f0 * tsite)
                break
            end
            ap += 5f0
            if ap > agmax
                tage = agmax; break
            end
        end
        sitage = tage; sitage <= 0f0 && (sitage = 1f0)
        return (sitage, sitht, agmax, htmax1)
    end
    toler = 2f0; ag = 2f0; incrng = 0; hguess = 0f0
    while true
        oldhg = hguess
        hguess = ws_htcalc(ifor, sindx, ispc, ag)
        if hguess >= 1f0
            if abs(hguess - h) <= toler || h < hguess
                return (ag, hguess, agmax, htmax1)
            end
            dd = hguess - oldhg
            (oldhg != 0f0 && dd >= 0.05f0) && (incrng = 1)
            (incrng == 1 && dd < 0.05f0) && return (ag, hguess, agmax, htmax1)
        end
        ag += 2f0
        ag > agmax && return (agmax, h, agmax, htmax1)
    end
end

"Alexander (1967) RM-32 breast-height-age site curve used by ws/htgf.f CASE(21) (base-age-100 ES/AF)."
@inline _ws_alexander(tsite::Float32, age::Float32)::Float32 =
    (2.75780f0 * fpow(tsite, 0.83312f0)) * fpow(1f0 - fexp(-0.015701f0 * age), 22.71944f0 * fpow(tsite, -0.63557f0)) + 4.5f0

"ws/htgf.f SIZCAP (col 4) compliance on one record: HT+HTG ≤ cap, floor 0.1."
@inline function _ws_sizcap(htg::Float32, h::Float32, cap::Float32)::Float32
    if cap > 0f0 && h + htg > cap
        htg = cap - h; htg < 0.1f0 && (htg = 0.1f0)
    end
    return htg
end

# ws/htgf.f ENTRY HTCONS — HTCON(ISPC) site intercept (CASE DEFAULT species only; 0 for surrogate species).
function ws_htcons!(s::StandState)
    c = s.calib; p = s.plot
    elev = p.elevation; slope = p.slope
    itlat = round(Int, p.latitude)
    ilat = itlat <= 35 ? 1 : itlat == 36 ? 2 : itlat == 37 ? 3 : itlat == 38 ? 4 : 5
    @inbounds for isp in 1:43
        if isp in WS_HTCON_ZERO
            c.htg_cor[isp] = 0f0
        else
            sitear = p.sp_site_index[isp]
            c.htg_cor[isp] = WS_HGLAT2[isp, ilat] + WS_HGEL2[isp] * elev +
                             WS_HGELQ2[isp] * elev * elev + WS_HGSL2[isp] * slope +
                             WS_HGSI[isp] * sitear
        end
    end
    return s
end

function height_growth!(s::StandState, ::WestSierra; scale::Float32 = 1.0f0)
    p, t, c = s.plot, s.trees, s.calib
    ba = p.basal_area
    alba = ba > 0f0 ? flog(ba) : 0f0
    ctl = s.control
    cor2on = ctl.htg_cor2_on
    sd = s.coef.species; ifor = Int(p.forest_idx); avh = p.avg_height; pccf = s.density.point_ccf
    dgsd = ctl.dg_sd
    # ws/htgf.f is SPECIES-MAJOR (DO 40 ISPC … I=IND1(I3)); GB's ZZRAN draws must follow that order.
    @inbounds for i in species_major_order(s)
        t.ht_growth[i] = 0f0; t.temhtg[i] = -1f0
        t.tpa[i] <= 0f0 && continue
        isp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]; dg = t.diam_growth[i]
        icr = Float32(t.crown_pct[i])
        htcon = c.htg_cor[isp]
        # XHT2 = HCOR2 (READCORH), default 1 (rides in multiplicatively — NOT into HTCON since non-log-linear).
        xht2 = (cor2on && ctl.htg_cor2[isp] > 0f0) ? ctl.htg_cor2[isp] : 1f0
        xht = 1f0                                                 # XHMULT growth multiplier (no MULTS kw)

        sindx = p.sp_site_index[isp]
        cap = ctl.sp_size_cap[isp, 4]
        if isp == 21
            # ── htgf CASE(21) GB — UT-variant Alexander even-aged curve (value later overwritten by REGENT, whose
            # XMIN/XMAX = 99/199 cover every GB; the ZZRAN draw still consumes the stream). BAU≡0 and AGERNG≡0 in
            # WS (grinit, no BADIST) ⇒ RATIO = 1 and the uneven-aged blend never fires. ISTAGF≡0 (no keyword).
            htg = 0f0
            if !(d < 0.5f0 || h <= 4.5f0)
                bark = ws_bratio(sd, isp, d)
                tsite = sindx < 20f0 ? 20f0 : sindx
                ap = t.birth_age[i] - (4.5f0 / (-0.22f0 + 0.0155f0 * tsite)); ap < 1f0 && (ap = 1f0)
                agetem = ap < 30f0 ? 30f0 : ap                     # TSITE is reset to the raw SITEAR here
                hhe1 = _ws_alexander(sindx, agetem); ap < agetem && (hhe1 = ((hhe1 - 4.5f0) / agetem) * ap + 4.5f0)
                agefut = ap + 10f0
                agetem = agefut < 30f0 ? 30f0 : agefut
                hhe2 = _ws_alexander(sindx, agetem)
                agefut < agetem && (hhe2 = ((hhe2 - 4.5f0) / agetem) * agefut + 4.5f0)
                htg = hhe2 - hhe1                                  # ×RATIO(=1)·ADJUST(=1)
                zzran = 0f0
                if dgsd > 0f0
                    while true
                        zzran = bachlo(s.rng, 0f0, 1f0)
                        (zzran > dgsd || zzran < -dgsd) || break
                    end
                end
                htg += zzran * 0.1f0
                htg < 0.1f0 && (htg = 0.1f0)
            end
            htg = htg * scale * xht * xht2                          # label 201
            htg *= mis_hg_mult(s, i)                                # ws/htgf.f:460 MISHGF
            t.ht_growth[i] = _ws_sizcap(htg, h, cap)
            continue
        elseif isp == 41
            # ── htgf CASE(41) MC — SO-variant Curtis potential ×(0.25·Hoerl CR + 0.75·Chapman-Richards RELHT).
            # FAITHFUL QUIRK: the normal path ends at label 161 without SCALE/XHT/SIZCAP (only the H≥HTMAX path
            # applies SCALE·XHT·XHT2), and HTMAX(41)=20 ft caps H+HTG.
            sitage, sitht, agmax, htmax = ws_findag(41, h, sindx, ifor, d, 0f0)
            if h >= htmax
                t.ht_growth[i] = 0.1f0 * scale * xht * xht2
                continue
            end
            local pothtg::Float32
            if sitage > agmax
                pothtg = 0.10f0
            else
                agp10 = sitage + 10f0
                hguess = (sindx - 4.5f0) / (0.6192f0 - 5.3394f0 / (sindx - 4.5f0) +
                         240.29f0 * fpow(agp10, -1.4f0) + (3368.9f0 / (sindx - 4.5f0)) * fpow(agp10, -1.4f0))
                pothtg = (hguess + 4.5f0) - sitht
            end
            relht = avh > 0f0 ? h / avh : 0f0; relht > 1.5f0 && (relht = 1.5f0)
            xcr = Float32(t.crown_pct[i]) / 100f0
            hgmdcr = (100f0 * fpow(xcr, 3f0)) * fexp(-5f0 * xcr); hgmdcr > 1f0 && (hgmdcr = 1f0)
            fctrkx = fpow(1f0 / 0.10f0, 1.10f0 - 1f0) - 1f0
            fctrrb = -1f0 * (15f0 / (1f0 - (-1.45f0)))
            fctrxb = fpow(relht, 1f0 - (-1.45f0)) - fpow(0f0, 1f0 - (-1.45f0))
            fctrm = -1f0 / (1.10f0 - 1f0)
            hgmdrh = 1f0 * fpow(1f0 + fctrkx * fexp(fctrrb * fctrxb), fctrm)
            htgmod = 0.25f0 * hgmdcr + 0.75f0 * hgmdrh
            htgmod >= 2f0 && (htgmod = 2f0); htgmod <= 0f0 && (htgmod = 0.1f0)
            htg = pothtg * htgmod
            (h + htg > htmax) && (htg = htmax - h)
            htg < 0.1f0 && (htg = 0.1f0)
            t.ht_growth[i] = htg
            continue
        elseif isp in WS_HTG_CA
            # ── htgf CASE(4,9:10,12,14:17,19:20,23,25:27) — CA-variant equations.
            local htg::Float32
            if isp == 4 || isp == 23
                # RW/GS Castle LN(HI) on the 10-yr outside-bark DG, bounded between 217 and 380 ft.
                brat = ws_bratio(sd, isp, d)
                dg10 = dg / brat; h < 4.5f0 && (dg10 = 0.1f0)
                lthtg = fexp(1.412947f0 - 0.000204f0 * d * d + 0.31971f0 * flog(d) + 0.394005f0 * flog(sindx) +
                            0.399888f0 * flog(dg10) - 0.451708f0 * flog(h))
                hgbnd = (h >= 217f0 && h < 380f0) ? max(1f0 - ((h - 217f0) / (380f0 - 217f0)), 0.1f0) :
                        (h < 217f0 ? 1f0 : 0.1f0)
                htg = lthtg * hgbnd
            else
                # FINDAG → Dunning/Levitan potential (ws_htcalc) × Ritchie-Hann XMOD.
                sitage, sitht, agmax, _ = ws_findag(isp, h, sindx, ifor, d, 0f0)
                pothtg = sitage > agmax ? 0.10f0 : ws_htcalc(ifor, sindx, isp, sitage + 10f0) - sitht
                cratio = Float32(t.crown_pct[i]) / 100f0
                relht = avh > 0f0 ? h / avh : 1f0              # H/AVH (AVH=0 ⇒ +Inf ⇒ capped to 1)
                relht > 1f0 && (relht = 1f0)
                pt = Int(t.plot_id[i])
                ((1 <= pt <= length(pccf)) ? pccf[pt] : 0f0) < 100f0 && (relht = 1f0)
                crmod = 1f0 - fexp(-4.26558f0 * cratio)
                rhmod = fexp(2.54119f0 * (fpow(relht, 0.250537f0) - 1f0))
                htg = pothtg * (1.016605f0 * crmod * rhmod)
            end
            htg < 0.1f0 && (htg = 0.1f0)
            htg = scale * xht * htg * xht2
            htg *= mis_hg_mult(s, i)                                # ws/htgf.f:784 MISHGF
            t.ht_growth[i] = _ws_sizcap(htg, h, cap)
            continue
        end

        # ── WS-native CASE DEFAULT linear regression (ws/htgf.f:822-831) ──
        pct = t.crown_ratio[i]
        bal = ((100f0 - pct) / 100f0) * ba
        bal <= 0f0 && (bal = 0.001f0)
        bai = fpow(d + dg, 2.0f0) - d * d                        # ws/htgf.f:828,873 ((DBH+DG)**2.0) = powf
        htg = htcon + WS_HGDG2[isp] * dg + WS_HGRDG2[isp] * sqrt(dg) + WS_HGBA2[isp] * ba +
              WS_HGBAI2[isp] * bai + WS_HGLBA2[isp] * alba + WS_HGBLT2[isp] * bal +
              WS_HGBAD2[isp] * bal / d + WS_HGCR2[isp] * icr + WS_HGDSQ[isp] * d * d
        htg *= _ws_htg_mult(isp)
        htg < 0.5f0 && (htg = 0.5f0)                              # temporary DSQ-neg trap (Dixon 8-18-93)
        (dg < 1f0 && d > 30f0) && (htg *= dg)                     # small-DG large-tree damp
        if isp == 42 || isp == 6                                  # JP/OS BAI-based MAXHTG cap
            maxhtg = -2.16f0 + 4.22f0 * flog(bai)
            htg > maxhtg && (htg = maxhtg)
        end
        # HTMAX cap (oaks 28:33,40,43 use (D+DG+1)², else (D+DG+1)); floor 0.1.
        htmax = (isp in (28,29,30,31,32,33,40,43)) ?
            fexp(WS_MXHTG1[isp] + WS_MXHTG2[isp] / fpow(d + dg + 1f0, 2.0f0)) + 4.5f0 :   # ws/htgf.f:889 **2.0 = powf
            fexp(WS_MXHTG1[isp] + WS_MXHTG2[isp] / (d + dg + 1f0)) + 4.5f0
        (h + htg) > htmax && (htg = htmax - h)
        htg < 0.1f0 && (htg = 0.1f0)
        htg = htg * scale * xht * xht2                            # × SCALE·XHT·XHT2
        htg *= mis_hg_mult(s, i)                                  # ws/htgf.f:914 MISHGF (before TEMHTG)
        t.temhtg[i] = htg                                         # htgf.f:918 TEMHTG → the tripled copies (ws_triple_htg!)
        # SIZCAP: HT+HTG ≤ species size cap (col 4).
        sizcap = ctl.sp_size_cap[isp, 4]
        if sizcap > 0f0 && (h + htg) > sizcap
            htg = sizcap - h
            htg < 0.1f0 && (htg = 0.1f0)
        end
        t.ht_growth[i] = htg
    end
    return s
end

"""
    ws_triple_htg!(s, stash)

ws/htgf.f:921-956 (CASE DEFAULT, LTRIP): each tripled copy's large-tree HTG is the central's scaled TEMHTG
(pre-SIZCAP) times the copy's spread DG over the central's, HTG(ITFN)=TEMHTG·DG(ITFN)/DGI (DGI=MAX(DG(I),0.01)),
then re-capped at label 131 by HTMAX from the COPY's DG (EXP(MXHTG1+MXHTG2/(DBH+DG(ITFN)+1))+4.5, oaks squared;
floor 0.1) and SIZCAP. height_growth! only grew the central (copy_tree! handed both copies its HTG). The other
WS branches (GB, MC, CA-surrogate) give the copies TEMHTG = the central's value and are left to copy_tree!.
Deterministic (no draws). Runs after height_growth!, before small_tree_growth! (HTGF → REGENT).
"""
function ws_triple_htg!(s::StandState, stash)
    (stash === nothing || isempty(stash.dgU)) && return s
    t, ctl = s.trees, s.control
    @inbounds for i in 1:stash.nlive
        tem = t.temhtg[i]; tem < 0f0 && continue
        isp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        dgi = t.diam_growth[i]; dgi < 0.01f0 && (dgi = 0.01f0)
        cap = ctl.sp_size_cap[isp, 4]
        oak = isp in (28, 29, 30, 31, 32, 33, 40, 43)
        for l in 1:2
            dgc = l == 1 ? stash.dgU[i] : stash.dgL[i]
            hc = tem * dgc / dgi
            htmax = oak ? fexp(WS_MXHTG1[isp] + WS_MXHTG2[isp] / fpow(d + dgc + 1f0, 2.0f0)) + 4.5f0 :
                          fexp(WS_MXHTG1[isp] + WS_MXHTG2[isp] / (d + dgc + 1f0)) + 4.5f0
            (h + hc) > htmax && (hc = htmax - h)
            hc < 0.1f0 && (hc = 0.1f0)
            (cap > 0f0 && (h + hc) > cap) && (hc = max(cap - h, 0.1f0))
            l == 1 ? (stash.htgU[i] = hc) : (stash.htgL[i] = hc)
        end
        stash.htg_copy[i] = true
    end
    return s
end
