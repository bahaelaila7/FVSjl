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

# ws/htgf.f ENTRY HTCONS — HTCON(ISPC) site intercept (CASE DEFAULT species only; 0 for surrogate species).
function ws_htcons!(s::StandState)
    c = s.calib; p = s.plot
    elev = p.elevation; slope = p.slope
    itlat = round(Int, p.latitude)
    ilat = itlat <= 35 ? 1 : itlat == 36 ? 2 : itlat == 37 ? 3 : itlat == 38 ? 4 : 5
    @inbounds for isp in 1:43
        if isp in WS_HTG_SURR
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
    alba = ba > 0f0 ? log(ba) : 0f0
    ctl = s.control
    cor2on = ctl.htg_cor2_on
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0f0 && continue
        isp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]; dg = t.diam_growth[i]
        icr = Float32(t.crown_pct[i])
        htcon = c.htg_cor[isp]
        # XHT2 = HCOR2 (READCORH), default 1 (rides in multiplicatively — NOT into HTCON since non-log-linear).
        xht2 = (cor2on && ctl.htg_cor2[isp] > 0f0) ? ctl.htg_cor2[isp] : 1f0
        xht = 1f0                                                 # XHMULT growth multiplier (no MULTS kw)

        if isp in WS_HTG_SURR
            # chunk 4b: GB(21)/MC(41)/CA-surrogate/RW-GS(4,23) potential-height branches (need ws_findag).
            # Absent from wst01. Faithful minimum: HTG=0.1·SCALE·XHT·XHT2 (ws/htgf.f floors HTG≥0.1).
            htg = 0.1f0 * scale * xht * xht2
            t.ht_growth[i] = htg
            continue
        end

        # ── WS-native CASE DEFAULT linear regression (ws/htgf.f:822-831) ──
        pct = t.crown_ratio[i]
        bal = ((100f0 - pct) / 100f0) * ba
        bal <= 0f0 && (bal = 0.001f0)
        bai = (d + dg) * (d + dg) - d * d
        htg = htcon + WS_HGDG2[isp] * dg + WS_HGRDG2[isp] * sqrt(dg) + WS_HGBA2[isp] * ba +
              WS_HGBAI2[isp] * bai + WS_HGLBA2[isp] * alba + WS_HGBLT2[isp] * bal +
              WS_HGBAD2[isp] * bal / d + WS_HGCR2[isp] * icr + WS_HGDSQ[isp] * d * d
        htg *= _ws_htg_mult(isp)
        htg < 0.5f0 && (htg = 0.5f0)                              # temporary DSQ-neg trap (Dixon 8-18-93)
        (dg < 1f0 && d > 30f0) && (htg *= dg)                     # small-DG large-tree damp
        if isp == 42 || isp == 6                                  # JP/OS BAI-based MAXHTG cap
            maxhtg = -2.16f0 + 4.22f0 * log(bai)
            htg > maxhtg && (htg = maxhtg)
        end
        # HTMAX cap (oaks 28:33,40,43 use (D+DG+1)², else (D+DG+1)); floor 0.1.
        htmax = (isp in (28,29,30,31,32,33,40,43)) ?
            exp(WS_MXHTG1[isp] + WS_MXHTG2[isp] / ((d + dg + 1f0)^2)) + 4.5f0 :
            exp(WS_MXHTG1[isp] + WS_MXHTG2[isp] / (d + dg + 1f0)) + 4.5f0
        (h + htg) > htmax && (htg = htmax - h)
        htg < 0.1f0 && (htg = 0.1f0)
        htg = htg * scale * xht * xht2                            # × SCALE·XHT·XHT2 (MISHGF=1)
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
