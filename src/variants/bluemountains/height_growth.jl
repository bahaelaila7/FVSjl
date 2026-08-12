# =============================================================================
# height_growth.jl (bluemountains) — BM large-tree height growth (bm/htgf.f + bm/findag.f). Chunk 4.
#
# Normal species (all except WJ/WB/LM/AS = 6,11,12,15): FINDAG(H)→SITAGE,SITHT,HTMAX,AGMAX; then
#   H≥HTMAX  → HTG = (AASM + BASM·RELSI)·XHT·SCALE·exp(HTCON),  RELSI=(SI−SLO)/(SHI−SLO)
#   SITAGE>AGMAX → POTHTG = 0.1  (PP/OS sp10/17: −1.31+0.05·SINDX)
#   else     → POTHTG = HTCALC(SINDX, sp, SITAGE+10) − SITHT
#   modifiers: HGMDCR=CRA·(CR)^CRB·exp(CRC·CR) [CRA100/CRB3/CRC−5, ≤1]; RELHT=H/AVH(≤1.5);
#     HGMDRH=RHK·(1+FCTRKX·exp(FCTRRB·FCTRXB))^FCTRM (gen. Chapman-Richards); HTGMOD=.25·HGMDCR+.75·HGMDRH ∈[.1,2]
#   HTG = POTHTG·HTGMOD; cap H+HTG≤HTMAX; HTG = SCALE·XHT·HTG·exp(HTCON).  HTCON = self-calib (htg_cor).
# WB/LM/AS (11,12,15) use Johnson's SBB (COF1/COF6) — deferred (not in bmt01). WJ (6) height from REGENT.
# =============================================================================

# bm/htgf.f DATA (species order WP WL DF GF MH WJ LP ES AF PP WB LM PY YC AS CW OS OH).
const BM_AASM  = Float32[1.5,2.0,0.4,2.1,0.0,0.0,1.5,1.5,1.5,1.3,0.0,0.0,0.0,0.0,0.0,0.0,1.3,0.0]
const BM_BASM  = Float32[0.003,0.0026,0.0080,0.0,0.0,0.0,0.0,0.0,0.0,0.002,0.0,0.0,0.0,0.0,0.0,0.0,0.002,0.0]
const BM_RHR   = Float32[15.0,12.0,15.0,20.0,20.0,0.0,12.0,16.0,16.0,13.0,0.0,0.0,20.0,16.0,0.0,12.0,13.0,15.0]
const BM_RHYXS = Float32[0.10,0.01,0.10,0.20,0.20,0.0,0.01,0.15,0.15,0.05,0.0,0.0,0.20,0.15,0.0,0.01,0.05,0.10]
const BM_RHM   = Float32[1.10,1.10,1.10,1.10,1.10,0.0,1.10,1.10,1.10,1.10,0.0,0.0,1.10,1.10,0.0,1.10,1.10,1.10]
const BM_RHB   = Float32[-1.45,-1.60,-1.45,-1.10,-1.10,0.0,-1.60,-1.20,-1.20,-1.60,0.0,0.0,-1.10,-1.20,0.0,-1.60,-1.60,-1.45]
const BM_HT_CRA = 100.0f0; const BM_HT_CRB = 3.0f0; const BM_HT_CRC = -5.0f0
const BM_HT_RHK = 1.0f0;   const BM_HT_RHXS = 0.0f0
# bm/findag.f DATA.
const BM_AGMAX = Float32[200,110,180,130,250,900,140,150,150,200,400,400,200,200,100,100,200,100]
const BM_AHMAX = Float32[2.3,12.86,-2.86,21.29,52.27,80.0,2.3,20.0,45.27,-5.0,85.0,85.0,50.0,100.0,75.0,125.0,-5.0,100.0]
const BM_BHMAX = Float32[2.39,1.32,1.54,1.24,1.14,0.0,1.75,1.1,1.24,1.30,0.0,0.0,0.0,0.0,0.0,0.0,1.30,0.0]

# bm/findag.f — invert the site height-age curve: age at which species sp reaches height h on site sindx.
# Returns (sitage, sitht, htmax, agmax).
function bm_findag(sp::Int, h::Float32, sindx::Float32)
    agmax = BM_AGMAX[sp]
    htmax = BM_AHMAX[sp] + BM_BHMAX[sp] * sindx
    sp == 5 && (htmax *= 3.281f0)
    if h >= htmax
        return (agmax + (h - htmax) / 0.10f0, h, htmax, agmax)
    end
    if sp == 15
        return ((h * 2.54f0 * 12.0f0 / 26.9825f0)^(1.0f0 / 1.1752f0), h, htmax, agmax)
    elseif sp == 6 || sp == 11 || sp == 12
        return (0f0, h, htmax, agmax)
    end
    ag = 2.0f0
    (sp == 10 || sp == 17) && (ag = 98.38f0 * exp(sindx * (-0.0422f0)) + 1.0f0)
    ag < 2.0f0 && (ag = 2.0f0)
    sp == 3 && (ag = 18.0f0)
    incrng = 0; hguess = 0f0
    while true
        oldhg = hguess
        hguess = bm_htcalc(sindx, sp, ag)
        if hguess < 1f0
            ag += 2f0
            if ag > agmax; return (agmax, h, htmax, agmax); end
            continue
        end
        if hguess <= 0f0
            ag += 2f0; continue
        end
        diff = abs(hguess - h)
        (diff <= 2.0f0 || h < hguess) && return (ag, hguess, htmax, agmax)
        d2 = hguess - oldhg
        (oldhg != 0f0 && d2 >= 0.05f0) && (incrng = 1)
        (incrng == 1 && d2 < 0.05f0) && return (ag, hguess, htmax, agmax)
        ag += 2f0
        if ag > agmax; return (agmax, h, htmax, agmax); end
    end
end

# bm/htgf.f Johnson-SBB height coeffs (Schreuder-Hafley SBB). COF1 = WB(11)/LM(12) (from TT/UT), COF6 = AS(15)
# (from UT). [9 coeffs × 3 crown classes], K from IICR=int(ICR/10+0.5): {1,2}→1 {3-7}→2 {8,9}→3. XI1=0.1, XI2=4.5.
# (Identical to EM _EM_COFLM/_EM_COFAS — the same shared Schreuder-Hafley tables.)
const _BM_HT_COF1 = Float32[
    37.0     45.0     45.0
    85.0    100.0     90.0
     1.77836  1.66674  1.64770
    -0.51147  0.25626  0.30546
     1.88795  1.45477  1.35015
     1.20654  1.11251  0.94823
     0.57697  0.67375  0.70453
     3.57635  2.17942  2.46480
     0.90283  0.88103  1.00316]
const _BM_HT_COF6 = Float32[
    30.0     30.0     35.0
    85.0     85.0     85.0
     2.00995  2.00995  1.80388
     0.03288  0.03288 -0.07682
     1.81059  1.81059  1.70032
     1.28612  1.28612  1.29148
     0.72051  0.72051  0.72343
     3.00551  3.00551  2.91519
     1.01433  1.01433  0.95244]

function height_growth!(s::StandState, ::BlueMountains; scale::Float32 = 1.0f0)
    p, t, c = s.plot, s.trees, s.calib
    site_lo = s.coef.species[:site_lo]; site_hi = s.coef.species[:site_hi]
    avh = p.avg_height
    # Johnson-SBB context (WB/LM/AS) + its young-tree accelerator (bm/htgf.f:365-390).
    relden = p.relative_density
    iage = Int(p.stand_age)
    iy1 = Int(s.control.cycle_year[1]); cur_year = current_cycle_year(s); icyc = Int(s.control.cycle)
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (d <= 0f0 || h <= 0f0) && continue
        sp == 6 && continue                                # WJ(6) height from REGENT
        if sp == 11 || sp == 12 || sp == 15                # WB/LM/AS — Johnson's SBB (bm/htgf.f CASE(11,12,15))
            htcon = exp(c.htg_cor[sp])
            iicr = trunc(Int, Float32(t.crown_pct[i]) / 10f0 + 0.5f0); iicr > 9 && (iicr = 9); iicr < 1 && (iicr = 1)
            k = iicr <= 2 ? 1 : (iicr <= 7 ? 2 : 3)
            cof = (sp == 11 || sp == 12) ? _BM_HT_COF1 : _BM_HT_COF6
            cof1 = cof[1,k]; cof2 = cof[2,k]; cof3 = cof[3,k]; cof4 = cof[4,k]; cof5 = cof[5,k]
            cof6 = cof[6,k]; cof7 = cof[7,k]; cof8 = cof[8,k]; cof9 = cof[9,k]
            # small / out-of-fitted-range ⇒ HTG=0.1 (bm/htgf.f:305,317-318); REGENT overrides small trees anyway.
            if h <= 4.5f0 || d <= 0.1f0 || (0.1f0 + cof1) <= d || (4.5f0 + cof2) <= h
                t.ht_growth[i] = 0.1f0 * scale * htcon; continue
            end
            y1 = (d - 0.1f0) / cof1; y2 = (h - 4.5f0) / cof2
            fby1 = log(y1 / (1f0 - y1)); fby2 = log(y2 / (1f0 - y2))
            z = (cof4 + cof6 * fby2 - cof7 * (cof3 + cof5 * fby1)) * (1f0 - cof7 * cof7)^(-0.5f0)
            # ZBIAS: AZBIAS=BZBIAS=0 for all BM SBB species (bm/htgf.f:145-147) ⇒ ZBIAS≡0, inert.
            if sp == 12 || sp == 15                         # LM/AS bias correction (NOT WB) — bm/htgf.f:357
                zadj = 0.1f0 - 0.10273f0 * z + 0.00273f0 * z * z
                zadj < 0f0 && (zadj = 0f0)
                z += zadj
            end
            # bm/htgf.f:365-390 young-small-tree HTG accelerator (Targhee). Cycle-1 only, young (10<IXAGE<40), D<9, Z<2.
            if iage != 0 && icyc == 0
                ixage = iage + cur_year - iy1
                if ixage < 40 && ixage > 10 && d < 9.0f0 && z <= 2.0f0
                    zadja = 0.3564f0 * t.diam_growth[i] * scale
                    closur = relden < 100.0f0 ? 1.0f0 : Float32(t.crown_ratio[i]) / 100.0f0
                    zadja *= closur
                    (iicr == 9 || iicr == 8) && (zadja *= 1.1f0)
                    z += zadja; z > 2.0f0 && (z = 2.0f0)
                end
            end
            bark = bm_bratio(s.coef.species, sp, d)
            dia = d + t.diam_growth[i] / bark
            if (0.1f0 + cof1) > dia
                psi = cof8 * ((dia - 0.1f0) / (0.1f0 + cof1 - dia))^cof9 *
                      exp(z * ((1f0 - cof7 * cof7)^0.5f0) / cof6)
                hnew = (psi / (1f0 + psi)) * cof2 + 4.5f0
                hnew < h && (hnew = h)
                htg = hnew - h; htg < 0.1f0 && (htg = 0.1f0)
                t.ht_growth[i] = htg * scale * htcon
            else
                t.ht_growth[i] = 0.1f0 * scale * htcon
            end
            continue
        end
        sindx = p.sp_site_index[sp]
        sitage, sitht, htmax, agmax = bm_findag(sp, h, sindx)
        htcon = exp(c.htg_cor[sp])
        if h >= htmax
            slo = site_lo[sp]; shi = site_hi[sp]
            si = sindx; si > shi && (si = shi); si <= slo && (si = slo + 0.5f0)
            relsi = (si - slo) / (shi - slo)
            htg = BM_AASM[sp] + BM_BASM[sp] * relsi
            htg < 0.1f0 && (htg = 0.1f0)
            t.ht_growth[i] = htg * scale * htcon
            continue
        end
        # potential height growth
        local pothtg::Float32
        if sitage > agmax
            pothtg = 0.1f0
            (sp == 10 || sp == 17) && (pothtg = max(-1.31f0 + 0.05f0 * sindx, 0.1f0))
        else
            pothtg = bm_htcalc(sindx, sp, sitage + 10f0) - sitht
            pothtg < 0.1f0 && (pothtg = 0.1f0)
        end
        # modifiers
        crf = Float32(t.crown_pct[i]) / 100f0
        hgmdcr = BM_HT_CRA * crf^BM_HT_CRB * exp(BM_HT_CRC * crf)
        hgmdcr > 1f0 && (hgmdcr = 1f0)
        relht = avh > 0f0 ? h / avh : 0f0
        relht > 1.5f0 && (relht = 1.5f0)
        rhx = relht
        fctrkx = (BM_HT_RHK / BM_RHYXS[sp])^(BM_RHM[sp] - 1f0) - 1f0
        fctrrb = -1f0 * (BM_RHR[sp] / (1f0 - BM_RHB[sp]))
        fctrxb = rhx^(1f0 - BM_RHB[sp]) - BM_HT_RHXS^(1f0 - BM_RHB[sp])
        fctrm  = -1f0 / (BM_RHM[sp] - 1f0)
        hgmdrh = BM_HT_RHK * (1f0 + fctrkx * exp(fctrrb * fctrxb))^fctrm
        htgmod = 0.25f0 * hgmdcr + 0.75f0 * hgmdrh
        htgmod >= 2f0 && (htgmod = 2f0)
        htgmod <= 0.1f0 && (htgmod = 0.1f0)
        htg = pothtg * htgmod
        # cap at HTMAX, then scale + calibration
        (h + htg > htmax) && (htg = htmax - h)
        htg < 0.1f0 && (htg = 0.1f0)
        t.ht_growth[i] = scale * htg * htcon
    end
    return s
end
