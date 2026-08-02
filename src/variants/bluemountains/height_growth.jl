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

function height_growth!(s::StandState, ::BlueMountains; scale::Float32 = 1.0f0)
    p, t, c = s.plot, s.trees, s.calib
    site_lo = s.coef.species[:site_lo]; site_hi = s.coef.species[:site_hi]
    avh = p.avg_height
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (d <= 0f0 || h <= 0f0) && continue
        # WJ(6) height from REGENT; WB/LM/AS(11,12,15) Johnson SBB (deferred — not in bmt01).
        (sp == 6 || sp == 11 || sp == 12 || sp == 15) && continue
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
