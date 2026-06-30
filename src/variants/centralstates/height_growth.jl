# =============================================================================
# height_growth.jl (centralstates) — CS periodic height growth (cs/htgf.f + htcalc.f)
#
# Same SHAPE as NE's height growth (the NC-128 Carmean/Hahn site-index curve via the
# shared HTCALC machinery, then a BAL competition modifier tempered by relative height
# and the OLDRN serial-correlation residual), with TWO CS specifics:
#   * the species→curve map is MAPCS (htcalc.f IVAR=2), not MAPNE — same LTBHEC curve
#     table (htcalc.f is byte-identical CS/NE; only the per-variant map differs).
#   * the BAL modifier is CS's own `cs_balmod` (cs/balmod.f: B3·(1−exp(−(B1/(BAL+TRBA)
#     +B2·D²)·√(1−min(BA,200)/210))), clamp ≥0.15) — NOT NE's simpler exp(−b3·BAL).
#
#   HTG1 = HTCALC(mode 0→9, IVAR=2, SI, H)            (10-yr NC-128 increment)
#   GMOD = cs_balmod(sp, D, BAL, BA)                  (cs/balmod.f)
#   GMOD = (1 − (1−GMOD)·(1−min(H/AVH,1)))·0.8        (cs/htgf.f:48)
#   HTG  = HTG1·(1+OLDRN)·GMOD                        (cs/htgf.f:49)
# =============================================================================

# MAPCS: CS species (1-96) → row in the shared LTBHEC NC-128 curve table (htcalc.f IVAR=2).
const MAPCS = let
    path = joinpath(CS_DATADIR, "htcalc_map.csv")
    m = zeros(Int, 96)
    for (k, line) in enumerate(eachline(path))
        k == 1 && continue
        f = split(strip(line), ',')
        m[parse(Int, f[1])] = parse(Int, f[3])
    end
    m
end

"HTCALC curve coefficients (B1,B2,B3,B4,B5,BH) for CS species `sp` via MAPCS (shared LTBHEC)."
@inline _cs_htcoef(sp::Integer) = @inbounds LTBHEC[MAPCS[sp]]

cs_htcalc_htmax(sp::Integer, si::Real) = (b = _cs_htcoef(sp); b[1] * Float32(si)^b[2])

"Tree age from current height (HTCALC mode 0). Returns 0 if at/above HTMAX."
function cs_htcalc_age(sp::Integer, si::Real, h::Real)
    b1, b2, b3, b4, b5, bh = _cs_htcoef(sp); sif = Float32(si)
    base = (Float32(h) - bh) / (b1 * sif^b2)
    base <= 0f0 && return 0f0
    return (1f0 / b3) * log(1f0 - base^(1f0 / (b4 * sif^b5)))
end

"10-yr height increment from starting age `aget` (HTCALC mode 9)."
function cs_htcalc_incr(sp::Integer, si::Real, aget::Real)
    b1, b2, b3, b4, b5, bh = _cs_htcoef(sp); sif = Float32(si); a = Float32(aget)
    hmax = b1 * sif^b2; ex = b4 * sif^b5
    h0  = bh + hmax * (1f0 - exp(b3 * a))^ex
    hp10 = bh + hmax * (1f0 - exp(b3 * (a + 10f0)))^ex
    return hp10 - h0
end

"CS BAL competition modifier for height growth (cs/balmod.f)."
@inline function cs_balmod(b1::Float32, b2::Float32, b3::Float32,
                           bal::Float32, ba::Float32, d::Float32)::Float32
    trba = d * d * 0.005454f0
    temba = ba > 200f0 ? 200f0 : ba
    part1 = -1f0 * (b1 / (bal + trba) + b2 * d * d)
    part2 = sqrt(1f0 - temba / 210f0)
    part3 = part1 * part2
    part4 = part3 > -85f0 ? exp(part3) : 0f0
    gmod = b3 * (1f0 - part4)
    gmod < 0.15f0 ? 0.15f0 : gmod
end

"""
    height_growth!(state, ::CentralStates; scale)

CS HTGF (cs/htgf.f): periodic height increment into `trees.ht_growth`. NC-128 curve
inversion (MAPCS) → 10-yr increment → cs_balmod·relative-height temper → (1+OLDRN)
serial correlation → XHMULT/SCALE/HCON, with the htcalc.f:389 near-HTMAX zero guard.
"""
function height_growth!(s::StandState, ::CentralStates; scale::Float32 = 1f0)
    p, t, c, sd = s.plot, s.trees, s.calib, s.coef.species
    b1 = sd[:balmod_b1]; b2 = sd[:balmod_b2]; b3 = sd[:balmod_b3]
    ba = p.basal_area
    avh = p.avg_height
    oldrn = t.old_random
    cur_year = current_cycle_year(s)
    @inbounds for i in 1:t.n
        sp = t.species[i]; hti = t.height[i]
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0f0 && continue
        si = p.sp_site_index[sp]
        htmax = cs_htcalc_htmax(sp, si)
        xht = active_multiplier(s.control, :htg, sp, cur_year)
        htcon = c.htg_cor[sp]
        if htmax - hti <= 1f0                       # htcalc.f:389 — at/near HTMAX
            htg1 = 0f0
        else
            aget = cs_htcalc_age(sp, si, hti)
            htg1 = cs_htcalc_incr(sp, si, aget)
        end
        bal = (1f0 - t.crown_ratio[i] / 100f0) * ba   # PCT = BA percentile (cs/htgf.f:71)
        gmod = cs_balmod(b1[sp], b2[sp], b3[sp], bal, ba, t.dbh[i])
        relht = avh > 0f0 ? min(hti / avh, 1f0) : 0f0
        gmod = (1f0 - (1f0 - gmod) * (1f0 - relht)) * 0.8f0
        htg = htg1 * (1f0 + oldrn[i]) * gmod
        htg < 0.1f0 && (htg = 0.1f0)
        htg = scale * xht * htg * exp(htcon)
        sc4 = s.control.sp_size_cap[sp, 4]
        (hti + htg) > sc4 && (htg = max(sc4 - hti, 0.1f0))
        t.ht_growth[i] = htg
    end
    return s
end
