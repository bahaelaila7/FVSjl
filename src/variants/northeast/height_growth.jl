# =============================================================================
# height_growth.jl (northeast) — NE periodic height growth (ne/htgf.f + htcalc.f)
#
# Structurally different from SN. The height/age curve is the NC-128
# (Carmean/Hahn/Jacobs) site-index family, indexed per species by MAPNE into the
# 127-row LTBHEC table (B1..B5 plus a breast-height offset BH):
#
#   HTMAX = B1·SI^B2
#   H(age) = BH + B1·SI^B2·(1 − exp(B3·age))^(B4·SI^B5)              (htcalc.f:404,412)
#   AGE(H) = (1/B3)·ln(1 − ((H−BH)/(B1·SI^B2))^(1/(B4·SI^B5)))       (htcalc.f:394)
#   HTG1   = H(age+10) − H(age)                                       (10-yr increment)
#
# The increment is then modulated NOT by SN's crown/relative-height blend but by
# the SAME BAL competition modifier used in diameter growth (ne/balmod.f), tempered
# by relative height, plus the DG serial-correlation random OLDRN (htgf.f:102-114):
#
#   GMOD   = BALMOD(sp, DBH)                       (= exp(−dg_b3·BAL), clamp ≥0.5)
#   RELHTA = AVH>0 ? min(HT/AVH, 1) : 0
#   GMOD   = (1 − (1−GMOD)·(1−RELHTA))·0.8
#   HTG    = HTG1·(1+OLDRN)·GMOD ; floor 0.1
#   HTG    = SCALE·XHT·HTG·exp(HTCON)             (SCALE=FINT/YR, YR=10; XHT=HTGMULT)
#   size-cap to SIZCAP[sp,4]
#
# HTCON (= ln(HCOR2) when the small-tree height calibration LHCOR2 is on; else 0)
# rides in via calib.htg_cor — currently 0 for NE (HCOR2 calibration not yet ported).
# =============================================================================

include("_htcalc_tables.jl")   # const MAPNE (108), LTBHEC (127×(B1..B5,BH))

"HTCALC curve coefficients (B1,B2,B3,B4,B5,BH) for NE species `sp` via MAPNE."
@inline _ne_htcoef(sp::Integer) = @inbounds LTBHEC[MAPNE[sp]]

ne_htcalc_htmax(sp::Integer, si::Real) = (b = _ne_htcoef(sp); b[1] * Float32(si)^b[2])

"Tree age from current height (HTCALC mode 0). Returns 0 if at/above HTMAX."
function ne_htcalc_age(sp::Integer, si::Real, h::Real)
    b1,b2,b3,b4,b5,bh = _ne_htcoef(sp); sif = Float32(si)
    base = (Float32(h) - bh) / (b1 * sif^b2)
    base <= 0f0 && return 0f0                       # tree below breast-height offset
    return (1f0 / b3) * log(1f0 - base^(1f0 / (b4 * sif^b5)))
end

"Height (ft) at `age` on the NC-128 curve (HTCALC mode 1)."
function ne_htcalc_height(sp::Integer, si::Real, age::Real)
    b1,b2,b3,b4,b5,bh = _ne_htcoef(sp); sif = Float32(si)
    return bh + b1 * sif^b2 * (1f0 - exp(b3 * Float32(age)))^(b4 * sif^b5)
end

"10-yr height increment from starting age `aget` (HTCALC mode 9)."
function ne_htcalc_incr(sp::Integer, si::Real, aget::Real, period::Real = 10f0)
    b1,b2,b3,b4,b5,bh = _ne_htcoef(sp); sif = Float32(si); a = Float32(aget)
    hmax = b1 * sif^b2; ex = b4 * sif^b5
    h0  = bh + hmax * (1f0 - exp(b3 * a))^ex
    # htcalc.f:413 evaluates the curve at AGET+YRS (the ACTUAL cycle length), NOT a fixed 10-yr
    # increment linearly scaled — the curve decelerates, so a 5-yr step is MORE than half a decade.
    # `period` defaults to 10 (the NE native cycle) so every existing caller is unchanged.
    hp5 = bh + hmax * (1f0 - exp(b3 * (a + Float32(period))))^ex
    return hp5 - h0
end

"""
    height_growth!(state, ::Northeast; scale)

NE HTGF: periodic height increment into `trees.ht_growth`. For each tree, invert the
NC-128 height/age curve, take the 10-yr increment, modulate by the BAL competition
modifier + relative height + the DG serial-correlation random OLDRN, then apply the
HTGMULT multiplier, SCALE (=FINT/YR), and HTCON. `scale = FINT/10` for NE.
"""
function height_growth!(s::StandState, ::Northeast; scale::Float32 = 1f0)
    p, t, c, sd = s.plot, s.trees, s.calib, s.coef.species
    b3_dg = sd[:dg_b3]                              # BALMOD uses the DG b3
    avh = p.avg_height
    oldrn = t.old_random
    ebau = zeros(Float32, 50); ne_badist!(ebau, s)  # cycle-start BAL basis (same as DGF)
    cur_year = current_cycle_year(s)
    @inbounds for i in 1:t.n
        sp = t.species[i]; hti = t.height[i]
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0f0 && continue
        si = p.sp_site_index[sp]
        htmax = ne_htcalc_htmax(sp, si)
        xht = active_multiplier(s.control, :htg, sp, cur_year)
        htcon = c.htg_cor[sp]
        # htcalc.f:389 — within 1 ft of HTMAX (or below BH): HTG1=0, falls to the 0.1 floor.
        if htmax - hti <= 1f0
            htg1 = 0f0; aget = 0f0
        else
            aget = ne_htcalc_age(sp, si, hti)
            # Evaluate the NC-128 curve over the ACTUAL cycle length (period = scale·YR = FINT),
            # per htcalc.f:413 (AGET+YRS). At the native 10-yr cycle scale=1 ⇒ period=10 (unchanged);
            # at a non-native 5-yr cycle this gives the true decelerated 5-yr step instead of ½·(10-yr),
            # which jl previously under-grew. The linear `scale` below is therefore dropped (the period
            # is now baked into the curve evaluation) — native stays bit-exact (scale was 1.0).
            htg1 = ne_htcalc_incr(sp, si, aget, scale * 10f0)
        end
        gmod = ne_balmod(b3_dg[sp], ebau, t.dbh[i])
        relht = avh > 0f0 ? min(hti / avh, 1f0) : 0f0
        gmod = (1f0 - (1f0 - gmod) * (1f0 - relht)) * 0.8f0
        htg = htg1 * (1f0 + oldrn[i]) * gmod
        htg < 0.1f0 && (htg = 0.1f0)
        htg = xht * htg * exp(htcon)
        # size cap (SIZCAP[sp,4], default 999 for NE — TREESZCP keyword only)
        sc4 = s.control.sp_size_cap[sp, 4]
        (hti + htg) > sc4 && (htg = max(sc4 - hti, 0.1f0))
        t.ht_growth[i] = htg
    end
    return s
end
