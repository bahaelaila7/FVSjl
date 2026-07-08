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
    # htcalc.f:394 (byte-identical CS/NE) — `(H-BH)/B1/SI**B2` and exponent `1./B4/SI**B5` are Fortran
    # LEFT-ASSOC sequential divisions ((a/b)/c), NOT a/(b·c). Match that op order (net01 stays bit-exact
    # because product==sequential for its inputs; this is defensive against inputs where they'd diverge).
    # htcalc.f transcendentals `SI**B2`, `**(...)`, `ALOG` are gfortran powf/logf — route via the FFI
    # companion (doctrine #8), not Julia's native `^`/`log` (openlibm), so the inverted age is bit-exact vs FVS.
    base = (Float32(h) - bh) / b1 / fpow(sif, b2)
    base <= 0f0 && return 0f0                       # tree below breast-height offset
    return (1f0 / b3) * flog(1f0 - fpow(base, 1f0 / b4 / fpow(sif, b5)))
end

"Height (ft) at `age` on the NC-128 curve (HTCALC mode 1)."
function ne_htcalc_height(sp::Integer, si::Real, age::Real)
    b1,b2,b3,b4,b5,bh = _ne_htcoef(sp); sif = Float32(si)
    return bh + b1 * fpow(sif, b2) * fpow(1f0 - fexp(b3 * Float32(age)), b4 * fpow(sif, b5))
end

"10-yr height increment from starting age `aget` (HTCALC mode 9)."
function ne_htcalc_incr(sp::Integer, si::Real, aget::Real)
    b1,b2,b3,b4,b5,bh = _ne_htcoef(sp); sif = Float32(si); a = Float32(aget)
    # SI**B2, **EX, EXP(B3*A) → FFI companion (gfortran) to match htcalc.f bit-for-bit.
    hmax = b1 * fpow(sif, b2); ex = b4 * fpow(sif, b5)
    h0  = bh + hmax * fpow(1f0 - fexp(b3 * a), ex)
    hp5 = bh + hmax * fpow(1f0 - fexp(b3 * (a + 10f0)), ex)
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
        # READCORH/REUSCORH: HTCON += ln(HCOR2) when LHCOR2 on (htgf.f:172-174, HTCONS entry).
        s.control.htg_cor2_on && s.control.htg_cor2[sp] > 0f0 && (htcon += log(s.control.htg_cor2[sp]))
        # htcalc.f:389 — within 1 ft of HTMAX (or below BH): HTG1=0, falls to the 0.1 floor.
        if htmax - hti <= 1f0
            htg1 = 0f0; aget = 0f0
        else
            aget = ne_htcalc_age(sp, si, hti)
            htg1 = ne_htcalc_incr(sp, si, aget)
        end
        gmod = ne_balmod(b3_dg[sp], ebau, t.dbh[i])
        relht = avh > 0f0 ? min(hti / avh, 1f0) : 0f0
        gmod = (1f0 - (1f0 - gmod) * (1f0 - relht)) * 0.8f0
        htg = htg1 * (1f0 + oldrn[i]) * gmod
        htg < 0.1f0 && (htg = 0.1f0)
        htg = scale * xht * htg * fexp(htcon)          # EXP(HTCON) htcalc.f → FFI companion (gfortran)
        # size cap (SIZCAP[sp,4], default 999 for NE — TREESZCP keyword only)
        sc4 = s.control.sp_size_cap[sp, 4]
        (hti + htg) > sc4 && (htg = max(sc4 - hti, 0.1f0))
        t.ht_growth[i] = htg
    end
    return s
end
