# =============================================================================
# height_growth.jl — Southern height growth (HTGF + HTCALC)
#
# Ported from: sn/htgf.f, sn/htcalc.f.
#
# HTCALC is a Chapman-Richards height/age curve: HTMAX = b1·SI^b2, and the height
# at an age is hb + HTMAX·(1−exp(b3·age))^(b4·SI^b5). HTGF inverts it to get the
# tree's age from its current height (mode 0), takes the 5-year increment from that
# age (mode 9), then scales it by crown-ratio and relative-height modifiers and the
# height calibration HTCON. Tables copied verbatim.
# =============================================================================

# HTCALC curve coefficients (ht_curve_b1..b5) and the relative-height modifier
# parameters (ht_relht_yxs/r/b) are loaded from data/southern/species_coefficients.csv.

const HTGF_RHM = 1.10f0; const HTGF_RHXS = 0.0f0; const HTGF_RHK = 1.0f0
const HTGF_CRA = 100.0f0; const HTGF_CRB = 3.0f0; const HTGF_CRC = -5.0f0
const HTGF_REGYR = 5.0f0

# b1..b5 for a species (HTCALC), with the yellow-poplar (sp 45, non-montane) special
# case. `bc` is a 5-tuple of the per-species ht_curve_b* columns.
@inline function _htcalc_coef(bc, sp::Integer, montane::Bool)
    if sp == 45 && !montane
        return (1.1789f0, 1.0f0, -0.0339f0, 0.8117f0, -0.0001f0)
    end
    (bc[1][sp], bc[2][sp], bc[3][sp], bc[4][sp], bc[5][sp])
end

htcalc_htmax(bc, sp::Integer, si::Real, montane::Bool=false) =
    (b = _htcalc_coef(bc, sp, montane); b[1] * Float32(si)^b[2])

"Height (ft) at a given age on the Chapman-Richards curve (HTCALC mode 1) — used by
ESSUBH to assign established-tree heights."
function htcalc_height(bc, sp::Integer, si::Real, age::Real, montane::Bool=false)
    b1,b2,b3,b4,b5 = _htcalc_coef(bc, sp, montane); sif = Float32(si)
    return (b1 * sif^b2) * (1f0 - exp(b3 * Float32(age)))^(b4 * sif^b5)
end

"Solve tree age from current height (HTCALC mode 0)."
function htcalc_age(bc, sp::Integer, si::Real, h::Real, montane::Bool=false)
    b1,b2,b3,b4,b5 = _htcalc_coef(bc, sp, montane); sif = Float32(si)
    ratio = Float32(h) / (b1 * sif^b2)
    ratio = clamp(ratio, 0f0, 1f0 - 1f-6)
    return (1f0 / b3) * log(1f0 - ratio^(1f0 / (b4 * sif^b5)))
end

"5-year height increment from a starting age (HTCALC mode 9)."
function htcalc_incr(bc, sp::Integer, si::Real, aget::Real, montane::Bool=false)
    b1,b2,b3,b4,b5 = _htcalc_coef(bc, sp, montane); sif = Float32(si); a = Float32(aget)
    hmax = b1 * sif^b2; ex = b4 * sif^b5
    h0  = hmax * (1f0 - exp(b3 * a))^ex
    hp5 = hmax * (1f0 - exp(b3 * (a + 5f0)))^ex
    return hp5 - h0
end

"""
    height_growth!(state, ::Southern)

Variant hook: periodic height growth into `trees.ht_growth` (HTGF). For each tree,
invert the height/age curve, take the 5-yr increment, apply crown-ratio + relative-
height modifiers and HTCON calibration. `scale = FINT/YR` (1 for snt01).
"""
function height_growth!(s::StandState, ::Southern; scale::Float32 = 1f0)
    p, t, c, sd = s.plot, s.trees, s.calib, s.coef.species
    bc = (sd[:ht_curve_b1], sd[:ht_curve_b2], sd[:ht_curve_b3], sd[:ht_curve_b4], sd[:ht_curve_b5])
    rhyxs_v = sd[:ht_relht_yxs]; rhr_v = sd[:ht_relht_r]; rhb_v = sd[:ht_relht_b]
    montane = !isempty(p.eco_unit) && p.eco_unit[1] == 'M'
    avh = p.avg_height
    # HTGMULT (MULTS kind 2): per-species height-growth multiplier — htgf.f:163,260 apply
    # XHT at HTG = HTG·XHT·SCALE·EXP(HTCON). cur_year = inventory year + cycle·period.
    cur_year = current_cycle_year(s)   # IY schedule (TIMEINT/CYCLEAT-aware)
    @inbounds for i in 1:t.n
        sp = t.species[i]; hti = t.height[i]
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0f0 && continue
        xht = active_multiplier(s.control, :htg, sp, cur_year)
        si = p.sp_site_index[sp]
        htmax = htcalc_htmax(bc, sp, si, montane)
        htcon = c.htg_cor[sp]
        if htmax - hti <= 1f0
            t.ht_growth[i] = 0.10f0 * xht * scale * exp(htcon)
            sc4 = s.control.sp_size_cap[sp, 4]
            (hti + t.ht_growth[i]) > sc4 && (t.ht_growth[i] = max(sc4 - hti, 0.1f0))
            continue
        end
        aget = htcalc_age(bc, sp, si, hti, montane)
        htg1 = htcalc_incr(bc, sp, si, aget, montane)
        relht = avh > 0f0 ? min(hti / avh, 1.5f0) : 0f0
        cr = Float32(t.crown_pct[i]) / 100f0
        hgmdcr = min(HTGF_CRA * cr^HTGF_CRB * exp(HTGF_CRC * cr), 1f0)
        rhyxs = rhyxs_v[sp]; rhr = rhr_v[sp]; rhb = rhb_v[sp]
        fctrkx = (HTGF_RHK / rhyxs)^(HTGF_RHM - 1f0) - 1f0
        fctrrb = -1f0 * (rhr / (1f0 - rhb))
        fctrxb = relht^(1f0 - rhb) - HTGF_RHXS^(1f0 - rhb)
        fctrm  = -1f0 / (HTGF_RHM - 1f0)
        hgmdrh = HTGF_RHK * (1f0 + fctrkx * exp(fctrrb * fctrxb))^fctrm
        htgmod = clamp(0.25f0 * hgmdcr + 0.75f0 * hgmdrh, 0.1f0, 2f0)
        htg = max(htg1 * htgmod, 0.1f0)
        t.ht_growth[i] = htg * xht * scale * exp(htcon)
        # htgf.f:286-288 — large-tree height cap (SIZCAP[4], set by TREESZCP). Default 999.
        # A tree already past the cap still crawls by the 0.1 floor (never shrinks).
        sc4 = s.control.sp_size_cap[sp, 4]
        (hti + t.ht_growth[i]) > sc4 && (t.ht_growth[i] = max(sc4 - hti, 0.1f0))
    end
    return s
end
