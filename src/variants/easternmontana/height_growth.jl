# =============================================================================
# height_growth.jl (easternmontana) — EM large-tree height growth (em/htgf.f + em/pothtg.f).
#
# Main Wykoff conifers (sp1-3,7-10,18): POTHTG (per-species site-index height-age curve → potential
# 10-yr height growth PHTG, dubs t.birth_age=EFAGE if ≤0) → RALPH correction → relative-height modifier.
#   RALPH: PHTG = PHTG·0.706·(1−e^(−10.19·CR))·(1−e^(−1.8158·DG))^0.944 + 0.0265·H
#   modifier: RLHTMD = exp(C3MOD·((H/AVH)^C4MOD − 1)); HTMOD=min(RLHTMD,1); HTG=max(PHTG·HTMOD, 0.1)
#   C3MOD=2.54119, C4MOD=0.250537 (em/htgf.f:130-131, fixed). SI50 for DF(sp3), SI100 for the others.
# LL(5)/LM(4)/aspen(12,17)/RM(6)/CO(11,13-16,19) special forms deferred (not in emt01).
# =============================================================================

const _EM_C3MOD = 2.54119f0
const _EM_C4MOD = 0.250537f0

# em/pothtg.f: potential 10-yr height growth for the main conifers. Sets t.birth_age[i]=EFAGE if ≤0.
@inline function em_pothtg(sp::Int, h::Float32, si::Float32, i::Int, t)::Float32
    if sp == 1 || sp == 2 || sp == 7 || sp == 18          # WB/WL/LP/OS — Alexander-Tackle-Dahms (SI100)
        ccf = 125f0
        a = 9.72443f0 - 0.00091f0 * si * ccf - h
        b = -0.23733f0 + 0.0149f0 * si
        c = 0.00160f0 - 0.00005f0 * si
        tem = b * b - 4f0 * a * c; tem < 0f0 && (tem = 0f0)
        efage = c != 0f0 ? (-b + sqrt(tem)) / (2f0 * c) : -a / b
        t.birth_age[i] <= 0f0 && (t.birth_age[i] = efage)
        e10 = efage + 10f0
        return -2.3733f0 + 0.00160f0 * (e10*e10 - efage*efage) + 0.1490f0 * si - 0.00005f0 * si * (e10*e10 - efage*efage)
    elseif sp == 3                                        # DF — Monserud 1985 (SI50)
        z1 = 0.3197f0; z2 = 1.0232f0
        temht = h - 4.5f0; temht == 1f0 && (temht = 1.1f0)
        temsi = si - 4.5f0
        term1 = (42.397f0 * temsi^z1) / temht - 1f0
        term1 <= 0f0 && return 0f0
        term1 = log(term1)
        efage = exp((term1 + z2 * log(temsi) - 9.7278f0) / (-1.2934f0))
        t.birth_age[i] <= 0f0 && (t.birth_age[i] = efage)
        e10 = efage + 10f0
        h10 = (42.397f0 * temsi^z1) / (1f0 + exp(9.7278f0 - 1.2934f0 * log(e10) - log(temsi) * z2))
        return h10 - temht
    elseif sp == 8 || sp == 9                             # ES/AF — Alexander look-alike (SI100)
        a = si; b = 0.931764f0; c = 0.01679f0; d = 0.302381f0
        h >= a && return 1f0
        term1 = log((1f0 - ((h / a)^(1f0 - d))) / b)
        efage = term1 / (-c)
        t.birth_age[i] <= 0f0 && (t.birth_age[i] = efage)
        e10 = efage + 10f0
        h10 = a * (1f0 - b * exp(-c * e10))^(1f0 / (1f0 - d))
        return h10 - h
    elseif sp == 10                                       # PP — Meyer look-alike (SI100)
        p1 = 3.635794f0; p2 = 0.916307f0; p3 = 6.09478f0; p4 = 0.96483f0; p5 = 0.277025f0
        term1 = ((p1 * si^p2) / h) - 1f0
        term1 <= 0f0 && return 0f0
        term1 = log(term1)
        efage = exp((term1 - p3 + p5 * log(si)) / (-p4))
        t.birth_age[i] <= 0f0 && (t.birth_age[i] = efage)
        e10 = efage + 10f0
        h10 = p1 * si^p2 / (1f0 + exp(p3 - p4 * log(e10) - p5 * log(si)))
        return h10 - h
    end
    return 0f0
end

function height_growth!(s::StandState, ::EasternMontana; scale::Float32 = 1.0f0)
    p, t = s.plot, s.trees
    avh = p.avg_height
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (d <= 0f0 || h <= 0f0) && continue
        if sp <= 3 || (7 <= sp <= 10) || sp == 18          # main Wykoff conifers
            h <= 4.5f0 && continue                          # → REGENT small-tree (chunk 6)
            cr = Float32(t.crown_pct[i]) / 100f0
            dg = t.diam_growth[i]
            si = p.sp_site_index[sp]                        # SI50 (DF) / SI100 (others), per site_setup!
            phtg = em_pothtg(sp, h, si, i, t)
            phtg = phtg * 0.706f0 * (1f0 - exp(-10.19f0 * cr)) * (1f0 - exp(-1.8158f0 * dg))^0.944f0 + 0.0265f0 * h
            rlhtmd = avh > 0f0 ? exp(_EM_C3MOD * ((h / avh)^_EM_C4MOD - 1f0)) : 1f0
            htmod = rlhtmd > 1f0 ? 1f0 : rlhtmd
            htg = phtg * htmod
            htg < 0.1f0 && (htg = 0.1f0)
            t.ht_growth[i] = htg * scale
        else
            error("EM height_growth! sp $sp (LL/LM/aspen/RM/CO) not yet ported")
        end
    end
    return s
end
