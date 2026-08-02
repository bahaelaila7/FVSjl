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

# em/htgf.f COFLM (LM, sp4) / COFAS (aspen+CO, sp11-17,19) Schreuder-Hafley SB height coeffs, by crown
# class K=1..3 (columns). Row order = COF1..COF9. K from IICR=int(ICR/10+0.5): {1,2}→1 {3-7}→2 {8,9}→3.
const _EM_COFLM = Float32[  # [9 coeffs × 3 crown classes]
    37.0    45.0    45.0
    85.0   100.0    90.0
    1.77836 1.66674 1.64770
   -0.51147 0.25626 0.30546
    1.88795 1.45477 1.35015
    1.20654 1.11251 0.94823
    0.57697 0.67375 0.70453
    3.57635 2.17942 2.46480
    0.90283 0.88103 1.00316]
const _EM_COFAS = Float32[
    30.0    30.0    35.0
    85.0    85.0    85.0
    2.00995 2.00995 1.80388
    0.03288 0.03288 -0.07682
    1.81059 1.81059 1.70032
    1.28612 1.28612 1.29148
    0.72051 0.72051 0.72343
    3.00551 3.00551 2.91519
    1.01433 1.01433 0.95244]

# em/htgf.f HTCONS: LL(sp5) height habitat constants. ITYPE(1..30) → IHT(1..8) via MAPHAB; per-IHT
# HGHC (intercept), HGLDD (ln DG coef), HGH2 (H² coef). HTCON(5)=HGHC(IHT)−0.5478 (others 0).
const _EM_HT_MAPHAB = Int[1,1,2,2,2,2,2,2,2,3,3,4,5,6,7,7,7,7,4,4,1,4,4,8,8,8,1,1,1,1]
const _EM_HGHC  = Float32[2.03035, 1.72222, 1.19728, 1.81759, 2.14781, 1.76998, 2.21104, 1.74090]
const _EM_HGLDD = Float32[0.62144, 1.02372, 0.85493, 0.75756, 0.46238, 0.49643, 0.37042, 0.34003]
const _EM_HGH2  = Float32[-13.358f-5, -3.809f-5, -3.715f-5, -2.607f-5, -5.200f-5, -1.605f-5, -3.631f-5, -4.460f-5]

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
    p, t, c = s.plot, s.trees, s.calib
    avh = p.avg_height
    # LL(sp5) habitat-dependent height constants (em/htgf.f HTCONS), resolved once per stand.
    itype = Int(p.habitat_input); (itype < 1 || itype > 30) && (itype = 1)
    iht = _EM_HT_MAPHAB[itype]
    ll_h2cof = _EM_HGH2[iht]; ll_hdgcof = _EM_HGLDD[iht]
    ll_htcon = _EM_HGHC[iht] - 0.5478f0
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
        elseif sp == 6
            # RM juniper: em/htgf.f:231 GO TO 30 — ALL height growth comes from REGENT. Leave 0.
            continue
        elseif sp == 5
            # LL (NI form, em/htgf.f:222-226): CON = HTCON + H2COF·H² − 0.1997·lnD + 0.23315·lnH.
            h <= 4.5f0 && continue
            dg = t.diam_growth[i]
            con = ll_htcon + ll_h2cof * h * h - 0.1997f0 * log(d) + 0.23315f0 * log(h)
            htg = dg > 0f0 ? exp(con + ll_hdgcof * log(dg)) + 0.4809f0 : 0.1f0
            htg < 0.1f0 && (htg = 0.1f0)
            t.ht_growth[i] = htg * scale
        else
            # LM(4)/CO(11,13-16,19)/aspen(12,17) — Schreuder-Hafley SB height (em/htgf.f:236-358). COFLM for
            # sp4, COFAS otherwise. (Young-tree accelerator at :311 is dead code — IAGE is never set ⇒ 0.)
            iicr = trunc(Int, Float32(t.crown_pct[i]) / 10f0 + 0.5f0); iicr > 9 && (iicr = 9); iicr < 1 && (iicr = 1)
            k = iicr <= 2 ? 1 : (iicr <= 7 ? 2 : 3)
            cof = sp == 4 ? _EM_COFLM : _EM_COFAS
            cof1 = cof[1,k]; cof2 = cof[2,k]; cof3 = cof[3,k]; cof4 = cof[4,k]; cof5 = cof[5,k]
            cof6 = cof[6,k]; cof7 = cof[7,k]; cof8 = cof[8,k]; cof9 = cof[9,k]
            # bounds: outside the fitted SB range ⇒ HTG=0.1 (from regent for these small trees)
            if h <= 4.5f0 || (0.1f0 + cof1) <= d || (4.5f0 + cof2) <= h || d <= 0.1f0
                t.ht_growth[i] = 0.1f0 * scale; continue
            end
            temd = d <= 0.2f0 ? 0.2f0 : d
            y1 = (temd - 0.1f0) / cof1; y2 = (h - 4.5f0) / cof2
            fby1 = log(y1 / (1f0 - y1)); fby2 = log(y2 / (1f0 - y2))
            z = (cof4 + cof6 * fby2 - cof7 * (cof3 + cof5 * fby1)) * (1f0 - cof7 * cof7)^(-0.5f0)
            if sp != 4
                zadj = 0.1f0 - 0.10273f0 * z + 0.00273f0 * z * z
                zadj < 0f0 && (zadj = 0f0)
                z += zadj
            end
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
            dia = d + t.diam_growth[i] / bark
            if (0.1f0 + cof1) > dia
                psi = cof8 * ((dia - 0.1f0) / (0.1f0 + cof1 - dia))^cof9 *
                      exp(z * ((1f0 - cof7 * cof7)^0.5f0) / cof6)
                hnew = (psi / (1f0 + psi)) * cof2 + 4.5f0
                hnew < h && (hnew = h)
                t.ht_growth[i] = (hnew - h) * scale
            else
                t.ht_growth[i] = 0.1f0 * scale
            end
        end
    end
    return s
end
