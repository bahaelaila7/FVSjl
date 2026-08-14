# =============================================================================
# small_tree_growth.jl (southcentraloregon) — SO small-tree growth (so/regent.f + so/smhtgf.f). Chunk 6.
#
# STEP 1 (this file): so_smhtgf — the per-species small-tree 5/10-yr HEIGHT increment (so/smhtgf.f SELECT
# CASE(ISPC)). Height-first small-tree model: SMHTGF gives HHT, regent then converts to DBH for trees
# crossing 4.5' and blends with the large-tree HTG via XWT. Most species are a linear-SI increment
# ((a+b·SI)/(c+d·SI))·DTIME·factor; WP(1)/RC(18) are the Chapman-Richards age-inversion form; WJ(11) a
# SI-ratio form; RF(9) the CA Ritchie&Hann DOMHTGR·SMHMOD; WO(27) the CA blackoak BAL form. MODE 0=from
# ESSUBH (estab), 1=from REGENT (invert age from H for the nonlinear species). MEASURED vs FVSso_g16 smhtgf.
#
# STEP 2 (TODO — the regent DBH driver): small_tree_growth!(::SouthCentralOregon) mirroring the CA driver,
# using so/regent.f's OWN height-DBH regression (HTT1/HTT2 per-crown-group, HT1/HT2, AA, IABFLG from
# so/blkdat.f — NOT the chunk-4a Curtis so_htdbh) + XMAX/XMIN/DIAM/REGYR. Needed to unblock grow_cycle!.
# =============================================================================

# so/smhtgf.f species case sets.
const SO_SMHT_WFGRP = (4, 12)      # WF/GF (·1.2)
const SO_SMHT_PYGRP = (20, 21, 23, 25, 26, 28, 29, 30, 31, 33)  # WC PY-form

# so/smhtgf.f — 5/10-yr small-tree height increment (ft). `cr1` = ICR (integer percent, 0..100).
function so_smhtgf(sp::Int, d::Float32, h::Float32, cr1::Float32, dtime::Float32, mode::Int;
                   si::Float32, ba::Float32, pct::Float32, avh::Float32)::Float32
    S = si
    if sp == 1 || sp == 18                                  # WP / RC — Chapman-Richards age inversion
        c1, c2, c3, c4 = sp == 1 ? (0.375045f0, 0.92503f0, -0.020796f0, 2.48811f0) :
                                   (0.752842f0, 1.0f0, -0.0174f0, 1.4711f0)
        effage = mode == 1 ? log((1f0 - (c1 / S * h)^(1f0 / c4)) / c2) / c3 : 0f0
        agepdt = effage + dtime
        hht1 = (S / c1) * (1f0 - c2 * exp(c3 * effage))^c4
        hht2 = (S / c1) * (1f0 - c2 * exp(c3 * agepdt))^c4
        return hht2 - hht1
    elseif sp == 2 || sp == 10                              # SP / PP
        return ((-1.0f0 + 0.32857f0 * S) / (28.0f0 - 0.042857f0 * S)) * dtime
    elseif sp == 3 || sp == 32                              # DF / OS
        return ((2.0f0 + 0.420f0 * S) / (28.5f0 - 0.05f0 * S)) * dtime * 1.1f0
    elseif sp in SO_SMHT_WFGRP                              # WF / GF
        return ((4.2435f0 + 0.1510f0 * S) / (19.0184f0 - 0.0570f0 * S)) * dtime * 1.2f0
    elseif sp == 5                                          # MH (metric → ft)
        return ((0.965758f0 + 0.082969f0 * S) / (55.249612f0 - 1.288852f0 * S)) * dtime * 3.280833f0 * 1.60f0
    elseif sp == 6                                          # IC
        return ((4.2435f0 + 0.1510f0 * S) / (19.0184f0 - 0.0570f0 * S)) * dtime * 1.3f0
    elseif sp == 7 || sp == 16                              # LP / WB
        h0 = 0.02008805f0 * S * dtime
        return sp == 16 ? 1.6f0 * h0 : h0
    elseif sp == 8                                          # ES
        return ((0.09211f0 + 0.208517f0 * S) / (43.358f0 - 0.168166f0 * S)) * dtime * 1.35f0
    elseif sp == 9                                          # RF — CA Ritchie&Hann
        relht = avh <= 0f0 ? 1f0 : h / avh
        relht > 1.05f0 && (relht = 1.05f0)
        domhtgr = 5f0 * (2.2227f0 + 0.4314f0 * S) / (29.0f0 - 0.05f0 * S)
        cr = cr1 / 100f0
        crmod = 1.0f0 - exp(-4.26558f0 * cr)
        rhmod = exp(2.54119f0 * (relht^0.250537f0 - 1.0f0))
        return domhtgr * 1.016605f0 * crmod * rhmod
    elseif sp == 11                                        # WJ — SI-ratio (SI clamp [SLO+.5, SHI])
        sj = S; sj > 75f0 && (sj = 75f0); sj <= 5f0 && (sj = 5.5f0)
        return (sj / 5.0f0) * (sj * 1.5f0 - h) / (sj * 1.5f0)
    elseif sp == 13                                        # AF
        return ((6.0f0 + 0.14f0 * S) / (33.882f0 - 0.06588f0 * S)) * dtime
    elseif sp == 14                                        # SF (EC form)
        return ((-0.6667f0 + 0.4333f0 * S) / (28.5f0 - 0.05f0 * S)) * dtime
    elseif sp == 15                                        # NF (WC form)
        return ((11.26677f0 + 0.12027f0 * S) / (27.93806f0 - 0.02873f0 * S)) * dtime
    elseif sp == 17                                        # WL (EC form)
        return ((-3.9725f0 + 0.50995f0 * S) / (28.1168f0 - 0.05661f0 * S)) * dtime
    elseif sp == 19                                        # WH (WC form)
        return ((-5.74874f0 + 0.54576f0 * S) / (26.15767f0 - 0.03596f0 * S)) * dtime
    elseif sp in SO_SMHT_PYGRP                             # PY + WC catch-all
        return ((1.47043f0 + 0.23317f0 * S) / (31.56252f0 - 0.05586f0 * S)) * dtime
    elseif sp == 22                                        # RA (WC form)
        return (-0.007025f0 + 0.056794f0 * S) * dtime * 1.20f0
    elseif sp == 27                                        # WO — CA blackoak (BAL form)
        bal = ((100.0f0 - pct) / 100.0f0) * ba
        factor = 0.80f0 + 0.004f0 * (S - 50.0f0)
        bal < 5.0f0 && (bal = 5.0f0)
        return exp(3.817f0 - 0.7829f0 * log(bal)) * factor
    end
    return 0f0
end
