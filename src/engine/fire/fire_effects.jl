# =============================================================================
# fire/fire_effects.jl — fire-caused tree mortality (FFE chunk F6 core)
#
# Ported from: bin/FVSsn_buildDir/fmeff.f (FMEFF mortality) + fmbrkt.f (FMBRKT bark).
#
# The probability that a tree is killed by a fire of a given flame length. This is the
# `.sum`-affecting part of the FFE — it reduces tree TPA (the kill applied by FMBURN).
# Two model forms (fmeff.f):
#   - SN oaks / hickory / red-maple / black-gum use a Regelbrugge & Smith (1994)
#     logistic in DBH and char height (= 0.7 × flame length);
#   - all other species use the Reinhardt crown-scorch + bark-thickness logistic.
#
# Given a flame length and the crown-volume-scorched fraction (which the fire-behavior
# chunk F5 supplies via scorch height), `fire_tree_mortality` is self-contained and
# unit-tested here; wiring it to the burn loop is part of F5/F6 integration.
# =============================================================================

# Bark-thickness multipliers B1 indexed by EQNUM 1..39 (fmbrkt.f); single bark = DBH·B1.
const _FM_BARK_B1 = Float32[
    0.019, 0.022, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.030, 0.031,
    0.032, 0.033, 0.034, 0.035, 0.036, 0.037, 0.038, 0.039, 0.040, 0.041,
    0.042, 0.043, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049, 0.050, 0.052,
    0.055, 0.057, 0.059, 0.060, 0.062, 0.063, 0.068, 0.072, 0.081]

# SN Regelbrugge-Smith fire-mortality logistic coefficients by mortality group 1..5 (fmeff.f).
const _FM_MORTB0 = Float32[1.0229, 0.1683, 1.2165, 0.8221, 2.7750]
const _FM_MORTB1 = Float32[-0.2646, -0.1332, -0.4758, -0.4098, -1.1224]
const _FM_MORTB2 = Float32[2.6232, 3.4152, 6.0415, 8.4682, 2.8312]

"""
    fire_bark_thickness(coef, sp, dbh) -> Float32

Single bark thickness (in) for fire mortality (FMBRKT, fmbrkt.f): `DBH · B1[EQNUM[sp]]`,
with shortleaf pine (sp 5) using the Harmon (1984) quadratic.
"""
@inline function fire_bark_thickness(coef::SpeciesCoefficients, sp::Integer, dbh::Float32)::Float32
    if sp == 5                                          # shortleaf pine (Harmon 1984)
        b = (0.07f0 + 0.09f0 * dbh * 2.54f0 - 0.0001f0 * (dbh * 2.54f0)^2) / 2.54f0
        return max(0f0, b)
    end
    return dbh * _FM_BARK_B1[Int(coef_col(coef, :bark_eqnum)[sp])]
end

"SN fire-mortality group (1–6) for species `sp` (fmeff.f:208-221)."
@inline function fire_mortality_group(sp::Integer)::Int
    sp == 63 || sp == 74               ? 1 :   # white oak, chestnut oak
    sp == 64 || sp == 75 || sp == 78   ? 2 :   # scarlet, black, northern red oak
    sp == 27                           ? 3 :   # hickory
    sp == 20                           ? 4 :   # red maple
    sp == 54                           ? 5 :   # black gum
                                         6     # everything else (Reinhardt)
end

"""
    scorch_height(byram, atemp, fwind) -> Float32

Van Wagner crown scorch height (ft) from Byram fireline intensity `byram`
(BTU/ft/min), air temperature `atemp` (°F), and wind speed `fwind` (fmburn.f:470).
"""
@inline function scorch_height(byram::Float32, atemp::Float32, fwind::Float32)::Float32
    b = byram / 60f0                                    # BTU/ft/min → BTU/ft/sec
    return (63f0 / (140f0 - atemp)) * (b^(7f0 / 6f0) / sqrt(b + fwind^3))
end

"""
    crown_volume_scorched(sch, ht, crown_pct) -> Float32

Percent crown volume scorched (CSV) for a tree of height `ht` (ft) and crown ratio
`crown_pct` (%) under scorch height `sch` (ft) (FMEFF, fmeff.f:170-186). A tree with
no live crown is treated as fully scorched.
"""
@inline function crown_volume_scorched(sch::Float32, ht::Float32, crown_pct::Integer)::Float32
    crl = ht * (Float32(crown_pct) / 100f0)            # crown length
    crl > 0f0 || return 100f0
    sl = sch - (ht - crl)                              # scorch length within the crown
    sl < 0f0 && (sl = 0f0)
    sl > crl && (sl = crl)
    return 100f0 * (sl * (2f0 * crl - sl) / (crl * crl))
end

"""
    fire_tree_mortality(coef, sp, dbh, flame, csv) -> Float32

Probability (0–1) of fire-caused mortality (FMEFF, fmeff.f) for a tree of species `sp`,
DBH `dbh` (in), under flame length `flame` (ft) with crown-volume-scorched `csv` (%).
SN groups 1–5 (oaks/hickory/red-maple/black-gum) use the Regelbrugge-Smith DBH + char-
height logistic (char height = 0.7·flame); group 6 uses the Reinhardt bark-thickness +
crown-scorch logistic.
"""
function fire_tree_mortality(coef::SpeciesCoefficients, sp::Integer, dbh::Float32,
                             flame::Float32, csv::Float32)::Float32
    g = fire_mortality_group(sp)
    if 1 <= g <= 5
        charht = flame * 0.7f0                          # max (uphill) char height
        xm = -(_FM_MORTB0[g] + _FM_MORTB1[g] * dbh * 2.54f0 + _FM_MORTB2[g] * charht / 3.28f0)
        mnmort = log(1f0 / 0.000001f0 - 1f0)            # guard against exp() overflow
        return xm >= mnmort ? 0f0 : 1f0 / (1f0 + exp(xm))
    else                                                # Reinhardt crown-scorch + bark
        bt = fire_bark_thickness(coef, sp, dbh)
        xm = exp(-1.941f0 + 6.316f0 * (1f0 - exp(-bt)) - 0.000535f0 * csv * csv)
        return 1f0 / (1f0 + xm)
    end
end

"""
    fire_mortality_adjust(pmort, sp, dbh, burnseas) -> Float32

SN species/size/season adjustments to the base fire mortality `pmort` (FMEFF,
fmeff.f:281-300): small maples (<4") and very small hardwoods (≤1") are fully killed;
in an early-season burn (`burnseas` ≤ 2) hardwood mortality is reduced — oaks ≥2.5"
halved, other hardwoods ×0.8. Applied in Fortran order (the season reduction can act on
the maple's 1.0 before the ≤1" rule restores it).
"""
@inline function fire_mortality_adjust(pmort::Float32, sp::Integer, dbh::Float32, burnseas::Integer)::Float32
    p = pmort
    (sp == 18 || sp == 19 || sp == 26 || sp == 27 || sp == 51 || sp == 52) && dbh < 4f0 && (p = 1f0)
    if burnseas <= 2 && sp > 14
        p = (30 <= sp <= 36) ? (dbh >= 2.5f0 ? p / 2f0 : p * 0.8f0) : p * 0.8f0
    end
    (sp > 14 && dbh <= 1f0) && (p = 1f0)
    return p
end
