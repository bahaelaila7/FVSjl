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
@inline function fire_bark_thickness(coef::SpeciesCoefficients, sp::Integer, dbh::Float32,
                                     variant::AbstractVariant = Southern())::Float32
    # shortleaf pine uses the Harmon (1984) quadratic — species 5 in SN, species 3 in CS (cs/fmbrkt.f).
    slpine = variant isa CentralStates ? 3 : 5
    if sp == slpine
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

"CS fire-mortality group (1–6) for species `sp` (cs/fmeff.f:223-236)."
@inline function cs_fire_mortality_group(sp::Integer)::Int
    sp == 47 || sp == 59               ? 1 :   # white oak, chestnut oak
    48 <= sp <= 51                     ? 2 :   # scarlet/black/northern+southern red oak
    14 <= sp <= 23                     ? 3 :   # hickories
    sp == 29                           ? 4 :   # red maple
    sp == 11 || sp == 13               ? 5 :   # black & swamp tupelo
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
                             flame::Float32, csv::Float32, variant::AbstractVariant = Southern())::Float32
    # The SN/CS Regelbrugge-Smith species groups (1-5) are gated `IF VARACD .EQ. 'SN'/'CS'`
    # (fmeff.f:196); NE skips them and uses the base Reinhardt logistic for every species.
    g = variant isa Northeast ? 6 :
        variant isa CentralStates ? cs_fire_mortality_group(sp) : fire_mortality_group(sp)
    if 1 <= g <= 5
        charht = flame * 0.7f0                          # max (uphill) char height
        xm = -(_FM_MORTB0[g] + _FM_MORTB1[g] * dbh * 2.54f0 + _FM_MORTB2[g] * charht / 3.28f0)
        mnmort = log(1f0 / 0.000001f0 - 1f0)            # guard against exp() overflow
        return xm >= mnmort ? 0f0 : 1f0 / (1f0 + exp(xm))
    else                                                # Reinhardt crown-scorch + bark
        bt = fire_bark_thickness(coef, sp, dbh, variant)
        xm = exp(-1.941f0 + 6.316f0 * (1f0 - exp(-bt)) - 0.000535f0 * csv * csv)
        return 1f0 / (1f0 + xm)
    end
end

"""
    fire_mortality_adjust(pmort, sp, dbh, burnseas) -> Float32

Variant-specific post-logistic mortality adjustments (FMEFF, fmeff.f:278-326). **For SN
(and CS) these do not apply** — the maple-<4″, hardwood-≤1″ and early-season (`burnseas`
≤ 2) reductions are gated by `IF (VARACD .EQ. 'LS'/'ON'/'NE')` and are skipped in the SN
variant. The only universal SN rule is `dbh ≤ 1″ & csv > 50% ⇒ 1.0` (fmeff.f:330), which
needs the scorch volume and is applied in `fmburn!`. This is a no-op kept as the seam
where a future LS/NE/ON port would re-introduce those branches.
"""
@inline function fire_mortality_adjust(pmort::Float32, sp::Integer, dbh::Float32, burnseas::Integer,
                                       variant::AbstractVariant = Southern())::Float32
    variant isa Northeast || return pmort                # SN/CS: no post-logistic adjustment
    # NE dormant-season (BURNSEAS≤2) reductions (ne/fmeff.f:304-326):
    (burnseas <= 2 && sp <= 25) && (pmort /= 2f0)        # conifers ×½ before greenup
    sp == 1 && (pmort = max(0.7f0, pmort))               # balsam fir floor 70%
    ((26 <= sp <= 29) || (99 <= sp <= 100)) && dbh < 4f0 && (pmort = 1f0)   # small maples die
    if burnseas <= 2 && sp > 25                          # hardwoods before greenup
        if (55 <= sp <= 70) || sp == 89                  # oaks — especially resistant
            pmort = dbh >= 2.5f0 ? pmort / 2f0 : pmort * 0.8f0
        else
            pmort *= 0.8f0
        end
    end
    (sp > 25 && dbh <= 1f0) && (pmort = 1f0)             # hardwoods ≤1″ die
    return pmort
end
