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
# CR fire bark: DBH·B1[sp] with a CR-SPECIFIC per-species B1 (FOFEM V5.0, cr/fmbrkt.f) — NOT the eastern
# EQNUM-indexed table (CR has no bark_eqnum). Thin CR bark (0.025-0.063) drives the fireline-intensity kill.
const _CR_FM_BARK_B1 = Float32[
    0.041, 0.041, 0.063, 0.046, 0.048, 0.040, 0.035, 0.063, 0.030, 0.030,
    0.028, 0.030, 0.063, 0.030, 0.035, 0.025, 0.031, 0.036, 0.025, 0.044,
    0.038, 0.038, 0.045, 0.045, 0.045, 0.045, 0.045, 0.044, 0.025, 0.025,
    0.025, 0.025, 0.030, 0.030, 0.030, 0.063, 0.030, 0.038]

# N-Rockies western fire bark: plain DBH·B1[sp], per-species B1 (ie/fmbrkt.f, no special cases). KT's
# 11 species are the first 11 of IE's 23 (KT⊂IE), so KT reuses this table truncated. EM/BM/UT/TT/CI add
# their own western fmbrkt tables as their FFE lands (see docs/EXTENSIONS_ROLLOUT_PLAN.md).
const _IE_FM_BARK_B1 = Float32[
    0.035, 0.063, 0.063, 0.046, 0.040, 0.035, 0.028, 0.036, 0.041, 0.063, 0.040,
    0.030, 0.030, 0.050, 0.030, 0.025, 0.025, 0.044, 0.038, 0.040, 0.027, 0.026, 0.040]
# BM fire bark (bm/fmbrkt.f, FOFEM V5.0 Reinhardt) — plain DBH·B1[sp], 18 species.
const _BM_FM_BARK_B1 = Float32[
    0.035, 0.063, 0.063, 0.046, 0.040, 0.025, 0.028, 0.036, 0.041, 0.063,
    0.030, 0.030, 0.025, 0.022, 0.044, 0.044, 0.063, 0.044]
# NC (Klamath) fire bark-thickness B1 per species (nc/fmbrkt.f, FOFEM v5.0). NC previously fell to the SN
# _FM_BARK_B1[bark_eqnum=1]=0.019 for EVERY species — ~3.3× too thin for the conifers (DF 0.063) ⇒ cambium
# over-kill of the 10-20" overstory once the crown-fire flame/scorch was corrected. `bt = DBH·B1` (fmbrkt.f).
const _NC_FM_BARK_B1 = Float32[
    0.063, 0.072, 0.063, 0.048, 0.060, 0.060, 0.030, 0.052, 0.039, 0.063, 0.052, 0.081]
# WS (WestSierra) fire bark-thickness B1 per species (ws/fmbrkt.f, FOFEM v5.0), 43 species. Same as NC's fix:
# WS previously fell to the SN _FM_BARK_B1[bark_eqnum=1]=0.019 for EVERY species — ~2-4× too thin for the WS
# conifers (SP 0.072/DF 0.063/WF 0.048/RF 0.039) ⇒ cambium over-kill of the large overstory (wst01+SIMFIRE
# post-fire BA kill 142 vs oracle 49) once the fuel-model selection was corrected. `bt = DBH·B1` (fmbrkt.f).
const _WS_FM_BARK_B1 = Float32[
    0.072, 0.063, 0.048, 0.081, 0.060, 0.068, 0.039, 0.063, 0.028, 0.030,
    0.035, 0.030, 0.047, 0.030, 0.030, 0.063, 0.030, 0.030, 0.033, 0.030,
    0.030, 0.063, 0.081, 0.040, 0.025, 0.025, 0.025, 0.050, 0.024, 0.033,
    0.030, 0.043, 0.034, 0.052, 0.045, 0.044, 0.026, 0.060, 0.062, 0.024,
    0.044, 0.028, 0.030]
# CA (CentralCalifornia) fire bark-thickness B1 per species (ca/fmbrkt.f, FOFEM v5.0), 50 species. Same class
# as WS/NC: CA fell to the SN _FM_BARK_B1[bark_eqnum=1]=0.019 default ⇒ fire mortality wrong. `bt = DBH·B1`.
const _CA_FM_BARK_B1 = Float32[
    0.081, 0.06, 0.035, 0.048, 0.039, 0.039, 0.063, 0.035, 0.04, 0.03,
    0.03, 0.028, 0.063, 0.03, 0.068, 0.072, 0.035, 0.063, 0.03, 0.033,
    0.025, 0.025, 0.081, 0.025, 0.063, 0.05, 0.024, 0.033, 0.059, 0.029,
    0.03, 0.043, 0.034, 0.024, 0.036, 0.026, 0.06, 0.045, 0.062, 0.042,
    0.041, 0.052, 0.033, 0.044, 0.044, 0.041, 0.025, 0.026, 0.03, 0.081]

# WC (WestCascades) fire bark-thickness B1 per species (wc/fmbrkt.f, FOFEM v5.0 Reinhardt), 39 species.
# Same class as WS/NC/CA: without this WC would fall to the SN _FM_BARK_B1[bark_eqnum=1]=0.019 default for
# every species ⇒ fire mortality wrong. `bt = DBH·B1` (fmbrkt.f). Species 6/38 are blank slots (B1=0).
const _WC_FM_BARK_B1 = Float32[
    0.047, 0.048, 0.046, 0.041, 0.039, 0.000, 0.045, 0.022, 0.081, 0.036,
    0.028, 0.068, 0.072, 0.035, 0.063, 0.063, 0.081, 0.035, 0.040, 0.040,
    0.024, 0.026, 0.060, 0.027, 0.045, 0.044, 0.044, 0.029, 0.025, 0.050,
    0.030, 0.030, 0.025, 0.062, 0.038, 0.062, 0.041, 0.000, 0.044]

# PN (PacificNorthwest) fire bark-thickness B1 per species (pn/fmbrkt.f, FOFEM v5.0 Reinhardt), 39 sp.
# IDENTICAL to WC except slot 6 = 0.027 (SS; WC blank -0.000). `bt = DBH·B1` (fmbrkt.f).
const _PN_FM_BARK_B1 = Float32[
    0.047, 0.048, 0.046, 0.041, 0.039, 0.027, 0.045, 0.022, 0.081, 0.036,
    0.028, 0.068, 0.072, 0.035, 0.063, 0.063, 0.081, 0.035, 0.040, 0.040,
    0.024, 0.026, 0.060, 0.027, 0.045, 0.044, 0.044, 0.029, 0.025, 0.050,
    0.030, 0.030, 0.025, 0.062, 0.038, 0.062, 0.041, 0.000, 0.044]

# EC (EastCascades) fire bark-thickness B1 per species (ec/fmbrkt.f, FOFEM v5.0 Reinhardt), 32 species.
# Same class as WS/NC/CA/WC/PN: without this EC would fall to the SN _FM_BARK_B1[bark_eqnum=1]=0.019 default
# for every species (~2-3× too thin for the conifers, DF 0.063) ⇒ fire mortality wrong. `bt = DBH·B1`.
const _EC_FM_BARK_B1 = Float32[
    0.035, 0.063, 0.063, 0.047, 0.035, 0.046, 0.028, 0.036, 0.041, 0.063,
    0.040, 0.040, 0.025, 0.030, 0.045, 0.046, 0.050, 0.022, 0.025, 0.024,
    0.024, 0.026, 0.027, 0.045, 0.062, 0.044, 0.044, 0.029, 0.062, 0.041,
    0.040, 0.044]

# SO (SouthCentralOregon) fire bark-thickness B1 per species (so/fmbrkt.f, FOFEM v5.0 Reinhardt), 33 species.
# Same class as WS/NC/CA/WC/PN/EC: without this SO would fall to the SN _FM_BARK_B1[bark_eqnum=1]=0.019 default
# for every species (~3× too thin for the conifers, DF 0.063) ⇒ fire mortality wrong. `bt = DBH·B1`.
const _SO_FM_BARK_B1 = Float32[
    0.035, 0.072, 0.063, 0.048, 0.040, 0.060, 0.028, 0.036, 0.039, 0.063,
    0.025, 0.046, 0.041, 0.047, 0.045, 0.030, 0.063, 0.035, 0.040, 0.025,
    0.062, 0.026, 0.024, 0.044, 0.044, 0.062, 0.029, 0.041, 0.045, 0.044,
    0.044, 0.063, 0.044]

@inline function fire_bark_thickness(coef::SpeciesCoefficients, sp::Integer, dbh::Float32,
                                     variant::AbstractVariant = Southern())::Float32
    variant isa CentralRockies && return dbh * _CR_FM_BARK_B1[Int(sp)]   # cr/fmbrkt.f
    # IE (23 sp) and KT (11 sp = first 11 of IE) share ie/fmbrkt.f's plain DBH·B1 form.
    (variant isa InlandEmpire || variant isa Kootenai) && return dbh * _IE_FM_BARK_B1[Int(sp)]
    variant isa BlueMountains && return dbh * _BM_FM_BARK_B1[Int(sp)]   # bm/fmbrkt.f
    variant isa Klamath && return dbh * _NC_FM_BARK_B1[Int(sp)]         # nc/fmbrkt.f
    variant isa WestSierra && return dbh * _WS_FM_BARK_B1[Int(sp)]      # ws/fmbrkt.f
    variant isa CentralCalifornia && return dbh * _CA_FM_BARK_B1[Int(sp)]  # ca/fmbrkt.f
    variant isa WestCascades && return dbh * _WC_FM_BARK_B1[Int(sp)]       # wc/fmbrkt.f
    variant isa PacificNorthwest && return dbh * _PN_FM_BARK_B1[Int(sp)]   # pn/fmbrkt.f
    variant isa EastCascades && return dbh * _EC_FM_BARK_B1[Int(sp)]       # ec/fmbrkt.f
    variant isa SouthCentralOregon && return dbh * _SO_FM_BARK_B1[Int(sp)] # so/fmbrkt.f
    # Shortleaf pine uses the Harmon (1984) quadratic INSTEAD of the B1 table — but ONLY in the variants
    # where it is a species: SN sp5 (sn/fmbrkt.f:126) and CS sp3 (cs/fmbrkt.f:133). NE and LS have NO such
    # special case (their fmbrkt.f is a plain DBH·B1[EQNUM] for every species) and sp5 there is NOT shortleaf
    # (LS sp5 = eastern white pine, EQNUM 24 → B1 .045) — applying Harmon there over-thickens the bark ~2× and
    # under-kills the tree in fire. (Was gated `: 5` for all non-CS, silently mis-barking NE sp5 too.)
    slpine = variant isa Southern ? 5 : variant isa CentralStates ? 3 : 0
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
    # fmburn.f:471-472 `SCH=(63/(140-ATEMP))*(BYRAM**(7/6)/(BYRAM+FWIND**3.0)**0.5)`. ALL THREE powers are FVS
    # `**` (gfortran powf): `**(7/6)`, `FWIND**3.0` (a REAL exponent ⇒ powf, NOT x*x*x — gfortran lowers them
    # differently ~26% of the time!), and the outer `**0.5` (powf, NOT sqrt — differ ~0.07%). Route all three
    # through the companion (doctrine #8) to match FVS bit-exactly. (The old fwind^3 + sqrt caused the scorch Δ0.002.)
    return (63f0 / (140f0 - atemp)) * (fpow(b, 7f0 / 6f0) / fpow(b + fpow(fwind, 3f0), 0.5f0))
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
    # The Regelbrugge-Smith groups are gated `IF VARACD .EQ. 'SN'/'CS'` — NE and LS/ON skip them
    # (ls/fmeff.f:196) and use the base Reinhardt logistic (group 6) for every species.
    # CR gates the Regelbrugge-Smith groups (1-5) to VARACD=='SN'/'CS' ONLY (cr/fmeff.f:196), so CR — like
    # NE/LS/ON — uses the base Reinhardt crown-scorch+bark logistic (group 6) for EVERY species. jl previously
    # fell CR through to fire_mortality_group(sp) (the SN species map), which mis-assigned CR sp20 (aspen) to
    # SN group 4 (red maple Regelbrugge-Smith) and sp27 to group 3 ⇒ the wrong DBH+char-height logistic
    # UNDER-killed large aspen (crt01 FFE 10-20" kill 9 vs live 17 ⇒ +9% surviving BA).
    # NC (like NE/LS/CR/BM) gates the Regelbrugge-Smith groups (1-5) to VARACD=='SN'/'CS' ONLY (nc/fmeff.f:196),
    # so it uses the base Reinhardt crown-scorch+bark logistic (group 6) for EVERY species. Without Klamath here it
    # fell through to fire_mortality_group(sp) (the SN species map) ⇒ NC large trees mis-assigned to Regelbrugge-
    # Smith DBH+char-height logistics ⇒ over-killed the 10-20" class (jl 41/58 vs live 19/58) once the flame/scorch
    # was corrected (the surface-only low scorch had masked it — the two-canceling-bugs case).
    # WS (like NE/LS/CR/BM/NC) gates the Regelbrugge-Smith groups (1-5) to VARACD=='SN'/'CS' ONLY
    # (ws/fmeff.f:196), so it uses the base Reinhardt crown-scorch+bark logistic (group 6) for EVERY species.
    # Without WestSierra here it fell through to fire_mortality_group(sp) (the SN species map) ⇒ WS species
    # mis-assigned to SN Regelbrugge-Smith DBH+char-height logistics ⇒ fire OVER-kill (wst01+SIMFIRE post-fire
    # TPA 112 vs oracle 218) once the fuel-model selection was corrected (the low-scorch surface fire exposed it).
    # WC (like NE/LS/CR/BM/NC/WS/CA) gates the Regelbrugge-Smith groups (1-5) to VARACD=='SN'/'CS' ONLY
    # (wc/fmeff.f), so it uses the base Reinhardt crown-scorch+bark logistic (group 6) for EVERY species.
    g = (variant isa Northeast || variant isa LakeStates || variant isa CentralRockies || variant isa BlueMountains || variant isa Klamath || variant isa WestSierra || variant isa CentralCalifornia || variant isa WestCascades || variant isa PacificNorthwest || variant isa EastCascades || variant isa SouthCentralOregon) ? 6 :
        variant isa CentralStates ? cs_fire_mortality_group(sp) : fire_mortality_group(sp)
    if 1 <= g <= 5
        charht = flame * 0.7f0                          # max (uphill) char height
        xm = -(_FM_MORTB0[g] + _FM_MORTB1[g] * dbh * 2.54f0 + _FM_MORTB2[g] * charht / 3.28f0)
        mnmort = flog(1f0 / 0.000001f0 - 1f0)           # guard against exp() overflow (ALOG → FFI companion)
        # fmeff.f mortality logistic EXP → FFI companion (gfortran), not native openlibm (doctrine #8): the
        # per-tree fire-kill probability feeds the RANN gate, so a libm-exp ULP flips a near-threshold kill.
        return xm >= mnmort ? 0f0 : 1f0 / (1f0 + fexp(xm))
    else                                                # Reinhardt crown-scorch + bark
        bt = fire_bark_thickness(coef, sp, dbh, variant)
        xm = fexp(-1.941f0 + 6.316f0 * (1f0 - fexp(-bt)) - 0.000535f0 * csv * csv)
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
    if variant isa LakeStates
        # LS/ON dormant-season reductions (ls/fmeff.f:278-300). BURNSEAS≤2 = before greenup.
        (burnseas <= 2 && sp <= 14) && (pmort /= 2f0)            # conifers ×½ before greenup
        sp == 8 && (pmort = max(0.7f0, pmort))                   # balsam fir floor 70%
        (sp in (18,19,26,27,51,52) && dbh < 4f0) && (pmort = 1f0)  # certain maples <4″ die
        if burnseas <= 2 && sp > 14                              # hardwoods before greenup
            if 30 <= sp <= 36                                    # oaks — especially resistant
                pmort = dbh >= 2.5f0 ? pmort / 2f0 : pmort * 0.8f0
            else
                pmort *= 0.8f0
            end
        end
        (sp > 14 && dbh <= 1f0) && (pmort = 1f0)                 # hardwoods ≤1″ die
        return pmort
    end
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
