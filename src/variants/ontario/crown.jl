# =============================================================================
# ontario/crown.jl — Ontario crown-ratio dub/update (canada/on/crown.f).
#
# ON's crown.f IS the shared TWIGS NC-125 model (GTR NC-125) used by NE/CS/LS:
#     CRNEW = 10·( BCR1/(1 + BCR2·BA) + BCR3·(1 − exp(BCR4·D)) )        (crown.f:180)
# with the same ±1%/yr change limit, CRNMLT window, crown-length cap, top-kill reduction and
# 1/10/95 bounds. BCR1..BCR4 are per species (data/ontario/species_coefficients.csv crown_bcr1..4,
# transcribed from canada/on/crown.f DATA via a dump of FVSon_g16). BCR4's sign is baked into the
# constant (ON uses exp(BCR4·D) directly), matching the shared kernel's `fexp(bcr4·d)`.
#
# VALIDATED: the LSTART crown dub reproduces FVSon_g16's cyc0 ICR 8/8 given the oracle's crown-time
# stand BA (scratchpad/on/check_crown_formula.jl). End-to-end cyc0 crown bit-exactness is gated on
# the metric plot/TPA expansion (stand BA/QMD), the documented next density blocker.
# =============================================================================

crown_ratio_update!(s::StandState, ::Ontario; kwargs...) = _twigs_crown_update!(s; kwargs...)

# =============================================================================
# ON open-grown crown WIDTH (canada/on/cwcalc.f, IWHO=1) → CCF (canada/on/ccfcal.f,
# the LS-form CCFCAL the ON buildDir links). cwcalc.f remaps ISPC → a US 2-char code
# (JSP2, blkdat-style DATA) and evaluates the eastern-US Bechtold/Ek crown-width
# library — the SAME four families `_cw_eval` already implements — then ccfcal does
# `CCFT = 0.001803·CW²` (D>0.1) or 0.001, with CR fixed at 90 and IWHO=1. The generic
# `stand_ccf` `else` path fed ON's raw code2 ("PW","SW",…) which never keys the
# US-coded library ⇒ cw→0.5 ⇒ CCF≈15; this remap fixes it.
#
# ON_JSP2 is cwcalc.f's ISPC→US-code table (lines 74-81), verbatim. ON_CW_OPEN is the
# IWHO=1 CWEQ each JSP2 code selects (cwcalc.f SELECT CASE). Only the species an ON
# stand exercises are transcribed here; species whose cwcalc CASE is commented-out (CW
# stays 0) — and any not yet transcribed — return 0, faithful to cwcalc's default.
# =============================================================================

const ON_JSP2 = String[
    "JP","SC","RN","RP","WP","WS","NS","BF","BS","TA","WC","EH","OS","RC","BA","GA",
    "CW","SV","RM","BC","AE","RL","RE","YB","BW","SM","BM","AB","WA","WO","SW","BR",
    "CK","RO","BO","NP","BH","PH","SH","BT","QA","BP","PB","CH","BN","WN","HH","BK",
    "NC","BE","ST","MM","AH","AC","HB","DW","HT","AP","BG","SY","PR","CC","PL","WI",
    "BL","DM","SS","MA","JP","WP","WS","BS",
]

# US 2-char code → open-grown (IWHO=1) CWEQ (cwcalc.f SELECT CASE, IWHO=1 branch).
const ON_CW_OPEN = Dict{String,String}(
    "JP" => "10503",   # jack pine  — Ek 1974
    "WP" => "12903",   # e. white pine — Ek 1974
    "WS" => "09403",   # white spruce — Ek 1974
    "BS" => "09503",   # black spruce — Ek 1974
    "BF" => "01203",   # balsam fir  — Ek 1974
    "SM" => "31803",   # sugar maple — Ek 1974
    "WC" => "24101",   # n. white cedar — Bechtold m2 (no IWHO branch: forest eqn used)
    "AB" => "53101",   # american beech — Bechtold m3 (no IWHO branch: forest eqn used)
)

# CWEQ → equation (cwcalc.f coefficients). ek: a+b·D^power (thr 3", clamp max_cw).
# bechtold: a+b·Dc+c·Dc²+cr_coef·CR+hi_coef·HI (thr 5", Dc=min(D,dbh_cap), clamp).
const ON_CW_EQS = Dict{String,CrownWidthEq}(
    "10503" => CrownWidthEq(:ek, 0.2990f0, 5.6440f0, 0f0, 0f0, 0f0, 0.6036f0, 0f0, 30f0),
    "12903" => CrownWidthEq(:ek, 1.6200f0, 3.1970f0, 0f0, 0f0, 0f0, 0.7981f0, 0f0, 58f0),
    "09403" => CrownWidthEq(:ek, 3.5940f0, 1.9630f0, 0f0, 0f0, 0f0, 0.8820f0, 0f0, 37f0),
    "09503" => CrownWidthEq(:ek, 3.6550f0, 1.3980f0, 0f0, 0f0, 0f0, 1.0000f0, 0f0, 27f0),
    "01203" => CrownWidthEq(:ek, 0.3270f0, 5.1160f0, 0f0, 0f0, 0f0, 0.5035f0, 0f0, 34f0),
    "31803" => CrownWidthEq(:ek, 0.8680f0, 4.1500f0, 0f0, 0f0, 0f0, 0.7514f0, 0f0, 54f0),
    # bechtold linear (no inner D cap): dbh_cap huge so min(D,cap) is inert.
    "24101" => CrownWidthEq(:bechtold, -0.0634f0, 0.7057f0, 0f0, 0.0837f0,  0f0,     0f0, 1f9, 27f0),
    "53101" => CrownWidthEq(:bechtold,  3.9361f0, 1.1500f0, 0f0, 0.1237f0, -0.0691f0, 0f0, 1f9, 80f0),
)

"""
    on_open_crown_width(ispc, d, lat, long, elev) -> Float32

Open-grown crown width (ft) for ON species index `ispc`, DBH `d` (in), evaluated at the
fixed CR=90 that ccfcal.f passes (IWHO=1). Species with no active cwcalc.f equation → 0.
"""
@inline function on_open_crown_width(ispc::Integer, d::Real, lat::Real, long::Real, elev::Real)::Float32
    (1 <= ispc <= length(ON_JSP2)) || return 0f0
    eqn = get(ON_CW_OPEN, ON_JSP2[ispc], "")
    eqn == "" && return 0f0
    e = ON_CW_EQS[eqn]
    cw = _cw_eval(e, Float32(d), 90f0, hopkins_index(lat, long, elev))
    cw < 0.5f0 && (cw = 0.5f0)      # cwcalc.f final clamp (inert for ont01's large trees)
    cw > 99.9f0 && (cw = 99.9f0)
    return cw
end

"Per-tree CCF for ON (ccfcal.f): 0.001803·CW² for D>0.1, else 0.001 (×P applied by caller)."
@inline function on_tree_ccf(ispc::Integer, d::Real, lat::Real, long::Real, elev::Real)::Float32
    df = Float32(d)
    df <= 0.1f0 && return 0.001f0
    cw = on_open_crown_width(ispc, df, lat, long, elev)
    return 0.001803f0 * cw * cw
end
