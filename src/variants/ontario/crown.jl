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
# ON species JSP2 US-code → open-grown CWEQ (cwcalc.f IWHO=1 SELECT CASE, all ON species).
const ON_CW_OPEN = Dict{String,String}(
    "AB" => "53101",
    "AC" => "40703",
    "AE" => "97203",
    "AH" => "39101",
    "AP" => "76102",
    "BA" => "54301",
    "BC" => "76203",
    "BE" => "31301",
    "BF" => "01203",
    "BG" => "69301",
    "BH" => "40201",
    "BK" => "90101",
    "BL" => "97203",
    "BM" => "31803",
    "BN" => "60201",
    "BO" => "83704",
    "BP" => "74101",
    "BR" => "82303",
    "BS" => "09503",
    "BT" => "74301",
    "BW" => "95101",
    "CC" => "76102",
    "CH" => "60201",
    "CK" => "82601",
    "CW" => "74203",
    "DM" => "97203",
    "DW" => "49101",
    "EH" => "26101",
    "GA" => "54403",
    "HB" => "46201",
    "HH" => "70101",
    "HT" => "49101",
    "JP" => "10503",
    "MA" => "55201",
    "MM" => "31301",
    "NC" => "49101",
    "NP" => "80901",
    "NS" => "09104",
    "OS" => "06801",
    "PB" => "37503",
    "PH" => "40301",
    "PL" => "76102",
    "PR" => "76102",
    "QA" => "74603",
    "RC" => "06801",
    "RE" => "97203",
    "RL" => "97501",
    "RM" => "31603",
    "RN" => "12503",
    "RO" => "83303",
    "RP" => "12503",
    "SC" => "13001",
    "SH" => "40703",
    "SM" => "31803",
    "SS" => "93101",
    "ST" => "31301",
    "SV" => "31701",
    "SW" => "80204",
    "SY" => "73101",
    "TA" => "07103",
    "WA" => "54101",
    "WC" => "24101",
    "WI" => "97203",
    "WN" => "60201",
    "WO" => "80204",
    "WP" => "12903",
    "WS" => "09403",
    "YB" => "37101",
)


# CWEQ → equation (cwcalc.f coefficients). ek: a+b·D^power (thr 3", clamp max_cw).
# bechtold: a+b·Dc+c·Dc²+cr_coef·CR+hi_coef·HI (thr 5", Dc=min(D,dbh_cap), clamp).
# ON open-grown crown-width equations (cwcalc.f IWHO=1), coefficients from the eastern
# crown-width library (data/lakestates/crown_width_equations.csv — same cwcalc.f source).
const ON_CW_EQS = Dict{String,CrownWidthEq}(
    "01203" => CrownWidthEq(:ek, 0.327f0, 5.116f0, 0f0, 0f0, 0f0, 0.5035f0, 999f0, 34.0f0),
    "06801" => CrownWidthEq(:bechtold, 1.2359f0, 1.2962f0, 0.0f0, 0.0545f0, 0.0f0, 0f0, 999.0f0, 33.0f0),
    "07103" => CrownWidthEq(:ek, 2.205f0, 3.475f0, 0f0, 0f0, 0f0, 0.7506f0, 999f0, 29.0f0),
    "09104" => CrownWidthEq(:ek, 5.057f0, 1.1313f0, 0f0, 0f0, 0f0, 1.0f0, 999f0, 47.0f0),
    "09403" => CrownWidthEq(:ek, 3.594f0, 1.963f0, 0f0, 0f0, 0f0, 0.882f0, 999f0, 37.0f0),
    "09503" => CrownWidthEq(:ek, 3.655f0, 1.398f0, 0f0, 0f0, 0f0, 1.0f0, 999f0, 27.0f0),
    "10503" => CrownWidthEq(:ek, 0.299f0, 5.644f0, 0f0, 0f0, 0f0, 0.6036f0, 999f0, 30.0f0),
    "12503" => CrownWidthEq(:ek, 4.233f0, 1.462f0, 0f0, 0f0, 0f0, 1.0f0, 999f0, 39.0f0),
    "12903" => CrownWidthEq(:ek, 1.62f0, 3.197f0, 0f0, 0f0, 0f0, 0.7981f0, 999f0, 58.0f0),
    "13001" => CrownWidthEq(:bechtold, 3.5522f0, 0.6742f0, 0.0f0, 0.0985f0, 0.0f0, 0f0, 999.0f0, 27.0f0),
    "24101" => CrownWidthEq(:bechtold, -0.0634f0, 0.7057f0, 0.0f0, 0.0837f0, 0.0f0, 0f0, 999.0f0, 27.0f0),
    "26101" => CrownWidthEq(:bechtold, 6.1924f0, 1.4491f0, -0.0178f0, 0.0f0, -0.0341f0, 0f0, 40.0f0, 0.0f0),
    "31301" => CrownWidthEq(:bechtold, 6.4741f0, 1.0778f0, 0.0f0, 0.0719f0, -0.0637f0, 0f0, 999.0f0, 57.0f0),
    "31603" => CrownWidthEq(:ek, 0.0f0, 4.776f0, 0f0, 0f0, 0f0, 0.7656f0, 999f0, 55.0f0),
    "31701" => CrownWidthEq(:bechtold, 3.3576f0, 1.1312f0, 0.0f0, 0.1011f0, -0.173f0, 0f0, 999.0f0, 45.0f0),
    "31803" => CrownWidthEq(:ek, 0.868f0, 4.15f0, 0f0, 0f0, 0f0, 0.7514f0, 999f0, 54.0f0),
    "37101" => CrownWidthEq(:bechtold, -1.1151f0, 2.2888f0, -0.0493f0, 0.0985f0, -0.0396f0, 0f0, 24.0f0, 0.0f0),
    "37503" => CrownWidthEq(:ek, 3.639f0, 1.953f0, 0f0, 0f0, 0f0, 1.0f0, 999f0, 42.0f0),
    "39101" => CrownWidthEq(:bechtold, 0.9219f0, 1.6303f0, 0.0f0, 0.115f0, -0.1113f0, 0f0, 999.0f0, 42.0f0),
    "40201" => CrownWidthEq(:bechtold, 8.0118f0, 1.4212f0, 0.0f0, 0.0f0, 0.0f0, 0f0, 999.0f0, 41.0f0),
    "40301" => CrownWidthEq(:bechtold, 3.9234f0, 1.522f0, 0.0f0, 0.0405f0, 0.0f0, 0f0, 999.0f0, 53.0f0),
    "40703" => CrownWidthEq(:ek, 2.36f0, 3.548f0, 0f0, 0f0, 0f0, 0.7986f0, 999f0, 54.0f0),
    "46201" => CrownWidthEq(:bechtold, 7.1043f0, 1.3041f0, 0.0f0, 0.0456f0, 0.0f0, 0f0, 999.0f0, 51.0f0),
    "49101" => CrownWidthEq(:bechtold, 2.9646f0, 1.9917f0, 0.0f0, 0.0707f0, 0.0f0, 0f0, 999.0f0, 36.0f0),
    "53101" => CrownWidthEq(:bechtold, 3.9361f0, 1.15f0, 0.0f0, 0.1237f0, -0.0691f0, 0f0, 999.0f0, 80.0f0),
    "54101" => CrownWidthEq(:bechtold, 1.7625f0, 1.3413f0, 0.0f0, 0.0957f0, 0.0f0, 0f0, 999.0f0, 62.0f0),
    "54301" => CrownWidthEq(:bechtold, 5.2824f0, 1.1184f0, 0.0f0, 0.0f0, 0.0f0, 0f0, 999.0f0, 34.0f0),
    "54403" => CrownWidthEq(:ek, 0.0f0, 4.755f0, 0f0, 0f0, 0f0, 0.7381f0, 999f0, 61.0f0),
    "55201" => CrownWidthEq(:bechtold, 4.1971f0, 1.5567f0, 0.0f0, 0.088f0, 0.0f0, 0f0, 999.0f0, 46.0f0),
    "60201" => CrownWidthEq(:bechtold, 3.6031f0, 1.1472f0, 0.0f0, 0.1224f0, 0.0f0, 0f0, 999.0f0, 37.0f0),
    "69301" => CrownWidthEq(:bechtold, 5.5037f0, 1.0567f0, 0.0f0, 0.088f0, 0.061f0, 0f0, 999.0f0, 50.0f0),
    "70101" => CrownWidthEq(:bechtold, 7.8084f0, 0.8129f0, 0.0f0, 0.0941f0, -0.0817f0, 0f0, 999.0f0, 39.0f0),
    "73101" => CrownWidthEq(:bechtold, -1.3973f0, 1.3756f0, 0.0f0, 0.1835f0, 0.0f0, 0f0, 999.0f0, 66.0f0),
    "74101" => CrownWidthEq(:bechtold, 6.2498f0, 0.8655f0, 0.0f0, 0.0f0, 0.0f0, 0f0, 999.0f0, 25.0f0),
    "74203" => CrownWidthEq(:ek, 2.934f0, 2.538f0, 0f0, 0f0, 0f0, 0.8617f0, 999f0, 80.0f0),
    "74301" => CrownWidthEq(:bechtold, 0.6847f0, 1.105f0, 0.0f0, 0.142f0, -0.0265f0, 0f0, 999.0f0, 43.0f0),
    "74603" => CrownWidthEq(:ek, 4.203f0, 2.129f0, 0f0, 0f0, 0f0, 1.0f0, 999f0, 43.0f0),
    "76102" => CrownWidthEq(:braggm, 4.102718f0, 1.396006f0, 0f0, 0f0, 0f0, 1.077474f0, 999f0, 52.0f0),
    "76203" => CrownWidthEq(:ek, 0.621f0, 7.059f0, 0f0, 0f0, 0f0, 0.5441f0, 999f0, 52.0f0),
    "80204" => CrownWidthEq(:ek, 1.8f0, 1.883f0, 0f0, 0f0, 0f0, 1.0f0, 999f0, 69.0f0),
    "80901" => CrownWidthEq(:bechtold, 4.8935f0, 1.6069f0, 0.0f0, 0.0f0, 0.0f0, 0f0, 999.0f0, 44.0f0),
    "82303" => CrownWidthEq(:ek, 0.942f0, 3.539f0, 0f0, 0f0, 0f0, 0.7952f0, 999f0, 78.0f0),
    "82601" => CrownWidthEq(:bechtold, 0.5189f0, 1.4134f0, 0.0f0, 0.1365f0, -0.0806f0, 0f0, 999.0f0, 45.0f0),
    "83303" => CrownWidthEq(:ek, 2.85f0, 3.782f0, 0f0, 0f0, 0f0, 0.7968f0, 999f0, 82.0f0),
    "83704" => CrownWidthEq(:ek, 4.51f0, 1.67f0, 0f0, 0f0, 0f0, 1.0f0, 999f0, 52.0f0),
    "90101" => CrownWidthEq(:bechtold, 3.0012f0, 0.8165f0, 0.0f0, 0.1395f0, 0.0f0, 0f0, 999.0f0, 48.0f0),
    "93101" => CrownWidthEq(:bechtold, 4.6311f0, 1.0108f0, 0.0f0, 0.0564f0, 0.0f0, 0f0, 999.0f0, 29.0f0),
    "95101" => CrownWidthEq(:bechtold, 1.6871f0, 1.211f0, 0.0f0, 0.1194f0, -0.0264f0, 0f0, 999.0f0, 61.0f0),
    "97203" => CrownWidthEq(:ek, 2.829f0, 3.456f0, 0f0, 0f0, 0f0, 0.8575f0, 999f0, 72.0f0),
    "97501" => CrownWidthEq(:bechtold, 9.0023f0, 1.3933f0, 0.0f0, 0.0f0, -0.0785f0, 0f0, 999.0f0, 49.0f0),
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
