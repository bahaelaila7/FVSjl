# =============================================================================
# BM forest-grown crown width CRWDTH — base/cwidth.f → cwcalc.f (IWHO=0), faithful.
#
# cwcalc.f is one national routine: BM picks the 5-char equation code from BMMAP (cwcalc.f:102-106) and, because
# every BM location is a Region-6 forest (forkod.f KODFOR ∈ {604,607,614,616} after the 619→616 / 8117→614 remap),
# applies the R6 forest bias factor BF(KODFOR, FIASP=CWEQN(1:3)) (cwcalc.f:468-876) to the leading coefficient of
# the Crookston-R6 model-2 ('…05') forms. The equation forms themselves are the shared `_cwcalc_national` library.
#
# Replaces the old `bm_cwcalc(…; forest_bf=true)` TreeList path, which re-used CR's equation carriers (nearest-genus
# stand-ins for BM-unique WJ/PY/YC/CW/OH codes) and baked in only forest 614's BF.
# =============================================================================

# cwcalc.f:102-106 BMMAP — WP WL DF GF MH WJ LP ES AF PP WB LM PY YC AS CW OS OH
const _BM_CWMAP = ("11905", "07303", "20205", "01703", "26403", "06405",
                   "10805", "09305", "01905", "12205", "10105", "11301",
                   "23104", "04205", "74605", "74705", "12205", "31206")

# cwcalc.f SELECT CASE (KODFOR) → SELECT CASE (FIASP) BF table, the four BM forests (604 MALHEUR :554, 607 OCHOCO
# :613, 614 UMATILLA :699, 616 WALLOWA-WHITMAN :751). Any (forest, FIASP) pair not listed keeps BF=1.0.
const _BM_CW_R6BF = Dict{Int,Dict{String,Float32}}(
    604 => Dict("019" => 1.110f0, "073" => 0.818f0, "108" => 1.196f0, "093" => 1.121f0, "119" => 1.081f0,
                "202" => 1.058f0),
    607 => Dict("019" => 1.110f0, "073" => 0.879f0, "093" => 1.169f0, "108" => 1.196f0, "202" => 1.055f0),
    614 => Dict("017" => 1.076f0, "019" => 1.110f0, "073" => 0.907f0, "108" => 1.244f0, "093" => 1.137f0,
                "117" => 1.097f0, "119" => 1.128f0, "122" => 1.035f0, "202" => 1.055f0, "242" => 1.055f0,
                "263" => 1.106f0),
    616 => Dict("073" => 0.818f0, "108" => 1.114f0, "093" => 1.070f0, "264" => 1.077f0))

# cwcalc.f:477 — BF only for 601 ≤ KODFOR < 1000 (all BM forests); KODFOR is the forkod-remapped location code.
@inline function bm_cw_bf(kodfor::Int, fiasp::AbstractString)::Float32
    (kodfor < 601 || kodfor >= 1000) && return 1f0
    tbl = get(_BM_CW_R6BF, kodfor, nothing)
    tbl === nothing && return 1f0
    return get(tbl, fiasp, 1f0)
end

"""
    bm_crwdth(sp, d, h, cr, barea, el, hi, kodfor) -> Float32

BM CRWDTH (cwcalc.f IWHO=0): BMMAP equation code + R6 forest BF for the forkod-remapped `kodfor`, evaluated by the
national cwcalc library (which applies the final [0.5, 99.9] clamp). `cr` is the crown-ratio PERCENT FVS passes
(ICR for live records; the CRDUM=1 dummy for new sprouts, esuckr.f:314).
"""
function bm_crwdth(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32, el::Float32, hi::Float32,
                   kodfor::Int)::Float32
    (1 <= sp <= 18) || return 0f0
    eqn = _BM_CWMAP[sp]
    return _cwcalc_national(eqn, d, h, cr, barea, el, hi; bf = bm_cw_bf(kodfor, eqn[1:3]))
end
