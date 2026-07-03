# =============================================================================
# r9vol_gevorkiantz.jl — the R9 Gevorkiantz volume model ('900DVEE', r9vol.f R9VOL)
#
# Ported from: bin/FVScs_buildDir/{r9vol.f R9VOL, r9init.f R9_MHTS + R9INIT}.
#
# This is the volume model FVS selects for CS species whose cubic-method code
# METHC==5 — set by the VOLUME keyword's field 7, e.g. `VOLUME  0  All … 5`. jl's
# default CS path (`r9clark_cubic`, the Clark taper) is METHC 6/9; DVEE is a
# distinct FORM-FACTOR model. Without it, a `VOLUME …5` stand (e.g. cst01_method5's
# BARE-GROUND-PLANT of shortleaf pine / bitternut hickory) diverges up to ~40% Mcuft.
#
# Formulas VALIDATED bit-exact vs a live R9VOL debug-stamp (see DIVERGENCE_FIX_CAMPAIGN
# D35): total cubic and the pulpwood merch cubic reproduce live exactly; the pulp-log
# count HT2PRD comes from R9_MHTS (a per-species Chapman-Richards height model whose
# B(1..6) coefficients were extracted into data/centralstates/dvee_r9_height_coef.csv).
#
# SCOPE (current): the pulpwood path (PROD='02', DBH < board-min ~9") — total cubic
# VOL(1) + pulp merch cubic VOL(4). Board feet VOL(2) is 0 across the sub-9" range and
# the sawtimber VOL(4) branch (PROD='01', r9vol.f '912' section with per-species CFs)
# is NOT yet ported — extend when a DVEE stand grows stems past the board minimum.
# =============================================================================

# R9INIT height coefficients B(1..6) per (FIA species, region∈{"LS","CS","NE"}).
const _DVEE_HTCOEF = let
    path = joinpath(@__DIR__, "..", "..", "data", "centralstates", "dvee_r9_height_coef.csv")
    d = Dict{Tuple{Int,String},NTuple{6,Float32}}()
    for (k, line) in enumerate(eachline(path))
        k == 1 && continue
        f = split(strip(line), ',')
        (isempty(f[1]) || length(f) < 8) && continue
        g(i) = parse(Float32, f[i])
        d[(parse(Int, f[1]), String(f[2]))] = (g(3), g(4), g(5), g(6), g(7), g(8))
    end
    d
end

"IFORST (CS forest = KODFOR−900) → R9_MHTS region flag: LS / CS / NE (r9init.f:13-24)."
function _dvee_region(iforst::Int)::String
    iforst in (2, 3, 4, 6, 7, 9, 10) && return "LS"
    iforst in (5, 8, 12)             && return "CS"
    return "NE"
end

"""
    _r9_mhts_ht2prd(fia, dbh, httot, si, ba, region; topd=4) -> Int

R9_MHTS (r9init.f): the pulpwood log count HT2PRD from total height. Estimates total
height ESTTHT and the height to the `topd`-inch pulp top ESTCHT via the per-species
Chapman-Richards curve, scales the bole by HTTOT/ESTTHT, and floors to 8⅓-ft logs.
Returns 0 when the species has no coefficients.
"""
function _r9_mhts_ht2prd(fia::Int, dbh::Float32, httot::Float32, si::Int, ba::Int,
                         region::String; topd::Float32 = 4f0)::Int
    b = get(_DVEE_HTCOEF, (fia, region), nothing)
    b === nothing && return 0
    b1, b2, b3, b4, b5, b6 = b
    ba <= 0 && (ba = 90)
    sif = Float32(si); baf = Float32(ba)
    base = b1 * (1f0 - exp(-b2 * dbh))^b3 * sif^b4
    esttht = 4.5f0 + base * (1.00001f0)^b5 * baf^b6           # FACTOR=0
    httot <= 0f0 && (httot = esttht)
    factor = min(topd / dbh, 1f0)
    estcht = 4.5f0 + base * (1.00001f0 - factor)^b5 * baf^b6
    pulbol = estcht * (httot / esttht)
    return trunc(Int, pulbol / 8.333333f0)
end

"""
    r9vol_gevorkiantz(fia, dbh, httot, si, ba, iforst) -> (tcf, mcf, scf, bf)

The '900DVEE' (Gevorkiantz, r9vol.f R9VOL) cubic volumes for a pulpwood-sized stem.
`tcf` = form-factor total cubic; `mcf` = pulp merch cubic to a 4" top. `scf`/`bf`
(sawtimber cubic / board feet) are returned 0 — the sawtimber & board branches are
not yet ported (0 across the sub-board-min range this model currently covers).
"""
function r9vol_gevorkiantz(fia::Int, dbh::Float32, httot::Float32, si::Int, ba::Int,
                           iforst::Int)
    tcf = httot > 0f0 ? 0.42f0 * Float32(pi) * dbh * dbh * httot / 576f0 : 0f0
    region = _dvee_region(iforst)
    h2 = _r9_mhts_ht2prd(fia, dbh, httot, si, ba, region)
    mcf = 0f0
    if h2 > 0
        h2f = Float32(h2)
        mcf = 0.001f0 * dbh * dbh * (1.9f0 + 0.01f0 * dbh) *
              (0.208f0 * h2f - 0.009984f0 * h2f * h2f + 0.04f0 / h2f) * 79f0
    end
    return (tcf, mcf, 0f0, 0f0)
end
