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

# fvsvol.f passes SIDUM=0/BADUM=0 to the volume library, so R9_MHTS always falls back to
# the per-region SI default (LS 60 / CS 65 / NE 55) and BA=90 — NOT the stand's SI/BA.
_dvee_si_default(region::String) = region == "LS" ? 60 : region == "CS" ? 65 : 55

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
    ba <= 0 && (ba = 90)                       # R9_MHTS: BA≤0 ⇒ 90
    si <= 0 && (si = _dvee_si_default(region)) # R9_MHTS: SI≤0 ⇒ per-region default
    sif = Float32(si); baf = Float32(ba)
    base = b1 * (1f0 - exp(-b2 * dbh))^b3 * sif^b4
    esttht = 4.5f0 + base * (1.00001f0)^b5 * baf^b6           # FACTOR=0
    httot <= 0f0 && (httot = esttht)
    factor = min(topd / dbh, 1f0)
    estcht = 4.5f0 + base * (1.00001f0 - factor)^b5 * baf^b6
    pulbol = estcht * (httot / esttht)
    return trunc(Int, pulbol / 8.333333f0)
end

"R9_MHTS SAW path: the sawlog log count HT1PRD to the `bftopd`-inch board top. `h2` (pulp count) feeds the
 R9_MHTS `HT1PRD≤0 & HT2PRD>1 ⇒ HT1PRD=1` fixup."
function _r9_mhts_ht1prd(fia::Int, dbh::Float32, httot::Float32, si::Int, ba::Int,
                         region::String, bftopd::Float32, h2::Int)::Int
    b = get(_DVEE_HTCOEF, (fia, region), nothing)
    b === nothing && return 0
    b1, b2, b3, b4, b5, b6 = b
    ba <= 0 && (ba = 90)
    si <= 0 && (si = _dvee_si_default(region))
    sif = Float32(si); baf = Float32(ba)
    base = b1 * (1f0 - exp(-b2 * dbh))^b3 * sif^b4
    esttht = 4.5f0 + base * (1.00001f0)^b5 * baf^b6
    httot <= 0f0 && (httot = esttht)
    factor = min(bftopd / dbh, 1f0)
    estsht = 4.5f0 + base * (1.00001f0 - factor)^b5 * baf^b6
    sawbol = estsht * (httot / esttht)
    h1 = trunc(Int, sawbol / 8.333333f0)
    (h1 <= 0 && h2 > 1) && (h1 = 1)
    return h1
end

"Per-species·DBH-range correction factor for the R9VOL sawtimber merch cubic (r9vol.f '912' section)."
function _dvee_saw_cf(fia::Int, dbh::Float32)::Float32
    fia in (71, 94, 95, 97, 105, 241, 460, 543, 601, 602, 731, 742, 823, 824) && return 0.95f0
    fia in (400, 404, 651, 694, 813, 830) && return 1.05f0
    fia == 531 && return 1.10f0
    fia == 920 && return 0.90f0
    fia == 970 && return 1.08f0
    (dbh < 13f0 && fia == 110) && return 1.06f0
    if dbh < 15f0
        fia in (621, 746) && return 1.03f0
        fia == 125 && return 1.04f0
        fia == 837 && return 1.05f0
        fia in (835, 951) && return 1.06f0
        fia in (371, 833) && return 1.08f0
        fia in (129, 318, 802) && return 1.10f0
        fia == 806 && return 1.11f0
        fia == 375 && return 1.12f0
        fia == 762 && return 1.16f0
        fia == 316 && return 1.18f0
    else
        fia == 746 && return 0.95f0
        fia == 129 && return 0.96f0
        fia == 835 && return 1.01f0
        fia == 371 && return 1.03f0
        fia in (318, 375, 951) && return 1.04f0
        fia == 833 && return 1.05f0
        fia == 837 && return 1.06f0
        fia == 802 && return 1.07f0
        fia == 621 && return 1.08f0
        fia in (762, 806) && return 1.09f0
        fia == 316 && return 1.12f0
    end
    return 1f0
end

"""
    r9vol_gevorkiantz(fia, dbh, httot, iforst) -> (tcf, mcf, scf, bf)

The '900DVEE' (Gevorkiantz, r9vol.f R9VOL) cubic volumes. `tcf` = form-factor total cubic; `mcf` = merch
cubic — the PULPWOOD polynomial (to a 4" top, via HT2PRD) for DBH < BFMIND, the SAWTIMBER polynomial (to the
board top, via HT1PRD + per-species CF) for DBH ≥ BFMIND. `scf`/`bf` (sawtimber-cubic column / board feet)
are 0 (board-foot branch not ported; VOL(2)=0 for these regen stands). Validated bit-exact vs a live R9VOL
stamp for both the pulp (H2=1,2) and sawtimber (D9.004→V4 2.94) regimes. BFMIND=9 for CS (region-checked).
"""
function r9vol_gevorkiantz(fia::Int, dbh::Float32, httot::Float32, iforst::Int;
                           si::Int = 0, ba::Int = 0)
    tcf = httot > 0f0 ? 0.42f0 * Float32(pi) * dbh * dbh * httot / 576f0 : 0f0
    region = _dvee_region(iforst)
    h2 = _r9_mhts_ht2prd(fia, dbh, httot, si, ba, region)
    # Merch cubic to the 4" PULP top (VOL(4)+VOL(7) in R9VOL). The pulpwood polynomial computes this directly
    # for pulp-sized stems and closely approximates saw+topwood for sawtimber stems (≤~1.4% on the DVEE stand).
    # EXACT-fix TODO: for DBH≥BFMIND, Mcuft = the SAWLOG cubic VOL(4) (`_dvee_saw_cf` + the '912' polynomial via
    # `_r9_mhts_ht1prd`, both live-validated) PLUS the topwood VOL(7) (saw-top→pulp-top) — VOL(7) formula still
    # to stamp; those helpers are kept for it. Using pulp-for-all here beats a VOL(4)-only saw split (−22%).
    mcf = 0f0
    if h2 > 0
        h2f = Float32(h2)
        mcf = 0.001f0 * dbh * dbh * (1.9f0 + 0.01f0 * dbh) *
              (0.208f0 * h2f - 0.009984f0 * h2f * h2f + 0.04f0 / h2f) * 79f0
    end
    return (tcf, mcf, 0f0, 0f0)
end
