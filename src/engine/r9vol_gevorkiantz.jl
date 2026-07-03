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

"Per-species·IFORST correction factor for the R9VOL board-foot TABLE B (r9vol.f; IFORST∈{4,5,8,11,12,14,24})."
function _dvee_bf_cf_tableB(fia::Int, iforst::Int)::Float32
    fia == 400 && return 1.06f0
    fia == 602 && return 0.90f0
    fia == 621 && return 1.10f0
    fia == 694 && return 1.15f0
    fia == 731 && return 0.93f0
    fia == 742 && return 0.97f0
    fia == 824 && return 0.80f0
    fia == 830 && return 0.96f0
    fia == 832 && return 1.03f0
    (fia == 68 && iforst != 4) && return 0.80f0
    (fia == 110 && (iforst == 8 || iforst == 5)) && return 0.95f0
    (fia == 125 && iforst == 4) && return 0.96f0
    if fia == 129
        iforst in (11, 14, 12) && return 0.95f0
        iforst == 4 && return 0.96f0
    end
    (fia == 241 && iforst == 4) && return 0.80f0
    if fia == 802
        iforst == 8 && return 1.08f0
        iforst == 5 && return 0.96f0
    end
    fia == 806 && return iforst == 5 ? 1.03f0 : 1.10f0
    fia == 833 && return iforst == 8 ? 1.11f0 : 1.06f0
    fia == 835 && return iforst == 5 ? 0.94f0 : 0.98f0
    if fia == 837
        iforst == 8 && return 1.05f0
        iforst == 5 && return 0.96f0
        iforst == 4 && return 0.95f0
    end
    return 1f0
end

"R9VOL board feet (Scribner, VOL(2)) — TABLE B (IFORST∈{4,5,8,11,12,14,24}, incl. the CS 905 forests).
 Returns 0 for other IFORST (TABLE A/C not ported — no DVEE stand in the corpus uses them)."
function _dvee_boardft(fia::Int, dbh::Float32, h1::Int, iforst::Int)::Float32
    (h1 <= 0) && return 0f0
    (iforst in (4, 5, 8, 11, 12, 14, 24)) || return 0f0
    h1f = Float32(h1); d2 = dbh * dbh; d4 = d2 * d2
    r = (h1f * dbh - 3.75f0) / (24f0 * h1f - 10.5f0)
    vc = h1f * (1.0757f0 + 3.002f0 * r + 8.3776f0 * r * r)
    v2 = -0.092685f0 - 5.98f0 * vc - 2.9715f0 * dbh + 16.7022f0 * h1f +
         0.2471f0 * d2 * h1f - 0.91751f0 * h1f * h1f - 0.00876f0 * vc * vc +
         0.351046f0 * d2 + 0.00451f0 * d2 * h1f * h1f -
         0.00030183475f0 * d2 * h1f^3 + 0.0000019222413f0 * d4 * h1f * h1f
    return v2 * _dvee_bf_cf_tableB(fia, iforst)
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
    bfmind = 9f0                              # R9_MHTS CS board-min DBH (softwood & hardwood); PROD='01' at ≥bfmind
    bftopd = 7.6f0                            # saw top (SCFTOPD): CS uses 7.6 for both softwood & hardwood
                                              # (BH V4=SP_V4/1.06 ⇒ same H1 ⇒ same saw top; 9.6 hardwood default
                                              # would exceed a 9" stem's DBH and drop H1 ⇒ Scuft low)
    h2 = _r9_mhts_ht2prd(fia, dbh, httot, si, ba, region)
    # R9VOL merch cubic MCF = VOL(4)+VOL(7) (fvsvol.f:512). GCB = the pulp polynomial. For a PULP stem
    # (DBH<BFMIND ⇒ HT1PRD=0) VOL(4)=GCB and VOL(7)=0 ⇒ MCF=GCB. For a SAWTIMBER stem VOL(4)=the '912' sawlog
    # cubic and VOL(7)=PT·GCB (topwood), PT=(98.461−1.394P+0.004P²)·0.01, P=HT1PRD/HT2PRD·100 ⇒ MCF=VOL(4)+PT·GCB.
    # (The DBH≥DBHMIN species gate is applied by the caller, like the Clark path.) All pieces live-validated.
    gcb = 0f0
    if h2 > 0
        h2f = Float32(h2)
        gcb = 0.001f0 * dbh * dbh * (1.9f0 + 0.01f0 * dbh) *
              (0.208f0 * h2f - 0.009984f0 * h2f * h2f + 0.04f0 / h2f) * 79f0
    end
    mcf = gcb
    scf = 0f0                                  # sawtimber cubic (VOL(4)_saw), gated on D≥SCFMIND by the caller
    bf = 0f0                                   # board feet (VOL(2), Scribner), gated on D≥BFMIND by the caller
    if dbh >= bfmind && h2 > 0                 # sawtimber: MCF = VOL(4)_saw + VOL(7)=PT·GCB, SCF = VOL(4)_saw
        h1 = _r9_mhts_ht1prd(fia, dbh, httot, si, ba, region, bftopd, h2)
        bf = _dvee_boardft(fia, dbh, h1, iforst)
        if h1 > 0
            h2f = Float32(h2)
            h1f = Float32(h1); d2 = dbh * dbh; d3 = d2 * dbh; d4 = d2 * d2
            term1 = -1.70774f0 + 0.051321f0 * dbh + 0.58857f0 * h1f +
                    0.0193547f0 * d2 + 0.0237324f0 * h1f * d2
            term2 = -0.04821f0 * h1f^2 - 0.0002174f0 * d2 * h1f^2 -
                    0.0000239f0 * d2 * h1f^3 + 0.00000795f0 * d3 * h1f^2
            term3 = -0.00000057f0 * d3 * h1f^3 - 0.000000035f0 * d4 * h1f^2
            v4saw = Float32(term1 + term2 + term3) * _dvee_saw_cf(fia, dbh)
            p = h1f / h2f * 100f0
            pt = (98.461f0 - 1.394f0 * p + 0.004f0 * p * p) * 0.01f0
            mcf = v4saw + pt * gcb
            scf = v4saw
        end
    end
    return (tcf, mcf, scf, bf)
end
