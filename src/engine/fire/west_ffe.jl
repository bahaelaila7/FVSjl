# =============================================================================
# fire/west_ffe.jl — the western (non-eastern) FFE volume basis for snags and the live-carbon stem.
#
# {v}/fmsvol.f (identical in every western variant, :146-151): FMSVOL/FMSVL2 call NATCRS on the tree's DBH and
# height with XHT=−1 for an intact stem ⇒ LTKIL=.FALSE. ⇒ no CFTOPK, and return VOL2HT = MAX(X,TCF)
# (X = 0.005454154·H) — or MCF when LMERCH (the merch-carbon pool, fmcrbout.f:127 LVWEST). A snag that has lost
# height (FMSNGHT) is volumed as its death-form tree (DBHS, HTDEAD) trimmed by CFTOPK at IHT=INT(HTIH·100).
# Only CS/LS/NE/SN take the MAX(X,MCF)/SCF branch.
#
# `ffe_west_nocut` gives each western variant's NATCRS (TCF, MCF) on an arbitrary (d, h) with no top-kill — the same
# kernels as its compute_volumes (FW2 / MATW r4vol / DVE / Behre / NVB), minus the CFTOPK broken-top trim — plus the
# FMSVOL bark ratio BRATIO(JS,D,H) and whether that equation family takes the CFTOPK trim (`trim`: the DVE woodland
# families in CI/TT/UT/EM are kept untrimmed as in their compute_volumes; CR trims every family — fmsvol.f calls
# NATCRS, which returns CTKFLG=.TRUE. (fvsvol.f:535), and FMSVOL has no vols.f:191 NVB exemption).
# `nothing` for a variant not yet on this layer (its own snag/live paths stay in force).
# =============================================================================

"Western variants whose FFE snag bole / live-carbon stem use `ffe_west_nocut` (fmsvol.f non-eastern branch)."
_ffe_west_vol(v) = v isa InlandEmpire || v isa Kootenai || v isa CentralIdaho || v isa Teton || v isa Utah ||
                   v isa EasternMontana || v isa CentralRockies

"""
    ffe_west_nocut(s, sp, d, h) -> (tcf, mcf, bark, trim) | nothing

NATCRS total/merch cubic for species `sp` at DBH `d`, height `h`, with NO top-kill (fmsvol.f XHT=−1), for the
variants in `_ffe_west_vol`. MCF is gated at the variant's DBHMIN like its compute_volumes; `bark` is the variant
BRATIO that FMSVOL hands CFTOPK and `trim` whether that family takes the CFTOPK trim (`ffe_west_snag_vol_at`).
"""
function ffe_west_nocut(s::StandState, sp::Int, d::Float32, h::Float32)
    v = s.variant
    _ffe_west_vol(v) || return nothing
    (d < 1f0 || h <= 0f0 || sp < 1 || sp > length(s.species.vol_eq)) && return (0f0, 0f0, 1f0, false)
    eq = s.species.vol_eq[sp]; se = strip(eq); mdl = length(se) >= 7 ? se[4:6] : "   "
    if v isa EasternMontana
        tcf, mcf = em_nocut_cuft(s, sp, d, h)
        return (tcf, mcf, em_bratio(sp, d), startswith(eq, "I") || mdl == "FW2")
    elseif v isa CentralRockies                                # cr: NVEL DVE / NVB / FW2 (compute_volumes_cr!)
        tcf, mcf = cr_nocut_cuft(s, sp, d, h)
        return (tcf, mcf, cr_bratio(s.coef.species, sp, d, Int(s.plot.model_type)), true)
    elseif v isa Kootenai                                      # kt: FW2 only (compute_volumes_kt!)
        bark = bark_ratio(s.calib.bark_a, s.calib.bark_b, sp, d)
        startswith(eq, "I") || return (0f0, 0f0, bark, false)
        w = cr_fw2_vol(eq, d, h; bark = bark, topd = 4.5f0, bftopd = 4.5f0, stump = 1f0, iregn = 1)
        return (max(w[1], 0f0), d >= (sp == 7 ? 6f0 : 7f0) ? max(w[4] + w[7], 0f0) : 0f0, bark, true)
    elseif v isa InlandEmpire                                  # ie: region-6 Behre / FW2 / region-1-2 DVE
        bark = ie_bratio(sp, d); dbhmin = sp == 7 ? 6f0 : 7f0
        if occursin("BEH", eq)
            tcf, mcf, _, _ = ie_behre_vol(sp, Int(s.plot.forest_idx), d, h, bark)
            return (max(tcf, 0f0), d >= dbhmin ? max(mcf, 0f0) : 0f0, bark, true)
        end
        w = startswith(eq, "I") ?
            cr_fw2_vol(eq, d, h; bark = bark, topd = 4.5f0, bftopd = 4.5f0, stump = 1f0, iregn = 1,
                       sf_hs = true, ht2td = zeros(Float32, 2)) :
            ie_dve_vol(eq, d, h, true)
        return (max(w[1], 0f0), d >= dbhmin ? max(w[4] + w[7], 0f0) : 0f0, bark, true)
    end
    # CI / TT / UT (compute_volumes_{ci,tt,ut}!): MATW r4vol, FW2 (region 4), DVE woodland (not trimmed)
    dbhmin = sp == 7 ? 7f0 : 8f0
    bark = v isa CentralIdaho ? ci_bratio(s.coef.species, sp, d) :
           v isa Teton        ? tt_bratio(sp, d) :
                                bark_ratio(s.calib.bark_a, s.calib.bark_b, sp, d)
    if mdl == "MAT"
        mtopp = v isa Teton ? 6f0 : 6f0 * bark                # TT: MTOPS=TOPD (no bark), CI/UT: TOPD·BARK
        tcf, mcf = r4vol_volumes(eq, d, h, mtopp, 0f0)
        return (max(tcf, 0f0), d >= dbhmin ? max(mcf, 0f0) : 0f0, bark, true)
    elseif mdl == "FW2" && !(v isa Teton)
        w = cr_fw2_vol(eq, d, h; bark = bark, topd = 6f0, bftopd = 6f0, stump = 1f0, iregn = 4)
        return (max(w[1], 0f0), d >= dbhmin ? max(w[4] + w[7], 0f0) : 0f0, bark, true)
    end
    vol1 = (v isa Utah && se[1] == '3') ? r3d2hv_vol1(eq, d, h) : r4d2h_vol1(eq, d, h)
    return (max(vol1, 0f0), d >= (v isa Teton ? 8f0 : dbhmin) ? max(vol1, 0f0) : 0f0, bark, false)
end

"FMSVOL snag bole for an intact stem: MAX(0.005454154·H, TCF) (fmsvol.f:150). `nothing` off the western layer."
function ffe_west_snag_bole(s::StandState, sp::Int, d::Float32, h::Float32)
    w = ffe_west_nocut(s, sp, d, h); w === nothing && return nothing
    h <= 0f0 && return 0f0
    return max(0.005454154f0 * h, w[1])
end

"""
    ffe_west_snag_vol_at(s, sp, d, htd, htcur) -> Float32 | nothing

FMSVOL(II, XHT=HTIH): the snag's death-form tree (DBHS, HTDEAD) through NATCRS, then CFTOPK at IHT=INT(XHT·100)
when the snag has lost height (the fat lower bole, not a short tree), VOL2HT = MAX(0.005454154·HTDEAD, TCF).
"""
function ffe_west_snag_vol_at(s::StandState, sp::Int, d::Float32, htd::Float32, htcur::Float32)
    w = ffe_west_nocut(s, sp, d, htd); w === nothing && return nothing
    htd <= 0f0 && return 0f0
    tcf, mcf, bark, trim = w
    if htcur < htd && tcf > 0f0 && trim
        tcf, _ = cr_cftopk(tcf, mcf, d, htd, tcf, bark, unsafe_trunc(Int, htcur * 100f0), 1f0, 4.5f0)
    end
    return max(0.005454154f0 * htd, tcf)
end

# ---------------------------------------------------------------------------------------------------------------
# Snag height loss (FMSNGHT) parameters. {v}/fmvinit.f HTX(I,1:4) per species (hard first-50% / hard 50-95% / soft
# first-50% / soft 50-95%; TT sets 1:2 and copies to 3:4, fmvinit.f:411-412), extracted from the species SELECT CASE.
# HTR1 = 0.0228, HTR2 = 0.01 in all six; HTXSFT = 2 (IE/KT/CI/TT), 10 (UT, CR — cr/fmvinit.f:132). jl seeded HTX only
# for NE/LS/EM, so these variants' snags never lost height (live IE S248112 LP input snags shrink every year).
const _IE_FM_HTX = NTuple{4,Float32}[(0.9f0, 0.9f0, 0.9f0, 0.9f0), (0.9f0, 0.9f0, 0.9f0, 0.9f0), (0.9f0, 0.9f0, 0.9f0, 0.9f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0)]
const _KT_FM_HTX = NTuple{4,Float32}[(0.9f0, 0.9f0, 0.9f0, 0.9f0), (0.9f0, 0.9f0, 0.9f0, 0.9f0), (0.9f0, 0.9f0, 0.9f0, 0.9f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0)]
const _CI_FM_HTX = NTuple{4,Float32}[(0.4f0, 0.4f0, 0.4f0, 0.4f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.5f0, 1.5f0, 1.5f0, 1.5f0), (0.9f0, 0.9f0, 0.9f0, 0.9f0), (0.3f0, 0.3f0, 0.3f0, 0.3f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (1.5f0, 1.5f0, 1.5f0, 1.5f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0)]
const _TT_FM_HTX = NTuple{4,Float32}[(0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (1.01f0, 4.61f0, 1.01f0, 4.61f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0)]
const _UT_FM_HTX = NTuple{4,Float32}[(0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.9f0, 0.9f0, 0.9f0, 0.9f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0)]

const _CR_FM_HTX = NTuple{4,Float32}[(1.494f0, 1.494f0, 1.494f0, 1.494f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (0.9f0, 0.9f0, 0.9f0, 0.9f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.9f0, 0.9f0, 0.9f0, 0.9f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0)]

"{v}/fmvinit.f HTX table for the western layer (EM keeps its own EM_FM_HTX seeding), else `nothing`."
_ffe_west_htx(v) = v isa InlandEmpire ? _IE_FM_HTX : v isa Kootenai ? _KT_FM_HTX : v isa CentralIdaho ? _CI_FM_HTX :
                   v isa Teton ? _TT_FM_HTX : v isa Utah ? _UT_FM_HTX : v isa CentralRockies ? _CR_FM_HTX : nothing
_snag_htr1(::InlandEmpire) = 0.0228f0
_snag_htr1(::Kootenai) = 0.0228f0
_snag_htr1(::CentralIdaho) = 0.0228f0
_snag_htr1(::Teton) = 0.0228f0
_snag_htr1(::Utah) = 0.0228f0
_snag_htr1(::CentralRockies) = 0.0228f0
"HTXSFT (soft-snag height-loss multiplier, {v}/fmvinit.f): UT/CR 10, the default 2 elsewhere."
_snag_htxsft(v) = (v isa Utah || v isa CentralRockies) ? 10f0 : 2f0

"""
Snag fall FMSFALL ({v}/fmsfall.f) for the western layer: the linear small-snag fall runs below 18" (IE/KT/CI/TT/UT;
EM 12"), with no redcedar exception, and some species fall slower once ≥ 18": BASE = MAX(0.01, BASE·0.32) for
CI ksp 2/8 (ci/fmsfall.f:39), TT 3/8, UT 3/5/8, CR (CASE DEFAULT) 3/17-19 (the shared tt/ut/cr fmsfall.f LDFSP).
Returns (linear_max, slow).
"""
function _ffe_west_fall(v, ksp::Integer)
    v isa EasternMontana && return (12f0, false)
    slow = v isa CentralIdaho ? (ksp == 2 || ksp == 8) :
           v isa Teton        ? (ksp == 3 || ksp == 8) :
           v isa Utah         ? (ksp == 3 || ksp == 5 || ksp == 8) :
           v isa CentralRockies ? (ksp == 3 || 17 <= ksp <= 19) : false
    return (18f0, slow)
end
