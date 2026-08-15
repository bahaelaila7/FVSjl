# =============================================================================
# fire/fuel_model.jl — dynamic surface fuel model construction (FFE chunk F5b)
#
# Ported from: bin/FVSsn_buildDir/fmcfmd2.f (FMCFMD3, the SN custom model) +
# fmgfmv.f (the dead-herb moisture split) + the FMINIT USAV/UBD/CANMHT defaults.
#
# The SN FFE does not use a static fuel model — it builds a *custom* one each fire from
# the stand's own fuels: the down-wood pools (FireState.cwd, F3), the live herb/shrub
# load (FireState.flive, F3), and the understory crown biomass (crown_biomass, F2). The
# result (loads, SAV, depth, moisture of extinction) is exactly what the Rothermel model
# (F5) consumes — so this is the keystone that ties F2+F3 into the fire behavior and
# makes the crown-biomass chunk a live input rather than an inert one.
# =============================================================================

const _FM_USAV   = (2000f0, 1800f0, 1500f0)   # USAV: dead-1hr / live-herb / live-woody SAV (fminit.f:826)
const _FM_UBD    = (0.10f0, 0.75f0)           # UBD: fuelbed bulk-density bounds (fminit.f:829)
const _FM_CANMHT = 6.0f0                       # CANMHT: understory height threshold, ft (fminit.f:147)
const _TONS_TO_LBFT2 = 0.04591f0              # tons/acre → lb/ft²

@inline _fm_algslp2(x, x1, x2, y1, y2) =       # 2-point clamped linear interpolation (ALGSLP)
    x < x1 ? y1 : x >= x2 ? y2 : y1 + (y2 - y1) / (x2 - x1) * (x - x1)

"""
    standard_fuel_model(coef, model) -> (load, sav, depth, mext)

Rothermel inputs for one standard fire-behavior fuel model (1–13, the Anderson models,
fminit.f / `fire_fuel_models.csv`). `load[2,4]`/`sav[2,4]` are loads (lb/ft²) and surface-
area-to-volume by [1=dead/2=live, class]; the 10-hr / 100-hr / dead-herb / live-herb SAVs
take the FFE constant defaults (109 / 30 / 1500 / 1500). `depth` is bed depth (ft), `mext`
the dead moisture of extinction.
"""
# The fuel model for `model`, honoring a DEFULMOD custom override (FireState.defulmod) if one was defined,
# else the standard table (standard_fuel_model). Returns (load[2,4], sav[2,4], depth, mext).
@inline function fuel_model_resolved(s::StandState, model::Integer)
    fs = s.fire
    if fs !== nothing && !isempty(fs.defulmod)
        ov = get(fs.defulmod, Int32(model), nothing)
        ov !== nothing && return ov
    end
    return standard_fuel_model(s.coef, model)
end

function standard_fuel_model(coef::SpeciesCoefficients, model::Integer)
    m = @view coef.ffe_fuel_models[model, :]   # [sav_1hr, sav_lwoody, l_1hr,l_10,l_100,l_lwoody,l_lherb, depth, mext]
    load = zeros(Float32, 2, 4); sav = zeros(Float32, 2, 4)
    load[1, 1] = m[3]; load[1, 2] = m[4]; load[1, 3] = m[5]      # dead 1/10/100-hr
    load[2, 1] = m[6]; load[2, 2] = m[7]                          # live woody / herb
    sav[1, 1] = m[1]; sav[1, 2] = 109f0; sav[1, 3] = 30f0; sav[1, 4] = 1500f0
    # live-herb SAV: the FVS SURFVL(I,2,2) default is 1500, overridden per-model for the extended
    # Scott-Burgan grass/timber-understory models (fminit.f). Read from the trailing sav_lherb column
    # when present (LS); the standard 13-model CSVs omit it ⇒ the 1500 default (bit-identical).
    sav[2, 1] = m[2]; sav[2, 2] = length(m) >= 10 ? m[10] : 1500f0
    return (load, sav, m[8], m[9])
end

# ----------------------------------------------------------------------------
# FMCFMD/FMDYN — weighted *standard* fuel-model selection (FFE chunk F4-select).
#
# The SN FFE's real surface-fire path is NOT the custom dynamic model above: FMBURN
# calls FMCFMD, which (a) picks a set of candidate standard fuel models from the FFE
# forest type + fuel moisture, then (b) hands them to FMDYN, which places the stand's
# (SMALL, LARGE) down-wood point in the 2-D fuel space and weights the candidate models
# by inverse perpendicular distance to each model's iso-line. FMFINT then runs Rothermel
# on each weighted model and sums the weighted flame (F4-weight). build_dynamic_fuel_model
# is only used when FUELMODL/IFLOGIC forces a static/custom model — not for snt01 stand 4.
# ----------------------------------------------------------------------------

# XPTS (fmcfmd.f:79): per-model iso-line as (SMALL-intercept, LARGE-intercept), tons/ac. PER-VARIANT (the
# model-10 intercept + candidate-set differ): SN has 14 models w/ model-10 (10,30); NE has 13 w/ model-10 (15,30).
const _FMD_XPTS = Float32[
    5  15;  5  15;  5  15;  5  15;  5  15;  5  15;  5  15;  # models 1–7
    5  15;  5  15;                                          # models 8–9
   10  30; 15  30; 30  60; 45 100; 30  60]                  # models 10–14
const _FMD_XPTS_NE = Float32[                               # ne/fmcfmd.f:78-91 (model 10 = (15,30); no model 14)
    5  15;  5  15;  5  15;  5  15;  5  15;  5  15;  5  15;  # models 1–7
    5  15;  5  15;                                          # models 8–9
   15  30; 15  30; 30  60; 45 100;  0   0]                  # models 10–13 (14 = degenerate/unused)
const _FMD_XPTS_CR = Float32[                               # cr/fmcfmd.f:119-131 — model 10 = (15,30); ICLSS=12
    5  15;  5  15;  5  15;  5  15;  5  15;  5  15;  5  15;  # models 1–7
    5  15;  5  15;                                          # models 8–9
   15  30; 15  30; 30  60;                                  # models 10–12 (10 shares iso-line with 11)
    0   0;  0   0]                                          # 13,14 unused in CR (degenerate ⇒ _fmdyn skips)
const _FMD_XPTS_IE = Float32[                               # ie/fmcfmd.f:22-36 — ICLSS=14; model 10 = (15,30),
    5  15;  5  15;  5  15;  5  15;  5  15;  5  15;  5  15;  # models 1–7    # 14 = (30,60) (shares w/ 12)
    5  15;  5  15;                                          # models 8–9
   15  30; 15  30; 30  60; 45 100; 30  60]                  # models 10–14 (IE has all 14; = kt/em/bm/ci)
fmd_xpts(::Northeast) = _FMD_XPTS_NE
fmd_xpts(::CentralRockies) = _FMD_XPTS_CR
# IE-family western XPTS (ie/fmcfmd.f). KT identical; EM/BM/CI share the XPTS breakpoints (candidate-set
# selection via NIFMHAB/IDRY differs per variant — ported with each variant's fmcfmd).
fmd_xpts(::InlandEmpire) = _FMD_XPTS_IE
fmd_xpts(::Kootenai) = _FMD_XPTS_IE
fmd_xpts(::EasternMontana) = _FMD_XPTS_IE   # em/fmcfmd.f XPTS verified identical to ie
fmd_xpts(::CentralIdaho) = _FMD_XPTS_IE     # ci/fmcfmd.f XPTS verified identical to ie
fmd_xpts(::BlueMountains) = _FMD_XPTS_IE     # bm/fmcfmd.f ICLSS=14, XPTS identical to ie
fmd_xpts(::Teton) = _FMD_XPTS_CR            # TT/UT ICLSS=12 (models 1-12), same XPTS breakpoints as CR
fmd_xpts(::Utah) = _FMD_XPTS_CR
fmd_xpts(::AbstractVariant) = _FMD_XPTS
const _FMD_ICLSS = 14

# LS (ls/fmcfmd.f) uses ICLSS=22 fuel-model classes: the standard 13 plus the extended
# Minnesota models 105/142/143/146/161/162/164/186/189 (classes 14–22), mapped to their
# printed model numbers by IPTR. XPTS (fmcfmd.f:104) iso-line intercepts (SMALL, LARGE):
# model 10/11 = (15,30), 12 = (30,60), 13 = (45,100), everything else (5,15).
const _FMD_IPTR_LS = Int[1,2,3,4,5,6,7,8,9,10,11,12,13,105,142,143,146,161,162,164,186,189]
const _FMD_XPTS_LS = Float32[
    5  15;  5  15;  5  15;  5  15;  5  15;  5  15;  5  15;  # models 1–7
    5  15;  5  15;                                          # models 8–9
   15  30; 15  30; 30  60; 45 100;                          # models 10–13
    5  15;  5  15;  5  15;  5  15;  5  15;  5  15;  5  15;  # 105,142,143,146,161,162,164
    5  15;  5  15]                                          # 186,189
const _FMD_ICLSS_LS = 22
const _FMD_MXFMOD = 5     # MXFMOD (FMPARM.F77)

"SMALL/LARGE down-wood loads (tons/ac): classes 1–3 + litter(10) are SMALL, 4–9 are LARGE (fmtret.f:382)."
function _small_large_fuel(fs)
    small = 0f0; large = 0f0
    @inbounds for k in 1:2, l in 1:4
        small += fs.cwd[1, k, l] + fs.cwd[2, k, l] + fs.cwd[3, k, l] + fs.cwd[10, k, l]
        for j in 4:9
            large += fs.cwd[j, k, l]
        end
    end
    return small, large
end

"""
    select_fuel_models(s, mois) -> Vector{Tuple{Int,Float32}}

The SN weighted standard fuel models (FMCFMD + FMDYN) for the current stand: candidate
models chosen from the FFE forest type (`ffe_forest_type`/FMSNFT) and dead 100-hr fuel
moisture `mois[1,4]`, then weighted by the (SMALL, LARGE) down-wood point's inverse
distance to each model's iso-line. Returns up to `MXFMOD` (model, weight) pairs whose
weights sum to 1. This is the input FMFINT integrates over for the surface fire.
"""
function select_fuel_models(s::StandState, mois::AbstractMatrix{Float32}; fire_basis::Bool = false)
    # FUELMODL (fmusrfm.f → fmcfmd.f:113 IF(LUSRFM)RETURN): a forced fuel-model list for this cycle skips
    # the FMCFMD auto-selection entirely.
    if s.fire !== nothing && !isempty(s.fire.fuelmodl)
        yr = Int(current_cycle_year(s)); fvscyc = Int(s.control.cycle) + 1
        for (date, pairs) in s.fire.fuelmodl
            (Int(date) == yr || (0 < Int(date) < 1000 && Int(date) == fvscyc)) || continue
            return [(Int(m), w) for (m, w) in pairs]
        end
    end
    eqwt = zeros(Float32, _FMD_ICLSS)
    iffeft = ffe_forest_type(s)
    # An actual SIMFIRE burns on the start-of-cycle + 1-annual-step down wood (FVS interleaves the annual
    # fuel loop with FMBURN); the PotFire report and the no-stash case use the live (period-end) cwd.
    sm, lg = (fire_basis && s.fire.fire_smlg[1] >= 0f0) ? s.fire.fire_smlg : _small_large_fuel(s.fire)
    m14 = mois[1, 4]                                   # dead 100-hr (3+") moisture

    # NE FMCFMD (ne/fmcfmd.f:120-148) is a stripped-down selection (NO forest-type/moisture logic): base
    # model 9 + the natural-fuel candidates 10/12/13, resolved by the (SMALL,LARGE) point. (Model 11 = the
    # 5-yr post-activity fuel-jump candidate, AFWT/SLCHNG — deferred; net01 takes the natural-fuel path.)
    if s.variant isa Northeast
        eqwt[9] = 1f0; eqwt[10] = 1f0; eqwt[12] = 1f0; eqwt[13] = 1f0
        return _fmdyn(sm, lg, eqwt, fmd_xpts(s.variant))
    end

    # LS uses a full cover-type × PERCOV × season selection (ls/fmcfmd.f), NOT the SN forest-type path.
    if s.variant isa LakeStates
        return ls_select_fuel_models(s, mois, sm, lg)
    end

    # CR / TT / UT share cr/fmcfmd.f's cover-type algorithm (VARACD-branched, parameterized in cr_select).
    if s.variant isa CentralRockies || s.variant isa Teton || s.variant isa Utah
        return cr_select_fuel_models(s, mois, sm, lg)
    end

    # IE-family (ie/fmcfmd.f) — the simplest western selection: MAPDRY habitat→dryness → base model +
    # natural fuels {10,12,13}. KT shares IE's fmcfmd + MAPDRY (verified identical).
    if s.variant isa InlandEmpire || s.variant isa Kootenai
        return ie_select_fuel_models(s, mois, sm, lg)
    end

    # EM (em/fmcfmd.f) — same structure as IE but the two candidate models come from the habitat table
    # EMMD: (M1,M2)=MD1/MD2[IEMTYP], weighted by PERCOV, + natural fuels {10,12,13}.
    if s.variant isa EasternMontana
        return em_select_fuel_models(s, mois, sm, lg)
    end

    # CI (ci/fmcfmd.f) — cover-type-group (ICT=MAPPVG[ICINDX]) with a grand-fir-understory sub-model.
    if s.variant isa CentralIdaho
        return ci_select_fuel_models(s, mois, sm, lg)
    end

    # BM (bm/fmcfmd.f) — DF/PP-fraction weighted between a size-class-distribution path (BMSTAGE canopy
    # stratification + BMSZCLS Monte-Carlo size classes) and a PERCOV-band path.
    if s.variant isa BlueMountains
        return bm_select_fuel_models(s, mois, sm, lg)
    end

    # NC (nc/fmcfmd.f + cwhr.f) — California CWHR size×density structural-stage classification.
    if s.variant isa Klamath
        return nc_select_fuel_models(s, mois, sm, lg)
    end

    # WS (ws/fmcfmd.f + cwhr.f) — California CWHR, WS's own CWHRFMD 12×18 + 43-species forest-type + models 25/26.
    if s.variant isa WestSierra
        return ws_select_fuel_models(s, mois, sm, lg)
    end

    # CA (ca/fmcfmd.f + cwhr.f) — California CWHR, CA CWHRFMD 11×18 + 50-species forest-type; reuses _ca_cwhr.
    if s.variant isa CentralCalifornia
        return ca_select_fuel_models(s, mois, sm, lg)
    end

    # --- SN candidate-model selection (fmcfmd.f:131) ---
    if iffeft in (1, 2, 3)                             # hardwood / hwd-pine / pine-hwd
        if sm > 6f0
            eqwt[5] = 1f0
        else
            moiswt8 = 0f0; moiswt9 = 0f0
            if m14 <= 0.15f0
                moiswt9 = 1f0
            elseif m14 > 0.25f0
                moiswt8 = 1f0
            else
                moiswt9 = 1f0 - (m14 - 0.15f0) / 0.1f0
                moiswt8 = 1f0 - (0.25f0 - m14) / 0.1f0
            end
            if sm <= 4f0
                eqwt[8] = moiswt8; eqwt[9] = moiswt9
            else
                eqwt[8] = (1f0 - (sm - 4f0) / 2f0) * moiswt8
                eqwt[9] = (1f0 - (sm - 4f0) / 2f0) * moiswt9
                eqwt[5] = 1f0 - (6f0 - sm) / 2f0
            end
        end
    elseif iffeft in (4, 8)                            # pine / saint francis
        if m14 <= 0.15f0
            eqwt[9] = 1f0
        elseif m14 > 0.25f0
            eqwt[8] = 1f0
        else
            eqwt[9] = 1f0 - (m14 - 0.15f0) / 0.1f0
            eqwt[8] = 1f0 - (0.25f0 - m14) / 0.1f0
        end
    elseif iffeft in (5, 6)                            # pine bluestem / oak savannah
        eqwt[2] = 1f0
    elseif iffeft == 7                                 # eastern redcedar
        rcht = 0f0; rctpa = 0f0
        t = s.trees
        @inbounds for i in 1:t.n
            if Int(t.species[i]) == 2                  # redcedar
                rcht += t.height[i]; rctpa += t.tpa[i]
            end
        end
        rcht = rctpa > 0f0 ? rcht / rctpa : 0f0
        if rcht > 7.5f0
            eqwt[4] = 1f0
        elseif rcht <= 4.5f0
            eqwt[6] = 1f0
        else
            eqwt[6] = 1f0 - (rcht - 4.5f0) / 3f0
            eqwt[4] = 1f0 - (7.5f0 - rcht) / 3f0
        end
    elseif iffeft == 9
        eqwt[6] = 1f0
    end
    # models 10 & 12 are always candidates for natural fuels (fmcfmd.f:202)
    eqwt[10] = 1f0; eqwt[12] = 1f0

    return _fmdyn(sm, lg, eqwt, fmd_xpts(s.variant))
end

# FINDMOD (fmcfmd.f:641): the class index whose IPTR entry is model `jmod`; 8 if not found.
@inline _ls_findmod(jmod::Integer) = (k = findfirst(==(jmod), _FMD_IPTR_LS); k === nothing ? 8 : k)

"""
    ls_select_fuel_models(s, mois, sm, lg) -> Vector{Tuple{Int,Float32}}

The Lake States weighted standard fuel models (ls/fmcfmd.f + FMDYN). LS selects candidate
models from the stand's COVER-TYPE metagroup (jack pine / northern hardwood / red-white
pine / mixed wood / oak / aspen-birch, by basal area), then a per-cover-type branch keyed
on percent cover, arithmetic mean height, burn season, understory composition, and midflame
wind; a post-selection block adds the natural-fuel candidates 10/12/13 (cover-type-limited).
`(sm, lg)` is the (SMALL, LARGE) down-wood point FMDYN resolves against.

Activity-fuel (the 5-yr post-harvest model-11 fuel jump, SLCHNG/LATFUEL/AFWT) is DEFERRED
— natural stands take LATFUEL=false ⇒ AFWT=0 (same deferral as the SN/NE model-11 path).
"""
function ls_select_fuel_models(s::StandState, mois::AbstractMatrix{Float32}, sm::Float32, lg::Float32)
    t = s.trees; coef = s.coef; fs = s.fire
    eqwt = zeros(Float32, _FMD_ICLSS_LS)
    percov  = fs.percov
    burnseas = Int(fs.burnseas)
    fwind   = fs.swind * fire_wind_reduction(percov)      # 20-ft → midflame wind
    itype   = Int(s.plot.forest_idx)                      # ITYPE (/PLOT/) — FDc/FDn habitat checks
    ldry    = false                                       # DROUGHT (IDRYB..IDRYE) — none in these stands
    # ATAVH (arithmetic mean tree height) drives the jack-pine height thresholds
    hsum = 0f0; hden = 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        hsum += t.height[i] * t.tpa[i]; hden += t.tpa[i]
    end
    fmavh = hden > 0f0 ? hsum / hden : 0f0

    # crown area (ft²/ac) helper for a record: forest-grown crown width (CWCALC iwho=0)
    carea(i) = begin
        sp2 = s.species.code2[Int(t.species[i])]
        cw = crown_width(coef, sp2, t.dbh[i], t.height[i], Float32(t.crown_pct[i]), 0,
                         s.plot.latitude, s.plot.longitude, s.plot.elevation)
        3.1415927f0 * cw * cw / 4f0 * t.tpa[i]
    end

    # understory flags (fmcfmd.f:167-190): balsam-fir (sp8) / small-tree / conifer (≤14) /
    # balsam-fir-or-white-pine 1–3″ counts, each ≥500 TPA turns its flag on.
    bfc = 0f0; smc = 0f0; cnc = 0f0; bfwp = 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        sp = Int(t.species[i]); d = t.dbh[i]; p = t.tpa[i]
        (sp == 8 && 1f0 <= d <= 3f0) && (bfc += p)
        (d <= 2f0) && (smc += p)
        (sp <= 14 && 1f0 <= d <= 3f0) && (cnc += p)
        ((sp == 8 || sp == 5) && 1f0 <= d <= 3f0) && (bfwp += p)
    end
    lbfunder = bfc >= 500f0; lsmtrees = smc >= 500f0
    lcnunder = cnc >= 500f0; lbfwpund = bfwp >= 500f0

    # LT3: dead surface fuel 0–3″ (down-wood size classes 1–3, tons/ac)
    lt3 = 0f0
    @inbounds for isz in 1:3, k in 1:2, l in 1:4
        lt3 += fs.cwd[isz, k, l]
    end

    # ---- cover-type metagroup basal areas (fmcfmd.f:215-253) ----
    JPCT=1; NHCT=2; RPCT=3; MWCT=4; OACT=5; ABCT=7
    ctba = zeros(Float32, 7); stndba = 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        sp = Int(t.species[i]); d = t.dbh[i]
        x = t.tpa[i] * d * d * 0.0054542f0
        if sp == 1
            ctba[JPCT] += x; stndba += x
        elseif sp in (12,18,19,25,26,27,28)
            ctba[NHCT] += x; stndba += x
        elseif sp in (3,4,5)
            ctba[RPCT] += x; stndba += x
        elseif sp in (24,40,41,43)
            ctba[MWCT] += x; ctba[ABCT] += x; stndba += x
        elseif sp in (6,7,8,9)
            (d >= 5f0 ? ctba[MWCT] += x : ctba[ABCT] += x); stndba += x
        elseif sp in (30,31,32,33,35,36)
            ctba[OACT] += x; stndba += x
        elseif sp == 34
            ctba[NHCT] += x; ctba[OACT] += x; stndba += x
        end
    end

    # dominant metagroup ICT; combined oak+pine can override to oak or pine (fmcfmd.f:261-278)
    ict = 0
    if t.n > 0 && stndba > 0.001f0
        bamost = 0f0
        for i in 1:7
            ctba[i] > bamost && (bamost = ctba[i]; ict = i)
        end
        poct = ctba[OACT] + ctba[RPCT]
        if poct > bamost
            ict = ctba[RPCT] > ctba[OACT] ? RPCT : OACT
        end
    else
        ict = fs.covtyp_ict > 0 ? Int(fs.covtyp_ict) : RPCT   # OLDICT, FMVINIT default RPCT
    end
    fs.covtyp_ict = Int32(ict)

    # ---- detailed low-fuel-model selection per cover type (fmcfmd.f:280-537) ----
    if ict == JPCT
        if percov <= 70f0
            if fmavh <= 25f0
                eqwt[4] = 1f0
            elseif lbfunder
                burnseas <= 2 ? (eqwt[10] = 1f0) : (eqwt[_ls_findmod(162)] = 1f0)
            elseif itype == 6 || itype == 7
                eqwt[2] = 1f0
            else
                burnseas <= 2 ? (eqwt[10] = 1f0) : (eqwt[_ls_findmod(161)] = 1f0)
            end
        else
            if fmavh <= 15f0
                eqwt[4] = 1f0
            elseif fwind <= 4f0
                eqwt[8] = 1f0
            else
                eqwt[10] = 1f0
            end
        end
    elseif ict == NHCT
        mapbasba = 0f0; oakba = 0f0; asbiba = 0f0; hembary = 0f0; hemcra = 0f0; sba = 0f0
        @inbounds for i in 1:t.n
            t.tpa[i] > 0f0 || continue
            sp = Int(t.species[i]); d = t.dbh[i]; x = t.tpa[i]*d*d*0.0054542f0
            if sp in (18,19,25,26,27); mapbasba += x; sba += x
            elseif sp in (30,31,32,33,34,35,36); oakba += x; sba += x
            elseif sp in (24,40,41,43); asbiba += x; sba += x
            elseif sp == 12; hembary += x; sba += x; hemcra += carea(i)
            elseif sp == 28; sba += x
            end
        end
        hemcov = (1f0 - exp(-hemcra/43560f0)) * 100f0
        if hemcov >= 30f0 || asbiba > 0f0
            eqwt[8] = 1f0
        elseif sba > 0f0
            if mapbasba/sba >= 0.5f0
                eqwt[_ls_findmod(186)] = 1f0
            else
                burnseas == 4 ? (eqwt[9] = 1f0) : (eqwt[8] = 1f0)
            end
        else
            eqwt[8] = 1f0
        end
    elseif ict == RPCT
        sba = 0f0; rpba = 0f0; wpba = 0f0; pinecra = 0f0; hardba = 0f0
        @inbounds for i in 1:t.n
            t.tpa[i] > 0f0 || continue
            sp = Int(t.species[i]); d = t.dbh[i]; x = t.tpa[i]*d*d*0.0054542f0
            sba += x
            (3 <= sp <= 5) && (pinecra += carea(i))
            if sp in (3,4); rpba += x
            elseif sp == 5; wpba += x
            elseif sp >= 15; hardba += x
            end
        end
        pinecov = (1f0 - exp(-pinecra/43560f0)) * 100f0
        if lt3 >= 5f0
            eqwt[10] = 1f0
        elseif percov <= 0f0
            eqwt[8] = 1f0
        elseif pinecov/percov < 0.5f0 && hardba > 0f0
            eqwt[8] = 1f0
        elseif itype == 4 || itype == 10
            ldry ? (eqwt[_ls_findmod(146)] = 1f0) : (eqwt[_ls_findmod(143)] = 1f0)
        elseif lbfunder || lbfwpund
            fwind <= 4f0 ? (eqwt[10] = 1f0) : (eqwt[_ls_findmod(146)] = 1f0)
        elseif percov >= 50f0
            eqwt[9] = 1f0
        elseif percov <= 30f0 && rpba > wpba
            eqwt[2] = 1f0
        else
            eqwt[9] = 1f0
        end
    elseif ict == MWCT
        sba = 0f0; birba = 0f0; concra = 0f0; overcra = 0f0
        @inbounds for i in 1:t.n
            t.tpa[i] > 0f0 || continue
            sp = Int(t.species[i]); d = t.dbh[i]; x = t.tpa[i]*d*d*0.0054542f0
            (d >= 5f0) && (overcra += carea(i))
            if sp in (24,43); birba += x; sba += x
            elseif 1 <= sp <= 14; sba += x; (d >= 5f0) && (concra += carea(i))
            elseif sp in (40,41); sba += x
            end
        end
        overcov = (1f0 - exp(-overcra/43560f0)) * 100f0
        concov  = (1f0 - exp(-concra/43560f0)) * 100f0
        if sba <= 0f0
            eqwt[8] = 1f0
        elseif birba/sba >= 0.5f0
            eqwt[9] = 1f0
        elseif overcov > 0f0
            concov/overcov >= 0.30f0 ? (eqwt[10] = 1f0) : (eqwt[8] = 1f0)
        else
            eqwt[8] = 1f0
        end
    elseif ict == OACT
        sba = 0f0; oakba = 0f0; overcra = 0f0
        @inbounds for i in 1:t.n
            t.tpa[i] > 0f0 || continue
            sp = Int(t.species[i]); d = t.dbh[i]; x = t.tpa[i]*d*d*0.0054542f0
            sba += x
            (d >= 5f0) && (overcra += carea(i))
            (sp == 30 || sp == 35) && (oakba += x)
        end
        overcov = (1f0 - exp(-overcra/43560f0)) * 100f0
        if overcov >= 45f0
            if burnseas == 1
                mois[1,1] < 0.08f0 ? (eqwt[_ls_findmod(186)] = 1f0) : (eqwt[8] = 1f0)
            elseif sba <= 0f0
                eqwt[9] = 1f0
            elseif oakba/sba >= 0.30f0
                eqwt[_ls_findmod(189)] = 1f0
            else
                eqwt[9] = 1f0
            end
        elseif overcov >= 15f0
            lsmtrees ? (eqwt[_ls_findmod(142)] = 1f0) : (eqwt[2] = 1f0)
        else
            lsmtrees ? (eqwt[_ls_findmod(142)] = 1f0) : (eqwt[_ls_findmod(105)] = 1f0)
        end
    elseif ict == ABCT
        aspba = 0f0; birba = 0f0
        @inbounds for i in 1:t.n
            t.tpa[i] > 0f0 || continue
            sp = Int(t.species[i]); d = t.dbh[i]; x = t.tpa[i]*d*d*0.0054542f0
            (sp == 40 || sp == 41) && (aspba += x)
            (sp == 24 || sp == 43) && (birba += x)
        end
        if lt3 >= 5f0
            eqwt[10] = 1f0
        elseif !lcnunder
            birba > aspba ? (eqwt[9] = 1f0) : (eqwt[8] = 1f0)
        elseif fwind > 4f0
            eqwt[_ls_findmod(164)] = 1f0
        else
            birba > aspba ? (eqwt[9] = 1f0) : (eqwt[8] = 1f0)
        end
    end

    # ---- post-selection natural-fuel candidates (fmcfmd.f:541-620) ----
    # Activity fuel (model 11) deferred: LATFUEL=false ⇒ AFWT=0 for all natural stands.
    if ict in (JPCT, MWCT, OACT)
        eqwt[10] = 1f0; eqwt[12] = 1f0; eqwt[13] = 1f0
        ict == OACT && (eqwt[12] = 0f0; eqwt[13] = 0f0)   # oak never gets 12 or 13
    elseif ict == NHCT
        eqwt[10] = 0f0; eqwt[11] = 0f0; eqwt[12] = 0f0; eqwt[13] = 0f0
    elseif ict == RPCT
        eqwt[10] = 1f0; eqwt[11] = 0f0; eqwt[12] = 1f0; eqwt[13] = 1f0
    elseif ict == ABCT
        eqwt[10] = 1f0; eqwt[11] = 0f0; eqwt[12] = 0f0; eqwt[13] = 0f0
    end

    return _fmdyn(sm, lg, eqwt, _FMD_XPTS_LS; iptr = _FMD_IPTR_LS)
end

# CR species → fuel-model cover-type metagroup (cr/fmcfmd.f:281-315). 1=OBCT oak-brush,
# 2=PJCT pinyon-juniper, 3=PPCT ponderosa, 4=WSCT white-spruce, 5=SFCT spruce-fir,
# 6=LPCT lodgepole, 7=MCCT mixed-conifer, 8=ASCT aspen.
@inline function _cr_fm_covtype(sp::Int)::Int
    (23 <= sp <= 27)                        && return 1   # OBCT oak brush
    (sp == 12 || sp == 16 || 29 <= sp <= 35) && return 2  # PJCT pinyon-juniper
    (sp == 13 || sp == 36)                  && return 3   # PPCT ponderosa, chihuahua pine
    (sp == 19)                              && return 4   # WSCT white spruce
    (sp == 1 || sp == 17 || sp == 18)       && return 5   # SFCT spruce-fir
    (sp == 11)                              && return 6   # LPCT lodgepole
    (sp in (2, 3, 4, 5, 6, 7, 8, 9, 10, 14, 15, 37)) && return 7  # MCCT mixed conifer
    (sp in (20, 21, 22, 28, 38))            && return 8   # ASCT aspen
    return 7                                              # default mixed conifer
end

# TT/UT share cr's fmcfmd algorithm (VARACD-branched); only the species→covtype map + a few species-role
# indices differ (tt/ut fmcfmd.f SELECT CASE(ISP)). Groups: OBCT1 PJCT2 PPCT3 WSCT4 SFCT5 LPCT6 MCCT7 ASCT8.
@inline function _tt_fm_covtype(sp::Int)::Int   # TT (18 sp): no OBCT/WSCT
    (sp == 4 || sp == 11 || sp == 12)       && return 2   # PJCT pinyon-juniper (PM,UJ,RM)
    (sp == 10)                              && return 3   # PPCT ponderosa
    (sp == 5 || sp == 8 || sp == 9)         && return 5   # SFCT spruce-fir (BS,ES,AF)
    (sp == 7)                               && return 6   # LPCT lodgepole
    (sp in (1, 2, 3, 17))                   && return 7   # MCCT mixed conifer (WB,LM,DF,OS)
    (sp in (6, 13, 14, 15, 16, 18))         && return 8   # ASCT aspen/hardwoods
    return 7
end
@inline function _ut_fm_covtype(sp::Int)::Int   # UT (24 sp): +OBCT oak
    (sp == 13)                              && return 1   # OBCT oak brush
    (sp in (11, 12, 14, 15, 16))            && return 2   # PJCT pinyon-juniper
    (sp == 10)                              && return 3   # PPCT ponderosa
    (sp == 5 || sp == 8 || sp == 9)         && return 5   # SFCT spruce-fir
    (sp == 7)                               && return 6   # LPCT lodgepole
    (sp in (1, 2, 3, 4, 17, 23))            && return 7   # MCCT mixed conifer
    (sp in (6, 18, 19, 20, 21, 22, 24))     && return 8   # ASCT aspen/hardwoods
    return 7
end
# Western fmcfmd per-variant role params (cr_select_fuel_models generalization).
@inline _fm_covtype(v::AbstractVariant, sp::Int) =
    v isa Teton ? _tt_fm_covtype(sp) : v isa Utah ? _ut_fm_covtype(sp) : _cr_fm_covtype(sp)
@inline _fm_ppct_sp(v::AbstractVariant) = (v isa Teton || v isa Utah) ? 10 : 13   # ponderosa index (lppdom)
# ASCT lcundr crown-cover EXCLUSION set (tt/ut/cr fmcfmd ASCT SELECT CASE(ISP)): species skipped in the CVR10 sum.
@inline _fm_asct_excl(v::AbstractVariant, sp::Int) =
    v isa Teton ? (sp in (4, 6, 11, 12, 13, 14, 15, 16, 18)) :
    v isa Utah  ? (sp == 6 || (11 <= sp <= 16) || (18 <= sp <= 22) || sp == 24) :
                  (sp == 12 || sp == 16 || (20 <= sp <= 35) || sp == 38)

"""
    cr_select_fuel_models(s, mois, sm, lg) -> [(model, weight)]

CR FFE fuel-model selection (cr/fmcfmd.f `CASE('CR')`). CR is COVER-TYPE based (like the
western TT/UT variants), NOT forest-type like SN: accumulate BA into 8 metagroups, pick the
dominant cover type (>50% BA, else mixed conifer), then per-cover-type rules build the
candidate `eqwt`, plus the always-added natural-fuel candidates 10 & 12; `_fmdyn` resolves
by the actual (SMALL,LARGE) fuel load. IN PROGRESS: MCCT (mixed conifer) fully ported +
validated (crt01); the biomass-heavy cover types (OB/PP/WS/SF/LP/AS) and the exact FMSSTAGE
structure class (`sstage.f`) are follow-ups — see docs/CR_VARIANT_PORT_AUDIT.md.
"""
function cr_select_fuel_models(s::StandState, mois::AbstractMatrix{Float32}, sm::Float32, lg::Float32)
    t = s.trees; fs = s.fire
    percov = fs.percov
    fwind = fs.swind * fire_wind_reduction(percov)
    ldry = false                                   # DROUGHT (IDRYB..IDRYE) — none in crt01
    NSP = nspecies(s.variant)                       # CR=38, UT=24, TT=18 (fmcfmd shared, VARACD-branched)
    usht = 0.5f0 * stand_top_height(s)             # UNDERSTORY = HT ≤ 0.5·FMAVH (top-40 ht) (fmcfmd.f:182)
    ctba = zeros(Float32, 8)                        # per-cover-type BA (CTBA)
    usba = zeros(Float32, 8)                        # per-cover-type UNDERSTORY BA (USBA, HT≤USHT; for LCUNDR)
    fmtba = zeros(Float32, NSP)                     # per-species BA (FMTBA — LPPDOM)
    sumtpa = 0f0; sumd = 0f0; sumd2 = 0f0           # ΣFMPROB, Σ(FMPROB·DBH), Σ(FMPROB·DBH²) — avg DBH / QMD
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        sp = Int(t.species[i]); d = t.dbh[i]
        x = t.tpa[i] * d * d * 0.0054542f0
        ct = _fm_covtype(s.variant, sp)
        ctba[ct] += x
        t.height[i] <= usht && (usba[ct] += x)     # understory BA (fmcfmd.f:203-231)
        (1 <= sp <= NSP) && (fmtba[sp] += x)
        sumtpa += t.tpa[i]; sumd += t.tpa[i] * d; sumd2 += t.tpa[i] * d * d
    end
    # LCUNDR: coniferous understory present (fmcfmd.f:368) — PPCT+WSCT+SFCT+LPCT+MCCT understory BA > 1
    lcundr = (usba[3] + usba[4] + usba[5] + usba[6] + usba[7]) > 1f0
    avgdbh = sumtpa > 1f-6 ? sumd / sumtpa : 0f0
    qmd = sumtpa > 1f-6 ? sqrt(sumd2 / sumtpa) : 0f0
    stndba = sum(ctba)
    # dominant cover-type metagroup: first > 50% BA, else mixed conifer (MCCT=7) (fmcfmd.f:331-343)
    ict = 7
    if t.n > 0 && stndba > 0.001f0
        for i in 1:8
            if ctba[i] / stndba > 0.5f0
                ict = i; break
            end
        end
    end
    # LPPDOM: is ponderosa the single highest-BA species? (fmcfmd.f:811-822; CR J=13, UT/TT J=10)
    pp_sp = _fm_ppct_sp(s.variant)
    lppdom = true
    @inbounds for i in 1:NSP
        (i != pp_sp && fmtba[i] > fmtba[pp_sp]) && (lppdom = false)
    end
    # IFMST structure class = FMSSTAGE (sstage.f) — already ported + validated bit-exact vs the SSTAGE
    # "Structural statistics" report as `structure_class`. Called with the CR fmcfmd params
    # (fmcfmd.f:388-398): GAPPCT=20, SSDBH=5, SAWDBH=18 (12 for lodgepole ICT=6), CCMIN=5, TPAMIN=200,
    # PCTSMX=30 — distinct from the STRCLASS-report defaults (gappct=30, sawdbh=25). (The crown-width fix
    # in fmcba/structure_class — cr_cwcalc, not the 0.5-default crown_width — makes this return the right
    # class, e.g. crt01 MCCT → IFMST 3 → model 10.)
    sawdbh = ict == 6 ? 12f0 : 18f0
    ifmst = structure_class(s; thresh = (20f0, 5f0, sawdbh, 5f0, 200f0, 30f0)).class
    imodty = Int(s.plot.model_type)
    # USCT = dominant UNDERSTORY cover-type group (fmcfmd.f:348-364): argmax(USBA), else ICT.
    usct = ict; ybest = -1f0; usbatot = 0f0
    for i in 1:8; usbatot += usba[i]; end
    if usbatot > 0f0
        for i in 1:8
            usba[i] > ybest && (usct = i; ybest = usba[i])
        end
    end
    eqwt = zeros(Float32, _FMD_ICLSS)
    # GOTO-111 loopback (fmcfmd.f:407): some cover-type rules reassign ICT and re-dispatch.
    # eqwt is NOT reset on loopback (the triggering branch sets nothing before GOTO 111), so a
    # loopback is simply re-running the dispatch with the new ICT. Bounded pass count = backstop.
    for _pass in 1:5
    redo_ict = 0                                   # >0 ⇒ loop back and re-dispatch with this ICT
    if ict == 1                                    # OBCT oak-brush (fmcfmd.f:413-518, CR)
        if fwind <= 7f0
            eqwt[8] = 1f0                           # low wind ⇒ model 8 (the common case)
        elseif ctba[1] <= 0f0                       # CTBA(OBCT)==0 ⇒ model 5
            eqwt[5] = 1f0
        else
            p2t = _FM_P2T; v2t = coef_col(s.coef, :v2t)
            # X = BA(FMPROB)-weighted avg height of oak species (23–27). NOTE: live uses a buggy IND1(J)
            # index (species-sort-order dependent) that in dominant-oak stands reads oak trees anyway; we
            # compute the INTENDED avg oak height. In mixed stands with Y>50 (the only X-dependent path)
            # this can diverge from live's buggy X — a documented bounded residual (rare dead-oak case).
            xh = 0f0; psum = 0f0; bl = 0f0; bd = 0f0
            @inbounds for i in 1:t.n
                t.tpa[i] > 0f0 || continue
                spi = Int(t.species[i]); dd = t.dbh[i]; hh2 = t.height[i]; pr = t.tpa[i]
                (23 <= spi <= 27) && (xh += hh2 * pr; psum += pr)
                xv = crown_biomass(s, spi, dd, hh2, Int(t.crown_pct[i]))   # BL over ALL trees (NO USHT filter)
                bl += xv[1] * pr * p2t
                for j in 2:6
                    bl += (xv[j] + t.ffe_oldcrw[j - 1, i]) * pr * p2t
                end
                bl += pr * v2t[spi] * p2t * cr_snag_bole_cuft(s, spi, dd, hh2)  # bole; v2t RAW ⇒ ·P2T
            end
            x = psum > 1f-6 ? xh / psum : 0f0
            sn = fs.snags
            @inbounds for k in eachindex(sn.sp)                            # snags (NO height filter for OBCT)
                bd += sn.fallvol[k] * (sn.den_hard[k] + sn.den_soft[k])
            end
            bd += lg + sm                                                   # + LARGE + SMALL down-wood
            y = (bd + bl) > 1f-6 ? 100f0 * bd / (bd + bl) : 0f0
            if x <= 2f0 || y <= 50f0
                eqwt[5] = 1f0
            elseif x > 6f0 && y > 50f0
                eqwt[4] = 1f0
            else                                                           # 2<X≤6, Y>50: model 4/5 blend
                w4 = _fm_algslp2(x, 2f0, 6f0, 0f0, 1f0); eqwt[4] = w4; eqwt[5] = 1f0 - w4
            end
        end
    elseif ict == 2                                # PJCT pinyon-juniper (fmcfmd.f:530-543, CR) — pure PERCOV
        if percov <= 25f0
            eqwt[2] = 1f0
        elseif percov <= 35f0
            eqwt[2] = 1f0 - (percov - 25f0) / 10f0; eqwt[5] = 1f0 - (35f0 - percov) / 10f0
        elseif percov <= 45f0
            eqwt[5] = 1f0
        elseif percov <= 55f0
            eqwt[5] = 1f0 - (percov - 45f0) / 10f0; eqwt[6] = 1f0 - (55f0 - percov) / 10f0
        else
            eqwt[6] = 1f0
        end
    elseif ict == 4                                # WSCT white spruce (fmcfmd.f:671-690, CR)
        if percov <= 40f0
            eqwt[2] = 1f0
        else
            (ifmst >= 4 && avgdbh > 12f0) ? (eqwt[10] = 1f0) : (eqwt[8] = 1f0)
        end
    elseif ict == 5                                # SFCT spruce-fir (fmcfmd.f:695-727, CR)
        if ifmst == 0
            eqwt[2] = 1f0
        elseif ifmst == 1
            qmd > 1f0 ? (eqwt[5] = 1f0) : (eqwt[2] = 1f0)
        elseif ifmst == 2
            eqwt[8] = 1f0
        else                                       # IFMST 3–6
            eqwt[10] = 1f0
        end
    elseif ict == 6                                # LPCT lodgepole (fmcfmd.f:742-759, CR)
        if ifmst == 0 || ifmst == 1 || ifmst == 5
            eqwt[5] = 1f0
        elseif ifmst == 2
            eqwt[8] = 1f0
        else                                       # IFMST 3,4,6
            eqwt[10] = 1f0
        end
    elseif ict == 7                                # MCCT mixed conifer (fmcfmd.f:827-869, CR)
        if lppdom
            eqwt[9] = 1f0
        elseif ifmst == 0
            ldry ? (eqwt[1] = 1f0) : (eqwt[2] = 1f0)
        elseif ifmst == 1
            ldry ? (eqwt[6] = 1f0) : (eqwt[5] = 1f0)
        elseif ifmst == 2
            ldry ? (eqwt[6] = 1f0) : (eqwt[8] = 1f0)
        elseif ifmst == 3 || ifmst == 4 || ifmst == 6
            eqwt[10] = 1f0
        elseif ifmst == 5
            if percov >= 55f0
                eqwt[8] = 1f0
            elseif percov < 45f0
                eqwt[2] = 1f0
            else
                eqwt[2] = 1f0 - (percov - 45f0) / 10f0
                eqwt[8] = 1f0 - (55f0 - percov) / 10f0
            end
        end
    elseif ict == 8                                # ASCT aspen (fmcfmd.f:936-1020, CR)
        if ctba[8] / max(1f-3, stndba) > 0.80f0    # aspen-DOMINANT (>80% BA)
            # (VARACD≠UT ∧ ≠TT ∧ IMODTY==1) → model 2 (CR SW-mixed), ELSE model 5 (fmcfmd.f ASCT)
            (s.variant isa CentralRockies && imodty == 1) ? (eqwt[2] = 1f0) : (eqwt[5] = 1f0)
        elseif lcundr                              # conifer understory present ⇒ COVOLP crown cover CVR10
            # Sum crown area (sq ft/ac) of trees HT>10 EXCLUDING sp {12,16,20:35,38} (fmcfmd.f:988-1000),
            # then COVOLP (canopy cover %). CVR10>40 ⇒ model 8, else model 2.
            cccoef = Float64(s.control.cc_coef)
            _cr_ba = s.plot.basal_area; _cr_el = s.plot.elevation
            _cr_hi = _cr_hopkins(s.plot.latitude, s.plot.longitude, s.plot.elevation)
            area = 0.0
            @inbounds for i in 1:t.n
                (t.tpa[i] > 0f0 && t.height[i] > 10f0) || continue
                spi = Int(t.species[i])
                _fm_asct_excl(s.variant, spi) && continue    # per-variant CVR10 exclusion set
                cw = cr_cwcalc(spi, t.dbh[i], t.height[i], Float32(t.crown_pct[i]), _cr_ba, _cr_el, _cr_hi)
                area += Float64(cw)^2 * Float64(t.tpa[i]) * 0.785398
            end
            pccu = cccoef * (area / 43560.0)
            cvr10 = pccu > 5.0 ? 100.0 : (1.0 - exp(-pccu)) * 100.0
            eqwt[cvr10 > 40.0 ? 8 : 2] = 1f0
        elseif (ctba[3] + ctba[4] + ctba[5] + ctba[6] + ctba[7]) > 1f0   # PPCT+WSCT+SFCT+LPCT+MCCT
            redo_ict = 7                           # ⇒ reprocess as MCCT (GOTO 111)
        elseif ctba[1] > 1f0
            redo_ict = 1                           # ⇒ reprocess as OBCT
        elseif ctba[2] > 1f0
            redo_ict = 2                           # ⇒ reprocess as PJCT
        end
    elseif ict == 3                                # PPCT ponderosa (fmcfmd.f:562-665, CR)
        if percov > 60f0                           # PERCOV>60 branch (fmcfmd.f:643-665) — no understory biomass
            if fwind > 7f0
                if ctba[2] / max(1f-3, stndba) > 0.2f0    # PJCT pinyon-juniper component
                    eqwt[ldry ? 6 : 5] = 1f0
                elseif ifmst == 6
                    eqwt[lcundr ? 5 : 2] = 1f0
                else
                    eqwt[9] = 1f0
                end
            else
                eqwt[9] = 1f0
            end
        else                                       # PERCOV≤60 branch (fmcfmd.f:571-631) — understory BL/BD biomass
            p2t = _FM_P2T; v2t = coef_col(s.coef, :v2t)
            bl = 0f0; bd = 0f0
            @inbounds for i in 1:t.n
                (t.tpa[i] > 0f0 && t.height[i] <= usht) || continue
                spi = Int(t.species[i]); dd = t.dbh[i]; hh2 = t.height[i]
                xv = crown_biomass(s, spi, dd, hh2, Int(t.crown_pct[i]))   # (0..5) crown biomass (pounds)
                bl += xv[1] * t.tpa[i] * p2t                                # foliage (current crown only)
                for j in 2:6
                    bl += (xv[j] + t.ffe_oldcrw[j - 1, i]) * t.tpa[i] * p2t # wood j: current + OLDCRW(j-1)
                end
                bl += t.tpa[i] * v2t[spi] * p2t * cr_snag_bole_cuft(s, spi, dd, hh2)  # bole; v2t is RAW ⇒ ·P2T
            end
            sn = fs.snags
            @inbounds for k in eachindex(sn.sp)
                den = sn.den_hard[k] + sn.den_soft[k]
                (den > 0f0 && sn.height[k] <= usht) || continue
                bd += sn.fallvol[k] * den                              # snag biomass (fallvol=total-cuft·V2T)
            end
            y = (bd + bl) > 1f-6 ? 100f0 * bd / (bd + bl) : 0f0
            if (bd + bl) > 0f0
                if fwind > 7f0
                    if y <= 50f0
                        (usct == 1 || usct == 2) ? (redo_ict = usct) : (eqwt[5] = 1f0)  # OBCT/PJCT loopback
                    else
                        eqwt[6] = 1f0
                    end
                else
                    eqwt[5] = 1f0
                end
            else
                eqwt[2] = 1f0
            end
        end
    end
    redo_ict == 0 && break
    ict = redo_ict
    end                                            # end GOTO-111 loopback loop
    # ALL CR cover-type fuel-model rules are now ported: OBCT/PJCT/PPCT(both PERCOV branches)/WSCT/SFCT/
    # LPCT/MCCT/ASCT(dominant + LCUNDR-COVOLP + MCCT/OBCT/PJCT loopback), plus the GOTO-111 loopback infra
    # (PPCT FWIND>7 Y≤50 → USCT OBCT/PJCT; ASCT → MCCT/OBCT/PJCT). Bounded residuals (documented): the OBCT
    # FWIND>7 buggy-IND1(J) avg-oak-height X (only affects the rare Y>50 mixed-stand blend; dominant-oak
    # agrees), and BD snag CURRENT-broken-height (FMSVOL to HTIH/HTIS vs jl total fallvol — model-inert on
    # validated stands). SFCT/WSCT residuals near the model-8-vs-10 boundary track the F3 down-wood-fuel
    # magnitude (jl cwd pools), not these rules.
    # Always-added natural-fuel candidates (fmcfmd.f:1045-1046; AFWT=0 with no recent harvest).
    eqwt[10] = 1f0
    eqwt[12] = 1f0
    return _fmdyn(sm, lg, eqwt, fmd_xpts(s.variant))
end

# MAPDRY (ie/fmcba.f:82) habitat KODTYP → IDRY class (1=dry grassy, 2=dry shrubby); absent ⇒ 0 (other).
# IE and KT are bit-identical (verified). EM/BM/CI pass their own MAPDRY when their FFE lands.
const _IE_MAPDRY = Dict{Int,Int}(
    130 => 1, 140 => 1, 210 => 1, 220 => 1, 230 => 1,                       # dry grassy (PIPO/PSME grass)
    161 => 2, 170 => 2, 171 => 2, 172 => 2, 180 => 2, 181 => 2, 182 => 2,   # dry shrubby (PIPO/PSME shrub)
    310 => 2, 311 => 2, 312 => 2, 313 => 2)

"""IE-family FMCFMD candidate selection (ie/fmcfmd.f): MAPDRY habitat→dryness class → a base fuel model
weighted by canopy cover, plus the always-added natural-fuel candidates {10,12,13}. FMDYN then resolves
the (SMALL,LARGE) down-wood point among the weighted models. Activity fuels (models 11/14 via
AFWT/SLCHNG/LATFUEL) are DEFERRED — natural stands take LATFUEL=false ⇒ AFWT=0 (same as SN/NE/CR)."""
function ie_select_fuel_models(s::StandState, mois::AbstractMatrix{Float32}, sm::Float32, lg::Float32)
    percov = s.fire.percov
    eqwt = zeros(Float32, _FMD_ICLSS)
    idry = get(_IE_MAPDRY, Int(s.plot.habitat_code), 0)          # FMKOD → IDRY (fmcba.f MAPDRY)
    # Canopy-cover split: ALGSLP(PERCOV, X=[30,50], Y=[0,1]) — high-cover weight → model 9 (fmcfmd.f:61-77).
    wt9 = percov <= 30f0 ? 0f0 : percov >= 50f0 ? 1f0 : (percov - 30f0) / 20f0
    wt_lo = 1f0 - wt9
    if idry == 1                                                 # dry grassy: model 1 (low cover) + 9 (high)
        eqwt[1] = wt_lo; eqwt[9] = wt9
    elseif idry == 2                                             # dry shrubby: model 2 + 9
        eqwt[2] = wt_lo; eqwt[9] = wt9
    else                                                         # all other habitats: model 8
        eqwt[8] = 1f0
    end
    # Always-added natural-fuel candidates (ie/fmcfmd.f:95-98; AFWT=0, no recent harvest ⇒ EQWT 10/12=1-0).
    eqwt[10] = 1f0; eqwt[12] = 1f0; eqwt[13] = 1f0
    return _fmdyn(sm, lg, eqwt, fmd_xpts(s.variant))
end

"""EM FMCFMD candidate selection (em/fmcfmd.f): the two candidate models come from the habitat table
EMMD ((M1,M2)=MD1/MD2[IEMTYP]), split by PERCOV via ALGSLP([30,50]) — WT1(1)→M1, WT1(2)→M2 — plus the
natural-fuel candidates {10,12,13}. IEMTYP (the EM habitat subscript 1..122) is recomputed from the stand's
habitat code via the growth-port em_habtyp. Activity fuels (11/14) deferred like IE (AFWT=0 natural path)."""
function em_select_fuel_models(s::StandState, mois::AbstractMatrix{Float32}, sm::Float32, lg::Float32)
    percov = s.fire.percov
    eqwt = zeros(Float32, _FMD_ICLSS)
    iemtyp = em_habtyp(Int(s.plot.habitat_code))[1]              # EM habitat subscript (em/habtyp.f JTYPE bucket)
    m1, m2 = em_md_models(iemtyp)                               # EMMD: (MD1,MD2)[iemtyp]
    wt2 = percov <= 30f0 ? 0f0 : percov >= 50f0 ? 1f0 : (percov - 30f0) / 20f0   # ALGSLP(PERCOV,[30,50],[0,1])
    eqwt[m1] += 1f0 - wt2                                        # WT1(1) → M1 (low cover)
    eqwt[m2] += wt2                                              # WT1(2) → M2 (high cover)
    eqwt[10] = 1f0; eqwt[12] = 1f0; eqwt[13] = 1f0             # natural fuels ASSIGNED (overwrite M1/M2 if 10/12/13)
    return _fmdyn(sm, lg, eqwt, fmd_xpts(s.variant))
end

"""CI FMCFMD candidate selection (ci/fmcfmd.f) — cover-type-group (ICT=MAPPVG[ICINDX], 1..11) based, with a
grand-fir-understory sub-model for ICT 5:6. ICINDX is stashed in p.habitat_input by the CI site setup.
PRLONG = BA-fraction in long-needle pines (sp 1,10). Activity fuels (11/14) deferred (AFWT=0 natural path)."""
function ci_select_fuel_models(s::StandState, mois::AbstractMatrix{Float32}, sm::Float32, lg::Float32)
    t = s.trees; coef = s.coef; percov = s.fire.percov
    eqwt = zeros(Float32, _FMD_ICLSS)
    icindx = Int(s.plot.habitat_input); ict = ci_pvg(icindx)
    # PRLONG: BA-fraction in sp {1,10} (long-needle pines). FMTBA = per-species BA (ci/fmcba.f uses 0.0054542·D²·TPA).
    prlong = 0f0; stndba = 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        x = t.tpa[i] * t.dbh[i]^2 * 0.0054542f0; stndba += x
        (t.species[i] == 1 || t.species[i] == 10) && (prlong += x)
    end
    prlong = (prlong > 0.01f0 && stndba > 0.01f0) ? prlong / stndba : 0f0
    alg(v, x1, x2) = _cr_algslp2(Float32(v), Float32(x1), Float32(x2), 0f0, 1f0)   # ALGSLP(v,[x1,x2],[0,1])
    if 1 <= ict <= 4
        k = ict == 1 ? 1 : ict == 2 ? (ci_s9b(icindx) == 1 ? 5 : 2) : ict == 3 ? 5 : 2
        w2 = alg(percov, 30, 50); w1 = 1f0 - w2
        w1 > 0f0 && (eqwt[k] += w1)
        w2 > 0f0 && (eqwt[9] += w2 * prlong; eqwt[8] += w2 * (1f0 - prlong))
    elseif ict == 5 || ict == 6
        k = ict == 5 ? 2 : 5
        # grand-fir (sp 4) sapling (DBH≤3) understory: avg crown ratio CRGF (%) + crown-cover CCGF (%).
        crgf = 0f0; fmtpa = 0f0; totcra = 0f0
        @inbounds for i in 1:t.n
            (t.tpa[i] > 0f0 && t.species[i] == 4 && t.dbh[i] <= 3f0) || continue
            crgf += t.tpa[i] * Float32(t.crown_pct[i]); fmtpa += t.tpa[i]
            cw = crown_width(coef, s.species.code2[4], t.dbh[i], t.height[i], Float32(t.crown_pct[i]), 0,
                             s.plot.latitude, s.plot.longitude, s.plot.elevation)
            totcra += 3.1415927f0 * cw * cw / 4f0 * t.tpa[i]
        end
        crgf = fmtpa > 0f0 ? crgf / fmtpa : 0f0
        lcrgf = crgf >= 75f0
        ccgf = 100f0 * (1f0 - exp(-totcra / 43560f0))
        if lcrgf
            w2 = alg(ccgf, 50, 70); w1 = 1f0 - w2       # WT1 on CCGF
            if w1 > 0f0
                u2 = alg(percov, 40, 60); u1 = 1f0 - u2  # WT2 on PERCOV
                u1 > 0f0 && (eqwt[k] += w1 * u1)
                u2 > 0f0 && (eqwt[9] += w1 * u2 * prlong; eqwt[8] += w1 * u2 * (1f0 - prlong))
            end
            if w2 > 0f0
                if w2 == 1f0
                    eqwt[5] += w2
                else
                    u2 = alg(percov, 40, 60); u1 = 1f0 - u2
                    u1 > 0f0 && (eqwt[5] += w2 * u1)
                    u2 > 0f0 && (eqwt[9] += w2 * u2 * prlong; eqwt[8] += w2 * u2 * (1f0 - prlong))
                end
            end
        else                                              # no significant long-crown GF understory
            w2 = alg(percov, 40, 60); w1 = 1f0 - w2
            w1 > 0f0 && (eqwt[k] += w1)
            w2 > 0f0 && (eqwt[9] += w2 * prlong; eqwt[8] += w2 * (1f0 - prlong))
        end
    else                                                  # ICT 7:11
        eqwt[8] = 1f0
    end
    eqwt[10] = 1f0; eqwt[12] = 1f0; eqwt[13] = 1f0        # natural fuels (AFWT=0 natural path)
    return _fmdyn(sm, lg, eqwt, fmd_xpts(s.variant))
end

# COVOLP (covolp.f): canopy cover % from a set of per-tree crown areas CRAREA over idx[lo:hi].
@inline function _bm_covolp(idx::Vector{Int}, lo::Int, hi::Int, wk6::Vector{Float32}, cccoef::Float32)::Float32
    hi < lo && return 0f0
    ssum = 0f0
    @inbounds for ii in lo:hi; ssum += wk6[idx[ii]]; end
    pccu = cccoef * (ssum / 43560f0)
    return pccu > 5f0 ? 100f0 : (1f0 - exp(-pccu)) * 100f0
end

"""
    bm_stage(s) -> (cova, covb, la)

BMSTAGE (bm/fmcfmd.f:337): a stripped SSTAGE that stratifies the stand into at most two canopy layers by
the largest height gap (>= max(10ft, 30% of the taller tree), ladder trees < 2 TPA absorbed), returns the
upper/lower stratum canopy cover (COVA/COVB, %) via COVOLP and a per-tree upper-layer membership flag LA.
CRWDTH is the western forest-grown crown width (bm_cwcalc -> cr_cwcalc library).
"""
function bm_stage(s::StandState)
    t = s.trees; n = t.n
    cccoef = Float32(s.control.cc_coef)
    la = falses(n)
    wk6 = zeros(Float32, n)
    _ba = s.plot.basal_area; _el = s.plot.elevation
    _hi = _cr_hopkins(s.plot.latitude, s.plot.longitude, s.plot.elevation)
    cwof(i) = bm_cwcalc(Int(t.species[i]), t.dbh[i], t.height[i], Float32(t.crown_pct[i]), _ba, _el, _hi)
    idx = Int[]; sprob = 0f0
    @inbounds for i in 1:n
        sprob += t.tpa[i]
        t.tpa[i] > 0.00001f0 && push!(idx, i)
    end
    ntrees = length(idx)
    ntrees == 0 && return (0f0, 0f0, la)
    ht = t.height
    crs1 = 0f0; crs2 = 0f0
    is1i1 = 1; is1i2 = ntrees; is2i1 = 0; is2i2 = ntrees
    if ntrees <= 1
        i = idx[1]; w = cwof(i)
        crs1 = w * w * t.tpa[i] * 0.785398f0 / 43560f0
        is1i1 = 1; is1i2 = 1
    else
        rdpsrt!(ntrees, ht, idx, false)                     # HT descending; idx[1] = tallest
        @inbounds for ii in 1:ntrees
            i = idx[ii]; w = cwof(i)
            wk6[i] = w * w * t.tpa[i] * 0.785398f0          # crown area (sqft/ac)
        end
        diff1 = -1f20; id1i1 = 0; id1i2 = 0
        iilg = 1; ilarge = idx[iilg]; sumprb = 0f0
        @inbounds for ii in 2:ntrees
            ismall = idx[ii]
            x = max(10f0, ht[ilarge] * 30f0 * 0.01f0)
            if ht[ismall] < ht[ilarge] - x
                if t.tpa[ismall] + sumprb < 2f0
                    sumprb += t.tpa[ismall]
                else
                    dff = ht[ilarge] - ht[ismall]
                    if dff > diff1
                        diff1 = dff; id1i1 = iilg; id1i2 = ii
                    end
                    ilarge = ismall; iilg = ii; sumprb = 0f0
                end
            else
                if t.tpa[ismall] + sumprb < 2f0
                    sumprb += t.tpa[ismall]
                else
                    ilarge = ismall; iilg = ii; sumprb = 0f0
                end
            end
        end
        nstr = 1
        if id1i1 > 0
            nstr = 2; is1i2 = id1i1; is2i1 = id1i2; is2i2 = ntrees
        end
        is1i2 = max(is1i2, is2i1 - 1)
        crs1 = _bm_covolp(idx, is1i1, is1i2, wk6, cccoef)
        is1ok = crs1 > 5f0 ? 1 : 0
        is2ok = 0
        if nstr >= 2
            crs2 = _bm_covolp(idx, is2i1, is2i2, wk6, cccoef)
            is2ok = crs2 > 5f0 ? 1 : 0
        end
        nstr = is1ok + is2ok
        if nstr == 0 && sprob >= 200f0
            crs2 = 0f0; is1i1 = 1; is1i2 = ntrees
            crs1 = _bm_covolp(idx, is1i1, is1i2, wk6, cccoef)
        end
    end
    @inbounds for ii in is1i1:is1i2
        (1 <= ii <= ntrees) && (la[idx[ii]] = true)
    end
    return (crs1, crs2, la)
end

# BMSZCLS complexity index per composite size category (bm/fmcfmd.f C(13)).
const _BM_SZ_CC = (1,1,1,1,1,1,1, 2,2,2, 3,3,2)

"""
    bm_szcls(s, la) -> Vector{Float32}(13)

BMSZCLS (bm/fmcfmd.f:549): BA-weighted dominant size-class distribution over the 13 BM size categories,
Monte-Carlo averaged over 50 passes of +/-20% jittered DBH (upper-layer trees LA only). Each pass bins the
jittered BA into 7 basic + 6 composite categories, picks the max-BA category (ties -> lowest complexity),
and adds 1/50 to that category's WDOM. RANN is bracketed by RANNGET/RANNPUT so it consumes ZERO net RNG.
"""
function bm_szcls(s::StandState, la::AbstractVector{Bool})
    t = s.trees; n = t.n
    wdom = zeros(Float32, 13)
    pass = 50; jitter = 0.2f0
    saveso = rannget(s.rng)
    for _j in 1:pass
        sa = zeros(Float32, 13)
        @inbounds for i in 1:n
            la[i] || continue
            xran = rann!(s.rng)
            dt = t.dbh[i] * (1f0 + jitter * ((xran * 2f0) - 1f0))
            ba = dt * dt * t.tpa[i]
            if dt < 1f0;        sa[1] += ba
            elseif dt < 5f0;    sa[2] += ba
            elseif dt < 9f0;    sa[3] += ba
            elseif dt < 15f0;   sa[4] += ba
            elseif dt < 21f0;   sa[5] += ba
            elseif dt < 32f0;   sa[6] += ba
            else                sa[7] += ba
            end
        end
        sa[8]  = sa[1] + sa[2]
        sa[9]  = sa[2] + sa[3]
        sa[10] = sa[4] + sa[5]
        sa[11] = sa[3] + sa[10]
        sa[12] = sa[10] + sa[6]
        sa[13] = sa[6] + sa[7]
        icls = 0; ccls = 4; swt = -1f0
        @inbounds for i in 1:13
            if sa[i] == swt
                _BM_SZ_CC[i] < ccls && (ccls = _BM_SZ_CC[i]; icls = i)
            elseif sa[i] > swt
                swt = sa[i]; ccls = _BM_SZ_CC[i]; icls = i
            end
        end
        icls >= 1 && (wdom[icls] += 1f0 / Float32(pass))
    end
    rannput!(s.rng, saveso)
    return wdom
end

const _BM_SCLAB = (1f0,2f0,3f0,4f0,5f0,6f0,6.5f0,7f0,7.5f0,8f0,9f0,10f0,11f0)

"""BM FMCFMD candidate selection (bm/fmcfmd.f). WT1 splits the stand by max(PRDF,PRPP) (DF sp3 / PP sp10
BA fraction) over ALGSLP([0.40,0.60]). WT1(1) (neither dominates) drives a size-class-distribution path
(BMSTAGE + BMSZCLS) weighting models 5/8/10 by size class x PERCOV x stratum cover; WT1(2) (DF/PP dominant)
drives a PERCOV-band path (models 1/2/9/10 + a PRPP/(PRPP+PRDF) split). Natural fuels {10,12,13} always
ASSIGNED last (overwriting any 10 accumulation). Activity fuels (11/14) deferred (AFWT=0)."""
function bm_select_fuel_models(s::StandState, mois::AbstractMatrix{Float32}, sm::Float32, lg::Float32)
    t = s.trees; fs = s.fire
    percov = fs.percov
    fmtba = zeros(Float32, 18); stndba = 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        sp = Int(t.species[i]); (1 <= sp <= 18) || continue
        x = t.tpa[i] * t.dbh[i] * t.dbh[i] * 0.0054542f0
        fmtba[sp] += x; stndba += x
    end
    prdf = stndba > 0.01f0 ? fmtba[3] / stndba : 0f0
    prpp = stndba > 0.01f0 ? fmtba[10] / stndba : 0f0
    eqwt = zeros(Float32, _FMD_ICLSS)
    wt1b = _fm_algslp2(max(prdf, prpp), 0.40f0, 0.60f0, 0f0, 1f0)   # WT1(2)
    wt1a = 1f0 - wt1b                                               # WT1(1)

    if wt1a > 0f0                                                   # CASE 1: neither DF nor PP dominates
        cova, covb, la = bm_stage(s)
        wd = bm_szcls(s, la)
        @inbounds for k in 1:13
            wd[k] > 0f0 || continue
            szcls = _BM_SCLAB[k]
            if szcls <= 3f0
                eqwt[5] += wt1a * wd[k]
            elseif szcls < 7f0
                w2b = _fm_algslp2(percov, 25f0, 35f0, 0f0, 1f0); w2a = 1f0 - w2b
                w2a > 0f0 && (eqwt[5] += wt1a * w2a * wd[k])
                w2b > 0f0 && (eqwt[8] += wt1a * w2b * wd[k])
            else                                                    # szcls > 7
                w2b = _fm_algslp2(percov, 25f0, 35f0, 0f0, 1f0); w2a = 1f0 - w2b
                w2a > 0f0 && (eqwt[5] += wt1a * w2a * wd[k])
                if w2b > 0f0
                    w3b = _fm_algslp2(cova, 25f0, 35f0, 0f0, 1f0); w3a = 1f0 - w3b
                    w4b = _fm_algslp2(covb, 10f0, 20f0, 0f0, 1f0); w4a = 1f0 - w4b
                    eqwt[8]  += wt1a * w2b * w3a * w4a * wd[k]
                    eqwt[8]  += wt1a * w2b * w3a * w4b * wd[k]
                    eqwt[8]  += wt1a * w2b * w3b * w4a * wd[k]
                    eqwt[10] += wt1a * w2b * w3b * w4b * wd[k]
                end
            end
        end
    end

    if wt1b > 0f0                                                  # CASE 2: DF or PP dominates BA
        wt2 = zeros(Float32, 4)
        if percov < 15f0
            v = _fm_algslp2(percov, 5f0, 15f0, 0f0, 1f0); wt2[2] = v; wt2[1] = 1f0 - v
        elseif percov < 35f0
            v = _fm_algslp2(percov, 25f0, 35f0, 0f0, 1f0); wt2[3] = v; wt2[2] = 1f0 - v
        else
            v = _fm_algslp2(percov, 45f0, 55f0, 0f0, 1f0); wt2[4] = v; wt2[3] = 1f0 - v
        end
        wt2[1] > 0f0 && (eqwt[1] += wt1b * wt2[1])
        wt2[2] > 0f0 && (eqwt[2] += wt1b * wt2[2])
        if wt2[3] > 0f0
            denom = prpp + prdf
            w3b = _fm_algslp2(denom > 0f0 ? prpp / denom : 0f0, 0.40f0, 0.60f0, 0f0, 1f0)
            w3a = 1f0 - w3b
            w3a > 0f0 && (eqwt[10] += wt1b * wt2[3] * w3a)
            w3b > 0f0 && (eqwt[2]  += wt1b * wt2[3] * w3b)
        end
        wt2[4] > 0f0 && (eqwt[9] += wt1b * wt2[4])
    end

    # Natural-fuel candidates always assigned (bm/fmcfmd.f:318-320; AFWT=0 => 10/12=1, 13=1).
    eqwt[10] = 1f0; eqwt[12] = 1f0; eqwt[13] = 1f0
    return _fmdyn(sm, lg, eqwt, fmd_xpts(s.variant))
end

"""
    _fmdyn(sm, lg, eqwt) -> Vector{Tuple{Int,Float32}}

FMDYN (fmdyn.f): resolve candidate fuel models (`eqwt[i] > 0`) into weighted models by
the inverse perpendicular distance from the point `(sm, lg)` to each model's iso-line.
All SN iso-lines are sloped (ITYP≡0), so only the sloped-line geometry is ported.
Collinear candidates (identical XPTS — e.g. the litter models 1–9) share their bracket's
weight in proportion to `eqwt`. Returns up to `MXFMOD` (model, weight) pairs summing to 1.
"""
function _fmdyn(sm::Float32, lg::Float32, eqwt::Vector{Float32}, xpts::AbstractMatrix{Float32} = _FMD_XPTS;
               iptr::Union{Nothing,Vector{Int}} = nothing)
    ic = length(eqwt); mx = _FMD_MXFMOD
    out = Tuple{Int,Float32}[]
    (sm < 0f0 || lg < 0f0) && return out

    lok = falses(ic)
    for i in 1:ic
        eqwt[i] > 0f0 && (lok[i] = true)
    end
    # unset candidates with a zero intercept (degenerate line) — none in the SN table
    for i in 1:ic
        if lok[i] && (xpts[i, 1] == 0f0 || xpts[i, 2] == 0f0)
            lok[i] = false; eqwt[i] = 0f0
        end
    end
    # EQMOD: tag each candidate with the first candidate sharing its iso-line (collinear)
    eqmod = zeros(Int, ic)
    for i in 1:ic
        lok[i] || continue
        for j in i:ic
            if eqmod[j] == 0 && lok[j] && xpts[i, 1] == xpts[j, 1] && xpts[i, 2] == xpts[j, 2]
                eqmod[j] = i
            end
        end
    end
    # rescale each collinear group's eqwt to sum to 1
    for i in 1:ic
        lok[i] || continue
        xwt = 0f0
        for j in i:ic
            (lok[j] && eqmod[j] == i) && (xwt += eqwt[j])
        end
        if xwt > 1f-6
            for j in i:ic
                (lok[j] && eqmod[j] == i) && (eqwt[j] /= xwt)
            end
        end
    end
    # XD/YD: signed distance from the point to where each line crosses its LARGE / SMALL
    xd = zeros(Float32, ic); yd = zeros(Float32, ic)
    for i in 1:ic
        lok[i] || continue
        m1 = xpts[i, 2] / (-xpts[i, 1]); b1 = xpts[i, 2]
        xd[i] = (lg - b1) / m1 - sm
        yd[i] = (m1 * sm + b1) - lg
    end
    # nearest left/right (xd) and below/above (yd) candidate lines
    nbr = zeros(Int, 4)
    prv = Float32[-9.99f30, 9.99f30, -9.99f30, 9.99f30]
    for i in 1:ic
        lok[i] || continue
        if xd[i] < 0f0 && xd[i] > prv[1]
            prv[1] = xd[i]; nbr[1] = i
        elseif xd[i] >= 0f0 && xd[i] < prv[2]
            prv[2] = xd[i]; nbr[2] = i
        end
        if yd[i] < 0f0 && yd[i] > prv[3]
            prv[3] = yd[i]; nbr[3] = i
        elseif yd[i] >= 0f0 && yd[i] < prv[4]
            prv[4] = yd[i]; nbr[4] = i
        end
    end
    # perpendicular distance from the point to each neighbor line
    wt = zeros(Float32, 4)
    for k in 1:4
        i = nbr[k]; (i == 0 || !lok[i]) && continue
        m1 = xpts[i, 2] / (-xpts[i, 1]); b1 = xpts[i, 2]
        m2 = -(1f0 / m1); b2 = lg - m2 * sm
        npt1 = (b2 - b1) / (m1 - m2); npt2 = m2 * npt1 + b2
        wt[k] = sqrt((lg - npt2)^2 + (sm - npt1)^2)
    end
    # merge duplicate neighbors, accumulating distance into fmod/fwt
    fmod = zeros(Int, mx); fwt = zeros(Float32, mx); k2 = 0
    for i in 1:4
        nbr[i] == 0 && continue
        k2 += 1; k = k2; found = false
        for j in 1:i
            if fmod[j] == nbr[i]
                k = j; found = true; break
            end
        end
        !found && k <= mx && (fmod[k] = nbr[i])
        k <= mx && (fwt[k] += wt[i])
    end
    # weight by inverse distance, normalize
    xwt = 0f0
    for i in 1:mx
        fmod[i] == 0 && continue
        fwt[i] = 1f0 / (fwt[i] + 1f-6); xwt += fwt[i]
    end
    for i in 1:mx
        fmod[i] != 0 && (fwt[i] /= xwt)
    end
    # compact nonzero entries to the top
    k = 0
    for i in 1:mx
        fmod[i] == 0 && continue
        k += 1
        if i != k && k <= mx
            fmod[k] = fmod[i]; fwt[k] = fwt[i]; fmod[i] = 0; fwt[i] = 0f0
        end
    end
    # split each bracket's weight among its collinear models in proportion to eqwt
    fmod2 = zeros(Int, mx); fwt2 = zeros(Float32, mx); k = 1
    for i in 1:mx
        fmod[i] == 0 && continue
        ii = fmod[i]; xw = fwt[i]
        if eqmod[ii] == 0 && k <= mx
            fwt2[k] = fwt[i]; fmod2[k] = fmod[i]; k += 1
        else
            for j in 1:ic
                if eqmod[j] == eqmod[ii] && eqwt[j] > 0f0 && k <= mx
                    fwt2[k] += eqwt[j] * xw; fmod2[k] = j; k += 1
                end
            end
        end
    end
    # merge any duplicate models that arose from collinear additions
    fmod = zeros(Int, mx); fwt = zeros(Float32, mx); k2 = 0
    for i in 1:mx
        fmod2[i] == 0 && continue
        k2 += 1; k = k2
        for j in 1:i
            if fmod[j] == fmod2[i]
                k = j; break
            end
        end
        k <= mx && fmod[k] == 0 && (fmod[k] = fmod2[i])
        k <= mx && (fwt[k] += fwt2[i])
    end
    for i in 1:mx
        fmod[i] != 0 && push!(out, (iptr === nothing ? fmod[i] : iptr[fmod[i]], fwt[i]))
    end
    return out
end

"""
    build_dynamic_fuel_model(s, mois) -> (load, sav, depth, mext)

Construct the SN dynamic surface fuel model (FMCFMD3, fmcfmd2.f) from the stand's fuel
state. `load[2,4]`/`sav[2,4]` are loads (lb/ft²) and surface-area-to-volume by
[1=dead/2=live, class], `depth` the fuel-bed depth (ft), `mext` the dead moisture of
extinction. Dead loads come from the down-wood pools (`fire.cwd`, with the 1-hr class =
0–.25" + litter), the live-woody load from the understory crown biomass (foliage +
½·fine for trees ≤ CANMHT) plus the live shrub, and the live-herb load from `fire.flive`.
A moisture-dependent share of the live herb is moved to a dead-herb class (fmgfmv.f).
`mois` is the fuel-moisture matrix from `fuel_moisture`.
"""
function build_dynamic_fuel_model(s::StandState, mois::AbstractMatrix{Float32})
    fs = s.fire; t = s.trees
    # down-wood pools → load by size class (tons/ac → lb/ft²)
    currcwd = zeros(Float32, 11)
    @inbounds for j in 1:11, k in 1:2, l in 1:4
        currcwd[j] += fs.cwd[j, k, l] * _TONS_TO_LBFT2
    end
    herb = fs.flive[1] * _TONS_TO_LBFT2
    # understory live-woody load: crown foliage + ½ of the 0–.25" crown for trees ≤ CANMHT
    woody = 0f0
    @inbounds for i in 1:t.n
        (t.tpa[i] > 0f0 && t.height[i] <= _FM_CANMHT) || continue
        xv = crown_biomass(s, t.species[i], t.dbh[i], t.height[i], Int(t.crown_pct[i]))
        woody += (xv[1] + 0.5f0 * xv[2]) * t.tpa[i] * _FM_P2T    # ×P2T undoes crown_biomass's /P2T
    end
    woody = (woody + fs.flive[2]) * _TONS_TO_LBFT2

    load = zeros(Float32, 2, 4); sav = zeros(Float32, 2, 4)
    load[1, 1] = max(0f0, currcwd[1] + currcwd[10])  # 1-hr: 0–.25" + litter
    load[1, 2] = max(0f0, currcwd[2])                # 10-hr
    load[1, 3] = max(0f0, currcwd[3])                # 100-hr
    load[2, 1] = max(0f0, woody)                     # live woody
    load[2, 2] = max(0f0, herb)                      # live herb
    sav[1, 1] = _FM_USAV[1]; sav[1, 2] = 109f0; sav[1, 3] = 30f0
    sav[2, 1] = _FM_USAV[3]; sav[2, 2] = _FM_USAV[2]

    # dead-herb split: when herb is dry (moisture < 1.2) move part of the live herb into a
    # dead-herb class, which takes the live-herb SAV (fmgfmv.f:79/88-97).
    if load[2, 2] > 0f0 && mois[2, 2] < 1.2f0
        wt = _fm_algslp2(mois[2, 2], 0.30f0, 1.2f0, 0f0, 1f0)
        load[1, 4] = (1f0 - wt) * load[2, 2]; sav[1, 4] = sav[2, 2]
        load[2, 2] = wt * load[2, 2]
    end

    # fuel-bed depth from a load-weighted bulk density; moisture of extinction (fmcfmd2.f:582)
    fdfl = currcwd[1] + currcwd[10]
    ffl  = fdfl + herb + woody
    wf   = ffl > 0f0 ? fdfl / ffl : 0f0
    bdavg = _FM_UBD[1] + wf * (_FM_UBD[2] - _FM_UBD[1])
    depth = bdavg > 0f0 ? (ffl + currcwd[2] + currcwd[3]) / bdavg : 0f0
    mext  = (12f0 + 480f0 * bdavg / 32f0) / 100f0
    return (load, sav, depth, mext)
end
