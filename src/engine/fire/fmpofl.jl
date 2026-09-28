# =============================================================================
# fire/fmpofl.jl — the Potential Fire Report (FMPOFL, fire/vbase/fmpofl.f) for every FFE variant
#
# FMPOFL (called from FMMAIN every FFE year after FMBURN/FMPOCR, fmmain.f:188-194) burns two hypothetical
# fires — SEVERE (FMOIS=1, PREWND(1)/POTEMP(1)) and MODERATE (FMOIS=3, PREWND(2)/POTEMP(2)) — without touching
# the stand, and reports (DBSFMPF / DBSFMPFC):
#   surface + total flame, fire type, torching probability, torching/crowning index (FMCFIR, not SN/CS),
#   crown base height ACTCBH + crown bulk density CBD (FMPOCR), potential mortality (FMEFF ICALL=1: % of basal
#   area and the cubic volume killed), potential PM2.5 smoke (FMCONS ICALL=1) and the fuel models/weights.
# The routines it calls are mirrored here in their FMPOFL shape (FMFINT's weighted intermediates, FMCFIR,
# FMEFF's ICALL=1 accumulation, FMCONS's ICALL=1 smoke, FMPOFL_FMPTRH); they reuse the shared primitives
# (rothermel_surface_fire, fmgfmv, fire_tree_mortality, …) and never mutate the stand.
# =============================================================================

# fmpofl.f:118-119 SWIND = INT(PREWND(IDPL)); POTF* keyword overrides (fs.params.potf) replace the defaults.
const _POFL_FOLMC = 100f0                         # fminit.f:150 FOLMC (foliar moisture, %)

"""
    _pofl_fmfint(s, models, mois, fwind, dpmod) -> NamedTuple

FMFINT's weighted pass over the selected fuel models (fmfint.f:73-419, the FTYP≠2 / ICALL=2 branch): Byram
intensity, flame from the WEIGHTED Byram, and the weighted Rothermel intermediates SSIGMA/SRHOBQ/SXIR/SIRXI/
SPHIS/SFRATE/SCBE, plus HPA = SXIR·384/SSIGMA — each accumulated `x·FWT` in model order, in REAL*4.
"""
function _pofl_fmfint(s::StandState, models, mois::AbstractMatrix{Float32}, fwind::Float32, dpmod::Float32)
    byram = 0f0; ssig = 0f0; srho = 0f0; sxir = 0f0; sirxi = 0f0; sphis = 0f0; sfrate = 0f0; scbe = 0f0
    for (fm, w) in models
        (fm == 0 || w <= 0f0) && continue                              # fmfint.f:89
        load, sav, depth, mext = fmgfmv(s, fm, mois)
        depth *= dpmod
        r = rothermel_surface_fire(load, sav, depth, mext, mois; wind = fwind, slope_tan = s.plot.slope)
        byram += r.byram * w; ssig += r.sigma * w; srho += r.rhobqig * w; sxir += r.xir * w
        sirxi += r.xio * w; sphis += r.phis * w; sfrate += r.spread * w; scbe += r.scbe * w
    end
    flame = 0.45f0 * fpow(byram / 60f0, 0.46f0)                         # fmfint.f:810
    hpa = ssig != 0f0 ? sxir * 384f0 / ssig : 0f0
    return (; byram, flame, ssigma = ssig, srhobq = srho, sxir, sirxi, sphis, sfrate, scbe, hpa)
end

# FMFINT(…, FTYP=2, …, ICALL=2): only the weighted spread SFRATE(2) is used by FMCFIR's bisection/passive steps.
_pofl_sfrate(s::StandState, models, mois, fwind::Float32, dpmod::Float32) =
    _pofl_fmfint(s, models, mois, fwind, dpmod).sfrate

"""
    _pofl_fmcfir(s, models, mois, f, swind, wmult, cbd, actcbh, dpmod) -> NamedTuple

FMCFIR (fire/base/fmcfir.f, byte-identical in every variant build) for one FMPOFL scenario: the crowning index
OACT1 from the FM10 reference model at the fixed midflame wind SWIND·0.4, the torching index OINIT1 (analytic
start, then fmcfir.f's 1000-step bisection of the weighted surface spread against RINIT1), the fire type and the
crown fraction burned CRBURN / final spread RFINAL. `f` is the scenario's FMFINT (slot FMOIS) and `swind` the
INTEGER 20-ft wind.
"""
function _pofl_fmcfir(s::StandState, models, mois::AbstractMatrix{Float32}, f, swind::Int, wmult::Float32,
                      cbd::Float32, actcbh::Integer, dpmod::Float32)
    cbd <= 0f0 && return (; cftype = "SURFACE", crburn = 0f0, oinit = -1f0, oact = -1f0, rfinal = f.sfrate)
    init1 = 0f0
    if actcbh != -1
        init1 = (460f0 + 25.9f0 * _POFL_FOLMC) * 0.001333f0                          # fmcfir.f:94
        init1 = fpow(init1 * Float32(actcbh), 3f0 / 2f0)
    end
    # FM10 at FWIND = SWIND·0.4 (fmcfir.f:122-143: FMFINT(…,2,HPA2,1) ⇒ FMGFMV(IYR,10), no interpolation)
    load, sav, depth, mext = _fm10(s)
    r10 = rothermel_surface_fire(load, sav, depth, mext, mois; wind = Float32(swind) * 0.4f0, slope_tan = s.plot.slope)
    b = 0.02526f0 * fpow(f.ssigma, 0.54f0)                                            # fmcfir.f:146
    oact = 0f0
    if r10.xio < 1f0                                                                   # `IF (SIRXI(2) .LT. 00001)` — integer 1
        oact = -1f0
    else
        oact = ((2.95f0 * r10.rhobqig / (r10.xio * cbd)) - r10.phis - 1f0) / 0.001612f0
        oact = oact > 0f0 ? fpow(oact, 0.7f0) * 0.01137f0 / 0.4f0 : 0f0
    end
    ract = 3.34f0 * r10.spread
    rinit1 = f.hpa > 0f0 ? 60f0 * init1 / f.hpa : 0f0
    oinit = -1f0
    if actcbh != -1 && f.hpa > 0f0
        oinit = ((60f0 * init1 * f.srhobq / (f.hpa * f.sirxi)) - f.sphis - 1f0) / f.scbe  # fmcfir.f:180
        oinit = oinit > 0f0 ? fpow(oinit, 1f0 / b) * 0.01137f0 / wmult : 0f0
        sfr = _pofl_sfrate(s, models, mois, 999f0 * wmult, dpmod)
        if sfr < rinit1
            oinit = 999f0
        else
            oinit >= 999f0 && (oinit = 999f0)
            diff = _pofl_sfrate(s, models, mois, oinit * wmult, dpmod) - rinit1
            if !(diff <= 0.001f0 && diff >= -0.001f0)
                boundl = 0f0; boundu = 0f0
                for i in 1:1000                                                        # fmcfir.f DO 200
                    if i == 1
                        diff > 0.001f0 ? (boundl = 0f0; boundu = oinit) : (boundl = oinit; boundu = 999f0)
                    end
                    oinit = (boundu + boundl) / 2f0
                    diff = _pofl_sfrate(s, models, mois, oinit * wmult, dpmod) - rinit1
                    (diff <= 0.001f0 && diff >= -0.001f0) && break
                    diff > 0.001f0 && (boundu = oinit)
                    diff < -0.001f0 && (boundl = oinit)
                    boundu <= 0f0 && break
                    if boundu > 0f0 && boundu < 0.00000000001f0
                        oinit = 0f0; break
                    end
                end
            end
        end
        oinit = min(oinit, 999f0)
    end
    sw = Float32(swind)
    cftype = "SURFACE"; crburn = 0f0; rfinal = f.sfrate
    if oinit > sw
        oact > sw ? (cftype = "SURFACE"; crburn = 0f0; rfinal = f.sfrate) :
                    (cftype = "COND_CRN"; crburn = 1f0; rfinal = ract)
    elseif oact > sw
        cftype = "PASSIVE"
    else
        cftype = "ACTIVE"; crburn = 1f0; rfinal = ract
    end
    if oinit == -1f0 || oact == -1f0
        cftype = "SURFACE"; crburn = 0f0; rfinal = f.sfrate
    end
    if cftype == "PASSIVE"                                                             # fmcfir.f:347-359
        sfr2 = _pofl_sfrate(s, models, mois, oact * wmult, dpmod)
        cfb = (f.sfrate - rinit1) / (sfr2 - rinit1)
        rfinal = f.sfrate + cfb * (ract - f.sfrate)
        crburn = min(cfb, 1f0)
    end
    return (; cftype, crburn, oinit, oact, rfinal)
end

"""
    _pofl_fmeff(s, flame, sch, crburn, burnseas, psburn, year, cyclen) -> (pomort, pvolkl, bcrown)

FMEFF with ICALL=1 (fmeff.f): per record a RANN draw against PSBURN (the stream is restored, RANNGET/RANNPUT), the
FMEFF mortality PMORT for burned records, POMORT = PMORT·FMPROB plus the crown-fire share CRBURN·(FMPROB−POMORT),
accumulated as BA (DBH/24)² and volume (MCFV for CS/LS/NE/SN, CFV elsewhere). Returns POMORT = BAMORT/TOTBA, the
volume killed PVOLKL and the potential crown material burned PBRNCR (BCROWN, tons/ac) that FMCONS smokes.
"""
function _pofl_fmeff(s::StandState, flame::Float32, sch::Float32, crburn::Float32, burnseas::Integer,
                     psburn::Float32, year::Integer, cyclen::Real)
    fs = s.fire; t = s.trees; coef = s.coef
    bcrown = 0f0
    if crburn > 0f0                                                                    # fmeff.f:118-130
        @inbounds for isz in axes(fs.cwd2b, 2), idc in axes(fs.cwd2b, 1), itm in axes(fs.cwd2b, 3)
            bcrown += crburn * fs.cwd2b[idc, isz, itm] * psburn / 100f0 * _FM_P2T
            bcrown += crburn * fs.cwd2b2[idc, isz, itm] * psburn / 100f0 * _FM_P2T
        end
    end
    saved = rannget(s.rng)                                                             # fmeff.f:143
    bamort = 0f0; totba = 0f0; pvolkl = 0f0; pomort = 0f0
    merch = _fm_volkill_merch(s.variant)
    @inbounds for i in 1:t.n
        pmort = 0f0
        lpsburn = !(rann!(s.rng) * 100f0 > psburn)
        fmprob = t.tpa[i]; d = t.dbh[i]
        if lpsburn && fmprob > 0f0
            csv = crown_volume_scorched(sch, t.height[i], Int(t.crown_pct[i]))
            sp = Int(t.species[i])
            pmort = fire_tree_mortality(coef, sp, d, flame, csv, s.variant)
            pmort = fire_mortality_adjust(pmort, sp, d, burnseas, s.variant)
            (d <= 1f0 && csv > 50f0) && (pmort = 1f0)
            pmort *= active_fmort_mult(s.control, sp, year, d)
            pmort = clamp(pmort, 0f0, 1f0)
            bcrown += _fm_bcrown(s, i, crburn, sch, cyclen, true)
        end
        pomort = pmort * fmprob                                                        # fmeff.f:556-564
        lpsburn && (pomort = pomort + crburn * (fmprob - pomort))
        bamort = bamort + pomort * ((d / 24f0) * (d / 24f0))
        totba = totba + fmprob * ((d / 24f0) * (d / 24f0))
        pvolkl = pvolkl + (merch ? t.merch_cuft_vol[i] : t.cuft_vol[i]) * pomort
    end
    rannput!(s.rng, saved)
    totba != 0f0 && (pomort = bamort / totba)                                          # fmeff.f:603
    return (pomort, pvolkl, bcrown)
end

# FMCONS ICALL=1 (fmcons.f): the PM2.5 a fire would emit — consumed down wood × EMMFAC + live herb/shrub + PBRNCR,
# without consuming anything (the pools are restored). Returns PSMOKE·P2T (tons/ac), as DBSFMPF reports it.
function _pofl_smoke(s::StandState, mois::AbstractMatrix{Float32}, psburn::Float32, pbrncr::Float32)
    fs = s.fire
    cwd0 = copy(fs.cwd)
    c = fire_consumption!(fs, mois; psburn = psburn, burncr = pbrncr)
    copyto!(fs.cwd, cwd0)
    return c.smoke[1]
end

# FMPOFL_NPROB (fmpofl.f, Adams 1969 Algorithm 39): returns Q, the upper-tail probability (FMPOFL calls it as
# NPROB(Z,Q,PT,PDF), so PT receives the routine's Q).
function _pofl_nprob_upper(z::Float64)::Float64
    zabs = abs(z)
    zabs > 12.7 && return z < 0 ? 1.0 : 0.0
    y = 0.5 * z * z
    pdf = exp(-y) * 0.398942280385
    q = zabs > 1.28 ?
        pdf / (zabs - 3.8052e-8 + 1.00000615302 / (zabs + 3.98064794e-4 + 1.98615381364 / (zabs - 0.151679116635 +
               5.29330324926 / (zabs + 4.8385912808 - 15.1508972451 / (zabs + 0.742380924027 + 30.789933034 /
               (zabs + 3.99019417011)))))) :
        0.5 - zabs * (0.398942280444 - 0.399903438504 * y / (y + 5.75885480458 - 29.8213557808 /
               (y + 2.62433121679 + 48.6959930692 / (y + 5.92885724438))))
    return z < 0 ? 1.0 - q : q
end

"""
    _pofl_fmptrh(s, year, flm1, flm2) -> (ptr1, ptr2)

FMPOFL_FMPTRH (fmpofl.f:446-651): the probability of torching for the severe and moderate SURFACE flames. Crown
base heights sorted DESCENDING (RDPSRT), MOD(IYR,10) discarded draws, 30 virtual 0.025-ac plots (one RANN per
record in sorted order), the lowest ignitable crown per plot, then the NPROB upper tail on the log scale
(σ .25) averaged over the reps in mixed REAL/DOUBLE as FVS sums it. The RNG state is restored (RANNGET/RANNPUT).
"""
function _pofl_fmptrh(s::StandState, year::Integer, flm1::Float32, flm2::Float32, stash = nothing)
    prb, ht, icr = _pofl_fmmain_records(s.trees, stash)
    mxi = length(prb)
    (mxi <= 0 || (flm1 <= 0f0 && flm2 <= 0f0)) && return (0f0, 0f0)
    cbh = Float32[ht[i] * (1f0 - (Float32(icr[i]) * 0.01f0)) for i in 1:mxi]
    indx = zeros(Int, max(mxi, 30))
    rdpsrt!(mxi, cbh, view(indx, 1:mxi), true)
    saved = rannget(s.rng)
    for _ in 1:mod(Int(year), 10); rann!(s.rng); end
    avht = 0f0; ssum = 0f0
    @inbounds for i in 1:mxi                                                           # record order
        p = prb[i]
        ssum + p > 40f0 && (p = 40f0 - ssum)
        ssum = ssum + p
        avht = avht + ht[i] * p
        ssum >= 40f0 && break
    end
    ssum > 0f0 && (avht = avht / ssum)
    crit = max(5f0, min(0.5f0 * avht, 50f0))
    mincb = fill(-1f0, 30); yes = zeros(Int, mxi)
    for irep in 1:30
        itop = false; nyes = 0
        @inbounds for ii in 1:mxi
            i = indx[ii]
            ran = rann!(s.rng)
            if prb[i] > 1000f0 || ran > exp(-prb[i] * 0.025f0)
                nyes += 1; yes[nyes] = i
                ht[i] >= crit && (itop = true)
            end
        end
        mincb[irep] = -1f0
        itop || continue
        @inbounds for ii in nyes:-1:1
            i = yes[ii]
            if ht[i] >= crit
                mincb[irep] = cbh[i]; break
            else
                if ii > 1
                    mxht = ht[i]
                    for jj in (ii - 1):-1:1
                        j = yes[jj]
                        if mxht * 1.25f0 > cbh[j]
                            mxht < ht[j] && (mxht = ht[j])
                            if mxht >= crit
                                mincb[irep] = cbh[i]; break
                            end
                        end
                    end
                end
                mincb[irep] > -1f0 && break
            end
        end
    end
    sel = [i for i in 1:30 if mincb[i] != -1f0]
    lmin = Float32[log(mincb[i]) for i in sel]
    p = 1f0 / 30f0
    ptr(flm) = begin
        acc = 0f0
        if flm > 0.0001f0
            mxnt = log(fpow(flm / 0.0775f0, 1.45f0) / 30.5f0)
            for lm in lmin
                z = Float64(lm - mxnt) / 0.25
                pt = _pofl_nprob_upper(z)
                pt < 1e-7 && (pt = 0.0)
                acc = Float32(Float64(acc) + pt * Float64(p))
            end
        end
        acc
    end
    r1 = ptr(flm1); r2 = ptr(flm2)
    rannput!(s.rng, saved)
    return (r1, r2)
end

# The FMMAIN-time record list (FMPROB, HT, ICR in storage order). In a tripling cycle FVS has already run TRIPLE
# (grincr.f:543) before GRADD calls FMMAIN, while jl's non-fire path triples after the FMMAIN seam: rebuild the
# tripled view from the stash exactly as triple_records! lays it out (originals ×0.60, then per original the
# upper ×0.25 / lower ×0.15 copies at nlive+2i−1 / nlive+2i, same HT, crown from crU/crL when the variant sets it).
function _pofl_fmmain_records(t, stash)
    n = t.n
    if stash === nothing
        return (t.tpa[1:n], t.height[1:n], Int[t.crown_pct[i] for i in 1:n])
    end
    nl = stash.nlive
    prb = Vector{Float32}(undef, 3nl); ht = Vector{Float32}(undef, 3nl); icr = Vector{Int}(undef, 3nl)
    hascr = hasproperty(stash, :crU)
    @inbounds for i in 1:nl
        prb[i] = t.tpa[i] * 0.60f0; ht[i] = t.height[i]; icr[i] = t.crown_pct[i]
        u = nl + 2i - 1; l = nl + 2i
        prb[u] = t.tpa[i] * 0.25f0; prb[l] = t.tpa[i] * 0.15f0
        ht[u] = t.height[i]; ht[l] = t.height[i]
        icr[u] = (hascr && stash.crU[i] >= 0) ? stash.crU[i] : t.crown_pct[i]
        icr[l] = (hascr && stash.crL[i] >= 0) ? stash.crL[i] : t.crown_pct[i]
    end
    return (prb, ht, icr)
end

# SN/CS skip FMCFIR (CRBURN=0, OINIT1=OACT1=−1); CR/CS/LS/SN/TT/UT re-select the fuel models under each scenario's
# moisture (fmpofl.f:133-139), every other variant keeps the year's FMCFMD3 selection (moisture-independent).
_pofl_east(v) = v isa Southern || v isa CentralStates
_pofl_reselect(v) = v isa CentralRockies || v isa CentralStates || v isa LakeStates || v isa Southern ||
                    v isa Teton || v isa Utah

"""
    fmpofl_report(s, year; cyclen) -> NamedTuple

One FMPOFL year: the DBSFMPF row values (surface/total flame, fire types, PTORCH, OINIT1(1)/OACT1(1), ACTCBH,
CBD, INT(POKILL·100), INT(POVOLK), PSMOKE·P2T, the severe (SFMOD/SFWT) and current (FMOD/FWT) fuel models with
weights INT(W·100+.5)) and the DBSFMPFC conditions (wind PREWND, INT(POTEMP), 100·MOIS) per scenario.
`ptorch=false` leaves PTORCH to the FMMAIN-time seam (`_pofl_fmptrh` needs FVS's FMMAIN RNG state and record list).
`fire_basis`: a fire burned this FMMAIN year — FMCFMD3 re-selects the models after FMBURN (fmmain.f:189) but on
FMTRET's year-start SMALL/LARGE (fmtret.f:371-387, computed before the fire consumed the down wood).
"""
function fmpofl_report(s::StandState, year::Integer; cyclen::Real = 5, fire_basis::Bool = false, ptorch::Bool = true)
    fs = s.fire
    (fs === nothing || !fs.active) && return nothing
    wmult = fire_wind_reduction(fs.percov)
    cf = canopy_bulk_density(s)
    dpmod = _fueltret_dpmod(s, Int(year))
    env = potfire_env(s.variant)                     # (PREWND(1), POTEMP(1), PREWND(2), POTEMP(2))
    east = _pofl_east(s.variant)
    base_models = nothing
    sc = Vector{Any}(undef, 2)
    for (k, fmois) in enumerate((1, 3))
        pc = fs.params.potf[k]
        mois = pc.mois[1] >= 0f0 ? _moisture_matrix(pc.mois) : fuel_moisture(fmois, s.variant)
        prewnd = pc.wind >= 0f0 ? pc.wind : env[2k - 1]
        potemp = pc.temp >= 0f0 ? pc.temp : env[2k]
        burnseas = pc.season > 0 ? Int(pc.season) : 1
        psburn = pc.pab >= 0f0 ? pc.pab : 100f0
        swind = unsafe_trunc(Int, prewnd)                                             # SWIND = INT(PREWND)
        fwind = Float32(swind) * wmult
        models = (_pofl_reselect(s.variant) || base_models === nothing) ? select_fuel_models(s, mois; fire_basis = fire_basis) : base_models
        base_models === nothing && (base_models = models)
        f = _pofl_fmfint(s, models, mois, fwind, dpmod)
        surf = f.flame; pflam = f.flame
        by = f.byram / 60f0
        sch = (63f0 / (140f0 - potemp)) * (fpow(by, 7f0 / 6f0) / fpow(by + fpow(fwind, 3f0), 0.5f0))
        cfir = east ? (; cftype = "SURFACE", crburn = 0f0, oinit = -1f0, oact = -1f0, rfinal = f.sfrate) :
                      _pofl_fmcfir(s, models, mois, f, swind, wmult, cf.cbd, cf.actcbh, dpmod)
        if cfir.crburn > 0f0                                                           # fmpofl.f:188-196
            finten = (f.hpa + cf.tcload * 7744.8f0 * cfir.crburn) * cfir.rfinal / 60f0
            flb = 0.45f0 * fpow(finten, 0.46f0)
            flt = 0.2f0 * fpow(finten, 0.667f0)
            pflam = flb + cfir.crburn * (flt - flb)
            sch = (63f0 / (140f0 - potemp)) * (fpow(finten, 7f0 / 6f0) / fpow(finten + fpow(fwind, 3f0), 0.5f0))
        end
        pomort, povolk, pbrncr = _pofl_fmeff(s, pflam, sch, cfir.crburn, burnseas, psburn, year, cyclen)
        smoke = _pofl_smoke(s, mois, psburn, pbrncr)
        sc[k] = (; surf, pflam, cftype = cfir.cftype, oinit = cfir.oinit, oact = cfir.oact, pomort, povolk, smoke,
                 models = collect(models), prewnd, potemp, mois)
    end
    pt = ptorch ? _pofl_fmptrh(s, year, sc[1].surf, sc[2].surf) : (0f0, 0f0)
    mw(m) = (ntuple(i -> i <= length(m) ? Int(m[i][1]) : 0, 4),
             ntuple(i -> i <= length(m) ? Float64(unsafe_trunc(Int, m[i][2] * 100f0 + 0.5f0)) : 0.0, 4))
    smod, swt = mw(sc[1].models); fmod, fwt = mw(sc[2].models)
    return (; surf_sev = sc[1].surf, surf_mod = sc[2].surf, tot_sev = sc[1].pflam, tot_mod = sc[2].pflam,
            type_sev = sc[1].cftype, type_mod = sc[2].cftype, ptorch_sev = pt[1], ptorch_mod = pt[2],
            torch_index = sc[1].oinit, crown_index = sc[1].oact, canopy_ht = Int(cf.actcbh), cbd = cf.cbd,
            mort_ba_sev = unsafe_trunc(Int, sc[1].pomort * 100f0), mort_ba_mod = unsafe_trunc(Int, sc[2].pomort * 100f0),
            mort_vol_sev = unsafe_trunc(Int, sc[1].povolk), mort_vol_mod = unsafe_trunc(Int, sc[2].povolk),
            smoke_sev = sc[1].smoke, smoke_mod = sc[2].smoke, smod, swt, fmod, fwt,
            cond = ntuple(k -> (; wind = sc[k].prewnd, temp = unsafe_trunc(Int, sc[k].potemp),
                                mois = (100f0 * sc[k].mois[1, 1], 100f0 * sc[k].mois[1, 2], 100f0 * sc[k].mois[1, 3],
                                        100f0 * sc[k].mois[1, 4], 100f0 * sc[k].mois[1, 5], 100f0 * sc[k].mois[2, 1],
                                        100f0 * sc[k].mois[2, 2])), 2))
end
