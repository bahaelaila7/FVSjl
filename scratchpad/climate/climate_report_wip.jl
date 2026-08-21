# WIP climate_report + _climate_potestab (validated: Viability/SiteMult/MxDenMult bit-exact,
# BA/TPA/GrowthMult cornered; mort1/mort2 need SPCALIB presence-calibration, potestab needs AutoEstb source)
"""
    climate_report(s; report_year, fint) -> Vector of per-species NamedTuples

The Climate-FVS Viability-and-Effects report (clauestb.f:178-207 → DBSCLSUM / FVS_Climate table). One entry per
species that passes the importance filter `SPIMP>0.05 || SPVIAB>0.4` (with a viability column). Computed at the
POST-growth stand state, sampling climate at `report_year + fint/2` — the same THISYR = IY(ICYC)+FINT/2 midpoint
clgmult/clmorts/clauestb use (clgmult.f:81). FVS labels the row with the cycle-START year but the data reflects
the post-growth tree list + the mid-cycle climate — so a caller collects this right AFTER `grow_cycle!`, passing
the pre-advance cycle year as `report_year`. Fields (dbsclsum.f column order): sp (1-based), viab (raw SPVIAB),
ba, tpa, mort1 (FYRMORT viability rate), mort2 (transfer-distance DMORT rate), gmult (Σtreemult·prob/Σprob),
sitgm (xgsite^clgrowmult), mxden, potestab. No-op (empty) when climate is inactive.
"""
function climate_report(s::StandState; report_year::Real, fint::Real)
    c = s.climate
    (c === nothing || !c.active) && return NamedTuple[]
    cd = c.data; ix = c.indices; t = s.trees; ns = length(c.plant_symbols)
    ty = Float32(report_year) + Float32(fint) / 2f0
    fi = Float32(fint)
    A(sym, yr) = algslp(yr, cd.years, view(cd.attrs, :, ix[sym]))
    smi(yr) = (g = A(:gsp, yr); g > 0f0 ? A(:dd5, yr) / g : 0f0)
    have_grow = !(ix[:mtcm] == 0 || ix[:mmin] == 0 || ix[:dd0] == 0 || ix[:d100] == 0 ||
                  ix[:dd5] == 0 || ix[:gsp] == 0)
    xgsite = ix[:pSite] > 0 ? clim_xgsite(A(:pSite, ty), A(:pSite, Float32(c.inv_year))) : 1f0
    mtcm_now = have_grow ? A(:mtcm, ty) : 0f0; mmin_now = have_grow ? A(:mmin, ty) : 0f0
    smi_now = have_grow ? smi(ty) : 0f0
    # per-species accumulators (clauestb.f:178-192 SPBA/SPTPA + clgmult.f:243 SPGMULT weighting)
    spba = zeros(Float32, ns); sptpa = zeros(Float32, ns)
    gm_num = zeros(Float32, ns); gm_wt = zeros(Float32, ns)
    spviab = ones(Float32, ns); vscore = ones(Float32, ns)
    @inbounds for sp in 1:ns
        spviab[sp], vscore[sp] = species_vscore(cd, c.plant_symbols[sp], ty)
    end
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        sp = Int(t.species[i]); (sp < 1 || sp > ns) && continue
        pr = t.tpa[i]
        spba[sp]  += d * d * pr * 0.005454154f0
        sptpa[sp] += pr
        if have_grow
            by = ty - t.birth_age[i]
            xdf = leites_xdf(mtcm_now, A(:mtcm, by))
            xwl = leites_xwl(mmin_now, A(:mmin, by), A(:dd0, by))
            xpp = leites_xpp(smi_now, smi(by), A(:d100, by))
            xr  = clim_xrelgr(c.plant_symbols[sp], xdf, xpp, xwl)
            _, tm = clim_treemult(xgsite, xr, vscore[sp], c.growmult[sp])
            gm_num[sp] += tm * pr; gm_wt[sp] += pr
        end
    end
    # per-species FYRMORT viability mortality rate (clmorts.f:79-126 SPMORT1)
    mort1 = zeros(Float32, ns)
    @inbounds for sp in 1:ns
        _, mort1[sp] = clim_mort_rates(clim_survival(spviab[sp]), fi, c.mortmult[sp])
    end
    # MXDENMLT (clmaxden.f:48 — 1 at ICYC≤1) — report uses the same mid-cycle sample. Collected right after
    # grow_cycle! advanced s.control.cycle to the just-grown cycle's index (jl cycle k → FVS ICYC = k+1 = cycle).
    icyc = max(1, Int(s.control.cycle))
    mxden = 1f0
    if icyc > 1
        clmx = 1f0; @inbounds for e in c.mxden; e[1] <= icyc && (clmx = e[2]); end
        mxden = clim_maxden_mult(s, ty, clmx)
    end
    # POTESTAB per species (clauestb.f:80-113): the auto-establishment TPA scored + scaled + normalized.
    potestab = _climate_potestab(s, c, cd, ty, fint, icyc)
    # importance filter SPIMP (clauestb.f:180-192): normalized (BA-share + TPA-share)
    tba = sum(spba); ttpa = sum(sptpa)
    out = NamedTuple[]
    spimp = zeros(Float32, ns)
    @inbounds for sp in 1:ns
        b = tba > 0f0 ? spba[sp] / tba : 0f0
        p = ttpa > 0f0 ? sptpa[sp] / ttpa : 0f0
        spimp[sp] = b + p
    end
    si = sum(spimp); si > 0f0 && (spimp ./= si)
    @inbounds for sp in 1:ns
        (spimp[sp] > 0.05f0 || spviab[sp] > 0.4f0) && (findfirst(==(c.plant_symbols[sp]), cd.labels) !== nothing) || continue
        gm = gm_wt[sp] > 0f0 ? gm_num[sp] / gm_wt[sp] : 1f0
        push!(out, (sp = sp, viab = spviab[sp], ba = spba[sp], tpa = sptpa[sp],
                    mort1 = mort1[sp], mort2 = mort1[sp], gmult = gm,
                    sitgm = xgsite^c.growmult[sp], mxden = mxden, potestab = potestab[sp]))
    end
    return out
end

# POTESTAB(I) per species (clauestb.f:80-113): raw viability → top-nespecies scaling → normalize → PTREES·TTOADD·frac.
# Mirrors clim_autoestb! but RETURNS the per-species TPA instead of scheduling. Empty AutoEstb ⇒ all zero.
function _climate_potestab(s::StandState, c, cd, ty::Real, fint::Real, icyc::Integer)::Vector{Float32}
    ns = length(c.plant_symbols); pot = zeros(Float32, ns)
    isempty(c.autoestb) && return pot
    active = nothing
    @inbounds for e in c.autoestb; e[1] <= icyc && (active = e); end
    active === nothing && return pot
    aestock = active[2]; aesntrees = active[3]; nespecies = active[4]
    t = s.trees
    xmax = stand_sdimax(s)
    if icyc > 1
        clmx = 1f0; @inbounds for e in c.mxden; e[1] <= icyc && (clmx = e[2]); end
        xmax *= clim_maxden_mult(s, Float32(ty), clmx)
    end
    rmsqd = max(5f0, stand_qmd(s))
    tmaxtrs = (xmax / 0.02483133f0) * rmsqd^(-1.605f0)
    tprob = 0f0; @inbounds for i in 1:t.n; tprob += t.tpa[i]; end
    ptrees = tmaxtrs > 1f0 ? clamp(2f0 - 4f0 * (tprob / tmaxtrs), 0f0, 1f0) : 1f0
    (ptrees * aesntrees > 0f0) || return pot
    sc = zeros(Float32, ns)
    @inbounds for sp in 1:ns; sc[sp] = species_vscore(cd, c.plant_symbols[sp], Float32(ty))[1]; end
    order = sortperm(sc; rev = true)
    nspec = 0; for sp in order; sc[sp] < 0.4f0 && break; nspec += 1; end
    nspec == 0 && return pot
    nspec > nespecies && (nspec = round(Int, nespecies))
    top = order[1:nspec]
    @inbounds for sp in top; sc[sp] = clamp(-1f0 + 2.5f0 * sc[sp], 0f0, 1f0); end
    ttoadd = sc[top[1]] > 0.8f0 ? aesntrees : aesntrees * sc[top[1]]
    ssum = 0f0; @inbounds for sp in top; ssum += sc[sp]; end
    ssum > 0.001f0 || return pot
    tprob > tmaxtrs * aestock * 0.01f0 && return pot
    @inbounds for sp in top
        xx = ptrees * ttoadd * (sc[sp] / ssum)
        xx <= 1f0 && (xx = 0f0)
        pot[sp] = xx
    end
    return pot
end

