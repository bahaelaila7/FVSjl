# =============================================================================
# climate.jl — Climate-FVS extension (clinit/clin/clgmult/clmorts/clmaxden/clputget/clauestb).
#
# CHUNK-1a (LANDED, validated): the CLIMDATA reader — `parse_climdata` ports clin.f:108-239's
# inline table parse; `resolve_climate_indices` ports clin.f:210-227's attribute-column lookup.
# Validated vs live FVSie_clean tests/FVSie/climate.key (CGCM3_A2): 81 attributes × 4 years
# (1990/2030/2060/2090), junk/wrong-scenario records skipped, matching the live echo
# "NUMBER OF ATTRIBUTES=81; NUMBER OF YEARS=4". Oracle staged at /workspace/.iework/climate/.
#
# NOT YET WIRED into keyword dispatch or growth (chunk-1b = clgmult VSCORE/PS viability +
# TREEMULT=1+(PS-1)·CLGROWMULT applied to DDS; see docs/WESTERN_ROLLOUT_RECONCILIATION_2026-08-07.md).
# This file is `include`d but nothing calls it yet ⇒ zero effect on existing runs.
# =============================================================================

const MXCLYEARS = 4       # CLIMATE.F77 PARAMETER (MXCLYEARS=4) — max climate-data years (rows of ATTRS)
const MXCLATTRS = 130     # CLIMATE.F77 PARAMETER (MXCLATTRS=130) — max attributes (columns of ATTRS)

"Parsed CLIMDATA for one stand+scenario: attribute labels, the years present, and the
year×attribute value matrix (clin.f COMMON /CLIMATE/ ATTRS(MXCLYEARS,MXCLATTRS) + YEARS)."
struct ClimateData
    labels::Vector{String}          # ATTR_LABELS (dd5/mat/mtcm/… + per-species viability + DE*)
    years::Vector{Int}              # YEARS(1:NYEARS)
    attrs::Matrix{Float32}          # ATTRS[year, attr]
end

"""
Per-stand Climate-FVS state (held by `StandState.climate`, `nothing` until a CLIMATE keyword fires).
`active` = LCLIMATE (NATTRS>0 && NYEARS>0). `growmult`/`mortmult` are the per-species CLGROWMULT /
CLMRTMLT1 keyword weights (GrowMult/MortMult), default 1. `plant_symbols` maps species index → PLANTS
symbol for viability-column lookup. `indices` caches the named-attribute columns (resolve_climate_indices).
"""
# Per-variant PLANTS symbols (PLNJSP) — species index → USDA PLANTS code, used to locate each
# species' climate-viability column and drive the clgmult XRELGR dispatch. IE: blkdat.f:186.
const _IE_PLNJSP = String[
    "PIMO3","LAOC","PSME","ABGR","TSHE","THPL","PICO","PIEN","ABLA","PIPO",
    "TSME","PIAL","PIFL2","LALY","PIMO","JUSC2","TABR2","POTR5","POPUL","ACGL","BEPA","2TB","2TN"]

"PLANTS symbols per species for the Climate-FVS viability lookup (empty ⇒ variant not climate-wired)."
climate_plant_symbols(::AbstractVariant) = String[]

mutable struct ClimateState <: AbstractClimateState
    active::Bool
    data::ClimateData
    indices::Dict{Symbol,Int}
    plant_symbols::Vector{String}
    growmult::Vector{Float32}
    mortmult::Vector{Float32}
    inv_year::Int
    # GrowMult/MortMult keyword weights, cycle-SCHEDULED (clin.f opt5/opt3 → OPADD; clgmult.f:44/clmorts.f:41
    # OPGET applies the value when ICYC reaches the scheduled cycle, then it persists). Each event =
    # (cycle, sp, value); sp=0 ⇒ all species. Applied per-cycle by `apply_climate_schedule!`.
    grow_events::Vector{Tuple{Int,Int,Float32}}
    mort_events::Vector{Tuple{Int,Int,Float32}}
    # AutoEstb (climate auto-establishment, clauestb.f) events = (cycle, aestock%, aesntrees, nespecies).
    # A RECURRING activity: once icyc ≥ its cycle it fires EVERY cycle (OPINCR), scheduling NATURAL regen.
    # The latest event with cycle ≤ icyc supplies the active params.
    autoestb::Vector{Tuple{Int,Float32,Float32,Int}}
end

"""
    parse_climdata(lines, nplt, climname) -> ClimateData

Port of the CLIMDATA inline read (clin.f:108-239). `lines` are the records after the
`ClimData`/scenario/`*` header (i.e. starting at the column-label line). Line 1 = header
`Stand_ID,Scenario,Year,<labels…>` (first 3 dropped, trailing blanks trimmed → NATTRS).
Subsequent rows `Stand_ID,Scenario,Year,<vals…>` are kept only when `Stand_ID==nplt` AND
`Scenario==climname` (clin.f:142-143 skip otherwise — e.g. JUNK/other-stand rows). A new year
beyond `MXCLYEARS` raises the "TOO MANY YEARS" error (clin.f:157). `-999` terminates.
"""
function parse_climdata(lines::AbstractVector{<:AbstractString}, nplt::AbstractString, climname::AbstractString)::ClimateData
    isempty(lines) && error("CLIMDATA: empty block")
    hdr = split(strip(lines[1]), ',')
    length(hdr) < 4 && error("CLIMDATA: header has no attribute labels")
    labels = String[strip(x) for x in hdr[4:end]]
    while !isempty(labels) && isempty(last(labels)); pop!(labels); end   # clin.f:132-136
    nattrs = length(labels)
    years = Int[]; rows = Vector{Float32}[]
    @inbounds for i in 2:length(lines)
        s = strip(lines[i])
        (isempty(s) || startswith(s, '*')) && continue
        s == "-999" && break
        f = split(s, ',')
        length(f) < 3 + nattrs && continue
        (strip(f[1]) == nplt && strip(f[2]) == climname) || continue     # clin.f:142-143
        yr = parse(Int, strip(f[3]))
        yr in years && continue
        length(years) >= MXCLYEARS && error("TOO MANY YEARS IN CLIMATE DATA (>$MXCLYEARS)")  # clin.f:157
        push!(years, yr)
        push!(rows, Float32[parse(Float32, strip(f[3 + k])) for k in 1:nattrs])
    end
    attrs = isempty(rows) ? zeros(Float32, 0, nattrs) : reduce(vcat, (permutedims(r) for r in rows))
    return ClimateData(labels, years, attrs)
end

"""
    algslp(xx, x, y) -> Float32

Port of `algslp.f` — piecewise-linear interpolation of the series `(x, y)` at `xx`, with FLAT
extrapolation: `xx < x[1] → y[1]`; `xx ≥ x[end] → y[end]`; else linearly interpolate between the
bracketing knots. `x` must be ascending (the climate `YEARS`). Used by clgmult/clmorts to sample
a climate attribute's time-series at the inventory / current / birth year.
"""
function algslp(xx::Real, x::AbstractVector{<:Real}, y::AbstractVector{<:Real})::Float32
    n = length(x)
    x1 = Float32(xx)
    x1 < x[1] && return Float32(y[1])
    x1 >= x[n] && return Float32(y[n])
    @inbounds for i in 1:(n - 1)
        if x1 < x[i + 1]
            return Float32(y[i] + ((y[i + 1] - y[i]) / (x[i + 1] - x[i])) * (x1 - x[i]))
        end
    end
    return Float32(y[n])
end

"""
    vscore_transform(spviab) -> Float32

Port of clgmult.f:118-123 — turn a raw species viability score (0..1, the interpolated
per-species climate-viability attribute) into the growth-viability score VSCORE:
`>0.5 → 1.0`; else `-0.66666667 + 3.3333333·spviab`; floored at 0.2.
"""
@inline function vscore_transform(spviab::Real)::Float32
    s = Float32(spviab)
    v = s > 0.5f0 ? 1.0f0 : (-0.66666667f0 + s * 3.3333333f0)
    return v < 0.2f0 ? 0.2f0 : v
end

"""
    species_vscore(cd, plant_symbol, thisyr) -> (spviab, vscore)

Port of clgmult.f:112-127 for one species: locate the species' viability column in `cd`
by its PLANTS symbol (INDXSPECIES = the label matching `plant_symbol`, clin.f:241-248),
interpolate it to `thisyr` (= IY(ICYC)+FINT/2), and apply [`vscore_transform`](@ref).
Returns `(1.0, 1.0)` when the species has no viability column (clgmult leaves VSCORE=1).
"""
function species_vscore(cd::ClimateData, plant_symbol::AbstractString, thisyr::Real)
    col = findfirst(==(plant_symbol), cd.labels)
    col === nothing && return (1.0f0, 1.0f0)
    spviab = algslp(thisyr, cd.years, view(cd.attrs, :, col))
    return (spviab, vscore_transform(spviab))
end

# --- Leites et al. (Ecol. Appl. 22(1):154-165) transfer-distance relative-growth models ---
# clgmult.f:161-197. Each returns GROW_TD1/GROW_TD0 with both terms floored at 0.05 before the ratio.
# TD = "now" minus "birth" climate (time substituted for space per the FVS note).

"XDF — Douglas-fir/lodgepole growth ratio from MTCM transfer distance (clgmult.f:161-167)."
@inline function leites_xdf(mtcm_now::Real, mtcm_birth::Real)::Float32
    td = Float32(mtcm_now) - Float32(mtcm_birth); mb = Float32(mtcm_birth)
    td0 = 373.97f0 + 38.52f0 * mb
    td1 = 373.97f0 + 6.799f0 * td - 3.726f0 * td * td + 38.52f0 * mb - 3.602f0 * mb * td
    max(td1, 0.05f0) / max(td0, 0.05f0)
end

"XWL — larch/spruce/hemlock growth ratio from MMIN transfer distance + DD0 (clgmult.f:176-182)."
@inline function leites_xwl(mmin_now::Real, mmin_birth::Real, dd0_birth::Real)::Float32
    td = Float32(mmin_now) - Float32(mmin_birth); db = Float32(dd0_birth)
    td0 = 542.20f0 - 0.1468f0 * db
    td1 = 542.20f0 + 17.50f0 * td - 1.215f0 * td * td - 0.1468f0 * db - 0.0187f0 * td * db
    max(td1, 0.05f0) / max(td0, 0.05f0)
end

"XPP — ponderosa growth ratio from SMI(=dd5/gsp) transfer distance + D100 (clgmult.f:191-197)."
@inline function leites_xpp(smi_now::Real, smi_birth::Real, d100_birth::Real)::Float32
    td = Float32(smi_now) - Float32(smi_birth); db = Float32(d100_birth)
    td0 = 551.20221f0 - 2.02135f0 * db
    td1 = 551.20221f0 - 14.88483f0 * td - 0.58027f0 * td * td - 2.02135f0 * db + 0.15582f0 * td * db
    max(td1, 0.05f0) / max(td0, 0.05f0)
end

"""
    clim_xgsite(psite_now, psite_invyr) -> Float32

Site-productivity growth adjustment from the pSite attribute (clgmult.f:90-106): 2.0 if the
inventory-year pSite ≤ .001; 1.0 if unchanged; else 1.5819767·(1−exp(−pSite_now/pSite_inv)),
snapped to 1.0 within .01.
"""
@inline function clim_xgsite(psite_now::Real, psite_invyr::Real)::Float32
    pi = Float32(psite_invyr); pn = Float32(psite_now)
    pi <= 0.001f0 && return 2.0f0
    abs(pn - pi) < 0.001f0 && return 1.0f0
    x = 1.5819767f0 * (1.0f0 - exp(-(pn / pi)))
    abs(x - 1.0f0) < 0.01f0 ? 1.0f0 : x
end

"""
    clim_xrelgr(plant_symbol, xdf, xpp, xwl) -> Float32

Per-species relative-growth dispatch (clgmult.f:199-217): PSME/PICO→XDF, PIPO→XPP,
LAOC/PIMO3/PIEN/TSHE→XWL, else the half-weighted mean of the three; snapped to 1.0 within .005.
"""
@inline function clim_xrelgr(sym::AbstractString, xdf::Real, xpp::Real, xwl::Real)::Float32
    x = (sym == "PSME" || sym == "PICO") ? Float32(xdf) :
        sym == "PIPO"                    ? Float32(xpp) :
        (sym == "LAOC" || sym == "PIMO3" || sym == "PIEN" || sym == "TSHE") ? Float32(xwl) :
        1.0f0 + (((Float32(xdf) + Float32(xpp) + Float32(xwl)) / 3.0f0 - 1.0f0) * 0.5f0)
    abs(x - 1.0f0) < 0.005f0 ? 1.0f0 : x
end

"""
    clim_treemult(xgsite, xrelgr, vscore, clgrowmult) -> (ps, treemult)

Combine into the per-tree growth multiplier (clgmult.f:222-226): PS = min(xgsite, xrelgr, vscore)
unless PS>0.99 (then max of the three), capped at 3; TREEMULT = 1+(PS−1)·CLGROWMULT, floored at 0.
`clgrowmult` is the species' GrowMult keyword weight (default 1 when GrowMult sets ALL/species).
"""
@inline function clim_treemult(xgsite::Real, xrelgr::Real, vscore::Real, clgrowmult::Real)
    ps = min(Float32(xgsite), Float32(xrelgr), Float32(vscore))
    ps > 0.99f0 && (ps = max(Float32(xgsite), Float32(xrelgr), Float32(vscore)))
    ps > 3.0f0 && (ps = 3.0f0)
    tm = 1.0f0 + (ps - 1.0f0) * Float32(clgrowmult)
    return ps, (tm < 0.0f0 ? 0.0f0 : tm)
end

"""
    apply_climate_dds!(s, wk2, thisyr)

Apply the Climate-FVS growth multiplier to the large-tree DDS (`wk2` = ln(DDS)), mirroring
dgdriv.f:217 `DDS = EXP(WK2)·WK4` where WK4 = clgmult's per-tree TREEMULT — so here `wk2[i] +=
log(treemult[i])`. Per cycle: XGSITE + per-species VSCORE at `thisyr`; per tree: BIRTHYR = thisyr −
birth_age, the Leites XDF/XWL/XPP, XRELGR (by PLNJSP), and `clim_treemult`. No-op when climate is
inactive or the required attribute columns are absent.
"""
function apply_climate_dds!(s::StandState, wk2::AbstractVector{Float32}, thisyr::Real)
    c = s.climate
    (c === nothing || !c.active) && return s
    cd = c.data; ix = c.indices; t = s.trees
    (ix[:mtcm] == 0 || ix[:mmin] == 0 || ix[:dd0] == 0 || ix[:d100] == 0 ||
     ix[:dd5] == 0 || ix[:gsp] == 0) && return s
    ty = Float32(thisyr)
    A(sym, yr) = algslp(yr, cd.years, view(cd.attrs, :, ix[sym]))
    smi(yr) = (g = A(:gsp, yr); g > 0f0 ? A(:dd5, yr) / g : 0f0)
    xgsite = ix[:pSite] > 0 ? clim_xgsite(A(:pSite, ty), A(:pSite, Float32(c.inv_year))) : 1f0
    mtcm_now = A(:mtcm, ty); mmin_now = A(:mmin, ty); smi_now = smi(ty)
    ns = length(c.plant_symbols)
    vscore = ones(Float32, ns)
    @inbounds for sp in 1:ns
        _, vscore[sp] = species_vscore(cd, c.plant_symbols[sp], ty)
    end
    @inbounds for i in 1:t.n
        t.dbh[i] <= 0f0 && continue
        sp = Int(t.species[i]); (sp < 1 || sp > ns) && continue
        birthyr = ty - t.birth_age[i]
        xdf = leites_xdf(mtcm_now, A(:mtcm, birthyr))
        xwl = leites_xwl(mmin_now, A(:mmin, birthyr), A(:dd0, birthyr))
        xpp = leites_xpp(smi_now, smi(birthyr), A(:d100, birthyr))
        xr  = clim_xrelgr(c.plant_symbols[sp], xdf, xpp, xwl)
        _, tm = clim_treemult(xgsite, xr, vscore[sp], c.growmult[sp])
        tm > 0f0 && (wk2[i] += log(tm))
    end
    return s
end

# --- Climate mortality (clmorts.f) — viability → survival → mortality rate ---
# VALIDATED 8/8 IE species vs live FVSie_g16 clmorts debug (cyc1: WH XV.0585→X0→MORT1; RC .489→.9633→.0367;
# LP mult=.5). This is the base viability-mortality path (SPMORT1/FYRMORT). NOT YET ported: the SPCALIB first-
# cycle presence-calibration branch (clmorts.f:92-98) + the SPMORT2 transfer-distance DMORT (clmorts.f:128-223,
# uses the DE* climate-distance attributes) — chunk-1c.2.
const _CLM_VS = Float32[0.2f0, 0.5f0]   # clmorts.f DATA VS/.2,.5/  (viability knots)
const _CLM_SR = Float32[0.0f0, 1.0f0]   # clmorts.f DATA SR/0.,1./  (survival-rate knots)

"""
    clim_survival(xv) -> Float32

10-yr climate survival from a species' raw viability score `xv` (clmorts.f:93): interpolate `xv`
against the VS→SR curve, i.e. `xv<0.2 → 0`, `xv≥0.5 → 1`, linear between. (The SPCALIB-calibrated
branch, clmorts.f:95-97, is a first-cycle refinement — chunk C, needs the presence-calibration state.)
"""
@inline clim_survival(xv::Real)::Float32 = algslp(xv, _CLM_VS, _CLM_SR)

"""
    clim_mort_rates(x, fint, mult) -> (spmort1, fyrmort)

Climate mortality from 10-yr survival `x` (clmorts.f:103-120): `spmort1 = (1−x)·mult` (10-yr, for
reporting); `fyrmort = (1 − x^(fint/10))·mult` (the applied fint-yr rate; `x→0` when `x≤1e-5`).
`mult` = the species' MortMult keyword weight `CLMRTMLT1` (default 1).
"""
function clim_mort_rates(x::Real, fint::Real, mult::Real)
    x1 = Float32(x); m = Float32(mult)
    spmort1 = (1f0 - x1) * m
    xf = x1 > 1f-5 ? clamp(exp(log(x1) / 10f0)^Float32(fint), 0f0, 1f0) : 0f0
    return spmort1, (1f0 - xf) * m
end

"""
    apply_climate_mort!(s, killed, thisyr, fint)

Apply Climate-FVS mortality (clmorts.f:240-259): per tree the applied mortality rate = MAX(base FVS rate,
climate rate), where the climate rate is the viability-driven `fyrmort` (clim_survival→clim_mort_rates).
If the climate rate exceeds the base, `killed[i] = tpa[i]·fyrmort[sp]` (clmorts.f:259 WK2=PROB·DMORT).
No-op when climate is inactive. (The SPMORT2 transfer-distance DMORT path, clmorts.f:145-237, needs the
DE* attributes + LDMORT gate — chunk-M2; the viability path drives the warming-scenario stand die-off.)
"""
function apply_climate_mort!(s::StandState, killed::AbstractVector{Float32}, thisyr::Real, fint::Real)
    c = s.climate
    (c === nothing || !c.active) && return s
    cd = c.data; ix = c.indices; t = s.trees; ns = length(c.plant_symbols)
    ty = Float32(thisyr); fi = Float32(fint)
    # (1) Per-species viability FYRMORT (clmorts.f:79-126) — the SPMORT1/FYRMORT loop.
    fy = zeros(Float32, ns)
    @inbounds for sp in 1:ns
        xv, _ = species_vscore(cd, c.plant_symbols[sp], ty)   # raw viability at THISYR
        _, fy[sp] = clim_mort_rates(clim_survival(xv), fi, c.mortmult[sp])
    end
    # (2) Climate-transfer-distance DMORT (clmorts.f:133-237). LDMORT gate = all DE* + the
    # climate columns present (clmorts.f:133-138). CTHISYR = current-year climate metrics.
    A(sym, yr) = algslp(yr, cd.years, view(cd.attrs, :, ix[sym]))
    de(sym) = (j = ix[sym]; j > 0 ? cd.attrs[1, j] : 0f0)     # DE threshold = row-1 value (constant)
    ldmort = ix[:DEmtwm] > 0 && ix[:DEmtcm] > 0 && ix[:DEdd5] > 0 && ix[:DEsdi] > 0 &&
             ix[:DEdd0] > 0 && ix[:DEpdd5] > 0 && ix[:mtwm] > 0 && ix[:mtcm] > 0 &&
             ix[:dd5] > 0 && ix[:dd0] > 0 && ix[:map] > 0 && ix[:gsp] > 0 && ix[:gsdd5] > 0
    if ldmort
        ct1 = A(:mtwm, ty); ct2 = A(:mtcm, ty); ct3 = A(:dd5, ty)
        ct4 = sqrt(A(:gsdd5, ty)) / A(:gsp, ty)               # sdi = sqrt(gsdd5)/gsp
        ct5 = A(:dd0, ty); ct6 = A(:map, ty) * ct3 / 1000f0   # mapdd5 = map·dd5/1000
        de1 = de(:DEmtwm); de2 = de(:DEmtcm); de3 = de(:DEdd5)
        de4 = de(:DEsdi); de5 = de(:DEdd0); de6 = de(:DEpdd5)
        clmrtmlt2 = 1f0                                        # CLMRTMLT2 (MortMult 2nd param); default 1
        @inbounds for i in 1:t.n
            pr = t.tpa[i]; pr <= 0f0 && continue
            sp = Int(t.species[i]); (sp < 1 || sp > ns) && continue
            by = ty - t.birth_age[i]                          # BIRTHYR = THISYR - ABIRTH
            # DTV = CTHISYR - CBIRTH, each normalized by its DE* threshold (1.0 if threshold ≤0)
            d1 = de1 > 0f0 ? (ct1 - A(:mtwm, by)) / de1 : 1f0
            d2 = de2 > 0f0 ? (ct2 - A(:mtcm, by)) / de2 : 1f0
            d3 = de3 > 0f0 ? (ct3 - A(:dd5, by)) / de3 : 1f0
            d4 = de4 > 0f0 ? (ct4 - sqrt(A(:gsdd5, by)) / A(:gsp, by)) / de4 : 1f0
            d5 = de5 > 0f0 ? (ct5 - A(:dd0, by)) / de5 : 1f0
            d6 = de6 > 0f0 ? (ct6 - A(:map, by) * A(:dd5, by) / 1000f0) / de6 : 1f0
            dm = (d1 + d2 + d3 + d4 + d5 + d6) / 6f0 - 1.1f0
            dm = clamp(dm, 0f0, 5.9f0)
            dm = 0.9f0 * (1f0 - exp(-dm^2.5f0))               # transfer-distance mortality curve
            surv = 1f0 - dm                                    # → survival, then FINT-yr rate
            surv = surv > 1f-5 ? clamp(exp(log(surv) / 10f0)^fi, 0f0, 1f0) : 0f0
            dmr = (1f0 - surv) * clmrtmlt2                     # back to mortality rate ·CLMRTMLT2
            rate = max(fy[sp], dmr)                            # clmorts.f:258 DMORT = max(FYRMORT, DMORT)
            rate > killed[i] / pr && (killed[i] = pr * rate)   # clmorts.f:259 WK2 = PROB·DMORT
        end
    else
        @inbounds for i in 1:t.n                               # viability-only fallback (no DE* columns)
            pr = t.tpa[i]; pr <= 0f0 && continue
            sp = Int(t.species[i]); (sp < 1 || sp > ns) && continue
            fy[sp] > killed[i] / pr && (killed[i] = pr * fy[sp])
        end
    end
    return s
end

"""
    apply_climate_schedule!(s, icyc)

Realize the cycle-SCHEDULED GrowMult/MortMult keyword weights for cycle `icyc` (FVS OPADD/OPGET: an event
scheduled for cycle N takes effect once ICYC≥N and persists). Resets `growmult`/`mortmult` to 1 then applies
every event with `cycle ≤ icyc` in keyword-file order (later entries override earlier — matches a specific
species overriding a preceding `All`). No-op when climate is inactive or no events were parsed.
"""
function apply_climate_schedule!(s::StandState, icyc::Integer)
    c = s.climate
    (c === nothing || !c.active) && return s
    ns = length(c.growmult)
    if !isempty(c.grow_events)
        fill!(c.growmult, 1f0)
        @inbounds for (cyc, sp, val) in c.grow_events
            cyc <= icyc || continue
            sp == 0 ? fill!(c.growmult, val) : (1 <= sp <= ns && (c.growmult[sp] = val))
        end
    end
    if !isempty(c.mort_events)
        fill!(c.mortmult, 1f0)
        @inbounds for (cyc, sp, val) in c.mort_events
            cyc <= icyc || continue
            sp == 0 ? fill!(c.mortmult, val) : (1 <= sp <= ns && (c.mortmult[sp] = val))
        end
    end
    return s
end

"""
    clim_autoestb!(s, icyc, fint)

Climate auto-establishment (clauestb.f): when an AutoEstb activity is active (icyc ≥ its cycle), schedule
NATURAL(431) regen for next cycle. TMAXTRS = (XMAX/0.02483133)·max(5,QMD)^−1.605 (max trees at current QMD);
PTREES = clamp(2−4·TPROB/TMAXTRS, 0,1) (low stocking ⇒ establish more); per-species viability POTESTAB (top
`nespecies` with viab≥0.4), scaled XX=clamp(−1+2.5·viab,0,1), normalized, trees = PTREES·TTOADD·frac (dropped
if ≤1). Skips if stocking TPROB > TMAXTRS·AESTOCK·0.01. Trees enter via jl's existing `establish!` (params =
[sp, tpa, 100%surv, age0, ht0, shade0], clauestb.f:226-232). No-op unless a CLIMATE block parsed AutoEstb.
"""
function clim_autoestb!(s::StandState, icyc::Integer, fint::Real)
    c = s.climate
    (c === nothing || !c.active || isempty(c.autoestb)) && return s
    active = nothing
    @inbounds for e in c.autoestb; e[1] <= icyc && (active = e); end   # latest event with cycle ≤ icyc
    active === nothing && return s
    aestock = active[2]; aesntrees = active[3]; nespecies = active[4]
    t = s.trees; cd = c.data; ns = length(c.plant_symbols)
    xmax = stand_sdimax(s)                                             # SDICAL XMAX (pre-CLMAXDEN)
    rmsqd = max(5f0, stand_qmd(s))
    tmaxtrs = (xmax / 0.02483133f0) * rmsqd^(-1.605f0)                 # 0.02483133 = 10^-1.605
    tprob = 0f0; @inbounds for i in 1:t.n; tprob += t.tpa[i]; end
    ptrees = tmaxtrs > 1f0 ? clamp(2f0 - 4f0 * (tprob / tmaxtrs), 0f0, 1f0) : 1f0
    (ptrees * aesntrees > 0f0) || return s
    ty = Float32(current_cycle_year(s)) + Float32(fint) / 2f0
    pot = zeros(Float32, ns)
    @inbounds for sp in 1:ns; pot[sp] = species_vscore(cd, c.plant_symbols[sp], ty)[1]; end   # raw viability
    order = sortperm(pot; rev = true)                                 # RDPSRT desc (ties rare among viabilities)
    nspec = 0
    for sp in order; pot[sp] < 0.4f0 && break; nspec += 1; end
    nspec == 0 && return s
    nspec > nespecies && (nspec = round(Int, nespecies))
    top = order[1:nspec]
    @inbounds for sp in top; pot[sp] = clamp(-1f0 + 2.5f0 * pot[sp], 0f0, 1f0); end
    ttoadd = pot[top[1]] > 0.8f0 ? aesntrees : aesntrees * pot[top[1]]
    ssum = 0f0; @inbounds for sp in top; ssum += pot[sp]; end
    ssum > 0.001f0 || return s
    tprob > tmaxtrs * aestock * 0.01f0 && return s                    # stocking gate (clauestb.f:219)
    nextyr = Int(current_cycle_year(s)) + round(Int, fint)
    @inbounds for sp in top
        trees = ptrees * ttoadd * (pot[sp] / ssum)
        trees <= 1f0 && continue
        push!(s.control.schedule, ScheduledActivity(nextyr, Int32(431),
                                                    (Float32(sp), trees, 100f0, 0f0, 0f0, 0f0)))
        s.estab.active = true    # activate the ESTAB packet so establish! processes the scheduled NATURAL regen
    end
    return s
end

"""
    resolve_climate_indices(labels) -> Dict{Symbol,Int}

Port of clin.f:210-227 — locate the column of each named climate attribute used by the
viability/growth models (dd5/mat/map/mtcm/mtwm/gsp/d100/mmin/dd0/gsdd5/pSite and the
species-deviation DEmtwm/DEmtcm/DEdd5/DEsdi/DEdd0/DEpdd5). Returns 0 for any absent label.
"""
function resolve_climate_indices(labels::AbstractVector{<:AbstractString})::Dict{Symbol,Int}
    want = (:dd5, :mat, :map, :mtcm, :mtwm, :gsp, :d100, :mmin, :dd0, :gsdd5, :pSite,
            :DEmtwm, :DEmtcm, :DEdd5, :DEsdi, :DEdd0, :DEpdd5)
    idx = Dict{Symbol,Int}(k => 0 for k in want)
    @inbounds for (i, lab) in enumerate(labels)
        s = Symbol(strip(lab))
        haskey(idx, s) && (idx[s] = i)
    end
    return idx
end
