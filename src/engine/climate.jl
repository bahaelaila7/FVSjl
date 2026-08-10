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
