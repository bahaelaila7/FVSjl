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
