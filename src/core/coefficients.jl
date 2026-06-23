# =============================================================================
# coefficients.jl (core) — variant-agnostic coefficient container + CSV loader
#
# All species coefficients are loaded from human-readable CSVs (one directory per
# variant) instead of being baked into source — compact code, reviewable data, and
# a new variant is just a new data directory (requirement #10). One-time parse at
# startup, cached per variant (see each variant's `coefficients(::Variant)`).
#
# `species[:descriptive_name]` is a per-species column keyed by its CSV header, so a
# new coefficient is just a new CSV column — no struct edits.
# =============================================================================

"""
    CrownWidthEq

One crown-width equation: a `family` (`:bechtold`/`:bragg`/`:ek`/`:smith`) plus its
coefficients. Evaluated by `crown_width` (engine/crown_width.jl). Data-only here.
"""
struct CrownWidthEq
    family::Symbol
    a::Float32; b::Float32; c::Float32     # intercept, D coef, D²(or Smith D_cm²) coef
    cr_coef::Float32; hi_coef::Float32     # crown-ratio and Hopkins-index coefs (bechtold)
    power::Float32                         # D exponent (bragg/ek)
    dbh_cap::Float32                       # D cap before the quadratic (bechtold)
    max_cw::Float32                        # upper clamp (0 ⇒ none)
end

"""
    SpeciesCoefficients

Loaded reference data for one variant: every numeric per-species column (keyed by
its descriptive CSV header) plus the structural lookup tables. Immutable, shared
across stands/threads.
"""
struct SpeciesCoefficients
    species::Dict{Symbol,Vector{Float32}}        # numeric per-species columns
    code_alpha::Vector{String}
    code_fia::Vector{String}
    code_plants::Vector{String}
    translation::Vector{NTuple{7,String}}        # alpha,fia,plants,cs,ls,ne,sn crosswalk
    site_species_index::Vector{Int32}            # ISNSIS
    site_master_group::Vector{Int32}             # ISNGRP
    master_group_rep::Vector{Int32}              # MGSISP
    forest_location::Dict{Int,NTuple{3,Float32}} # forest → (lat, long, elev)
    valid_habitat::Vector{Int32}
    crown_species::Dict{String,Tuple{String,String}} # 2-char code → (eqn_default, eqn_open)
    crown_eqs::Dict{String,CrownWidthEq}             # equation number → coefficients
    fia_group::Dict{Int,Int}                         # FIA code → ITG stocking group (TAB3 col1)
    fia_stock_eq::Dict{Int,Int}                      # FIA code → stocking-equation index (TAB3 col2)
    stock_b0::Vector{Float32}                        # 36 stocking-equation intercepts (TAB2)
    stock_b1::Vector{Float32}                        # 36 stocking-equation exponents  (TAB2)
    forest_type_codes::Vector{Int32}                 # valid FIA forest-type codes (FTYPE)
    ffe_fuel_dead::Matrix{Float32}                   # FFE dead surface fuels [9 forest types × 11 size classes] (FUINI)
    ffe_fuel_live::Matrix{Float32}                   # FFE live herb/shrub fuels [4 live types × 2] (FULIV)
    ffe_fuel_models::Matrix{Float32}                 # standard fire-behavior fuel models [13 models × params] (fminit.f)
end

"Per-species coefficient column by its descriptive name."
@inline coef_col(c::SpeciesCoefficients, name::Symbol) = c.species[name]

function _read_csv(path::AbstractString)
    lines = readlines(path)
    header = Symbol.(strip.(split(lines[1], ',')))
    rows = [String.(strip.(split(l, ','))) for l in @view(lines[2:end]) if !isempty(strip(l))]
    return header, rows
end

"""
    load_species_coefficients(datadir) -> SpeciesCoefficients

Read a variant's CSV files from `datadir`. Every non-code column of
`species_coefficients.csv` becomes a per-species Float32 vector keyed by its header.
"""
function load_species_coefficients(datadir::AbstractString)
    hdr, rows = _read_csv(joinpath(datadir, "species_coefficients.csv"))
    nsp = length(rows)
    ci = Dict(h => i for (i, h) in enumerate(hdr))
    alpha  = [rpad(rows[i][ci[:code_alpha]], 4)  for i in 1:nsp]
    fia    = [rows[i][ci[:code_fia]]              for i in 1:nsp]
    plants = [rpad(rows[i][ci[:code_plants]], 6) for i in 1:nsp]
    codecols = (:species_index, :code_alpha, :code_fia, :code_plants)
    species = Dict{Symbol,Vector{Float32}}()
    for (j, h) in enumerate(hdr)
        h in codecols && continue
        species[h] = Float32[parse(Float32, rows[i][j]) for i in 1:nsp]
    end

    _, trows = _read_csv(joinpath(datadir, "species_translation.csv"))
    translation = NTuple{7,String}[(r[1], r[2], r[3], r[4], r[5], r[6], r[7]) for r in trows]

    _, ss = _read_csv(joinpath(datadir, "site_species.csv"))
    site_idx = Int32[parse(Int32, r[1]) for r in ss]
    site_grp = Int32[parse(Int32, r[2]) for r in ss]
    _, mg = _read_csv(joinpath(datadir, "site_master_group.csv"))
    grp_rep = Int32[parse(Int32, r[2]) for r in mg]

    _, fl = _read_csv(joinpath(datadir, "forest_locations.csv"))
    forests = Dict{Int,NTuple{3,Float32}}(
        parse(Int, r[1]) => (parse(Float32, r[2]), parse(Float32, r[3]), parse(Float32, r[4]))
        for r in fl)

    _, hc = _read_csv(joinpath(datadir, "valid_habitat_codes.csv"))
    habitat = zeros(Int32, 122)
    for (k, r) in enumerate(hc); k <= 122 && (habitat[k] = parse(Int32, r[1])); end

    _, cws = _read_csv(joinpath(datadir, "crown_width_species.csv"))
    crown_species = Dict{String,Tuple{String,String}}(r[1] => (r[2], r[3]) for r in cws)
    ceh, cwe = _read_csv(joinpath(datadir, "crown_width_equations.csv"))
    cidx = Dict(h => i for (i, h) in enumerate(ceh))
    g(r, name) = parse(Float32, r[cidx[name]])
    crown_eqs = Dict{String,CrownWidthEq}(
        r[cidx[:equation]] => CrownWidthEq(Symbol(r[cidx[:family]]),
            g(r, :a), g(r, :b), g(r, :c), g(r, :cr_coef), g(r, :hi_coef),
            g(r, :power), g(r, :dbh_cap), g(r, :max_cw))
        for r in cwe)

    _, sc = _read_csv(joinpath(datadir, "stocking_coeffs.csv"))
    stock_b0 = zeros(Float32, 36); stock_b1 = zeros(Float32, 36)
    for r in sc; e = parse(Int, r[1]); stock_b0[e] = parse(Float32, r[2]); stock_b1[e] = parse(Float32, r[3]); end
    _, fm = _read_csv(joinpath(datadir, "fia_stocking_map.csv"))
    fia_group = Dict{Int,Int}(parse(Int, r[1]) => parse(Int, r[2]) for r in fm)
    fia_stock_eq = Dict{Int,Int}(parse(Int, r[1]) => parse(Int, r[3]) for r in fm)
    _, ftc = _read_csv(joinpath(datadir, "forest_type_codes.csv"))
    forest_type_codes = Int32[parse(Int32, r[1]) for r in ftc]

    # FFE surface fuel-loading tables (FMCBA): dead fuels by forest type × size class,
    # live herb/shrub by live-fuel type. Numeric columns only (skip the id + name columns).
    hd, dr = _read_csv(joinpath(datadir, "fire_fuel_dead.csv"))
    ffe_fuel_dead = Float32[parse(Float32, dr[i][j]) for i in 1:length(dr), j in 3:length(hd)]
    hl, lr = _read_csv(joinpath(datadir, "fire_fuel_live.csv"))
    ffe_fuel_live = Float32[parse(Float32, lr[i][j]) for i in 1:length(lr), j in 3:length(hl)]
    # standard fire-behavior fuel models (Anderson 13): numeric columns (skip model id + name)
    hm, mr = _read_csv(joinpath(datadir, "fire_fuel_models.csv"))
    ffe_fuel_models = Float32[parse(Float32, mr[i][j]) for i in 1:length(mr), j in 3:length(hm)]

    # merchantability specs (MRULES) + HTDBH height-dub coeffs → per-species columns
    for fname in ("merch_specs.csv", "htdbh_coeffs.csv", "crown_ratio_coeffs.csv",
                  "sprout_essprt.csv", "sprout_htdbh_wykoff.csv", "fire_biomass.csv",
                  "fire_species_props.csv")
        h2, rows2 = _read_csv(joinpath(datadir, fname))
        for (j, h) in enumerate(h2)
            h == :species_index && continue
            species[h] = Float32[parse(Float32, rows2[i][j]) for i in 1:length(rows2)]
        end
    end

    return SpeciesCoefficients(species, alpha, fia, plants, translation,
                               site_idx, site_grp, grp_rep, forests, habitat,
                               crown_species, crown_eqs,
                               fia_group, fia_stock_eq, stock_b0, stock_b1, forest_type_codes,
                               ffe_fuel_dead, ffe_fuel_live, ffe_fuel_models)
end

"Cached coefficient load, one per variant (filled by each variant's `coefficients`)."
const _COEF_CACHE = Dict{String,SpeciesCoefficients}()
