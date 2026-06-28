# =============================================================================
# init.jl — stand initialization driver (INITRE entry)
#
# Ported from: base/fvs.f (the INITRE call) + base/initre.f setup.
#
# Brings a fresh `StandState` up to the "ready to simulate" point for ONE stand:
# apply variant BLOCK DATA defaults, then process the keyword stream (which loads
# trees, sets plot design, cycles, site, etc.) until PROCESS/STOP/EOF.
#
# A keyword file may hold several stands separated by PROCESS; `initialize!` does
# one stand and returns how it terminated, so a driver can loop for the rest.
# =============================================================================

"Strip a trailing `.key`/`.KEY` (or any extension) to get the run's base path."
function strip_key_ext(keypath::AbstractString)
    lk = findfirst(".k", keypath)
    lk === nothing && (lk = findfirst(".K", keypath))
    return lk === nothing ? keypath : keypath[1:first(lk)-1]
end

# Open a keyword source as a `KeywordReader` + the run's base path. A `.yaml`/`.yml`
# input (positional OR structured form) is parsed to records and rendered to `.key`
# text so the engine reads it through the same KEYRDR path; `.key`/anything else is
# read verbatim. The base path (used to locate the companion `.tre`/`.csv`) is the
# input minus its extension.
function _keyword_reader(keypath::AbstractString)
    ext = lowercase(splitext(keypath)[2])
    if ext == ".yaml" || ext == ".yml"
        text = keyfile_string(read_keywords_yaml(keypath))
        return KeywordReader(IOBuffer(text)), first(splitext(keypath))
    end
    return open(io -> KeywordReader(io), keypath), strip_key_ext(keypath)
end

# Resolve the variant to run a stand file as: an explicit caller argument always wins;
# otherwise a YAML's top-level `variant:` field selects it; otherwise default to Southern
# (a `.key` has no variant — it's the stock-FVS binary choice). Keeps `.key` runs and
# explicit-variant callers behaving exactly as before.
function _resolve_variant(keypath::AbstractString, variant::Union{AbstractVariant,Nothing})
    variant === nothing || return variant
    ext = lowercase(splitext(keypath)[2])
    if ext == ".yaml" || ext == ".yml"
        vc = yaml_variant_code(keypath)
        vc === nothing || return variant_from_code(vc)
    end
    return Southern()
end

"""
    initialize!(state, kr, base_path) -> Symbol

Initialize one stand from an already-open keyword reader. Applies BLOCK DATA
defaults then processes keywords until PROCESS/STOP/EOF (the returned reason).
"""
function initialize!(s::StandState, kr::KeywordReader, base_path::AbstractString;
                     inherited_format::AbstractString = "")
    load_species_coefficients!(s, s.variant)      # BLOCK DATA (species, TREFMT, RNG seed)
    # TREFMT persists across stands in FVS (it lives in COMMON; INITRE never resets it —
    # only the TREEFMT keyword changes it). A 2nd+ stand with no TREEFMT keyword inherits
    # the previous stand's format, so re-applying the BLOCK DATA default here would break
    # it. Restore the inherited format; a TREEFMT keyword in this stand still overrides.
    isempty(inherited_format) || (s.control.tree_format = inherited_format)
    ranseed!(s.rng, false, s.rng.ss)              # INITRE: RANSED(false,...) → reset to seed
    s.plot.gross_space = -1f0                      # GRINIT reset (sn/grinit.f:156)
    @inbounds for i in 1:MAXSP                      # GRINIT size-cap defaults (sn/grinit.f:62)
        s.control.sp_size_cap[i, 1] = 999f0
        s.control.sp_size_cap[i, 2] = 1f0
        s.control.sp_size_cap[i, 3] = 0f0
        s.control.sp_size_cap[i, 4] = 999f0
    end
    reason = process_keywords!(s, kr, base_path)
    finalize_design!(s)                            # INITRE end: PI:=IPTINV, GROSPC
    site_setup!(s, s.variant)                      # SITSET: fan site index to all species (variant-specific)
    return reason
end

"""
    initialize(keypath; variant=Southern(), faithful=true) -> (state, reason)

Convenience: open `keypath`, build a fresh state, and initialize the FIRST stand.
"""
function initialize(keypath::AbstractString; variant::AbstractVariant = Southern(),
                    faithful::Bool = true)
    s = StandState(variant; faithful = faithful)
    kr, base = _keyword_reader(keypath)            # accepts .key OR .yaml/.yml (structured or positional)
    reason = initialize!(s, kr, base)
    return s, reason
end

"""
    each_stand(keypath; variant=Southern(), faithful=true) -> Vector{StandState}

Initialize EVERY stand in a multi-stand keyword file (stands are separated by
`PROCESS` and the run ends at `STOP`/EOF). FVS re-runs INITRE per stand — each stand
gets a fresh `StandState` (so `ITRN` resets) — but the tree-record format (`TREFMT`)
persists in COMMON across stands, so it is carried forward here; a `TREEFMT` keyword
inside a later stand still overrides it. Returns the per-stand initialized states
(cyc0-ready; run `notre!`/`setup_growth!`/`compute_volumes!` to project each).
"""
function each_stand(keypath::AbstractString;
                    variant::Union{AbstractVariant,Nothing} = nothing,
                    faithful::Bool = true)
    variant = _resolve_variant(keypath, variant)   # explicit arg wins; else YAML `variant:`; else SN
    kr, base = _keyword_reader(keypath)
    stands = StandState[]
    fmt = ""                                       # TREFMT carried across stands
    while true
        s = StandState(variant; faithful = faithful)
        reason = initialize!(s, kr, base; inherited_format = fmt)
        fmt = s.control.tree_format                # may have been set by a TREEFMT keyword
        # A bare STOP/EOF after the last stand's PROCESS is the run terminator, not a
        # stand: no STDINFO ran, no trees, no establishment. Don't emit a phantom stand.
        real = s.plot.user_forest_code != 0 || s.trees.n > 0 || s.estab.active
        real && push!(stands, s)
        reason in (:stop, :eof) && break
    end
    return stands
end
