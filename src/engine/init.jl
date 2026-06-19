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

"""
    initialize!(state, kr, base_path) -> Symbol

Initialize one stand from an already-open keyword reader. Applies BLOCK DATA
defaults then processes keywords until PROCESS/STOP/EOF (the returned reason).
"""
function initialize!(s::StandState, kr::KeywordReader, base_path::AbstractString)
    load_species_coefficients!(s, s.variant)      # BLOCK DATA (species, TREFMT, RNG seed)
    ranseed!(s.rng, false, s.rng.ss)              # INITRE: RANSED(false,...) → reset to seed
    s.plot.gross_space = -1f0                      # GRINIT reset (sn/grinit.f:156)
    reason = process_keywords!(s, kr, base_path)
    finalize_design!(s)                            # INITRE end: PI:=IPTINV, GROSPC
    return reason
end

"""
    initialize(keypath; variant=Southern(), faithful=true) -> (state, reason)

Convenience: open `keypath`, build a fresh state, and initialize the FIRST stand.
"""
function initialize(keypath::AbstractString; variant::AbstractVariant = Southern(),
                    faithful::Bool = true)
    s = StandState(variant; faithful = faithful)
    base = strip_key_ext(keypath)
    reason = open(keypath) do io
        initialize!(s, KeywordReader(io), base)
    end
    return s, reason
end
