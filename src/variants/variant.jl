# =============================================================================
# variant.jl — the variant abstraction
#
# Replaces Fortran's file-override mechanism (where sn/dgf.f shadows base/dgf.f at
# link time) with Julia multiple dispatch. The base `engine/` code is written once
# against generic functions; each variant provides specialized methods, dispatched
# on a zero-size singleton so the compiler devirtualizes them (keeping the hotpath
# allocation-free and trim-friendly — requirement #10).
#
# Adding a variant later (CS/NE/LS) = define a new `<: AbstractVariant` singleton
# and a method for each hook below; the engine is untouched.
# =============================================================================

"""
    AbstractVariant

Base type for all FVS geographic variants. Concrete variants are singletons
(e.g. `Southern()`), carried as the type parameter of `StandState`.
"""
abstract type AbstractVariant end

"""Two-letter FVS variant code (e.g. \"SN\"). Each variant must define this."""
function variant_code end

"""
    variant_from_code(code) -> AbstractVariant

Resolve an FVS variant designator to its singleton — the inverse of `variant_code`.
Accepts the 2-letter code or the full name (case-insensitive): `SN`/`Southern`,
`NE`/`Northeast`. This is what a YAML stand file's `variant:` field maps through, so a
config selects its own model. (`Southern`/`Northeast` are defined later in the load order;
this resolves them at call time.)
"""
function variant_from_code(code::AbstractString)
    c = uppercase(strip(String(code)))
    (c == "SN" || c == "SOUTHERN")  && return Southern()
    (c == "NE" || c == "NORTHEAST") && return Northeast()
    error("unknown FVS variant '$code' (supported: SN = Southern, NE = Northeast)")
end

"""Number of species in the variant (SN=90, NE=108). Array capacity is MAXSP (the max)."""
function nspecies end

# ---------------------------------------------------------------------------
# Variant hooks. These are the generic functions the base engine calls; each
# variant supplies methods. Declared here (no methods yet) so the engine can
# reference them; methods land in src/variants/<variant>/ during later chunks.
#
# Convention: every hook takes the `StandState` first; mutating hooks end in `!`.
# ---------------------------------------------------------------------------
function load_species_coefficients! end   # blkdat.f  — fill Coefficients for the variant
function diameter_growth! end             # dgf.f / dgdriv.f
function height_growth! end               # htgf.f
function height_from_dbh end              # htcalc.f / htdbh.f
function crown_ratio! end                 # crown.f
function mortality! end                   # morts.f
function regenerate! end                  # cratet.f / regent.f
function site_setup! end                  # sitset.f
function form_class end                   # formcl.f
function bark_ratio end                   # bratio.f
function max_sdi end                      # sdimax / dgf coefficients
