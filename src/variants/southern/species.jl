# =============================================================================
# species.jl — Southern variant species tables & BLOCK DATA defaults
#
# Ported from: sn/blkdat.f  (BLOCK DATA BLKDAT)
#
# The 90 Southern species and their identity codes (alpha / FIA SPCD / PLANTS),
# plus the diameter-growth regression standard errors and the BLOCK DATA scalar
# defaults (RNG seed, default tree-record FORMAT, etc.). These immutable tables
# live as `const` arrays in the variant; per-stand state is filled by
# `load_species_coefficients!` / `init_blockdata!` below.
# =============================================================================

"Number of Southern variant species."
const SN_NSPECIES = 90

# Species identity codes (alpha/FIA/PLANTS), DG regression standard errors (SIGMAR =
# dg_resid_sd), and the valid-habitat list now live in the variant's CSV data — see
# data/southern/{species_coefficients,valid_habitat_codes}.csv. They are read once at
# startup into the cached `SpeciesCoefficients` and reached through `s.coef`.

# BLOCK DATA scalar defaults (sn/blkdat.f). The Fortran unit numbers are kept as
# named fields rather than a global unit table (see io layer). The RNG main stream
# is seeded to 55329 here (blkdat sets S0=SS=55329), matching the establishment seed.
const SN_RNG_SEED = 55329.0f0
const SN_REGEN_BARK = 2.999f0      # REGNBK

"""
    init_blockdata!(s::StandState, ::Southern)

Apply the Southern BLOCK DATA defaults to a fresh state: species identity tables,
default tree FORMAT, RNG seed, and the handful of scalar defaults set in blkdat.f.
Called once at the start of stand initialization.
"""
function init_blockdata!(s::StandState, ::Southern)
    sd = s.species
    alpha = s.coef.code_alpha; fia = s.coef.code_fia; plants = s.coef.code_plants
    @inbounds for i in 1:SN_NSPECIES
        sd.alpha[i]  = alpha[i]
        sd.fia[i]    = fia[i]
        sd.plants[i] = plants[i]
        code = rstrip(alpha[i])
        sd.class_codes[i, 1] = code * "1"
        sd.class_codes[i, 2] = code * "2"
        sd.class_codes[i, 3] = code * "3"
        sd.code2[i] = String(rstrip(first(sd.class_codes[i, 1], 2)))  # pre-rstripped crown-width key
    end
    hab = s.coef.valid_habitat
    copyto!(s.plot.valid_habitat, 1, hab, 1, min(length(hab), length(s.plot.valid_habitat)))

    s.control.tree_format = DEFAULT_TREE_FORMAT
    s.control.year = 5.0f0                                # YR default cycle length
    s.control.zeide_sdi = true                            # SN uses Zeide/Reineke SDI

    # RNG: both streams seeded 55329 (blkdat S0/SS + ESBLKD ESS0/ESSS)
    s.rng.s0 = Float64(SN_RNG_SEED); s.rng.ss = SN_RNG_SEED
    return s
end

"""
    load_species_coefficients!(s, ::Southern)

Variant hook: load species data into the stand state. For now this is the
identity/BLOCK-DATA load; growth coefficients (DGF/HTGF tables) are attached in C3.
"""
load_species_coefficients!(s::StandState, v::Southern) = init_blockdata!(s, v)
