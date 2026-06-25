# =============================================================================
# species.jl (northeast) — NE species block-data init (ne/blkdat.f)
#
# Mirrors the Southern init_blockdata!: copy the species codes into the stand,
# set the habitat list, default tree format, cycle length, SDI method, and seed
# both RNG streams. NE specifics: 108 species, Zeide SDI (ne/grinit.f:129
# LZEIDE=.TRUE.), RNG seed 55329 (ne/blkdat.f S0/SS, same as SN).
# =============================================================================

const NE_RNG_SEED = 55329.0f0   # ne/blkdat.f DATA S0/55329D0, SS/55329.

function init_blockdata!(s::StandState, v::Northeast)
    sd = s.species
    alpha = s.coef.code_alpha; fia = s.coef.code_fia; plants = s.coef.code_plants
    @inbounds for i in 1:nspecies(v)
        sd.alpha[i]  = alpha[i]
        sd.fia[i]    = fia[i]
        sd.plants[i] = plants[i]
        code = rstrip(alpha[i])
        sd.class_codes[i, 1] = code * "1"
        sd.class_codes[i, 2] = code * "2"
        sd.class_codes[i, 3] = code * "3"
    end
    hab = s.coef.valid_habitat
    copyto!(s.plot.valid_habitat, 1, hab, 1, min(length(hab), length(s.plot.valid_habitat)))

    s.control.tree_format = DEFAULT_TREE_FORMAT
    s.control.year = 5.0f0                                # YR default cycle length
    s.control.zeide_sdi = true                            # NE uses Zeide SDI (ne/grinit.f:129)

    # RNG: both streams seeded 55329 (blkdat S0/SS + ESBLKD ESS0/ESSS)
    s.rng.s0 = Float64(NE_RNG_SEED); s.rng.ss = NE_RNG_SEED
    return s
end

"""
    load_species_coefficients!(s, ::Northeast)

Variant hook: load NE species data into the stand state (the BLOCK-DATA load).
Growth/volume coefficients are attached as their subsystems are ported.
"""
load_species_coefficients!(s::StandState, v::Northeast) = init_blockdata!(s, v)

# Species-translation crosswalk: NE's target code is column 6 (alpha,fia,plants,cs,ls,NE,sn).
spctrn_column(::Northeast) = 6
# Catch-all "other" species for an unmatched code: OH (other hardwood, index 97). TODO verify vs ne/spctrn.f.
other_species(::Northeast) = Int32(97)
