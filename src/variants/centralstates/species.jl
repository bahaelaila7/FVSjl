# =============================================================================
# species.jl (centralstates) — CS species block-data init (cs/blkdat.f)
#
# Mirrors the SN/NE init_blockdata!: copy the 96 CS species codes into the stand,
# default tree format, cycle length, SDI method, and seed both RNG streams. CS
# specifics: 96 species, Zeide SDI (cs/grinit.f LZEIDE=.TRUE.), RNG seed 55329
# (cs/blkdat.f S0/SS, same as SN/NE), YR=10 (10-yr cycle).
# =============================================================================

const CS_RNG_SEED = 55329.0f0   # cs/blkdat.f DATA S0/55329D0, SS/55329.

function init_blockdata!(s::StandState, v::CentralStates)
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
        sd.code2[i] = String(rstrip(first(sd.class_codes[i, 1], 2)))  # pre-rstripped crown-width key
    end
    hab = s.coef.valid_habitat
    copyto!(s.plot.valid_habitat, 1, hab, 1, min(length(hab), length(s.plot.valid_habitat)))

    s.control.tree_format = DEFAULT_TREE_FORMAT
    s.control.year = 10.0f0                               # CS YR default cycle length (cs/grinit.f FINT=10)
    s.control.zeide_sdi = true                            # CS uses Zeide SDI (cs/grinit.f LZEIDE=.TRUE.)

    # RNG: both streams seeded 55329 (blkdat S0/SS + ESBLKD ESS0/ESSS)
    s.rng.s0 = Float64(CS_RNG_SEED); s.rng.ss = CS_RNG_SEED
    return s
end

"""
    load_species_coefficients!(s, ::CentralStates)

Variant hook: load CS species data into the stand state (the BLOCK-DATA load).
Growth/volume coefficients attach as their subsystems are ported.
"""
load_species_coefficients!(s::StandState, v::CentralStates) = init_blockdata!(s, v)

# Species-translation crosswalk: CS's target code is column 4 (alpha,fia,plants,cs,ls,ne,sn)
# — cf. cs/spctrn.f CASE('CS') → ASPT(I,4).
spctrn_column(::CentralStates) = 4
# Unmatched/non-commercial species code → CS index 85 (cs/spctrn.f: VAR=='CS' ⇒ ISPC1=85).
other_species(::CentralStates) = Int32(85)
