# =============================================================================
# species.jl (lakestates) — LS species block-data init (ls/blkdat.f)
#
# Mirrors the SN/NE/CS init_blockdata!: copy the 68 LS species codes into the stand,
# default tree format, cycle length, SDI method, and seed both RNG streams. LS
# specifics: 68 species, Zeide SDI (ls/grinit.f LZEIDE=.TRUE.), RNG seed 55329
# (ls/blkdat.f S0/SS, same as SN/NE/CS), YR=10 (10-yr cycle, like NE/CS).
# =============================================================================

const LS_RNG_SEED = 55329.0f0   # ls/blkdat.f DATA S0/55329D0, SS/55329.

function init_blockdata!(s::StandState, v::LakeStates)
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
    s.control.year = 10.0f0                               # LS YR default cycle length (ls/blkdat.f:71)
    s.control.zeide_sdi = true                            # LS uses Zeide SDI (ls/grinit.f:126 LZEIDE=.TRUE.)

    # RNG: both streams seeded 55329 (blkdat S0/SS + ESBLKD ESS0/ESSS)
    s.rng.s0 = Float64(LS_RNG_SEED); s.rng.ss = LS_RNG_SEED
    return s
end

"""
    load_species_coefficients!(s, ::LakeStates)

Variant hook: load LS species data into the stand state (the BLOCK-DATA load).
Growth/volume coefficients attach as their subsystems are ported.
"""
load_species_coefficients!(s::StandState, v::LakeStates) = init_blockdata!(s, v)

# Species-translation crosswalk: LS's target code is column 5 (alpha,fia,plants,cs,ls,ne,sn)
# — cf. ls/spctrn.f CASE('LS') → ASPT(I,5).
spctrn_column(::LakeStates) = 5
# Unmatched/non-commercial species code → LS index 49 (OH, Other Hardwood; ls/spctrn.f: VAR=='LS' ⇒ ISPC1=49).
other_species(::LakeStates) = Int32(49)
