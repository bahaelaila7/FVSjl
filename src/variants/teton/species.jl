# =============================================================================
# species.jl (teton) — TT species block-data init (tt/blkdat.f + tt/grinit.f)
#
# Copies the 18 TT species codes into the stand, sets the default tree format,
# 10-yr cycle, and seeds both RNG streams (55329, tt/blkdat.f DATA S0/55329D0/,
# SS/55329.). TT's grinit (tt/grinit.f) differs from EM/KT in ONE key way:
# LZEIDE=.TRUE. — TT uses ZEIDE SDI (like CR), not Stage. Otherwise DGSD=2.0,
# IFINT=10, IFINTH=5, LHTDRG=.TRUE. (except sp13 BI / sp16 MC = .FALSE.).
# TT-specific coefficients (bark, DG, crown, mortality, site, habitat) land in
# their own chunks (data/teton/).
# =============================================================================

const TT_RNG_SEED = 55329.0f0   # tt/blkdat.f DATA S0/55329D0/, SS/55329.

function init_blockdata!(s::StandState, v::Teton)
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
        sd.code2[i] = String(rstrip(first(sd.class_codes[i, 1], 2)))
    end
    hab = s.coef.valid_habitat
    copyto!(s.plot.valid_habitat, 1, hab, 1, min(length(hab), length(s.plot.valid_habitat)))

    s.control.tree_format = DEFAULT_TREE_FORMAT
    s.control.year = 10.0f0            # TT YR default cycle length (tt IFINT=10)
    s.control.growth_fint = 10.0f0     # TT FINT default = 10 (tt/grinit.f), like CR/KT/EM
    s.control.zeide_sdi = true         # TT uses ZEIDE SDI (tt/grinit.f LZEIDE=.TRUE.) — like CR, UNLIKE EM/KT
    s.rng.s0 = Float64(TT_RNG_SEED); s.rng.ss = TT_RNG_SEED   # both streams (tt/blkdat.f S0/SS)
    fill!(s.control.ht_drag_sp, true)  # LHTDRG default .TRUE. (tt/grinit.f)
    s.control.ht_drag_sp[13] = false   # sp13 BI (bigtooth maple) — tt/grinit.f LHTDRG(13)=.FALSE.
    s.control.ht_drag_sp[16] = false   # sp16 MC (curlleaf mtn-mahogany) — tt/grinit.f LHTDRG(16)=.FALSE.
    s.control.dg_sd = 2.0f0            # DGSD default (tt/grinit.f) — enables DG serial-corr + ZZRAN
    return s
end

load_species_coefficients!(s::StandState, v::Teton) = init_blockdata!(s, v)

# Species crosswalk (SPCTRN, tt spctrn.f: shared western ASPT crosswalk). Only used for FOREIGN
# input codes; ttt01 uses native alpha codes (direct class_codes match). The translation table +
# exact ASPT target column land with the crosswalk sub-chunk (data/teton/species_translation.csv).
# TODO(chunk 1b): verify the TT ASPT column vs live before committing a crosswalk. Provisional.
spctrn_column(::Teton) = 4

# Catch-all for a code absent from the crosswalk: TT SPTRN default ISPC1=MAXSP=18 (OH).
other_species(::Teton) = Int32(18)
