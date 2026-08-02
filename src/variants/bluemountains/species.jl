# =============================================================================
# species.jl (bluemountains) — BM species block-data init (bm/blkdat.f + bm/grinit.f)
#
# Copies the 18 BM species codes into the stand, sets the default tree format, 10-yr cycle,
# seeds both RNG streams (55329). BM grinit: LZEIDE=.FALSE. (Stage SDI), DGSD=1.5, FINT=10,
# LHTDRG default .TRUE.. BM-specific coefficients land in their own chunks (data/bluemountains/).
# =============================================================================

const BM_RNG_SEED = 55329.0f0   # bm/blkdat.f:187 DATA S0/55329D0/,SS/55329./

function init_blockdata!(s::StandState, v::BlueMountains)
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
    s.control.year = 10.0f0            # BM YR default cycle length (bm FINT=10)
    s.control.growth_fint = 10.0f0     # BM FINT default = 10 (bm/grinit.f:203)
    s.control.zeide_sdi = false        # BM uses STAGE SDI (bm/grinit.f:157 LZEIDE=.FALSE.) — like EM/KT
    s.rng.s0 = Float64(BM_RNG_SEED); s.rng.ss = BM_RNG_SEED
    fill!(s.control.ht_drag_sp, true)  # LHTDRG default .TRUE.
    s.control.dg_sd = 1.5f0            # DGSD default (bm/grinit.f:200) — NOTE 1.5 not 2.0
    return s
end

load_species_coefficients!(s::StandState, v::BlueMountains) = init_blockdata!(s, v)

# Species crosswalk (SPCTRN, bm/spctrn.f: shared western ASPT, BM target = FOREST-DEPENDENT column 4..13;
# default col 4). NOTE (chunk 0): data/bluemountains/species_translation.csv carries EM's rows as a
# PLACEHOLDER — native BM stands resolve by direct alpha match; real BM ASPT extraction is a species TODO.
spctrn_column(::BlueMountains) = 4

# Catch-all for a code absent from the crosswalk: BM SPTRN default ISPC1=MAXSP=18 (OH).
other_species(::BlueMountains) = Int32(18)
