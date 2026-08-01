# =============================================================================
# species.jl (inlandempire) — IE species block-data init (ie/blkdat.f + ie/grinit.f).
# Mirrors KT: 23 species codes into the stand, default tree format, 10-yr cycle, seed 55329.
# IE grinit defaults IDENTICAL to KT (MEASURED ie/grinit.f): LHTDRG=.TRUE.(:105), LZEIDE=.FALSE.(:130,
# STAGE SDI), DGSD=2.0(:175), FINT=10(:178), IFINT=10(:193); RNG seed 55329 (ie/blkdat.f:237).
# =============================================================================

const IE_RNG_SEED = 55329.0f0   # ie/blkdat.f:237 DATA S0/55329D0/, SS/55329./

function init_blockdata!(s::StandState, v::InlandEmpire)
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
    s.control.year = 10.0f0            # IE YR default cycle length (IFINT=10)
    s.control.growth_fint = 10.0f0     # IE FINT default = 10 (ie/grinit.f:178)
    s.control.zeide_sdi = false        # IE uses STAGE SDI (ie/grinit.f:130 LZEIDE=.FALSE.)
    s.rng.s0 = Float64(IE_RNG_SEED); s.rng.ss = IE_RNG_SEED   # both streams (ie/blkdat.f S0/SS)
    fill!(s.control.ht_drag_sp, true)  # LHTDRG default .TRUE. (ie/grinit.f:105)
    s.control.dg_sd = 2.0f0            # DGSD default (ie/grinit.f:175)
    return s
end

load_species_coefficients!(s::StandState, v::InlandEmpire) = init_blockdata!(s, v)

# SPCTRN crosswalk target column (shared western ASPT). IE CSV mirrors KT's 7-col layout ⇒ col 4.
spctrn_column(::InlandEmpire) = 4
other_species(::InlandEmpire) = Int32(23)   # IE "OS" (other softwood) = species 23
