# =============================================================================
# species.jl (klamath) — NC species block-data init (nc/blkdat.f + nc/grinit.f). Chunk 1.
# 12 species into the stand; MEASURED nc/grinit.f defaults: LZEIDE=.TRUE. (Zeide SDI, like CI/UT/TT);
# DGSD=2.0; IFINT=10; LHTDRG default .FALSE. (nc/grinit.f:102); RNG seed 55329 (nc/blkdat.f S0/SS).
# =============================================================================

const NC_RNG_SEED = 55329.0f0   # nc/blkdat.f DATA S0/55329D0/,SS/55329./

function init_blockdata!(s::StandState, v::Klamath)
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
    s.control.year = 5.0f0             # NC YR default cycle length = 5 (nct01 live runs 5-yr; NC growth models
    s.control.growth_fint = 5.0f0      # are 5-yr-based: dgf TDDS/2, REGYR=5). IFINT=10 is the calibration interval.
    s.control.zeide_sdi = true         # NC uses ZEIDE SDI (nc/grinit.f LZEIDE=.TRUE.)
    s.rng.s0 = Float64(NC_RNG_SEED); s.rng.ss = NC_RNG_SEED
    fill!(s.control.ht_drag_sp, false) # nc/grinit.f:102 LHTDRG default .FALSE. (verify per-species exceptions ch4)
    s.control.dg_sd = 2.0f0            # DGSD default (nc/grinit.f DGSD=2.0)
    s.control.dg_stddev_bound = 2.0f0  # set BOTH (the BM/CI DGSD field-disconnect lesson)
    return s
end

load_species_coefficients!(s::StandState, v::Klamath) = init_blockdata!(s, v)

# SPCTRN crosswalk target column — NC CSV mirrors the 7-col layout (target_nc ×4) ⇒ col 4.
spctrn_column(::Klamath) = 4
other_species(::Klamath) = Int32(1)   # NC "OS" (other softwood) = species 1
