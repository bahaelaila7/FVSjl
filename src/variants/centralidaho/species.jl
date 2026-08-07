# =============================================================================
# species.jl (centralidaho) — CI species block-data init (ci/blkdat.f + ci/grinit.f).
# Mirrors IE (19 species into the stand, default tree format, 10-yr cycle, seed 55329).
# CI grinit defaults (MEASURED ci/grinit.f): LHTDRG=.TRUE. all but sp15(MC); LZEIDE=.TRUE.
# (ZEIDE SDI — differs from IE's STAGE); DGSD=1.7; IFINT=10. RNG seed 55329 (ci/blkdat.f S0/SS).
# =============================================================================

const CI_RNG_SEED = 55329.0f0   # ci/blkdat.f DATA S0/55329D0/, SS/55329./

function init_blockdata!(s::StandState, v::CentralIdaho)
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
    s.control.year = 10.0f0            # CI YR default cycle length (IFINT=10)
    s.control.growth_fint = 10.0f0     # CI FINT default = 10 (ci/grinit.f)
    s.control.zeide_sdi = true         # CI uses ZEIDE SDI (ci/grinit.f LZEIDE=.TRUE.) ← differs from IE
    s.rng.s0 = Float64(CI_RNG_SEED); s.rng.ss = CI_RNG_SEED
    fill!(s.control.ht_drag_sp, true)  # LHTDRG default .TRUE. (ci/grinit.f)
    s.control.ht_drag_sp[15] = false   # ci/grinit.f LHTDRG(15)=.FALSE. (MC)
    s.control.dg_sd = 1.7f0            # DGSD default (ci/grinit.f DGSD=1.7)
    s.control.dg_stddev_bound = 1.7f0  # ci/grinit.f DGSD=1.7 also bounds the large-tree DGSCOR reject/clamp
                                       # (was SN default 2.0) — same disconnected-field bug as BM #140
    return s
end

load_species_coefficients!(s::StandState, v::CentralIdaho) = init_blockdata!(s, v)

# SPCTRN crosswalk target column — CI CSV mirrors the 7-col layout (target_ci ×4) ⇒ col 4.
spctrn_column(::CentralIdaho) = 4
other_species(::CentralIdaho) = Int32(18)   # CI "OS" (other softwood) = species 18
