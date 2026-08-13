# =============================================================================
# species.jl (oregoncoast) — OC species block-data init (oc/blkdat.f + oc/grinit.f).
#
# Copies the 50 OC species codes into the stand and sets the ORGANON-style control defaults
# (MEASURED, oc/grinit.f): 5-yr cycle (FINT=5), Stage SDI (LZEIDE=.FALSE.), DGSD=0 (no DG serial
# correlation — ORGANON is self-contained), LHTDRG=.FALSE. all species. RNG seed 55329 (shared,
# oc/blkdat.f). ORGANON coefficients (DG/HTG/CR/mortality/volume) live in the ORGANON subsystem
# (organon/ + vorganon/, UNPORTED) — this scaffold carries only the species codes.
# =============================================================================

const OC_RNG_SEED = 55329.0f0   # oc/blkdat.f DATA S0/SS (shared across all variants)

function init_blockdata!(s::StandState, v::OregonCoast)
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
    s.control.year = 5.0f0             # OC YR default cycle length (oc FINT=5) — ORGANON 5-yr step
    s.control.growth_fint = 5.0f0      # OC FINT default = 5 (oc/grinit.f:175)
    s.control.zeide_sdi = false        # OC uses STAGE SDI (oc/grinit.f:129 LZEIDE=.FALSE.)
    s.rng.s0 = Float64(OC_RNG_SEED); s.rng.ss = OC_RNG_SEED
    fill!(s.control.ht_drag_sp, false) # LHTDRG all .FALSE. (oc/grinit.f:105) — ht from ORGANON
    s.control.dg_sd = 0.0f0            # DGSD=0 (oc/grinit.f:172) — no DG serial-corr
    return s
end

load_species_coefficients!(s::StandState, v::OregonCoast) = init_blockdata!(s, v)

# Species crosswalk: OC SPCTRN column (oc/spctrn.f). Placeholder col 4 (mirrors 7-col layout);
# native OC stands resolve by direct alpha match. Real ASPT extraction is a species-chunk TODO.
spctrn_column(::OregonCoast) = 4

# Catch-all for a code absent from the crosswalk: OC SPCTRN default (last species OH-ish). The
# live .out remaps unknown conifers into ORGANON species; scaffold defaults to OH (sp49).
other_species(::OregonCoast) = Int32(49)
