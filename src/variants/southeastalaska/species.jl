# =============================================================================
# species.jl (southeastalaska) — AK species block-data init (ak/blkdat.f + ak/grinit.f)
#
# Copies the 23 AK species codes into the stand, sets the default tree format, 10-yr
# cycle, and seeds both RNG streams (55329, shared blkdat). AK grinit (ak/grinit.f):
# LZEIDE=.TRUE. (Zeide SDI — like CR/UT), DGSD=2.0, FINT=10, FINTH=5, FINTM=5,
# LHTDRG(I)=.FALSE. (default — UNLIKE UT). AK-specific coefficients (bark, DG, permafrost,
# crown, mortality, site, volume) land in their own chunks.
# =============================================================================

const AK_RNG_SEED = 55329.0f0   # shared blkdat DATA S0/55329D0/,SS/55329./ (all variants)

function init_blockdata!(s::StandState, v::SoutheastAlaska)
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
    s.control.year = 10.0f0            # AK YR default cycle length (ak FINT=10)
    s.control.growth_fint = 10.0f0     # AK FINT default = 10 (ak/grinit.f:179)
    s.control.zeide_sdi = true         # AK uses ZEIDE SDI (ak/grinit.f:135 LZEIDE=.TRUE.) — like CR/UT
    s.rng.s0 = Float64(AK_RNG_SEED); s.rng.ss = AK_RNG_SEED   # both streams (shared blkdat S0/SS)
    fill!(s.control.ht_drag_sp, false) # LHTDRG(I)=.FALSE. default (ak/grinit.f:133) — UNLIKE UT
    s.control.dg_sd = 2.0f0            # DGSD default (ak/grinit.f:176) — enables DG serial-corr + ZZRAN
    return s
end

load_species_coefficients!(s::StandState, v::SoutheastAlaska) = init_blockdata!(s, v)

# Species crosswalk (SPCTRN, ak spctrn.f: shared western ASPT(442,21), AK target = column 4).
spctrn_column(::SoutheastAlaska) = 4

# Catch-all for a code absent from the crosswalk: AK default = MAXSP=23 (OH other hardwood).
other_species(::SoutheastAlaska) = Int32(23)
