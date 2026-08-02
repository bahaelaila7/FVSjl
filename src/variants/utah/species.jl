# =============================================================================
# species.jl (utah) — UT species block-data init (ut/blkdat.f + ut/grinit.f)
#
# Copies the 24 UT species codes into the stand, sets the default tree format,
# 10-yr cycle, and seeds both RNG streams (55329, same as CR/KT/EM/eastern). UT
# grinit (ut/grinit.f): LZEIDE=.TRUE. (Zeide SDI — like CR, UNLIKE KT/EM), DGSD=2.0,
# FINT=10, FINTH=5, LHTDRG=.TRUE. except MC(20)/BI(21). UT-specific coefficients (bark,
# DG, crown, mortality, site, habitat) land in their own chunks (data/utah/).
# =============================================================================

const UT_RNG_SEED = 55329.0f0   # ut/blkdat.f:201 DATA S0/55329D0/,SS/55329./

function init_blockdata!(s::StandState, v::Utah)
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
    s.control.year = 10.0f0            # UT YR default cycle length (ut FINT=10)
    s.control.growth_fint = 10.0f0     # UT FINT default = 10 (ut/grinit.f:177)
    s.control.zeide_sdi = true         # UT uses ZEIDE SDI (ut/grinit.f:133 LZEIDE=.TRUE.) — like CR
    s.rng.s0 = Float64(UT_RNG_SEED); s.rng.ss = UT_RNG_SEED   # both streams (ut/blkdat.f S0/SS)
    fill!(s.control.ht_drag_sp, true)  # LHTDRG default .TRUE. (ut/grinit.f:103)
    s.control.ht_drag_sp[20] = false   # MC (ut/grinit.f:125 LHTDRG(20)=.FALSE.)
    s.control.ht_drag_sp[21] = false   # BI (ut/grinit.f:126 LHTDRG(21)=.FALSE.)
    s.control.dg_sd = 2.0f0            # DGSD default (ut/grinit.f:174) — enables DG serial-corr + ZZRAN
    return s
end

load_species_coefficients!(s::StandState, v::Utah) = init_blockdata!(s, v)

# Species crosswalk (SPCTRN, ut/spctrn.f: shared western ASPT(442,21), UT target = column 4).
# In the 7-col CSV layout (code_alpha,code_fia,code_plants,target×4) the UT target lands in col 4.
# NOTE (chunk 0): data/utah/species_translation.csv currently carries EM's target values as a
# PLACEHOLDER — native UT stands resolve by direct alpha match (resolve_species), so this only
# affects non-native input codes; real UT ASPT col-4 extraction is a species-chunk TODO.
spctrn_column(::Utah) = 4

# Catch-all for a code absent from the crosswalk: UT SPCTRN default ISPC1=MAXSP=24 (OH).
other_species(::Utah) = Int32(24)
