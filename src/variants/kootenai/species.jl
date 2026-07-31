# =============================================================================
# species.jl (kootenai) — KT species block-data init (kt/blkdat.f + kt/grinit.f)
#
# Copies the 11 KT species codes into the stand, sets the default tree format,
# 10-yr cycle, and seeds both RNG streams (55329, same as CR/eastern). KT-specific
# coefficients (bark, DG, crown, mortality, site, habitat) land in their own chunks.
#
# KT vs CR grinit differences (MEASURED, kt/grinit.f): LZEIDE=.FALSE. (KT uses STAGE
# SDI, not CR's Zeide) — the SDI density-mortality threshold path differs. Shared with
# CR: FINT=10 (grinit.f:173), LHTDRG=.TRUE. (:101), DGSD=2.0 (:170), IFINT=10/IFINTH=5.
# =============================================================================

const KT_RNG_SEED = 55329.0f0   # kt/blkdat.f:196 DATA S0/55329D0/, SS/55329.

function init_blockdata!(s::StandState, v::Kootenai)
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
    s.control.year = 10.0f0            # KT YR default cycle length (kt IFINT=10)
    s.control.growth_fint = 10.0f0     # KT FINT default = 10 (kt/grinit.f:173), like CR
    s.control.zeide_sdi = false        # KT uses STAGE SDI (kt/grinit.f:125 LZEIDE=.FALSE.) — UNLIKE CR
    s.rng.s0 = Float64(KT_RNG_SEED); s.rng.ss = KT_RNG_SEED   # both streams (kt/blkdat.f S0/SS)
    fill!(s.control.ht_drag_sp, true)  # LHTDRG default .TRUE. (kt/grinit.f:101)
    s.control.dg_sd = 2.0f0            # DGSD default (kt/grinit.f:170) — enables DG serial-corr + ZZRAN
    return s
end

load_species_coefficients!(s::StandState, v::Kootenai) = init_blockdata!(s, v)

# Species crosswalk (SPCTRN, kt/spctrn.f: shared western ASPT(442,21), KT target = column 12).
# In the 7-col CSV layout (code_alpha,code_fia,code_plants,target×4) the KT target lands in col 4.
spctrn_column(::Kootenai) = 4

# Catch-all for a code absent from the 442-row table: KT "OT" (other tree) = species 11 (code_fia 999).
other_species(::Kootenai) = Int32(11)
