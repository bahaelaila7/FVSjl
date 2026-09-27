# =============================================================================
# species.jl (britishcolumbia) — BC species block-data init (bc/blkdat.f + bc/grinit.f).
# Mirrors IE (15 species, default tree format, 10-yr cycle, seed 55329). BC grinit (MEASURED):
# LHTDRG=.TRUE. all; LZEIDE=.FALSE. (STAGE SDI, like IE); DGSD=2.0; IFINT=10. Bark constant
# (BARK2=0, IMAP=2 ⇒ BRATIO=BARK1) — shared linear bark_ratio handles it (bark_a=0, bark_b=BARK1).
# =============================================================================

const BC_RNG_SEED = 55329.0f0

# bc/blkdat.f:37-39 TREFMT — the metric-variant record layout (F5.1 DBH, F4.1 DG, 2F5.1 HT/THT), the same as ON.
# The shared SN default (F4.1/F3.1/2F3.0) misreads a BC .tre record: HT lands blank and the real height is read
# as the top-kill height THT, so every tree became top-killed at its own height (ITRUNC=HT).
const BC_TREE_FORMAT =
    "(I4,T1,I7,F6.0,I1,A3,F5.1,F4.1,2F5.1,F5.1,I1,3(I2,I2),2I1,I2,2I3,2I1,F3.0)"

function init_blockdata!(s::StandState, v::BritishColumbia)
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

    s.control.tree_format = BC_TREE_FORMAT
    s.control.year = 10.0f0
    s.control.growth_fint = 10.0f0
    s.control.zeide_sdi = false        # BC uses STAGE SDI (bc/grinit.f LZEIDE=.FALSE.)
    s.rng.s0 = Float64(BC_RNG_SEED); s.rng.ss = BC_RNG_SEED
    fill!(s.control.ht_drag_sp, true)  # LHTDRG all .TRUE.
    s.control.dg_sd = 2.0f0            # DGSD (bc/grinit.f)
    return s
end

load_species_coefficients!(s::StandState, v::BritishColumbia) = init_blockdata!(s, v)

spctrn_column(::BritishColumbia) = 4        # BC translation CSV mirrors 7-col layout ⇒ col 4
other_species(::BritishColumbia) = Int32(14)   # FVS maps unrecognized codes → OC (sp14): oracle .out
                                                # "INPUT SPECIES CODE (BC) WAS SET TO (OC)". NOT OH(15).
