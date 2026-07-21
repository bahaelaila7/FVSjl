# =============================================================================
# species.jl (centralrockies) — CR species block-data init (cr/blkdat.f + grinit.f)
#
# Copies the 38 CR species codes into the stand, sets the default tree format,
# 10-yr cycle, Zeide SDI (cr/grinit.f:134 LZEIDE=.TRUE.), and seeds both RNG
# streams (55329, same as the eastern variants). CR-specific coefficients (bark,
# DG, crown, mortality, habitat) land in their own chunks.
# =============================================================================

const CR_RNG_SEED = 55329.0f0   # cr/blkdat.f DATA S0/55329D0, SS/55329.

function init_blockdata!(s::StandState, v::CentralRockies)
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
    s.control.year = 10.0f0            # CR YR default cycle length (cr/blkdat.f:126)
    s.control.zeide_sdi = true         # CR uses Zeide SDI (cr/grinit.f:134 LZEIDE=.TRUE.)
    s.rng.s0 = Float64(CR_RNG_SEED); s.rng.ss = CR_RNG_SEED   # both streams (cr/blkdat.f S0/SS)
    fill!(s.control.ht_drag_sp, true)  # LHTDRG default .TRUE. (cr/grinit.f:110) — enables cratet AA HT-DBH fit
    return s
end

load_species_coefficients!(s::StandState, v::CentralRockies) = init_blockdata!(s, v)

# Species crosswalk (SPCTRN, bin/FVScr_buildDir/spctrn.f ASPT(442,21)): the shared western table's
# CR target is column J=8; in the 7-col eastern CSV layout it lands in the first target slot (col 4).
spctrn_column(::CentralRockies) = 4
# Catch-all for a code absent from the 442-row table (rare — the table is FIA/PLANTS-comprehensive):
# other-hardwood OH (sp 38). Softwoods that miss resolve to OS via the table itself.
other_species(::CentralRockies) = Int32(38)
