# =============================================================================
# species.jl (easternmontana) — EM species block-data init (em/blkdat.f + em/grinit.f)
#
# Copies the 19 EM species codes into the stand, sets the default tree format,
# 10-yr cycle, and seeds both RNG streams (55329, same as CR/KT/eastern). EM's
# grinit is IDENTICAL to KT (em/grinit.f): LZEIDE=.FALSE. (Stage SDI, not Zeide),
# DGSD=2.0, IFINT=10, IFINTH=5, LHTDRG=.TRUE.. EM-specific coefficients (bark, DG,
# crown, mortality, site, habitat) land in their own chunks (data/easternmontana/).
# =============================================================================

const EM_RNG_SEED = 55329.0f0   # em/blkdat.f:274 DATA S0/55329D0/, SS/55329.

# EM VOLEQDEF (em/sitset.f VEQNNC) — dumped from live FVSem em_AS.out "NVEL EQUATION NUMBERS" table.
# Conifers (1-5,7-10,18) = Region-1 Flewelling FW2 (note the volume-FIA ≠ species-FIA mapping: LM→073,
# LL/AF→019). Non-conifers (RM/GA/AS/CW/BA/PW/NC/PB/OH) = DVEW (Chojnacky/Gevorkiantz woodland, region 1/2).
const EM_VOL_EQ = String[
    "I00FW2W012", "I00FW2W073", "I00FW2W202", "I00FW2W073", "I00FW2W019", "102DVEW106",
    "I00FW2W108", "I00FW2W093", "I00FW2W019", "I00FW2W122", "101DVEW740", "102DVEW746",
    "102DVEW740", "101DVEW740", "102DVEW740", "102DVEW740", "101DVEW375", "I00FW2W260",
    "200DVEW746"]

function init_blockdata!(s::StandState, v::EasternMontana)
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
    s.control.year = 10.0f0            # EM YR default cycle length (em IFINT=10)
    s.control.growth_fint = 10.0f0     # EM FINT default = 10 (em/grinit.f:206), like CR/KT
    s.control.zeide_sdi = false        # EM uses STAGE SDI (em/grinit.f:143 LZEIDE=.FALSE.) — like KT
    s.rng.s0 = Float64(EM_RNG_SEED); s.rng.ss = EM_RNG_SEED   # both streams (em/blkdat.f S0/SS)
    fill!(s.control.ht_drag_sp, true)  # LHTDRG default .TRUE. (em/grinit.f:118)
    s.control.dg_sd = 2.0f0            # DGSD default (em/grinit.f:188) — enables DG serial-corr + ZZRAN
    return s
end

load_species_coefficients!(s::StandState, v::EasternMontana) = init_blockdata!(s, v)

# Species crosswalk (SPCTRN, em spctrn.f: shared western ASPT(442,21), EM target = column 10).
# In the 7-col CSV layout (code_alpha,code_fia,code_plants,target×4) the EM target lands in col 4.
spctrn_column(::EasternMontana) = 4

# Catch-all for a code absent from the 442-row table: EM SPTRN default ISPC1=MAXSP=19 (OH).
other_species(::EasternMontana) = Int32(19)
