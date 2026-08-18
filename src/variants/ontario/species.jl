# =============================================================================
# species.jl (ontario) — ON species block-data init (canada/on/blkdat.f + grinit.f).
# 72 species; alpha JSP + FIAJSP + PLNJSP crosswalks (blkdat.f:104-140). Class codes are
# alpha·"1"/"2"/"3" (blkdat.f NSP). grinit.f (MEASURED): LHTDRG=.FALSE. all species
# (grinit.f:159), LZEIDE=.FALSE. (grinit.f:183 — Reineke/Stage SDI branch), DGSD=2.0
# (grinit.f:219), IFINT=10 / YR=10 (grinit.f:238, blkdat.f:46), seed S0=SS=55329 (blkdat.f:306).
# =============================================================================

const ON_RNG_SEED = 55329.0f0

# ON tree-record FORMAT (canada/on/blkdat.f:42-44 TREFMT) — WIDER than the SN default
# (F5.1 DBH / 2F5.1 HT&HTG vs the SN F4.1/2F3.0), so an inline TREEDATA record's DBH and
# height land in the correct columns. Without this the shared DEFAULT_TREE_FORMAT misreads
# ON .tre records (heights → 0, spurious top-kill).
const ON_TREE_FORMAT =
    "(I4,T1,I7,F6.0,I1,A3,F5.1,F4.1,2F5.1,F5.1,I1,3(I2,I2),2I1,I2,2I3,2I1,F3.0)"

function init_blockdata!(s::StandState, v::Ontario)
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

    s.control.tree_format = ON_TREE_FORMAT
    s.control.year = 10.0f0             # IFINT (grinit.f:238)
    s.control.growth_fint = 10.0f0
    s.control.zeide_sdi = false         # ON LZEIDE=.FALSE. (grinit.f:183) — Reineke/Stage SDI
    s.rng.s0 = Float64(ON_RNG_SEED); s.rng.ss = ON_RNG_SEED
    fill!(s.control.ht_drag_sp, false)  # LHTDRG all .FALSE. (grinit.f:159)
    s.control.dg_sd = 2.0f0             # DGSD (grinit.f:219)
    return s
end

load_species_coefficients!(s::StandState, v::Ontario) = init_blockdata!(s, v)

# ON translation targets are a scoped follow-on (the copied CSV carries SN targets in col 7).
# ont01 (and every alpha/FIA/PLANTS-coded ON stand) DIRECT-matches the variant's own codes in
# resolve_species BEFORE the SPCTRN crosswalk, so this column is unused for native-coded input.
spctrn_column(::Ontario) = 7
other_species(::Ontario) = Int32(58)   # DF (sp58, "other hardwood"/generic) — placeholder; ON SPCTRN follow-on
