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
    # canada/on/grinit.f never assigns CCCOEF/CCCOEF2 (ls/grinit.f:264-265 and bc/grinit.f:275-276 set 1.0), so the
    # COMMON keeps 0: COVOLP's PCCU = CCCOEF·ΣCRAREA/43560 = 0 ⇒ every FVS_StrClass_Metric cover is 0 (live ON goldens),
    # and cuts.f:850-851 CCC = 0 for the THINCC CRADIF.
    s.control.cc_coef = 0f0; s.control.cc_coef2 = 0f0
    return s
end

load_species_coefficients!(s::StandState, v::Ontario) = init_blockdata!(s, v)

# ON translation targets are a scoped follow-on (the copied CSV carries SN targets in col 7).
# ont01 (and every alpha/FIA/PLANTS-coded ON stand) DIRECT-matches the variant's own codes in
# resolve_species BEFORE the SPCTRN crosswalk, so this column is unused for native-coded input.
spctrn_column(::Ontario) = 7
other_species(::Ontario) = Int32(49)   # NC — canada/on/spctrn.f CASE DEFAULT ! non-commercial hardwood

# canada/on/spctrn.f: ON does NOT use the shared 562-row ASPT crosswalk — a code that matched none of the variant's own
# JSP/FIAJSP/PLNJSP lists (intree.f:240-263) goes to `SELECT CASE (SPCIN(:2))` over the ON alpha codes, and anything
# else (e.g. FIA 356 serviceberry, not in FIAJSP) to CASE DEFAULT = 49 (NC). The shared crosswalk sent FIA 356 to ML.
const _ON_SPCTRN = Dict{String,Int32}(
    "PJ"=>1, "PS"=>2, "RN"=>3, "RP"=>4, "PW"=>5, "SW"=>6, "SN"=>7, "BF"=>8, "SB"=>9, "TA"=>10, "CE"=>11, "HE"=>12,
    "SO"=>13, "CR"=>14, "AB"=>15, "AR"=>16, "CW"=>17, "MV"=>18, "MR"=>19, "CB"=>20, "EW"=>21, "ES"=>22, "ER"=>23,
    "BY"=>24, "BD"=>25, "MH"=>26, "MB"=>27, "BE"=>28, "AW"=>29, "OW"=>30, "OP"=>31, "OB"=>32, "OC"=>33, "OR"=>34,
    "BO"=>35, "PN"=>36, "HB"=>37, "HP"=>38, "HU"=>39, "PG"=>40, "PT"=>41, "PB"=>42, "BW"=>43, "CH"=>44, "BT"=>45,
    "WB"=>46, "IW"=>47, "LB"=>48, "NC"=>49, "MM"=>50, "MS"=>51, "MT"=>52, "BB"=>53, "CA"=>54, "BH"=>55, "DF"=>56,
    "HT"=>57, "ML"=>58, "GB"=>59, "SY"=>60, "CP"=>61, "CC"=>62, "PL"=>63, "WI"=>64, "SS"=>67, "AM"=>68, "JP"=>69,
    "WP"=>70, "SP"=>71, "BP"=>72)
function translate_species(code::AbstractString, ::Ontario, sp::SpeciesData, coef::SpeciesCoefficients)
    c = rstrip(code)
    fmt = (!isempty(c) && all(isdigit, c)) ? Int32(2) : Int32(1)
    return (get(_ON_SPCTRN, first(c, 2), Int32(49)), fmt)
end
