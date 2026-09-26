# =============================================================================
# climate_plants.jl — Climate-FVS PLANTS symbols (PLNJSP) per western variant.
#
# Each list is the variant's <v>/blkdat.f `DATA PLNJSP` verbatim (species index → USDA PLANTS code). Climate-FVS
# (clin.f INDXSPECIES) locates each species' viability column by it; with no method the shared default is empty ⇒
# ns = 0 ⇒ every climate hook (clgmult WK4, clmorts, clmaxden, clauestb) is a no-op — which is how CLIMATE ran for
# these variants until 2026-09-26 (#252). BM/EM/IE keep their own definitions. Codes blank in the Fortran stay "".
# =============================================================================

# ut/blkdat.f PLNJSP (24 species)
const _UT_PLNJSP = String[
    "PIAL", "PIFL2", "PSME", "ABCO", "PIPU", "POTR5", "PICO", "PIEN", "ABLA", "PIPO", "PIED", "JUOC", "QUGA",
    "PIMO", "JUSC2", "JUOS", "PILO", "POAN3", "POFR2", "CELE3", "ACGR3", "ACNE2", "2TN", "2TB"
]
climate_plant_symbols(::Utah) = _UT_PLNJSP

# tt/blkdat.f PLNJSP (18 species)
const _TT_PLNJSP = String[
    "PIAL", "PIFL2", "PSME", "PIMO", "PIPU", "POTR5", "PICO", "PIEN", "ABLA", "PIPO", "JUOS", "JUSC2", "ACGR3",
    "ACGL", "POAN3", "CELE3", "2TN", "2TB"
]
climate_plant_symbols(::Teton) = _TT_PLNJSP

# ci/blkdat.f PLNJSP (19 species)
const _CI_PLNJSP = String[
    "PIMO3", "LAOC", "PSME", "ABGR", "TSHE", "THPL", "PICO", "PIEN", "ABLA", "PIPO", "PIAL", "TABR2", "POTR5",
    "JUOC", "CELE3", "PIFL2", "POBAT", "2TN", "2TB"
]
climate_plant_symbols(::CentralIdaho) = _CI_PLNJSP

# cr/blkdat.f PLNJSP (38 species)
const _CR_PLNJSP = String[
    "ABLA", "ABLAA", "PSME", "ABGR", "ABCO", "TSME", "THPL", "LAOC", "PIAR", "PIFL2", "PICO", "PIED", "PIPO",
    "PIAL", "PIST3", "JUOS", "PIPU", "PIEN", "PIGL", "POTR5", "POAN3", "PODEM", "QUGA", "QUAR", "QUEM",
    "QUMA2", "QUHY", "BEPA", "JUDE2", "JUSC2", "JUMO", "JUVI", "PIMO", "PIDI3", "PIMOF", "PILE", "2TN", "2TB"
]
climate_plant_symbols(::CentralRockies) = _CR_PLNJSP

# kt/blkdat.f PLNJSP (11 species)
const _KT_PLNJSP = String[
    "PIMO3", "LAOC", "PSME", "ABGR", "TSHE", "THPL", "PICO", "PIEN", "ABLA", "PIPO", "2TREE"
]
climate_plant_symbols(::Kootenai) = _KT_PLNJSP

# nc/blkdat.f PLNJSP (12 species)
const _NC_PLNJSP = String[
    "2TN", "PILA", "PSME", "ABCO", "ARME", "CADE27", "QUKE", "LIDE3", "ABMA", "PIPO", "2TB", "SESE3"
]
climate_plant_symbols(::Klamath) = _NC_PLNJSP

# pn/blkdat.f PLNJSP (39 species)
const _PN_PLNJSP = String[
    "ABAM", "ABCO", "ABGR", "ABLA", "ABMA", "PISI", "ABPR", "CANO9", "CADE27", "PIEN", "PICO", "PIJE", "PILA",
    "PIMO3", "PIPO", "PSME", "SESE3", "THPL", "TSHE", "TSME", "ACMA3", "ALRU2", "ALRH2", "BEPA", "CHCHC4",
    "POTR5", "POBAT", "QUGA4", "JUOC", "LALY", "PIAL", "PIAT", "TABR2", "CONU4", "CRATA", "PREM", "SALIX", "",
    "2TREE"
]
climate_plant_symbols(::PacificNorthwest) = _PN_PLNJSP

# wc/blkdat.f PLNJSP (39 species)
const _WC_PLNJSP = String[
    "ABAM", "ABCO", "ABGR", "ABLA", "ABMA", "", "ABPR", "CANO9", "CADE27", "PIEN", "PICO", "PIJE", "PILA",
    "PIMO3", "PIPO", "PSME", "SESE3", "THPL", "TSHE", "TSME", "ACMA3", "ALRU2", "ALRH2", "BEPA", "CHCHC4",
    "POTR5", "POBAT", "QUGA4", "JUOC", "LALY", "PIAL", "PIAT", "TABR2", "CONU4", "CRATA", "PREM", "SALIX", "",
    "2TREE"
]
climate_plant_symbols(::WestCascades) = _WC_PLNJSP

# ec/blkdat.f PLNJSP (32 species)
const _EC_PLNJSP = String[
    "PIMO3", "LAOC", "PSME", "ABAM", "THPL", "ABGR", "PICO", "PIEN", "ABLA", "PIPO", "TSHE", "TSME", "TABR2",
    "PIAL", "ABPR", "ABCO", "LALY", "CANO9", "JUOC", "ACMA3", "ACCI", "ALRU2", "BEPA", "CHCHC4", "CONU4",
    "POTR5", "POBAT", "QUGA4", "PRUNU", "SALIX", "2TN", "2TB"
]
climate_plant_symbols(::EastCascades) = _EC_PLNJSP

# so/blkdat.f PLNJSP (33 species)
const _SO_PLNJSP = String[
    "PIMO3", "PILA", "PSME", "ABCO", "TSME", "CADE27", "PICO", "PIEN", "ABSH", "PIPO", "JUOC", "ABGR", "ABLA",
    "ABAM", "ABPR", "PIAL", "LAOC", "THPL", "TSHE", "TABR2", "ALRH2", "ALRU2", "ACMA3", "POTR5", "POBAT",
    "PREM", "QUGA4", "SALIX", "CHCHC4", "CELE3", "CEMOG", "2TN", "2TB"
]
climate_plant_symbols(::SouthCentralOregon) = _SO_PLNJSP

# ca/blkdat.f PLNJSP (50 species)
const _CA_PLNJSP = String[
    "CHLA", "CADE27", "THPL", "ABCO", "ABMA", "ABSH", "PSME", "TSHE", "TSME", "PIAL", "PIAT", "PICO", "PICO3",
    "PIFL2", "PIJE", "PILA", "PIMO3", "PIPO", "PIRA2", "PISA2", "JUOC", "PIBR", "SEGI2", "TABR2", "2TN",
    "QUAG", "QUCH2", "QUDO", "QUEN", "QUGA4", "QUKE", "QULO", "QUWI2", "ACMA3", "AECA", "ALRU2", "ARME",
    "CHCHC4", "CONU4", "FRLA", "JUGLA", "LIDE3", "PLRA", "POTR5", "POBAT", "SALIX", "TOCA", "UMCA", "2TB",
    "SESE3"
]
climate_plant_symbols(::CentralCalifornia) = _CA_PLNJSP

# ws/blkdat.f PLNJSP (43 species)
const _WS_PLNJSP = String[
    "PILA", "PSME", "ABCO", "SEGI2", "CADE27", "PIJE", "ABMA", "PIPO", "PICO", "PIAL", "PIMO3", "PIMO", "ABAM",
    "PIAT", "PIBA", "PICO3", "PIFL2", "PIRA2", "PISA2", "PIWA", "PILO", "PSMA", "SESE3", "TSME", "JUOC",
    "JUOS", "JUCA7", "QUAG", "QUCH2", "QUDO", "QUKE", "QULO", "QUWI2", "LIDE3", "CHCHC4", "POTR5", "UMCA",
    "ARME", "CONU4", "ACMA3", "CELE3", "2TN", "2TB"
]
climate_plant_symbols(::WestSierra) = _WS_PLNJSP

# oc/blkdat.f PLNJSP (50 species)
const _OC_PLNJSP = String[
    "CHLA", "CADE27", "THPL", "ABGR", "ABMA", "ABSH", "PSME", "TSHE", "TSME", "PIAL", "PIAT", "PICO", "PICO3",
    "PIFL2", "PIJE", "PILA", "PIMO3", "PIPO", "PIRA2", "PISA2", "JUOC", "PIBR", "SEGI2", "TABR2", "2TN",
    "QUAG", "QUCH2", "QUDO", "QUEN", "QUGA4", "QUKE", "QULO", "QUWI2", "ACMA3", "AECA", "ALRU2", "ARME",
    "CHCHC4", "CONU4", "FRLA", "JUGLA", "LIDE3", "PLRA", "POTR5", "POBAT", "SALIX", "TOCA", "UMCA", "2TB",
    "SESE3"
]
climate_plant_symbols(::OregonCoast) = _OC_PLNJSP

# op/blkdat.f PLNJSP (39 species)
const _OP_PLNJSP = String[
    "ABAM", "ABCO", "ABGR", "ABLA", "ABMA", "PISI", "ABPR", "CANO9", "CADE27", "PIEN", "PICO", "PIJE", "PILA",
    "PIMO3", "PIPO", "PSME", "SESE3", "THPL", "TSHE", "TSME", "ACMA3", "ALRU2", "ARME", "LIDE3", "CHCHC4",
    "POTR5", "POBAT", "QUGA4", "JUOC", "LALY", "PIAL", "PIAT", "TABR2", "CONU4", "CRATA", "PREM", "SALIX", "",
    "2TREE"
]
climate_plant_symbols(::Olympic) = _OP_PLNJSP

# bc/blkdat.f PLNJSP (15 species)
const _BC_PLNJSP = String[
    "PIMO3", "LAOC", "PSME", "ABGR", "TSHE", "THPL", "PICO", "PIEN", "ABLA", "PIPO", "BEPA", "POTR5", "POBAT",
    "SOFT", "HARD"
]
climate_plant_symbols(::BritishColumbia) = _BC_PLNJSP

# Which builds link the real Climate-FVS (clin.f, clgmult.f, ...) rather than the exclim.f stubs
# (bin/FVS<v>_buildDir/exclim.f present ⇒ stubbed: ak, on, sn, cs, ls, ne).
climate_extension_linked(::AbstractVariant) = true
climate_extension_linked(::Union{SoutheastAlaska,Ontario,Southern,CentralStates,LakeStates,Northeast}) = false
