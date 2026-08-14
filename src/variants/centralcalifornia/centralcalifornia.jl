# =============================================================================
# centralcalifornia.jl — the CA (Inland California / CentralCalifornia) variant singleton + registration.
#
# Ported from: ca/*.f (FVS "CA"). CA is a westside Wykoff-DDS variant on the SAME shared engine as
# West Cascades (WC) / Pacific Northwest (PN) / East Cascades (EC) — 10-yr native period, IFINT=10/
# IFINTH=5, RNG seed 55329, DGSD=1.7 — BUT with two CA-specific differences from the R6 westside trio:
#   • ZEIDE SDI (ca/grinit.f LZEIDE=.TRUE.), unlike WC/PN/EC's Reineke/Stage SDI.
#   • MAXSP=50 with a GROUP-COMPRESSED dgf: 50 species map (MAPSPC) onto 13 diameter-growth equations;
#     every dgf coefficient array is length 13 and indexed by the growth-group JSPC (like WC/PN, unlike EC).
# Species JSP (50): PC IC RC WF RF SH DF WH MH WB KP LP CP LM JP SP WP PP MP GP WJ BR GS PY OS LO CY BL EO
#   WO BO VO IO BM BU RA MA GC DG FL WN TO SY AS CW WI CN CL OH RW
# Oracle: /workspace/.cawork/FVSca_clean (relink bin/FVSca_buildDir/*.o + isoc23 shim);
#   instrumentable /workspace/.cawork/FVSca_dump (dgf.f fort-unit DGCON/per-tree-DDS dump).
#   Canonical stand: cat01 (forest 610→IFOR 6, elev 35 (hundred-ft), slope 30%, aspect 315°).
# =============================================================================

"""
    CentralCalifornia <: AbstractVariant

The FVS Inland California variant (VARACD "CA", MAXSP = 50) — a westside Wykoff-DDS variant (Zeide SDI,
DGSD 1.7) whose 50 species map onto 13 group-compressed diameter-growth equations.
"""
struct CentralCalifornia <: AbstractVariant end

variant_code(::CentralCalifornia) = "CA"
nspecies(::CentralCalifornia) = 50
htg_period(::CentralCalifornia) = 10f0   # /CONTRL/ YR=10 (ca/grinit.f FINT=10), like WC/PN/EC
mort_ri_scale(::CentralCalifornia) = 0.5f0   # background half-rate (ca/morts.f), like the westside trio

const CA_RNG_SEED = 55329.0f0            # ca/blkdat.f DATA S0/55329D0/,SS/55329./

# ca/blkdat.f JSP/FIAJSP/PLNJSP (50 species).
const CA_JSP = String[
    "PC","IC","RC","WF","RF","SH","DF","WH","MH","WB","KP","LP","CP","LM","JP","SP","WP","PP","MP","GP",
    "WJ","BR","GS","PY","OS","LO","CY","BL","EO","WO","BO","VO","IO","BM","BU","RA","MA","GC","DG","FL",
    "WN","TO","SY","AS","CW","WI","CN","CL","OH","RW"]
const CA_FIA = String[
    "041","081","242","015","020","021","202","263","264","101","103","108","109","113","116","117","119","122","124","127",
    "064","092","212","231","299","801","805","807","811","815","818","821","839","312","333","351","361","431","492","542",
    "600","631","730","746","747","920","251","981","998","211"]
const CA_PLANTS = String[
    "CHLA","CADE27","THPL","ABCO","ABMA","ABSH","PSME","TSHE","TSME","PIAL","PIAT","PICO","PICO3","PIFL2","PIJE","PILA","PIMO3","PIPO","PIRA2","PISA2",
    "JUOC","PIBR","SEGI2","TABR2","2TN","QUAG","QUCH2","QUDO","QUEN","QUGA4","QUKE","QULO","QUWI2","ACMA3","AECA","ALRU2","ARME","CHCHC4","CONU4","FRLA",
    "JUGLA","LIDE3","PLRA","POTR5","POBAT","SALIX","TOCA","UMCA","2TB","SESE3"]

"""
    ca_grinit!(s)

Set CA's /CONTRL/ + SDI-family flags (ca/grinit.f + blkdat.f): ZEIDE SDI (LZEIDE=.TRUE.), DGSD 1.7,
10-yr native period, RNG seed 55329, LHTDRG default .FALSE.
"""
function ca_grinit!(s::StandState)
    s.control.year = 10.0f0
    s.control.growth_fint = 10.0f0
    # ca/grinit.f sets LZEIDE=.TRUE. + CALCSDI=' ', but ca/sitset.f:95 `IF(CALCSDI.EQ.' ')LZEIDE=.FALSE.`
    # RESETS it to false (no CALCSDI keyword) ⇒ the reported .sum SDI + Stage/Reineke mortality use REINEKE,
    # not Zeide. (The DGF's point relative density PRD uses point-Zeide SDICZ separately, unaffected.)
    s.control.zeide_sdi = false
    s.control.dg_sd = 1.7f0             # ca/grinit.f DGSD=1.7
    s.control.dg_stddev_bound = 1.7f0
    s.rng.s0 = Float64(CA_RNG_SEED); s.rng.ss = CA_RNG_SEED
    fill!(s.control.ht_drag_sp, false)  # ca/grinit.f LHTDRG default .FALSE.
    return s
end
