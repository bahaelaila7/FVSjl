# =============================================================================
# eastcascades.jl — the EC (East Cascades) variant singleton + registration.
#
# Ported from: ec/*.f (FVS "EC"). EC is a westside R6 Wykoff-DDS variant on the SAME shared engine as
# West Cascades (WC) / Pacific Northwest (PN) — Reineke/Stage SDI (LZEIDE=.FALSE.), DGSD=1.7,
# IFINT=10/IFINTH=5, native 10-yr period, RNG seed 55329 — but with a DIFFERENT 32-species table and its
# own coefficient DATA. Unlike PN (which shared WC's group-compressed dgf coefficients), EC's dgf.f is
# NOT group-compressed: every coefficient array is length MAXSP=32 and indexed directly by ISPC.
# Species JSP (32): WP WL DF SF RC GF LP ES AF PP WH MH PY WB NF WF LL YC WJ BM VN RA PB GC DG AS CW WO
#   PL WI OS OH
# Oracle: /workspace/.ecwork/FVSec_clean (relink bin/FVSec_buildDir/*.o + isoc23 shim; relink_ec.sh);
#   instrumentable /workspace/.ecwork/FVSec_g16. Canonical stand: ect01 (forest 608, habitat 12, S248112).
# =============================================================================

"""
    EastCascades <: AbstractVariant

The FVS East Cascades variant (VARACD "EC", MAXSP = 32) — a westside R6 Wykoff-DDS variant (Reineke SDI,
DGSD 1.7) with its own 32-species table and (uncompressed) per-species coefficient DATA.
"""
struct EastCascades <: AbstractVariant end

variant_code(::EastCascades) = "EC"
nspecies(::EastCascades) = 32
htg_period(::EastCascades) = 10f0   # /CONTRL/ YR=10 (ec/blkdat.f), like WC/PN
mort_ri_scale(::EastCascades) = 0.5f0   # ec/morts.f RI = 0.5·RI (background half-rate), like CR/NE

const EC_RNG_SEED = 55329.0f0           # ec/blkdat.f DATA S0/55329D0/,SS/55329./

# ec/blkdat.f JSP/FIAJSP/PLNJSP (32 species).
const EC_JSP = String[
    "WP","WL","DF","SF","RC","GF","LP","ES","AF","PP","WH","MH","PY","WB","NF","WF","LL","YC","WJ","BM",
    "VN","RA","PB","GC","DG","AS","CW","WO","PL","WI","OS","OH"]
const EC_FIA = String[
    "119","073","202","011","242","017","108","093","019","122","263","264","231","101","022","015","072","042","064","312",
    "324","351","375","431","492","746","747","815","760","920","299","998"]
const EC_PLANTS = String[
    "PIMO3","LAOC","PSME","ABAM","THPL","ABGR","PICO","PIEN","ABLA","PIPO","TSHE","TSME","TABR2","PIAL","ABPR","ABCO","LALY","CANO9","JUOC","ACMA3",
    "ACCI","ALRU2","BEPA","CHCHC4","CONU4","POTR5","POBAT","QUGA4","PRUNU","SALIX","2TN","2TB"]

"""
    ec_grinit!(s)

Set EC's /CONTRL/ + SDI-family flags (ec/grinit.f + blkdat.f): Reineke/Stage SDI, DGSD 1.7, 10-yr native
period, RNG seed 55329, LHTDRG default .FALSE.
"""
function ec_grinit!(s::StandState)
    s.control.year = 10.0f0
    s.control.growth_fint = 10.0f0
    s.control.zeide_sdi = false        # ec/grinit.f:187 LZEIDE=.FALSE. ⇒ Reineke/Stage SDI
    s.control.dg_sd = 1.7f0            # ec/grinit.f:230 DGSD=1.7
    s.control.dg_stddev_bound = 1.7f0
    s.rng.s0 = Float64(EC_RNG_SEED); s.rng.ss = EC_RNG_SEED
    fill!(s.control.ht_drag_sp, false) # ec/grinit.f LHTDRG default .FALSE.
    return s
end
