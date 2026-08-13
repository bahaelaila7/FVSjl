# =============================================================================
# pacificnorthwest.jl — the PN (Pacific Northwest) variant singleton + registration.
#
# Ported from: pn/*.f (FVS "PN"). PN is a NEAR-CLONE of West Cascades (WC) on the SAME shared
# westside Wykoff-DDS engine — MAXSP=39, Reineke/Stage SDI (LZEIDE=.FALSE.), DGSD=1.7, IFINT=10/
# IFINTH=5, RNG seed 55329 — differing only in the coefficient DATA (dgf/htdbh/crown/sichg/sitset/
# ecocls/habtyp/formcl/bratio/ccfcal/htcalc) and the species table (slot 6 = SS/Sitka spruce vs WC's
# blank). The PN functions mirror the WC ones with PN_* coefficients so WC stays provably inert.
# Species JSP (39): SF WF GF AF RF SS NF YC IC ES LP JP SP WP PP DF RW RC WH MH BM RA WA PB GC AS CW WO
#   WJ LL WB KP PY DG HT CH WI __ OT
# Oracle: /workspace/.pnwork/FVSpn_clean (relink bin/FVSpn_buildDir/*.o + isoc23 shim; relink_pn.sh);
#   instrumentable /workspace/.pnwork/FVSpn_g16. Canonical stand: pnt01 (= same trees as wct01, S248112).
# =============================================================================

"""
    PacificNorthwest <: AbstractVariant

The FVS Pacific Northwest variant (VARACD "PN", MAXSP = 39) — a westside R6 Wykoff-DDS near-clone of
West Cascades (Reineke SDI, DGSD 1.7). PORT IN PROGRESS (coefficient-swap on the shared WC engine).
"""
struct PacificNorthwest <: AbstractVariant end

variant_code(::PacificNorthwest) = "PN"
nspecies(::PacificNorthwest) = 39
htg_period(::PacificNorthwest) = 10f0   # /CONTRL/ YR=10 (pn/blkdat.f), like WC

const PN_RNG_SEED = 55329.0f0           # pn/blkdat.f DATA S0/55329D0/,SS/55329./

# pn/blkdat.f JSP/FIAJSP/PLNJSP (39 species; slot 6 = SS/Sitka, slot 38 blank).
const PN_JSP = String[
    "SF","WF","GF","AF","RF","SS","NF","YC","IC","ES","LP","JP","SP","WP","PP","DF","RW","RC","WH","MH",
    "BM","RA","WA","PB","GC","AS","CW","WO","WJ","LL","WB","KP","PY","DG","HT","CH","WI","__","OT"]
const PN_FIA = String[
    "011","015","017","019","020","098","022","042","081","093","108","116","117","119","122","202","211","242","263","264",
    "312","351","352","375","431","746","747","815","064","072","101","103","231","492","500","768","920","   ","999"]
const PN_PLANTS = String[
    "ABAM","ABCO","ABGR","ABLA","ABMA","PISI","ABPR","CANO9","CADE27","PIEN","PICO","PIJE","PILA","PIMO3","PIPO","PSME","SESE3","THPL","TSHE","TSME",
    "ACMA3","ALRU2","ALRH2","BEPA","CHCHC4","POTR5","POBAT","QUGA4","JUOC","LALY","PIAL","PIAT","TABR2","CONU4","CRATA","PREM","SALIX","","2TREE"]

"""
    pn_grinit!(s)

Set PN's /CONTRL/ + SDI-family flags (pn/grinit.f + blkdat.f), identical to WC: Reineke/Stage SDI,
DGSD 1.7, 10-yr native period, RNG seed 55329, LHTDRG default .FALSE.
"""
function pn_grinit!(s::StandState)
    s.control.year = 10.0f0
    s.control.growth_fint = 10.0f0
    s.control.zeide_sdi = false        # pn/grinit.f LZEIDE=.FALSE. ⇒ Reineke/Stage SDI (base self-thin)
    s.control.dg_sd = 1.7f0            # pn/grinit.f DGSD=1.7
    s.control.dg_stddev_bound = 1.7f0
    s.rng.s0 = Float64(PN_RNG_SEED); s.rng.ss = PN_RNG_SEED
    fill!(s.control.ht_drag_sp, false) # pn/grinit.f LHTDRG default .FALSE.
    return s
end
