# =============================================================================
# westcascades.jl — the WC (West Cascades) variant singleton + registration.
#
# Ported from: wc/*.f (FVS "WEST CASCADES"). WC is the westside R6 "Prognosis" family's
# STRUCTURAL ANCHOR (roadmap docs/WESTERN_UNPORTED_VARIANTS_ROADMAP.md): PN/EC/CA/WS borrow
# its htgf/regent/cratet/dgdriv. It is a western Wykoff-DDS variant, template = BM
# (formcl/htdbh/sichg/ecocls identical structure).
# Distinguishing infra (MEASURED, wc/grinit.f + wc/blkdat.f):
#   • VARACD='WC', MAXSP=39 (slot 6 blank), IFINT=10 / IFINTH=5.
#   • LZEIDE=.FALSE. ⇒ Reineke/Stage SDI (NOT Zeide — unlike NC/CI/UT/TT). ⇒ reuse the base self-thin.
#   • DGSD=1.7 ; LHTDRG default .FALSE. ; RNG seed 55329 (wc/blkdat.f S0/SS).
#   • wc/dgf.f DDS is the standard western Wykoff LN(DDS) with 19 SPECIES GROUPS (MAPSPC 39→19),
#     a single DEFAULT branch + RA (group 13) and RW (ISPC 17) special equations.
# Species JSP (39; slot 6 blank): SF WF GF AF RF __ NF YC IC ES LP JP SP WP PP DF RW RC WH MH BM
#     RA WA PB GC AS CW WO WJ LL WB KP PY DG HT CH WI __ OT
# Oracle: /workspace/.wcwork/FVSwc_clean (relink bin/FVSwc_buildDir/*.o + isoc23 shim; relink_wc.sh).
# Canonical test stand: wct01 (tests/FVSwc/wct01.{key,tre}; stand 248112, forest JFOR=6).
#
# PORT STATUS — chunk 0 foundation + chunk 3 large-tree DGF (VALIDATED at formula level vs
# FVSwc_clean DEBUG DGF: DEFAULT DDS 27/27 trees bit-exact, DGCON reconstruction bit-exact).
# Remaining chunks (species CSV / site / density / height / regent / crown / mortality / volume)
# are TODO — dispatching the un-ported hooks on WestCascades() errors loudly (doctrine #5).
# =============================================================================

"""
    WestCascades <: AbstractVariant

The FVS West Cascades variant (VARACD "WC", MAXSP = 39) — westside R6 Wykoff-DDS anchor,
Reineke SDI (BM template). PORT IN PROGRESS (chunk 0 foundation + chunk 3 DGF validated).
"""
struct WestCascades <: AbstractVariant end

variant_code(::WestCascades) = "WC"
nspecies(::WestCascades) = 39
htg_period(::WestCascades) = 5f0    # wc IFINTH=5 (height growth is a 5-yr rise); DG model is 5-yr (TDDS/2)

const WC_RNG_SEED = 55329.0f0       # wc/blkdat.f DATA S0/55329D0/,SS/55329./

# wc/blkdat.f JSP/FIAJSP/PLNJSP (39 species; slot 6 blank "__"/"999" placeholder).
const WC_JSP = String[
    "SF","WF","GF","AF","RF","__","NF","YC","IC","ES","LP","JP","SP","WP","PP","DF","RW","RC","WH","MH",
    "BM","RA","WA","PB","GC","AS","CW","WO","WJ","LL","WB","KP","PY","DG","HT","CH","WI","__","OT"]
const WC_FIA = String[
    "011","015","017","019","020","   ","022","042","081","093","108","116","117","119","122","202","211","242","263","264",
    "312","351","352","375","431","746","747","815","064","072","101","103","231","492","500","768","920","   ","999"]
const WC_PLANTS = String[
    "ABAM","ABCO","ABGR","ABLA","ABMA","","ABPR","CANO9","CADE27","PIEN","PICO","PIJE","PILA","PIMO3","PIPO","PSME","SESE3","THPL","TSHE","TSME",
    "ACMA3","ALRU2","ALRH2","BEPA","CHCHC4","POTR5","POBAT","QUGA4","JUOC","LALY","PIAL","PIAT","TABR2","CONU4","CRATA","PREM","SALIX","","2TREE"]

"""
    wc_grinit!(s)

Set WC's /CONTRL/ + SDI-family flags (wc/grinit.f + blkdat.f). Called from the block-data
init once a WC species table lands (chunk 1). Kept separate so the DGF chunk can be wired
and validated ahead of the full species-coefficient loader.
"""
function wc_grinit!(s::StandState)
    s.control.year = 5.0f0             # cycle default (wct01 live runs 5-yr; WC DG models are 5-yr)
    s.control.growth_fint = 5.0f0
    s.control.zeide_sdi = false        # wc/grinit.f LZEIDE=.FALSE. ⇒ Reineke/Stage SDI (base self-thin)
    s.control.dg_sd = 1.7f0            # wc/grinit.f DGSD=1.7
    s.control.dg_stddev_bound = 1.7f0  # set BOTH (the BM/CI DGSD field-disconnect lesson)
    s.rng.s0 = Float64(WC_RNG_SEED); s.rng.ss = WC_RNG_SEED
    fill!(s.control.ht_drag_sp, false) # wc/grinit.f:102 LHTDRG default .FALSE.
    return s
end
