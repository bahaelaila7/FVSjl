# =============================================================================
# southcentraloregon.jl — the SO (SORNEC / South-Central Oregon & NE California) variant singleton.
#
# Ported from: so/*.f (FVS "SO", grinit.f VARACD='SO', blkdat.f "SORNEC-33 variant"). SO is a westside
# Wykoff-DDS variant on the SAME shared engine as WC/PN/EC/CA — 10-yr native period (FINT=10, IFINT=10/
# IFINTH=5), RNG seed 55329 (same as CA), DGSD=2.0. MAXSP=33 with a PER-SPECIES dgf (NOT group-compressed
# like CA's 13): DGLD(33)/DGDBAL(33)/DGCR/DGCRSQ(33), DGFOR(6,33) location classes, DGDS(4,33) diam-sq.
#
# ⚠ SDI type is FOREST-DEPENDENT (so/grinit.f:126 LZEIDE=.TRUE., CALCSDI=' '; so/sitset.f:81 resets
#   LZEIDE=.FALSE. ONLY for R6 forests `IF(IFOR.LE.3 .OR. IFOR.EQ.10)`). So R6 forests → Reineke;
#   interior/other forests → Zeide. The IFOR-dependent reset belongs in the SO sitset chunk (chunk 2);
#   grinit here sets the LZEIDE=.TRUE. default. sot01 is R6 (→Reineke, matches its .sum SDI).
# Species JSP (33): WP SP DF WF MH IC LP ES SH PP WJ GF AF SF NF WB WL RC WH PY WA RA BM AS CW CH WO WI
#   GC MC MB OS OH
# Oracle: /workspace/.sowork/FVSso_clean (relink_so.sh: bin/FVSso_buildDir/*.o + isoc23 shim). Canonical
#   stand: sot01 (tests/FVSso/). Growth cols bit-exact vs sot01.sum.save (vol = NVEL lib-version delta).
# =============================================================================

"""
    SouthCentralOregon <: AbstractVariant

The FVS South-Central Oregon / NE California variant (VARACD "SO", SORNEC-33, MAXSP = 33) — a westside
Wykoff-DDS variant (DGSD 2.0, forest-dependent Zeide/Reineke SDI) with per-species diameter-growth eqs.
"""
struct SouthCentralOregon <: AbstractVariant end

variant_code(::SouthCentralOregon) = "SO"
nspecies(::SouthCentralOregon) = 33
htg_period(::SouthCentralOregon) = 10f0   # /CONTRL/ YR=10 (so/grinit.f FINT=10), like WC/PN/EC/CA
mort_ri_scale(::SouthCentralOregon) = 0.5f0   # background half-rate (so/morts.f), like the westside cluster

const SO_RNG_SEED = 55329.0f0            # so/blkdat.f DATA S0/55329D0/,SS/55329./ (same as CA)

# so/blkdat.f JSP/FIAJSP/PLNJSP (33 species, SORNEC order).
const SO_JSP = String[
    "WP","SP","DF","WF","MH","IC","LP","ES","SH","PP","WJ","GF","AF","SF","NF","WB","WL",
    "RC","WH","PY","WA","RA","BM","AS","CW","CH","WO","WI","GC","MC","MB","OS","OH"]
const SO_FIA = String[
    "119","117","202","015","264","081","108","093","021","122","064","017","019","011","022","101","073",
    "242","263","231","352","351","312","746","747","768","815","920","431","475","478","299","998"]
const SO_PLANTS = String[
    "PIMO3","PILA","PSME","ABCO","TSME","CADE27","PICO","PIEN","ABSH","PIPO","JUOC","ABGR","ABLA","ABAM","ABPR","PIAL","LAOC",
    "THPL","TSHE","TABR2","ALRH2","ALRU2","ACMA3","POTR5","POBAT","PREM","QUGA4","SALIX","CHCHC4","CELE3","CEMOG","2TN","2TB"]

"""
    so_grinit!(s)

Set SO's /CONTRL/ + SDI-family flags (so/grinit.f + blkdat.f): 10-yr native period, DGSD 2.0, RNG seed
55329, LHTDRG default .FALSE. LZEIDE=.TRUE. here (so/grinit.f:126); the SO sitset (chunk 2) resets it to
Reineke for R6 forests (IFOR≤3 or =10) per so/sitset.f:81.
"""
function so_grinit!(s::StandState)
    s.control.year = 10.0f0
    s.control.growth_fint = 10.0f0
    s.control.zeide_sdi = true          # so/grinit.f:126 LZEIDE=.TRUE. (chunk-2 sitset resets to Reineke for R6 IFOR)
    s.control.dg_sd = 2.0f0             # so/grinit.f:166 DGSD=2.0
    s.control.dg_stddev_bound = 2.0f0
    s.rng.s0 = Float64(SO_RNG_SEED); s.rng.ss = SO_RNG_SEED
    fill!(s.control.ht_drag_sp, false)  # so/grinit.f:102 LHTDRG default .FALSE.
    # so/grinit.f:222-223 default location TLAT=42/TLONG=121 (sot01 carries no LOCATION keyword). Feeds the
    # Hopkins bioclimatic index used by the Bechtold-2004 crown-width eqns (GC/MC/MB); without it jl's lat/lon=0
    # gives HI=-323.6 (clamped) vs the oracle's ~-4.4 ⇒ those 3 species' FVS_TreeList CrWidth were off 1-2.6.
    # Only default when unset (an explicit LOCATION / DB lat-lon wins, as initre.f overrides grinit's default).
    s.plot.latitude  == 0f0 && (s.plot.latitude  = 42f0)
    s.plot.longitude == 0f0 && (s.plot.longitude = 121f0)
    return s
end
