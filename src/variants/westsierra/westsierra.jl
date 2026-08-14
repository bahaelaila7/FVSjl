# =============================================================================
# westsierra.jl — the WS (Western Sierra Nevada) variant singleton.
#
# Ported from: ws/*.f (FVS "WS", grinit.f VARACD='WS'). WS is a westside/California Wykoff-DDS variant on
# the SAME shared engine as WC/PN/EC/CA/SO — 10-yr native period (FINT=10, IFINT=10/IFINTH=5), RNG seed
# 55329 (same as CA/SO), DGSD=2.0. MAXSP=43 California species (redwood/sequoia/tanoak/oaks). WS is a
# close CentralCalifornia SIBLING (R5 California) — expect heavy reuse of the CA machinery in later chunks.
#
# ⚠ SDI type = ZEIDE (ws/grinit.f:209 LZEIDE=.TRUE., CALCSDI=' '; NO IFOR-dependent reset in ws/sitset.f —
#   unlike SO which resets to Reineke for R6 forests). So WS stays Zeide throughout (unless CALCSDI keyword).
# Species JSP (43): SP DF WF GS IC JP RF PP LP WB WP PM SF KP FP CP LM MP GP WE GB BD RW MH WJ UJ CJ LO CY
#   BL BO VO IO TO GC AS CL MA DG BM MC OS OH  (GS=giant sequoia, RW=coast redwood, TO=tanoak, oaks).
# Oracle: /workspace/.wswork/FVSws_clean (relink_ws.sh: bin/FVSws_buildDir/*.o + isoc23 shim); FVSws_g16
#   (build_g16.sh) instrumentable. Canonical stand: wst01 (tests/FVSws/). Growth cols BIT-EXACT vs
#   wst01.sum.save (vol = NVEL lib-version delta). Harness live in /workspace/.wswork/run.
# =============================================================================

"""
    WestSierra <: AbstractVariant

The FVS Western Sierra Nevada variant (VARACD "WS", MAXSP = 43) — a California R5 westside Wykoff-DDS
variant (DGSD 2.0, Zeide SDI) with per-species diameter-growth eqs; a CentralCalifornia sibling.
"""
struct WestSierra <: AbstractVariant end

variant_code(::WestSierra) = "WS"
nspecies(::WestSierra) = 43
htg_period(::WestSierra) = 10f0   # /CONTRL/ YR=10 (ws/grinit.f FINT=10), like WC/PN/EC/CA/SO
mort_ri_scale(::WestSierra) = 0.5f0   # background half-rate (ws/morts.f:569 RI = 0.5·RI), like the westside cluster

const WS_RNG_SEED = 55329.0f0            # ws/blkdat.f DATA S0/55329D0/,SS/55329./ (same as CA/SO)

# ws/blkdat.f JSP/FIAJSP/PLNJSP (43 species, WS order).
const WS_JSP = String[
    "SP","DF","WF","GS","IC","JP","RF","PP","LP","WB","WP","PM","SF","KP","FP","CP","LM","MP","GP","WE",
    "GB","BD","RW","MH","WJ","UJ","CJ","LO","CY","BL","BO","VO","IO","TO","GC","AS","CL","MA","DG","BM",
    "MC","OS","OH"]
const WS_FIA = String[
    "117","202","015","212","081","116","020","122","108","101","119","133","011","103","104","109","113",
    "124","127","137","142","201","211","264","064","065","062","801","805","807","818","821","839","631",
    "431","746","981","361","492","312","475","299","998"]
const WS_PLANTS = String[
    "PILA","PSME","ABCO","SEGI2","CADE27","PIJE","ABMA","PIPO","PICO","PIAL","PIMO3","PIMO","ABAM","PIAT",
    "PIBA","PICO3","PIFL2","PIRA2","PISA2","PIWA","PILO","PSMA","SESE3","TSME","JUOC","JUOS","JUCA7","QUAG",
    "QUCH2","QUDO","QUKE","QULO","QUWI2","LIDE3","CHCHC4","POTR5","UMCA","ARME","CONU4","ACMA3","CELE3",
    "2TN","2TB"]

"""
    ws_grinit!(s)

Set WS's /CONTRL/ + SDI-family flags (ws/grinit.f + blkdat.f): 10-yr native period, DGSD 2.0, RNG seed
55329, LHTDRG default .FALSE., LZEIDE=.TRUE. (Zeide SDI — no IFOR reset, unlike SO).
"""
function ws_grinit!(s::StandState)
    s.control.year = 10.0f0
    s.control.growth_fint = 10.0f0
    s.control.zeide_sdi = true          # ws/grinit.f:209 LZEIDE=.TRUE. (no sitset reset — WS stays Zeide)
    s.control.dg_sd = 2.0f0             # ws/grinit.f:251 DGSD=2.0
    s.control.dg_stddev_bound = 2.0f0
    s.rng.s0 = Float64(WS_RNG_SEED); s.rng.ss = WS_RNG_SEED
    fill!(s.control.ht_drag_sp, false)  # ws/grinit.f LHTDRG default .FALSE.
    return s
end
