# =============================================================================
# volume_equations.jl — per-species NVEL volume-equation identifiers (VOLEQDEF)
#
# Ported from: base/volstubs.jl (_r8_ceqn / VOLEQDEF) + sn/sitset.jl (forst/dist
# from KODFOR). Each SN species gets an R8 New Clark equation string
# "8<geoa>1CLKE<vol-species-code>" (e.g. sugar maple FIA 318 → "841CLKE318"); the
# geographic-area digit `geoa` comes from the national-forest number + district,
# and the volume species code from the SNFIA→SNSP crosswalk. The string drives
# `_R8CLARK_VOL` (engine/r8clark_vol.jl). For SN, METHC=METHB=6 → all species use
# this path, and the cubic and board equations are identical.
# =============================================================================

# FIA species code → volume-equation species code (binary-search SNFIA, take SNSP).
const _VOL_SNFIA = Int32[
   10, 57, 90,100,107,110,111,115,121,123,
  126,128,129,130,131,132,197,221,222,260,
  261,299,300,311,313,314,316,317,318,330,
  370,372,391,400,404,450,460,471,491,500,
  521,531,540,541,543,544,545,546,550,552,
  555,580,591,601,602,611,621,650,651,652,
  653,654,660,680,691,693,694,701,711,721,
  731,740,742,743,762,800,802,804,806,812,
  813,817,819,820,822,823,824,825,826,827,
  828,830,831,832,833,834,835,837,838,901,
  920,930,931,950,970,971,972,975,998,999]

const _VOL_SNSP = String[
  "261","100","115","100","107","110","111","115","121","123",
  "126","128","129","132","131","132","197","221","222","261",
  "261","132","300","500","313","314","316","317","318","330",
  "370","370","370","400","404","300","460","300","300","500",
  "521","531","541","541","300","544","545","546","550","500",
  "300","580","300","601","602","611","621","650","651","652",
  "653","300","300","300","691","693","694","500","711","300",
  "731","300","742","300","762","800","802","804","806","812",
  "813","817","800","820","822","823","800","825","826","827",
  "828","830","831","832","833","834","835","837","835","901",
  "920","930","300","950","970","970","970","970","300","300"]

# geographic-area digit from national-forest number + district (volstubs.f _r8_ceqn).
function _vol_geoa(fornum::Int, distnum::Int)
    fornum == 1            ? (distnum == 3 ? '1' : '4') :
    fornum in (2,4,8,60)   ? '3' :
    fornum == 3            ? (distnum == 8 ? '2' : '3') :
    fornum in (5,36)       ? '1' :
    fornum in (6,13)       ? '5' :
    fornum == 7            ? (distnum == 6 ? '7' : (distnum in (7,17) ? '4' : '5')) :
    fornum == 9            ? '6' :
    fornum == 10           ? (distnum == 7 ? '7' : '6') :
    fornum == 11           ? (distnum == 3 ? '1' : (distnum == 10 ? '2' : '3')) :
    fornum == 12           ? (distnum == 2 ? '3' : (distnum == 5 ? '1' : '2')) :
    '9'
end

"R8 cubic-equation string for an FIA species code (volstubs.f _r8_ceqn)."
function _r8_ceqn(forst::AbstractString, dist::AbstractString, spec::Int)
    fornum  = something(tryparse(Int, strip(forst)), 0)
    distnum = something(tryparse(Int, strip(dist)), 0)
    geoa = _vol_geoa(fornum, distnum)
    # binary search _VOL_SNFIA for spec
    first_, last_ = 1, 110; done = 0
    while true
        half = (last_ - first_ + 1) ÷ 2 + first_
        if _VOL_SNFIA[half] == spec
            done = half; break
        elseif first_ == last_
            done = spec < 300 ? 22 : 110; break
        elseif _VOL_SNFIA[half] < spec
            first_ = half
        else
            last_ = max(half - 1, first_)
        end
    end
    return "8" * string(geoa) * "1CLKE" * _VOL_SNSP[done]
end

"""
    setup_volume_equations!(state)

VOLEQDEF: assign each species its NVEL volume-equation id (`species.vol_eq`) from
the stand's national-forest code (KODFOR) and the species FIA code. Region-8 only.
"""
function setup_volume_equations!(s::StandState)
    kodfor = Int(s.plot.user_forest_code)
    iregn  = kodfor ÷ 10000
    iforst = kodfor ÷ 100 - iregn * 100
    intdist = kodfor - (kodfor ÷ 100) * 100
    forst = lpad(string(iforst), 2, '0')
    dist  = lpad(string(intdist), 2, '0')
    @inbounds for sp in 1:MAXSP
        ifia = something(tryparse(Int, strip(s.coef.code_fia[sp])), 0)
        s.species.vol_eq[sp] = (iregn == 8 && ifia > 0) ? _r8_ceqn(forst, dist, ifia) : "           "
    end
    return s
end
