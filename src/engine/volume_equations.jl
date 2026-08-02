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

# CR default NVEL volume-equation ids (VOLEQDEF for VAR='CR'), verified vs live crt01.out.
# NOTE: this is the forest-303 (crt01) assignment; VOLEQDEF is region/forest-keyed, so a full
# voleqdef.f port is needed for arbitrary CR forests. cuft==bdft eq for CR.
const _CR_VOLEQ = String[
    "300DVEW093",
    "NVB0000015",
    "300FW2W202",
    "300DVEW093",
    "NVB0000015",
    "301DVEW015",
    "301DVEW015",
    "301DVEW015",
    "300DVEW113",
    "300DVEW113",
    "301DVEW202",
    "300DVEW106",
    "300FW2W122",
    "300DVEW113",
    "NVBM240119",
    "300DVEW060",
    "300DVEW093",
    "NVBM330093",
    "300DVEW093",
    "NVB0000746",
    "300DVEW999",
    "300DVEW999",
    "300DVEW800",
    "300DVEW800",
    "300DVEW800",
    "300DVEW800",
    "300DVEW800",
    "300DVEW999",
    "300DVEW060",
    "300DVEW060",
    "300DVEW060",
    "300DVEW060",
    "300DVEW106",
    "300DVEW106",
    "300DVEW106",
    "300DVEW122",
    "300DVEW060",
    "300DVEW999"]

# CR forest-keyed VOLEQDEF table (data/centralrockies/volume_equations_by_forest.csv), dumped from the
# module-free NVEL voleqdef.o for VAR='CR' × all 29 CR forests × 38 species. Keyed (KODFOR, FIA)→veq.
const _CR_VEQ_BY_FOREST = Ref{Union{Nothing,Dict{Tuple{Int,Int},String}}}(nothing)
function _cr_veq_by_forest()
    _CR_VEQ_BY_FOREST[] !== nothing && return _CR_VEQ_BY_FOREST[]
    d = Dict{Tuple{Int,Int},String}()
    path = joinpath(@__DIR__, "..", "..", "data", "centralrockies", "volume_equations_by_forest.csv")
    for (li, line) in enumerate(eachline(path))
        li == 1 && continue
        f = split(line, ',')
        length(f) < 3 && continue
        kf = tryparse(Int, strip(f[1])); fia = tryparse(Int, strip(f[2]))
        (kf === nothing || fia === nothing) && continue
        veq = strip(f[3])
        (isempty(veq) || all(==('0'), veq)) && continue      # skip the "000…" no-eq rows
        d[(kf, fia)] = veq
    end
    _CR_VEQ_BY_FOREST[] = d
    return d
end

# KT per-species VOLUME FIA code (kt/sitset.f FIAJSP), validated vs live VEQNNC. sp11 OT = 260 (not code_fia 999).
const KT_VOL_FIA = Int[119, 73, 202, 17, 260, 242, 108, 93, 19, 122, 260]

function setup_volume_equations!(s::StandState)
    kodfor = Int(s.plot.user_forest_code)
    iregn  = kodfor ÷ 10000
    iforst = kodfor ÷ 100 - iregn * 100
    intdist = kodfor - (kodfor ÷ 100) * 100
    forst = lpad(string(iforst), 2, '0')
    dist  = lpad(string(intdist), 2, '0')
    cr_tbl = s.variant isa CentralRockies ? _cr_veq_by_forest() : nothing
    @inbounds for sp in 1:MAXSP
        ifia = something(tryparse(Int, strip(s.coef.code_fia[sp])), 0)
        if s.variant isa CentralRockies
            # Forest-keyed VOLEQDEF (KODFOR, FIA); fall back to the forest-303/crt01 default table.
            veq = get(cr_tbl, (kodfor, ifia), nothing)
            s.species.vol_eq[sp] = veq !== nothing ? veq :
                (sp <= length(_CR_VOLEQ) ? _CR_VOLEQ[sp] : "           ")
        elseif s.variant isa Kootenai
            # KT VOLEQDEF (kt/sitset.f): Region-1 Flewelling FW2, geocode "I00", per-species VOLUME FIA code
            # (FIAJSP, NOT code_fia — sp11 OT is 260 not 999). Validated vs live VEQNNC dump for all 11 species.
            vfia = KT_VOL_FIA[sp]
            s.species.vol_eq[sp] = "I00FW2W" * lpad(string(vfia), 3, '0')
        elseif s.variant isa InlandEmpire
            # IE VOLEQDEF (ie/sitset.f, VAR='IE'): FW2 (sp1-14,23) + DVE/Behre (sp15-22). Strings dumped from
            # live (forest 118); forest-keyed VOLEQDEF port needed for arbitrary IE forests (see volume.jl).
            s.species.vol_eq[sp] = sp <= length(IE_VOL_EQ) ? IE_VOL_EQ[sp] : "           "
        elseif s.variant isa Teton
            s.species.vol_eq[sp] = sp <= length(TT_VOL_EQ) ? TT_VOL_EQ[sp] : "          "
        elseif s.variant isa EasternMontana
            # EM VOLEQDEF (em/sitset.f, VAR='EM', IREGN=1): Region-1 Flewelling FW2, "I00FW2W<FIAJSP>" — the
            # conifers (emt01) use FW2. Hardwoods (GA/AS/CW/…) may use DVE/Behre (not emt01) — a forest-keyed
            # VOLEQDEF port covers those; here FW2 for all via the species FIA (code_fia = FIAJSP for EM).
            s.species.vol_eq[sp] = ifia > 0 ? "I00FW2W" * lpad(string(ifia), 3, '0') : "           "
        else
            s.species.vol_eq[sp] = (iregn == 8 && ifia > 0) ? _r8_ceqn(forst, dist, ifia) : "           "
        end
    end
    return s
end
