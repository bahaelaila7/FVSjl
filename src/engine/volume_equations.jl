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

# ---------------------------------------------------------------------------------
# BM (region-6 eastside) VOLEQDEF — faithful port of voleqdef.f R6_EQN eastside branch
# (VAR='BM'). Maps (national-forest number, FIA species) -> NVEL volume-equation id: an INGY
# Flewelling FW2 id (EQNUMI) where the per-forest table assigns one, otherwise region-6 Behre
# 616BEHW<fia>. BM's VEQNNC is FOREST-DEPENDENT (geosub 00/11/12/13/14 and FW2-vs-Behre both
# vary by forest), so a single hardcoded array was only correct for Umatilla (forest 14). This
# reproduces the live per-forest VEQNNC (validated: forest 14 == old BM_VOL_EQ; forest 4 Malheur
# == live dump WL->616BEHW073, PP->I12FW2W122, LP->I12FW2W108).
const _BM_R6_EQNUMI = String[
 "I00FW2W012","I00FW2W017","I00FW2W019","I00FW2W073","I00FW2W093","I00FW2W108",
 "I00FW2W119","I00FW2W122","I00FW2W202","I00FW2W242","I00FW2W260","I00FW2W263",
 "I11FW2W012","I11FW2W017","I11FW2W019","I11FW2W073","I11FW2W093","I11FW2W108",
 "I11FW2W119","I11FW2W122","I11FW2W202","I11FW2W242","I11FW2W260","I11FW2W263",
 "I12FW2W012","I12FW2W017","I12FW2W019","I12FW2W073","I12FW2W093","I12FW2W108",
 "I12FW2W119","I12FW2W122","I12FW2W202","I12FW2W242","I12FW2W260","I12FW2W263",
 "I13FW2W012","I13FW2W017","I13FW2W019","I13FW2W073","I13FW2W093","I13FW2W108",
 "I13FW2W119","I13FW2W122","I13FW2W202","I13FW2W242","I13FW2W260","I13FW2W263",
 "I14FW2W012","I14FW2W017","I14FW2W019","I14FW2W073","I14FW2W093","I14FW2W108",
 "I14FW2W119","I14FW2W122","I14FW2W202","I14FW2W242","I14FW2W260","I14FW2W263",
 "I00FW2W012","I00FW2W202","I00FW2W093","I00FW2W017","I00FW2W108","I00FW2W122",
 "I00FW2W019","I00FW2W073","I00FW2W119","I00FW2W242","I00FW2W260","I13FW2W202",
 "I11FW2W202","I12FW2W202","I14FW2W017","I13FW2W017","I14FW2W122","I12FW2W122",
 "I13FW2W122","I11FW2W122","I11FW2W017","I11FW2W073","I11FW2W242","I12FW2W017",
 "I12FW2W073","I13FW2W073","I14FW2W073"]

const _BM_R6_EQNUMF = String[
 "F00FW2W202","F00FW2W242","F00FW2W263","F01FW2W202","F01FW2W242","F01FW2W263",
 "F02FW2W202","F02FW2W242","F02FW2W263","F03FW2W202","F03FW2W242","F03FW2W263",
 "F04FW2W202","F04FW2W242","F04FW2W263","F05FW2W202","F05FW2W242","F05FW2W263",
 "F06FW2W202","F06FW2W242","F06FW2W263","F07FW2W202","F07FW2W242","F07FW2W263",
 "F08FW2W202","F08FW2W242","F08FW2W263","F03FW2W202","F01FW2W202","F02FW2W202",
 "F00FW2W202","F04FW2W202","F08FW2W202","F07FW2W202","F06FW2W202","F05FW2W202",
 "F00FW2W242","F00FW2W260","F01FW2W242","F01FW2W260","F02FW2W242","F02FW2W260",
 "F03FW2W242","F03FW2W260","F04FW2W242","F04FW2W260","F05FW2W260","F06FW2W260",
 "F07FW2W260","F08FW2W260"]

# FIA(53) sorted list used by R6_EQN's Behre fallback membership check (unknown -> 616BEHW000).
const _BM_R6_FIA = Set{Int}([11,15,17,19,20,21,22,42,64,66,72,73,81,93,98,101,103,106,108,113,
 116,117,119,122,202,211,231,242,263,264,290,299,312,321,351,352,361,375,431,475,478,492,
 500,631,740,746,747,768,815,818,920,998,999])

function _bm_r6_eqn(iforst::Int, idist::Int, fia::Int)::String
    donei = 0; donef = 0
    if iforst == 1                                    # Deschutes
        if fia in (11,15,17,21); donei = 2
        elseif fia == 73; donei = 16
        elseif fia == 108; donei = 18
        elseif fia == 122; donei = 20
        elseif fia == 202; donei = 21
        elseif fia == 81; donei = 22
        elseif fia == 117; donei = 44
        end
    elseif iforst == 2 || iforst == 20                # Fremont
        if fia == 15 || fia == 17; donei = 14
        elseif fia == 81; donei = 9
        elseif fia == 108; donei = 6
        elseif fia == 122; donei = 8
        elseif fia == 202; donei = 2
        end
    elseif iforst == 3                                # Gifford Pinchot
        if fia == 11; donei = 26
        elseif (fia == 263 || fia == 260) && idist == 3; donef = 3
        elseif fia == 202 && idist == 3; donei = 10
        end
    elseif iforst == 4                                # Malheur
        if fia == 108; donei = 30
        elseif fia == 17 || fia == 15; donei = 26
        elseif fia == 122; donei = 32
        elseif fia == 202; donei = 33
        end
    elseif iforst == 6                                # Mt Hood
        if fia == 11; donei = 26
        elseif fia == 17; donei = 38
        elseif fia == 93; donei = 17
        elseif fia == 108; donei = 18
        elseif fia == 122; donei = 32
        elseif fia == 263 || fia == 260; donei = 23
        elseif fia == 22; donei = 38
        elseif fia == 202; donef = (idist == 1 || idist == 6) ? 22 : 16
        end
    elseif iforst == 7                                # Ochoco
        if fia == 17 || fia == 15; donei = 26
        elseif fia == 73; donei = 28
        elseif fia == 122; donei = 32
        elseif fia == 202; donei = 33
        elseif fia == 108; donei = 30
        end
    elseif iforst == 14                               # Umatilla
        if fia == 17 || fia == 15; donei = 38
        elseif fia == 19; donei = 3
        elseif fia == 73; donei = 45
        elseif fia == 93; donei = 5
        elseif fia == 108; donei = 6
        elseif fia == 122; donei = 44
        elseif fia == 202; donei = 38
        end
    elseif iforst == 16                               # Wallowa-Whitman
        if fia == 17 || fia == 15; donei = 14
        elseif fia == 73; donei = 40
        elseif fia == 93; donei = 17
        elseif fia == 108; donei = 6
        elseif fia == 122; donei = 20
        elseif fia == 202; donei = 21
        end
    elseif iforst == 8 || iforst == 17                # Okanogan-Wenatchee
        if fia == 17; donei = 14
        elseif fia == 202; donei = 33
        elseif fia == 108; donei = 30
        elseif fia == 93; donei = 17
        elseif fia == 122 || fia == 73
            donei = (idist == 2 || idist == 3 || idist == 5 || idist == 7) ? 20 : 32
        end
    elseif iforst == 21                               # Colville
        if fia == 17 || fia == 260; donei = 14
        elseif fia == 19; donei = 21
        elseif fia == 73; donei = 16
        elseif fia == 93; donei = 41
        elseif fia == 108; donei = 18
        elseif fia == 119; donei = 7
        elseif fia == 122; donei = 32
        elseif fia == 202; donei = 21
        elseif fia == 242; donei = 22
        elseif fia == 263 || fia == 264; donei = 14
        end
    end
    donei > 0 && return _BM_R6_EQNUMI[donei]
    donef > 0 && return _BM_R6_EQNUMF[donef]
    return fia in _BM_R6_FIA ? "616BEHW" * lpad(string(fia), 3, '0') : "616BEHW000"
end

# NVEL voleqdef.f R4_EQN — Region-4 (Intermountain) default volume-equation assignment, keyed on FORNUM
# (forest number within region) + FIA species code. Faithful transcription incl. the per-species FORNUM
# branches and the binary-search fallback. Shared by the region-4 western variants (CI/UT/TT) whose forests
# span both the INGY (Flewelling FW2) and Matney sub-assignments. Reproduces the live VEQNNC dumps: CI cit01
# (forest 412 = FORNUM 12) and UT utt01 (forest 407 = FORNUM 7).
const _R4_EQNUM = String[
    "400MATW020","300DVEW060","300DVEW060","400DVEW066","400MATW073","400MATW081","300DVEW106","400MATW117",
    "400MATW117","400DVEW133","400DVEW998","400MATW081","400MATW015","401MATW015","400MATW108","400MATW108",
    "400MATW108","400MATW108","400DVEW475","400MATW746","400DVEW998","400MATW108","400MATW108","300DVEW800",
    "300DVEW800","400DVEW998","400DVEW998","I15FW2W017","401MATW015","400MATW015","I15FW2W017","405MATW019",
    "400MATW019","401DVEW065","400DVEW065","I15FW2W093","407FW2W093","400MATW093","401MATW108","400MATW108",
    "I15FW2W122","401MATW122","402MATW122","403MATW122","400MATW122","I15FW2W202","405MATW202","400MATW202",
    "400DVEW064","400DVEW106","407MATW093","401MATW202","407MATW093","401MATW202","401MATW202","407MATW093",
    "407MATW093"]
const _R4_FIA = Int[20,60,64,66,73,81,106,117,119,133,231,242,263,264,299,313,321,322,475,746,747,748,749,800,814,998,999]
const _R4_FIA_IDX = Dict(c => i for (i, c) in enumerate(_R4_FIA))

"NVEL voleqdef.f R4_EQN: (FORNUM, FIA species) → NVEL volume-equation id. `nothing` if unassigned (rare)."
function r4_voleq(fornum::Int, spec::Int)::Union{String,Nothing}
    E = _R4_EQNUM
    f2 = fornum == 2 || fornum == 12 || fornum == 13     # the INGY (Flewelling FW2) forests
    if spec == 15                                        # White fir
        return (f2 || fornum == 6) ? E[28] : (fornum == 9 || fornum == 17) ? E[29] : E[30]
    elseif spec == 17                                    # Grand fir
        return f2 ? E[31] : E[30]
    elseif spec == 19                                    # Subalpine fir
        return fornum == 5 ? E[32] : E[33]
    elseif spec == 64                                    # Utah/western juniper
        return E[49]
    elseif spec == 65                                    # Rocky Mtn juniper (UJ)
        return (fornum == 3 || fornum == 5 || fornum == 15 || fornum == 16) ? E[34] : E[35]
    elseif spec == 93 || spec == 96                      # Engelmann / blue spruce
        return f2 ? E[36] : fornum == 7 ? E[37] : fornum == 8 ? E[56] : E[38]
    elseif spec == 101 || spec == 108 || spec == 113 || spec == 142   # WB/LP/LM/bristlecone
        return (fornum == 9 || fornum == 17) ? E[39] : E[40]
    elseif spec == 106                                   # pinyon
        return (fornum in (1,4,7,8,9,10,17,18,19)) ? E[50] : nothing
    elseif spec == 122                                   # Ponderosa pine
        return f2 ? E[41] : fornum == 1 ? E[42] :
               (fornum in (7,8,10,18,19)) ? E[43] : (fornum == 9 || fornum == 17) ? E[44] : E[45]
    elseif spec == 202                                   # Douglas-fir
        return f2 ? E[46] : fornum == 5 ? E[47] : E[48]
    elseif spec == 998                                   # other hardwoods
        return (fornum in (2,6,12,13,14)) ? E[17] : E[21]
    end
    idx = get(_R4_FIA_IDX, spec, 0)                      # binary-search fallback (forest-independent)
    idx > 0 ? E[idx] : nothing
end

# (KODFOR → (region, FORNUM)) exactly as fvsvol.f decodes it (5-6 digit region-prefixed vs bare 3-digit).
@inline function _kodfor_region_fornum(kodfor::Int)
    if kodfor > 10000
        (kodfor ÷ 10000, kodfor ÷ 100 - (kodfor ÷ 10000) * 100)
    else
        (kodfor ÷ 100, kodfor - (kodfor ÷ 100) * 100)
    end
end

function setup_volume_equations!(s::StandState)
    kodfor = Int(s.plot.user_forest_code)
    iregn  = kodfor ÷ 10000
    iforst = kodfor ÷ 100 - iregn * 100
    intdist = kodfor - (kodfor ÷ 100) * 100
    forst = lpad(string(iforst), 2, '0')
    dist  = lpad(string(intdist), 2, '0')
    cr_tbl = s.variant isa CentralRockies ? _cr_veq_by_forest() : nothing
    r4reg, r4fornum = _kodfor_region_fornum(kodfor)       # for the region-4 R4_EQN western variants (UT/TT)
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
            # TT VOLEQDEF = R4_EQN (region-4, FORNUM-keyed via TT_JFOR[forest_idx]). Forest 405 (F5) selects
            # DF→405MATW202, AF→405MATW019 (confirmed live). The earlier apparent ".sum regression" from this
            # was NOT the equation — it was a latent BROKEN-TOP volume bug (r4_topkill, now fixed in
            # compute_volumes_tt!/ci!/ut!): a top-killed AF over-volumed with the correct 405 eq. With both
            # fixes, forest-405 TCuFt is bit-exact (residual = a small separate r4vol merch tail). The 405
            # SELECTION is validated (r4_voleq 555/555 vs authoritative R4_EQN).
            fidx_tt = Int(s.plot.forest_idx)
            fnum_tt = (1 <= fidx_tt <= length(TT_JFOR)) ? TT_JFOR[fidx_tt] % 100 : 3
            veq_tt = r4_voleq(fnum_tt, ifia)
            s.species.vol_eq[sp] = veq_tt !== nothing ? veq_tt : (sp <= length(TT_VOL_EQ) ? TT_VOL_EQ[sp] : "          ")
        elseif s.variant isa Utah
            # UT VOLEQDEF = R4_EQN (region-4, FORNUM-keyed). jl previously used ONLY the forest-407/utt01 table,
            # so WB/LP/LM/PP/WF/ES/BS volume was wrong on other UT forests (e.g. F9/F17 use 401MATW108/403MATW122
            # /401MATW015, F9/F10/F17 use 400MATW093 not 407FW2W093). Key on the ut_forkod!-CORRECTED forest —
            # UT_JFOR[forest_idx] — NOT the raw user_forest_code, because forkod remaps (Humboldt 409→Fishlake 408
            # ⇒ FORNUM 8, matching live's "FOREST-LOCATION CODE 408"). UT_VOL_EQ (= r4_voleq at F7) is the fallback.
            fidx = Int(s.plot.forest_idx)
            fnum = (1 <= fidx <= length(UT_JFOR)) ? UT_JFOR[fidx] % 100 : 7
            veq = r4_voleq(fnum, ifia)
            s.species.vol_eq[sp] = veq !== nothing ? veq :
                (sp <= length(UT_VOL_EQ) ? UT_VOL_EQ[sp] : "           ")
        elseif s.variant isa CentralIdaho
            # CI VOLEQDEF (NVEL voleqdef.f R4_EQN, FORNUM-keyed): 400MATW (R4 Matney) + I15FW2W (Flewelling,
            # DF/GF/ES/PP on INGY forests 402/412/413) + 400DVEW (woodland, PY/WJ/MC/CW). Forest-DEPENDENT: the
            # Salmon-Challis/Sawtooth forests (406/414 = forest_idx 3/6) use Matney for DF/GF/ES/PP — see
            # ci_vol_eq_table. forest_idx is set by ci_forkod! in site_setup! (runs before volume setup).
            citbl = ci_vol_eq_table(Int(s.plot.forest_idx))
            s.species.vol_eq[sp] = sp <= length(citbl) ? citbl[sp] : "           "
        elseif s.variant isa EasternMontana
            # EM VOLEQDEF (em/sitset.f VEQNNC, dumped from live): conifers = FW2 (volume-FIA ≠ species-FIA,
            # e.g. LM→073); non-conifers RM/GA/AS/CW/BA/PW/NC/PB/OH = DVEW woodland (region 1/2).
            s.species.vol_eq[sp] = sp <= length(EM_VOL_EQ) ? EM_VOL_EQ[sp] : "           "
        elseif s.variant isa BlueMountains
            # BM VOLEQDEF (voleqdef.f R6_EQN eastside) — FOREST-DEPENDENT: (forest, FIA)->veq. INGY FW2
            # (geosub varies by forest) where the forest table assigns one, else region-6 Behre 616BEHW.
            # Key off the forkod-REMAPPED forest (619 Whitman -> 616 W-W, 8117 -> 614), as live does.
            bmkf = bm_kodfor_remap(kodfor)
            s.species.vol_eq[sp] = _bm_r6_eqn(bmkf % 100, 0, ifia)
        elseif s.variant isa WestCascades
            # WC VOLEQDEF (voleqdef.f R6_EQN WESTSIDE branch, VAR='WC') — FORNUM=KODFOR%100. Willamette
            # (618→FORNUM 18) validated vs FVSwc_clean; westside Flewelling (DF/WH) + INGY (GF/NF/IC) + Behre.
            s.species.vol_eq[sp] = _wc_r6_eqn(kodfor % 100, ifia)
        elseif s.variant isa PacificNorthwest
            s.species.vol_eq[sp] = _pn_r6_eqn(kodfor % 100, ifia)
        else
            s.species.vol_eq[sp] = (iregn == 8 && ifia > 0) ? _r8_ceqn(forst, dist, ifia) : "           "
        end
    end
    return s
end
