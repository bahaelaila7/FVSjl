# =============================================================================
# site_index.jl (inlandempire) — IE habtyp/forkod/sitset (chunk 2).
#
# Produces the three inputs the DG chunk (dgf!/ie_dgcons!) reads:
#   ITYPE  (1..30) -> p.habitat_input   (ie/habtyp.f: JTYPE(95) bracket search -> KTYPE)
#   IFOR   (1..11) -> p.forest_idx       (ie/forkod.f: JFOR(15) match + remap)
#   SITEAR (per-sp) -> p.sp_site_index   (ie/sitset.f: MAPSIT(ITYPE,sp) default when 0)
# Also sets IGL (geo-location) -> p.geo_location and the SDImax defaults (BAMAXA/SDIDEF).
#
# habtyp NUMERIC path only (STDINFO habitat code): the PVREF plant-association-reference
# crosswalk + HBDECD Region-6 plant-community decode are deferred (fire only for CPVREF /
# R6 inputs). Validated: kod=570 -> ITYPE=17, forest 118 -> IFOR=11 (live FVSie iet01).
# =============================================================================

# ie/habtyp.f + ie/blkdat.f JTYPE(95): habitat-code brackets (identical to KT_JTYPE).
const IE_JTYPE = Int32[10,100,110,130,140,160,170,180,190,200,210,220,230,250,260,280,290,310,320,330,
                       340,350,360,370,380,400,410,420,430,440,450,460,470,480,500,501,502,505,506,510,
                       515,516,520,529,530,540,545,550,555,560,565,570,575,579,590,600,610,620,630,635,
                       640,650,660,670,675,680,685,690,700,701,710,720,730,740,750,770,780,790,800,810,
                       820,830,840,850,860,870,890,900,910,920,925,930,940,950,999]
# ie/habtyp.f KTYPE(95): bracket index -> ITYPE (1..30).
const IE_KTYPE = Int32[1,1,1,1,1, 2,2,2,2, 4, 1,1,1, 3,4,5,6,7,8,9,8,8,9,7,3, 10,10,10, 4,11,20,29,11,
                       11,13,14,17, 12,12,12,12, 13,13,13, 14,15,14,16,14,16, 17,17,17, 24,12,24,18,
                       19,21,19,20,20,21,22,19,23,19,24,27,25,25,26,27,22,24,27,24,
                       27,28,28,29,28,28, 29,29,29,29, 27,9,20,24,21,27,24,30]
# ie/habtyp.f MTYPE(30): representative habitat code per ITYPE (informational; KODTYP output).
const IE_MTYPE = Int32[130,170,250,260,280,290,310,320,330,420,470,510,520,530,540,550,570,610,620,640,
                       660,670,680,690,710,720,730,830,850,999]

# ie/forkod.f JFOR(15)/KFOR(15): national-forest code list + geographic-location class.
const IE_JFOR = Int[103,104,105,106,621,110,113,114,116,117,118,613,102,109,112]
const IE_KFOR = Int[1,1,3,2,1,1,1,1,1,3,2,1,1,1,1]

# ie/sitset.f BAMAXA(30) and MAPSS(30) (site-species by ITYPE).
const IE_BAMAXA = Float32[140,220,250,310,240,270,310,310,200,310,290,330,380,440,500,500,390,390,440,180,
                          290,400,350,390,260,300,220,220,160,300]
const IE_MAPSS = Int[10,10, 3,3,3,3,3,3,3, 8,8, 4,4, 6,6,6, 5, 9,9,9,9,9, 11, 9, 11, 9,9,9, 12, 3]

# ie/sitset.f MAPSIT(30,MAXSP): default per-species site index by ITYPE (0 = keep input SITEAR).
# Non-zero only for LM(13), PI(15), JU(16), PY(17), AS(18), CO(19), MM(20), PB(21), OH(22).
const _IE_MAPSIT_LM = Int32[25,29,35,35,32,34,35,32,27,39,37,32,36,41,41,41,43,38,39,32,36,34,34,32,36,36,36,26,22,26]
const _IE_MAPSIT_PI = Int32[7,9,13,12,10,12,12,11,8,14,11,13,16,15,15,15,16,14,14,11,12,13,13,11,13,13,13,8,6,8]
const _IE_MAPSIT_JU = Int32[6,8,10,10,9,10,10,9,7,11,9,11,12,12,12,12,13,11,11,9,10,10,10,9,10,10,10,7,6,7]
const _IE_MAPSIT_PY = Int32[25,29,35,35,33,34,35,32,27,39,37,32,36,41,41,41,43,38,39,32,36,34,34,32,36,36,36,26,22,26]
const _IE_MAPSIT_AS = Int32[36,43,51,50,44,49,48,46,41,55,46,52,59,58,58,58,60,54,55,46,49,52,52,46,52,52,52,38,33,38]
const _IE_MAPSIT_CO = Int32[44,57,76,62,72,72,71,67,55,86,67,80,95,93,93,93,98,84,86,67,73,80,80,67,80,80,80,49,36,49]
function _ie_build_mapsit()
    m = zeros(Int32, 30, 23)
    m[:,13] = _IE_MAPSIT_LM; m[:,15] = _IE_MAPSIT_PI; m[:,16] = _IE_MAPSIT_JU; m[:,17] = _IE_MAPSIT_PY
    m[:,18] = _IE_MAPSIT_AS; m[:,19] = _IE_MAPSIT_CO; m[:,20] = _IE_MAPSIT_AS  # MM == AS
    m[:,21] = _IE_MAPSIT_AS; m[:,22] = _IE_MAPSIT_CO                            # PB == AS, OH == CO
    return m
end
const IE_MAPSIT = _ie_build_mapsit()

"ie/habtyp.f numeric path: habitat code KODTYP -> ITYPE (1..30). JTYPE(95) bracket search, ITYPE=KTYPE(K-1)."
function ie_habtyp(kodtyp::Integer)
    (kodtyp < 10 || kodtyp > 999) && return 0        # ERRGRO(14); ITYPE unchanged (caller keeps prior)
    kk = findfirst(k -> kodtyp < IE_JTYPE[k], 1:95)
    k = kk === nothing ? 96 : kk
    return Int(IE_KTYPE[k-1])
end

"ie/forkod.f: KODFOR -> (IFOR 1..11, IGL). Reservation special-cases + JFOR(15) match + IFOR remap 12/13/14/15."
function ie_forkod!(p)
    kodfor = Int(p.user_forest_code)
    ifor = 0; useigl = true
    # reservation pseudo-codes -> national forest IFOR (forkod.f SELECT CASE KODFOR)
    resmap = Dict(8106=>5, 8107=>10, 8109=>7, 8131=>5, 8132=>7, 8133=>6, 8137=>11)
    if haskey(resmap, kodfor)
        ifor = resmap[kodfor]
    else
        idx = findfirst(==(kodfor), IE_JFOR)
        if idx === nothing
            useigl = false; ifor = 1               # not found (ERRGRO 3); default to a safe subscript
        else
            ifor = idx
        end
    end
    # IFOR remap for the JFOR entries outside 1..11 (forkod.f second SELECT)
    ifor == 12 && (ifor = 7)                        # Kaniksu 613 -> 113
    ifor == 13 && (ifor = 1)                        # Beaverhead 102 -> Bitterroot 103
    ifor == 14 && (ifor = 1)                        # Deerlodge 109 -> Bitterroot 103
    ifor == 15 && (ifor = 9)                        # Helena 112 -> Lolo 116
    p.forest_idx = Int32(ifor)
    useigl && (p.geo_location = Int32(IE_KFOR[ifor]))
    return ifor
end

"ie/sitset.f: fill SITEAR (site index by species) from MAPSIT(ITYPE,sp) when 0, and the SDImax defaults."
function ie_sitset!(s::StandState, itype::Int)
    p = s.plot
    (itype < 1 || itype > 30) && (itype = 1)
    @inbounds for sp in 1:23
        if p.sp_site_index[sp] == 0f0
            p.sp_site_index[sp] = Float32(IE_MAPSIT[itype, sp])
        end
    end
    bamax = s.control.ba_max
    bamax <= 0f0 && (bamax = IE_BAMAXA[itype])
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    @inbounds for sp in 1:23
        p.sp_sdi_def[sp] <= 0f0 && bamax > 0f0 && (p.sp_sdi_def[sp] = bamax / (0.5454154f0 * pmsdiu))
    end
    return s
end

function ie_site_index_setup!(s::StandState)
    p = s.plot
    ie_forkod!(p)                                    # IFOR -> p.forest_idx, IGL -> p.geo_location
    kodtyp_in = Int(p.habitat_code)
    itype = kodtyp_in > 0 ? ie_habtyp(kodtyp_in) : Int(p.habitat_input)
    itype < 1 && (itype = Int(p.habitat_input))      # habtyp error path keeps prior ITYPE
    (itype < 1 || itype > 30) && (itype = 1)
    p.habitat_input = Int32(itype)                   # ITYPE for dgf!/ie_dgcons! (MAPHAB/MAPCCF)
    ie_sitset!(s, itype)
    return s
end

site_setup!(s::StandState, ::InlandEmpire) = ie_site_index_setup!(s)
