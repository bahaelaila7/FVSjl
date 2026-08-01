# =============================================================================
# site_index.jl (easternmontana) — EM habtyp (em/habtyp.f) + sitset (em/sitset.f), chunk 2.
#
# Produces the inputs later chunks read:
#   ITYPE  (1..30) -> p.habitat_input   (em/habtyp.f: JTYPE(118) bucket -> IEMTYP -> NIHMAP)
#   IEMTYP (1..118)                      (used for the MAPSDI SDI-group lookup)
#   SITEAR (per-sp) -> p.sp_site_index   (em/sitset.f: MAPSIT(ITYPE, site-group(sp)))
#   ISISP  -> p.site_species             (MAPSS(ITYPE))
# SDImax defaults (SDIDEF): EM-SPECIFIC — SDICON(MAPSDI(IEMTYP)) when BAMAX<=0 (NOT BAMAXA-derived
# like IE); BAMAX defaults to BAMAXA(ITYPE) AFTER the SDIDEF loop. All tables extracted from
# em/blkdat.f + em/habtyp.f + em/sitset.f (tools/easternmontana/extract_{habitat,sitset}.py).
#
# VALIDATED vs live FVSem emt01 (STDINFO hab 260 -> IEMTYP 29, ITYPE 4): SDIDEF=SDICON(3)=696 for
# all 19 species (= live "SDI MAX" table); BAMAX=BAMAXA(4)=310; ISISP=MAPSS(4)=3 (DF); DF SI=MAPSIT(4,8)=51.
# =============================================================================

# em/blkdat.f JTYPE(118): sorted valid habitat codes (bucket search key).
const EM_JTYPE = Int32[
    10, 65, 70, 74, 79, 91, 92, 93, 95, 100, 110, 120, 130, 140, 141, 161, 170, 171, 172, 180,
    181, 182, 200, 210, 220, 221, 230, 250, 260, 261, 262, 280, 281, 282, 283, 290, 291, 292, 293, 310,
    311, 312, 313, 315, 320, 321, 322, 323, 330, 331, 332, 340, 350, 360, 370, 371, 400, 410, 430, 440,
    450, 460, 461, 470, 480, 591, 610, 620, 624, 625, 630, 632, 640, 641, 642, 650, 651, 653, 654, 655,
    660, 661, 662, 663, 670, 674, 690, 691, 692, 700, 710, 720, 730, 731, 732, 733, 740, 750, 751, 770,
    780, 790, 791, 792, 810, 820, 830, 832, 850, 860, 870, 900, 910, 920, 930, 940, 950, 999]
# em/habtyp.f NIHMAP(118): bucket index IEMTYP -> ITYPE (1..30 NI habitat types).
const EM_NIHMAP = Int32[
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
    2, 2, 4, 1, 1, 1, 1, 3, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7,
    7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 8, 8, 9, 7, 7, 10, 10, 4, 11,
    20, 29, 29, 11, 11, 12, 18, 19, 19, 19, 21, 21, 20, 20, 20, 20, 20, 20, 20, 20,
    21, 21, 21, 21, 22, 22, 24, 24, 24, 27, 25, 26, 27, 27, 27, 27, 22, 24, 24, 27,
    24, 27, 27, 27, 28, 29, 28, 28, 29, 29, 29, 27, 9, 20, 21, 27, 24, 30]
# em/sitset.f BAMAXA(30) BA-max/ITYPE, SDICON(9) SDImax/group, MAPSS(30) site-species/ITYPE.
const EM_BAMAXA = Float32[140, 220, 250, 310, 240, 270, 310, 310, 200, 310, 290, 330, 380, 440, 500, 500, 390, 390, 440, 180, 290, 400, 350, 390, 260, 300, 220, 220, 160, 300]
const EM_SDICON = Float32[467, 634, 696, 768, 775, 751, 707, 661, 635]
const EM_MAPSS  = Int[10, 10, 3, 3, 3, 3, 3, 3, 3, 8, 8, 3, 3, 3, 3, 3, 3, 9, 9, 9, 9, 9, 3, 9, 3, 9, 9, 9, 3, 3]
# em/sitset.f MAPSDI(122): IEMTYP -> SDICON index (1..9).
const EM_MAPSDI = Int[
    2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 1, 1, 1, 1, 1, 4, 3, 3, 3, 5, 5, 5, 5, 4, 4, 4, 4, 3,
    1, 3, 3, 3, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 7, 7, 3, 7,
    7, 9, 9, 4, 4, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    4, 4, 4, 4, 8, 8, 5, 5, 5, 7, 7, 5, 6, 6, 6, 7, 8, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 7, 7, 7, 4, 6, 9, 5, 5, 5,
    5, 5]
# em/sitset.f MAPSIT(30,11): site index by ITYPE (row) x site-species-group (col). Col-major DATA
# transcribed row-major here; EM_MAPSIT[itype, group].
const EM_MAPSIT = Int32[
    25 6 36 44 36 35 46 26 46 46 46;
    29 8 43 57 43 51 68 44 68 68 68;
    35 10 51 76 51 57 76 49 76 76 76;
    35 10 50 62 50 60 80 51 80 78 80;
    32 9 44 72 44 48 64 41 75 64 64;
    34 10 49 72 49 55 73 47 73 73 73;
    35 10 48 71 48 45 67 39 72 61 61;
    32 9 46 67 46 38 67 33 70 51 51;
    27 7 41 55 41 44 59 38 63 59 59;
    39 11 55 86 55 60 80 51 80 80 80;
    37 9 46 67 46 64 85 42 80 85 85;
    32 11 52 80 52 47 62 40 62 68 62;
    36 12 59 95 59 47 62 40 62 68 62;
    41 12 58 93 58 47 62 40 62 68 62;
    41 12 58 93 58 47 62 40 62 68 62;
    41 12 58 93 58 47 62 40 62 68 62;
    43 13 60 98 60 47 62 40 62 68 62;
    38 11 54 84 54 47 62 40 62 68 62;
    39 11 55 86 55 47 62 40 62 68 62;
    32 9 46 67 46 53 70 46 74 75 70;
    36 10 49 73 49 50 67 43 76 68 67;
    34 10 52 80 52 57 76 49 76 78 76;
    34 10 52 80 52 49 65 42 68 80 76;
    32 9 46 67 46 57 76 42 68 80 76;
    36 10 52 80 52 49 65 42 68 80 76;
    36 10 52 80 52 49 65 42 65 67 65;
    36 10 52 80 52 47 62 30 35 64 62;
    26 7 38 49 38 41 55 36 51 53 55;
    22 6 33 36 33 26 35 23 35 45 35;
    26 7 38 49 38 60 80 51 80 78 80;
]

# em/forkod.f JFOR(6)/KFOR(6): national-forest location codes + geographic-location class (all 1).
const EM_JFOR = Int[102, 108, 109, 111, 112, 115]
const EM_KFOR = Int[1, 1, 1, 1, 1, 1]

# em/forkod.f: translate the user forest location code KODFOR → IFOR (1..6, the JFOR subscript that
# MAPLOC/MAPDSQ index in the DG DGCONS) + IGL=KFOR[IFOR] (geographic location, regen/htdbh). Reservation
# pseudo-codes 7xxx map to IFOR directly; else match KODFOR against JFOR (not found → ERRGRO(3), IFOR=1).
function em_forkod!(p)
    kodfor = Int(p.user_forest_code)
    ifor = 1; useigl = true
    if kodfor in (7101, 7102, 7103, 7107, 7108, 7109, 7302, 7305)
        ifor = 2
    elseif kodfor in (7301, 7303, 7304, 7307)
        ifor = 6
    else
        idx = findfirst(==(kodfor), EM_JFOR)
        idx === nothing ? (useigl = false; ifor = 1) : (ifor = idx)
    end
    p.forest_idx = Int32(ifor)
    useigl && (p.geo_location = Int32(EM_KFOR[ifor]))
    return ifor
end

# em/habtyp.f: bucket KODTYP into JTYPE -> IEMTYP (largest idx with JTYPE<=KODTYP), ITYPE=NIHMAP(IEMTYP).
# Returns (iemtyp, itype). KODTYP<JTYPE[1] or <=0 -> defaults (iemtyp 1).
function em_habtyp(kodtyp::Integer)
    kodtyp <= 0 && return (1, Int(EM_NIHMAP[1]))
    jj = findfirst(k -> kodtyp < EM_JTYPE[k], 1:118)
    iemtyp = jj === nothing ? 118 : jj - 1
    iemtyp < 1 && (iemtyp = 1)
    return (iemtyp, Int(EM_NIHMAP[iemtyp]))
end

# em/sitset.f species -> MAPSIT site-species-group (1..11); else fixed SI 70.
@inline function _em_site_group(i::Int)
    (i == 4 || i == 5)               && return 1   # LM, LL
    i == 6                           && return 2   # RM
    i == 12                          && return 3   # AS
    (i == 11 || (13 <= i <= 16) || i == 19) && return 4  # GA,CW,BA,PW,NC,OH
    i == 17                          && return 5   # PB
    (i == 1 || i == 2 || i == 18)    && return 6   # WB, WL, OS
    i == 10                          && return 7   # PP
    i == 3                           && return 8   # DF
    i == 7                           && return 9   # LP
    i == 8                           && return 10  # ES
    i == 9                           && return 11  # AF
    return 0                                       # else -> fixed 70.
end

# em/sitset.f: fill SITEAR (site index by species) + ISISP + SDImax defaults.
function em_sitset!(s::StandState, itype::Int, iemtyp::Int)
    p = s.plot
    (itype < 1 || itype > 30) && (itype = 1)
    p.site_species == 0 && (p.site_species = Int32(EM_MAPSS[itype]))   # ISISP = MAPSS(ITYPE)
    @inbounds for i in 1:nspecies(s.variant)
        if p.sp_site_index[i] == 0f0
            g = _em_site_group(i)
            p.sp_site_index[i] = g == 0 ? 70f0 : Float32(EM_MAPSIT[itype, g])
        end
    end
    # SDIDEF: EM uses SDICON(MAPSDI(IEMTYP)) when keyword BAMAX<=0 (NOT BAMAXA-derived). BAMAX
    # is applied AFTER this loop, so during it BAMAX is only the keyword value.
    isdi = (iemtyp > 0) ? EM_MAPSDI[iemtyp] : 5
    (isdi < 1 || isdi > 9) && (isdi = 5)
    bamax_kw = s.control.ba_max
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    @inbounds for i in 1:nspecies(s.variant)
        if p.sp_sdi_def[i] <= 0f0
            p.sp_sdi_def[i] = bamax_kw > 0f0 ? bamax_kw / (0.5454154f0 * pmsdiu) : Float32(EM_SDICON[isdi])
        end
    end
    s.control.ba_max <= 0f0 && (s.control.ba_max = EM_BAMAXA[itype])   # BAMAX = BAMAXA(ITYPE)
    return s
end

function em_site_index_setup!(s::StandState)
    p = s.plot
    em_forkod!(p)                                    # IFOR → p.forest_idx (DG MAPLOC/MAPDSQ), IGL → p.geo_location
    kodtyp_in = Int(p.habitat_code)
    iemtyp, itype = kodtyp_in > 0 ? em_habtyp(kodtyp_in) : (Int(p.habitat_input) > 0 ? (0, Int(p.habitat_input)) : (1, Int(EM_NIHMAP[1])))
    (itype < 1 || itype > 30) && (itype = 1)
    p.habitat_input = Int32(itype)
    em_sitset!(s, itype, iemtyp)
    return s
end

site_setup!(s::StandState, ::EasternMontana) = em_site_index_setup!(s)
