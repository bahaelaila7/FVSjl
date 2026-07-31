# =============================================================================
# habitat_tables.jl (kootenai) — KT habitat-type + SDImax lookup tables (kt/blkdat.f,
# kt/habtyp.f, kt/sitset.f). Faithful transcription (doctrine #2, measured from source);
# consumed by site_index.jl's habtyp 2-level mapping and sitset SDIDEF.
#
# The 2-level habtyp flow (kt/habtyp.f):
#   habitat code -> KODTYP; search KT_KOTHAB(175): first K with KODTYP < KOTHAB[K] -> KKTYPE = K-1
#     (KKTYPE is the "Kootenai habitat type" index 1..175; feeds DG MAPHAB(KKTYPE,ISPC)->DGHAB, chunk 3);
#   KODTYP2 = KT_KOTHAB[KKTYPE]; search KT_JTYPE(95): first K with KODTYP2 < JTYPE[K] -> ITYPE = KT_KTYPE[K-1]
#     (ITYPE is the "Inland Empire habitat type" index 1..30); rep code = KT_MTYPE[ITYPE].
# sitset (kt/sitset.f): BAMAX defaults to KT_BAMAXA[ITYPE]; SDIDEF = BAMAX/(0.5454154*(PMSDIU/100)).
# VALIDATED vs live FVSkt (stand 753200841290487): KOOTENAI=531, IE=530, SDI MAX=949 (ITYPE=14, BAMAXA[14]=440,
# PMSDIU=85 -> 949).
# =============================================================================

# KOTHAB(175): sorted Kootenai habitat codes (kt/blkdat.f:55)
const KT_KOTHAB = Int32[
    10,65,70,74,79,91,92,93,95,100,110,120,130,140,141,
    161,170,171,172,180,181,182,200,210,220,221,230,250,
    260,261,262,263,280,281,282,283,290,291,292,293,310,
    311,312,313,315,320,321,322,323,324,330,331,332,340,
    350,360,370,371,400,410,420,421,422,430,440,450,460,
    461,470,480,500,501,502,505,507,508,510,511,515,520,
    521,522,523,524,529,530,531,532,533,534,540,542,545,
    550,565,570,571,572,573,574,575,577,578,579,580,585,
    590,591,592,610,620,621,622,623,624,625,630,632,635,
    636,640,641,642,650,651,652,653,654,655,660,661,662,
    663,670,671,672,673,674,680,690,691,692,693,700,710,
    712,720,730,731,732,733,734,740,750,751,761,770,780,
    790,791,792,810,820,830,831,832,850,860,870,900,910,
    920,930,940,950]

# JTYPE(95): sorted Inland Empire habitat codes (kt/blkdat.f:160; 68 real values then 27*0 padding = 95)
const KT_JTYPE = Int32[
    10,100,110,130,140,160,170,180,190,200,
    210,220,230,250,260,280,290,310,320,330,
    340,350,360,370,380,400,410,420,430,440,
    450,460,470,480,500,501,502,505,506,510,
    515,516,520,529,530,540,545,550,555,560,
    565,570,575,579,590,600,610,620,630,635,
    640,650,660,670,675,680,685,690,700,701,
    710,720,730,740,750,770,780,790,800,810,
    820,830,840,850,860,870,890,900,910,920,
    925,930,940,950,999]   # note: source pads to 95 with trailing 0s; the search never reaches them

# KTYPE(95): habitat-code-slot -> ITYPE (1..30) (kt/habtyp.f DATA KTYPE)
const KT_KTYPE = Int32[
    1,1,1,1,1, 2,2,2,2, 4, 1,1,1, 3,4,5,6,7,8,9,8,8,9,7,3, 10,10,10, 4,11,20,29,11,
    11,13,14,17, 12,12,12,12, 13,13,13, 14,15,14,16,14,16, 17,17,17, 24,12,24,18,
    19,21,19,20,20,21,22,19,23,19,24,27,25,25,26,27,22,24,27,24,
    27,28,28,29,28,28, 29,29,29,29, 27,9,20,24,21,27,24,30]

# MTYPE(30): ITYPE -> IE representative habitat code (kt/habtyp.f DATA MTYPE)
const KT_MTYPE = Int32[
    130,170,250,260,280,290,310,320,330,420,
    470,510,520,530,540,550,570,610,620,640,
    660,670,680,690,710,720,730,830,850,999]

# BAMAXA(30): ITYPE -> max BA (kt/sitset.f:32); SDIDEF = BAMAX/(0.5454154*(PMSDIU/100))
const KT_BAMAXA = Float32[
    140,220,250,310,240,270,310,310,200,310,
    290,330,380,440,500,500,390,390,440,180,
    290,400,350,390,260,300,220,220,160,300]

# KT_BKRAT(11): bark ratio per species (kt/blkdat.f BKRAT) — for the later bark chunk (bratio).
const KT_BKRAT = Float32[0.964,0.851,0.867,0.915,0.934,0.950,0.969,0.956,0.937,0.890,0.934]
