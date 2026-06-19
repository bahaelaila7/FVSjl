# =============================================================================
# species.jl — Southern variant species tables & BLOCK DATA defaults
#
# Ported from: sn/blkdat.f  (BLOCK DATA BLKDAT)
#
# The 90 Southern species and their identity codes (alpha / FIA SPCD / PLANTS),
# plus the diameter-growth regression standard errors and the BLOCK DATA scalar
# defaults (RNG seed, default tree-record FORMAT, etc.). These immutable tables
# live as `const` arrays in the variant; per-stand state is filled by
# `load_species_coefficients!` / `init_blockdata!` below.
# =============================================================================

"Number of Southern variant species."
const SN_NSPECIES = 90

"4-char alpha species codes (JSP)."
const SN_ALPHA = String[
    "FR  ","JU  ","PI  ","PU  ","SP  ","SA  ","SR  ","LL  ","TM  ","PP  ",
    "PD  ","WP  ","LP  ","VP  ","BY  ","PC  ","HM  ","FM  ","BE  ","RM  ",
    "SV  ","SM  ","BU  ","BB  ","SB  ","AH  ","HI  ","CA  ","HB  ","RD  ",
    "DW  ","PS  ","AB  ","AS  ","WA  ","BA  ","GA  ","HL  ","LB  ","HA  ",
    "HY  ","BN  ","WN  ","SU  ","YP  ","MG  ","CT  ","MS  ","MV  ","ML  ",
    "AP  ","MB  ","WT  ","BG  ","TS  ","HH  ","SD  ","RA  ","SY  ","CW  ",
    "BT  ","BC  ","WO  ","SO  ","SK  ","CB  ","TO  ","LK  ","OV  ","BJ  ",
    "SN  ","CK  ","WK  ","CO  ","RO  ","QS  ","PO  ","BO  ","LO  ","BK  ",
    "WI  ","SS  ","BD  ","EL  ","WE  ","AE  ","RL  ","OS  ","OH  ","OT  ",
]

"FIA SPCD codes as text (FIAJSP)."
const SN_FIA = String[
    "010","057","090","107","110","111","115","121","123","126",
    "128","129","131","132","221","222","260","311","313","316",
    "317","318","330","370","372","391","400","450","460","471",
    "491","521","531","540","541","543","544","552","555","580",
    "591","601","602","611","621","650","651","652","653","654",
    "660","680","691","693","694","701","711","721","731","740",
    "743","762","802","806","812","813","819","820","822","824",
    "825","826","827","832","833","834","835","837","838","901",
    "920","931","950","970","971","972","975","299","998","999",
]

"PLANTS symbols (PLNJSP)."
const SN_PLANTS = String[
    "ABIES ","JUNIP ","PICEA ","PICL  ","PIEC2 ","PIEL  ","PIGL2 ",
    "PIPA2 ","PIPU5 ","PIRI  ","PISE  ","PIST  ","PITA  ","PIVI2 ",
    "TADI2 ","TAAS  ","TSUGA ","ACBA3 ","ACNE2 ","ACRU  ","ACSA2 ",
    "ACSA3 ","AESCU ","BETUL ","BELE  ","CACA18","CARYA ","CATAL ",
    "CELTI ","CECA4 ","COFL2 ","DIVI5 ","FAGR  ","FRAXI ","FRAM2 ",
    "FRNI  ","FRPE  ","GLTR  ","GOLA  ","HALES ","ILOP  ","JUCI  ",
    "JUNI  ","LIST2 ","LITU  ","MAGNO ","MAAC  ","MAGR4 ","MAVI2 ",
    "MAMA2 ","MALUS ","MORUS ","NYAQ2 ","NYSY  ","NYBI  ","OSVI  ",
    "OXAR  ","PEBO  ","PLOC  ","POPUL ","POGR4 ","PRSE2 ","QUAL  ",
    "QUCO2 ","QUFA  ","QUPA5 ","QULA2 ","QULA3 ","QULY  ","QUMA3 ",
    "QUMI  ","QUMU  ","QUNI  ","QUPR2 ","QURU  ","QUSH  ","QUST  ",
    "QUVE  ","QUVI  ","ROPS  ","SALIX ","SAAL5 ","TILIA ","ULMUS ",
    "ULAL  ","ULAM  ","ULRU  ","2TN   ","2TB   ","2TREE ",
]

"Diameter-growth regression standard errors per species (SIGMAR)."
const SN_SIGMAR = Float32[
    0.451100, 0.529700, 0.451100, 0.542800, 0.498700, 0.525100, 0.436700,
    0.441000, 0.469300, 0.552500, 0.592100, 0.493700, 0.468700, 0.469300,
    0.551100, 0.626700, 0.451100, 0.563200, 0.560800, 0.593000, 0.593000,
    0.475500, 0.537300, 0.569600, 0.569600, 0.603200, 0.499300, 0.440100,
    0.527600, 0.545300, 0.538200, 0.516000, 0.480500, 0.595800, 0.422800,
    0.485600, 0.485600, 0.468200, 0.590800, 0.599600, 0.546900, 0.571000,
    0.571000, 0.577900, 0.518100, 0.572700, 0.512600, 0.570300, 0.572700,
    0.570300, 0.505600, 0.505600, 0.569800, 0.535600, 0.588800, 0.577300,
    0.504700, 0.568700, 0.557000, 0.440100, 0.440100, 0.578100, 0.440700,
    0.382700, 0.430500, 0.400000, 0.495700, 0.558600, 0.465900, 0.464900,
    0.491300, 0.505600, 0.466400, 0.429300, 0.404800, 0.407400, 0.485300,
    0.421900, 0.663500, 0.518700, 0.450800, 0.450400, 0.549600, 0.644700,
    0.528600, 0.536400, 0.535000, 0.529700, 0.577300, 0.557600,
]

"Valid habitat type codes (JTYPE), 122 slots (95 codes + 27 zero pads)."
const SN_VALID_HABITAT = Int32[
     10,100,110,130,140,160,170,180,190,200,
    210,220,230,250,260,280,290,310,320,330,
    340,350,360,370,380,400,410,420,430,440,
    450,460,470,480,500,501,502,505,506,510,
    515,516,520,529,530,540,545,550,555,560,
    565,570,575,579,590,600,610,620,630,635,
    640,650,660,670,675,680,685,690,700,701,
    710,720,730,740,750,770,780,790,800,810,
    820,830,840,850,860,870,890,900,910,920,
    925,930,940,950,999,
    (zeros(Int32, 27))...,
]

# BLOCK DATA scalar defaults (sn/blkdat.f). The Fortran unit numbers are kept as
# named fields rather than a global unit table (see io layer). The RNG main stream
# is seeded to 55329 here (blkdat sets S0=SS=55329), matching the establishment seed.
const SN_RNG_SEED = 55329.0f0
const SN_REGEN_BARK = 2.999f0      # REGNBK

"""
    init_blockdata!(s::StandState, ::Southern)

Apply the Southern BLOCK DATA defaults to a fresh state: species identity tables,
default tree FORMAT, RNG seed, and the handful of scalar defaults set in blkdat.f.
Called once at the start of stand initialization.
"""
function init_blockdata!(s::StandState, ::Southern)
    sd = s.species
    @inbounds for i in 1:SN_NSPECIES
        sd.alpha[i]  = SN_ALPHA[i]
        sd.fia[i]    = SN_FIA[i]
        sd.plants[i] = SN_PLANTS[i]
        code = rstrip(SN_ALPHA[i])
        sd.class_codes[i, 1] = code * "1"
        sd.class_codes[i, 2] = code * "2"
        sd.class_codes[i, 3] = code * "3"
    end
    copyto!(s.plot.valid_habitat, 1, SN_VALID_HABITAT, 1, length(SN_VALID_HABITAT))

    s.control.tree_format = DEFAULT_TREE_FORMAT
    s.control.year = 5.0f0                                # YR default cycle length

    # RNG: both streams seeded 55329 (blkdat S0/SS + ESBLKD ESS0/ESSS)
    s.rng.s0 = Float64(SN_RNG_SEED); s.rng.ss = SN_RNG_SEED
    return s
end

"""
    load_species_coefficients!(s, ::Southern)

Variant hook: load species data into the stand state. For now this is the
identity/BLOCK-DATA load; growth coefficients (DGF/HTGF tables) are attached in C3.
"""
load_species_coefficients!(s::StandState, v::Southern) = init_blockdata!(s, v)
