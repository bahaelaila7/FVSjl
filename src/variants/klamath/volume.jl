# =============================================================================
# volume.jl (klamath) — NC volume (chunk 8). NVEL Region-5/6, dispatched by VEQNNC geocode.
#
# ⚠ NC's VEQNNC is FOREST/REGION-DEPENDENT (sitset.f → VOLEQDEF(VAR='NC',IREGN=KODFOR/100,FORST)):
#   • REGION-5 forests (Klamath/Six Rivers/Trinity, IFOR 1-3) — the NC_VOL_EQ table below:
#       8 conifers + redwood = "500WO2W<fia>"  → Region-5 Wensel-Krumland R5TAP taper (r5tap.f)  [nc_wo2w_vol]
#       4 hardwoods MA/BO/TO/OH = "500DVEW<fia>" → Region-5 DVE California hardwood D²H (r5harv.f) [nc_r5harv_vol]
#   • REGION-6 forest SISKIYOU (IFOR 4, forest 611) — the NC_R6_VOL_EQ table further down: westside Flewelling
#       F06FW2W202 (DF) + INGY I00FW2W093/073 (WF/PP) + Region-6 Behre 616BEHW<fia> (all else). See that block.
# nct01 = 100% conifers (SP/DF/WF/RF) on a Region-5 forest ⇒ exercises only the R5 WO2W path; the DVEW hardwood
# path is ported faithfully from r5harv.f but not reachable on nct01 (no hardwood stand).
#
# Merch: NC grinit.f leaves TOPD/DBHMIN/STMP at 0/0/1 ⇒ the NVEL equation-default merch specs apply
# (MTOPP resolved by the volume driver — see nc_wo2w_vol / the merch-top plumbing).
# =============================================================================

# NC volume equations — species index 1..12 (OS SP DF WF MA IC BO TO RF PP OH RW).
# CONFIRMED bit-exact vs live FVSnc (nct01.out:68-70). Source: voleqdef.f R5_EQN FIA binary-search table.
const NC_VOL_EQ = String[
    "500WO2W108",  # 1  OS  other softwoods → lodgepole taper
    "500WO2W117",  # 2  SP  sugar pine
    "500WO2W202",  # 3  DF  Douglas-fir
    "500WO2W015",  # 4  WF  white fir
    "500DVEW361",  # 5  MA  Pacific madrone   (DVE hardwood, SPEC 10)
    "500WO2W081",  # 6  IC  incense cedar
    "500DVEW818",  # 7  BO  California black oak (DVE hardwood, SPEC 3)
    "500DVEW631",  # 8  TO  tanoak            (DVE hardwood, SPEC 12)
    "500WO2W020",  # 9  RF  California red fir
    "500WO2W122",  # 10 PP  ponderosa pine
    "500DVEW981",  # 11 OH  other hardwoods → laurel (DVE hardwood, SPEC 9)
    "500WO2W211",  # 12 RW  redwood (West-side taper, NOT DVE)
]

# ---------------------------------------------------------------------------
# DVEW: R5HARV — Region-5 California hardwood direct-volume estimator (r5harv.f).
# From PNW-414 (Pillsbury & Kirkley): total/wood/saw-log cubic for 15 CA hardwoods.
# VOLEQ(8:10) → SPEC index (r5harv.f:116-152). NC uses only the misc-hardwood power-law branch
# (SPEC 2-13,15); juniper(0)/red-alder(1)/giant-sequoia(14) branches ported for completeness/reuse.
# ---------------------------------------------------------------------------

# COFA/COFB/COFC(15,4): CV4 / CV8 / CVT power-law coefficients (a·D^b·H^c·IV^d, IV=10). r5harv.f:36-97.
const NC_R5_COFA = (
    (0.0,0.0,0.0,0.0),                                # 1  red alder (SPEC 1 uses own branch)
    (0.0034214162,2.35347,0.69586,0.0),               # 2  bigleaf maple 76,95
    (0.0036795695,2.12635,0.83339,0.0),               # 3  California black oak 81  ← NC BO
    (0.0042324071,2.53987,0.50591,0.0),               # 4  blue oak 88
    (0.0031670596,2.32519,0.74348,0.0),               # 5  canyon live oak 84
    (0.0055212937,2.07202,0.77467,0.0),               # 6  giant chinkapin 93
    (0.0024574847,2.53284,0.60764,0.0),               # 7  coast live oak 82,96,98
    (0.0041192264,2.14915,0.77843,0.0),               # 8  interior live oak 85
    (0.0016380753,2.05910,1.05293,0.0),               # 9  California laurel 91  ← NC OH(981)
    (0.0025616425,1.99295,1.01532,0.0),               # 10 Pacific madrone 94    ← NC MA(361)
    (0.0024277027,2.25575,0.87108,0.0),               # 11 Oregon white oak 86
    (0.000577497, 2.19576,1.14078,0.0),               # 12 tanoak 87,72,73,75    ← NC TO(631)
    (0.0009684363,2.39565,0.98878,0.0),               # 13 California white oak 83
    (0.0,0.0,0.0,0.0),                                # 14 (giant sequoia — own branch)
    (0.0053866353,2.61268,0.31103,0.0),               # 15 Engelmann oak 79
)
const NC_R5_COFB = (
    (0.0,0.0,0.0,0.0),
    (0.0004236332,2.10316,1.08584,0.40017),           # 2
    (0.0012478663,2.68099,0.42441,0.28385),           # 3  BO
    (0.0036912408,1.79732,0.838884,0.15958),          # 4
    (0.0006540144,2.24437,0.81358,0.43381),           # 5
    (0.0018985111,2.38285,0.77105,0.0),               # 6
    (0.0006540144,2.24437,0.81358,0.43381),           # 7
    (0.0006540144,2.24437,0.81358,0.43381),           # 8
    (0.0007741517,2.23009,1.037,  0.0),               # 9  OH
    (0.000618153, 1.72635,1.26462,0.37867),           # 10 MA
    (0.0008281647,2.10651,0.91215,0.32652),           # 11
    (0.0002526443,2.30949,1.21069,0.0),               # 12 TO
    (0.0001880044,1.87346,1.62443,0.0),               # 13
    (0.0,0.0,0.0,0.0),                                # 14
    (0.0001880044,1.87346,1.62443,0.0),               # 15 (white oak COFB per 2018/10/30 fix)
)
const NC_R5_COFC = (
    (0.0,0.0,0.0,0.0),
    (0.0101786350,2.22462,0.57561,0.0),               # 2
    (0.0070538108,1.97437,0.85034,0.0),               # 3  BO
    (0.0125103008,2.33089,0.46100,0.0),               # 4
    (0.0097438611,2.20527,0.61190,0.0),               # 5
    (0.0120372263,2.02232,0.68638,0.0),               # 6
    (0.0065261029,2.31958,0.62528,0.0),               # 7
    (0.0136818837,2.02989,0.63257,0.0),               # 8
    (0.0057821322,1.94553,0.88389,0.0),               # 9  OH
    (0.0067322665,1.96628,0.83458,0.0),               # 10 MA
    (0.0072695058,2.14321,0.74220,0.0),               # 11
    (0.0058870024,1.94165,0.86562,0.0),               # 12 TO
    (0.0042870077,2.33631,0.74872,0.0),               # 13
    (0.0,0.0,0.0,0.0),                                # 14
    (0.0191453191,2.40248,0.28060,0.0),               # 15
)
const NC_R5_SEQB = (0.001682608,1.755956,1.490641)    # giant sequoia VOL(2)
const NC_R5_SEQC = (0.002438339,1.694874,1.098957)    # giant sequoia VOL(4)

"VOLEQ(8:10) FIA code → R5HARV SPEC index (r5harv.f:116-152). -1 = unsupported."
function _nc_r5_spec(fia::AbstractString)::Int
    (fia == "060" || fia == "064") && return 0        # juniper
    fia == "351" && return 1                          # red alder
    fia == "312" && return 2                          # bigleaf maple
    fia == "818" && return 3                          # CA black oak
    fia == "807" && return 4                          # blue oak
    fia == "805" && return 5                          # canyon live oak
    fia == "431" && return 6                          # giant chinkapin
    fia == "801" && return 7                          # coast live oak
    fia == "839" && return 8                          # interior live oak
    fia == "981" && return 9                          # laurel (other hardwoods)
    fia == "361" && return 10                         # Pacific madrone
    fia == "815" && return 11                         # Oregon white oak
    fia == "631" && return 12                         # tanoak
    fia == "821" && return 13                         # CA white oak
    fia == "212" && return 14                         # giant sequoia
    fia == "811" && return 15                         # Engelmann oak
    return -1
end

"R5HARV DVE cubic+board (r5harv.f) as DVEST returns it (dvest.f:170 `VOL(2)=ANINT(VOL(2))` — board feet to the
nearest whole foot; VOL(1)/VOL(4) are NOT rounded since the 2025/05/07 NVEL change). Returns (tcuft VOL(1),
MERCH-cuft = VOL(4)+VOL(7), scribner-bf VOL(2)) — the FVS .sum aggregates MCF = VOL(4)+VOL(7) (fvsvol.f:512) and
BDFT = VOL(2). `mtopp` = merch top DIB (TOPC): [3,5)→CV4, [5,7)→CV6, [7,9]→CV8, <3→CVT (whole stem).
Precision is Fortran's: COFA/COFB/COFC/COEFSEQ* are REAL*8 arrays DATA-initialised from default-REAL literals
(⇒ Float64(Float32(lit))), each power product is evaluated in double and stored to a REAL CVx; everything after
(CV6, BA, TARIF, the Scribner ratios, juniper, red alder) is single precision with glibc powf/expf/log10f.
⚠ These are the FULL (unbroken-height) volumes — the driver applies CFTOPK/BFTOPK (r4_topkill) for broken tops
exactly as fvsvol.f:193/391 does for METHC=6 (both WO2W conifers AND DVE hardwoods go through NATCRS→CFTOPK)."
function nc_r5harv_vol(voleq::AbstractString, d::Float32, h::Float32, mtopp::Float32)
    length(voleq) < 10 && return (0f0, 0f0, 0f0)
    d < 1f0 && return (0f0, 0f0, 0f0)                  # r5harv.f:84 ERRFLAG=3
    spec = _nc_r5_spec(voleq[8:10])
    spec < 0 && return (0f0, 0f0, 0f0)
    D = d; H = h; TOPC = mtopp
    vol1 = 0f0; vol2 = 0f0; vol4 = 0f0; vol7 = 0f0
    if spec == 0                                       # juniper (060/064), r5harv.f:158-178
        if D < 5f0 || H < 10f0
            cvts = 0.00272708f0 * D * D * H; v = 0f0
        else
            f = 0.307f0 + 0.00086f0*H - 0.0037f0*D*H / (H - 4.5f0)
            ba = 0.005454154f0 * D * D
            cvts = ba * f * H * fpowi(H / (H - 4.5f0), 2)
            v = TOPC > 0f0 ? (cvts + 3.48f0) / (1.18052f0 + 0.32736f0 * fexp(-0.1f0 * D)) - 2.948f0 : cvts
        end
        vol4 = round(v*10f0 + 0.5f0, RoundNearestTiesAway) / 10.0f0     # VOL(4)=ANINT(V*10+0.5)/10.0
        vol1 = cvts
    elseif spec == 1                                   # red alder (351), r5harv.f:180-285
        ba = D * D * 0.005454154f0
        term1 = (1.033f0 * (1.0f0 + 1.382937f0 * fexp(-4.015292f0 * (D / 10.0f0)))) * (ba + 0.087266f0) - 0.174533f0
        dlog = _r5_log10(D); hlog = _r5_log10(H)
        cvts = fpow(10.0f0, -2.672775f0 + 1.920617f0 * dlog + 1.074024f0 * hlog)
        tarif = (cvts * 0.912733f0) / term1
        tarif <= 0f0 && (tarif = 0.01f0)
        cv4 = tarif * (ba - 0.087266f0) / 0.912733f0
        cv8 = cv4 * (0.983f0 - 0.983f0 * fpow(0.65f0, D - 8.6f0))
        cvt = tarif * (0.9679f0 - 0.1051f0 * fpow(0.5523f0, D - 1.5f0)) * term1 / 0.912733f0
        cv6 = cv4 * (0.993f0 - 0.993f0 * fpow(0.62f0, D - 6.0f0))
        vol4 = TOPC >= 3f0 && TOPC < 5f0 ? cv4 : TOPC >= 5f0 && TOPC < 7f0 ? cv6 :
               TOPC >= 7f0 && TOPC <= 9f0 ? cv8 : TOPC < 3f0 ? cvt : 0f0
        vol1 = cvt
        vol4 <= 0f0 && (vol4 = 0f0)
        if D >= 7.0f0
            b4 = tarif / 0.912733f0
            balog = _r5_log10(b4)
            rs616l = 0.174439f0 + 0.117594f0 * dlog * balog - 8.210585f0 / fpowi(D, 2) + 0.236693f0 * balog -
                     0.00001345f0 * fpowi(b4, 2) - 0.00001937f0 * fpowi(D, 2)
            sv616 = fpow(10.0f0, rs616l) * cv6
            sv816 = (0.99f0 - 0.58f0 * fpow(0.484f0, D - 9.5f0)) * sv616
            vol2 = TOPC >= 5f0 && TOPC < 7f0 ? sv616 : TOPC >= 7f0 && TOPC <= 9f0 ? sv816 : 0f0
            vol7 = TOPC <= 0f0 ? 0f0 : TOPC == 6f0 ? cv4 - cv6 : TOPC == 8f0 ? cv4 - cv8 : 0f0
        end
    elseif spec == 14                                  # giant sequoia, r5harv.f:293-296
        b = NC_R5_SEQB; c = NC_R5_SEQC
        vol2 = Float32(_r8(b[1]) * dpow(Float64(D), _r8(b[2])) * dpow(Float64(H), _r8(b[3])))
        vol4 = Float32(_r8(c[1]) * dpow(Float64(D), _r8(c[2])) * dpow(Float64(H), _r8(c[3])))
        vol1 = vol4
    else                                               # misc hardwoods, r5harv.f:300-408
        a = NC_R5_COFA[spec]; b = NC_R5_COFB[spec]; cc = NC_R5_COFC[spec]
        cv4 = _r5_pw(a, D, H); cv8 = _r5_pw(b, D, H)
        cv6 = (cv4 > 0f0 && cv8 > 0f0) ? cv4 - ((cv4 - cv8) * 0.4f0) : 0f0
        cvt = _r5_pw(cc, D, H)
        cuftgros = TOPC >= 3f0 && TOPC < 5f0 ? cv4 : TOPC >= 5f0 && TOPC < 7f0 ? cv6 :
                   TOPC >= 7f0 && TOPC <= 9f0 ? cv8 : TOPC < 3f0 ? cvt : 0f0
        vol4 = cuftgros; vol7 = cv4 - cuftgros; vol1 = cvt
        if D >= 5.0f0                                  # SCRIBNER VOL(2), TOPB = MTOPP
            if D < 11.0f0
                vol2 = TOPC >= 5f0 && TOPC < 7f0 ? cv6 * 4f0 : TOPC >= 7f0 && TOPC <= 9f0 ? cv8 * 4f0 : 0f0
            else
                ba = fpowi(D, 2) * 0.005454154f0
                tarif = (cv8 * 0.912733f0) / ((0.983f0 - 0.983f0 * fpow(0.65f0, D - 8.6f0)) * (ba - 0.087266f0))
                tarif <= 0f0 && (tarif = 0.01f0)
                b4 = tarif / 0.912733f0
                dlog = _r5_log10(D); balog = _r5_log10(b4)
                rs616l = 0.174439f0 + 0.117594f0 * dlog * balog - 8.210585f0 / fpowi(D, 2) + 0.236693f0 * balog -
                         0.00001345f0 * fpowi(b4, 2) - 0.00001937f0 * fpowi(D, 2)
                sv616 = fpow(10.0f0, rs616l) * cv6
                sv816 = (0.99f0 - 0.58f0 * fpow(0.484f0, D - 9.5f0)) * sv616
                vol2 = TOPC >= 5f0 && TOPC < 7f0 ? sv616 : TOPC >= 7f0 && TOPC <= 9f0 ? sv816 : 0f0
            end
        end
    end
    if mtopp > d                                       # r5harv.f:410-412 top>DBH ⇒ no merch/board (every branch)
        vol2 = 0f0; vol4 = 0f0
    end
    vol2 = round(vol2, RoundNearestTiesAway)           # dvest.f:170 VOL(2) = ANINT(VOL(2))
    vol7 < 0f0 && (vol7 = 0f0)                         # volinit.f:559 IF(VOL(7).LT.0.0) VOL(7)=0.0
    return (vol1, max(vol4 + vol7, 0f0), max(vol2, 0f0))
end

# REAL*8 DATA from a default-REAL literal (gfortran: the constant is REAL(4), widened on assignment).
@inline _r8(x::Float64) = Float64(Float32(x))
# COFx(1)·DBHOB**COFx(2)·HTTOT**COFx(3)·IV**COFx(4) in REAL*8 (REAL**REAL*8 ⇒ glibc pow), stored to a REAL CVx.
@inline _r5_pw(c, D::Float32, H::Float32) =
    Float32(_r8(c[1]) * dpow(Float64(D), _r8(c[2])) * dpow(Float64(H), _r8(c[3])) * dpow(10.0, _r8(c[4])))
# gfortran LOG10 of a REAL → glibc log10f
const _r5_log10 = log10f   # glibc (FMath single libm binding)

# ---------------------------------------------------------------------------
# WO2W: R5TAP — Wensel & Krumland (Region-5 California) profile taper (r5tap.f).
# NOT Flewelling: WO2W routes through the generic PROFILE driver (VOLEQ(4:4)='W' ⇒ FWINIT skipped,
# JSP=0), and TAPERMODEL calls R5TAP for the inside-bark DIB at any height. TCUBIC (total cubic VOL1)
# and the GETDIB→Smalian merch loop (VOL4) are the SAME integration the FW2 profile helpers implement
# (profile.f is the shared source of both) ⇒ reuse _fw2_tcubic / _fw2_merch_cuft / _fw2_board with a
# dibat closure = nc_r5tap_dib. R5 merch defaults (mrules.f REGN 5): MTOPP=6, stump=1, MERCHL=8,
# MAXLEN=16, MINLEN=2, TRIM=0.5, EVOD=2, OPT=22 — the last five == _NVB_R3 / _cr_merch_*(iregn=4).
# ---------------------------------------------------------------------------

# R5WKC(5,9)=C1..C5, R5WKB(2,9)=B1,B2 (r5tap.f:16-55). VOLEQ(8:10) FIA → SP 1..9 (r5tap.f:94-116).
const NC_R5WKC = (
    (0.84292,0.97062,-0.38163,-0.0074002,0.0),          # 1 DF 202
    (0.87278,1.26066,-1.91214,0.020445,0.0),            # 2 PP 122
    (0.90051,0.91588,-0.92964,0.0077119,-0.0011019),    # 3 SP 117
    (0.86039,1.45196,-2.42273,-0.15848,0.036947),       # 4 WF 015  (TERM2 clamp ≥ -1)
    (0.87927,0.91350,-0.56617,-0.014480,0.0037262),     # 5 RF 020
    (1.0,0.31550,-0.34316,0.0,-0.00039283),             # 6 IC 081
    (0.82932,1.50831,-4.08016,0.047053,0.0),            # 7 JP 116
    (1.0,0.84257,-0.98434,0.0,0.0),                     # 8 LP 108
    (0.955,0.387,-0.362,-0.00581,0.00122),              # 9 RW 211
)
const NC_R5WKB = (
    (0.1420,0.04302),(0.1031,0.03068),(0.0743,0.02936),(0.0844,0.03320),(0.1105,0.05061),
    (0.1177,0.03894),(0.1472,0.03880),(0.0147,0.03223),(0.153,0.035),
)

"VOLEQ(8:10) FIA code → R5TAP SP index 1..9 (r5tap.f:94-116). 0 = unsupported."
@inline function _nc_r5tap_sp(fia::AbstractString)::Int
    fia == "202" ? 1 : fia == "122" ? 2 : fia == "117" ? 3 : fia == "015" ? 4 :
    fia == "020" ? 5 : fia == "081" ? 6 : fia == "116" ? 7 : fia == "108" ? 8 :
    fia == "211" ? 9 : 0
end

"R5TAP (r5tap.f:88-153): inside-bark DIB (in) at height `htup` on a stem of DBH `dbh`, total height
`totht`. REAL*8 kernel → Float32. `sp` = 1..9 (NC_R5WKC/B index). DM=0 (unused stem-form term)."
function nc_r5tap_dib(sp::Int, dbh::Float32, totht::Float32, htup::Float32)::Float32
    (sp < 1 || sp > 9) && return 0f0
    DBH = Float64(dbh); TOTHT = Float64(totht); HTUP = Float64(htup)
    HTUP > TOTHT && return 0f0
    c = NC_R5WKC[sp]; b = NC_R5WKB[sp]
    local dibcor::Float64
    if HTUP >= 4.499
        term1 = c[1]
        term2 = c[3] + c[4] * DBH + c[5] * TOTHT
        (sp == 4 && term2 > -1.0) && (term2 = -1.0)               # white-fir clamp (r5tap.f:130)
        term3 = ((HTUP - 1.0) / (TOTHT - 1.0))^c[2]
        term4 = log(1.0 - term3 * (1.0 - exp(c[1] / term2)))      # DM/(DBH·term2)=0
        dibcor = DBH * (term1 - term2 * term4)
    else
        dibcor = (1.0 - b[1]) * DBH * exp(b[2] * (4.5 - HTUP))
    end
    return Float32(dibcor)
end

"MERLEN (profile.f:1011-1052, non-Flewelling branch): merch length stump→`mtopp` top. Binary search in
0.1-ft steps for the highest height where the tenth-inch-truncated DIB ≥ mtopp (TOP1=ANINT((mtopp+.005)·10))."
@inline function nc_merlen(dibat, mtopp::Float32, httot::Float32, stump::Float32)::Float32
    top1 = round(Int, (mtopp + 0.005f0) * 10f0)              # ANINT
    first = 1; last = trunc(Int, httot + 0.5f0) * 10
    @inbounds while first != last
        half = (first + last + 1) ÷ 2
        dibt = trunc(Int, (dibat(Float32(half) / 10f0) + 0.005f0) * 10f0)
        top1 <= dibt ? (first = half) : (last = half - 1)
    end
    max(Float32(first) / 10f0 - stump, 0f0)
end

"WO2W merch cubic. `topwood=false` (the .sum path) returns primary VOL(4) [stump→MTOPP=6] only — the FVS
summary volume call runs with SPFLG=0 ⇒ VOL(7)=0 ⇒ MCF=VOL(4) (fvsvol.f:512). `topwood=true` also adds the
secondary VOL(7) [primary-top→MTOPS=4] (profile.f:684-819) for the harvest/product reports that set SPFLG=1.
Both products bucked in 1-inch DIB classes (GETDIB), each log 0.1-rounded via the exact MERLEN length; the
topwood butt = last primary log's small end, its logs restarting above LENMS. Measured: on nct01 the .sum
MCuFt matches VOL(4)-only (427 vs oracle 449, the residual = the TPA/expansion normalization also seen in
TCuFt −3.6%); VOL(4)+VOL(7) overshoots to 524 ⇒ the .sum is SPFLG=0."
function nc_wo2w_merch(dibat, h::Float32; mtopp::Float32 = 6f0, mtops::Float32 = 4f0, stump::Float32 = 1f0,
                       minlen::Float32 = 2f0, merchl::Float32 = 8f0, topwood::Bool = false)::Float32
    # --- primary product, stump → MTOPP (6") ---
    lmp = nc_merlen(dibat, mtopp, h, stump)
    lmp < merchl && return 0f0
    nsp = _nvb_numlog(22, 2, lmp, 16f0, minlen, 0.5f0)
    nsp == 0 && return 0f0
    loglen_p, nsp = _nvb_segmnt(22, 2, lmp, 16f0, minlen, 0.5f0, nsp)
    dibl = _fw2_dclass(dibat(4.5f0)); ht2 = stump; vol4 = 0f0; lenms = 0f0
    @inbounds for i in 1:nsp
        ht2 += 0.5f0 + loglen_p[i]; lenms += 0.5f0 + loglen_p[i]
        dib = dibat(ht2); (i == nsp && dib < mtopp) && (dib = mtopp)
        dibs = _fw2_dclass(dib)
        vol4 += floor(0.00272708f0 * (dibl * dibl + dibs * dibs) * loglen_p[i] * 10f0 + 0.5f0) / 10f0
        dibl = dibs
    end
    topwood || return vol4
    # --- secondary topwood, MTOPP → MTOPS (4"); SPFLG=1 harvest-report path only ---
    tw = nc_merlen(dibat, mtops, h, stump) - lenms
    vol7 = 0f0
    if tw >= minlen
        nst = _nvb_numlog(22, 2, tw, 16f0, minlen, 0.5f0)
        if nst > 0
            loglen_t, nst = _nvb_segmnt(22, 2, tw, 16f0, minlen, 0.5f0, nst)
            dblt = dibl; ht2 = stump + lenms                 # butt = last primary small end
            @inbounds for i in 1:nst
                ht2 += 0.5f0 + loglen_t[i]
                dib = dibat(ht2); (i == nst && dib < mtops) && (dib = mtops)
                dibs = _fw2_dclass(dib)
                vol7 += floor(0.00272708f0 * (dblt * dblt + dibs * dibs) * loglen_t[i] * 10f0 + 0.5f0) / 10f0
                dblt = dibs
            end
        end
    end
    return vol4 + vol7
end

"WO2W per-tree volume (R5TAP profile). Returns (tcuft VOL1, merch-cuft VOL4, scribner VOL2).
★ VOL(1) TOTAL CUBIC is BIT-EXACT vs live (validated per-tree vs a standalone R5TAP driver + the FVSnc
treelist: SP D11.5 17.6==17.6, DF D12.7 20.9==20.9, big SP D34.6 277.4≈277.3). Total cubic integrates the
raw DIB (TCUBIC stump-cyl + 4-ft Smalian + tip); the R5TAP DIBs match live to 4 decimals (9.9332/8.4859/…).
MERCH (VOL4) here = primary product to MTOPP=6\" IB (MERLEN length == live). ⚠ The FVS .sum 'Merch CuFt'
is MCF = VOL(4)+VOL(7) [primary-to-6\" + topwood 6\"→4\"], gated on D≥DBHMIN(ISPC), and board uses
BFTOPD·BARK tops (fvsvol.f:337-512). VOL(7) topwood + the DBHMIN gate + the bark-adjusted board recompute
are the remaining merch-driver pieces (measured, not yet ported) — see docs/NC_VARIANT_PORT_AUDIT.md.
Needs CUTFLG/CUPFLG=1 (NC .sum ⇒ on). ERRFLAG paths: d<1 or h<5 ⇒ 0."
function nc_wo2w_vol(voleq::AbstractString, d::Float32, h::Float32; mtopp::Float32 = 6.0f0,
                     bftop::Float32 = mtopp, ht2td = nothing)
    # `ht2td` (optional 2-slot buffer) ← [HT1PRD cubic call, HT1PRD board call]: profile.f MERLEN non-Flewelling
    # branch, HT1PRD = LMERCH + STUMP with LMERCH = FIRST/10 − STUMP floored at 0 (profile.f:1011-1052, :342).
    ht2td === nothing || (ht2td[1] = 0f0; ht2td[2] = 0f0)
    length(voleq) < 10 && return (0f0, 0f0, 0f0)
    (d < 1f0 || h < 5f0) && return (0f0, 0f0, 0f0)
    sp = _nc_r5tap_sp(voleq[8:10])
    sp == 0 && return (0f0, 0f0, 0f0)
    dibat = ht -> nc_r5tap_dib(sp, d, h, Float32(ht))
    stump = 1.0f0; minl = 2.0f0; merl = 8.0f0                   # mrules.f REGN 5 defaults
    tcf = _nint(_fw2_tcubic(dibat, h) * 10.0f0) * 1f-1          # profile.f:293 VOL(1)=NINT(TCVOL*10.0)*1E-1 (×REAL 0.1, not /10)
    mcf = nc_wo2w_merch(dibat, h; mtopp = mtopp, stump = stump, minlen = minl, merchl = merl)  # VOL(4) (SPFLG=0)
    # The board pass's merch height is profile.f MERLEN (tenth-inch-truncated 0.1-ft search, nc_merlen), exactly
    # as the cubic pass uses — NOT _fw2_hs's diameter-tolerance bisection, which lands a hair off (96.97 vs 97.0,
    # 33.08 vs 32.9 ft) and shifts a 16-ft log boundary ⇒ single-tree ±10-bf Scribner steps. Verified against the
    # standalone VOLINITNVB oracle driver on WS/NC trees (SP 450, DF 470, OS 40 — the old path gave 440/460/50).
    bf  = _fw2_board(dibat, h, bftop, stump, minl, merl;
                     hs_solver = top -> nc_merlen(dibat, top, h, stump) + stump)   # VOL(2) Scribner
    if ht2td !== nothing
        ht2td[1] = nc_merlen(dibat, mtopp, h, stump) + stump
        ht2td[2] = bftop == mtopp ? ht2td[1] : nc_merlen(dibat, bftop, h, stump) + stump
    end
    return (max(tcf, 0f0), max(mcf, 0f0), max(bf, 0f0))
end

# ---------------------------------------------------------------------------
# NC REGION-6 (Siskiyou IFOR=4) volume — NVEL VOLEQDEF R6 branch, NOT R5.
#
# ROOT CAUSE of the large-tree WO2W total-cubic over-prediction: NC's VEQNNC is forest/region-dependent
# (sitset.f → VOLEQDEF(VAR='NC',IREGN=KODFOR/100,FORST)). The 500WO2W/500DVEW table above is the REGION-5
# assignment (forests Klamath/SixRivers/Trinity, IFOR 1-3). SISKIYOU (IFOR=4, forest 611, a Region-6 forest —
# and its reservation crosswalks 8103/8105) instead gets the REGION-6 equations: westside Flewelling
# F06FW2W202 (DF), INGY I00FW2W093/073 (WF/PP), and Region-6 Behre 616BEHW<fia> for everything else.
# jl previously hardcoded the R5 table for ALL forests, so on Siskiyou stands large sound DF was computed with
# the R5 Wensel-Krumland R5TAP taper (a much fatter stem) instead of the R6 Flewelling profile — +30% TCuFt at
# D=54. VEQNNC table confirmed bit-exact vs FVSnc_g16 forest 611 (.out VOLEQ table); per-tree TVOL1 bit-exact
# (D=54.4/H=221 DF: 827.4 == 827.4 vs the R5TAP 1075.3). The R6 kernels are shared NVEL and already ported for
# WC (westside Flewelling wc_fw2_westside_vol / f_west.f SHP_W3) and BM/WC (Behre bm_r6vol3/r6dibs/r6vol1).
# ---------------------------------------------------------------------------
const NC_R6_VOL_EQ = String[
    "616BEHW299",  # 1  OS  → Behre (other softwood → 299)
    "616BEHW117",  # 2  SP  → Behre
    "F06FW2W202",  # 3  DF  → westside Flewelling (F-model)
    "I00FW2W093",  # 4  WF  → INGY (Engelmann-spruce profile 093)
    "616BEHW361",  # 5  MA  → Behre
    "616BEHW081",  # 6  IC  → Behre
    "616BEHW818",  # 7  BO  → Behre
    "616BEHW631",  # 8  TO  → Behre
    "616BEHW020",  # 9  RF  → Behre
    "I00FW2W073",  # 10 PP  → INGY (western-larch profile 073)
    "616BEHW998",  # 11 OH  → Behre
    "616BEHW211",  # 12 RW  → Behre
]

# SISKFC(sp, IFCDBH) form class (formcl.f:35-40, IFOR 4). Column = DBH class (1..5), row = species 1..12.
# IFCDBH = INT((D-1)/10+1), clamp ≥1, D>40.9 → 5. Used only by the Region-6 Behre path.
const NC_SISKFC = (
    (91f0, 84f0, 79f0, 78f0, 78f0),  # 1  OS
    (96f0, 91f0, 85f0, 83f0, 82f0),  # 2  SP
    (90f0, 86f0, 81f0, 80f0, 80f0),  # 3  DF
    (98f0, 90f0, 86f0, 85f0, 85f0),  # 4  WF
    (98f0, 88f0, 84f0, 81f0, 80f0),  # 5  MA
    (89f0, 89f0, 77f0, 73f0, 72f0),  # 6  IC
    (98f0, 98f0, 98f0, 98f0, 98f0),  # 7  BO
    (91f0, 91f0, 82f0, 80f0, 79f0),  # 8  TO
    (92f0, 83f0, 80f0, 80f0, 79f0),  # 9  RF
    (93f0, 89f0, 83f0, 81f0, 80f0),  # 10 PP
    (95f0, 86f0, 78f0, 76f0, 75f0),  # 11 OH
    (82f0, 82f0, 79f0, 78f0, 78f0),  # 12 RW
)

@inline function nc_siskfc(sp::Int, d::Float32)::Int
    (sp < 1 || sp > 12) && return 80
    ifc = Int(floor((d - 1f0) / 10f0 + 1f0))
    ifc < 1 && (ifc = 1)
    d > 40.9f0 && (ifc = 5)
    Int(NC_SISKFC[sp][ifc])
end

# NC Region-6 Behre (616BEHW) per-tree volume. Total cubic via R6VOL3 (fvsvol.f→r6vol.f:100-114, ZONE=1,
# FC_HT=17.3 short-tree cylinder); merch VOL(4)/board VOL(2) via R6DIBS/R6VOL1 with NC TOPD=6.0 (·BARK IB top,
# sitset.f DEFAULT). Reuses the shared bm_r6vol3/r6dibs/r6vol1 kernels + NC SISKFC form class.
# ⚠ The cylinder cutoff is r6vol.f's FC_HT=17.3 (ZONE 1), NOT profile.f/profile2.f's 16.3 — NC's total cubic is
# driven by R6VOL (r6vol.f), not PROFILE. A tree with TTH∈(16.3,17.3] is SHORTER than the 17.3-ft butt log, so
# routing it through the full R6VOL3 profile makes HTUP=TTH−17.3<0 and the taper loop runs away (verified: a
# D=4.5/H=17 chinquapin gave 26.7 cuft vs the oracle's cylinder 0.76). WC/BM already guard at 17.3; this was an
# NC-only mis-port to the profile2 cutoff. (fixed 2026-09-05: NC Siskiyou short-tree TCuFt over-prediction bug.)
function nc_behre_vol(sp::Int, d::Float32, h::Float32, bark::Float32; topd::Float32 = 6.0f0)
    fc = nc_siskfc(sp, d)
    dbtbh = d * (1f0 - bark); dbhib = d - dbtbh
    vol2 = 0f0; vol4 = 0f0
    v1 = if h <= 17.3f0
        0.00272708f0 * (dbhib * dbhib) * h            # r6vol.f:105 short-tree cylinder (TTH ≤ FC_HT=17.3, ZONE 1)
    else
        v = bm_r6vol3(d, dbtbh, fc, h, 1)
        mtopp = topd * bark                          # TOPD·BARK inside-bark top (Siskiyou TOPD=4.5)
        xlogs, ld1 = bm_r6dibs(d, fc, mtopp, h)
        lv1, lv4 = bm_r6vol1(d, fc, xlogs, ld1)
        nlog = Int(floor(xlogs)); nacc = (xlogs - nlog) > 0f0 ? nlog + 1 : nlog
        for k in 1:nacc
            vol2 += bm_anint(lv1[k])
            vol4 += bm_anint(lv4[k] * 10f0) / 10f0
        end
        v
    end
    return (max(v1, 0f0), max(vol4, 0f0), max(vol2, 0f0))
end

# ---------------------------------------------------------------------------
# NC REGION-7 (HOOPA IFOR=5, forest 705) volume — NVEL BLMVOL (Behre-hyperbola BLM taper, blmtap.f/blmvol.f).
#
# ROOT CAUSE of the +16-22% TCuFt / −8-23% BdFt on the 57 LOC-705 stands: NC's VEQNNC is forest/region-
# dependent (sitset.f → VOLEQDEF(VAR='NC',IREGN=KODFOR/100,FORST)). Forest 705 = HOOPA, IREGN=7, FORST='05'
# ⇒ VOLEQDEF returns B00BEHW<fia> (DF = B01BEHW202) — the BLM Oregon volume routines (volinit.f:366 → BLMVOL
# → blmtap.f BEHRE'S HYPERBOLA), NOT the Region-5 WO2W/DVE table (IFOR 1-3) nor the Region-6 Behre (IFOR 4).
# jl previously had no IFOR-5 branch ⇒ fell through to the R5 WO2W taper (fatter stem, R5 board rules) ⇒ the
# one-directional cyc0 volume error. FORMCL (formcl.f) leaves FC=80 for IFOR 5 (neither SISKFC-4 nor BLM712-7).
# Merch specs for IFOR 5,7 (sitset.f CASE(5,7)): DBHMIN=9.0, TOPD=BFTOPD=5.0 (inside-bark top = 5.0·BARK).
#
# The driver (fvsvol.f NATCRS, IREGN=7) calls BLMVOL twice: (1) cubic — MTOPP=TOPD·BARK=5·BARK, BFPFLG=0,
# CUPFLG=1 ⇒ VOL(1) total cubic (BLMTCUB) + VOL(4) merch cubic; (2) board — MTOPP=BFTOPD·BARK=5·BARK,
# BFPFLG=1 ⇒ VOL(2) Scribner. Both MTOPP identical (TOPD==BFTOPD==5) ⇒ one pass yields VOL(1)/VOL(4)/VOL(2).
# STUMP is hard-wired to 1.0 inside BLMVOL. Broken tops trimmed by CFTOPK/BFTOPK (r4_topkill), METHC=6.
# The BLMVOL kernel itself is the shared engine `blm_vol` (src/engine/blm_vol.jl).
# ---------------------------------------------------------------------------

# VEQNNC/VEQNNB for forest 705 (VOLEQDEF VAR=NC IREGN=7 FORST=05). CONFIRMED bit-exact vs FVSnc_g16 705 .out.
const NC_R7_VOL_EQ = String[
    "B00BEHW999",  # 1  OS
    "B00BEHW117",  # 2  SP
    "B01BEHW202",  # 3  DF
    "B00BEHW015",  # 4  WF
    "B00BEHW361",  # 5  MA
    "B00BEHW081",  # 6  IC
    "B00BEHW800",  # 7  BO
    "B00BEHW631",  # 8  TO
    "B00BEHW021",  # 9  RF
    "B00BEHW122",  # 10 PP
    "B00BEHW999",  # 11 OH
    "B00BEHW211",  # 12 RW
]

"BLMVOL (blmvol.f) per-tree BLM Oregon volume for forest 705 (Behre-hyperbola taper, FC=80). Returns
(tcuft VOL(1), merch-cuft VOL(4), scribner-bf VOL(2)). `bark`=DIB/DOB (nc_bratio) ⇒ MTOPP=5·BARK inside-bark
merch/board top. The FULL (unbroken NORMHT) volumes; the driver applies r4_topkill (CFTOPK/BFTOPK) for
broken tops. VOL(1) is UNROUNDED (BLMTCUB, unlike profile.f TCUBIC's NINT·10/10). Small-tree cylinder for
TTH≤17.8 or when the DIB at 17.3 ft < MTOPP."
# nc/formcl.f BLM712 — Coos Bay BLM (IFOR 7) per-species form class; every other non-Siskiyou forest = 80.
const NC_BLM712_FC = Float32[74, 76, 76, 78, 72, 66, 72, 74, 78, 80, 70, 75]

function nc_blmvol_vol(sp::Int, d::Float32, h::Float32, bark::Float32; fclass::Float32 = 80f0, ifor::Int = 5,
                       ht_out = nothing)
    ht_out === nothing || (ht_out[] = 0f0)
    (d < 1f0 || h <= 0f0) && return (0f0, 0f0, 0f0)    # DBHOB<=0 (errflag3) / HTTOT<=0 (errflag4) ⇒ all zero
    # 712 Coos Bay (IFOR 7): VOLEQDEF gives DF B02BEHW202 (705 Hoopa: B01) ⇒ BLMTAPEQ PROFILE 2 / TAPEQU 2.
    eq = (ifor == 7 && sp == 3) ? "B02BEHW202" : NC_R7_VOL_EQ[sp]
    # MTOPP = TOPD(=5.0)·BARK for both calls (TOPD==BFTOPD) ⇒ one BLMVOL pass yields VOL(1)/VOL(4)/VOL(2)
    # (and the one HT1PRD, blmvol.f:427, both calls return).
    v1, v2, v4 = blm_vol(eq, 5.0f0 * bark, h, d, trunc(Int, fclass); bfpflg = true, ht_out = ht_out)
    return (max(v1, 0f0), max(v4, 0f0), max(v2, 0f0))
end

# ---------------------------------------------------------------------------
# Driver — dispatch each species' VEQNNC by model type (WO2W vs DVEW).
# ---------------------------------------------------------------------------
function compute_volumes_nc!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; sd = s.coef.species; c = s.control
    ncmerch = (stmp = c.sp_stump_ht, topd = c.sp_top_diam, scfstmp = c.sp_scf_stump,
               scftop = c.sp_scf_topd, bftopd = c.sp_bf_topd, bfstmp = c.sp_bf_stump)
    isr6 = Int(s.plot.forest_idx) == 4       # SISKIYOU (IFOR 4, forest 611) = Region-6 VOLEQDEF branch
    isr7 = Int(s.plot.forest_idx) in (5, 7)  # Region 7 (705 Hoopa IFOR 5, 712 BLM Coos Bay IFOR 7) = BLMVOL branch (voleqdef R7 table)
    # vols.f:86-90 zeroes HT2TD for every record; NATCRS (fvsvol.f:337-339 / 484-487) stores the NVEL HT1PRD
    # (BFMIND = DBHMIN = 9 for NC, nc/grinit.f).
    fill!(t.merch_top_cf, 0f0); fill!(t.merch_top_bf, 0f0)
    htb = zeros(Float32, 2); ht7 = Ref(0f0)
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        htb[1] = 0f0; htb[2] = 0f0
        if d < 1f0 || sp < 1 || sp > 12
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        dbhmin0 = 9.0f0
        hv0 = (t.trunc[i] > 0 && t.norm_ht[i] > 0) ? Float32(t.norm_ht[i]) / 100f0 : h
        # vols.f:150 BARK=BRATIO(ISPC,D,H) is taken at the START-of-cycle DBH (before D=D+DG/BARK) — the stashed
        # vol_bark; the grown-DBH bark only at cycle 0 / for dead records. It sets the merch/board tops (TOPD·BARK).
        bark = (i <= t.n && t.vol_bark[i] > 0f0) ? t.vol_bark[i] : nc_bratio(sd[:bark1][sp], sd[:bark2][sp], Int(sd[:bark_imap][sp]), d)
        if isr7
            # Region-7 (Hoopa 705): NVEL BLMVOL (Behre-hyperbola BLM taper), FC=80, TOPD=5.0.
            tcf7, mcf7, bf7 = nc_blmvol_vol(sp, d, hv0, bark;
                                            fclass = Int(s.plot.forest_idx) == 7 ? NC_BLM712_FC[sp] : 80f0, ifor = Int(s.plot.forest_idx),
                                            ht_out = ht7)
            d >= dbhmin0 && (t.merch_top_cf[i] = ht7[]; t.merch_top_bf[i] = ht7[])
            tcf7, mcf7, bf7 = r4_topkill(t, i, sp, d, hv0, bark, tcf7, mcf7, bf7, ncmerch, _R4_TOPD5)
            t.cuft_vol[i] = max(tcf7, 0f0)
            t.merch_cuft_vol[i] = d >= dbhmin0 ? max(mcf7, 0f0) : 0f0
            t.saw_cuft_vol[i] = 0f0
            t.bdft_vol[i] = d >= dbhmin0 ? max(bf7, 0f0) : 0f0
            continue
        end
        if isr6
            # Region-6 (Siskiyou): F06 westside Flewelling (DF) / I00 INGY (WF,PP) / 616 Behre (rest).
            eq6 = NC_R6_VOL_EQ[sp]; m6 = eq6[4:6]
            local tcf6::Float32, mcf6::Float32, bf6::Float32
            if m6 == "FW2" && (eq6[1] == 'F' || eq6[1] == 'f')
                tcf6, mcf6, bf6 = wc_fw2_westside_vol(eq6, d, hv0, bark; topd = c.sp_top_diam[sp], bftopd = c.sp_bf_topd[sp],
                                                      ht2td = htb)
            elseif m6 == "FW2"                              # INGY (I00) — reuse cr_fw2_vol
                v = cr_fw2_vol(eq6, d, hv0; bark = bark, topd = c.sp_top_diam[sp], bftopd = c.sp_bf_topd[sp], stump = 1f0, iregn = 6,
                               board_cor = 'N', merch_opt = 23, sf_hs = true, ht2td = htb)
                tcf6 = max(v[1], 0f0); mcf6 = max(v[4] + v[7], 0f0); bf6 = max(v[2], 0f0)
            else                                            # 616BEHW Behre
                tcf6, mcf6, bf6 = nc_behre_vol(sp, d, hv0, bark; topd = c.sp_top_diam[sp])
                htb[1] = r6vol_ht1prd(d, nc_siskfc(sp, d), c.sp_top_diam[sp] * bark, hv0, d - d * (1f0 - bark)); htb[2] = htb[1]
            end
            d >= dbhmin0 && (t.merch_top_cf[i] = htb[1]; t.merch_top_bf[i] = htb[2])
            tcf6, mcf6, bf6 = r4_topkill(t, i, sp, d, hv0, bark, tcf6, mcf6, bf6, ncmerch, _BM_TOPD45)
            t.cuft_vol[i] = max(tcf6, 0f0)
            t.merch_cuft_vol[i] = d >= dbhmin0 ? max(mcf6, 0f0) : 0f0
            t.saw_cuft_vol[i] = 0f0
            t.bdft_vol[i] = d >= dbhmin0 ? max(bf6, 0f0) : 0f0
            continue
        end
        eq = NC_VOL_EQ[sp]; mdl = eq[4:6]
        # nc/sitset.f:196-224 forest-default merch specs (IFOR 1 = Klamath 505 = DEFAULT case):
        # DBHMIN=9.0, TOPD=BFTOPD=6.0. Merch/board are ZEROED for D < DBHMIN (fvsvol.f:337,512 gate).
        dbhmin = 9.0f0
        # Top-killed: the full cubic uses the NORMAL height (norm_ht, cratet.f nc_htdbh_h SISKIY dub),
        # then CFTOPK/BFTOPK (vols.f:193) trims it back to the standing break (t.trunc/100) via the Behre
        # taper. WITHOUT the trim, a broken-top redwood's full-NORMHT cubic over-counts; WITHOUT the NORMHT
        # (using the recorded broken height) it under-counts 16-31%. Both together match the oracle.
        hv = (t.trunc[i] > 0 && t.norm_ht[i] > 0) ? Float32(t.norm_ht[i]) / 100f0 : h
        # fvsvol.f (the oracle buildDir) passes western tops INSIDE bark: cubic TOPDIAM=TOPD·BARK, board
        # MTOPP=BFTOPD·BARK. With a fixed 6" IB top jl lost ~11% MCuFt / ~20% BdFt (measured per tree vs live
        # FVSnc_g16 DEBUG FVSVOL: PP D11.5 TOPDIAM 5.148 = 6·0.858).
        topib = c.sp_top_diam[sp] * bark; bftib = c.sp_bf_topd[sp] * bark
        if mdl == "WO2"
            tcf, mcf, bf = nc_wo2w_vol(eq, d, hv; mtopp = topib, bftop = bftib, ht2td = htb)
            d >= dbhmin && (t.merch_top_cf[i] = htb[1]; t.merch_top_bf[i] = htb[2])
            # Broken/killed-top reduction (fvsvol.f CFTOPK/BFTOPK) — WO2W conifer/redwood path only; the DVE
            # California-hardwood path skips CFTOPK (fvsvol.f). NC merch TOPD=6.0 (nc/sitset.f DEFAULT).
            tcf, mcf, bf = r4_topkill(t, i, sp, d, hv, bark, tcf, mcf, bf, ncmerch, _R4_TOPD6)
            t.cuft_vol[i] = tcf
            t.merch_cuft_vol[i] = d >= dbhmin ? mcf : 0f0
            t.saw_cuft_vol[i] = 0f0
            t.bdft_vol[i] = d >= dbhmin ? bf : 0f0
        else                                       # DVE — California hardwood D²H (r5harv.f), MTOPP=6
            # Full (unbroken-height NORMHT) cubic/merch/board, then CFTOPK/BFTOPK trim for broken tops:
            # fvsvol.f METHC=6 → NATCRS sets VMAX=TCF & CTKFLG, then vols.f:193/391 CFTOPK/BFTOPK trim TCF/MCF
            # & BDFT — the SAME top-kill path as the WO2W conifers (both are NVEL method 6). DVE is NOT skipped.
            tcf, mcf, bf = nc_r5harv_vol(eq, d, hv, topib)                       # CF pass, TOPD·BARK
            bftib != topib && (bf = nc_r5harv_vol(eq, d, hv, bftib)[3])           # BF pass, BFTOPD·BARK
            tcf, mcf, bf = r4_topkill(t, i, sp, d, hv, bark, tcf, mcf, bf, ncmerch, _R4_TOPD6)
            t.cuft_vol[i] = tcf
            t.merch_cuft_vol[i] = d >= dbhmin ? mcf : 0f0
            t.saw_cuft_vol[i] = 0f0
            t.bdft_vol[i] = d >= dbhmin ? bf : 0f0
        end
    end
    return s
end

# ---------------------------------------------------------------------------
# Shared Region-5 NVEL per-tree volume for the western variants' R5 forests (voleqdef.f R5_EQN: the equation
# depends only on species + variant, not the forest). fvsvol.f passes the tops INSIDE bark (TOPD·BARK /
# BFTOPD·BARK, BARK = start-of-cycle BRATIO) and mrules.f REGN 5 (stump 1, MINLEN 2, MERCHL 8, COR 'Y',
# OPT 22). Returns (tcuft VOL1, merch cuft VOL4 (SPFLG=0), scribner VOL2) before the CFTOPK/BFTOPK trim.
function nvel_r5_vol(eq::AbstractString, d::Float32, hv::Float32, bark::Float32, topd::Float32, bftopd::Float32;
                     ht2td = nothing)
    # `ht2td` ← [cubic, board] HT1PRD (WO2W MERLEN, INGY SF_HS); R5HARV (DVE) sets none ⇒ 0.
    ht2td === nothing || (ht2td[1] = 0f0; ht2td[2] = 0f0)
    m = eq[4:6]
    if m == "WO2"
        return nc_wo2w_vol(eq, d, hv; mtopp = topd * bark, bftop = bftopd * bark, ht2td = ht2td)
    elseif m == "DVE"
        tcf, mcf, bf = nc_r5harv_vol(eq, d, hv, topd * bark)
        bftopd != topd && (bf = nc_r5harv_vol(eq, d, hv, bftopd * bark)[3])
        return (tcf, mcf, bf)
    else                                            # INGY FW2 (e.g. GF I15FW2W017) — cr_fw2_vol applies ·bark itself
        v = cr_fw2_vol(eq, d, hv; bark = bark, topd = topd, bftopd = bftopd, stump = 1f0,
                       iregn = 5, board_cor = 'Y', merch_opt = 22, sf_hs = true, ht2td = ht2td)
        return (max(v[1], 0f0), max(v[4], 0f0), max(v[2], 0f0))
    end
end

# ---------------------------------------------------------------------------
# nc_snag_bole_cuft — the FFE snag-bole TOTAL cubic (FMSVOL→TCF) for NC. NC's `vol_eq` is EMPTY (it uses
# the NVEL WO2W/DVE models in NC_VOL_EQ, not an R8-Clark string), so the shared _R8CLARK_VOL snag path
# returns 0 ⇒ every NC snag bole collapses to the tiny cone floor ⇒ snag falldown adds ~nothing to the
# >3" down-wood pool (the pool shrinks instead of growing). Mirror `cr_snag_bole_cuft`: return the total
# cubic VOL(1) from NC's own volume model (FVS FMSVOL fmsvol.f:153 VOL2HT=MAX(X,TCF) for the western
# NVEL variants ⇒ bole==fall==TCF). Used by _snag_merch_cuft_on + the input/SNAGINIT/fire snag paths.
function nc_snag_bole_cuft(s::StandState, sp::Int, d::Float32, h::Float32)::Float32
    (d < 1f0 || h <= 0f0 || sp < 1 || sp > 12) && return 0f0
    if Int(s.plot.forest_idx) in (5, 7)            # HOOPA (Region-7 705): BLMVOL Behre total cubic VOL(1)
        bark = nc_bratio(s.coef.species[:bark1][sp], s.coef.species[:bark2][sp],
                         Int(s.coef.species[:bark_imap][sp]), d)
        return max(nc_blmvol_vol(sp, d, h, bark;
                                 fclass = Int(s.plot.forest_idx) == 7 ? NC_BLM712_FC[sp] : 80f0, ifor = Int(s.plot.forest_idx))[1], 0f0)
    end
    if Int(s.plot.forest_idx) == 4              # SISKIYOU (Region-6): same VEQNNC dispatch as compute_volumes_nc!
        eq6 = NC_R6_VOL_EQ[sp]; m6 = eq6[4:6]
        bark = nc_bratio(s.coef.species[:bark1][sp], s.coef.species[:bark2][sp],
                         Int(s.coef.species[:bark_imap][sp]), d)
        if m6 == "FW2" && (eq6[1] == 'F' || eq6[1] == 'f')
            return max(wc_fw2_westside_vol(eq6, d, h, bark; topd = 6.0f0, bftopd = 6.0f0)[1], 0f0)
        elseif m6 == "FW2"
            return max(cr_fw2_vol(eq6, d, h; bark = bark, topd = 6.0f0, bftopd = 6.0f0, stump = 1f0, iregn = 6, board_cor = 'N', merch_opt = 23)[1], 0f0)
        else
            return max(nc_behre_vol(sp, d, h, bark)[1], 0f0)
        end
    end
    eq = NC_VOL_EQ[sp]
    length(eq) < 6 && return 0f0
    v = eq[4:6] == "WO2" ? nc_wo2w_vol(eq, d, h) : nc_r5harv_vol(eq, d, h, 6.0f0)
    return max(v[1], 0f0)                       # VOL(1) total cubic
end
