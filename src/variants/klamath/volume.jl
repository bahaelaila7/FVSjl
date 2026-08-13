# =============================================================================
# volume.jl (klamath) — NC volume (chunk 8). NVEL Region-5/6, dispatched by VEQNNC geocode.
#
# NC's 12 species (VOLEQ confirmed BIT-EXACT vs live oracle nct01.out:68-70, from voleqdef.f R5_EQN):
#   8 conifers + redwood = "500WO2W<fia>"  → Region-6 West-side Flewelling taper (r6vol.f)  [nc_wo2w_vol]
#   4 hardwoods MA/BO/TO/OH = "500DVEW<fia>" → Region-5 DVE California hardwood D²H (r5harv.f) [nc_r5harv_vol]
# nct01 = 100% conifers (SP/DF/WF/RF) ⇒ only the WO2W path is exercised for nct01 validation; the DVEW
# hardwood path below is ported faithfully from r5harv.f but is NOT reachable on nct01 (no hardwood stand).
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

"R5HARV DVE cubic (r5harv.f). Returns (tcuft VOL(1), merch-cuft VOL(4), topwood VOL(7)). `mtopp` = merch
top DIB (TOPC): [3,5)→CV4, [5,7)→CV6, [7,9]→CV8, <3→CVT (whole stem). Misc-hardwood power-law branch
(NC's MA/BO/TO/OH); juniper & giant-sequoia special branches included. Board/Intl (VOL(2)/VOL(10)) omitted
(NC .sum reports cubic only for DVE; board handled by the driver if ever needed)."
function nc_r5harv_vol(voleq::AbstractString, d::Float32, h::Float32, mtopp::Float32)
    length(voleq) < 10 && return (0f0, 0f0, 0f0)
    d < 1f0 && return (0f0, 0f0, 0f0)
    spec = _nc_r5_spec(voleq[8:10])
    spec < 0 && return (0f0, 0f0, 0f0)
    D = Float64(d); H = Float64(h); TOPC = Float64(mtopp); IV = 10.0
    if spec == 0                                       # juniper (060/064)
        if D < 5 || H < 10
            cvts = 0.00272708 * D * D * H; v = 0.0
        else
            f = 0.307 + 0.00086*H - 0.0037*D*H/(H - 4.5)
            ba = 0.005454154*D*D
            cvts = ba*f*H*(H/(H - 4.5))^2
            v = TOPC > 0 ? (cvts + 3.48)/(1.18052 + 0.32736*exp(-0.1*D)) - 2.948 : cvts
        end
        vol4 = round(v*10 + 0.5)/10.0
        return (Float32(cvts), Float32(max(vol4,0.0)), 0f0)
    elseif spec == 14                                  # giant sequoia
        vol4 = NC_R5_SEQC[1]*D^NC_R5_SEQC[2]*H^NC_R5_SEQC[3]
        return (Float32(vol4), Float32(vol4), 0f0)     # VOL(1)=VOL(4)
    end
    # SPEC 1 red alder shares the misc formula's shape via its own TARIF branch in Fortran; NC never uses
    # 351, so the misc power-law covers every NC DVE species (SPEC 3/9/10/12). Faithful misc branch:
    a = NC_R5_COFA[spec]; b = NC_R5_COFB[spec]; cc = NC_R5_COFC[spec]
    cv4 = a[1]*D^a[2]*H^a[3]*IV^a[4]
    cv8 = b[1]*D^b[2]*H^b[3]*IV^b[4]
    cv6 = (cv4 > 0 && cv8 > 0) ? cv4 - (cv4 - cv8)*0.4 : 0.0
    cvt = cc[1]*D^cc[2]*H^cc[3]*IV^cc[4]
    cuftgros = TOPC >= 3 && TOPC < 5 ? cv4 :
               TOPC >= 5 && TOPC < 7 ? cv6 :
               TOPC >= 7 && TOPC <= 9 ? cv8 :
               TOPC < 3 ? cvt : 0.0
    topwood = cv4 - cuftgros
    mtopp > d && (cuftgros = 0.0)                      # r5harv.f:412 top>DBH ⇒ no merch
    return (Float32(cvt), Float32(max(cuftgros,0.0)), Float32(max(topwood,0.0)))
end

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
function nc_wo2w_vol(voleq::AbstractString, d::Float32, h::Float32)
    length(voleq) < 10 && return (0f0, 0f0, 0f0)
    (d < 1f0 || h < 5f0) && return (0f0, 0f0, 0f0)
    sp = _nc_r5tap_sp(voleq[8:10])
    sp == 0 && return (0f0, 0f0, 0f0)
    dibat = ht -> nc_r5tap_dib(sp, d, h, Float32(ht))
    mtopp = 6.0f0; stump = 1.0f0; minl = 2.0f0; merl = 8.0f0    # mrules.f REGN 5 defaults
    tcf = Float32(round(_fw2_tcubic(dibat, h) * 10.0)) / 10.0f0  # VOL(1)=NINT(TCVOL*10)/10
    mcf = nc_wo2w_merch(dibat, h; mtopp = mtopp, stump = stump, minlen = minl, merchl = merl)  # VOL4+VOL7
    bf  = _fw2_board(dibat, h, mtopp, stump, minl, merl)         # VOL(2) Scribner (primary to 6")
    return (max(tcf, 0f0), max(mcf, 0f0), max(bf, 0f0))
end

# ---------------------------------------------------------------------------
# Driver — dispatch each species' VEQNNC by model type (WO2W vs DVEW).
# ---------------------------------------------------------------------------
function compute_volumes_nc!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0 || sp < 1 || sp > 12
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = NC_VOL_EQ[sp]; mdl = eq[4:6]
        # nc/sitset.f:196-224 forest-default merch specs (IFOR 1 = Klamath 505 = DEFAULT case):
        # DBHMIN=9.0, TOPD=BFTOPD=6.0. Merch/board are ZEROED for D < DBHMIN (fvsvol.f:337,512 gate).
        dbhmin = 9.0f0
        # Top-killed: full cubic uses the NORMAL height, then the profile naturally truncates at the break.
        hv = (t.trunc[i] > 0 && t.norm_ht[i] > 0) ? Float32(t.norm_ht[i]) / 100f0 : h
        if mdl == "WO2"
            tcf, mcf, bf = nc_wo2w_vol(eq, d, hv)
            t.cuft_vol[i] = tcf
            t.merch_cuft_vol[i] = d >= dbhmin ? mcf : 0f0
            t.saw_cuft_vol[i] = 0f0
            t.bdft_vol[i] = d >= dbhmin ? bf : 0f0
        else                                       # DVE — California hardwood D²H (r5harv.f), MTOPP=6
            tcf, mcf, _ = nc_r5harv_vol(eq, d, hv, 6.0f0)
            t.cuft_vol[i] = tcf
            t.merch_cuft_vol[i] = d >= dbhmin ? mcf : 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0   # DVE board deferred (no hardwood on nct01)
        end
    end
    return s
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
    eq = NC_VOL_EQ[sp]
    length(eq) < 6 && return 0f0
    v = eq[4:6] == "WO2" ? nc_wo2w_vol(eq, d, h) : nc_r5harv_vol(eq, d, h, 6.0f0)
    return max(v[1], 0f0)                       # VOL(1) total cubic
end
