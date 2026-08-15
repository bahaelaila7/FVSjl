# =============================================================================
# volume.jl (southeastalaska) — AK Region-10 volume (ak/sitset.f VOLEQDEF, IREGN=10). Chunk 8.
#
# AK assigns NVEL equations via VOLEQDEF(VAR='AK',IREGN=10): three families (see
# data/southeastalaska/volume_coefficients.jl). akt01 is 100% F32 (A00F32W###), the Region-10
# Flewelling 2-point profile that dispatches through the SAME NVEL PROFILE routine as the INGY/CR
# FW2 path jl already ports (cr_fw2_vol / _fw2_* kernels) — only the coefficient table (SHP_AK) and
# the bark→DBHIB conversion (FDBT_AK) are AK-specific. So this reuses the shared kernels
# (_fw2_shp_core / _fw2_sf_taper / _fw2_sf_yhat / _fw2_tcubic / _fw2_hs / _nvb_numlog / _nvb_segmnt /
# _scrib) with AK coefficients + Region-10 merch rules (mrules.f REGN.EQ.10).
#
# Merch rules (mrules.f REGN.EQ.10): COR='Y' EVOD=2 MAXLEN=16 MINLEN=8 OPT=23 STUMP=1 MTOPP=6 MTOPS=4
# TRIM=0.5 MERCHL=8 MINBFD=1. Merch top (MERLEN F-branch, profile.f:1004: SF_HS with DS=MTOPP) is the
# INSIDE-bark profile diameter = MTOPP directly (the F/Flewelling profile is inside bark — no ·bark).
#
# DVE (woodland → DVEST/R10D2H) and CUR (AD/RA → A32CURW351 R10TAP profile) are ported below.
# akt01 is 100% F32, so cyc0 is a full F32-path validation; DVE/CUR validated on synthetic all-species
# stands. VALIDATED BIT-EXACT vs live FVSak_clean akt01 TREELIST cyc0 (total + merch cubic + board).
# =============================================================================

# AK species (1..23) → internal Flewelling JSP (fwinit.f GEOCODE 'A'); 0 = not an F32 species (DVE/CUR).
# SF AF YC TA WS LS BE SS LP RC WH MH OS AD RA PB AB BA AS CW WI SU OH
const _AK_VOL_JSP = Int[34, 34, 31, 0, 0, 0, 0, 33, 34, 32, 34, 34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# AK species (1..23) → R10D2H equation group = the voleq(8:10) FIA code that VOLEQDEF assigns
# (MEASURED from the live FVSak_clean VOLEQHEAD table: TA/WS/LS/BE/OS=094, PB/AB/AS=375,
# BA/CW/WI/SU/OH=747). 0 = F32 (conifers) or CUR (AD/RA → A32CURW351, follow-on). These are the
# `A00DVEW###` woodland/interior species that dispatch grossvol MDL='DVE' → DVEST → R10D2H.
# SF AF YC  TA WS LS BE  SS LP RC WH MH  OS  AD RA  PB  AB  BA  AS  CW  WI  SU  OH
const _AK_DVE_GRP = Int[0, 0, 0, 94, 94, 94, 94, 0, 0, 0, 0, 0, 94, 0, 0, 375, 375, 747, 375, 747, 747, 747, 747]

# JSP → F-coefficient column (31=YC/F1, 32=RC/F2, 33=spruce/F3, 34=spruce+hemlock share/F3).
@inline function _ak_vol_fcoef(jsp::Int)
    jsp == 31 && return _AK_VOL_F1
    jsp == 32 && return _AK_VOL_F2
    return _AK_VOL_F3          # jsp 33, 34
end

# NOTE on bark: FVS does NOT use the NVEL FDBT_AK bark for the F32 profile — the volume driver passes
# DBT_USER = DBHOB−DBHOB·BRATIO(ISPC) from the AK VARIANT bark (ak/bratio.f, = ak_bratio), so SF_SHP's
# `IF(DBT_USER.LE.0.)FDBT_AK` branch is bypassed. MEASURED bit-exact vs live SF2PT DBH_IB: DBHIB =
# D·ak_bratio(sp,d) for YC/SS/LP/WH/MH/RC (e.g. SS D5.8 → 5.16169; FDBT_AK would give 5.517, wrong).
# So the profile is scaled to — and the merch tops converted with — the AK variant bark, not FDBT_AK.

# AK merch standards (setcubicdflts.f, VARACD 'AK', AKMERCHCAT by KODFOR). akt01 KODFOR=1005 → CAT 3:
#   DBHMIN=9, TOPD=7 (cubic top), SCFMIND=9, SCFTOPD=7 (board top), STMP=1, SCFSTMP=1.
# The cubic/board merch tops are OUTSIDE-bark diameters converted to the inside-bark profile via ·bark
# (bark = DBHIB/DBHOB from FDBT_AK), matching CI's validated FW2 merch (cr_fw2_vol topd·bark). Bucking
# uses the Region-10 NUMLOG/SEGMNT rules (mrules.f REGN.EQ.10).
const _AK_VOL_DBHMIN = 9.0f0   # min DBH for cubic merch (cfvol.f:150 D<DBHMIN ⇒ not merch)
const _AK_VOL_TOPD   = 7.0f0   # cubic merch top diameter (outside bark)
const _AK_VOL_SCFMIND = 9.0f0  # min DBH for board feet
const _AK_VOL_SCFTOPD = 7.0f0  # board merch top diameter (outside bark)
const _AK_VOL_STUMP = 1.0f0
const _AK_VOL_MINLEN = 8.0f0
const _AK_VOL_MAXLEN = 16.0f0
const _AK_VOL_MERCHL = 8.0f0
const _AK_VOL_TRIM = 0.5f0
const _AK_VOL_OPT = 23
const _AK_VOL_EVOD = 2

# AK merch bucking (profile.f MERLEN→NUMLOG/SEGMNT→GETDIB). Returns (VOL4 merch cubic, VOL2 Scribner
# board, or (0,0) if not merch). `mtop` = inside-bark merch top (= TOPD·bark). Faithful to profile.f:
#  • F3 MINLEN reset (line 367-372): N16SEG=INT(LMERCH/(MAXLEN+TRIM)); if odd, MINLEN=2 (else 8).
#  • GETDIB: per-log endpoint DIB is raw dibat; the top log's DIB is forced ≥ MTOPP; LOGDIA(,1)=inch
#    class (INT, +1 if frac>0.501) for cubic; LOGDIA(,2)=raw (INT-floored) for the 32-ft board.
#  • VOL(4) cubic (line 573-587): butt = dib-class at 4.5; Σ ANINT(.00272708·(DIBL²+DIBS²)·LEN·10)/10.
#  • VOL(2) board (line 464-542, REGN 10 32-ft rule): pair 16-ft logs into 32-ft, DIB=INT(raw small-end
#    of the pair's top), LEN=LOGLEN(2k)+LOGLEN(2k-1); odd top log kept as a 16-ft; Σ SCRIB(DIB,LEN)·10.
function _ak_buck(dibat, h::Float32, mtop::Float32)
    hs = _fw2_hs(dibat, mtop, h)
    lmerch = hs - _AK_VOL_STUMP
    lmerch < _AK_VOL_MERCHL && return (0f0, 0f0)
    n16 = trunc(Int, lmerch / (_AK_VOL_MAXLEN + _AK_VOL_TRIM))   # F3 reset (pre-SEGMNT LMERCH)
    minlen = isodd(n16) ? 2.0f0 : _AK_VOL_MINLEN
    numseg = _nvb_numlog(_AK_VOL_OPT, _AK_VOL_EVOD, lmerch, _AK_VOL_MAXLEN, minlen, _AK_VOL_TRIM)
    numseg == 0 && return (0f0, 0f0)
    loglen, numseg = _nvb_segmnt(_AK_VOL_OPT, _AK_VOL_EVOD, lmerch, _AK_VOL_MAXLEN, minlen, _AK_VOL_TRIM, numseg)
    numseg == 0 && return (0f0, 0f0)
    # per-log endpoint raw DIBs (GETDIB TAPERMODEL); butt = dib at 4.5; top log forced ≥ MTOPP
    rawdib = zeros(Float32, numseg + 1)
    rawdib[1] = dibat(4.5f0)
    ht2 = _AK_VOL_STUMP
    @inbounds for i in 1:numseg
        ht2 += _AK_VOL_TRIM + loglen[i]
        rawdib[i + 1] = dibat(ht2)
    end
    rawdib[numseg + 1] < mtop && (rawdib[numseg + 1] = mtop)
    # VOL(4) cubic: inch-class DIBs, per-log Smalian, 0.1-rounded
    dibl = _fw2_dclass(rawdib[1]); vol4 = 0f0
    @inbounds for i in 1:numseg
        dibs = _fw2_dclass(rawdib[i + 1])
        logv = 0.00272708f0 * (dibl * dibl + dibs * dibs) * loglen[i]
        vol4 += floor(logv * 10f0 + 0.5f0) / 10f0
        dibl = dibs
    end
    # VOL(2) board: 32-ft-log Scribner (pairs of 16-ft logs); DIB=INT(raw), LEN=pair sum
    vol2 = 0f0
    @inbounds for i in 2:2:numseg
        dib = Float32(floor(rawdib[i + 1]))                 # INT(LOGDIA(I+1,2))
        len = loglen[i] + loglen[i - 1]
        vol2 += _scrib(dib, len, 'Y') * 10f0
    end
    if isodd(numseg)                                        # top 16-ft log
        dib = Float32(floor(rawdib[numseg + 1]))
        vol2 += _scrib(dib, loglen[numseg], 'Y') * 10f0
    end
    return (vol4, vol2)
end

# Per-tree AK F32 volume (2-point Flewelling, inside-bark profile). Returns (total_cuft, merch_cuft, bdft).
function _ak_f32_vol(sp::Int, jsp::Int, d::Float32, h::Float32)
    (d < 1f0 || h <= 4.5f0) && return (0f0, 0f0, 0f0)   # need H>4.5 for the breast-height scaling
    fcoef = _ak_vol_fcoef(jsp)
    rflw, rhfw = _fw2_shp_core(fcoef, d, h, false)      # SHP_AK ≡ SHP_OT kernel; no lodgepole branch
    tapcoe = _fw2_sf_taper(rhfw, rflw)
    bark = ak_bratio(sp, d)                             # AK variant bark (ak/bratio.f) = DBT_USER path
    dbhib = d * bark                                    # inside-bark DBH (bit-exact vs live SF2PT DBH_IB)
    yhat_bh = _fw2_sf_yhat(4.5f0 / h, tapcoe, rhfw, rflw, 1.0f0)
    yhat_bh == 0f0 && return (0f0, 0f0, 0f0)
    f = dbhib / yhat_bh
    dibat = ht -> _fw2_sf_yhat(ht / h, tapcoe, rhfw, rflw, f)   # inside-bark section diameter
    stump_dib = h <= 15f0 ? _fw2_fwsmall(jsp, h, dibat(1.0f0), dbhib) : -1f0
    tcf = Float32(round(_fw2_tcubic(dibat, h; stump_dib = stump_dib) * 10.0f0)) / 10.0f0
    # Merch cubic (VOL4) + Scribner board (VOL2) share ONE bucking to the MTOPP=TOPD·bark top
    # (profile.f: board reuses the cubic logs). FVS-side gates: DBHMIN=9 (cubic), SCFMIND=9 (board).
    mcf, bf = _ak_buck(dibat, h, _AK_VOL_TOPD * bark)
    d < _AK_VOL_DBHMIN && (mcf = 0f0)
    d < _AK_VOL_SCFMIND && (bf = 0f0)
    return (tcf, mcf, bf)
end

# R10D2H (ak/r10d2h.f) — Region-10 direct D²H volume estimators (Larsen & Winterberger, PNW-RN-478/495;
# board Scribner from PNW-RN-495 eq.6). Coastal AK path EQN='00' (VOLEQ 'A00DVEW###', the only variant
# used per the VOLEQHEAD crosswalk — A01 interior-NW/A02 statewide are not assigned to any AK forest).
# DBHOB is the tree's OUTSIDE-bark dbh — no bark conversion (unlike the F32 profile). Returns
# (total_cuft VOL(1), merch_cuft VOL(4), bdft VOL(2)). dvest.f rounds ONLY VOL(2)=ANINT (2025/05/07:
# VOL(1)/VOL(4) left unrounded for biomass); r10d2h clamps every VOL<0 → 0.
#
# The FVS driver (fvsvol.f NATCRS→VOLINIT) applies the AK merch STANDARDS after r10d2h returns:
# TCF=VOL(1) unconditionally, but **MCF=VOL(4) only if DBHOB≥DBHMIN and BBFV=VOL(2) only if DBHOB≥BFMIND**
# (fvsvol.f:512/519). For KODFOR 1005 (AKMERCHCAT 3, setcubicdflts.f) DBHMIN=BFMIND=SCFMIND=9 — MEASURED
# from the live NATCRS dump (DBHMIN=9.0 BFMIND=9.0). So a D<9 tree contributes total-only. (r10d2h's own
# D>4 / D>6 inner gates still hold; the 9" driver gate dominates.)
const _AK_DVE_DBHMIN = 9f0   # DBHMIN (merch cubic) = BFMIND (board) for AK cat-3 (KODFOR 1005)
function _ak_dve_vol(grp::Int, d::Float32, h::Float32)
    (d <= 1f0 || h <= 0f0) && return (0f0, 0f0, 0f0)   # r10d2h ERRFLAG 3 (DBH≤1) / 4 (HT≤0)
    d2h = d * d * h
    v1 = 0f0; v4 = 0f0; v2 = 0f0
    if grp == 94                       # A00DVEW094: TA/WS/LS/BE/OS
        v1 = 0.65559f0 + 0.00191f0 * d2h
        (d < 6f0) && (v1 -= 0.65559f0 * (1f0 - (d / 6f0)^3))   # small-tree (DBH<6) correction
        d > 4f0 && (v4 = -0.21849f0 + 0.00189f0 * d2h)
        d > 6f0 && (v2 = 0.000136f0 * d2h^1.40338f0)
    elseif grp == 375                  # A00DVEW375: PB/AB/AS
        v1 = 0.64456f0 + 0.00206f0 * d2h
        (d < 6f0) && (v1 -= 0.64456f0 * (1f0 - (d / 6f0)^3))
        d > 4f0 && (v4 = -0.7126f0 + 0.00211f0 * d2h)
        d > 6f0 && (v2 = 0.000081f0 * d2h^1.48459f0)
    elseif grp == 747                  # A00DVEW747: BA/CW/WI/SU/OH
        v1 = 0.9864f0 + 0.00181f0 * d2h
        (d < 6f0) && (v1 -= 0.9864f0 * (1f0 - (d / 6f0)^3))
        d > 4f0 && (v4 = -1.39764f0 + 0.00188f0 * d2h)
        d > 6f0 && (v2 = -28.0674f0 + 0.00937f0 * d2h)
    end
    v1 = max(v1, 0f0)                                          # r10d2h VOL(1)<0 → 0
    v4 = d >= _AK_DVE_DBHMIN ? max(v4, 0f0) : 0f0              # driver DBHMIN gate on merch cubic
    v2 = d >= _AK_DVE_DBHMIN ? Float32(round(max(v2, 0f0))) : 0f0  # driver BFMIND gate + VOL(2)=ANINT
    return (v1, v4, v2)
end

# =============================================================================
# CUR (AD/RA → A32CURW351) — Region-10 red-alder taper (R10TAP DVREDA) via the shared PROFILE
# machinery (grossvol MDL='CUR' → PROFILE). Both hardwoods use the SAME NVEL equation A32CURW351
# (MEASURED from the live FVSak_clean VOLEQHEAD table). MDL='CUR' with VOLEQ(1:3)='A32' dispatches
# to PROFILE (not R10VOL) only for the large-tree branch; the FVS driver (volinit.f:466) routes
# DBHOB<9 (REGN 10) to R10VOL instead. So the CUR path is TWO routines:
#   • D≥9  → PROFILE: total cubic = TCUBIC 4-ft Smalian of the R10TAP profile; merch length from
#            R10HTS (Newton solve DVREDA(RH)=TOP²/D²); A32 32-ft-log bucking (VOL4 cubic, VOL2 board).
#   • D<9  → R10VOL small-tree: FSTGRO (D≤3.5 or H<18) / SECGRO cubic estimators, total cubic only.
# DBHOB is OUTSIDE-bark; R10TAP red-alder taper carries NO bark (BK=0 for RA), so profile diameters
# are used directly. Merch top MTOPP = BFTOPD·bark = 7·ak_bratio(sp,D) (fvsvol.f:382). Merch/board
# gated to D≥DBHMIN=9 by the driver. VALIDATED per-tree bit-exact vs live FVSak_clean on all-AD/all-RA
# synthetic stands (akAD/akRA .trl cyc0): merch cubic + Scribner board bit-exact for ALL trees; total
# cubic bit-exact except the damage-97 broken-top tree (0.1 cuft — the same synthetic broken-top height
# artifact documented for the DVE chunk) and Float32 0.1-cuft summation straddles on tall dead trees.
# =============================================================================

# DVREDA (r10tap.f/r10hts.f red-alder statement function): the inside-bark taper RATIO at relative
# height RH (before √·D). Fortran evaluates it in single precision (all operands REAL*4, result stored
# into a REAL variable); computing the powers in Float64 then narrowing to Float32 reproduces it to
# well within the 0.1-cuft output precision.
@inline function _ak_cur_dvreda(rh::Float32, rh32::Float32, rh40::Float32, h::Float32, d::Float32)::Float32
    r = Float64(rh); dd = Float64(d); hh = Float64(h)
    r15 = r^1.5; r3 = r^3.0
    v = 0.91274 * r15 -
        1.9758 * (r15 - r3) * (dd * 1e-2) +
        8.2375 * (r15 - r3) * hh * 1e-3 -
        4.964 * (r15 - Float64(rh32)^32.0) * (hh * dd) * 1e-5 +
        3.773 * (r15 - Float64(rh32)^32.0) * (hh^0.5) * 1e-3 -
        7.417 * (r15 - Float64(rh40)^40.0) * (hh^2.0) * 1e-6
    return Float32(v)
end

# R10TAP red-alder inside-bark section diameter at height `htup` (RH clamp per r10tap.f:169-196).
@inline function _ak_cur_dib(d::Float32, h::Float32, htup::Float32)::Float32
    rh = (h - htup) / (h - 4.5f0)
    rh <= 0f0 && return 0f0
    rh32 = rh; rh40 = rh
    if rh < 0.078f0
        rh40 = 0.15f0; rh32 = 0.078f0
    elseif rh < 0.15f0
        rh40 = 0.15f0
    end
    d2 = _ak_cur_dvreda(rh, rh32, rh40, h, d)
    d2 < 0f0 && (d2 = 0f0)
    return sqrt(d2) * d
end

# RH clamp for the R10HTS Newton iteration (RA rules). `rh13` selects the RH40=0.13 quirk that the
# initial-RH clamp uses (r10hts.f:1863) — RXL and inner-loop clamps use 0.15 (rh13=false).
@inline function _ak_cur_clamp(rh::Float32, rh13::Bool)
    rh <= 0f0 && return (0f0, 0f0, 0f0)
    rh < 0.078f0 && return (rh, 0.078f0, 0.15f0)
    rh < 0.15f0 && return (rh, rh, rh13 ? 0.13f0 : 0.15f0)
    return (rh, rh, rh)
end

# R10HTS (profile.f:1642, RA branch, total height given): Newton solve for the relative height RH where
# the taper ratio DVREDA(RH) equals TOP²/D², then merch length = (H − RH·(H−4.5)) − STUMP.
function _ak_cur_lmerch(d::Float32, h::Float32, top::Float32, stump::Float32)::Float32
    limd = top
    xll = 1f0 - (2f0/3f0) * (limd / d)
    hh = h * xll
    rh = (h - hh) / (h - 4.5f0)
    _, rh32, rh40 = _ak_cur_clamp(rh, true)
    dst = (limd * limd) / (d * d)
    dslo = dst - 0.0001f0; dshi = dst + 0.0001f0
    ds = _ak_cur_dvreda(rh, rh32, rh40, h, d); ds < 0f0 && (ds = 0f0)
    rxl = 0.9f0 * rh
    _, rh32x, rh40x = _ak_cur_clamp(rxl, false)
    dxl = _ak_cur_dvreda(rxl, rh32x, rh40x, h, d); dxl < 0f0 && (dxl = 0f0)
    taper = (ds - dxl) / (0.1f0 * rh)
    @inbounds for _ in 1:10
        (ds > dslo && ds < dshi) && break
        rh = taper != 0f0 ? rh + (dst - ds) / taper : rh
        rh, rh32, rh40 = _ak_cur_clamp(rh, false)
        ds = _ak_cur_dvreda(rh, rh32, rh40, h, d); ds < 0f0 && (ds = 0f0)
    end
    lm = h - rh * (h - 4.5f0) - stump
    return lm < 0f0 ? 0f0 : lm
end

# R10VOL small-tree cubic (r10vol.f FSTGRO/SECGRO), the CUR/DEM D<9 (REGN 10) path — total cubic only.
function _ak_cur_smalltree(d::Float32, h::Float32)::Float32
    if d <= 3.5f0 || h < 18f0                 # FSTGRO
        h <= 4.5f0 && return 0f0
        if h <= 18f0
            t1 = (h - 0.9f0) / (h - 4.5f0); t1 = t1 * t1
            t2 = t1 * (h - 0.9f0) / (h - 4.5f0)
            form = 0.406098f0 * t1 - 0.0762998f0 * d * t2 + 0.00262615f0 * d * h * t2
        else
            form = 0.480961f0 + 42.46542f0/(h*h) - 10.99643f0*d/(h*h) - 0.107809f0*d/h - 0.00409083f0*d
        end
        vn = 0.005454154f0 * form * d * d * h
        return vn < 0f0 ? 0f0 : vn
    else                                       # SECGRO (D>3.5 and H≥18, up to D<9/H≤40)
        return exp(-5.577f0 + 1.9067f0 * log(d) + 0.9416f0 * log(h))
    end
end

# A32 bucking (profile.f, REGN 10 32-ft-log rule) for the CUR profile — same structure as the F32
# `_ak_buck` but with the R10HTS-derived merch length (`lmerch`) instead of the SF_HS bisection.
function _ak_cur_buck(dibat, lmerch::Float32, mtop::Float32, stump::Float32)
    lmerch < _AK_VOL_MERCHL && return (0f0, 0f0)
    n16 = trunc(Int, lmerch / (_AK_VOL_MAXLEN + _AK_VOL_TRIM))    # F3/A32 MINLEN reset
    minlen = isodd(n16) ? 2.0f0 : _AK_VOL_MINLEN
    numseg = _nvb_numlog(_AK_VOL_OPT, _AK_VOL_EVOD, lmerch, _AK_VOL_MAXLEN, minlen, _AK_VOL_TRIM)
    numseg == 0 && return (0f0, 0f0)
    loglen, numseg = _nvb_segmnt(_AK_VOL_OPT, _AK_VOL_EVOD, lmerch, _AK_VOL_MAXLEN, minlen, _AK_VOL_TRIM, numseg)
    numseg == 0 && return (0f0, 0f0)
    rawdib = zeros(Float32, numseg + 1)
    rawdib[1] = dibat(4.5f0)
    ht2 = stump
    @inbounds for i in 1:numseg
        ht2 += _AK_VOL_TRIM + loglen[i]
        rawdib[i + 1] = dibat(ht2)
    end
    rawdib[numseg + 1] < mtop && (rawdib[numseg + 1] = mtop)
    dibl = _fw2_dclass(rawdib[1]); vol4 = 0f0
    @inbounds for i in 1:numseg
        dibs = _fw2_dclass(rawdib[i + 1])
        logv = 0.00272708f0 * (dibl * dibl + dibs * dibs) * loglen[i]
        vol4 += floor(logv * 10f0 + 0.5f0) / 10f0
        dibl = dibs
    end
    vol2 = 0f0
    @inbounds for i in 2:2:numseg
        dib = Float32(floor(rawdib[i + 1]))
        len = loglen[i] + loglen[i - 1]
        vol2 += _scrib(dib, len, 'Y') * 10f0
    end
    if isodd(numseg)
        dib = Float32(floor(rawdib[numseg + 1]))
        vol2 += _scrib(dib, loglen[numseg], 'Y') * 10f0
    end
    return (vol4, vol2)
end

# Per-tree AK CUR volume (AD/RA). Returns (total_cuft, merch_cuft, bdft).
function _ak_cur_vol(sp::Int, d::Float32, h::Float32)
    (d < 1f0 || h < 5f0) && return (0f0, 0f0, 0f0)          # PROFILE DBHOB<1 / HTTOT<5 guard
    if d < _AK_VOL_DBHMIN                                   # driver routes D<9 (REGN 10) to R10VOL
        return (Float32(round(_ak_cur_smalltree(d, h) * 10f0)) / 10f0, 0f0, 0f0)
    end
    bark = ak_bratio(sp, d)
    mtop = _AK_VOL_TOPD * bark                              # MTOPP = BFTOPD·bark = 7·ak_bratio
    stump = d > 36f0 ? d / 36f0 : _AK_VOL_STUMP             # DEM/CUR stump fix (profile.f:1145)
    dibat = ht -> _ak_cur_dib(d, h, ht)
    tcf = Float32(round(_fw2_tcubic(dibat, h) * 10f0)) / 10f0
    lmerch = _ak_cur_lmerch(d, h, mtop, stump)
    mcf, bf = _ak_cur_buck(dibat, lmerch, mtop, stump)
    return (tcf, mcf, bf)
end

function compute_volumes_ak!(s::StandState)
    t = s.trees
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        jsp = (1 <= sp <= 23) ? _AK_VOL_JSP[sp] : 0
        grp = (1 <= sp <= 23) ? _AK_DVE_GRP[sp] : 0
        iscur = (sp == 14 || sp == 15)                        # AD/RA → A32CURW351
        if d < 1f0 || (jsp == 0 && grp == 0 && !iscur)
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        tcf, mcf, bf = jsp != 0 ? _ak_f32_vol(sp, jsp, d, Float32(h)) :
                       iscur    ? _ak_cur_vol(sp, d, Float32(h)) :
                                  _ak_dve_vol(grp, d, Float32(h))
        t.cuft_vol[i] = max(tcf, 0f0)
        t.merch_cuft_vol[i] = max(mcf, 0f0)
        t.saw_cuft_vol[i] = 0f0
        t.bdft_vol[i] = max(bf, 0f0)
    end
    return s
end
