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
function _ak_buck(dibat, h::Float32, mtop::Float32; hs_solver = nothing)
    hs = hs_solver === nothing ? _fw2_hs(dibat, mtop, h) : hs_solver(mtop)
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

# f_alaska.f SHP_AK: the SHP_OT model form (same U1..U9 → R1..R5/A3/RHI1/RHI2/RHC tail, no lodgepole branch), but
# DMEDIAN and DFORM are declared REAL*4 there (REAL*8 in SHP_OT): DMEDIAN = F10·(HT−4.5)**(F11+F12·HT) is rounded to
# single, DFORM = DBHIB/DMEDIAN − 1.0D0 is a SINGLE division widened for the double subtract, then rounded to single.
# Through the shared double kernel RFLW(2)/RFLW(4) came out 1 ULP off (RC 10.5"/50': 3F164684/3F473A3F vs live
# 3F164685/3F473A40) ⇒ TAPCOE/F off ⇒ SF_HS merch-top height (FVS_TreeList Ht2TD) 24.78582 vs 24.785955.
function _ak_shp_core(f, d::Float32, h::Float32)
    @inline ff(i) = f[i - 9]
    H = Float64(h); lnH = Float64(flog(h))
    dmedian = Float32(ff(10) * dpow(Float64(h - 4.5f0), ff(11) + ff(12) * H))
    dform = Float64(Float32(Float64(d / dmedian) - 1.0))
    D = Float64(d)
    u7 = ff(13) + ff(14) * lnH + ff(15) * dform
    u9t = clamp(ff(18) + ff(19) * lnH + ff(20) * dform, -7.0, 7.0)
    u9 = ff(16) * dexp(u9t) / (1.0 + dexp(u9t))
    u8 = ff(21) + ff(22) * H + ff(23) * lnH + ff(24) * dform
    u1 = ff(25) + ff(26) * lnH + ff(27) * dform + ff(28) * dform * lnH
    u2 = ff(29) + ff(30) * dform + ff(31) * lnH + ff(32) * dform * lnH + ff(33) * D
    u3 = ff(34) + ff(35) * dform + ff(36) * lnH + ff(37) * lnH * dform
    u4 = ff(38) + ff(39) * dform + ff(40) * lnH + ff(41) * D
    u5 = ff(42) + ff(43) * lnH
    u6 = ff(45) + ff(46) * dform + ff(47) * lnH
    u1 = clamp(u1, -7.0, 7.0); u2 = clamp(u2, -7.0, 7.0)
    u3 = clamp(u3, -7.0, 7.0); u4 = clamp(u4, -7.0, 7.0)
    u5 = clamp(u5, -7.0, 7.1)
    u6 = u6 < Float64(1.005f0) ? 1.005 : (u6 > 10.0 ? 10.0 : u6)
    u7 = clamp(u7, -7.0, 7.0)
    u8 > Float64(0.99f0) && (u8 = Float64(0.99f0))
    u9 = u9 > 0.3 ? 0.3 : (u9 < 0.0 ? 0.0 : u9)
    r1 = dexp(u1) / (1.0 + dexp(u1)); r2 = dexp(u2) / (1.0 + dexp(u2))
    r3 = dexp(u3) / (1.0 + dexp(u3)); r4 = dexp(u4) / (1.0 + dexp(u4))
    r5 = u5 <= 7.0 ? 0.5 + 0.5 * dexp(u5) / (1.0 + dexp(u5)) : 1.0
    rhi1 = Float32(dexp(u7) / (1.0 + dexp(u7))); rhi1 > 0.5f0 && (rhi1 = 0.5f0)
    rhlongi = Float32(u9)
    rhi2 = rhi1 + rhlongi
    rhc = Float32(u8)
    rhc < rhi2 + 0.01f0 && (rhc = min(rhi2 + 0.01f0, (rhi2 + 1.0f0) / 2.0f0))
    return (Float32(r1), Float32(r2), Float32(r3), Float32(r4), Float32(r5), Float32(u6)), (rhi1, rhi2, rhc, rhlongi)
end

# Per-tree AK F32 volume (2-point Flewelling, inside-bark profile). Returns (total_cuft, merch_cuft, bdft).
# `m` = AK merch standards (ak_merch): cubic top TOPD / board top BFTOPD (outside bark, ·bark inside the
# kernel), DBHMIN / BFMIND gates. fvsvol.f makes the cubic VOLINIT call with MTOPP=TOPD·BARK and, for
# D>=BFMIND, a separate board call with MTOPP=BFTOPD·BARK.
function _ak_f32_vol(sp::Int, jsp::Int, d::Float32, h::Float32, m = _AK_MERCH_CAT3, bk::Float32 = -1f0;
                     ht2td = nothing)
    (d < 1f0 || h <= 4.5f0) && return (0f0, 0f0, 0f0)   # need H>4.5 for the breast-height scaling
    fcoef = _ak_vol_fcoef(jsp)
    rflw, rhfw = _ak_shp_core(fcoef, d, h)              # f_alaska.f SHP_AK (SHP_OT form, REAL*4 DMEDIAN/DFORM)
    tapcoe = _fw2_sf_taper(rhfw, rflw)
    bark = bk > 0f0 ? bk : ak_bratio(sp, d)             # AK variant bark (ak/bratio.f) = DBT_USER path (vols.f BARK)
    # sf_shp.f:30 DBHIB = DBHOB − DBTBH with DBTBH = DBT_USER = D·(1−BARK) (fvsvol.f DBTBH) — not D·BARK: the two
    # round differently (FIA 38178692010760 2069 RC rec 8: TCuFt 28.6 live vs 28.7 through D·BARK).
    dbhib = d - d * (1f0 - bark)
    # sf_2pt.f:63-66 F = DBH_IB / SF_YHAT(4.5/TOTALH) through sf_yhat.f's own precision (as cr_fw2_vol)
    yhat_bh = _fw2_sf_yhat_f(4.5f0 / h, tapcoe, rhfw, rflw, 1.0f0, h, false)[1]
    yhat_bh == 0f0 && return (0f0, 0f0, 0f0)
    f = dbhib / yhat_bh
    # inside-bark section diameter = TAPERMODEL → SF_DS(INEEDSL=0) → SF_YHAT at sf_yhat.f's OWN precision (REAL X /
    # REAL sub-terms, REAL*8 Y) — `_fw2_sf_yhat_f`, not the all-Float64 kernel: FIA 38178692010760 2069 RC 13.9"×65'
    # TCVOL live 28.649929 (41E5330E) = _f, the Float64 kernel gave 28.650032 ⇒ NINT(·10) 28.7 vs 28.6.
    dibat = ht -> _fw2_sf_yhat_f(ht / h, tapcoe, rhfw, rflw, f, h, false)[1]
    stump_dib = h <= 15f0 ? _fw2_fwsmall(jsp, h, dibat(1.0f0), dbhib) : -1f0
    tcf = _nint(_fw2_tcubic(dibat, h; stump_dib = stump_dib) * 10.0f0) * 1f-1   # profile.f:293 NINT(TCVOL*10.0)*1E-1
    # Merch cubic (VOL4) + Scribner board (VOL2) share ONE bucking to the MTOPP=TOPD·bark top
    # (profile.f: board reuses the cubic logs). FVS-side gates: DBHMIN=9 (cubic), SCFMIND=9 (board).
    # MERLEN (profile.f:1004, VOLEQ(4:4)='F'): merch-top height from SF_HS — the Newton/bisection solver itself,
    # not a diameter-tolerance bisection (the two straddle an even-foot SEGMNT boundary, e.g. YC 22.5"×59').
    hs_solver = top -> _fw2_sf_hs(tapcoe, rhfw, rflw, f, h, top)
    # HT1PRD (fvsvol.f:338/485 → HT2TD): profile.f:335-342 MERLEN F-branch LMERCH = HS−STUMP (≥0, profile.f:1052),
    # HT1PRD = LMERCH+STUMP — of the cubic call (MTOPP=TOPD·BARK) and of the board call (MTOPP=BFTOPD·BARK).
    if ht2td !== nothing
        _ht1prd(top) = (lm = hs_solver(top) - _AK_VOL_STUMP; lm < 0f0 && (lm = 0f0); lm + _AK_VOL_STUMP)
        ht2td[1] = _ht1prd(m.topd * bark); ht2td[2] = _ht1prd(m.bftopd * bark)
    end
    mcf, bf = _ak_buck(dibat, h, m.topd * bark; hs_solver = hs_solver)
    m.bftopd != m.topd && (bf = _ak_buck(dibat, h, m.bftopd * bark; hs_solver = hs_solver)[2])
    d < m.dbhmin && (mcf = 0f0)
    d < m.bfmind && (bf = 0f0)
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
# (fvsvol.f:512/519), with the per-forest ak/sitset.f standards (`ak_merch`: DBHMIN 9/9/6/5 for MERCHCAT
# 3/4/1/2, BFMIND 9). r10d2h's own D>4 / D>6 inner gates still hold underneath.
function _ak_dve_vol(grp::Int, d::Float32, h::Float32, m = _AK_MERCH_CAT3)
    (d <= 1f0 || h <= 0f0) && return (0f0, 0f0, 0f0)   # r10d2h ERRFLAG 3 (DBH≤1) / 4 (HT≤0)
    d2h = d * d * h
    v1 = 0f0; v4 = 0f0; v2 = 0f0
    if grp == 94                       # A00DVEW094: TA/WS/LS/BE/OS
        v1 = 0.65559f0 + 0.00191f0 * d2h
        (d < 6f0) && (v1 -= 0.65559f0 * (1f0 - (d / 6f0)^3))   # small-tree (DBH<6) correction
        d > 4f0 && (v4 = -0.21849f0 + 0.00189f0 * d2h)
        d > 6f0 && (v2 = 0.000136f0 * fpow(d2h, 1.40338f0))   # D2H**1.40338 = powf
    elseif grp == 375                  # A00DVEW375: PB/AB/AS
        v1 = 0.64456f0 + 0.00206f0 * d2h
        (d < 6f0) && (v1 -= 0.64456f0 * (1f0 - (d / 6f0)^3))
        d > 4f0 && (v4 = -0.7126f0 + 0.00211f0 * d2h)
        d > 6f0 && (v2 = 0.000081f0 * fpow(d2h, 1.48459f0))
    elseif grp == 747                  # A00DVEW747: BA/CW/WI/SU/OH
        v1 = 0.9864f0 + 0.00181f0 * d2h
        (d < 6f0) && (v1 -= 0.9864f0 * (1f0 - (d / 6f0)^3))
        d > 4f0 && (v4 = -1.39764f0 + 0.00188f0 * d2h)
        d > 6f0 && (v2 = -28.0674f0 + 0.00937f0 * d2h)
    end
    v1 = max(v1, 0f0)                                          # r10d2h VOL(1)<0 → 0
    v4 = d >= m.dbhmin ? max(v4, 0f0) : 0f0                     # driver DBHMIN gate on merch cubic
    v2 = d >= m.bfmind ? round(max(v2, 0f0), RoundNearestTiesAway) : 0f0   # BFMIND gate + dvest.f VOL(2)=ANINT
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
# height RH (before √·D), evaluated in REAL*4 exactly as written (`_r10tap_dvreda`, with the `**` REAL
# powers as powf).
@inline _ak_cur_dvreda(rh::Float32, rh32::Float32, rh40::Float32, h::Float32, d::Float32)::Float32 =
    _r10tap_dvreda(rh, rh32, rh40, h, d)

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
    return fpow(d2, 0.5f0) * d                              # (D2)**.5*D — powf, not sqrt
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

# R10VOL small-tree branch (r10vol.f:~60-80): taken by volinit.f for every A01/A02 (DEM) tree and for CUR trees
# with (HTTYPE F and HTTOT<=40) or DBHOB<9 (REGN 10). FSTGRO (D<=3.5 or H<18) sets VOL(1) only; SECGRO sets
# VOL(1)=VOL(4)=ANINT(CUBVOL*10.0)/10.0 (the driver's DBHMIN gate then decides MCF). REAL arithmetic as written.
function _ak_r10vol_small(d::Float32, h::Float32)
    local cub::Float32
    if d <= 3.5f0 || h < 18f0                 # FSTGRO
        if h <= 4.5f0
            cub = 0f0
        else
            if h <= 18f0
                t1 = ((h - 0.9f0) * (h - 0.9f0)) / ((h - 4.5f0) * (h - 4.5f0))
                t2 = t1 * (h - 0.9f0) / (h - 4.5f0)
                form = 0.406098f0 * t1 - 0.0762998f0 * d * t2 + 0.00262615f0 * d * h * t2
            else
                form = 0.480961f0 + 42.46542f0 / (h * h) - 10.99643f0 * d / (h * h) - 0.107809f0 * d / h - 0.00409083f0 * d
            end
            cub = 0.005454154f0 * form * d * d * h
            cub < 0f0 && (cub = 0f0)
        end
        v1 = round(cub * 10.0f0, RoundNearestTiesAway) / 10.0f0        # VOL(1) = ANINT(CUBVOL*10.0)/10.0 (ties away)
        return (max(v1, 0f0), 0f0)
    else                                       # SECGRO
        cub = fexp(-5.577f0 + 1.9067f0 * flog(d) + 0.9416f0 * flog(h))
        v = round(cub * 10.0f0, RoundNearestTiesAway) / 10.0f0
        return (max(v, 0f0), max(v, 0f0))       # VOL(1) = VOL(4)
    end
end
@inline _ak_r10_is_small(d::Float32, h::Float32) = h <= 40f0 || d < 9.0f0   # volinit.f R10VOL routing (HTTYPE 'F')

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
function _ak_cur_vol(sp::Int, d::Float32, h::Float32, m = _AK_MERCH_CAT3, bk::Float32 = -1f0; ht2td = nothing)
    (d < 1f0 || h < 5f0) && return (0f0, 0f0, 0f0)          # PROFILE DBHOB<1 / HTTOT<5 guard
    if _ak_r10_is_small(d, h)                              # volinit.f routes H<=40 or D<9 (REGN 10) to R10VOL
        v1, v4 = _ak_r10vol_small(d, h)
        return (v1, d >= m.dbhmin ? v4 : 0f0, 0f0)
    end
    bark = bk > 0f0 ? bk : ak_bratio(sp, d)
    stump = d > 36f0 ? d / 36f0 : _AK_VOL_STUMP             # DEM/CUR stump fix (profile.f:1145)
    dibat = ht -> _ak_cur_dib(d, h, ht)
    tcf = _nint(_fw2_tcubic(dibat, h) * 10f0) * 1f-1                   # profile.f:293 NINT(TCVOL*10.0)*1E-1
    mtop = m.topd * bark                                    # cubic call MTOPP = TOPD·BARK
    if ht2td !== nothing                                    # profile.f:337-342 HT1PRD = LMERCH+STUMP (cubic, board)
        _ht1prd(top) = (lm = _ak_cur_lmerch(d, h, top, stump); lm < 0f0 && (lm = 0f0); lm + stump)
        ht2td[1] = _ht1prd(mtop); ht2td[2] = _ht1prd(m.bftopd * bark)
    end
    mcf, bf = _ak_cur_buck(dibat, _ak_cur_lmerch(d, h, mtop, stump), mtop, stump)
    if m.bftopd != m.topd                                   # board call MTOPP = BFTOPD·BARK
        btop = m.bftopd * bark
        bf = _ak_cur_buck(dibat, _ak_cur_lmerch(d, h, btop, stump), btop, stump)[2]
    end
    d < m.dbhmin && (mcf = 0f0)
    d < m.bfmind && (bf = 0f0)
    return (tcf, mcf, bf)
end

# =============================================================================
# DEM (A01DEMW000 — Demars/Schmidt Chugach equation; FORST '04' SF/AF/LP/RC/MH). volinit.f sends every A01/A02
# equation to R10VOL: the small-tree FSTGRO/SECGRO branch (`_ak_r10vol_small`) when HTTOT<=40 or DBHOB<9, else
# R10VOLO — R10TAPO (Demars taper, height-type 'F') locates the merch top and the 16.3-ft log breaks, VOL(4) is the
# Smalian sum over ANINT log diameters, VOL(2) the Chugach 16-ft Scribner table (BFCHUG, 'A01'), and VOL(1) is
# R10TCO's 4-ft Smalian integral of R10TAP (unrounded). All REAL*4 as written; `**` with a REAL exponent is powf.
# =============================================================================

# r10tap.f: CHARACTER*2 ISP is a saved local (-fno-automatic) assigned only for TAPEQU(8:10) '042'/'242'/'098'/
# '351'. The DEM equation ('…000') never sets it, so a DEM tree's R10TAP runs whichever branch the LAST R10TAP call
# left: blank at process start (⇒ BB·DD2MI Demars taper) or 'RA' after any large AD/RA (A32CURW351 PROFILE) tree
# (⇒ the red-alder DVREDA taper with BK=0). Held in Control.ak_r10tap_ra (`ra` below); run_keyfile carries it.

# DD2MI (r10tapo.f:409 FUNCTION / r10tap.f statement function — same expression) and BB.
@inline function _r10_dd2mi(rh::Float32, r32::Float32, d::Float32, h::Float32)::Float32
    p15 = fpow(rh, 1.5f0)
    a = ((-(0.0052554f0 * h)) + 0.000034947f0 * (h * h)) + (0.104477f0 * h) / d
    b = (7.76807f0 / (d * d) - 0.0000094852f0 * (h * h)) - (0.011351f0 * h) / d
    return (p15 + a * (p15 - fpowi(rh, 3))) + b * (p15 - fpowi(r32, 32))
end
@inline _r10_bb(d::Float32, h::Float32) = (0.8467f0 + 0.0009144f0 * d) + 0.0003568f0 * h

# r10tap.f DVREDA statement function (REAL*4 operands; 10**(-2.) etc. are folded constants).
@inline function _r10tap_dvreda(rh::Float32, rh32::Float32, rh40::Float32, h::Float32, d::Float32)::Float32
    p15 = fpow(rh, 1.5f0); p3 = fpow(rh, 3f0)
    r32 = fpow(rh32, 32f0); r40 = fpow(rh40, 40f0)
    t1 = 0.91274f0 * p15
    t2 = (1.9758f0 * (p15 - p3)) * (d * 0.01f0)
    t3 = ((8.2375f0 * (p15 - p3)) * h) * 0.001f0
    t4 = ((4.964f0 * (p15 - r32)) * (h * d)) * 1f-5
    t5 = ((3.773f0 * (p15 - r32)) * fpow(h, 0.5f0)) * 0.001f0
    t6 = ((7.417f0 * (p15 - r40)) * (h * h)) * 1f-6
    return ((((t1 - t2) + t3) - t4) + t5) - t6
end

# R10TAP (r10tap.f) for a DEM equation: inside-bark diameter at HTUP; branch chosen by the saved ISP.
function _r10tap_dem(d::Float32, h::Float32, htup::Float32, ra::Bool)::Float32
    bk = ra ? 0f0 : _r10_bb(d, h)
    rh = (h - htup) / (h - 4.5f0)
    rh <= 0f0 && return 0f0
    if rh < 0.078f0
        rh40 = 0.15f0; rh32 = 0.078f0
    elseif rh < 0.15f0
        rh40 = 0.15f0; rh32 = rh
    else
        rh40 = rh; rh32 = rh
    end
    d2 = ra ? _r10tap_dvreda(rh, rh32, rh40, h, d) : _r10_dd2mi(rh, rh32, d, h)
    d2 < 0f0 && (d2 = 0f0)
    return ra ? fpow(d2, 0.5f0) * d : fpow(d2 * bk, 0.5f0) * d
end

# R10TCO (r10volo.f:312): 1-ft stump cylinder + 4-ft Smalian sections + Smalian tip, from R10TAP.
function _r10tco(d::Float32, h::Float32, ra::Bool)::Float32
    htloop = trunc(Int, (h - 1f0) / 4f0)
    ht2 = 1f0
    d2 = _r10tap_dem(d, h, ht2, ra)
    r = d2 / 2f0
    tcvol = ((3.1416f0 * r) * r) / 144f0
    for _ in 1:htloop
        d2old = d2
        ht2 += 4f0
        d2 = _r10tap_dem(d, h, ht2, ra)
        tcvol += (0.00272708f0 * (d2old * d2old + d2 * d2)) * 4f0
    end
    (h - ht2) > 0f0 && (tcvol += (0.00272708f0 * (d2 * d2)) * (h - ht2))
    return tcvol
end

@inline _r10_rhclamp(rh::Float32) = rh <= 0f0 ? (0f0, 0f0) : rh < 0.078f0 ? (rh, 0.078f0) : (rh, rh)

# R10TAPO (r10tapo.f, LTYPE 'F' branch, TOPDIA(1)=MTOPP) + R10VOLO's LOGDIA/LOGSPAN fill. Returns
# (nlog, logdia (ANINT, 1..nlog+1, top clamped ≥ MTOPP), logspan (1..nlog), CL).
function _r10tapo_logs(d::Float32, h::Float32, mtop::Float32)
    st = d > 36f0 ? d / 36f0 : 1f0
    bk = _r10_bb(d, h)
    dst = (mtop * mtop) / (bk * (d * d))
    dslo = dst - 0.00001f0; dshi = dst + 0.00001f0
    xll = 1f0 - (2f0 / 3f0) * (mtop / d)
    hh = h * xll
    rh, r32 = _r10_rhclamp((h - hh) / (h - 4.5f0))
    ds = _r10_dd2mi(rh, r32, d, h)
    rxl, r32x = _r10_rhclamp(0.9f0 * rh)
    dxl = _r10_dd2mi(rxl, r32x, d, h)
    taper = (ds - dxl) / (0.1f0 * rh)
    for _ in 1:10
        (ds > dslo && ds < dshi) && break
        rh, r32 = _r10_rhclamp(rh + (dst - ds) / taper)
        ds = _r10_dd2mi(rh, r32, d, h)
    end
    hh1 = h - rh * (h - 4.5f0)
    nn = 21
    for n in 1:21
        if 16.3f0 * Float32(n) + st >= hh1
            nn = n; break
        end
    end
    kbd = nn - 1
    cl = (hh1 - 16.3f0 * Float32(kbd)) - st
    sh = Vector{Float32}(undef, nn + 1)
    sh[1] = 4.5f0
    for kb in 1:kbd
        sh[kb + 1] = st + Float32(kb) * 16.3f0
    end
    sh[kbd + 2] = hh1
    logdia = Vector{Float32}(undef, nn + 1)
    for jj in 1:(nn + 1)
        rhj, r32j = _r10_rhclamp((h - sh[jj]) / (h - 4.5f0))
        dib = _r10_dd2mi(rhj, r32j, d, h)
        dib < 0f0 && (dib = 0f0)
        logdia[jj] = round(fpow(dib * bk, 0.5f0) * d, RoundNearestTiesAway)
    end
    sh[1] = 1f0                                            # HTSEG(1)=1.0; D>36 stump correction (08/19/2021)
    d > 36f0 && (sh[1] = sh[2] - 16.3f0)
    logspan = Float32[(sh[i + 1] - sh[i]) - 0.3f0 for i in 1:nn]
    logdia[nn + 1] < mtop && (logdia[nn + 1] = mtop)
    return nn, logdia, logspan, cl
end

# r10volo.f TAB16: (bd.ft.)/10 Scribner for 16-ft logs, 5" to 70" small-end diameter.
const _R10_TAB16 = Float32[2, 2, 3, 3, 4, 6, 7, 8, 10, 11, 14, 16, 18, 21, 24, 28, 30, 33, 38, 40, 46, 50, 55, 58,
    61, 66, 71, 74, 78, 80, 88, 92, 103, 107, 112, 120, 127, 134, 140, 148, 152, 159, 166, 173, 180, 187, 195, 202,
    210, 218, 227, 235, 244, 252, 261, 270, 280, 289, 299, 309, 319, 329, 339, 350, 361, 372]

# R10VOLO board, EQNUM 'A01' ⇒ VOL(2)=BFCHUG: per 16.3-ft log, 16-ft Scribner interpolated on the small-end
# diameter, the top log replaced by 20·CL/16.3, each rounded to whole board feet.
function _r10_bfchug(nlog::Int, logdia::Vector{Float32}, cl::Float32)::Float32
    bfchug = 0f0
    for i in 1:nlog
        di = logdia[i + 1]
        di1 = trunc(Int, di)
        (5 <= di1 <= 69) || error("R10VOLO TAB16 index out of range (DI=$di)")
        w2 = 10f0 * (di - Float32(di1)); w1 = 10f0 - w2
        bf1 = w1 * _R10_TAB16[di1 - 4] + w2 * _R10_TAB16[di1 - 3]
        i == nlog && (bf1 = (20f0 * cl) / 16.3f0)
        bf1 = Float32(trunc(Int, bf1 + 0.5f0))               # JBF = BF1+0.5; later ITEMP=BF1+0.5 (already whole)
        bfchug += round(bf1 * 10f0, RoundNearestTiesAway) / 10f0
    end
    return bfchug
end

# Per-tree DEM volume. Returns (total_cuft, merch_cuft, bdft) after the driver's DBHMIN/BFMIND gates.
function _ak_dem_vol(sp::Int, d::Float32, h::Float32, m, ra::Bool, bk::Float32 = -1f0)
    (d < 1f0 || h <= 0f0) && return (0f0, 0f0, 0f0)
    if _ak_r10_is_small(d, h)                              # r10vol.f small-tree branch
        v1, v4 = _ak_r10vol_small(d, h)
        return (v1, d >= m.dbhmin ? v4 : 0f0, 0f0)
    end
    bark = bk > 0f0 ? bk : ak_bratio(sp, d)
    nlog, logdia, logspan, _ = _r10tapo_logs(d, h, m.topd * bark)   # cubic call: MTOPP = TOPD·BARK
    vol4 = 0f0; dibl = logdia[1]
    for i in 1:nlog
        dibs = logdia[i + 1]
        logv = (0.00272708f0 * (dibl * dibl + dibs * dibs)) * logspan[i]
        vol4 += round(logv * 10f0, RoundNearestTiesAway) / 10f0
        dibl = dibs
    end
    vol4 = max(vol4, 0f0)
    vol1 = _r10tco(d, h, ra)                               # D>=9 and H>40 here
    vol2 = 0f0
    if d >= m.bfmind                                       # board call: MTOPP = BFTOPD·BARK
        nb, ldb, _, clb = _r10tapo_logs(d, h, m.bftopd * bark)
        vol2 = _r10_bfchug(nb, ldb, clb)
    end
    return (max(vol1, 0f0), d >= m.dbhmin ? vol4 : 0f0, vol2)
end

# ak/sitset.f merch standards: MERCHCDS maps KODFOR → MERCHCAT (a code not in the table ⇒ 3), and DBHMIN/BFMIND/
# TOPD/BFTOPD/SCFMIND/SCFTOPD are filled by category where still ≤ 0 (grinit.f zeroes them; stumps 1).
#   cat 1 (713):              DBHMIN 6, TOPD 4, BFTOPD 6        cat 2 (720, 7400-7408): DBHMIN 5, TOPD 4, BFTOPD 6
#   cat 3 (1005, 703, 81xx):  DBHMIN 9, TOPD 7, BFTOPD 7        cat 4 (1004):           DBHMIN 9, TOPD 6, BFTOPD 6
# BFMIND = SCFMIND = 9 in every category; SCFTOPD = BFTOPD. jl hard-coded category 3 for every forest.
const _AK_MERCHCDS = Dict(713 => 1, 720 => 2, 7400 => 2, 7401 => 2, 7402 => 2, 7403 => 2, 7404 => 2, 7405 => 2,
                          7406 => 2, 7407 => 2, 7408 => 2, 703 => 3, 1005 => 3, 8134 => 3, 8135 => 3, 8112 => 3,
                          1004 => 4)
@inline function _ak_merch_cat(cat::Int)
    dbhmin = cat == 1 ? 6f0 : cat == 2 ? 5f0 : 9f0
    topd   = (cat == 1 || cat == 2) ? 4f0 : cat == 4 ? 6f0 : 7f0
    bftopd = cat == 3 ? 7f0 : 6f0
    return (dbhmin = dbhmin, topd = topd, bfmind = 9f0, bftopd = bftopd, scfmind = 9f0, scftopd = bftopd, stump = 1f0)
end
const _AK_MERCH_CAT3 = _ak_merch_cat(3)
ak_merch_cat(kodfor::Integer) = get(_AK_MERCHCDS, Int(kodfor), 3)

# The per-stand standards as the volume kernels see them: the Control merch arrays (init_merch_standards! fills
# them from ak/sitset.f; VOLUME/BFVOLUME keywords override them there).
function ak_merch(s::StandState, sp::Int)
    s.control.merch_init || init_merch_standards!(s)
    c = s.control
    return (dbhmin = c.sp_dbh_min[sp], topd = c.sp_top_diam[sp], bfmind = c.sp_bf_dbhmin[sp],
            bftopd = c.sp_bf_topd[sp], scfmind = c.sp_scf_dbhmin[sp], scftopd = c.sp_scf_topd[sp],
            stump = c.sp_stump_ht[sp], bfstump = c.sp_bf_stump[sp])
end
# AK merch-cubic min DBH (DBHMIN) for species `sp` (FMCROWE's DBHMIN(SPIYV)).
ak_merch_dbhmin(s::StandState, sp::Int)::Float32 = ak_merch(s, sp).dbhmin

# ak/sitset.f: KODFOR 713, 720 and the 74xx reservations (and 1004 itself, INTFOR 4) run VOLEQDEF with FORST='04'
# (Chugach), everything else '05' (Tongass). Live FVSak_g16 equation tables: on FORST '04' SF/AF/LP/RC/MH take
# A01DEMW000 (R10VOL Demars) and YC A00DVEW094; on '05' those are A00F32W260/242/042. The rest are identical.
@inline _ak_forst04(kodfor::Integer) = kodfor == 1004 || kodfor == 713 || kodfor == 720 || kodfor ÷ 100 == 74
@inline _ak_is_dem04(sp::Int) = sp == 1 || sp == 2 || sp == 9 || sp == 10 || sp == 12

# Per-tree AK NVEL volume (TCF, MCF, BF) at DBH `d`, height `h` — the NATCRS result before any broken-top trim
# (also FFE's FMSVL2 volume for FMCROWE's DBHMIN tree).
function ak_tree_vol(s::StandState, sp::Int, d::Float32, h::Float32, bk::Float32 = -1f0; ht2td = nothing)
    (1 <= sp <= 23) || return (0f0, 0f0, 0f0)
    m = ak_merch(s, sp)
    if _ak_forst04(Int(s.plot.user_forest_code))
        _ak_is_dem04(sp) && return _ak_dem_vol(sp, d, h, m, s.control.ak_r10tap_ra, bk)
        sp == 3 && return _ak_dve_vol(94, d, h, m)          # YC → A00DVEW094 on FORST '04'
    end
    jsp = _AK_VOL_JSP[sp]; grp = _AK_DVE_GRP[sp]
    iscur = (sp == 14 || sp == 15)                            # AD/RA → A32CURW351
    (d < 1f0 || (jsp == 0 && grp == 0 && !iscur)) && return (0f0, 0f0, 0f0)
    # a large CUR tree goes through PROFILE → R10TAP('…351'), which leaves the saved ISP = 'RA'
    iscur && d >= 1f0 && h >= 5f0 && !_ak_r10_is_small(d, h) && (s.control.ak_r10tap_ra = true)
    return jsp != 0 ? _ak_f32_vol(sp, jsp, d, h, m, bk; ht2td = ht2td) :
           iscur    ? _ak_cur_vol(sp, d, h, m, bk; ht2td = ht2td) :
                      _ak_dve_vol(grp, d, h, m)
end

function compute_volumes_ak!(s::StandState)
    t = s.trees
    htb = zeros(Float32, 2)
    @inbounds for i in 1:(t.n + t.ndead)
        # vols.f:88 zeroes HT2TD for every record, then vols.f:110 `IF(P.LE.0.0) GO TO 200` skips a record MORTS zeroed
        # this cycle — its CFV/MCFV/BFV keep the previous VOLS values (FIA 38178692010760 2049: 15 WH seedlings killed
        # at D≈1.4" report TCuFt 0.0, their pre-growth volume; jl re-volumed them at 0.1).
        t.merch_top_cf[i] = 0f0; t.merch_top_bf[i] = 0f0
        (i <= t.n && t.tpa[i] <= 0f0) && continue
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        # ak/vols.f:145-146 (== ca/vols.f): a top-killed tree (H≥4.5, ITRUNC>0) is volumed at its NORMAL height
        # NORMHT, then vols.f:193/391 trim it back to the break with CFTOPK/BFTOPK (NATCRS sets CTKFLG=BTKFLG=.TRUE.
        # for every NVEL family — F32, DVE and CUR alike; VMAX=BFMAX=TCF). jl volumed the broken stem at its
        # standing height with no trim.
        tkill = h >= 4.5f0 && t.trunc[i] > 0
        hv = tkill ? Float32(t.norm_ht[i]) / 100f0 : Float32(h)
        # vols.f:150-151: BARK=BRATIO(ISPC,D,H) is taken on the START-of-cycle DBH, THEN D=D+DG/BARK; NVEL gets
        # DBTBH=D·(1−BARK) (fvsvol.f:131/159/337) and the merch tops TOPD·BARK with that start bark. jl volumes the
        # grown DBH, so use the BRATIO(D_start) the DBH update stashed (vol_bark); cycle 0 (LSTART) has none ⇒ current.
        bk = (i <= t.n && t.vol_bark[i] > 0f0) ? t.vol_bark[i] : -1f0
        # HT2TD (vols.f:88 zeroes both for every record; fvsvol.f:337-339/484-487 then set them for the live IT>0
        # records): CF = the cubic call's HT1PRD if D≥DBHMIN, BF = the board call's if D≥BFMIND. Only the PROFILE
        # families return one (F32, large CUR); R10D2H (DVE), R10VOLO (DEM) and the small-tree R10VOL leave it 0.
        htb[1] = 0f0; htb[2] = 0f0
        tcf, mcf, bf = ak_tree_vol(s, sp, d, hv, bk; ht2td = htb)
        if 1 <= sp <= 23                                     # live (IT=I) and the cycle-0 dead records alike
            m = ak_merch(s, sp)
            t.merch_top_cf[i] = d >= m.dbhmin ? htb[1] : 0f0
            t.merch_top_bf[i] = d >= m.bfmind ? htb[2] : 0f0
        end
        if tkill && tcf > 0f0 && 1 <= sp <= 23
            bark = bk > 0f0 ? bk : ak_bratio(sp, d)          # vols.f:150 BARK=BRATIO(ISPC,D,H)
            vmax = tcf
            m = ak_merch(s, sp)                              # CFTOPK: STMP/TOPD; BFTOPK: BFSTMP/BFTOPD (per forest)
            tcf, mcf = cr_cftopk(tcf, mcf, d, hv, vmax, bark, Int(t.trunc[i]), m.stump, m.topd)
            bf = cr_bftopk(bf, d, hv, vmax, bark, Int(t.trunc[i]), m.bfstump, m.bftopd)
        end
        t.cuft_vol[i] = max(tcf, 0f0)
        t.merch_cuft_vol[i] = max(mcf, 0f0)
        t.saw_cuft_vol[i] = 0f0
        t.bdft_vol[i] = max(bf, 0f0)
    end
    return s
end
