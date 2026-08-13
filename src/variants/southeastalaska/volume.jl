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
# STUB for DVE (woodland TA/WS/LS/BE/OS/PB/AB/BA/AS/CW/WI/SU/OH → DVEST) and CUR (AD/RA → R10 CUR
# profile): those families are not yet ported (0 volume). akt01 contains none of them, so cyc0 is a
# full F32-path validation; DVE/CUR are a follow-on chunk.
# VALIDATED BIT-EXACT vs live FVSak_clean akt01 TREELIST cyc0 (total + merch cubic + board).
# =============================================================================

# AK species (1..23) → internal Flewelling JSP (fwinit.f GEOCODE 'A'); 0 = not an F32 species (DVE/CUR).
# SF AF YC TA WS LS BE SS LP RC WH MH OS AD RA PB AB BA AS CW WI SU OH
const _AK_VOL_JSP = Int[34, 34, 31, 0, 0, 0, 0, 33, 34, 32, 34, 34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

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

function compute_volumes_ak!(s::StandState)
    t = s.trees
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        jsp = (1 <= sp <= 23) ? _AK_VOL_JSP[sp] : 0
        if jsp == 0 || d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        tcf, mcf, bf = _ak_f32_vol(sp, jsp, d, Float32(h))
        t.cuft_vol[i] = max(tcf, 0f0)
        t.merch_cuft_vol[i] = max(mcf, 0f0)
        t.saw_cuft_vol[i] = 0f0
        t.bdft_vol[i] = max(bf, 0f0)
    end
    return s
end
