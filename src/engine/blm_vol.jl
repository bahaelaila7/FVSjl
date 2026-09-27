# =============================================================================
# blm_vol.jl — NVEL BLM (Oregon Bureau of Land Management) volume: BLMVOL + BLMTAP + helpers.
#
# Ported from the NVEL volume library shipped in every FVS build (bin/FVSwc_buildDir/blmvol.f, blmtap.f).
# VOLINIT routes every VOLEQ beginning with 'B' here (volinit.f:247); WC's BLM forests 708-711 and PN's
# 708/709/712 carry B00BEHW<fia> / B01BEHW202 / B02BEHW202 (R7_EQN). FVS calls it through NATCRS
# (fvsvol.f) with HTTYPE='F', HT1PRD=0 — so only the total-height (TLH=0) branches are live:
#   • cubic call: MTOPP = TOPD·BARK (inside bark), BFPFLG=0, CUPFLG=1  → VOL(1) total, VOL(4) merch
#   • board call (D ≥ BFMIND): MTOPP = BFTOPD·BARK, BFPFLG=1, CUPFLG=1  → VOL(2) Scribner
# All arithmetic is REAL (Float32), statement order as written.
# =============================================================================

# BLMTHT(4,10) (blmtap.f DATA): taper A = B0 + B1·DBH + B2·HT + B3·DBH·HT, one column per PROFILE.
const _BLM_THT = Float32[
    0.6448   -0.00196  0.0       0.0;
    0.6096   -0.00196  0.0       0.0;
    0.31385   0.0      0.002985 -0.00003386;
    0.4779    0.0      0.0       0.0;
    0.5455   -0.00196  0.0       0.0;
    0.45648   0.00289  0.0       0.0;
    0.6014    0.0      0.0       0.0;
    0.54568   0.0      0.0       0.00000546;
    0.4606    0.0      0.0       0.0;
    0.6200    0.0      0.0       0.0]           # row = PROFILE 1..10, col = coefficient 1..4

"""BLMTAP (blmtap.f:4-85), total-height branch (TLH=0): inside-bark diameter at height `htup` on a stem of
height `httot` (the caller's TTH = HT+1.5), butt-log length `xlen`, form-class diameter `d17`."""
function _blm_tap(dbhob::Float32, httot::Float32, htup::Float32, d17::Float32, xlen::Float32, profile::Int)::Float32
    hbutt = httot - (xlen + 1.5f0)
    hbutt <= 0f0 && return 0f0
    htdib = httot - htup
    t = @view _BLM_THT[profile, :]
    a = t[1] + t[2] * dbhob + t[3] * httot + t[4] * dbhob * httot
    b = 1.0f0 - a
    return d17 * ((htdib / hbutt) / (a * (htdib / hbutt) + b))
end

"BLMTAPEQ (blmvol.f:937): taper-equation number TAPEQU and profile for a BLM VOLEQ (ITAPER=0)."
function _blm_tapeq(voleq::AbstractString)::Tuple{Int,Int}
    pre = voleq[1:3]; sp = voleq[8:10]
    sp == "202" && pre == "B01" && return (1, 1)
    sp == "202" && pre == "B02" && return (2, 2)
    sp == "202" && pre == "B03" && return (3, 10)
    sp == "202" && pre == "B04" && return (4, 10)
    sp == "211" && return (5, 10)
    sp == "202" && pre == "B05" && return (6, 10)
    sp == "122" && pre == "B01" && return (10, 3)
    sp == "122" && pre == "B00" && return (11, 3)
    sp == "116" && return (12, 10)
    sp == "117" && return (13, 4)
    sp == "119" && return (14, 5)
    sp == "108" && return (15, 10)
    sp == "231" && return (20, 10)
    sp == "631" && return (21, 10)
    sp == "351" && return (22, 10)
    sp == "998" && return (23, 10)
    sp == "312" && return (24, 10)
    sp == "361" && return (25, 10)
    sp == "431" && return (26, 10)
    sp == "542" && return (27, 10)
    sp == "747" && return (28, 10)
    sp == "800" && return (29, 10)
    sp == "015" && pre == "B01" && return (30, 6)
    sp == "015" && (pre == "B00" || pre == "B02") && return (31, 10)
    sp == "021" && return (32, 7)
    sp == "017" && return (33, 6)
    sp == "011" && return (34, 7)
    sp == "022" && return (35, 7)
    sp == "093" && return (41, 10)
    sp == "098" && return (42, 10)
    (sp == "260" || sp == "263") && return (48, 8)
    sp == "081" && return (51, 9)
    sp == "042" && return (52, 10)
    sp == "041" && return (53, 10)
    sp == "242" && return (54, 9)
    sp == "073" && return (55, 9)
    return (56, 10)
end

"""DOUBLE_BARK (blmvol.f:838): DBH inside bark for taper equation `tapequ`. The taper numbers with no case
(4, 6, 10, 30) leave DBHIB unset in the Fortran; none is reachable from the WC/PN BLM equation tables."""
function _blm_double_bark(tapequ::Int, dbhob::Float32)::Float32
    tq = tapequ
    if tq in (1, 2, 3, 5, 35)
        return 0.903563f0 * fpow(dbhob, 0.989388f0)
    elseif tq == 11 || tq == 12
        return 0.809427f0 * fpow(dbhob, 1.016866f0)
    elseif tq == 13 || tq == 14
        return 0.859045f0 * dbhob
    elseif tq == 15
        return dbhob - (0.3147f0 + 0.0274f0 * dbhob)
    elseif tq == 20 || tq == 25
        return -0.03425f0 + 0.98155f0 * dbhob
    elseif tq == 21
        return -4.36852f0 + 0.95354f0 * dbhob + 0.18307f0 * 4.5f0
    elseif tq in (22, 23, 24, 26, 27)
        return 0.39534f0 + 0.90182f0 * dbhob     # blmvol.f:880 "0.3953 4" — fixed-form blanks are ignored
    elseif tq == 28 || tq == 29
        return -0.78034f0 + 0.95956f0 * dbhob
    elseif tq == 31 || tq == 33
        return 0.904973f0 * dbhob
    elseif tq == 32 || tq == 34
        return 0.86951f0 * fpow(dbhob, 1.00983f0)
    elseif tq == 41 || tq == 42
        return dbhob - (0.2113f0 + 0.0445f0 * dbhob)
    elseif tq == 48 || tq == 56
        return dbhob / 1.071f0
    elseif tq == 52 || tq == 54
        return dbhob / 1.053f0
    elseif tq == 51 || tq == 53
        return 0.837291f0 * dbhob
    elseif tq == 55
        return dbhob - (0.1231f0 + 0.1306f0 * dbhob)
    end
    return 0f0
end

"BLMTCUB (blmvol.f:648): total cubic by 4-ft Smalian sections from a 1-ft stump disc to TTH."
function _blm_tcub(profile::Int, dbhob::Float32, tth::Float32, d17::Float32, htlog::Float32)::Float32
    tcvol = 0f0
    tth > 0f0 || return tcvol
    htloop = trunc(Int, (tth + 0.5f0 - 1.0f0) / 4.0f0)
    hgt2 = 1.0f0
    d2 = _blm_tap(dbhob, tth, hgt2, d17, htlog, profile)
    r = d2 / 2.0f0
    tcvol = (3.1416f0 * r * r) / 144.0f0
    for _ in 1:htloop
        d2old = d2
        hgt2 = hgt2 + 4.0f0
        d2 = _blm_tap(dbhob, tth, hgt2, d17, htlog, profile)
        tcvol = tcvol + 0.00272708f0 * (d2old * d2old + d2 * d2) * 4.0f0
    end
    if (tth - hgt2) > 0f0
        tcvol = tcvol + 0.00272708f0 * (d2 * d2) * (tth - hgt2)
    end
    return tcvol
end

"BLMMLEN (blmvol.f:710): merchantable length to the TOP diameter by 0.1-ft bisection."
function _blm_mlen(profile::Int, tth::Float32, dbhob::Float32, d17::Float32, stump::Float32, top::Float32,
                   htlog::Float32)::Float32
    top1 = trunc(Int, top * 10.0f0)                   # TOP1 = AINT(TOP*10.0), an INTEGER
    first = 1
    last = trunc(Int, tth + 0.5f0) * 10
    toplop = last
    for _ in 1:toplop
        first == last && break
        half = (first + last + 1) ÷ 2
        hgt2 = Float32(half) / 10.0f0
        d2 = _blm_tap(dbhob, tth, hgt2, d17, htlog, profile)
        d2 = Float32(trunc((d2 + 0.005f0) * 10.0f0))
        if Float32(top1) <= d2
            first = half
        else
            last = half - 1
        end
    end
    lmerch = Float32(first) / 10.0f0 - stump
    lmerch < 0f0 && (lmerch = 0f0)
    return lmerch
end

"""
    blm_vol(voleq, mtopp, httot, dbhob, fclass; bfpflg, cupflg) -> (vol1, vol2, vol4)

BLMVOL (blmvol.f:6) with HTTYPE='F': total cubic VOL(1), Scribner board VOL(2) (BFPFLG=1) and merch cubic
VOL(4) for a BLM VOLEQ. `mtopp` is the inside-bark merch top, `fclass` the FVS form class.
"""
function blm_vol(voleq::AbstractString, mtopp::Float32, httot::Float32, dbhob::Float32, fclass::Int;
                 bfpflg::Bool = false, cupflg::Bool = true)
    vol1 = 0f0; vol2 = 0f0; vol4 = 0f0
    fclass <= 0 && return (vol1, vol2, vol4)                       # ERRFLAG 2
    httot <= 0f0 && return (vol1, vol2, vol4)                      # ERRFLAG 4
    tth = httot + 1.5f0
    tapequ, profile = _blm_tapeq(voleq)
    dbhib = _blm_double_bark(tapequ, dbhob)
    dbhib <= 0.0001f0 && return (vol1, vol2, vol4)                 # ERRFLAG 14
    mtopp <= 0f0 && (mtopp = Float32(round((0.184f0 * dbhob) + 2.24f0, RoundNearestTiesAway)))
    if tth <= 17.8f0
        return (0.00272708f0 * (dbhib * dbhib) * tth, vol2, vol4)
    end
    smd_17 = Float32(trunc(sqrt(dbhib * dbhib - (dbhib * dbhib) * 17.3f0 / tth) + 0.5f0))
    smd_17 < mtopp && return (0.00272708f0 * (dbhib * dbhib) * tth, vol2, vol4)
    htlog = 16.3f0
    if cupflg
        d17 = Float32(round((dbhob * Float32(fclass)) / 100.0f0, RoundNearestTiesAway))
        vol1 = _blm_tcub(profile, dbhob, tth, d17, htlog)
    end
    # log bucking: OPT 23, EVOD 2, 16-ft logs (8-ft min), 1-ft stump, 0.3-ft trim, 8-ft minimum merch length
    evod = 2; maxlen = 16.0f0; minlen = 8.0f0; opt = 23; stump = 1.0f0; trim = 0.3f0; merchl = 8.0f0
    d17 = Float32(round((dbhob * Float32(fclass)) / 100f0, RoundNearestTiesAway))
    lmerch = _blm_mlen(profile, tth, dbhob, d17, stump, mtopp, htlog)
    lmerch < merchl && return (vol1, vol2, vol4)
    (lmerch / (maxlen + trim)) > 20f0 && return (vol1, vol2, vol4)   # ERRFLAG 12
    numseg = _nvb_numlog(opt, evod, lmerch, maxlen, minlen, trim)
    loglen, numseg = _nvb_segmnt(opt, evod, lmerch, maxlen, minlen, trim, numseg)
    # BLMGDIB (blmvol.f:768): small-end diameters along the bucked logs; the top log never below TOP
    logdia = zeros(Float32, 21)
    if numseg > 0
        logdia[1] = dbhib                                          # STUMP ≤ 4.5
        hgt2 = stump
        for i in 1:numseg
            hgt2 = hgt2 + trim + loglen[i]
            logdia[i + 1] = _blm_tap(dbhob, tth, hgt2, d17, htlog, profile)
        end
        logdia[numseg + 1] < mtopp && (logdia[numseg + 1] = mtopp)
    end
    if bfpflg
        for i in 1:numseg
            dib = Float32(round(logdia[i + 1], RoundNearestTiesAway))
            logv = _scrib(dib, loglen[i], 'N')
            vol2 = vol2 + Float32(round(Float32(round(logv, RoundNearestTiesAway)), RoundNearestTiesAway))
        end
        vol2 < 0f0 && (vol2 = 0f0)
    end
    if cupflg
        dibl = Float32(round(dbhib, RoundNearestTiesAway))
        for i in 1:numseg
            dibs = Float32(round(logdia[i + 1], RoundNearestTiesAway))
            logv = 0.00272708f0 * (dibl * dibl + dibs * dibs) * loglen[i]
            vol4 = vol4 + Float32(round(logv * 10.0f0, RoundNearestTiesAway)) / 10.0f0
            dibl = dibs
        end
        vol4 < 0f0 && (vol4 = 0f0)
    end
    return (vol1, vol2, vol4)
end

"""NATCRS (fvsvol.f) for a BLM VOLEQ: the cubic VOLINIT call (MTOPP = TOPD·BARK → VOL(1), VOL(4)) and, for
D ≥ BFMIND, the board call (MTOPP = BFTOPD·BARK → VOL(2)). FORMCL supplies the form class `fc`.
Returns (TCF, MCF, BBFV) before the driver's DBHMIN/BFMIND gates and the broken-top trim."""
function _blm_natcrs(eq::AbstractString, fc::Int, d::Float32, h::Float32, bark::Float32,
                     topd::Float32, bftopd::Float32, bfmind::Float32)
    v1, _, v4 = blm_vol(eq, topd * bark, h, d, fc; bfpflg = false, cupflg = true)
    bf = 0f0
    if d >= bfmind
        _, bf, _ = blm_vol(eq, bftopd * bark, h, d, fc; bfpflg = true, cupflg = true)
    end
    return (max(v1, 0f0), max(v4, 0f0), max(bf, 0f0))
end
