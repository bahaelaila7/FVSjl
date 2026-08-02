# =============================================================================
# r4vol.jl — Region-4 Matney volume (NVEL R4VOL / R4MATTAPER). Used by TT (MATW eq).
#
# Ported from bin/FVStt_buildDir/r4vol.f. VOLEQ "400MATW<code>" (MDL=MAT) → the Matney
# taper: STUMPD/BUTTCF/CF0/B from CFCOEF(20,7) (indexed by II from VOLEQ code+geocode),
# then a Smalian log-by-log cubic integration over 16.5-ft logs. Gross merch cubic (CFGRS).
# =============================================================================

const TT_R4C1 = Float32[0.5563, 0.4496, 0.2359, 0.8435, 1.1152, 0.5405, 1.2552, 0.4302, 0.4727, 2.0461, 0.6294, 1.2467, 0.5334, 0.5301, 0.8315, 0.5028, 0.7572, 0.631, 0.485, 0.6076]
const TT_R4C2 = Float32[-0.0636, -0.3956, -0.7237, -0.28, -0.3678, -0.5908, -0.1727, -0.5229, -0.0508, -0.508, -0.0939, -0.1615, -0.0071, 0.2327, 0.0967, -0.0942, -0.4527, -0.1026, -0.3636, -0.173]
const TT_R4C3 = Float32[0.99, 1.2756, 1.6343, 1.0264, 0.9828, 1.3374, 0.9435, 1.3824, 1.0615, 0.9113, 0.9886, 0.9013, 1.0127, 0.8748, 0.8293, 1.06, 1.1502, 1.023, 1.2324, 1.0585]
const TT_R4C4 = Float32[-0.2, -0.041, -0.041, -0.041, -0.191, -0.191, -0.131, -0.131, -0.159, -0.041, -0.488, -0.2, -0.2, -0.365, -0.365, -0.365, -0.143, -0.143, 0.0, -0.153]
const TT_R4C5 = Float32[0.964, 0.884, 0.884, 0.884, 0.943, 0.943, 0.886, 0.886, 0.832, 0.884, 0.894, 0.964, 0.964, 0.887, 0.887, 0.887, 0.933, 0.933, 0.887, 0.883]
const TT_R4C6 = Float32[-0.201, -0.041, -0.041, -0.041, -0.192, -0.192, -0.159, -0.159, -0.055, -0.041, -0.445, -0.201, -0.201, -0.367, -0.367, -0.367, -0.144, -0.144, 0.0, -0.143]
const TT_R4C7 = Float32[0.968, 0.888, 0.888, 0.888, 0.947, 0.947, 0.891, 0.891, 0.837, 0.888, 0.897, 0.968, 0.968, 0.891, 0.891, 0.891, 0.937, 0.937, 0.893, 0.886]

@inline _fint(x) = trunc(Int, x)      # Fortran INT (truncate toward zero)

# r4vol.f R4MATTAPER — II index from VOLEQ(8:10)=code, VOLEQ(1:3)=geocode.
function _tt_r4_ii(code::AbstractString, geo::AbstractString)
    code == "746" && return 1
    code == "202" && return geo == "400" ? 2 : geo == "405" ? 3 : geo == "401" ? 4 : 0
    code == "019" && return geo == "400" ? 5 : geo == "405" ? 6 : 0
    code == "015" && return geo == "400" ? 7 : geo == "401" ? 8 : 0
    code == "081" && return 9
    code == "073" && return 10
    code == "122" && return geo == "403" ? 11 : geo == "401" ? 14 : geo == "402" ? 15 : geo == "400" ? 16 : 0
    code == "108" && return geo == "400" ? 12 : geo == "401" ? 13 : 0
    code == "093" && return geo == "400" ? 17 : geo == "407" ? 18 : 0
    code == "020" && return 19
    code == "117" && return 20
    return 0
end

"Region-4 Matney gross merch cubic-foot volume (r4vol.f, M=3 path). Returns CFGRS."
r4vol_cubic(voleq::AbstractString, dbhob::Float32, httot::Float32, mtopp::Float32, ht1prd::Float32)::Float32 =
    r4vol_volumes(voleq, dbhob, httot, mtopp, ht1prd)[2]

"Region-4 Matney volumes (r4vol.f): returns (CF0 total-stem cubic [VOL(1)], CFGRS gross merch cubic [VOL(4)])."
function r4vol_volumes(voleq::AbstractString, dbhob::Float32, httot::Float32, mtopp::Float32, ht1prd::Float32)::Tuple{Float32,Float32}
    (dbhob < 1f0 || httot <= 4.5f0) && return (0f0, 0f0)
    tht = httot - 1f0
    tht <= 5f0 && return (dbhob * dbhob * httot * 0.00272708f0, 0f0)   # small tree (THT≤5): total=Smalian, merch=0 (r4vol.f:169-172 only sets VOL(1))
    ii = _tt_r4_ii(strip(voleq)[8:10], strip(voleq)[1:3])
    ii == 0 && return (0f0, 0f0)
    ht67 = TT_R4C1[ii] * dbhob^TT_R4C2[ii] * tht^TT_R4C3[ii]
    buttcf = TT_R4C5[ii] * dbhob + TT_R4C4[ii]
    stumpd = sqrt(buttcf * buttcf * tht / (tht - 4f0))
    d67 = TT_R4C7[ii] * dbhob * (2f0 / 3f0) + TT_R4C6[ii]
    cf0 = 0.002727f0 * (ht67 * stumpd * stumpd + d67 * d67 * tht)
    f = cf0 / (0.005454f0 * stumpd * stumpd * tht)
    b = (1f0 - f) / (2f0 * f)
    trm = 0.5f0
    topdia = mtopp <= 0f0 ? 1f0 : mtopp
    topdia >= buttcf && return (cf0, 0f0)                           # M=3: TOPDIA≥BUTTCF → no merch vol
    dratio = topdia / stumpd; dratio <= 0f0 && (dratio = 0.0001f0)
    merlen = tht - tht * dratio^(1f0 / b)
    if ht1prd > 0f0 && ht1prd < merlen
        merlen = ht1prd
    end
    merlen < 2.5f0 && return (cf0, 0f0)
    totlgs = merlen / 16.5f0
    toplen = (totlgs - _fint(totlgs)) * 16.5f0
    if toplen < 2.5f0
        toplen = totlgs >= 1f0 ? 16.5f0 : 0f0
        totlgs = Float32(_fint(totlgs))
    else
        set = false
        for L in 4:2:16
            if toplen < Float32(L) + trm
                toplen = Float32(L - 2) + trm; set = true; break
            end
        end
        set || (toplen = 16.5f0)
    end
    totlgs <= 0f0 && return (cf0, 0f0)
    totlgs = toplen < 16.5f0 ? Float32(_fint(totlgs) + 1) : Float32(_fint(totlgs))
    merlen = (totlgs - 1f0) * 16.5f0 + toplen
    numlgs = _fint(totlgs)
    @inline dsm_at(hcut) = Float32(_fint(stumpd * ((tht - hcut) / tht)^b + 0.499f0))
    @inline logcf(dlg, dsm, len) = Float32(_fint(0.002727f0 * (dlg * dlg + dsm * dsm) * len * 10f0 + 0.499f0)) / 10f0
    cfgrs = 0f0
    dsm_prev = 0f0
    if numlgs > 1
        dlg1 = Float32(_fint(buttcf + 0.499f0))
        dsm1 = dsm_at(16.5f0)
        cfgrs += logcf(dlg1, dsm1, 16f0)
        dsm_prev = dsm1
        for num in 2:(numlgs - 1)
            dlg = dsm_prev
            dsm = dsm_at(16.5f0 * num)
            cfgrs += logcf(dlg, dsm, 16f0)
            dsm_prev = dsm
        end
    end
    # top log
    dlg_top = numlgs == 1 ? Float32(_fint(buttcf + 0.499f0)) : dsm_prev
    len_top = toplen - trm
    dsm_top = Float32(_fint(stumpd * ((tht - merlen) / tht)^b + 0.499f0))
    cfgrs += logcf(dlg_top, dsm_top, len_top)
    return (cf0, cfgrs)
end

# =============================================================================
# Region-4 Matney BOARD-foot (Scribner) — r4vol.f M=1 path. VOL(2)=BFGRS.
# =============================================================================
const TT_SCRIBC = Float32[0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0, 6.0, 7.0, 8.0, 8.0, 8.0, 9.0, 9.0, 10.0, 10.0, 11.0, 11.0, 12.0, 14.0, 14.0, 15.0, 16.0, 17.0, 18.0, 18.0, 19.0, 20.0, 21.0, 21.0, 23.0, 24.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 39.0, 40.0, 41.0, 42.0, 44.0, 45.0, 47.0, 48.0, 49.0, 51.0, 52.0, 54.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 6.0, 6.0, 7.0, 8.0, 8.0, 10.0, 10.0, 11.0, 13.0, 14.0, 15.0, 15.0, 16.0, 18.0, 18.0, 20.0, 20.0, 22.0, 23.0, 26.0, 27.0, 28.0, 30.0, 32.0, 33.0, 35.0, 37.0, 38.0, 39.0, 41.0, 43.0, 45.0, 47.0, 48.0, 50.0, 52.0, 54.0, 56.0, 59.0, 61.0, 63.0, 65.0, 67.0, 70.0, 72.0, 74.0, 77.0, 79.0, 82.0, 85.0, 87.0, 90.0, 93.0, 96.0, 98.0, 101.0, 104.0, 107.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 11.0, 12.0, 13.0, 14.0, 15.0, 17.0, 19.0, 21.0, 22.0, 23.0, 25.0, 27.0, 28.0, 29.0, 30.0, 33.0, 35.0, 39.0, 40.0, 42.0, 45.0, 48.0, 50.0, 52.0, 56.0, 57.0, 59.0, 62.0, 65.0, 67.0, 70.0, 73.0, 76.0, 79.0, 82.0, 85.0, 88.0, 91.0, 95.0, 98.0, 101.0, 105.0, 108.0, 112.0, 116.0, 119.0, 123.0, 127.0, 131.0, 135.0, 139.0, 144.0, 148.0, 152.0, 157.0, 161.0, 1.0, 1.0, 1.0, 2.0, 3.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 11.0, 12.0, 14.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 29.0, 31.0, 33.0, 36.0, 37.0, 39.0, 40.0, 44.0, 46.0, 51.0, 54.0, 56.0, 60.0, 64.0, 67.0, 70.0, 74.0, 76.0, 79.0, 83.0, 86.0, 90.0, 94.0, 97.0, 101.0, 105.0, 109.0, 113.0, 118.0, 122.0, 126.0, 131.0, 135.0, 140.0, 145.0, 149.0, 154.0, 159.0, 164.0, 170.0, 175.0, 180.0, 186.0, 192.0, 197.0, 203.0, 209.0, 215.0, 1.0, 1.0, 2.0, 3.0, 3.0, 4.0, 5.0, 6.0, 7.0, 9.0, 10.0, 12.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 29.0, 31.0, 34.0, 36.0, 38.0, 41.0, 44.0, 46.0, 49.0, 50.0, 55.0, 58.0, 64.0, 67.0, 70.0, 75.0, 79.0, 84.0, 87.0, 93.0, 95.0, 99.0, 104.0, 108.0, 112.0, 117.0, 122.0, 127.0, 132.0, 137.0, 142.0, 147.0, 152.0, 158.0, 163.0, 169.0, 175.0, 181.0, 187.0, 193.0, 199.0, 206.0, 212.0, 219.0, 226.0, 232.0, 240.0, 247.0, 254.0, 261.0, 269.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 6.0, 7.0, 9.0, 11.0, 12.0, 14.0, 16.0, 18.0, 21.0, 23.0, 25.0, 28.0, 30.0, 34.0, 37.0, 41.0, 44.0, 46.0, 49.0, 53.0, 55.0, 59.0, 60.0, 66.0, 69.0, 77.0, 80.0, 84.0, 90.0, 95.0, 101.0, 105.0, 111.0, 114.0, 119.0, 124.0, 130.0, 135.0, 140.0, 146.0, 152.0, 158.0, 164.0, 170.0, 176.0, 183.0, 189.0, 196.0, 203.0, 210.0, 217.0, 224.0, 232.0, 239.0, 247.0, 254.0, 262.0, 271.0, 279.0, 287.0, 296.0, 305.0, 314.0, 323.0, 1.0, 2.0, 2.0, 3.0, 4.0, 5.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0, 19.0, 21.0, 24.0, 27.0, 29.0, 33.0, 35.0, 40.0, 44.0, 48.0, 51.0, 53.0, 57.0, 62.0, 64.0, 69.0, 70.0, 77.0, 81.0, 90.0, 93.0, 98.0, 105.0, 111.0, 117.0, 122.0, 129.0, 133.0, 139.0, 145.0, 151.0, 157.0, 164.0, 170.0, 177.0, 184.0, 191.0, 198.0, 206.0, 213.0, 221.0, 229.0, 237.0, 245.0, 253.0, 261.0, 270.0, 279.0, 288.0, 297.0, 306.0, 316.0, 325.0, 335.0, 345.0, 356.0, 366.0, 377.0, 2.0, 3.0, 3.0, 4.0, 6.0, 7.0, 8.0, 10.0, 11.0, 14.0, 16.0, 18.0, 21.0, 24.0, 28.0, 30.0, 33.0, 38.0, 40.0, 46.0, 50.0, 55.0, 58.0, 61.0, 66.0, 71.0, 74.0, 78.0, 80.0, 88.0, 92.0, 103.0, 107.0, 112.0, 120.0, 127.0, 134.0, 140.0, 148.0, 152.0, 159.0, 166.0, 173.0, 180.0, 187.0, 195.0, 202.0, 210.0, 218.0, 227.0, 235.0, 244.0, 252.0, 261.0, 270.0, 280.0, 289.0, 299.0, 309.0, 319.0, 329.0, 339.0, 350.0, 361.0, 372.0, 383.0, 395.0, 406.0, 418.0, 430.0, 2.0, 3.0, 3.0, 4.0, 6.0, 8.0, 9.0, 11.0, 13.0, 16.0, 18.0, 21.0, 24.0, 27.0, 31.0, 34.0, 38.0, 42.0, 45.0, 52.0, 56.0, 62.0, 65.0, 68.0, 74.0, 80.0, 83.0, 88.0, 90.0, 98.0, 104.0, 116.0, 120.0, 126.0, 135.0, 143.0, 151.0, 157.0, 166.0, 171.0, 178.0, 186.0, 194.0, 202.0, 211.0, 219.0, 228.0, 237.0, 246.0, 255.0, 264.0, 274.0, 284.0, 294.0, 304.0, 315.0, 325.0, 336.0, 348.0, 358.0, 370.0, 381.0, 393.0, 406.0, 419.0, 430.0, 444.0, 457.0, 471.0, 484.0, 2.0, 3.0, 3.0, 4.0, 7.0, 8.0, 10.0, 12.0, 14.0, 18.0, 20.0, 23.0, 27.0, 30.0, 35.0, 38.0, 42.0, 47.0, 50.0, 57.0, 62.0, 68.0, 73.0, 76.0, 82.0, 89.0, 92.0, 98.0, 100.0, 109.0, 115.0, 129.0, 133.0, 140.0, 150.0, 159.0, 168.0, 174.0, 185.0, 190.0, 198.0, 207.0, 216.0, 225.0, 234.0, 243.0, 253.0, 263.0, 273.0, 283.0, 294.0, 304.0, 315.0, 327.0, 338.0, 350.0, 362.0, 373.0, 387.0, 398.0, 412.0, 423.0, 437.0, 452.0, 465.0, 478.0, 493.0, 508.0, 523.0, 538.0]  # (70,10) col-major: [(len-1)*70+(dib-1)]
@inline _scribc(dib::Int, len::Int) = (1 <= dib <= 70 && 1 <= len <= 10) ? TT_SCRIBC[(len-1)*70 + dib] : 0f0

"Region-4 Matney gross Scribner board-foot volume (r4vol.f M=1 path). Returns BFGRS."
function r4vol_board(voleq::AbstractString, dbhob::Float32, httot::Float32, mtopp::Float32, ht1prd::Float32)::Float32
    (dbhob < 1f0 || httot <= 4.5f0) && return 0f0
    tht = httot - 1f0
    tht <= 5f0 && return 0f0
    ii = _tt_r4_ii(strip(voleq)[8:10], strip(voleq)[1:3]); ii == 0 && return 0f0
    ht67 = TT_R4C1[ii] * dbhob^TT_R4C2[ii] * tht^TT_R4C3[ii]
    buttcf = TT_R4C5[ii] * dbhob + TT_R4C4[ii]
    stumpd = sqrt(buttcf * buttcf * tht / (tht - 4f0))
    d67 = TT_R4C7[ii] * dbhob * (2f0 / 3f0) + TT_R4C6[ii]
    cf0 = 0.002727f0 * (ht67 * stumpd * stumpd + d67 * d67 * tht)
    f = cf0 / (0.005454f0 * stumpd * stumpd * tht)
    b = (1f0 - f) / (2f0 * f)
    trm = 0.5f0
    topdia = mtopp < 6f0 ? 6f0 : mtopp                     # M=1: board top ≥ 6" (r4vol.f:212)
    topdia >= stumpd && return 0f0                          # M=1 gate vs STUMPD
    dratio = topdia / stumpd; dratio <= 0f0 && (dratio = 0.0001f0)
    merlen = tht - tht * dratio^(1f0 / b)
    (ht1prd > 0f0 && ht1prd < merlen) && (merlen = ht1prd)
    merlen < 2.5f0 && return 0f0
    totlgs = merlen / 16.5f0
    toplen = (totlgs - _fint(totlgs)) * 16.5f0
    if toplen < 2.5f0
        toplen = totlgs >= 1f0 ? 16.5f0 : 0f0; totlgs = Float32(_fint(totlgs))
    else
        set = false
        for L in 4:2:16
            if toplen < Float32(L) + trm; toplen = Float32(L - 2) + trm; set = true; break; end
        end
        set || (toplen = 16.5f0)
    end
    totlgs <= 0f0 && return 0f0
    totlgs = toplen < 16.5f0 ? Float32(_fint(totlgs) + 1) : Float32(_fint(totlgs))
    merlen = (totlgs - 1f0) * 16.5f0 + toplen
    numlgs = _fint(totlgs)
    @inline dsm_at(hcut) = Float32(_fint(stumpd * ((tht - hcut) / tht)^b + 0.499f0))
    volr4 = 0f0; dsm_prev = 0f0
    if numlgs > 1
        dsm1 = dsm_at(16.5f0)
        volr4 += _scribc(_fint(dsm1 - 5f0), 8)             # SCRIBC(INT(DSM-5), INT(16/2)=8)
        dsm_prev = dsm1
        for num in 2:(numlgs - 1)
            dsm = dsm_at(16.5f0 * num); volr4 += _scribc(_fint(dsm - 5f0), 8); dsm_prev = dsm
        end
    end
    len_top = toplen - trm
    dsm_top = Float32(_fint(stumpd * ((tht - merlen) / tht)^b + 0.499f0))
    volr4 += _scribc(_fint(dsm_top - 5f0), _fint(len_top / 2f0))
    return volr4 * 10f0                                     # BFGRS = VOLR4(1)·10
end
