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
