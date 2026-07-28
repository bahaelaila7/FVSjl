# =============================================================================
# cr_fw2_vol.jl — CR volume, FW2 method (NVEL Flewelling & Raynes stem profile).
#
# The 2-point (DBH + total height) Flewelling variable-form segmented taper, for the CR
# region-2/3 species (ponderosa/DF/LP/WF/ES/aspen — JSP 23-29). Chain (fwinit.f→profile.f):
#   FWINIT: VOLEQ → JSP (internal Flewelling species).  SF_SHP→SHP_OT (f_other.f): the F(:,jrsp)
#   coefficient table + DBH/HT → geometric form params RFLW(6)/RHFW(4).  SF_TAPER (sf_taper.f):
#   RHFW/RFLW → 12 taper-polynomial coefs TAPCOE.  Scale F = DBHIB / SF_YHAT(4.5/HT) so the profile
#   hits DBH at breast height.  SF_YHAT (sf_yhat.f): TAPCOE + relative height → dib.  TCUBIC
#   (profile.f:883): stump 1-ft cylinder + 4-ft Smalian sections up the stem; VOL(1)=NINT(TCVOL·10)/10.
#
# Precision: the Fortran does the form/taper math in DOUBLE, stores RFLW/RHFW/TAPCOE as REAL*4, and
# integrates TCUBIC in REAL(single) — mirrored here (Float64 kernels, Float32 at those boundaries).
# The 0.1-rounded VOL(1) gives slack. Only the module-free 2-point path is ported (no 3-pt/sf_zero,
# no upper-stem measurement — the standard FVS case). Small trees (HT≤15, FWSMALL) TODO.
# =============================================================================

include(joinpath(@__DIR__, "..", "..", "data", "centralrockies", "fw2", "fw2_coefs.jl"))
include(joinpath(@__DIR__, "..", "..", "data", "centralrockies", "fw2", "ingy_coefs.jl"))

"FWINIT (fwinit.f): map a FW2 VOLEQ to the internal Flewelling species JSP. CR uses GEOCODE 2/3
(regions 2/3). Returns 0 for an unrecognized/unsupported eq (⇒ no FW2 volume)."
function _fw2_jsp(voleq::AbstractString)
    length(voleq) < 10 && return 0
    geocode = voleq[1]; geosub = voleq[2:3]; spec = voleq[8:10]
    if geocode == '2'
        spec == "122" && return geosub == "03" ? 22 : 23
        spec == "108" && return 25
        spec == "202" && return 26
        spec == "015" && return 27
        spec == "746" && return 28
    elseif geocode == '3'
        if geosub == "00"
            spec == "122" && return 29
            spec == "202" && return 26
        elseif geosub == "01"
            spec == "122" && return 29
            spec == "108" && return 25
            spec == "202" && return 26
            spec == "015" && return 27
        end
    elseif geocode == '4'                         # region 4 (fwinit.f): GEOSUB 07 → R2 profiles + R4 bark
        if geosub == "07"
            spec == "093" && return 24            # Engelmann spruce (Dixie ES model)
            spec == "122" && return 23            # R2 ponderosa (with R4 bark)
        end
    elseif geocode == 'I' || geocode == 'i' || geocode == '1'   # INGY / east-side (fwinit.f)
        (spec == "202" || spec == "205" || spec == "204") && return 11   # Douglas-fir
        (spec == "073" || spec == "070") && return 12   # Western larch
        spec == "017" && return 13                # Grand fir
        spec == "122" && return 14                # Ponderosa pine
        spec == "108" && return 15                # Lodgepole pine
        (spec == "242" || spec == "240") && return 16   # Western red cedar
        (spec == "260" || spec == "263" || spec == "264") && return 17  # Mountain hemlock
        spec == "119" && return 18                # White pine
        (spec == "093" || spec == "090") && return 19   # Engelmann spruce
        spec == "019" && return 20                # Subalpine fir
        spec == "012" && return 21                # Balsam fir
    end
    return 0
end

"SHP_OT/SHP_C2 (f_other.f / f_ingy.f): Flewelling geometric form params from an F-coefficient row.
Shared by region-2/3 (SHP_OT) and INGY (SHP_C2) — same regression, `is_lp` picks the lodgepole U3 form.
Returns (RFLW::NTuple{6,Float32}, RHFW::NTuple{4,Float32}). Double-precision kernel, REAL*4 storage."
function _fw2_shp_core(f, d::Float32, h::Float32, is_lp::Bool)
    @inline ff(i) = f[i - 9]                  # f[k] == Fortran F(k+9)
    D = Float64(d); H = Float64(h); lnH = log(H)
    dmedian = ff(10) * (H - 4.5)^(ff(11) + ff(12) * H)
    dform = D / dmedian - 1.0
    u7 = ff(13) + ff(14) * lnH + ff(15) * dform
    u9t = ff(18) + ff(19) * lnH + ff(20) * dform
    u9t = clamp(u9t, -7.0, 7.0)
    u9 = ff(16) * exp(u9t) / (1.0 + exp(u9t))
    u8 = ff(21) + ff(22) * H + ff(23) * lnH + ff(24) * dform
    u1 = ff(25) + ff(26) * lnH + ff(27) * dform + ff(28) * dform * lnH
    u2 = ff(29) + ff(30) * dform + ff(31) * lnH + ff(32) * dform * lnH + ff(33) * D
    u3 = is_lp ?
        ff(34) + ff(35) * dform + ff(36) * (1.0 - exp(ff(37) * H)) :
        ff(34) + ff(35) * dform + ff(36) * lnH + ff(37) * lnH * dform
    u4 = ff(38) + ff(39) * dform + ff(40) * lnH + ff(41) * D
    u5 = ff(42) + ff(43) * lnH
    u6 = ff(45) + ff(46) * dform + ff(47) * lnH
    u1 = clamp(u1, -7.0, 7.0); u2 = clamp(u2, -7.0, 7.0)
    u3 = clamp(u3, -7.0, 7.0); u4 = clamp(u4, -7.0, 7.0)
    u5 = clamp(u5, -7.0, 7.1)
    u6 = u6 < 1.005 ? 1.005 : (u6 > 10.0 ? 10.0 : u6)
    u7 = clamp(u7, -7.0, 7.0)
    u8 > 0.99 && (u8 = 0.99)
    u9 = u9 > 0.3 ? 0.3 : (u9 < 0.0 ? 0.0 : u9)
    r1 = exp(u1) / (1.0 + exp(u1)); r2 = exp(u2) / (1.0 + exp(u2))
    r3 = exp(u3) / (1.0 + exp(u3)); r4 = exp(u4) / (1.0 + exp(u4))
    r5 = u5 <= 7.0 ? 0.5 + 0.5 * exp(u5) / (1.0 + exp(u5)) : 1.0
    a3 = u6
    rhi1 = exp(u7) / (1.0 + exp(u7)); rhi1 > 0.5 && (rhi1 = 0.5)
    rhlongi = u9
    rhi2 = rhi1 + rhlongi
    rhc = u8
    rhc < rhi2 + 0.01 && (rhc = min(rhi2 + 0.01, (rhi2 + 1.0) / 2.0))
    rflw = (Float32(r1), Float32(r2), Float32(r3), Float32(r4), Float32(r5), Float32(a3))
    rhfw = (Float32(rhi1), Float32(rhi2), Float32(rhc), Float32(rhlongi))
    return rflw, rhfw
end

@inline _fw2_is_ingy(jsp::Int) = 11 <= jsp <= 21

"Form params for a JSP: region-2/3 (SHP_OT, F=_FW2_F[jsp-22]) or INGY (SHP_C2, F=_FW2_F_INGY[jsp-10])."
function _fw2_shp(jsp::Int, d::Float32, h::Float32)
    if _fw2_is_ingy(jsp)
        return _fw2_shp_core(_FW2_F_INGY[jsp - 10], d, h, jsp == 15)   # INGY; JSP15 = lodgepole
    else
        return _fw2_shp_core(_FW2_F[jsp - 22], d, h, jsp == 25)        # region 2/3; JSP25 = R2 lodgepole
    end
end

"FDBT_C2 (f_ingy.f): double bark thickness at breast height for INGY JSP 11-21 (GEOSUB '00' path).
Used to get DBHIB = DBHOB - DBTBH — the INGY profile is calibrated inside bark."
function _fw2_fdbt_c2(jsp::Int, d::Float32, h::Float32)::Float32
    a = _FW2_FDBT_A[jsp - 10]                 # (a1,a2,a3,a4,a5,a6)
    D = Float64(d); H = Float64(h)
    a00 = a[1]                                # GEOSUB '00' ⇒ A00 = A(1,jspr)
    duse = D
    if a[3] < 0.0 && a[6] == 0.0
        dmax = -(a[2] + a[5] * H) / (2.0 * a[3]); duse = min(D, dmax)
    elseif a[3] < 0.0 && a[6] == -1.0
        dmin = -(a[2] + a[5] * H) / (2.0 * a[3]); duse = max(D, dmin)
    end
    y2 = a00 + a[2] * duse + a[3] * duse * duse + a[4] * H + a[5] * H * duse
    y2 = clamp(y2, -8.0, 8.0)
    ratio = exp(y2) / (1.0 + exp(y2))
    return Float32(ratio * D)
end

"SF_TAPER (sf_taper.f): RHFW/RFLW → the 12 taper-polynomial coefficients TAPCOE (REAL*4)."
function _fw2_sf_taper(rhfw, rflw)
    r1 = Float64(rflw[1]); r2 = Float64(rflw[2]); r3 = Float64(rflw[3])
    r4 = Float64(rflw[4]); r5 = Float64(rflw[5]); a3 = Float64(rflw[6])
    rhi1 = Float64(rhfw[1]); rhi2 = Float64(rhfw[2]); rhc = Float64(rhfw[3]); rhlongi = Float64(rhfw[4])
    k = 1.0
    yc = k * (1.0 - rhc)
    c2 = r5 * yc; c1 = 3.0 * (yc - c2); slope = -(3.0 - r5) * k / 2.0
    s1 = slope * (rhc - rhi2)
    yi_min = yc - s1 * (1.0 + 2.0 * r3) / 3.0
    yi_max = yc - s1 * (5.0 + 4.0 * r3) / 9.0
    yi2 = yi_min + r4 * (yi_max - yi_min)
    s0 = r3 * s1
    b1 = (6.0 * yc - 6.0 * yi2 - 2.0 * s0 - 4.0 * s1) / (-3.0 * yc + 3.0 * yi2 + 2.0 * s0 + s1)
    b2 = s1 * (1.0 - r3) / (0.5 - 1.0 / (b1 + 1.0))
    b4 = s0; b0 = yi2
    slope_rhi = r3 * s1 / (rhc - rhi2)
    yi1 = yi2 - slope_rhi * rhlongi
    if rhlongi > 0.0
        e2 = (yi2 - yi1) / rhlongi; e1 = yi1 - e2 * rhi1
    else
        e1 = yi2; e2 = 0.0
    end
    s3 = -slope_rhi * rhi1; k2 = s3 / r1
    f_a3 = 1.0 / (6.0 * a3 * a3) + log(1.0 - 1.0 / a3) + 1.0 / (3.0 * (a3 - 1.0)) + 2.0 / (3.0 * a3)
    g_a3 = (1.0 / (a3 - 1.0) - 1.0 / a3 - 1.0 / (a3 * a3) - 1.0 / (a3 - 1.0)^3) / f_a3
    yb_min = yi1 + (2.0 * s3 + k2) / 3.0 + (s3 - k2) * f_a3 /
             (1.0 / (a3 - 1.0) - 1.0 / a3 - 1.0 / (a3^2) - 1.0 / (a3 * a3 * a3))
    yb_max = yi1 + (2.0 * s3 + k2) / 3.0 + (s3 - k2) / g_a3
    yb = yb_min + r2 * (yb_max - yb_min)
    a0 = yi1
    a2 = (yb - yi1 - (2.0 * s3 + k2) / 3.0) / f_a3
    a1 = (k2 - s3 + a2 * (1.0 / (a3 - 1.0) - 1.0 / a3 - 1.0 / a3^2)) / 3.0
    a4 = s3
    return (Float32(a0), Float32(a1), Float32(a2), Float32(a4), Float32(b0), Float32(b1),
            Float32(b2), Float32(b4), Float32(c1), Float32(c2), Float32(e1), Float32(e2))
end

"SF_YHAT (sf_yhat.f, JSP≠22): profile dib at relative height rh. Coefs REAL*4, x/y REAL*8, result REAL*4."
function _fw2_sf_yhat(rh::Float32, tapcoe, rhfw, rflw, f::Float32)::Float32
    rh > 1.0f0 && return 0.0f0
    rh < 0.0f0 && return f
    r3 = rflw[3]; a3 = Float64(rflw[6])
    rhi1 = rhfw[1]; rhi2 = rhfw[2]; rhc = rhfw[3]; rhlongi = rhfw[4]
    a0 = Float64(tapcoe[1]); a1 = Float64(tapcoe[2]); a2 = Float64(tapcoe[3]); a4 = Float64(tapcoe[4])
    b0 = Float64(tapcoe[5]); b1 = Float64(tapcoe[6]); b2 = Float64(tapcoe[7]); b4 = Float64(tapcoe[8])
    c1 = Float64(tapcoe[9]); c2 = Float64(tapcoe[10])
    e1 = Float64(tapcoe[11]); e2 = Float64(tapcoe[12])
    R = Float64(rh)
    local y::Float64
    if rh >= rhc                                    # upper segment
        x = (1.0 - R) / (1.0 - Float64(rhc))
        y = x * (c2 + x * ((c1 / 2.0) - (c1 / 6.0) * x))
    elseif rh >= rhi2                               # middle segment
        x = (R - Float64(rhi2)) / (Float64(rhc) - Float64(rhi2))
        if x > 0.0
            sus2 = (b1 * log10(x) <= -20.0) ? 0.0 : x^b1
            y = b0 + x * (b4 + x * (-b2 / ((b1 + 1.0) * (b1 + 2.0)) * sus2 + b2 / 6.0 * x))
        else
            y = b0
        end
    elseif rhlongi > 0.0f0 && rh > rhi1             # straight segment
        y = e1 + e2 * R
    else                                            # lower segment
        x = (Float64(rhi1) - R) / Float64(rhi1)
        y = a0 + x * ((a4 + a2 / a3) + x * (a2 / (2.0 * a3 * a3) + a1 * x)) + a2 * log(1.0 - x / a3)
    end
    return Float32(Float64(f) * y)
end

# BRK_OT (f_other.f) bark model, JSPR=JSP-21 columns 1-8. Rows 1-6 = B2,B3,B5,C1,C2,C3 (the
# double-bark-thickness proportion PY; the D0/D1/D2 rows 7-9 fit DBHIB when DBTBH≤0 — unused here
# since cr_bratio always gives DBTBH=D·(1-bark)>0). JSPR 8 = R3 ponderosa; 2=San Juan, 4=R2 LP, etc.
const _FW2_BK = (
    (-0.310745, -5.267465, 4.056924, 1.1159037603, -0.082066096, 0.424652164),        # 1 Black Hills PP
    (12.88990159, 17.90876955, 0.05237532, 0.0145082565, 0.0027058753, 0.1022613683), # 2 San Juan PP
    (-0.8527459114, 0.9336248438, 0.2809059946, 0.0029765304, 0.0535125376, 0.3088136745), # 3 Dixie ES
    (-4.053133738, -3.891047743, 0.380866177, 0.7617392463, 0.2071264529, 0.6491861037),   # 4 R2 Lodgepole
    (2.0518, 0.5569, 0.2009, 1.1163, 0.2194, 0.2525),          # 5 R2 Douglas-fir
    (-0.3114, -0.5979, -0.00029, 0.4661, 0.2448, 0.1619),      # 6 R2 White fir
    (20.1848, 111.927, 0.0808502, 2.24620, 0.336067, 0.486677),# 7 R2 Aspen
    (2.67998, 3.49916, 0.111418, 2.28510, 0.374458, 0.450724), # 8 R3 Ponderosa pine
)

"BRK_OT (f_other.f): inside-bark dib at height ht2 from the profile's outside-bark diameter `dob`,
given double-bark-thickness `dbtbh` (>0 from cr_bratio ⇒ the model-DBHIB branch is skipped)."
@inline function _fw2_brk_ot(jsp::Int, dbhob::Float32, dob::Float32, ht2::Float32, dbtbh::Float32)::Float32
    bk = _FW2_BK[jsp - 21]
    b2 = bk[1]; b3 = bk[2]; b5 = bk[3]; c1 = bk[4]; c2 = bk[5]; c3 = bk[6]
    DOB = Float64(dob); DBT = Float64(dbtbh)
    dr = DOB > 0.0 ? DOB / Float64(dbhob) : 0.0
    local py::Float64
    if ht2 > 4.5f0
        py = dr > 0.01 ? (dr * ((b2 - 1.0) / (b2 - dr^b3)) - ((dr^b5 - 1.0) / DBT)) : 0.0
    elseif ht2 == 4.5f0
        py = 1.0
    else
        clx = c1 * (dr - 1.0)
        py = clx >= 0.0 ? 1.0 + clx^(c2 + c3 * DBT) : 1.0
    end
    dib = DOB - py * DBT
    return dib < 0.0 ? 0.0f0 : Float32(dib)
end

"FWSMALL (profile.f): corrected stump dib at 1 ft for small trees (HTTOT≤15), `dib_at_1`=profile dib
at 1 ft, `dbhib`=D·BARK. JSP14 (INGY PP) uses different height-ratio bounds."
function _fw2_fwsmall(jsp::Int, h::Float32, dib_at_1::Float32, dbhib::Float32)::Float32
    hr5 = 0.18f0; hr15lo = 0.11f0; hr15hi = 0.25f0
    if jsp == 14
        hr5 = 0.25f0; hr15lo = 0.21f0; hr15hi = 0.27f0
    end
    hratio = dib_at_1 / (h - 1.0f0)
    hr_min = hr5 + (h - 5.0f0) / 10.0f0 * (hr15lo - hr5)
    hr_max = hr5 + (h - 5.0f0) / 10.0f0 * (hr15hi - hr5)
    hratio = clamp(hratio, hr_min, hr_max)
    dib1 = hratio * (h - 1.0f0)
    dr_min = h < 15.0f0 ? 1.0f0 + 0.3f0 * (15.0f0 - h) / 10.0f0 : 1.0f0
    (dib1 / dbhib) < dr_min && (dib1 = dr_min * dbhib)
    return dib1
end

"TCUBIC (profile.f:883): total cubic via stump cylinder + 4-ft Smalian sections. `dibat(ht)` gives the
inside-bark section diameter (region-2/3: SF_YHAT then BRK_OT; INGY: SF_YHAT directly — already ib).
`stump_dib` (small trees, HTTOT≤15) overrides the 1-ft stump diameter with the FWSMALL value."
function _fw2_tcubic(dibat, h::Float32; stump_dib::Float32 = -1f0)::Float32
    htloop = trunc(Int, (h + 0.5f0 - 1.0f0) / 4.0f0)
    ht2 = 1.0f0
    dib = stump_dib >= 0f0 ? stump_dib : dibat(ht2)       # dib at 1 ft (FWSMALL for small trees)
    r = dib / 2.0f0
    tcvol = (3.1416f0 * r * r) / 144.0f0                  # 1-ft stump cylinder
    @inbounds for _ in 1:htloop
        d2old = dib
        ht2 += 4.0f0
        dib = dibat(ht2)
        tcvol += 0.00272708f0 * (d2old * d2old + dib * dib) * 4.0f0
    end
    if (h - ht2) > 0.0f0                                  # tip (dib→0)
        tcvol += 0.00272708f0 * (dib * dib) * (h - ht2)
    end
    return tcvol
end

"SF_HS surrogate: bisection for the height where the inside-bark dib == `topd`. SF_HS itself is a
Newton+bisection solver to TOL=0.0005·H — below the 1-ft NUMLOG rounding, so a tight bisection matches."
function _fw2_hs(dibat, topd::Float32, h::Float32)::Float32
    lo = 0f0; hi = h                           # dib decreases monotonically with height
    @inbounds for _ in 1:80
        mid = (lo + hi) / 2f0
        dd = dibat(mid)
        abs(dd - topd) < 0.0005f0 && return mid
        dd > topd ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2f0
end

@inline function _fw2_dclass(x::Float32)::Float32   # GETDIB diameter class: INT(x), +1 if frac>0.501
    c = floor(Int, x)
    return Float32((x - Float32(c)) > 0.501f0 ? c + 1 : c)
end

"FW2 merch cubic VOL(4) (profile.f MERLEN→NUMLOG/SEGMNT→GETDIB→CUPFLG loop): buck stump→merch-top,
sum per-log Smalian .00272708·(DIBL²+DIBS²)·LEN with inch-class DIBs (butt = dib at breast height),
each log 0.1-rounded. `mtop` = inside-bark merch top = TOPD·BARK."
function _fw2_merch_cuft(dibat, h::Float32, mtop::Float32, stump::Float32)::Float32
    hs = _fw2_hs(dibat, mtop, h)
    lmerch = hs - stump
    lmerch < _NVB_R3_MERCHL && return 0f0
    numseg = _nvb_numlog(_NVB_R3_OPT, _NVB_R3_EVOD, lmerch, _NVB_R3_MAXLEN, _NVB_R3_MINLEN, _NVB_R3_TRIM)
    numseg == 0 && return 0f0
    loglen, numseg = _nvb_segmnt(_NVB_R3_OPT, _NVB_R3_EVOD, lmerch, _NVB_R3_MAXLEN, _NVB_R3_MINLEN, _NVB_R3_TRIM, numseg)
    dibl = _fw2_dclass(dibat(4.5f0))           # butt log large end = DIB class at breast height
    ht2 = stump; vol4 = 0f0
    @inbounds for i in 1:numseg
        ht2 += _NVB_R3_TRIM + loglen[i]
        dib = dibat(ht2)
        (i == numseg && dib < mtop) && (dib = mtop)   # GETDIB forces the top log ≥ MTOPP
        dibs = _fw2_dclass(dib)
        logv = 0.00272708f0 * (dibl * dibl + dibs * dibs) * loglen[i]
        vol4 += floor(logv * 10f0 + 0.5f0) / 10f0
        dibl = dibs
    end
    return vol4
end

"FW2 Scribner board VOL(2) (profile.f BFPFLG loop): buck stump→board-top (BFTOPD·BARK) and sum
SCRIB(small-end inch-class dib, len)·10 per log. Same region-3 log-bucking as the cubic."
function _fw2_board(dibat, h::Float32, bftop::Float32, stump::Float32)::Float32
    hs = _fw2_hs(dibat, bftop, h)
    lmerch = hs - stump
    lmerch < _NVB_R3_MERCHL && return 0f0
    numseg = _nvb_numlog(_NVB_R3_OPT, _NVB_R3_EVOD, lmerch, _NVB_R3_MAXLEN, _NVB_R3_MINLEN, _NVB_R3_TRIM)
    numseg == 0 && return 0f0
    loglen, numseg = _nvb_segmnt(_NVB_R3_OPT, _NVB_R3_EVOD, lmerch, _NVB_R3_MAXLEN, _NVB_R3_MINLEN, _NVB_R3_TRIM, numseg)
    ht2 = stump; vol2 = 0f0
    @inbounds for i in 1:numseg
        ht2 += _NVB_R3_TRIM + loglen[i]
        dib = dibat(ht2)
        (i == numseg && dib < bftop) && (dib = bftop)
        dibs = _fw2_dclass(dib)                    # small-end inch class
        vol2 += _scrib(dibs, loglen[i], 'Y') * 10f0   # COR='Y' ⇒ ×10 (decimal-C → bdft)
    end
    return vol2
end

"CR FW2 per-tree volume (2-point Flewelling), region-2/3 (JSP 23-29) + INGY (JSP 11-21). `bark`=DIB/DOB
(cr_bratio); `topd`/`bftopd`=cubic/board top DOB. Returns a 15-vec: VOL[1]=total cubic, VOL[4]=merch
cubic, VOL[2]=Scribner board feet (all 0.1-rounded where applicable). Merch tops are inside bark (·BARK).

INGY vs region-2/3 differ in bark handling: region-2/3 profiles are OUTSIDE bark (calibrated to DBHOB),
so the section diameter is SF_YHAT reduced by BRK_OT (DBTBH=D·(1-bark)); INGY profiles are INSIDE bark
(calibrated to DBHIB=DBHOB-FDBT_C2), so the section diameter is SF_YHAT directly (BRK_UP only adds bark
for DOB, not needed for cubic)."
function cr_fw2_vol(voleq::AbstractString, d::Float32, h::Float32;
                    bark::Float32 = 1f0, topd::Float32 = 4f0, stump::Float32 = 1f0, bftopd::Float32 = 6f0)
    vol = zeros(Float32, 15)
    (d < 1f0 || h <= 5f0) && return vol
    jsp = _fw2_jsp(voleq)
    (_fw2_is_ingy(jsp) || (23 <= jsp <= 29)) || return vol   # supported 2-pt families
    ingy = _fw2_is_ingy(jsp)
    rflw, rhfw = _fw2_shp(jsp, d, h)
    tapcoe = _fw2_sf_taper(rhfw, rflw)
    # INGY: profile is inside bark, calibrated to DBHIB. SF_SHP uses the PASSED DBTBH (fvsvol DBTBH=D·(1-BARK))
    # when >0 — so DBHIB=D·BARK (cr_bratio), NOT FDBT_C2 (that's only the no-bark-input fallback). Region-2/3: DBHOB.
    dbhib = ingy ? d * bark : d
    yhat_bh = _fw2_sf_yhat(4.5f0 / h, tapcoe, rhfw, rflw, 1.0f0)
    yhat_bh == 0f0 && return vol
    f = dbhib / yhat_bh
    dbtbh = d * (1f0 - bark)
    dibat = ingy ? (ht -> _fw2_sf_yhat(ht / h, tapcoe, rhfw, rflw, f)) :
                   (ht -> _fw2_brk_ot(jsp, d, _fw2_sf_yhat(ht / h, tapcoe, rhfw, rflw, f), ht, dbtbh))
    # Small trees (HTTOT≤15): FWSMALL corrects the stump diameter; merch/board stay 0 (LMERCH<MERCHL).
    stump_dib = h <= 15f0 ? _fw2_fwsmall(jsp, h, dibat(1.0f0), d * bark) : -1f0
    vol[1] = Float32(round(_fw2_tcubic(dibat, h; stump_dib = stump_dib) * 10.0f0)) / 10.0f0    # NINT(TCVOL*10)/10
    vol[4] = _fw2_merch_cuft(dibat, h, topd * bark, stump)
    vol[2] = _fw2_board(dibat, h, bftopd * bark, stump)
    return vol
end
