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
    end
    return 0
end

"SHP_OT (f_other.f): Flewelling geometric form params for JSP 23-29 (jrsp=JSP-22). Returns
(RFLW::NTuple{6,Float32}, RHFW::NTuple{4,Float32}). Double-precision kernel, REAL*4 storage."
function _fw2_shp_ot(jsp::Int, d::Float32, h::Float32)
    jrsp = jsp - 22
    f = _FW2_F[jrsp]                          # f[k] == Fortran F(k+9, jrsp)
    @inline ff(i) = f[i - 9]
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
    u3 = jrsp == 3 ?
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

"TCUBIC (profile.f:883): total cubic via stump cylinder + 4-ft Smalian sections. Each section diameter
is the profile's outside-bark dib (SF_YHAT) reduced to inside bark by BRK_OT (TAPERMODEL flow)."
function _fw2_tcubic(jsp::Int, d::Float32, h::Float32, tapcoe, rhfw, rflw, f::Float32, dbtbh::Float32)::Float32
    @inline dibat(ht) = _fw2_brk_ot(jsp, d, _fw2_sf_yhat(ht / h, tapcoe, rhfw, rflw, f), ht, dbtbh)
    htloop = trunc(Int, (h + 0.5f0 - 1.0f0) / 4.0f0)
    ht2 = 1.0f0
    dib = dibat(ht2)                                      # dib at 1 ft
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

"CR FW2 per-tree volume. `bark`=DIB/DOB ratio (cr_bratio) → DBTBH=D·(1-bark) for the BRK_OT bark
reduction. Returns a 15-vec with VOL[1]=total cubic (0.1-rounded). MCF/board TODO."
function cr_fw2_vol(voleq::AbstractString, d::Float32, h::Float32; bark::Float32 = 1f0)
    vol = zeros(Float32, 15)
    (d < 1f0 || h <= 5f0) && return vol
    jsp = _fw2_jsp(voleq)
    (jsp < 23 || jsp > 29) && return vol      # only the module-free region-2/3 2-pt path (JSP 23-29)
    h <= 15f0 && return vol                   # FWSMALL small-tree path TODO
    rflw, rhfw = _fw2_shp_ot(jsp, d, h)
    tapcoe = _fw2_sf_taper(rhfw, rflw)
    dbhib = d                                 # JSP 22-30: profile calibrated to DBHOB (SF_SHP DBHIB=DBHOB)
    yhat_bh = _fw2_sf_yhat(4.5f0 / h, tapcoe, rhfw, rflw, 1.0f0)
    yhat_bh == 0f0 && return vol
    f = dbhib / yhat_bh
    dbtbh = d * (1f0 - bark)                   # double bark thickness (fvsvol.f: DBTBH=D*(1-BARK))
    tcvol = _fw2_tcubic(jsp, d, h, tapcoe, rhfw, rflw, f, dbtbh)
    vol[1] = Float32(round(tcvol * 10.0f0)) / 10.0f0      # NINT(TCVOL*10)/10
    return vol
end
