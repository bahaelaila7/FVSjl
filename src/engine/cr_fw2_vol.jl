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
    # SHP_OT / SHP_C2 precision (f_other.f / f_ingy.f, IMPLICIT DOUBLE PRECISION but DBH/HT are REAL*4):
    #  * log(HT) of a REAL*4 argument is the SINGLE-precision generic LOG (gfortran logf), widened afterwards;
    #  * (HT-4.5) is a REAL*4 subtraction (4.5 is a default-REAL constant) before promotion to the REAL*8 power.
    D = Float64(d); H = Float64(h); lnH = Float64(flog(h))
    dmedian = ff(10) * dpow(Float64(h - 4.5f0), ff(11) + ff(12) * H)
    dform = D / dmedian - 1.0
    u7 = ff(13) + ff(14) * lnH + ff(15) * dform
    u9t = ff(18) + ff(19) * lnH + ff(20) * dform
    u9t = clamp(u9t, -7.0, 7.0)
    u9 = ff(16) * dexp(u9t) / (1.0 + dexp(u9t))
    u8 = ff(21) + ff(22) * H + ff(23) * lnH + ff(24) * dform
    u1 = ff(25) + ff(26) * lnH + ff(27) * dform + ff(28) * dform * lnH
    u2 = ff(29) + ff(30) * dform + ff(31) * lnH + ff(32) * dform * lnH + ff(33) * D
    u3 = is_lp ?
        ff(34) + ff(35) * dform + ff(36) * (1.0 - dexp(ff(37) * H)) :
        ff(34) + ff(35) * dform + ff(36) * lnH + ff(37) * lnH * dform
    u4 = ff(38) + ff(39) * dform + ff(40) * lnH + ff(41) * D
    u5 = ff(42) + ff(43) * lnH
    u6 = ff(45) + ff(46) * dform + ff(47) * lnH
    u1 = clamp(u1, -7.0, 7.0); u2 = clamp(u2, -7.0, 7.0)
    u3 = clamp(u3, -7.0, 7.0); u4 = clamp(u4, -7.0, 7.0)
    u5 = clamp(u5, -7.0, 7.1)
    # `U6 .lt. 1.005` / `U8 .gt. 0.99` compare against default-REAL (single) constants widened to double; the
    # U8 reset value 0.99 is likewise the widened single constant.
    u6 = u6 < Float64(1.005f0) ? 1.005 : (u6 > 10.0 ? 10.0 : u6)
    u7 = clamp(u7, -7.0, 7.0)
    u8 > Float64(0.99f0) && (u8 = Float64(0.99f0))
    u9 = u9 > 0.3 ? 0.3 : (u9 < 0.0 ? 0.0 : u9)
    r1 = dexp(u1) / (1.0 + dexp(u1)); r2 = dexp(u2) / (1.0 + dexp(u2))
    r3 = dexp(u3) / (1.0 + dexp(u3)); r4 = dexp(u4) / (1.0 + dexp(u4))
    r5 = u5 <= 7.0 ? 0.5 + 0.5 * dexp(u5) / (1.0 + dexp(u5)) : 1.0
    a3 = u6
    # RHI1/RHLONGI/RHI2/RHC are REAL*4: each double is rounded on assignment, then RHI2=RHI1+RHLONGI and the RHC
    # test/min (RHI2+.01, (RHI2+1.)/2.0) are SINGLE-precision arithmetic on those rounded values.
    rhi1 = Float32(dexp(u7) / (1.0 + dexp(u7))); rhi1 > 0.5f0 && (rhi1 = 0.5f0)
    rhlongi = Float32(u9)
    rhi2 = rhi1 + rhlongi
    rhc = Float32(u8)
    rhc < rhi2 + 0.01f0 && (rhc = min(rhi2 + 0.01f0, (rhi2 + 1.0f0) / 2.0f0))
    rflw = (Float32(r1), Float32(r2), Float32(r3), Float32(r4), Float32(r5), Float32(a3))
    rhfw = (rhi1, rhi2, rhc, rhlongi)
    return rflw, rhfw
end

# SHP_BH (f_other.f:682) — Black Hills NF ponderosa (JSP=22) Flewelling shape: distinct hardcoded PP14
# coefficients (U7 uses H not lnH; U6 capped at 100 not 10), same U→R tail as SHP_OT. Separate routine
# because JRSP=JSP-22=0 has no F(:,0) column. Ported so Black Hills ponderosa stands (voleq 203FW2W..122,
# geosub 03) get volume instead of 0.
function _fw2_shp_bh(d::Float32, h::Float32)
    D = Float64(d); H = Float64(h); lnH = dlog(H)
    dmedian = 1.6802 * dpow(H - 4.5, 0.4085 + 0.00169 * H)
    dform = D / dmedian - 1.0
    u7 = -1.2726446 - 0.0048259438 * H
    u9 = 0.1821947
    u8 = 0.99
    u1 = -1.5505171 - 0.017174522 * H
    u2 = 0.27722769 - 0.21540189 * D
    u3 = 2.0426515 - 0.83434213 * lnH
    u4 = -7.0
    u5 = 7.7448837
    u6 = 1.376637 - 0.47598661 * dform
    u1 = clamp(u1, -7.0, 7.0); u2 = clamp(u2, -7.0, 7.0)
    u3 = clamp(u3, -7.0, 7.0); u4 = clamp(u4, -7.0, 7.0)
    u5 = u5 < -7.0 ? -7.0 : (u5 > 7.1 ? 7.1 : u5)
    u6 = u6 < 1.005 ? 1.005 : (u6 > 100.0 ? 100.0 : u6)   # SHP_BH caps U6 at 100 (SHP_OT at 10)
    u7 = clamp(u7, -7.0, 7.0)
    u8 > 0.99 && (u8 = 0.99)
    u9 = u9 > 0.3 ? 0.3 : (u9 < 0.0 ? 0.0 : u9)
    r1 = dexp(u1) / (1.0 + dexp(u1)); r2 = dexp(u2) / (1.0 + dexp(u2))
    r3 = dexp(u3) / (1.0 + dexp(u3)); r4 = dexp(u4) / (1.0 + dexp(u4))
    r5 = u5 <= 7.0 ? 0.5 + 0.5 * dexp(u5) / (1.0 + dexp(u5)) : 1.0
    a3 = u6
    rhi1 = dexp(u7) / (1.0 + dexp(u7)); rhi1 > 0.5 && (rhi1 = 0.5)
    rhlongi = u9; rhi2 = rhi1 + rhlongi; rhc = u8
    rhc < rhi2 + 0.01 && (rhc = min(rhi2 + 0.01, (rhi2 + 1.0) / 2.0))
    rflw = (Float32(r1), Float32(r2), Float32(r3), Float32(r4), Float32(r5), Float32(a3))
    rhfw = (Float32(rhi1), Float32(rhi2), Float32(rhc), Float32(rhlongi))
    return rflw, rhfw
end

@inline _fw2_is_ingy(jsp::Int) = 11 <= jsp <= 21

"SHP_C2 (f_ingy.f) GEOSUB subregion coefficient override: return the (possibly modified) INGY
F-coefficient row for `jsp`/`geosub`. geosub \"00\" (or an unmatched subregion) leaves the base row
unchanged; matched subregions replace individual F(25)/F(34)/F(38) coefficients (JSP 12/13/14), or the
whole column with F(:,23) for grand fir geosub 15/03. Faithful to SHP_C2 lines 276-351."
function _fw2_ingy_frow(jsp::Int, geosub::AbstractString)
    base = _FW2_F_INGY[jsp - 10]
    geosub == "00" && return base
    local ovr
    if jsp == 11
        (geosub == "15" || geosub == "03") || return base             # DF: geotemp forced '15'
        ovr = _FW2_DFSUB
    elseif jsp == 19
        geosub == "15" || return base                                 # ES
        ovr = _FW2_ESSUB
    elseif jsp == 13                                                  # grand/white fir
        (geosub == "15" || geosub == "03") && return _FW2_F_INGY_GF23 # whole-column swap, no GFSUB
        o = get(_FW2_GFSUB, geosub, nothing); o === nothing && return base
        ovr = o
    elseif jsp == 14
        o = get(_FW2_PPSUB, geosub, nothing); o === nothing && return base
        ovr = o
    elseif jsp == 12
        o = get(_FW2_WLSUB, geosub, nothing); o === nothing && return base
        ovr = o
    else
        return base                                                   # JSP 15-18,20,21: no subregion mods
    end
    frow = copy(base)
    @inbounds for (k, v) in ovr
        frow[k] = v
    end
    return frow
end

"Form params for a JSP: region-2/3 (SHP_OT, F=_FW2_F[jsp-22]) or INGY (SHP_C2, F=_FW2_F_INGY[jsp-10]
with the GEOSUB subregion override applied)."
function _fw2_shp(jsp::Int, geosub::AbstractString, d::Float32, h::Float32)
    if _fw2_is_ingy(jsp)
        return _fw2_shp_core(_fw2_ingy_frow(jsp, geosub), d, h, jsp == 15)  # INGY; JSP15 = lodgepole
    elseif jsp == 22
        return _fw2_shp_bh(d, h)                                       # Black Hills PP (SHP_BH)
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
    ratio = dexp(y2) / (1.0 + dexp(y2))
    return Float32(ratio * D)
end

"SF_TAPER (sf_taper.f): RHFW/RFLW → the 12 taper-polynomial coefficients TAPCOE (REAL*4)."
function _fw2_sf_taper(rhfw, rflw)
    # sf_taper.f: IMPLICIT DOUBLE PRECISION, R1..R5/A3 REAL*8 (widened from RFLW), but RHI1/RHI2/RHC/RHLONGI are
    # declared REAL*4 — so `1.0-RHC` and `RHC-RHI2` are SINGLE-precision subtractions (default-REAL constants),
    # widened only when they meet a REAL*8 operand. LOG and `**` on REAL*8 are glibc log/pow (one rounding).
    r1 = Float64(rflw[1]); r2 = Float64(rflw[2]); r3 = Float64(rflw[3])
    r4 = Float64(rflw[4]); r5 = Float64(rflw[5]); a3 = Float64(rflw[6])
    rhi1 = rhfw[1]; rhi2 = rhfw[2]; rhc = rhfw[3]; rhlongi = rhfw[4]            # REAL*4
    k = 1.0
    yc = k * Float64(1.0f0 - rhc)                                                 # :36 k*(1.0-rhc)
    c2 = r5 * yc; c1 = 3.0 * (yc - c2); slope = -(3.0 - r5) * k / 2.0
    s1 = slope * Float64(rhc - rhi2)                                              # :43
    yi_min = yc - s1 * (1.0 + 2.0 * r3) / 3.0
    yi_max = yc - s1 * (5.0 + 4.0 * r3) / 9.0
    yi2 = yi_min + r4 * (yi_max - yi_min)
    s0 = r3 * s1
    b1 = (6.0 * yc - 6.0 * yi2 - 2.0 * s0 - 4.0 * s1) / (-3.0 * yc + 3.0 * yi2 + 2.0 * s0 + s1)
    b2 = s1 * (1.0 - r3) / (0.5 - 1.0 / (b1 + 1.0))
    b4 = s0; b0 = yi2
    slope_rhi = r3 * s1 / Float64(rhc - rhi2)                                     # :54
    yi1 = yi2 - slope_rhi * Float64(rhlongi)
    if rhlongi > 0f0
        e2 = (yi2 - yi1) / Float64(rhlongi); e1 = yi1 - e2 * Float64(rhi1)
    else
        e1 = yi2; e2 = 0.0
    end
    s3 = -slope_rhi * Float64(rhi1); k2 = s3 / r1
    f_a3 = 1.0 / (6.0 * a3 * a3) + dlog(1.0 - 1.0 / a3) + 1.0 / (3.0 * (a3 - 1.0)) + 2.0 / (3.0 * a3)
    g_a3 = (1.0 / (a3 - 1.0) - 1.0 / a3 - 1.0 / (a3 * a3) - 1.0 / dpow(a3 - 1.0, 3.0)) / f_a3
    yb_min = yi1 + (2.0 * s3 + k2) / 3.0 + (s3 - k2) * f_a3 /
             (1.0 / (a3 - 1.0) - 1.0 / a3 - 1.0 / dpow(a3, 2.0) - 1.0 / (a3 * a3 * a3))
    yb_max = yi1 + (2.0 * s3 + k2) / 3.0 + (s3 - k2) / g_a3
    yb = yb_min + r2 * (yb_max - yb_min)
    a0 = yi1
    a2 = (yb - yi1 - (2.0 * s3 + k2) / 3.0) / f_a3
    a1 = (k2 - s3 + a2 * (1.0 / (a3 - 1.0) - 1.0 / a3 - 1.0 / dpow(a3, 2.0))) / 3.0
    a4 = s3
    return (Float32(a0), Float32(a1), Float32(a2), Float32(a4), Float32(b0), Float32(b1),
            Float32(b2), Float32(b4), Float32(c1), Float32(c2), Float32(e1), Float32(e2))
end

"""SF_YHAT (sf_yhat.f, JSP≠22): profile dib at relative height rh — the faithful mixed-precision kernel
`_fw2_sf_yhat_f`: X and the REAL coefficient sub-terms (C1/6, B2/((B1+1)(B1+2)), A4+A2/A3, A2/(2·A3²), E1+E2·RH)
are formed in REAL*4, Y in REAL*8, F·Y rounded to REAL. The all-Float64 kernel this replaced put GETDIB's log
small-end dibs 1 ULP off (MEASURED FVSie_g16 374547584489998 climate 2035 WL 12.05"×107.0': 83.5-ft dib 4.5010004 jl
vs 4.5009999 live, across GETDIB's 0.501 diameter-class cut ⇒ MCuFt 29.8/BdFt 154 vs live 29.3/148)."""
_fw2_sf_yhat(rh::Float32, tapcoe, rhfw, rflw, f::Float32)::Float32 =
    _fw2_sf_yhat_f(rh, tapcoe, rhfw, rflw, f, 1f0, false)[1]

"SF_YHAT with slope (sf_yhat.f, ineedsl=1, JSP≠22): returns (dib::Float32, dDIB/dH::Float32). The diameter is
identical to _fw2_sf_yhat (validated: max-diff 0 over rh∈[0,1]; slope matches the numerical derivative); the slope
adds the per-segment dy_dx (dd_dH = dy_dx·F/(RH_LENGTH·TOTALH), sign-flipped for segments 1/4). Middle-segment
slope uses SUS3=x^(b1+1) (≠ diameter's SUS2=x^b1).
WIP BUILDING BLOCK (currently UNUSED, inert) toward the deferred SF_HS-Newton _fw2_hs port (see audit 2026-07-30):
SF_DS (sf_ds.f, NEXTRA=0) = raw SF_YHAT (this slope), NO BRK_UP; SF_HS solves SF_DS==DIB_target on the SF_YHAT
profile with a HEIGHT-tol Newton (`H+=-(D-DIB)/SLOPE`, |ADJUST|≤TOL·TOTALH) — replacing jl's current diameter-tol
bisection (which stops ~0.1 ft early at the flat top ⇒ the MERCHL=10 cliff-flip). NOTE: the exact merch-top/bark convention SF_HS solves against is NOT yet pinned — live SF_DS shows the crossing
DIB≠the naive topd·bark, and hand-reconstruction of f/bark is unreliable; the focused port must instrument live
SF_DS/SF_HS end-to-end to fix mtop + the inside/outside-bark comparison before wiring the Newton."
function _fw2_sf_yhat_sl(rh::Float32, tapcoe, rhfw, rflw, f::Float32, totalh::Float32)
    rh > 1.0f0 && return (0.0f0, 0.0f0)
    rh < 0.0f0 && return (f, 0.0f0)
    a3 = Float64(rflw[6])
    rhi1 = rhfw[1]; rhi2 = rhfw[2]; rhc = rhfw[3]; rhlongi = rhfw[4]
    a0 = Float64(tapcoe[1]); a1 = Float64(tapcoe[2]); a2 = Float64(tapcoe[3]); a4 = Float64(tapcoe[4])
    b0 = Float64(tapcoe[5]); b1 = Float64(tapcoe[6]); b2 = Float64(tapcoe[7]); b4 = Float64(tapcoe[8])
    c1 = Float64(tapcoe[9]); c2 = Float64(tapcoe[10])
    e1 = Float64(tapcoe[11]); e2 = Float64(tapcoe[12])
    R = Float64(rh)
    local y::Float64, dy_dx::Float64, rh_length::Float64, iseg::Int
    if rh >= rhc                                    # upper (I_SEG=1)
        x = (1.0 - R) / (1.0 - Float64(rhc))
        y = x * (c2 + x * ((c1 / 2.0) - (c1 / 6.0) * x))
        dy_dx = c2 + x * (c1 - c1 / 2.0 * x)
        rh_length = 1.0 - Float64(rhc); iseg = 1
    elseif rh >= rhi2                               # middle (I_SEG=2)
        x = (R - Float64(rhi2)) / (Float64(rhc) - Float64(rhi2))
        if x > 0.0
            sus2 = (b1 * dlog10(x) <= -20.0) ? 0.0 : dpow(x, b1)
            y = b0 + x * (b4 + x * (-b2 / ((b1 + 1.0) * (b1 + 2.0)) * sus2 + b2 / 6.0 * x))
            sus3 = (b1 * dlog10(x) <= -20.0) ? 0.0 : dpow(x, b1 + 1.0)
            dy_dx = b4 - b2 / (b1 + 1.0) * sus3 + b2 / 2.0 * x * x
        else
            y = b0; dy_dx = b4
        end
        rh_length = Float64(rhc) - Float64(rhi2); iseg = 2
    elseif rhlongi > 0.0f0 && rh > rhi1             # straight (I_SEG=3)
        y = e1 + e2 * R
        dy_dx = e2; rh_length = 1.0; iseg = 3
    else                                            # lower (I_SEG=4)
        x = (Float64(rhi1) - R) / Float64(rhi1)
        y = a0 + x * ((a4 + a2 / a3) + x * (a2 / (2.0 * a3 * a3) + a1 * x)) + a2 * dlog(1.0 - x / a3)
        dy_dx = a4 + a2 / a3 + a2 / (a3 * a3) * x + 3.0 * a1 * x * x - a2 / (a3 - x)
        rh_length = Float64(rhi1); iseg = 4
    end
    dd_dh = dy_dx * Float64(f) / (rh_length * Float64(totalh))
    (iseg == 1 || iseg == 4) && (dd_dh = -dd_dh)
    return (Float32(Float64(f) * y), Float32(dd_dh))
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
    # f_other.f BRK_OT: BK is REAL*8 but DATA-initialized from default-REAL literals (Float32-rounded, then widened);
    # DR=DOB/DBHOB is a REAL quotient stored to REAL*8; PY is REAL*8; DBT=PY*DBTBH is stored to REAL (DBT) before
    # the REAL DIB=DOB-DBT. Each of those roundings is load-bearing (1-ULP Ht2TD on R2 aspen 200FW2W746).
    bk = _FW2_BK[jsp - 21]
    w(x) = Float64(Float32(x))
    b2 = w(bk[1]); b3 = w(bk[2]); b5 = w(bk[3]); c1 = w(bk[4]); c2 = w(bk[5]); c3 = w(bk[6])
    dr = dob > 0f0 ? Float64(dob / dbhob) : 0.0
    DBT = Float64(dbtbh)
    local py::Float64
    if ht2 > 4.5f0
        py = dr > 0.01 ? (dr * ((b2 - 1.0) / (b2 - dpow(dr, b3))) - ((dpow(dr, b5) - 1.0) / DBT)) : 0.0
    elseif ht2 == 4.5f0
        py = 1.0
    else
        clx = c1 * (dr - 1.0)
        py = clx >= 0.0 ? 1.0 + dpow(clx, c2 + c3 * DBT) : 1.0
    end
    dbt = Float32(py * DBT)
    dib = (dob > 0f0 ? dob : 0f0) - dbt
    return dib < 0f0 ? 0f0 : dib
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

"""
SF_YHAT with sf_yhat.f's own precision (JSP≠22), for the SF_HS solver: X is formed in REAL (single operands)
before widening to REAL*8; the REAL sub-terms (c1/2, c1/6, b2/((b1+1)(b1+2)), a4+a2/a3, a2/(2·a3²), 3·a1,
the straight segment e1+e2·rh) stay single; Y is REAL*8; DY_DX and the slope DD_DH are REAL. Returns
(dib, slope); slope is 0 unless `needsl`. `_fw2_sf_yhat` is this kernel (GETDIB, TCUBIC); `_fw2_sf_yhat_sl` keeps an all-Float64 kernel
(unused); SF_HS (`_fw2_sf_hs`) uses this one with the slope.
"""
function _fw2_sf_yhat_f(rh::Float32, tapcoe, rhfw, rflw, f::Float32, totalh::Float32, needsl::Bool)
    rh > 1f0 && return (0f0, 0f0)
    rh < 0f0 && return (f, 0f0)
    a3 = rflw[6]; rhi1 = rhfw[1]; rhi2 = rhfw[2]; rhc = rhfw[3]; rhlongi = rhfw[4]
    a0, a1, a2, a4, b0, b1, b2, b4, c1, c2, e1, e2 = tapcoe
    local y::Float64
    dydx = 0f0; rhl = 0f0; iseg = 0
    if rh >= rhc                                           # upper (I_SEG=1)
        x = Float64((1f0 - rh) / (1f0 - rhc))
        y = x * (Float64(c2) + x * (Float64(c1 / 2f0) - Float64(c1 / 6f0) * x))
        rhl = 1f0 - rhc; iseg = 1
        needsl && (dydx = Float32(Float64(c2) + x * (Float64(c1) - Float64(c1 / 2f0) * x)))
    elseif rh >= rhi2                                      # middle (I_SEG=2)
        x = Float64((rh - rhi2) / (rhc - rhi2))
        if x > 0.0
            sus2 = Float64(b1) * dlog10(x) <= -20.0 ? 0.0 : dpow(x, Float64(b1))
            y = Float64(b0) + x * (Float64(b4) + x * (-Float64(b2 / ((b1 + 1f0) * (b1 + 2f0))) * sus2 +
                                                     Float64(b2) / 6.0 * x))
        else
            y = Float64(b0)
        end
        if needsl
            rhl = rhc - rhi2; iseg = 2
            if x > 0.0
                sus3 = Float64(b1) * dlog10(x) <= -20.0 ? 0.0 : dpow(x, Float64(b1 + 1f0))
                dydx = Float32(Float64(b4) - Float64(b2) / (Float64(b1) + 1.0) * sus3 + Float64(b2) / 2.0 * x * x)
            else
                dydx = b4
            end
        end
    elseif rhlongi > 0f0 && rh > rhi1                      # straight (I_SEG=3)
        y = Float64(e1 + e2 * rh)
        needsl && (iseg = 3; dydx = e2; rhl = 1f0)
    else                                                   # lower (I_SEG=4)
        x = Float64((rhi1 - rh) / rhi1)
        y = Float64(a0) + x * (Float64(a4 + a2 / a3) + x * (Float64(a2 / (2f0 * a3 * a3)) + Float64(a1) * x)) +
            Float64(a2) * dlog(1.0 - x / Float64(a3))
        if needsl
            rhl = rhi1; iseg = 4
            dydx = Float32(Float64(a4 + a2 / a3) + Float64(a2 / (a3 * a3)) * x + Float64(3f0 * a1) * x * x -
                           Float64(a2) / (Float64(a3) - x))
        end
    end
    d = Float32(Float64(f) * y)
    needsl || return (d, 0f0)
    dd = dydx * f / (rhl * totalh)
    (iseg != 2 && iseg != 3) && (dd = -dd)
    return (d, dd)
end

"""
SF_HS (sf_hs.f) — the height at which the SF_YHAT profile reaches inside-bark diameter `dib`, for the
no-BRK_UP families (JSP∉22:30, NEXTRA=0 ⇒ SF_DS==SF_YHAT), e.g. the INGY I00/I13 FW2 equations.
Newton on the profile slope from a segment-based start; converged when |ADJUST|≤TOL·TOTALH AND |ERR|≤EPSILON;
re-started (IBREAK) on a positive slope; bisection fallback after 30 iterations. REAL*4 throughout, and the
Fortran's quirks kept (TOOLOW/TOOHIGH are assigned diameters in two places). Its answer can sit on the other
side of a NUMLOG/SEGMNT even-foot boundary from a diameter-tolerance bisection (BM bmt01 DF 12.7"×67':
bisection HS 53.99994 ⇒ 2-ft top log; SF_HS ⇒ 4-ft top log, MCF 21.1 / BdFt 107 as live).
"""
function _fw2_sf_hs(tapcoe, rhfw, rflw, f::Float32, totalh::Float32, dib::Float32; brk = nothing)::Float32
    epsilon = 0.001f0; tol = 0.0005f0
    # SF_DS (NEXTRA=0) = SF_YHAT; above the tip DIB=0, SLOPE=-1 (sf_ds.f). INEEDSL=0 calls return no slope.
    ds(h::Float32, sl::Bool = true) = h > totalh ? (0f0, -1f0) : _fw2_sf_yhat_f(h / totalh, tapcoe, rhfw, rflw, f, totalh, sl)
    # JSP 22-30 (region-2/3 outside-bark profiles): after every SF_DS, sf_hs.f calls BRK_UP(…,HTUP,D,DOB,DBT), which
    # (brk_up.f) sets DOB=D and overwrites D with BRK_OT's INSIDE-bark diameter at HTUP. `brk(htup, d)` does that
    # (identity for the INGY families). The HTUP sf_hs.f passes is kept verbatim: HI2/HI1/0.0 for the start guess,
    # HI2 inside the Newton loop (not the current H — FVS's own argument), HHIGH/HLOW/HTRY in the bisection.
    bu(htup::Float32, d::Float32) = brk === nothing ? d : brk(htup, d)
    rhi1 = rhfw[1]; rhi2 = rhfw[2]; rhlongi = rhfw[4]
    hi2 = rhi2 * totalh
    toohigh = totalh; toolow = 0f0
    di2 = bu(hi2, ds(hi2, false)[1])
    local rh::Float32
    if dib > di2
        toohigh = hi2
        start_found = false
        local di1::Float32
        if rhlongi > 0f0
            hi1 = rhi1 * totalh
            di1 = bu(hi1, ds(hi1, false)[1])
            if dib < di1
                toolow = di1
                rh = rhi2 - (rhi2 - rhi1) * (dib - di2) / (di1 - di2)
                start_found = true
            else
                toohigh = di1
            end
        else
            di1 = di2
        end
        if !start_found
            dbase = bu(0f0, ds(0f0, false)[1])
            dbase <= dib && return 0f0
            rz = fpow((dib - di1) / (dbase - di1), 0.25f0)
            rh = (1f0 - rz) * rhi1
        end
    else
        toolow = di2
        rz = 1f0 - fpow(dib / di2, 2f0)
        rh = rhi2 + (1f0 - rhi2) * rz
    end
    h = rh * totalh
    ibreak = 0
    while true                                         # label 110
        outcome = :exhausted                           # ITER > 30 ⇒ label 200
        for _ in 1:30                                  # label 20, ITER ≤ 30
            d, slope = ds(h)
            d = bu(hi2, d)                             # sf_hs.f: BRK_UP(…,HI2,D,…) — HI2, not H
            err = d - dib
            err < 0f0 && (toohigh = h)
            adjust = -err / slope
            h = h + adjust
            h > totalh && (h = (h - adjust + totalh) / 2f0)
            h < 0f0 && (h = (h - adjust) / 2f0)
            if !(abs(adjust) > tol * totalh || abs(err) > epsilon)
                if slope > 0f0 && ibreak < 2
                    ibreak += 1
                    h = ibreak == 1 ? 0.8f0 * h : h + (totalh - h) * 0.25f0
                    outcome = :restart
                else
                    outcome = :converged
                end
                break
            end
        end
        outcome === :converged && return h
        outcome === :restart && continue
        break
    end
    hhigh = toohigh; hlow = toolow                     # label 200: bisection fallback
    ehigh = bu(hhigh, ds(hhigh, false)[1]) - dib; elow = bu(hlow, ds(hlow, false)[1]) - dib
    elow * ehigh > 0f0 && return 0f0
    eps = epsilon * 2f0; iter = 0
    while true
        htry = (hhigh + hlow) / 2f0
        err = bu(htry, ds(htry, false)[1]) - dib
        abs(err) < eps && return htry
        iter > 40 && return htry
        iter += 1
        err > 0f0 ? (hlow = htry) : (hhigh = htry)
    end
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
function _fw2_merch_cuft(dibat, h::Float32, mtop::Float32, stump::Float32, minlen::Float32, merchl::Float32;
                         opt::Int = _NVB_R3_OPT, hs_solver = nothing, logs = nothing)::Float32
    hs = hs_solver === nothing ? _fw2_hs(dibat, mtop, h) : hs_solver(mtop)
    lmerch = hs - stump
    lmerch < merchl && return 0f0
    numseg = _nvb_numlog(opt, _NVB_R3_EVOD, lmerch, _NVB_R3_MAXLEN, minlen, _NVB_R3_TRIM)
    numseg == 0 && return 0f0
    loglen, numseg = _nvb_segmnt(opt, _NVB_R3_EVOD, lmerch, _NVB_R3_MAXLEN, minlen, _NVB_R3_TRIM, numseg)
    dibl = _fw2_dclass(dibat(4.5f0))           # butt log large end = DIB class at breast height
    ht2 = stump; vol4 = 0f0
    @inbounds for i in 1:numseg
        ht2 += _NVB_R3_TRIM + loglen[i]
        dib = dibat(ht2)
        (i == numseg && dib < mtop) && (dib = mtop)   # GETDIB forces the top log ≥ MTOPP
        dibs = _fw2_dclass(dib)
        logv = 0.00272708f0 * (dibl * dibl + dibs * dibs) * loglen[i]
        lv4 = floor(logv * 10f0 + 0.5f0) / 10f0
        vol4 += lv4
        logs === nothing || push!(logs, (dibs, lv4))      # ECVOL: (LOGDIA(I+1,1), LOGVOL(4,I))
        dibl = dibs
    end
    return vol4
end

"FW2 Scribner board VOL(2) (profile.f BFPFLG loop): buck stump→board-top (BFTOPD·BARK) and sum
per-log SCRIB. Same region-3 log-bucking as the cubic. `cor` = MRULES Scribner flag: 'Y' (region 2/3)
⇒ decimal-C ×10; 'N' (region 6/EC) ⇒ ANINT(raw board feet) per log (profile.f:441-446)."
function _fw2_board(dibat, h::Float32, bftop::Float32, stump::Float32, minlen::Float32, merchl::Float32;
                    cor::Char = 'Y', opt::Int = _NVB_R3_OPT, hs_solver = nothing, logs = nothing)::Float32
    hs = hs_solver === nothing ? _fw2_hs(dibat, bftop, h) : hs_solver(bftop)
    lmerch = hs - stump
    lmerch < merchl && return 0f0
    numseg = _nvb_numlog(opt, _NVB_R3_EVOD, lmerch, _NVB_R3_MAXLEN, minlen, _NVB_R3_TRIM)
    numseg == 0 && return 0f0
    loglen, numseg = _nvb_segmnt(opt, _NVB_R3_EVOD, lmerch, _NVB_R3_MAXLEN, minlen, _NVB_R3_TRIM, numseg)
    ht2 = stump; vol2 = 0f0
    @inbounds for i in 1:numseg
        ht2 += _NVB_R3_TRIM + loglen[i]
        dib = dibat(ht2)
        (i == numseg && dib < bftop) && (dib = bftop)
        dibs = _fw2_dclass(dib)                    # small-end inch class
        logv = _scrib(dibs, loglen[i], cor)
        lv1 = cor == 'Y' ? logv * 10f0 : _nint(logv)       # 'Y' ⇒ ×10 decimal-C; 'N' ⇒ ANINT(raw bdft)
        vol2 += lv1
        logs === nothing || push!(logs, (dibs, lv1))      # ECVOL: (LOGDIA(I+1,1), LOGVOL(1,I))
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
                    bark::Float32 = 1f0, topd::Float32 = 4f0, stump::Float32 = 1f0,
                    bftopd::Float32 = 6f0, iregn::Int = 3, board_cor::Char = 'Y',
                    merch_opt::Int = _NVB_R3_OPT, sf_hs::Bool = false,
                    log_bf = nothing, log_ft3 = nothing,   # optional ECVOL per-log capture (ECON units 4/5)
                    ht2td = nothing)   # optional 2-slot buffer ← [HT1PRD of the cubic call, HT1PRD of the board call]
    vol = zeros(Float32, 15)
    ht2td === nothing || (ht2td[1] = 0f0; ht2td[2] = 0f0)   # profile.f:117 early exits leave HT1PRD=0
    (d < 1f0 || h < 5f0) && return vol   # profile.f:117 HTTOT.LT.5 (strict; h==5.0 IS computed)
    jsp = _fw2_jsp(voleq)
    (_fw2_is_ingy(jsp) || (22 <= jsp <= 29)) || return vol   # supported 2-pt families (22 = Black Hills PP)
    ingy = _fw2_is_ingy(jsp)
    geosub = length(voleq) >= 3 ? voleq[2:3] : "00"    # SHP_C2 subregion (INGY only)
    rflw, rhfw = _fw2_shp(jsp, geosub, d, h)
    tapcoe = _fw2_sf_taper(rhfw, rflw)
    # INGY: profile is inside bark, calibrated to DBHIB. SF_SHP uses the PASSED DBTBH (fvsvol DBTBH=D·(1-BARK))
    # when >0 — so DBHIB=D·BARK (cr_bratio), NOT FDBT_C2 (that's only the no-bark-input fallback). Region-2/3: DBHOB.
    # sf_shp.f (JSP 11-21): DBHIB = DBHOB - DBTBH with the passed DBTBH = D·(1-BARK) (fvsvol.f) — formed as that
    # REAL*4 subtraction, not D·BARK (equal in exact arithmetic, 1 ULP apart in Float32 on many trees).
    dbtbh = d * (1f0 - bark)
    dbhib = ingy ? d - dbtbh : d
    # sf_2pt.f:63-66: F = DBH_IB / SF_YHAT(4.5/TOTALH) through sf_yhat.f itself (REAL*4 X, REAL*8 Y — the faithful
    # kernel _fw2_sf_yhat_f), which the SF_HS merch-top solve also uses.
    yhat_bh = _fw2_sf_yhat_f(4.5f0 / h, tapcoe, rhfw, rflw, 1.0f0, h, false)[1]
    yhat_bh == 0f0 && return vol
    f = dbhib / yhat_bh
    dibat = ingy ? (ht -> _fw2_sf_yhat(ht / h, tapcoe, rhfw, rflw, f)) :
                   (ht -> _fw2_brk_ot(jsp, d, _fw2_sf_yhat(ht / h, tapcoe, rhfw, rflw, f), ht, dbtbh))
    # Small trees (HTTOT≤15): FWSMALL corrects the stump diameter; merch/board stay 0 (LMERCH<MERCHL).
    stump_dib = h <= 15f0 ? _fw2_fwsmall(jsp, h, dibat(1.0f0), d * bark) : -1f0
    minl = _cr_merch_minlen(iregn); merl = _cr_merch_merchl(iregn)
    vol[1] = _nint(_fw2_tcubic(dibat, h; stump_dib = stump_dib) * 10.0f0) * 1f-1     # profile.f:293 VOL(1)=NINT(TCVOL*10.0)*1E-1 — MULTIPLY by REAL 0.1 (≠ /10: e.g. 1334·0.1f0=133.40000915 vs 133.4f0=133.39999390)
    # `sf_hs=true`: MERLEN's merch-top height from the faithful SF_HS Newton (INGY, and JSP 22-30 with BRK_UP);
    # otherwise the legacy diameter-tolerance bisection (kept for the callers not yet re-validated on it).
    hs_solver = !sf_hs ? nothing :
                ingy ? (top -> _fw2_sf_hs(tapcoe, rhfw, rflw, f, h, top)) :
                (top -> _fw2_sf_hs(tapcoe, rhfw, rflw, f, h, top;
                                   brk = (hu, dd) -> _fw2_brk_ot(jsp, d, dd, hu, dbtbh)))   # JSP 22-30 BRK_UP
    # HT1PRD (fvsvol.f:338/485 → HT2TD): profile.f:335-342 MERLEN, VOLEQ(4:4)='F' ⇒ LMERCH = HS−STUMP (sf_hs.f at
    # DS=TOP), clamped ≥0 (profile.f:1052), then HT1PRD = LMERCH+STUMP — both Float32. Same HS as the volume uses.
    if ht2td !== nothing
        _ht1prd(top) = (hs = hs_solver === nothing ? _fw2_hs(dibat, top, h) : hs_solver(top);
                        lm = hs - stump; lm < 0f0 && (lm = 0f0); lm + stump)
        ht2td[1] = _ht1prd(topd * bark); ht2td[2] = _ht1prd(bftopd * bark)
    end
    vol[4] = _fw2_merch_cuft(dibat, h, topd * bark, stump, minl, merl; opt = merch_opt, hs_solver = hs_solver,
                             logs = log_ft3)
    vol[2] = _fw2_board(dibat, h, bftopd * bark, stump, minl, merl; cor = board_cor, opt = merch_opt,
                        hs_solver = hs_solver, logs = log_bf)
    return vol
end
