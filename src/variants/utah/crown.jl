# =============================================================================
# crown.jl (utah) — per-tree CCF (ut/ccfcal.f MODE=1). Chunk 5 (partial: the CCF path
# that feeds dgf!/htgf RELDEN + PCCF; crown-ratio dub + crown-width come later).
#
# UT CCF is the direct per-species polynomial (NI INT-133 / Paine-Hann form), summed over
# the tree list = stand CCF (RELDEN) + per-point PCCF (like EM/KT/IE/TT, not the CR crown-
# width→area path):
#   large (D≥BREAK):  CCFT = RD1 + D·RD2 + D²·RD3
#   small (D>0.1):    sp 20/21 (MC/BI): D·(RD1+RD2+RD3);  else RDA·D^RDB
#   tiny  (D≤0.1):    sp 20/21: D·(RD1+RD2+RD3);          else 0.001
#   BREAK = 10 for sp 17:19,22 (GB/NC/FC/BE); else 1.
# Coefficients verbatim from ut/ccfcal.f DATA (species 1=WB..24=OH).
# =============================================================================

const UT_RD1 = Float32[0.01925, 0.01925, 0.11, 0.04, 0.03, 0.03, 0.01925, 0.03, 0.03, 0.03,
                       0.01925, 0.01925, 0.03, 0.01925, 0.01925, 0.01925, 0.01925, 0.03, 0.03, 0.0204,
                       0.0204, 0.03, 0.01925, 0.03]
const UT_RD2 = Float32[0.01676, 0.01676, 0.0333, 0.0270, 0.0173, 0.0238, 0.01676, 0.0173, 0.0216, 0.0180,
                       0.01676, 0.01676, 0.0215, 0.01676, 0.01676, 0.01676, 0.01676, 0.0215, 0.0215, 0.0246,
                       0.0246, 0.0215, 0.01676, 0.0215]
const UT_RD3 = Float32[0.00365, 0.00365, 0.00259, 0.00405, 0.00259, 0.00490, 0.00365, 0.00259, 0.00405, 0.00281,
                       0.00365, 0.00365, 0.00363, 0.00365, 0.00365, 0.00365, 0.00365, 0.00363, 0.00363, 0.0074,
                       0.0074, 0.00363, 0.00365, 0.00363]
const UT_RDA = Float32[0.009187, 0.009187, 0.017299, 0.015248, 0.007875, 0.008915, 0.009187, 0.007875, 0.011402, 0.007813,
                       0.009187, 0.009187, 0.011109, 0.009187, 0.009187, 0.009187, 0.009187, 0.011109, 0.011109, 0.0,
                       0.0, 0.011109, 0.009187, 0.011109]
const UT_RDB = Float32[1.7600, 1.7600, 1.5571, 1.7333, 1.7360, 1.7800, 1.7600, 1.7360, 1.7560, 1.7680,
                       1.7600, 1.7600, 1.7250, 1.7600, 1.7600, 1.7600, 1.7600, 1.7250, 1.7250, 0.0,
                       0.0, 1.7250, 1.7600, 1.7250]

# ut/ccfcal.f MODE=1: per-tree CCF (before the ×P expansion the caller applies).
# ut/crown.f WEIBULL crown-ratio change model (= TT/CR form). Coefficients verbatim ut/crown.f DATA.
const UT_WEIBA  = Float32[1,1,1,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0]
const UT_WEIBB0 = Float32[-0.82631,-0.82631,-0.24217,-0.89553,-0.90648,-0.08414,0.17162,-0.90648,-0.89553,-0.82631,0,0,0,0,0,0,0,0,0,-0.23830,-0.23830,0,-0.26595,0]
const UT_WEIBB1 = Float32[1.06217,1.06217,0.96529,1.07728,1.08122,1.14765,1.07338,1.08122,1.07728,1.06217,0,0,0,0,0,0,0,0,0,1.18016,1.18016,0,0.98326,0]
const UT_WEIBC0 = Float32[3.31429,3.31429,-7.94832,1.74621,3.48889,2.775,3.15,3.48889,1.74621,-1.02873,0,0,0,0,0,0,0,0,0,3.04,3.04,0,-1.60411,0]
const UT_WEIBC1 = Float32[0,0,1.93832,0.29052,0,0,0,0,0.29052,0.80143,0,0,0,0,0,0,0,0,0,0,0,0,1.60411,0]
const UT_CRC0   = Float32[6.19911,6.19911,7.46296,7.65751,6.81087,4.01678,6.00567,6.81087,7.65751,6.19911,0,0,0,0,0,0,0,0,0,4.62512,4.62512,0,7.92810,0]
const UT_CRC1   = Float32[-0.02216,-0.02216,-0.02944,-0.03513,-0.01037,-0.01516,-0.0352,-0.01037,-0.03513,-0.02216,0,0,0,0,0,0,0,0,0,-0.01604,-0.01604,0,-0.06298,0]

# Weibull species (nonzero coeffs) = 1-10,20,21,23; PJ/hardwoods (11-19,22,24) get crown from REGENT (skip).
@inline _ut_crown_weibull(sp::Int) = sp <= 10 || sp == 20 || sp == 21 || sp == 23

# ut/crown.f crown-ratio change (Weibull) — mirror TT crown_ratio_update! with UT coeffs + sp5/sp23 RELSDI cap.
function crown_ratio_update!(s::StandState, ::Utah; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    relden = p.relative_density; sdiac = crown_sdi
    key = Vector{Float32}(undef, n); idx = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        bk = bark_ratio(s.calib.bark_a, s.calib.bark_b, Int(t.species[i]), t.dbh[i])
        key[i] = t.dbh[i] + t.diam_growth[i] / bk; idx[i] = Int32(i)
    end
    _rdpsrt!(key, idx; lseq = false)
    isort = Vector{Int32}(undef, n)
    @inbounds for jj in 1:n; isort[idx[jj]] = Int32(n - jj + 1); end
    @inbounds for i in 1:n
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (lstart && t.crown_pct[i] > 0) && continue
        _ut_crown_weibull(sp) || continue                  # PJ/hardwoods → crown from REGENT
        icr = Int(t.crown_pct[i])
        relsdi = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0
        relsdi > 1.5f0 && (relsdi = 1.5f0)
        ((sp == 5 || sp == 23) && relsdi > 1f0) && (relsdi = 1f0)   # BS/OS cap (ut/crown.f:184)
        acrnew = UT_CRC0[sp] + UT_CRC1[sp] * relsdi * 100f0
        A = UT_WEIBA[sp]
        B = UT_WEIBB0[sp] + UT_WEIBB1[sp] * acrnew; B < 1f0 && (B = 1f0)
        C = UT_WEIBC0[sp] + UT_WEIBC1[sp] * acrnew; C < 2f0 && (C = 2f0)
        scale = 1f0 - 0.00167f0 * (relden - 100f0)
        scale > 1f0 && (scale = 1f0); scale < 0.30f0 && (scale = 0.30f0)
        x = d > 0f0 ? (Float32(isort[i]) / Float32(n)) * scale : 0.5f0 * scale
        x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
        crnew = (A + B * (-log(1f0 - x))^(1f0 / C)) * 10f0
        if !(lstart || icr == 0)
            chg = crnew - Float32(icr); pdifpy = chg / Float32(icr) / fint
            pdifpy > 0.01f0 && (chg = Float32(icr) * 0.01f0 * fint)
            pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
            crnew = Float32(icr) + chg
        end
        icri = trunc(Int, crnew + 0.5f0)
        if !(lstart || icr == 0)
            htg = t.ht_growth[i]; crln = h * Float32(icr) / 100f0
            crmax = (crln + htg) / (h + htg) * 100f0
            (icri < 10 && Float32(icri) <= crmax) && (icri = trunc(Int, crmax + 0.5f0))
            Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
        end
        icri < 0 && (icri = 0); icri > 100 && (icri = 100)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end

@inline function ut_tree_ccf(sp::Integer, d::Real)::Float32
    d <= 0f0 && return 0f0
    dd = Float32(d)
    brk = (17 <= sp <= 19 || sp == 22) ? 10f0 : 1f0
    if dd >= brk
        return UT_RD1[sp] + dd * UT_RD2[sp] + dd * dd * UT_RD3[sp]
    elseif dd > 0.1f0
        return (sp == 20 || sp == 21) ? dd * (UT_RD1[sp] + UT_RD2[sp] + UT_RD3[sp]) : UT_RDA[sp] * dd ^ UT_RDB[sp]
    else
        return (sp == 20 || sp == 21) ? dd * (UT_RD1[sp] + UT_RD2[sp] + UT_RD3[sp]) : 0.001f0
    end
end
