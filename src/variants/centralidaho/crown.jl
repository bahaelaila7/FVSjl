# =============================================================================
# crown.jl (centralidaho) — CI per-tree CCF (ci/ccfcal.f MODE=1). Chunk 5 (CCF portion).
# The stand CCF = Σ CCFT·TPA is the RELDEN the DG (dgf!) CONSPP term reads. Same western
# polynomial as KT/IE (Bush/Crookston): D≥10 → RD1+D·RD2+D²·RD3; D<10 → RDA·D^RDB. CI's RD*
# coefficients (ci/ccfcal.f DATA) match KT for the shared N-Rockies conifers.
# (Crown-WIDTH MODE=2 B1..B6 + crown-ratio model land separately with crown_ratio_update!.)
# =============================================================================

const CI_RD1 = Float32[0.03,0.02,0.11,0.04,0.03,0.03,0.01925,0.03,0.03,0.03,0.01925,0.01925,0.03,0.01925,0.0204,0.01925,0.03,0.03,0.03]
const CI_RD2 = Float32[0.0167,0.0148,0.0333,0.027,0.0215,0.0238,0.01676,0.0173,0.0216,0.018,0.01676,0.01676,0.0238,0.01676,0.0246,0.01676,0.0215,0.0215,0.0215]
const CI_RD3 = Float32[0.0023,0.00338,0.00259,0.00405,0.00363,0.0049,0.00365,0.00259,0.00405,0.00281,0.00365,0.00365,0.0049,0.00365,0.0074,0.00365,0.00363,0.00363,0.00363]
const CI_RDA = Float32[0.009884,0.007244,0.017299,0.015248,0.011109,0.008915,0.009187,0.007875,0.011402,0.007813,0.009187,0.009187,0.008915,0.009187,0.0,0.009187,0.011109,0.011109,0.011109]
const CI_RDB = Float32[1.6667,1.8182,1.5571,1.7333,1.725,1.78,1.76,1.736,1.756,1.768,1.76,1.76,1.78,1.76,0.0,1.76,1.725,1.725,1.725]

"""Per-tree CCF contribution (ci/ccfcal.f MODE=1), excluding the ×TPA factor."""
@inline function ci_tree_ccf(sp::Integer, d::Real)::Float32
    d <= 0f0 && return 0f0
    if d >= 10f0
        return CI_RD1[sp] + d * CI_RD2[sp] + d * d * CI_RD3[sp]
    else
        return CI_RDA[sp] * Float32(d) ^ CI_RDB[sp]
    end
end

# ci/crown.f Weibull crown ratio (chunk 5b). ACRNEW=C0+C1·RELSDI·100; A=WEIBA; B=WEIBB0+WEIBB1·ACRNEW;
# C=WEIBC0 (WEIBC1=0); crnew=(A+B·(−ln(1−x))^(1/C))·10, x=ISORT/N·SCALE. CRNMLT=1/DLOW=0/DHI=99 (no-op mods).
const CI_WEIBA = Float32[2.0,0.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,0.0,1.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0]
const CI_WEIBB0 = Float32[-2.12713,0.07609,-1.19297,-1.19297,0.06593,-1.38636,0.07609,-0.91567,-0.91567,0.24916,-0.82631,-0.82631,-0.08414,0.0,-0.2383,-0.82631,0.0,0.07609,0.0]
const CI_WEIBB1 = Float32[1.10526,1.10184,1.12928,1.12928,1.09624,1.16801,1.10184,1.06469,1.06469,1.04831,1.06217,1.06217,1.14765,0.0,1.18016,1.06217,0.0,1.10184,0.0]
const CI_WEIBC0 = Float32[2.77,3.01,3.42,3.42,3.71,3.02,3.01,3.5,3.5,4.36,3.31429,3.31429,2.775,0.0,3.04,3.31429,0.0,3.01,0.0]
const CI_CRC0 = Float32[7.16846,5.50719,5.52653,5.52653,6.61291,6.17373,5.50719,6.774,6.12779,6.41166,6.19911,6.19911,4.01678,0.0,4.62512,6.19911,0.0,7.238,0.0]
const CI_CRC1 = Float32[-0.02375,-0.01833,0.0,0.0,-0.02182,-0.01795,-0.01833,0.0,-0.01269,-0.02041,-0.02216,-0.02216,-0.01516,0.0,-0.01604,-0.02216,0.0,0.0,0.0]

function crown_ratio_update!(s::StandState, ::CentralIdaho; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    relden = p.relative_density; sdiac = crown_sdi
    sd = s.coef.species
    key = Vector{Float32}(undef, n); idx = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        bk = ci_bratio(sd, Int(t.species[i]), t.dbh[i])
        key[i] = t.dbh[i] + t.diam_growth[i] / bk; idx[i] = Int32(i)
    end
    _rdpsrt!(key, idx; lseq = false)
    isort = Vector{Int32}(undef, n)
    @inbounds for jj in 1:n; isort[idx[jj]] = Int32(n - jj + 1); end
    @inbounds for i in 1:n
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (sp < 1 || sp > 19) && continue
        (lstart && t.crown_pct[i] > 0) && continue
        (d < 1f0 && lstart) && continue
        icr = Int(t.crown_pct[i])
        relsdi = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0
        relsdi > 1.5f0 && (relsdi = 1.5f0)
        acrnew = CI_CRC0[sp] + CI_CRC1[sp] * relsdi * 100f0
        A = CI_WEIBA[sp]
        B = CI_WEIBB0[sp] + CI_WEIBB1[sp] * acrnew; B < 1f0 && (B = 1f0)
        C = CI_WEIBC0[sp]; C < 2f0 && (C = 2f0)
        scale = 1f0 - 0.00167f0 * (relden - 100f0)
        scale > 1f0 && (scale = 1f0); scale < 0.30f0 && (scale = 0.30f0)
        x = d > 0f0 ? (Float32(isort[i]) / Float32(n)) * scale : 0.5f0 * scale
        x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
        crnew = (A + B * (-log(1f0 - x))^(1f0 / C)) * 10f0
        if !(lstart || icr == 0)
            chg = crnew - Float32(icr); pdifpy = chg / Float32(icr) / fint
            pdifpy > 0.01f0 && (chg = Float32(icr) * 0.01f0 * fint)
            pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
            crnew = Float32(icr) + chg                      # CRNMLT=1, DLOW=0/DHI=99 ⇒ no-op branch
        end
        icri = trunc(Int, crnew + 0.5f0)
        if !(lstart || icr == 0)
            crln = h * Float32(icr) / 100f0; htg = t.ht_growth[i]
            crmax = (h + htg) > 0f0 ? (crln + htg) / (h + htg) * 100f0 : 100f0
            Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
        end
        icri < 1 && (icri = 1); icri > 95 && (icri = 95)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end
