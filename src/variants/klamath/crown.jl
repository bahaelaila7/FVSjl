# =============================================================================
# crown.jl (klamath) — NC crown ratio (nc/crown.f). Chunk 5.
# Weibull crown-ratio: ACRNEW=C0+C1·RELSDI·100; B=WEIBB0+WEIBB1·ACRNEW; C=WEIBC0+WEIBC1·ACRNEW; A=0.
# Per-tree X = (rank/N)·SCALE, SCALE=clamp(1.5−RELSDI,0.30,1.0); CR=(A+B·(−ln(1−X))^(1/C))·10.
# sp12 RW = logistic X = 1/(1+exp(−1.021064+0.309296·lnHDR+0.869720·PRD−0.116274·D/QMDPLT)); CR=X·100.
# 1%/yr change-limit vs old ICR (CRNMLT=1, DLOW=0/DHI=99 ⇒ no band effect). Floor [10,95] (redwood [5,95]).
# =============================================================================

const NC_WEIBB0 = Float32[0.52909,0.25115,0.52909,0.48464,0.08402,0.29964,0.06607,0.25667,0.16601,0.03685,0.25667,0.0]
const NC_WEIBB1 = Float32[1.00677,1.05987,1.00677,1.01272,1.10297,1.05398,1.10705,1.06474,1.08150,1.09499,1.06474,0.0]
const NC_WEIBC0 = Float32[-3.48211,0.33383,-3.48211,-2.78353,0.91078,-1.09270,2.04714,0.11729,0.91420,4.01340,0.11729,0.0]
const NC_WEIBC1 = Float32[1.38780,0.63833,1.38780,1.27283,0.45819,0.80687,0.15070,0.61681,0.45768,0.04946,0.61681,0.0]
const NC_CRC0 = Float32[7.48846,6.92893,7.48846,7.44422,3.64292,5.12357,6.82187,5.95912,6.14578,6.04928,5.95912,0.0]
const NC_CRC1 = Float32[-0.02899,-0.04053,-0.02899,-0.04779,-0.00317,-0.01042,-0.02247,-0.01812,-0.02781,-0.01091,-0.01812,0.0]

function crown_ratio_update!(s::StandState, ::Klamath; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    sd = s.coef.species
    ba = p.basal_area; sdiac = crown_sdi
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    qmd = stand_qmd(s)
    # rank by (D+G) descending → ISORT (largest = rank 1)
    key = Vector{Float32}(undef, n); idx = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        bk = bark_ratio(bark_a, bark_b, Int(t.species[i]), t.dbh[i])
        key[i] = t.dbh[i] + t.diam_growth[i] / bk; idx[i] = Int32(i)
    end
    _rdpsrt!(key, idx; lseq = false)
    isort = Vector{Int32}(undef, n)
    @inbounds for jj in 1:n; isort[idx[jj]] = Int32(n - jj + 1); end
    @inbounds for i in 1:n
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (sp < 1 || sp > 12) && continue
        (lstart && t.crown_pct[i] > 0) && continue
        (d < 1f0 && lstart) && continue
        icr = Int(t.crown_pct[i])
        relsdi = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0
        relsdi > 1.5f0 && (relsdi = 1.5f0)
        local crnew::Float32
        if sp == 12                                       # redwood — logistic
            hdr = d > 0f0 ? h / d : 1f0
            prd = relsdi
            qp = qmd > 0f0 ? d / qmd : 1f0
            xl = -1.021064f0 + 0.309296f0 * log(max(hdr, 1f-3)) + 0.869720f0 * prd - 0.116274f0 * qp
            x = 1f0 / (1f0 + exp(xl))
            x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
            crnew = x * 100f0
        else
            acrnew = NC_CRC0[sp] + NC_CRC1[sp] * relsdi * 100f0
            b = NC_WEIBB0[sp] + NC_WEIBB1[sp] * acrnew
            c = NC_WEIBC0[sp] + NC_WEIBC1[sp] * acrnew
            scale = 1.5f0 - relsdi
            scale > 1f0 && (scale = 1f0); scale < 0.30f0 && (scale = 0.30f0)
            x = d > 0f0 ? (Float32(isort[i]) / Float32(n)) * scale : 0.5f0 * scale
            x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
            crnew = (b * (-log(1f0 - x))^(1f0 / c)) * 10f0    # A=WEIBA=0
        end
        if !(lstart || icr == 0)
            chg = crnew - Float32(icr); pdifpy = chg / Float32(icr) / fint
            pdifpy > 0.01f0 && (chg = Float32(icr) * 0.01f0 * fint)
            pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
            crnew = Float32(icr) + chg                      # CRNMLT=1 ⇒ no band multiplier
        end
        icri = trunc(Int, crnew + 0.5f0)
        lo = sp == 12 ? 5 : 10
        icri > 95 && (icri = 95); icri < lo && (icri = lo)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end
