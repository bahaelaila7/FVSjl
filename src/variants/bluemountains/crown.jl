# =============================================================================
# crown.jl (bluemountains) — CCF / crown width (bm/ccfcal.f). Chunk 5 (partial: the CCF
# piece is a prerequisite for chunk-3 DG, since RELDEN = stand CCF feeds CONSPP).
#
# bm/ccfcal.f MODE=1 per-tree CCFT (Paine & Hann, Oregon RP46):
#   CASE(1:12,15,17):  D>=1  → RD1 + D·RD2 + D²·RD3 ; 0.1<D<1 → RDA·D^RDB ; D<=0.1 → 0.001
#   CASE(13,14,16,18): D<1   → D·(RD1+RD2+RD3)      ; D>=1    → RD1 + RD2·D + RD3·D²
# Stand CCF = Σ CCFT·TPA = RELDEN.  Coeffs from data/bluemountains/ccf_coeffs_bm.csv (verified vs source).
# =============================================================================

let
    path = joinpath(BM_DATADIR, "ccf_coeffs_bm.csv")
    rows = [split(strip(l), ',') for l in readlines(path)[2:end]]
    col(j) = Float32[parse(Float32, rows[sp][j]) for sp in 1:18]
    global const BM_RD1 = col(2); global const BM_RD2 = col(3); global const BM_RD3 = col(4)
    global const BM_RDA = col(5); global const BM_RDB = col(6)
end

@inline function bm_tree_ccf(sp::Integer, d::Real)::Float32
    d <= 0f0 && return 0f0
    dd = Float32(d)
    if sp == 13 || sp == 14 || sp == 16 || sp == 18
        return dd < 1f0 ? dd * (BM_RD1[sp] + BM_RD2[sp] + BM_RD3[sp]) :
                          BM_RD1[sp] + BM_RD2[sp] * dd + BM_RD3[sp] * dd * dd
    else
        return dd >= 1f0 ? BM_RD1[sp] + dd * BM_RD2[sp] + dd * dd * BM_RD3[sp] :
               dd > 0.1f0 ? BM_RDA[sp] * dd ^ BM_RDB[sp] : 0.001f0
    end
end

# --- bm/crown.f WEIBULL crown-ratio change model (= TT/UT/CR form). Coeffs from crown_coeffs_bm.csv. ---
let
    path = joinpath(BM_DATADIR, "crown_coeffs_bm.csv")
    rows = [split(strip(l), ',') for l in readlines(path)[2:end]]
    col(name) = Float32[parse(Float32, rows[sp][findfirst(==(name), split(strip(readlines(path)[1]), ','))]) for sp in 1:18]
    global const BM_WEIBA  = col("WEIBA");  global const BM_WEIBB0 = col("WEIBB0"); global const BM_WEIBB1 = col("WEIBB1")
    global const BM_WEIBC0 = col("WEIBC0"); global const BM_WEIBC1 = col("WEIBC1")
    global const BM_CRC0   = col("C0");     global const BM_CRC1   = col("C1")
    global const BM_CRNMLT = col("CRNMLT"); global const BM_CR_DLOW = col("DLOW"); global const BM_CR_DHI = col("DHI")
end

# bm/crown.f — Weibull crown-ratio (all species; small trees D<1 → REGENT). RELSDI=SDIAC/SDIDEF,
# ACRNEW=C0+C1·RELSDI·100, Weibull A/B/C (B<1→1, C<2→2), SCALE=1−0.00167·(RELDEN−100), rank-based X.
function crown_ratio_update!(s::StandState, ::BlueMountains; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    relden = p.relative_density; sdiac = crown_sdi
    key = Vector{Float32}(undef, n); idx = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        bk = bm_bratio(s.coef.species, Int(t.species[i]), t.dbh[i])
        key[i] = t.dbh[i] + t.diam_growth[i] / bk; idx[i] = Int32(i)
    end
    _rdpsrt!(key, idx; lseq = false)
    isort = Vector{Int32}(undef, n)
    @inbounds for jj in 1:n; isort[idx[jj]] = Int32(n - jj + 1); end
    @inbounds for i in 1:n
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (lstart && t.crown_pct[i] > 0) && continue
        (d < 1f0 && lstart) && continue                    # small trees → REGENT (bm/crown.f:237)
        icr = Int(t.crown_pct[i])
        relsdi = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0
        relsdi > 1.5f0 && (relsdi = 1.5f0)
        acrnew = BM_CRC0[sp] + BM_CRC1[sp] * relsdi * 100f0
        A = BM_WEIBA[sp]
        B = BM_WEIBB0[sp] + BM_WEIBB1[sp] * acrnew; B < 1f0 && (B = 1f0)
        C = BM_WEIBC0[sp] + BM_WEIBC1[sp] * acrnew; C < 2f0 && (C = 2f0)
        scale = 1f0 - 0.00167f0 * (relden - 100f0)
        scale > 1f0 && (scale = 1f0); scale < 0.30f0 && (scale = 0.30f0)
        x = d > 0f0 ? (Float32(isort[i]) / Float32(n)) * scale : 0.5f0 * scale
        x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
        crnew = (A + B * (-log(1f0 - x))^(1f0 / C)) * 10f0
        if !(lstart || icr == 0)
            chg = crnew - Float32(icr); pdifpy = chg / Float32(icr) / fint
            pdifpy > 0.01f0 && (chg = Float32(icr) * 0.01f0 * fint)
            pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
            crnew = (d >= BM_CR_DLOW[sp] && d <= BM_CR_DHI[sp]) ? Float32(icr) + chg * BM_CRNMLT[sp] :
                    Float32(icr) + chg
        end
        icri = trunc(Int, crnew + 0.5f0)
        if lstart || icr == 0
            (d >= BM_CR_DLOW[sp] && d <= BM_CR_DHI[sp]) && (icri = trunc(Int, Float32(icri) * BM_CRNMLT[sp]))
        else
            # CRMAX cap (bm/crown.f:301-314): crown length can't exceed the height-growth-adjusted max.
            crln = h * Float32(icr) / 100f0; htg = t.ht_growth[i]
            crmax = (crln + htg) / (h + htg) * 100f0
            Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
            (icri < 10 && BM_CRNMLT[sp] == 1f0) && (icri = trunc(Int, crmax + 0.5f0))
        end
        # final clamps (bm/crown.f:347-349)
        icri > 95 && (icri = 95)
        (icri < 10 && BM_CRNMLT[sp] == 1f0) && (icri = 10)
        icri < 1 && (icri = 1)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end
