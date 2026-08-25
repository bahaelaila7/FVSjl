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

# --- bm/dubscr.f DUBSCR — crown-ratio dub for read-inventory D<1 (missing CR) + regen inserts. ---
# Logistic (sp 1-12,15,17) / linear-rescale (sp 13,14,16,18) crown model. Coeffs verified vs bm/dubscr.f DATA.
const BM_BCR0  = Float32[-1.669490,-1.669490,-.426688,-.426688,-.426688,-2.19723,-1.669490,-.426688,-.426688,-1.669490,-1.66949,-1.66949,6.489813,7.558538,-.426688,5.0,-1.669490,5.0]
const BM_BCR1  = Float32[-.209765,-.209765,-.093105,-.093105,-.093105,0.0,-.209765,-.093105,-.093105,-.209765,-.209765,-.209765,0.0,0.0,-.093105,0.0,-.209765,0.0]
const BM_BCR2  = Float32[0.0,0.0,.022409,.022409,.022409,0.0,0.0,.022409,.022409,0.0,0.0,0.0,-0.029815,-0.015637,.022409,0.0,0.0,0.0]
const BM_BCR3  = Float32[.003359,.003359,.002633,.002633,.002633,0.0,.003359,.002633,.002633,.003359,.003359,.003359,-0.009276,-0.009064,.002633,0.0,.003359,0.0]
const BM_BCR5  = Float32[.011032,.011032,0.0,0.0,0.0,0.0,.011032,0.0,0.0,.011032,.011032,.011032,0.0,0.0,0.0,0.0,.011032,0.0]
const BM_BCR6  = Float32[0.0,0.0,-.045532,-.045532,-.045532,0.0,0.0,-.045532,-.045532,0.0,0.0,0.0,0.0,0.0,-.045532,0.0,0.0,0.0]
const BM_BCR8  = Float32[.017727,.017727,0.0,0.0,0.0,0.0,.017727,0.0,0.0,.017727,.017727,.017727,0.0,0.0,0.0,0.0,.017727,0.0]
const BM_BCR9  = Float32[-.000053,-.000053,.000022,.000022,.000022,0.0,-.000053,.000022,.000022,-.000053,-.000053,-.000053,0.0,0.0,.000022,0.0,-.000053,0.0]
const BM_BCR10 = Float32[.014098,.014098,-.013115,-.013115,-.013115,0.0,.014098,-.013115,-.013115,.014098,.014098,.014098,0.0,0.0,-.013115,0.0,.014098,0.0]
const BM_DUBSD = Float32[.5000,.5000,.6957,.6957,.6957,0.200,.6124,.6957,.6957,.4942,.5000,.5000,2.0426,1.9658,.9310,0.500,.4942,0.500]

# so/adjmai.f ADJMAI — adjusted MAI site-index polynomial (10 eq-groups keyed by FIA site-species code).
# maical.f: ISPNUM maps BM species idx 1-18 → FIA code for the site-index lookup; RMAI capped at 128.
const _BM_ISPNUM = Int32[119,117,202,15,264,101,108,93,21,122,101,101,101,42,746,746,122,746]
const _BM_ADJ_ISP = Int32[19,41,42,71,81,94,95,242,263,264,98,298, 101,103,108,109,120,124, 119, 116,117,122, 201,202, 73, 92,93, 11,15,17,19,20,21,22, 211,212]
const _BM_ADJ_IMAP = Int32[fill(Int32(1),12); fill(Int32(2),6); Int32(3); fill(Int32(4),3); fill(Int32(5),2); Int32(6); fill(Int32(7),2); fill(Int32(9),7); fill(Int32(10),2)]

@inline function _bm_adjmai_group(grp::Integer, si::Float32)::Float32
    a = 0f0
    if grp == 1;      si < 33f0 && return 0f0; a = -63.689706f0 + 1.9402941f0*si
    elseif grp == 2;  si < 11f0 && return 0f0; a = -12.0388f0 + 1.18672f0*si
    elseif grp == 3;  a = 5.972615f0 + 1.857675f0*si
    elseif grp == 4;  a = 2.305357f0 + 0.033890056f0*si + 0.0090108543f0*si*si
    elseif grp == 5;  si < 29f0 && return 0f0; a = -10.303313f0 + .032929911f0*si + .012207163f0*si*si - .00003543129f0*si*si*si
    elseif grp == 6;  si < 11f0 && return 0f0; a = -6.0892857f0 + .45178571f0*si + .014464286f0*si*si
    elseif grp == 7;  si < 10f0 && return 0f0; a = -18.4f0 + 1.92f0*si
    elseif grp == 8;  si < 32f0 && return 0f0; a = -53.892857f0 + 1.7178571f0*si
    elseif grp == 9;  a = -4.89001f0 + 311.29546f0*((exp((si/170f0-1f0)^3/0.343f0)-0.055f0)/0.95f0)
    elseif grp == 10; si < 62f0 && return 0f0; a = 157.94643f0 - 1.78125f0*si + .014330357f0*si*si
    end
    a < 0f0 && (a = 0f0)                                  # ADJMAI = ADJMAI*POINTS/10 (POINTS=10 ⇒ ×1), floored 0
    return a
end

# adjmai.f DO 2 I=1,32 — only the first 32 ISP entries are searched (a live quirk: codes 21,22,211,212 as SITE
# species fall through to illegal ⇒ ADJMAI=0). Replicated faithfully.
@inline function bm_adjmai(inspec::Integer, sindex::Real)::Float32
    si = Float32(sindex)
    inspec >= 300 && return _bm_adjmai_group(8, si)
    grp = 0
    @inbounds for i in 1:32
        if inspec == _BM_ADJ_ISP[i]; grp = Int(_BM_ADJ_IMAP[i]); break; end
    end
    grp == 0 && return 0f0
    return _bm_adjmai_group(grp, si)
end

# maical.f RMAI — stand-scalar adjusted MAI for the site species (default DF/idx3; SITEAR default 140), capped 128.
@inline function bm_rmai(p)::Float32
    isisp = Int(p.site_species); isisp == 0 && (isisp = 3)
    (isisp < 1 || isisp > 18) && (isisp = 3)
    sssi = p.sp_site_index[isisp]; sssi == 0f0 && (sssi = 140f0)
    rmai = bm_adjmai(_BM_ISPNUM[isisp], sssi)
    rmai > 128f0 && (rmai = 128f0)
    return rmai
end

@inline function bm_dubscr(rng, sp::Integer, d, h, ba, tpccf, avh, rmai)::Float32
    hf = Float32(h); hf <= 0f0 && (hf = 0.1f0)
    cr = BM_BCR2[sp]*hf + BM_BCR1[sp]*Float32(d) + BM_BCR5[sp]*Float32(tpccf) +
         BM_BCR6[sp]*(Float32(avh)/hf) + BM_BCR8[sp]*Float32(avh) + BM_BCR3[sp]*Float32(ba) +
         BM_BCR9[sp]*(Float32(ba)*Float32(tpccf)) + BM_BCR10[sp]*Float32(rmai) + BM_BCR0[sp]
    sd = BM_DUBSD[sp]
    fcr = 0f0
    while true                                             # dubscr.f label 10: reject |FCR|>SD
        fcr = bachlo(rng, 0f0, sd)
        abs(fcr) > sd && continue
        break
    end
    if sp == 13 || sp == 14 || sp == 16 || sp == 18        # CASE(13,14,16,18): linear rescale
        cr = cr + fcr
        cr = ((cr - 1f0)*10f0 + 1f0)/100f0
    else                                                   # CASE(1:12,15,17): logistic
        abs(cr + fcr) >= 86f0 && (cr = 86f0)               # faithful: sets +86 regardless of sign
        cr = 1f0/(1f0 + exp(cr + fcr))
    end
    cr < 0.05f0 && (cr = 0.05f0); cr > 0.95f0 && (cr = 0.95f0)
    return cr
end

# bm/crown.f — Weibull crown-ratio (all species; small trees D<1 → REGENT). RELSDI=SDIAC/SDIDEF,
# ACRNEW=C0+C1·RELSDI·100, Weibull A/B/C (B<1→1, C<2→2), SCALE=1−0.00167·(RELDEN−100), rank-based X.
function crown_ratio_update!(s::StandState, ::BlueMountains; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    relden = p.relative_density; sdiac = crown_sdi
    p_pccf = s.density.point_ccf
    rmai = lstart ? bm_rmai(p) : 0f0                        # maical RMAI (stand scalar) for the DUBSCR crown dub
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
        if d < 1f0 && lstart                               # bm/crown.f:336 label 58 — D<1 missing-CR at LSTART → DUBSCR
            pt = Int(t.plot_id[i])
            tpccf = (1 <= pt <= length(p_pccf)) ? p_pccf[pt] : 0f0
            cr = bm_dubscr(s.rng, sp, d, h, p.basal_area, tpccf, p.avg_height, rmai)
            icri = trunc(Int, cr*100f0 + 0.5f0)
            (d >= BM_CR_DLOW[sp] && d <= BM_CR_DHI[sp]) && (icri = trunc(Int, Float32(icri) * BM_CRNMLT[sp]))
            icri > 95 && (icri = 95)
            (icri < 10 && BM_CRNMLT[sp] == 1f0) && (icri = 10)
            icri < 1 && (icri = 1)
            t.crown_pct[i] = Int32(icri)
            continue
        end
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

# bm/cratet.f LSTART crown-init: DENSE runs over the FULL inventory (live + HISTORY 6-9 standing-dead records),
# and that dead-inclusive density (BA/AVH/point-CCF) is what CROWN→DUBSCR sees when dubbing the D<1 / missing-CR
# LIVE trees. jl partitions the dead into t.n+1:t.n+ndead, so it computes the dead-inclusive scalars by temporarily
# extending the live range, then restores live-only (the grow cycle recomputes density before use). Measured on
# 504443988126144: dead-inclusive BA 55.56 / AVH 85.07 / TPCCF 84.3 (vs live-only 35/45/63) — matches live DUBSCR.
function bm_crown_init_lstart!(s::StandState)
    t = s.trees
    nlive = t.n
    if t.ndead > 0
        t.n = nlive + t.ndead
        # #151: dense.f:83-87 — in the CRATET backdating DENSE, standing-dead records get WK3=DBH EXCEPT
        # IMC(I)==9 (HISTORY 8,9, older-dead) which LOAD DBH=0. Only HISTORY 6,7 (dead ≤5yr, IMC=7) keep their
        # DBH. So HISTORY 8,9 contribute 0 to the DBH-based density (BA/CCF/SDI) while their height still counts
        # toward AVH (stand_top_height, which live does NOT zero). Replicate by zeroing the 8/9 DBH for this pass.
        # AVHT40/DENSE top-height (dense.f:285-297) sums HT over the 40 largest-DBH TPA using IND, the
        # descending-REAL-DBH sort — NOT WK3. WK3 (with IMC9→0) drives only the BA/CCF/SDI accumulation.
        # So the dead HISTORY 8/9 heights DO enter AVH, ranked by their real DBH. Compute AVH from real DBH
        # FIRST (before the WK3-zeroing), then restore it after compute_density! overwrites it with the
        # zeroed-DBH sort. Without this, the zeroed dead sink below the live seedlings and their heights are
        # lost ⇒ DUBSCR sees AVH≈seedling-height instead of the dead-inclusive top height (measured on
        # 449747082489998: live DUBSCR AVH 67.34 vs jl 1.01 ⇒ seedling crowns dubbed 80 not the capped 95 ⇒
        # over-vigorous small-tree height/DBH growth, BA/SDI/CCF/QMD one-directionally high).
        avht_real = stand_top_height(s)    # real-DBH IND sort, real HT (AVHT40 over live + all dead records)
        saved = Tuple{Int,Float32}[]
        @inbounds for i in (nlive + 1):(nlive + t.ndead)
            (t.history[i] == 8 || t.history[i] == 9) || continue
            push!(saved, (i, t.dbh[i])); t.dbh[i] = 0f0
        end
        compute_density!(s)                # dead-inclusive BA / point-CCF (CRATET DENSE over all inv records)
        @inbounds for (i, d) in saved; t.dbh[i] = d; end
        s.plot.avg_height = avht_real      # AVHT40 top height from real DBH (dead heights included), not the WK3 sort
        t.n = nlive
    end
    crown_ratio_update!(s, s.variant; lstart = true)   # DUBSCR-dub live D<1 seedlings + Weibull-dub missing-CR overstory
    compute_density!(s)                    # restore live-only density so nothing downstream sees the dead-inclusive BA
    return s
end
