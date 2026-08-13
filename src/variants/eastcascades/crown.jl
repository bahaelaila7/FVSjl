# =============================================================================
# crown.jl (eastcascades) — EC crown ratio (ec/crown.f + ec/dubscr.f) + CCF (ec/ccfcal.f). Chunks 5 + density.
#
# Rank-based Weibull crown ratio, SAME formula as WC/PN (crown_ratio_update! mirrors the shared dgdriv
# framework: B floor 3, C floor 2, scale=1-0.00167(RELDEN-100), pdifpy ±1%/yr limiting, crmax cap), but
# EC's crown.f is NOT group-indexed — every WEIB/C0/C1 array is per-species (identity, 32-long). No
# redwood branch (EC sp 17 = subalpine larch, a normal Weibull species). DUBSCR (sub-1"/dead crown dub)
# is the EC bespoke 2-branch model (EC-form logistic {1:10,12,16,31}; WC-form linear {11,13:15,17:30,32}).
# CCF (ec/ccfcal.f) is per-species RD1/RD2/RD3 (D≥1") + RDA·D^RDB (0.1<D<1"), no INDCCF.
# =============================================================================

# ec/crown.f DATA — per-species (index = ISPC) Weibull crown coefficients.
const EC_CROWN_WEIBA  = Float32[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0]
const EC_CROWN_WEIBB0 = Float32[0.08106,0.00603,-0.28295,-0.09734,-0.01129,-0.09734,-0.00047,-0.15678,0.08247,0.08106,0.490848,-0.01129,0.196054,0.196054,-0.135807,-0.09734,0.196054,-0.811424,-0.238295,-0.818809,-0.818809,-1.112738,-0.238295,-0.238295,-0.238295,-0.238295,-0.238295,-0.238295,-0.238295,-0.238295,-0.01129,-0.238295]
const EC_CROWN_WEIBB1 = Float32[1.10253,1.12276,1.18232,1.14675,1.11665,1.14675,1.13172,1.14894,1.10804,1.10253,1.014138,1.11665,1.073909,1.073909,1.147712,1.14675,1.073909,1.056190,1.180163,1.054176,1.054176,1.123138,1.180163,1.180163,1.180163,1.180163,1.180163,1.180163,1.180163,1.180163,1.11665,1.180163]
const EC_CROWN_WEIBC0 = Float32[1.04477,2.73400,3.03400,2.71600,3.35500,2.71600,2.22700,3.05300,1.45931,1.04477,3.164558,3.35500,0.345647,0.345647,3.017494,2.71600,0.345647,-3.831124,3.044134,-2.366108,-2.366108,2.533158,3.044134,3.044134,3.044134,3.044134,3.044134,3.044134,3.044134,3.044134,3.35500,3.044134]
const EC_CROWN_WEIBC1 = Float32[0.42828,0,0,0,0,0,0,0,0.25495,0.42828,0,0,0.620145,0.620145,0,0,0.620145,1.401938,0,1.202413,1.202413,0,0,0,0,0,0,0,0,0,0,0]
const EC_CROWN_C0     = Float32[5.23986,4.98675,4.99727,4.79981,5.74915,4.79981,3.85379,6.04394,6.00795,5.23986,5.488532,5.74915,5.417431,5.417431,5.568864,4.79981,5.417431,5.200550,4.625125,4.420000,4.420000,4.120478,4.625125,4.625125,4.625125,4.625125,4.625125,4.625125,4.625125,4.625125,5.74915,4.625125]
const EC_CROWN_C1     = Float32[-0.02569,-0.02466,-0.01043,-0.00653,-0.01090,-0.00653,-0.00795,-0.01825,-0.02301,-0.02569,-0.007173,-0.01090,-0.011608,-0.011608,-0.021293,-0.00653,-0.011608,-0.014890,-0.016042,-0.010660,-0.010660,-0.006357,-0.016042,-0.016042,-0.016042,-0.016042,-0.016042,-0.016042,-0.016042,-0.016042,-0.01090,-0.016042]

# ec/dubscr.f DATA — per-species (no BCR4/BCR7).
const EC_DUB_BCR0  = Float32[-1.669490,-1.669490,-0.426688,-0.426688,-0.426688,-0.426688,-1.669490,-0.426688,-0.426688,-1.669490,7.558538,-2.19723,6.489813,6.489813,8.042774,-0.426688,6.489813,7.558538,9.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,-2.19723,5.0]
const EC_DUB_BCR1  = Float32[-0.209765,-0.209765,-0.093105,-0.093105,-0.093105,-0.093105,-0.209765,-0.093105,-0.093105,-0.209765,0,0,0,0,0,-0.093105,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
const EC_DUB_BCR2  = Float32[0,0,0.022409,0.022409,0.022409,0.022409,0,0.022409,0.022409,0,-0.015637,0,-0.029815,-0.029815,0.007198,0.022409,-0.029815,-0.015637,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
const EC_DUB_BCR3  = Float32[0.003359,0.003359,0.002633,0.002633,0.002633,0.002633,0.003359,0.002633,0.002633,0.003359,-0.009064,0,-0.009276,-0.009276,-0.016163,0.002633,-0.009276,-0.009064,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
const EC_DUB_BCR5  = Float32[0.011032,0.011032,0,0,0,0,0.011032,0,0,0.011032,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
const EC_DUB_BCR6  = Float32[0,0,-0.045532,-0.045532,-0.045532,-0.045532,0,-0.045532,-0.045532,0,0,0,0,0,0,-0.045532,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
const EC_DUB_BCR8  = Float32[0.017727,0.017727,0,0,0,0,0.017727,0,0,0.017727,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
const EC_DUB_BCR9  = Float32[-0.000053,-0.000053,0.000022,0.000022,0.000022,0.000022,-0.000053,0.000022,0.000022,-0.000053,0,0,0,0,0,0.000022,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
const EC_DUB_BCR10 = Float32[0.014098,0.014098,-0.013115,-0.013115,-0.013115,-0.013115,0.014098,-0.013115,-0.013115,0.014098,0,0,0,0,0,-0.013115,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
const EC_DUB_CRSD  = Float32[0.5000,0.5000,0.6957,0.6957,0.6957,0.9310,0.6124,0.6957,0.6957,0.4942,1.9658,0.2,2.0426,2.0426,1.3167,0.9310,2.0426,1.9658,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.2,0.5]

# ec/dubscr.f — 6-arg DUBSCR(ISPC,D,H,CR,TPCT,TPCCF); BA/AVH/RMAI from common.
@inline function _ec_dubscr(rng, sp::Int, d::Real, h::Real, ba::Real, avh::Real, rmai::Real, tpccf::Real)::Float32
    H = Float32(h); BA = Float32(ba); AVH = Float32(avh)
    if sp in EC_NATIVE_SP                              # EC-form logistic
        cr = EC_DUB_BCR0[sp] + EC_DUB_BCR1[sp]*Float32(d) + EC_DUB_BCR2[sp]*H + EC_DUB_BCR3[sp]*BA +
             EC_DUB_BCR5[sp]*Float32(tpccf) + EC_DUB_BCR6[sp]*(AVH/H) + EC_DUB_BCR8[sp]*AVH +
             EC_DUB_BCR9[sp]*(BA*Float32(tpccf)) + EC_DUB_BCR10[sp]*Float32(rmai)
        sd = EC_DUB_CRSD[sp]; fcr = 0f0
        while true; fcr = bachlo(rng, 0f0, sd); abs(fcr) > sd && continue; break; end
        abs(cr + fcr) >= 86f0 && (cr = 86f0)
        cr = 1f0 / (1f0 + exp(cr + fcr))
    else                                                # WC-form linear
        cr = EC_DUB_BCR0[sp] + EC_DUB_BCR2[sp]*H + EC_DUB_BCR3[sp]*BA
        sd = EC_DUB_CRSD[sp]; fcr = 0f0
        while true; fcr = bachlo(rng, 0f0, sd); abs(fcr) > sd && continue; break; end
        cr = ((cr + fcr - 1f0) * 10f0 + 1f0) / 100f0
    end
    cr > 0.95f0 && (cr = 0.95f0); cr < 0.05f0 && (cr = 0.05f0)
    return cr
end

function crown_ratio_update!(s::StandState, ::EastCascades; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    sd = s.coef.species
    relden = p.relative_density; sdiac = crown_sdi; ba = p.basal_area; avh = p.avg_height
    rmai = s.plot.mai_adj; dens = s.density
    key = Vector{Float32}(undef, n); idx = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        bk = wc_bratio(sd, Int(t.species[i]), t.dbh[i])
        key[i] = t.dbh[i] + t.diam_growth[i] / bk; idx[i] = Int32(i)
    end
    _rdpsrt!(key, idx; lseq = false)
    isort = Vector{Int32}(undef, n)
    @inbounds for jj in 1:n; isort[idx[jj]] = Int32(n - jj + 1); end
    @inbounds for i in 1:n
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (sp < 1 || sp > 32) && continue
        (lstart && t.crown_pct[i] > 0) && continue
        icr = Int(t.crown_pct[i])
        if d < 1f0 && lstart
            icr != 0 && continue
            pt_i = Int(t.plot_id[i])
            tpccf = (1 <= pt_i <= length(dens.point_ccf)) ? dens.point_ccf[pt_i] : relden
            cr = _ec_dubscr(s.rng, sp, d, h, ba, avh, rmai, tpccf)
            icri = trunc(Int, cr * 100f0 + 0.5f0)
            icri > 95 && (icri = 95); icri < 10 && (icri = 10); icri < 1 && (icri = 1)
            t.crown_pct[i] = Int32(icri); continue
        end
        relsdi = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0
        relsdi > 1.5f0 && (relsdi = 1.5f0)
        acrnew = EC_CROWN_C0[sp] + EC_CROWN_C1[sp] * relsdi * 100f0
        A = EC_CROWN_WEIBA[sp]
        B = EC_CROWN_WEIBB0[sp] + EC_CROWN_WEIBB1[sp] * acrnew; B < 3f0 && (B = 3f0)
        C = EC_CROWN_WEIBC0[sp] + EC_CROWN_WEIBC1[sp] * acrnew; C < 2f0 && (C = 2f0)
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
            (icri < 10) && (icri = trunc(Int, crmax + 0.5f0))
            Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
        end
        if lstart && t.trunc[i] != 0
            hn = Float32(t.norm_ht[i]) / 100f0; hd = hn - Float32(t.trunc[i]) / 100f0
            cl = (Float32(icri) / 100f0) * hn - hd
            icri = trunc(Int, (cl * 100f0 / hn) + 0.5f0)
        end
        icri > 95 && (icri = 95); icri < 10 && (icri = 10); icri < 1 && (icri = 1)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end

# ec/ccfcal.f — per-species CCF (no INDCCF). D≥1": RD1+D·RD2+D²·RD3; 0.1<D<1": RDA·D^RDB.
const EC_CCF_RD1 = Float32[0.03,0.02,0.0388,0.04,0.03,0.04,0.01925,0.03,0.03,0.0219,0.03758,0.03,0.0204,0.01925,0.02453,0.04,0.0194,0.0194,0.0194,0.0204,0.0204,0.03561,0.0204,0.0160,0.0204,0.0204,0.0204,0.0204,0.0204,0.0204,0.03,0.0204]
const EC_CCF_RD2 = Float32[0.0167,0.0148,0.0269,0.0270,0.0238,0.027,0.01676,0.0173,0.0216,0.0169,0.0233,0.0215,0.0246,0.0168,0.0115,0.027,0.0142,0.0142,0.0142,0.0246,0.0246,0.02731,0.0246,0.0167,0.0246,0.0246,0.0246,0.0246,0.0246,0.0246,0.0215,0.0246]
const EC_CCF_RD3 = Float32[0.00230,0.00338,0.00466,0.00405,0.00490,0.00405,0.00365,0.00259,0.00405,0.00325,0.00361,0.00363,0.0074,0.00365,0.00134,0.00405,0.00261,0.00261,0.00261,0.0074,0.0074,0.00524,0.0074,0.00434,0.0074,0.0074,0.0074,0.0074,0.0074,0.0074,0.00363,0.0074]
const EC_CCF_RDA = Float32[0.009884,0.007244,0.017299,0.015248,0.008915,0.015248,0.009187,0.007875,0.011402,0.007813,0,0.011109,0,0,0,0.015248,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.011109,0]
const EC_CCF_RDB = Float32[1.6667,1.8182,1.5571,1.7333,1.7800,1.7333,1.7600,1.7360,1.7560,1.7780,0,1.7250,0,0,0,1.7333,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.7250,0]

@inline function ec_tree_ccf(sp::Integer, d::Real)::Float32
    (sp < 1 || sp > 32) && return 0f0
    D = Float32(d)
    D >= 1.0f0 && return EC_CCF_RD1[sp] + EC_CCF_RD2[sp] * D + EC_CCF_RD3[sp] * D * D
    if sp in EC_NATIVE_SP                               # EC-form small-tree CCF
        D > 0.1f0 && return EC_CCF_RDA[sp] * D^EC_CCF_RDB[sp]
        return 0.001f0
    end
    return D * (EC_CCF_RD1[sp] + EC_CCF_RD2[sp] + EC_CCF_RD3[sp])  # WC-form small-tree CCF
end
