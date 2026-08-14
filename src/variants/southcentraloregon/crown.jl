# =============================================================================
# crown.jl (southcentraloregon) — SO crown ratio (so/crown.f + so/dubscr.f). Chunk 5.
#
# Rank-based Weibull crown ratio, SAME framework as EC/WC/PN (crown_ratio_update! mirrors ec_crown):
# per-species (identity, 33-long) WEIB/C0/C1; B floor 3, C floor 2; SCALE=1−0.00167(RELDEN−100) clamp
# [0.30,1]; rank X=ISORT/ITRN·SCALE clamp [.05,.95]; CRNEW=(A+B·(−ln(1−X))^(1/C))·10; ±1%/yr change limit;
# CRMAX crown-length cap; NORMHT topkill reduction. so/crown.f's CRNMLT(1.0)/DLOW(0)/DHI(99)/ICFLG(0) are
# NO-OP defaults ⇒ their DLOW≤D≤DHI CRNMLT branches are inert. DBH<1" at LSTART → so/dubscr.f: a single
# CR = BCR0 + BCR1·D + BCR2·H + BCR3·BA + BCR5·TPCCF + BCR6·(AVH/H) + BCR8·AVH + BCR9·(BA·TPCCF) + BCR10·RMAI,
# a random BACHLO(0,CRSD) draw, then a species-split transform (logistic {1:8,10:14,16:18,24,32}; else linear).
# MEASURED vs FVSso_g16 crown per-tree ICRI.
# =============================================================================

# so/crown.f DATA — per-species Weibull crown coefficients (index = ISPC).
const SO_CROWN_WEIBA  = Float32[2,2,1,0,1,1,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0]
const SO_CROWN_WEIBB0 = Float32[-2.12713,-2.27296,-1.19297,0.06593,-0.94138,-1.38636,0.07609,-0.91567,0.16601,0.24916,0.07609,0.06593,-0.91567,-0.09734,-0.13581,-0.82631,0.00603,-0.01129,0.49085,0.19605,-0.23830,-1.11274,-0.23830,-0.08414,-0.23830,-0.23830,0.06607,-0.23830,-0.23830,-0.23830,-0.23830,-1.19297,-0.23830]
const SO_CROWN_WEIBB1 = Float32[1.10526,1.12434,1.12928,1.09624,1.08256,1.16801,1.10184,1.06469,1.08150,1.04831,1.10184,1.09624,1.06469,1.14675,1.14771,1.06217,1.12276,1.11665,1.01414,1.07391,1.18016,1.12314,1.18016,1.14765,1.18016,1.18016,1.10705,1.18016,1.18016,1.18016,1.18016,1.12928,1.18016]
const SO_CROWN_WEIBC0 = Float32[2.77,3.34,3.42,3.71,3.47,3.02,3.01,3.50,0.9142,4.36,3.01,3.71,3.50,2.716,3.02,3.31429,2.734,3.355,3.16,0.35,3.04,2.53,3.04,2.775,3.04,3.04,2.04714,3.04,3.04,3.04,3.04,3.42,3.04]
const SO_CROWN_WEIBC1 = Float32[0,0,0,0,0,0,0,0,0.45768,0,0,0,0,0,0,0,0,0,0,0.62015,0,0,0,0,0,0,0.15070,0,0,0,0,0,0]
const SO_CROWN_C0     = Float32[7.16846,6.59700,5.52653,6.61291,7.45097,6.17373,5.50719,6.77400,6.14578,6.41166,7.23800,6.61291,6.12779,4.79981,5.56886,6.19911,4.98675,5.74915,5.48853,5.41743,4.62512,4.12048,4.62512,4.01678,4.62512,4.62512,6.82187,4.62512,4.62512,4.62512,4.62512,5.52653,4.62512]
const SO_CROWN_C1     = Float32[-0.02375,-0.01954,0.0,-0.02182,-0.02406,-0.01795,-0.01833,0.0,-0.02781,-0.02041,0.0,-0.02182,-0.01269,-0.00653,-0.02129,-0.02216,-0.02466,-0.01090,-0.00717,-0.01161,-0.01604,-0.00636,-0.01604,-0.01516,-0.01604,-0.01604,-0.02247,-0.01604,-0.01604,-0.01604,-0.01604,0.0,-0.01604]

# so/dubscr.f DATA — sub-1" crown dub (index = ISPC).
const SO_DUB_BCR0  = Float32[-1.66949,-1.669490,-0.426688,-0.426688,-0.426688,-0.426688,-1.669490,-0.426688,8.042774,-1.669490,-2.19723,-0.426688,-0.426688,-0.426688,8.042774,-1.66949,-1.669490,-0.426688,7.558538,6.489813,5.0,5.0,5.0,-0.426688,5.0,5.0,6.489813,5.0,5.0,5.0,5.0,-0.426688,5.0]
const SO_DUB_BCR1  = Float32[-0.209765,-0.209765,-0.093105,-0.093105,-0.093105,-0.093105,-0.209765,-0.093105,0.0,-0.209765,0.0,-0.093105,-0.093105,-0.093105,0.0,-0.209765,-0.209765,-0.093105,0.0,0.0,0.0,0.0,0.0,-0.093105,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.093105,0.0]
const SO_DUB_BCR2  = Float32[0.0,0.0,0.022409,0.022409,0.022409,0.022409,0.0,0.022409,0.007198,0.0,0.0,0.022409,0.022409,0.022409,0.007198,0.0,0.0,0.022409,-0.015637,-0.029815,0.0,0.0,0.0,0.022409,0.0,0.0,-0.029815,0.0,0.0,0.0,0.0,0.022409,0.0]
const SO_DUB_BCR3  = Float32[0.003359,0.003359,0.002633,0.002633,0.002633,0.002633,0.003359,0.002633,-0.016163,0.003359,0.0,0.002633,0.002633,0.002633,-0.016163,0.003359,0.003359,0.002633,-0.009064,-0.009276,0.0,0.0,0.0,0.002633,0.0,0.0,-0.009276,0.0,0.0,0.0,0.0,0.002633,0.0]
const SO_DUB_BCR5  = Float32[0.011032,0.011032,0.0,0.0,0.0,0.0,0.011032,0.0,0.0,0.011032,0.0,0.0,0.0,0.0,0.0,0.011032,0.011032,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const SO_DUB_BCR6  = Float32[0.0,0.0,-0.045532,-0.045532,-0.045532,-0.045532,0.0,-0.045532,0.0,0.0,0.0,-0.045532,-0.045532,-0.045532,0.0,0.0,0.0,-0.045532,0.0,0.0,0.0,0.0,0.0,-0.045532,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.045532,0.0]
const SO_DUB_BCR8  = Float32[0.017727,0.017727,0.0,0.0,0.0,0.0,0.017727,0.0,0.0,0.017727,0.0,0.0,0.0,0.0,0.0,0.017727,0.017727,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const SO_DUB_BCR9  = Float32[-0.000053,-0.000053,0.000022,0.000022,0.000022,0.000022,-0.000053,0.000022,0.0,-0.000053,0.0,0.000022,0.000022,0.000022,0.0,-0.000053,-0.000053,0.000022,0.0,0.0,0.0,0.0,0.0,0.000022,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.000022,0.0]
const SO_DUB_BCR10 = Float32[0.014098,0.014098,-0.013115,-0.013115,-0.013115,-0.013115,0.014098,-0.013115,0.0,0.014098,0.0,-0.013115,-0.013115,-0.013115,0.0,0.014098,0.014098,-0.013115,0.0,0.0,0.0,0.0,0.0,-0.013115,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.013115,0.0]
const SO_DUB_CRSD  = Float32[0.5000,0.5000,0.6957,0.6957,0.6957,0.9310,0.6124,0.6957,1.3167,0.4942,0.200,0.6957,0.6957,0.6957,1.3167,0.5000,0.5000,0.6957,1.9658,2.0426,0.500,0.500,0.500,0.9310,0.500,0.500,2.0426,0.500,0.500,0.500,0.500,0.6957,0.500]

# so/dubscr.f SELECT CASE — logistic species; the rest use the linear transform.
const SO_DUB_LOGISTIC = Set{Int}([1,2,3,4,5,6,7,8,10,11,12,13,14,16,17,18,24,32])

# so/dubscr.f — CR for a DBH<1" inventory seedling (BA/AVH/RMAI from the stand).
@inline function _so_dubscr(rng, sp::Int, d::Float32, h::Float32, ba::Float32, avh::Float32,
                            rmai::Float32, tpccf::Float32)::Float32
    cr = SO_DUB_BCR2[sp]*h + SO_DUB_BCR1[sp]*d + SO_DUB_BCR5[sp]*tpccf + SO_DUB_BCR6[sp]*(avh/h) +
         SO_DUB_BCR8[sp]*avh + SO_DUB_BCR3[sp]*ba + SO_DUB_BCR9[sp]*(ba*tpccf) + SO_DUB_BCR10[sp]*rmai +
         SO_DUB_BCR0[sp]
    sd = SO_DUB_CRSD[sp]; fcr = 0f0
    while true; fcr = bachlo(rng, 0f0, sd); abs(fcr) > sd && continue; break; end
    if sp in SO_DUB_LOGISTIC
        abs(cr + fcr) >= 86f0 && (cr = 86f0)
        cr = 1f0 / (1f0 + exp(cr + fcr))
    else
        cr = ((cr - 1f0) * 10f0 + 1f0) / 100f0
    end
    cr > 0.95f0 && (cr = 0.95f0); cr < 0.05f0 && (cr = 0.05f0)
    return cr
end

function crown_ratio_update!(s::StandState, ::SouthCentralOregon; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    sd = s.coef.species
    relden = p.relative_density; sdiac = crown_sdi; ba = p.basal_area; avh = p.avg_height
    rmai = p.mai_adj; dens = s.density
    # rank the trees by projected DBH (so/crown.f ISORT via RDPSRT on D+DG/BARK), descending
    key = Vector{Float32}(undef, n); idx = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        bk = so_bratio(sd, Int(t.species[i]), t.dbh[i])
        key[i] = t.dbh[i] + t.diam_growth[i] / bk; idx[i] = Int32(i)
    end
    _rdpsrt!(key, idx; lseq = false)
    isort = Vector{Int32}(undef, n)
    @inbounds for jj in 1:n; isort[idx[jj]] = Int32(n - jj + 1); end
    @inbounds for i in 1:n
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (sp < 1 || sp > 33) && continue
        (lstart && t.crown_pct[i] > 0) && continue
        icr = Int(t.crown_pct[i])
        if d < 1f0 && lstart                                  # so/crown.f:58 → DUBSCR
            icr != 0 && continue
            pt_i = Int(t.plot_id[i])
            tpccf = (1 <= pt_i <= length(dens.point_ccf)) ? dens.point_ccf[pt_i] : relden
            cr = _so_dubscr(s.rng, sp, d, h, ba, avh, rmai, tpccf)
            icri = trunc(Int, cr * 100f0 + 0.5f0)
            icri > 95 && (icri = 95); icri < 10 && (icri = 10); icri < 1 && (icri = 1)
            t.crown_pct[i] = Int32(icri); continue
        end
        relsdi = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0
        relsdi > 1.5f0 && (relsdi = 1.5f0)
        acrnew = SO_CROWN_C0[sp] + SO_CROWN_C1[sp] * relsdi * 100f0
        A = SO_CROWN_WEIBA[sp]
        B = SO_CROWN_WEIBB0[sp] + SO_CROWN_WEIBB1[sp] * acrnew; B < 3f0 && (B = 3f0)
        C = SO_CROWN_WEIBC0[sp] + SO_CROWN_WEIBC1[sp] * acrnew; C < 2f0 && (C = 2f0)
        scale = 1f0 - 0.00167f0 * (relden - 100f0)
        scale > 1f0 && (scale = 1f0); scale < 0.30f0 && (scale = 0.30f0)
        x = d > 0f0 ? (Float32(isort[i]) / Float32(n)) * scale : 0.5f0 * scale
        x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
        crnew = (A + B * (-log(1f0 - x))^(1f0 / C)) * 10f0
        if !(lstart || icr == 0)                              # crown CHANGE, ±1%/yr limit
            chg = crnew - Float32(icr); pdifpy = chg / Float32(icr) / fint
            pdifpy > 0.01f0 && (chg = Float32(icr) * 0.01f0 * fint)
            pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
            crnew = Float32(icr) + chg
        end
        icri = trunc(Int, crnew + 0.5f0)
        if !(lstart || icr == 0)                              # CRMAX crown-length cap
            htg = t.ht_growth[i]; crln = h * Float32(icr) / 100f0
            crmax = (crln + htg) / (h + htg) * 100f0
            Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
            (icri < 10) && (icri = trunc(Int, crmax + 0.5f0))
        end
        if lstart && t.trunc[i] != 0                          # topkill NORMHT reduction
            hn = Float32(t.norm_ht[i]) / 100f0; hd = hn - Float32(t.trunc[i]) / 100f0
            cl = (Float32(icri) / 100f0) * hn - hd
            icri = trunc(Int, (cl * 100f0 / hn) + 0.5f0)
        end
        icri > 95 && (icri = 95); icri < 10 && (icri = 10); icri < 1 && (icri = 1)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end
