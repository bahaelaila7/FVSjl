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

# nc/ccfcal.f MODE=1 per-tree CCF: RD1+D·RD2+D²·RD3 (D≥1); RDA·D^RDB (0.1<D<1); 0.001 (D≤0.1). (×P by caller.)
const NC_RD1 = Float32[.0388,.0392,.0388,.0690,.0212,.0194,.0204,.0356,.0172,.0219,.0356,.0388]
const NC_RD2 = Float32[.0269,.0180,.0269,.0225,.0167,.0142,.0246,.0273,.00877,.0169,.0273,.0269]
const NC_RD3 = Float32[.00466,.00207,.00466,.00183,.00330,.00261,.0074,.00524,.00112,.00325,.00524,.00466]
const NC_RDA = Float32[0.009884,0.007244,0.017299,0.015248,0.011109,0.008915,0.009187,0.007875,0.011402,0.007813,0.011109,0.017299]
const NC_RDB = Float32[1.6667,1.8182,1.5571,1.7333,1.7250,1.7800,1.7600,1.7360,1.7560,1.7780,1.7250,1.5571]

# nc/dubscr.f — dub crown ratio for sub-1" inventory trees with MISSING crown (crown.f:265,377). Logistic on
# BA / per-point PCCF / top-40 AVH / RMAI (=50, grinit.f:139); sp12 RW on HDR/PRD/D:QMDPLT. The FCR random
# perturbation (BACHLO, |FCR|≤CRSD) is a small logit-scale jitter: for the dense-seedling regime the
# deterministic logit is ≲−6 ⇒ CR clamps to 0.95 REGARDLESS of FCR, so the crown VALUE is RNG-independent and
# bit-exact without a draw (jl's index-order crown loop can't reproduce FVS's species-grouped RANN order anyway,
# and NOT drawing leaves jl's stream — and the gate — untouched). Rare non-clamped sub-1" trees are cornered on
# the FCR ULP. Was SKIPPED entirely (crown_pct stayed 0 ⇒ htgr5 CR²=0 ⇒ seedlings never crossed 4.5' ⇒ no
# diameter ⇒ self-thinning never fired: the extreme-dense under-kill).
const NC_BCR0  = Float32[-1.66949,-1.66949,-0.426688,-0.426688,-0.426688,-0.426688,-1.66949,-0.426688,-0.426688,-1.66949,-2.19723,0.0]
const NC_BCR1  = Float32[-0.209765,-0.209765,-0.093105,-0.093105,-0.093105,-0.093105,-0.209765,-0.093105,-0.093105,-0.209765,0.0,0.0]
const NC_BCR2  = Float32[0.0,0.0,0.022409,0.022409,0.022409,0.022409,0.0,0.022409,0.022409,0.0,0.0,0.0]
const NC_BCR3  = Float32[0.003359,0.003359,0.002633,0.002633,0.002633,0.002633,0.003359,0.002633,0.002633,0.003359,0.0,0.0]
const NC_BCR5  = Float32[0.011032,0.011032,0.0,0.0,0.0,0.0,0.011032,0.0,0.0,0.011032,0.0,0.0]
const NC_BCR6  = Float32[0.0,0.0,-0.045532,-0.045532,-0.045532,-0.045532,0.0,-0.045532,-0.045532,0.0,0.0,0.0]
const NC_BCR8  = Float32[0.017727,0.017727,0.0,0.0,0.0,0.0,0.017727,0.0,0.0,0.017727,0.0,0.0]
const NC_BCR9  = Float32[-0.000053,-0.000053,0.000022,0.000022,0.000022,0.000022,-0.000053,0.000022,0.000022,-0.000053,0.0,0.0]
const NC_BCR10 = Float32[0.014098,0.014098,-0.013115,-0.013115,-0.013115,-0.013115,0.014098,-0.013115,-0.013115,0.014098,0.0,0.0]

# nc/dubscr.f deterministic crown (FCR=0). ba/tpccf/avh stand+point context; prd/qmdplt only for RW (sp12).
# nc/dubscr.f DATA CRSD — the DUBSCR random-error SD per species.
const NC_CRSD = Float32[0.5, 0.5, 0.6957, 0.6957, 0.6957, 0.9310, 0.6124, 0.6957, 0.6957, 0.4942, 0.200, 0.1500]

"""nc/dubscr.f — crown for a d<1\" seedling or cycle-0 dead record with a missing crown. Includes the
rejection-bounded BACHLO(0,CRSD) draw that the old port omitted (it consumed no draw ⇒ NC's stream shifted on
every stand with a missing crown). RMAI=50 (nc/grinit.f:139)."""
@inline function nc_dubscr(rng, sp::Integer, d::Real, h::Real, ba::Real, tpccf::Real, avh::Real,
                           prd::Real, qmdplt::Real)::Float32
    hf = Float32(h); df = Float32(d)
    local cr::Float32
    if sp == 12                                          # redwood — logistic on HDR/PRD/D:QMDPLT
        hdr = (hf * 12f0) / df
        cr = -1.021064f0 + 0.309296f0 * log(hdr) + 0.869720f0 * Float32(prd) - 0.116274f0 * (df / Float32(qmdplt))
    else
        cr = NC_BCR0[sp] + NC_BCR1[sp]*df + NC_BCR2[sp]*hf + NC_BCR3[sp]*Float32(ba) +
             NC_BCR5[sp]*Float32(tpccf) + NC_BCR6[sp]*(Float32(avh)/hf) + NC_BCR8[sp]*Float32(avh) +
             NC_BCR9[sp]*(Float32(ba)*Float32(tpccf)) + NC_BCR10[sp]*50f0
    end
    sd = NC_CRSD[sp]; fcr = 0f0
    while true; fcr = bachlo(rng, 0f0, sd); abs(fcr) > sd && continue; break; end
    abs(cr + fcr) >= 86f0 && (cr = 86f0)                 # dubscr.f: IF(ABS(CR+FCR).GE.86.)CR=86.
    cr = 1f0 / (1f0 + exp(cr + fcr))
    cr > 0.95f0 && (cr = 0.95f0); cr < 0.05f0 && (cr = 0.05f0)
    return cr
end

@inline function nc_tree_ccf(sp::Integer, d::Real)::Float32
    dd = Float32(d)
    if dd >= 1.0f0
        return NC_RD1[sp] + dd * NC_RD2[sp] + dd * dd * NC_RD3[sp]
    elseif dd > 0.1f0
        return NC_RDA[sp] * dd^NC_RDB[sp]
    else
        return 0.001f0
    end
end

function crown_ratio_update!(s::StandState, ::Klamath; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    sd = s.coef.species
    ba = p.basal_area; sdiac = crown_sdi
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    # rank by (D+G) descending → ISORT (largest = rank 1)
    key = Vector{Float32}(undef, n); idx = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        bk = variant_bratio(s, t.species[i], t.dbh[i], t.height[i])   # crown.f ISORT key D+DG/BRATIO — shared variant bark
        key[i] = t.dbh[i] + t.diam_growth[i] / bk; idx[i] = Int32(i)
    end
    _rdpsrt!(key, idx; lseq = false)
    isort = Vector{Int32}(undef, n)
    @inbounds for jj in 1:n; isort[idx[jj]] = Int32(n - jj + 1); end
    # nc/crown.f:87-104 — SDICAL XMAXPT + per-point Zeide ZRD (live trees), computed once; PRD/QMDPLT per point.
    xmaxpt, zrd, _ = point_zeide!(s)
    dens = s.density
    _prd(pt) = (1 <= pt <= length(xmaxpt) && xmaxpt[pt] > 0f0) ? zrd[pt] / xmaxpt[pt] : 0f0
    function _qmdplt(pt)
        baplt = (1 <= pt <= length(dens.point_ba)) ? dens.point_ba[pt] : 0f0
        tpaplt = (1 <= pt <= length(dens.point_tpa)) ? dens.point_tpa[pt] : 0f0
        q = tpaplt > 0f0 ? sqrt((baplt / tpaplt) / 0.005454f0) : 1f0
        q <= 1f0 ? 1f0 : q
    end
    # crown.f DO 70 ISPC … I=IND1(I3): SPECIES-MAJOR — the DUBSCR/RANN draws follow this order, not storage.
    @inbounds for i in species_major_order(s)
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (sp < 1 || sp > 12) && continue
        (lstart && t.crown_pct[i] > 0) && continue
        if d < 1f0 && lstart                              # crown.f:265,377 — sub-1" missing-crown ⇒ DUBSCR
            pt = Int(t.plot_id[i])
            tpccf = (1 <= pt <= length(s.density.point_ccf)) ? s.density.point_ccf[pt] : 0f0
            ptba = (1 <= pt <= length(s.density.point_ba)) ? s.density.point_ba[pt] : ba
            pttpa = (1 <= pt <= length(s.density.point_tpa)) ? s.density.point_tpa[pt] : 0f0
            qmdplt = pttpa > 0f0 ? sqrt((ptba / pttpa) / 0.005454f0) : 1f0
            qmdplt < 1f0 && (qmdplt = 1f0)
            cr = nc_dubscr(s.rng, sp, d, h, ba, tpccf, p.avg_height, _prd(pt), qmdplt)   # PRD = ZRD/XMAXPT (point)
            icri = trunc(Int, cr * 100f0 + 0.5f0)
            lo = sp == 12 ? 5 : 10
            icri > 95 && (icri = 95); icri < lo && (icri = lo); icri < 1 && (icri = 1)
            t.crown_pct[i] = Int32(icri)
            continue
        end
        icr = Int(t.crown_pct[i])
        relsdi = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0
        relsdi > 1.5f0 && (relsdi = 1.5f0)
        local crnew::Float32
        if sp == 12                                       # redwood — logistic (nc/crown.f:233-277, all POINT terms)
            pt = Int(t.plot_id[i])
            hdr = (h * 12f0) / d
            xl = -1.021064f0 + 0.309296f0 * log(hdr) + 0.869720f0 * _prd(pt) - 0.116274f0 * (d / _qmdplt(pt))
            x = 1f0 / (1f0 + exp(xl))
            x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
            crnew = x * 100f0
        else
            acrnew = NC_CRC0[sp] + NC_CRC1[sp] * relsdi * 100f0
            b = NC_WEIBB0[sp] + NC_WEIBB1[sp] * acrnew
            c = NC_WEIBC0[sp] + NC_WEIBC1[sp] * acrnew
            scale = 1.5f0 - relsdi
            scale > 1f0 && (scale = 1f0); scale < 0.30f0 && (scale = 0.30f0)
            x = d > 0f0 ? (Float32(isort[i]) / Float32(n)) * scale : rann!(s.rng) * scale
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
    # nc/crown.f:404-452 DO 79 — cycle-0 dead-record crown dub (point PRD/QMDPLT, TPCCF = PCCF(ITRE(I))).
    lstart && dub_dead_crowns!(s) do i
        pt = Int(t.plot_id[i])
        tpccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0f0
        icri_round(nc_dubscr(s.rng, Int(t.species[i]), t.dbh[i], t.height[i], ba, tpccf, p.avg_height, _prd(pt), _qmdplt(pt)))
    end
    return s
end

# nc/cratet.f LSTART crown-init = the shared `crown_init_lstart_dead_inclusive!` (bluemountains/crown.jl): nc/cratet.f
# runs the same backdated (LBKDEN) dead-inclusive DENSE before CROWN. The former NC-local copy skipped the backdating.
