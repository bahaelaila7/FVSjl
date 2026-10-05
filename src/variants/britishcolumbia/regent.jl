# =============================================================================
# regent.jl (britishcolumbia) — BC small-tree growth engine hook (canada/bc/regent.f, V3).
#
# CHUNK 6. `small_tree_growth!(s, stash, ::BritishColumbia; fint)` overrides DG/HTG for small trees,
# mirroring KT (src/variants/kootenai/regent.jl:289-411) with BC swaps:
#   • height increment = bc_v3_sthg (VALIDATED regression) accumulated over subcycles;
#   • ZZRAN (V3): HTGR = max(0, HTGR1 + ZZRAN·ST_COEF%SD)  [cornered class, seed 55329];
#   • DBH-from-height = bc_st_dbh (POWER form HHT1·(H-4.5)^HHT2 + DADJ);
#   • per-subcycle density via bc_tree_ccf (chunk-5 CCF spine).
# ⚠ V2 regime + SBS logistic branch TODO. Bounds XMINV3/XMAXV3: sp14=2/4 confirmed, others fill(2/4) TODO.
# =============================================================================

"""Small-tree site setup (regent MORCON-analog): resolve ST_COEF ip + RHCON per species (V3)."""
function bc_regcons!(s::StandState)
    ctl = s.control
    zone, series = bc_stand_zone(s)
    nsp = nspecies(BritishColumbia())
    ip = zeros(Int, nsp); rhcon = zeros(Float32, nsp)
    @inbounds for sp in 1:nsp
        ip[sp] = bc_resolve_stcoef(sp, series, zone)
        # V3 RHCON (regent.f:2112): 0 base, ln(RCOR2) when small-tree HT calibrated (default RCOR2=1 ⇒ 0)
        rhcon[sp] = (ctl.regh_cor2_on && ctl.regh_cor2[sp] > 0f0) ? flog(ctl.regh_cor2[sp]) : 0f0
    end
    return ip, rhcon
end

"""BC `small_tree_growth!` — V3 small-tree height + DBH increment (regent.f main body)."""
function small_tree_growth!(s::StandState, stash, ::BritishColumbia; fint::Float32 = 10.0f0)
    p, t = s.plot, s.trees
    t.n == 0 && return s
    n = t.n
    zone, series = bc_stand_zone(s)
    ip, rhcon = bc_regcons!(s)
    ba = p.basal_area; relden = p.relative_density; avh = p.avg_height
    aspect = p.aspect; slope = p.slope; dgsd = s.control.dg_sd
    regyr = 5.0f0
    # V2 (LV2ATV) small-tree: CON = RHCON(habitat) + HCOR(htg_cor_small); exp-form HTGRL; XMAXV2 bounds;
    # multiplicative ZZRAN (HSIGMA). No ST_COEF/ip needed. (regent.f LV2ATV branches.)
    v2 = bc_lv2atv(zone)
    nsp = nspecies(BritishColumbia())
    con_v2 = zeros(Float32, nsp)
    if v2
        regch_v2 = bc_v2_regch(aspect, slope)
        @inbounds for sp in 1:nsp
            con_v2[sp] = bc_v2_rhcon(sp, regch_v2) + s.calib.htg_cor_small[sp]
        end
    end
    # subcycle count + lengths (regent.f:186-203, mirror KT)
    ntyr = Int(round(fint)); iyr = Int(regyr)
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1); nper < 1 && (nper = 1)
    kper = zeros(Int, nper); itot = ntyr; nn = nper
    @inbounds for i in 1:nper
        if nn == 1; kper[i] = itot; break; end
        kper[i] = itot ÷ nn; itot -= kper[i]; nn -= 1
    end
    # per-subcycle stand density from LARGE trees (regent.f:236-256, via bc_tree_ccf)
    banext = fill(ba, nper); rdnext = fill(relden, nper)
    if nper > 1
        @inbounds for i in 1:n
            d1 = t.dbh[i]; d1 < 3.0f0 && continue
            sp = Int(t.species[i]); pr = t.tpa[i]
            bark = bc_bratio(sp)
            d2 = d1 + t.diam_growth[i] / bark
            b1 = 0.005454154f0*d1*d1; b2 = 0.005454154f0*d2*d2
            # CCFCAL (regent.f:1173-1174) computes the per-tree CCF WITH the tree's TPA (C=ccf·P), so CI=P·Δccf/10
            # and CI/P·PN recovers Δccf/10·PN. Using P=1 here dropped the ×TPA ⇒ RDNEXT barely moved (70.4 vs the
            # oracle's −163.79) ⇒ the small-tree HTGRL CCF term collapsed ⇒ ~4× height under-growth. Use pr.
            c1 = bc_tree_ccf(sp, d1, pr); c2 = bc_tree_ccf(sp, d2, pr)
            bi = (b2 - b1) / 10f0; ci = (c2 - c1) / 10f0
            k = 0
            for j in 2:nper
                k += kper[j-1]; pn = pr * fpowi(0.985f0, k)       # 0.985**K = __powisf2
                rdnext[j] += k * ci / pr * pn
                banext[j] += k * bi * pn
            end
        end
    end
    delmax = (avh/36f0)*(0.01232f0*relden - 1.75f0); delmax > 0f0 && (delmax = 0f0)
    wk3 = Float32[t.height[i] for i in 1:n]
    # subcycle height accumulation (regent.f:1273-1458). Each subcycle ALSO adds the SMALL trees' own CCF/BA
    # contribution to the NEXT subcycle's density (regent.f:1450-1458) — omitting it collapsed RDNEXT (e.g.
    # 74 vs the oracle's −163.79 for garbage-height trees, whose power-form D1 from H1 dwarfs the actual DBH),
    # which killed the small-tree HTGRL CCF term ⇒ ~4× height under-growth. KY = cumulative years thru subcycle j.
    ky = 0
    # regent.f:1289-1321 DO 16 ISPC / DO 15 I3 / I=IND1(I3): the subcycle walks species-major IND1, and its RDNEXT/BANEXT
    # (J+1) accumulations are REAL*4 sums in that order (MEASURED FVSbc_dbfix SkyRanch-Control cycle-2 RDNEXT(2) 4084C978
    # live / 4084C842 in record order).
    sub_order = species_major_order(s)
    @inbounds for j in 1:nper
        baj = banext[j]; rdj = rdnext[j]; kpj = Float32(kper[j]); ky += kper[j]
        decay = fpowi(0.985f0, ky)                         # 0.985**KY = __powisf2
        for i in sub_order
            sp = Int(t.species[i]); d = t.dbh[i]
            xmx_sp = v2 ? BC_RG_V2_XMAX[sp] : BC_RG_XMAX[sp]
            (d >= xmx_sp || t.tpa[i] <= 0f0 || (!v2 && ip[sp] < 1)) && continue
            h1 = wk3[i]; pct = t.crown_ratio[i]; pr = t.tpa[i]
            bal = (1f0 - pct/100f0) * baj
            incr = v2 ? fexp(max(-40f0, con_v2[sp] + BC_RG_V2_RHLH[sp]*flog(h1) +
                                BC_RG_V2_RHCCF[sp]*rdj + BC_RG_V2_RHBAL[sp]*bal)) :
                        bc_v3_sthg(sp, ip[sp], h1, bal, rdj, rhcon[sp], aspect, slope)
            xrhgro = active_multiplier(s.control, :regh, sp, current_cycle_year(s))
            h2 = h1 + incr * (kpj/regyr) * xrhgro
            wk3[i] = h2
            # small-tree density contribution to RDNEXT(j+1)/BANEXT(j+1) (regent.f:1437-1458). Skip last subcycle
            # and D≥3in (large trees handled by the pre-loop projection); needs H2>4.5 for the power-form DBH.
            if j < nper && d < 3f0 && h2 > 4.5f0
                relh = abs(avh - 4.5f0) < 0.01f0 ? 0f0 : clamp((h1 - 4.5f0)/(avh - 4.5f0), 0f0, 1f0)
                dadj = delmax*relh*relh - 2f0*delmax*relh + 0.65f0
                d1pf = h1 > 4.5f0 ? BC_RG_HHT1[sp]*fpow(h1 - 4.5f0, BC_RG_HHT2[sp]) + dadj : BC_RG_DIAM[sp] + dadj
                d2pf = BC_RG_HHT1[sp]*fpow(h2 - 4.5f0, BC_RG_HHT2[sp]) + dadj
                xrdgro = active_multiplier(s.control, :regd, sp, current_cycle_year(s))
                dgj = (d2pf - d1pf) * xrdgro; dgj < 0f0 && (dgj = 0f0)
                d2 = d + dgj
                c1 = bc_tree_ccf(sp, d1pf, pr); c2 = bc_tree_ccf(sp, d2, pr)
                b1 = 0.005454154f0 * d * d
                rdnext[j+1] += Float32(ky) * (c2 - c1) / 10f0 * decay
                banext[j+1] += (0.005454154f0*d2*d2 - b1) * pr * decay
            end
        end
    end
    # final: HTGR1 + ZZRAN + blend + DBH-dub (regent.f:1490-1690), species order for RNG determinism.
    # TRIPLING (regent.f:1682-1685 L-loop): while tripling (stash≠nothing), REGENT draws a FRESH ZZRAN per
    # tripled record — central→trees, upper(l=1)/lower(l=2)→the stash (htgU/htgL + dgU/dgL for DBH-dub) —
    # so the 3 copies get 3 DIFFERENT random heights AND the RNG stream advances 3×/tree like live FVS. Without
    # it BC drew 1 ZZRAN/tree (27 vs live's 81 on all_BC_essf) ⇒ copies identical + every downstream draw desynced.
    nrec = stash !== nothing ? 3 : 1
    scale_rg = htg_period(s.variant) / fint                 # regent.f:1142 SCALE=YR/FINT
    order = species_major_order(s)   # IND1: SPESRT lineage order within a species (post-TRIPLE copy1, original, copy2)
    @inbounds for i in order
        sp = Int(t.species[i]); d = t.dbh[i]
        xmn = v2 ? BC_RG_V2_XMIN[sp] : BC_RG_XMIN[sp]; xmx = v2 ? BC_RG_V2_XMAX[sp] : BC_RG_XMAX[sp]
        (d >= xmx || t.tpa[i] <= 0f0 || (!v2 && ip[sp] < 1)) && continue
        h = t.height[i]
        htgr1 = wk3[i] - h; htgr1 < 0f0 && (htgr1 = 0f0)
        htg_large = t.ht_growth[i]      # large-tree HTG (height_growth!) for the XWT blend — read before l=0 overwrites
        cap = s.control.sp_size_cap[sp, 4]
        # regent.f:1508/1638: D starts as DBH(I) and, once a pass takes the D<3in branch, the crown block re-reads
        # D=DBH(K) — so the NEXT copy's XWT (:1556), its D<3 test (:1576) and its DDS use the previous pass's DBH(K):
        # the record's own DBH reset to 0.1+DIAM·.01+HK·.001 when HK<4.5 (:1593), else the DBH that DGDRIV wrote into
        # the copy slot (dgdriv.f:271/278 DBH(ITRIPU/L)=DBH(I)). MEASURED FVSbc_dbfix SkyRanch-Control cycle 1: a 6-cm
        # PL at 0.65 m (record 923) — copy 2823 gets XWT 0 after the original's reset, HTG 2.846 ft, keeps DBH 6 cm.
        dcur = d
        for l in 0:(nrec - 1)
            xwt = dcur <= xmn ? 0f0 : (dcur - xmn)/(xmx - xmn)
            small_d = dcur < 3.0f0
            zzran = 0f0
            if dgsd >= 1f0
                while true
                    zzran = bachlo(s.rng, 0f0, 1f0)
                    (zzran <= 1f0 && zzran >= -1.5f0) && break
                end
            end
            # V2: MULTIPLICATIVE error HTGR=HTGR1·exp(ZZRAN·HSIGMA) (regent.f:1544); V3: additive + ST_COEF.SD
            htgr = v2 ? htgr1 * fexp(zzran * BC_RG_V2_HSIGMA) :
                        max(0f0, htgr1 + zzran * BC_STCOEF[ip[sp]].SD)
            htg = htgr*(1f0 - xwt) + xwt*htg_large
            (h + htg > cap) && (htg = max(cap - h, 0.1f0))
            dg = 0f0; dbhk = -1f0
            if small_d                                            # DBH dub (regent.f:1571-1629)
                relh = (h - 4.5f0)/(avh - 4.5f0); relh = clamp(relh, 0f0, 1f0)
                dadj = delmax*relh*relh - 2f0*delmax*relh + 0.65f0
                d1 = bc_st_dbh(sp, h, dadj)                        # D1 (regent.f:1588-1589)
                hk = h + htg
                xrdgro = active_multiplier(s.control, :regd, sp, current_cycle_year(s))
                if hk < 4.5f0
                    # regent.f:1593-1595 DBH(K)=0.1+DIAM(ISPC)*.01+HK*0.001; DG(K)=0 — set directly, and the
                    # :1627 min-diameter floor sits in the ELSE branch, so it does not apply here.
                    dbhk = 0.1f0 + BC_RG_DIAM[sp] * 0.01f0 + hk * 0.001f0
                    dg = 0f0
                else
                    dk = BC_RG_HHT1[sp]*fpow(hk - 4.5f0, BC_RG_HHT2[sp]) + dadj   # regent.f:1597 (powf)
                    dk < BC_RG_DIAM[sp] && (dk = BC_RG_DIAM[sp])   # 1600 DIAM floor on DK
                    dk += hk * 0.001f0                             # 1601
                    # regent.f:1615-1625: DGK=(DK−D1)·XRDGRO ≥0; BARK=BRATIO; DG=DGK·BARK; DDS=DG·(2·BARK·D+DG)·SCALE;
                    # DG=SQRT((D·BARK)**2.0+DDS)−BARK·D with SCALE=YR/FINT (regent.f:1142): the small-tree DG on the YR=10
                    # basis that GRADD rescales (gradd.f:79-90), like the other _gradd_rescale variants.
                    dgk = (dk - d1) * xrdgro; dgk < 0f0 && (dgk = 0f0)
                    bark = bc_bratio(sp)
                    dg = dgk * bark
                    dds = dg * (2f0 * bark * dcur + dg) * scale_rg
                    dg = sqrt(fpow(dcur * bark, 2.0f0) + dds) - bark * dcur
                end
                dbhk < 0f0 && (d + dg) < BC_RG_DIAM[sp] && (dg = BC_RG_DIAM[sp] - d) # MIN-DIAMETER floor on DBH(K) (regent.f:1627-1629)
                dcur = dbhk >= 0f0 ? dbhk : d                     # regent.f:1638 D=DBH(K) (K's slot holds DBH(I) unless reset)
            end
            if l == 0
                t.ht_growth[i] = htg
                small_d && (t.diam_growth[i] = dg)
                dbhk >= 0f0 && (t.dbh[i] = dbhk)                   # the record's own DBH(K) (regent.f:1594)
            elseif l == 1
                stash.htgU[i] = htg; stash.is_small[i] = true
                small_d && (stash.dgU[i] = dg)                     # D<3in ⇒ regent DBH-dub; else keep driver DDS dgU
                small_d && (stash.dbhU[i] = dbhk >= 0f0 ? dbhk : d)
            else
                stash.htgL[i] = htg
                small_d && (stash.dgL[i] = dg)
                small_d && (stash.dbhL[i] = dbhk >= 0f0 ? dbhk : d)
            end
        end
    end
    return s
end

# --- V2 (LV2ATV) small-tree coefficients (regent.f DATA) — verified via instrument (RHCON(14)=1.090813). ---
# RHCON(sp) = REGCH + RHSC(sp) + RHHAB(IRHHAB,sp); REGCH = RHGL(IGL=2) + (RSAB0+RSAB1·cosA+RSAB2·sinA)·slope.
# IRHHAB = regent.f's OWN MAPHAB(ITYPE=4,sp) [4th distinct subsystem MAPHAB]. See v2_smalltree_data.txt.
const BC_RG_V2_RHLH  = Float32[0.4214,0.2716,0.3907,0.3487,0.3417,0.2354,0.5843,0.2827,0.374,0.4485,0.2354,0.2354,0.2354,0.3907,0.2354]
const BC_RG_V2_RHCCF = Float32[-0.00591,-0.00654,-0.00591,-0.00391,-0.00391,-0.00391,-0.00654,-0.00391,-0.00391,-0.00654,-0.00391,-0.00391,-0.00391,-0.00591,-0.00391]
const BC_RG_V2_RHBAL = Float32[-0.37199,-0.41532,-0.40043,-0.25355,-0.34693,-0.12013,-0.24172,-0.253,-0.22957,-0.47299,-0.25349,-0.25349,-0.25349,-0.40043,-0.25349]
const BC_RG_V2_RHSC  = Float32[1.47,1.6204,1.4932,0.9981,1.0202,0.8953,1.2336,1.0964,1.0667,1.7311,0.8953,0.8953,0.8953,1.4932,0.8953]
const BC_RG_V2_IRHHAB = Int[3,3,4,3,1,1,5,1,4,3,3,3,3,4,3]     # regent MAPHAB(ITYPE=4, sp)
const BC_RG_V2_RHHAB = ([  # [sp][1..6]
    Float32[-0.2146,-0.0941,-0.3141,0,0,0], Float32[-0.2146,-0.0941,-0.3296,0,0,0],
    Float32[-0.2146,-0.0941,-0.5401,-0.3948,0,0], Float32[-0.2146,-0.0941,-0.2776,0,0,0],
    Float32[-0.2146,-0.0941,0,0,0,0], Float32[-0.2146,-0.0941,0,0,0,0],
    Float32[-0.2146,-0.0941,-0.2484,-0.5134,-0.3495,0], Float32[-0.2146,-0.0941,-0.3431,0,0,0],
    Float32[-0.2146,-0.0941,-0.4916,-0.3582,0,0], Float32[-0.2146,-0.0941,-0.4345,0,0,0],
    Float32[-0.2146,-0.0941,-0.3738,0,0,0], Float32[-0.2146,-0.0941,-0.3738,0,0,0],
    Float32[-0.2146,-0.0941,-0.3738,0,0,0], Float32[-0.2146,-0.0941,-0.5401,-0.3948,0,0],
    Float32[-0.2146,-0.0941,-0.3738,0,0,0],
]...,)
const BC_RG_V2_XMAX = Float32[10,10,10,10,10,10,5,10,10,10,10,10,10,10,10]
const BC_RG_V2_XMIN = Float32[2,2,2,2,2,2,1,2,2,2,2,2,2,2,2]
const BC_RG_V2_RHGL = Float32[-0.2785,-0.0480,0.0]   # RHGL(IGL); IGL=2 (grinit.f:204)
const BC_RG_V2_RSAB = Float32[-0.10987,0.22157,-0.12432]   # RSAB0/1/2 (aspect/slope)
const BC_RG_V2_HSIGMA = 0.59f0

"""V2 stand-level REGCH = RHGL(2) + (RSAB0 + RSAB1·cosA + RSAB2·sinA)·slope (regent.f:2020)."""
bc_v2_regch(aspect::Real, slope::Real) = BC_RG_V2_RHGL[2] +
    (BC_RG_V2_RSAB[1] + BC_RG_V2_RSAB[2]*fcos(Float32(aspect)) + BC_RG_V2_RSAB[3]*fsin(Float32(aspect))) * Float32(slope)

"""V2 per-species RHCON = REGCH + RHSC(sp) + RHHAB(IRHHAB,sp) (regent.f:2024, NI case, no RCOR2)."""
bc_v2_rhcon(sp::Integer, regch::Real) =
    Float32(regch) + BC_RG_V2_RHSC[sp] + BC_RG_V2_RHHAB[sp][BC_RG_V2_IRHHAB[sp]]

"""
    bc_esgent!(s, nstart; fint, atba, atccf, atavh, ba_pre, relden_pre, pccf_pre)

strp/esgent.f → canada/bc/regent.f REGENT(LESTB=.TRUE., ITRNIN=nstart+1): grow the records ESTAB just added for the
rest of their birth cycle. NTYR = FINT−5 (LSKIPH when ≤0, regent.f:1107-1108); the subcycle densities interpolate from
TEMBA/TEMCCF (= ATBA/ATCCF, else the gradd.f:192 DENSE BA/RELDEN) toward that DENSE's BA/RELDEN starting 5 years in
(regent.f:1222-1233); DO 13 dubs each new record's crown 0.89722−0.0000461·PCCF + 0.07985·RAN (storage order,
regent.f:1213-1228); DELMAX from TEMCCF/TEMAHT (:1235-1243); the subcycle height and DO 30 ZZRAN/DBH passes run over
the SPESRT species-major order with XWT=0 and DBH=DG=DK (:1606-1609), no tripling, no DUBSCR. esgent.f then scales
HTG by WK4=HTIMLT, rescales DBH when WK4<1 and clamps HT at HHTMAX.
"""
function bc_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0, atba::Float32 = 0f0, atccf::Float32 = 0f0,
                    atavh::Float32 = 0f0, ba_pre::Float32 = 0f0, relden_pre::Float32 = 0f0,
                    pccf_pre::Vector{Float32} = Float32[])
    t = s.trees
    nstart >= t.n && return s
    zone, series = bc_stand_zone(s)
    ip, rhcon = bc_regcons!(s)
    v2 = bc_lv2atv(zone)
    nsp = nspecies(BritishColumbia())
    con_v2 = zeros(Float32, nsp)
    if v2
        regch_v2 = bc_v2_regch(s.plot.aspect, s.plot.slope)
        @inbounds for sp in 1:nsp; con_v2[sp] = bc_v2_rhcon(sp, regch_v2) + s.calib.htg_cor_small[sp]; end
    end
    regyr = 5.0f0; dgsd = s.control.dg_sd; yr_now = current_cycle_year(s)
    # regent.f:1106-1121 subcycles (NTYR = INT(FINT) − 5 under LESTB); ITOT ends as the last KPER.
    ntyr = trunc(Int, fint) - 5
    lskiph = ntyr <= 0
    iyr = Int(regyr)
    nper = lskiph ? 0 : ntyr ÷ iyr + (ntyr % iyr != 0 ? 1 : 0)
    kper = zeros(Int, max(nper, 1)); itot = ntyr; nn = nper
    @inbounds for i in 1:nper
        if nn == 1; break; end
        kper[i] = itot ÷ nn; itot -= kper[i]; nn -= 1
    end
    nper > 0 && (kper[nper] = itot)
    # regent.f:1125-1130 TEMBA/TEMCCF/TEMAHT; :1222-1233 statement 8-10 interpolation.
    ba = ba_pre; relden = relden_pre
    temba = atba > 0f0 ? atba : ba
    temccf = atccf > 0f0 ? atccf : relden
    temaht = atavh
    bayr = 0f0; ccfyr = 0f0
    if !lskiph
        bayr = (ba - temba) / Float32(itot); ccfyr = (relden - temccf) / Float32(itot)
    end
    banext = zeros(Float32, max(nper, 1)); rdnext = zeros(Float32, max(nper, 1))
    nyr = 5
    @inbounds for j in 1:nper
        rdnext[j] = temccf + Float32(nyr) * ccfyr
        banext[j] = temba + Float32(nyr) * bayr
        nyr += kper[j]
    end
    # DO 13 (storage order): crown dub of the new records; WK3=HT, WK5=DBH.
    @inbounds for i in (nstart + 1):t.n
        ipt = Int(t.plot_id[i])
        pccf = (1 <= ipt <= length(pccf_pre)) ? pccf_pre[ipt] : s.density.point_ccf[ipt]
        cr = 0.89722f0 - 0.0000461f0 * pccf
        ran = 0f0
        while true; ran = bachlo(s.rng, 0f0, 1f0); (-1f0 <= ran <= 1f0) && break; end
        cr = cr + 0.07985f0 * ran
        cr > 0.90f0 && (cr = 0.90f0); cr < 0.20f0 && (cr = 0.20f0)
        t.crown_pct[i] = Int32(trunc(Int, cr * 100f0 + 0.5f0))
    end
    wk3 = Float32[t.height[i] for i in 1:t.n]
    delmax = (temaht / 36f0) * (0.01232f0 * temccf - 1.75f0); delmax > 0f0 && (delmax = 0f0)
    order = sort(collect((nstart + 1):t.n); by = i -> (Int(t.species[i]), i))   # esgent.f SPESRT → IND1 (new records)
    ky = 0
    @inbounds for j in 1:nper
        baj = banext[j]; rdj = rdnext[j]; kpj = Float32(kper[j]); ky += kper[j]
        decay = fpowi(0.985f0, ky)
        for i in order
            sp = Int(t.species[i]); d = t.dbh[i]
            (!v2 && ip[sp] < 1) && continue                                   # regent.f:1300 STG%FIT
            xmx = v2 ? BC_RG_V2_XMAX[sp] : BC_RG_XMAX[sp]
            d >= xmx && continue
            pr = t.tpa[i]; h1 = wk3[i]
            bal = (1f0 - t.crown_ratio[i] / 100f0) * baj                       # PCT=0 on a new record (estab.f:640)
            incr = v2 ? fexp(max(-40f0, con_v2[sp] + BC_RG_V2_RHLH[sp]*flog(h1) +
                                BC_RG_V2_RHCCF[sp]*rdj + BC_RG_V2_RHBAL[sp]*bal)) :
                        bc_v3_sthg(sp, ip[sp], h1, bal, rdj, rhcon[sp], s.plot.aspect, s.plot.slope)
            xrhgro = active_multiplier(s.control, :regh, sp, yr_now)
            h2 = h1 + incr * (kpj / regyr) * xrhgro
            wk3[i] = h2
            if j < nper && d < 3f0 && h2 > 4.5f0                               # regent.f:1437-1458
                relh = abs(temaht - 4.5f0) < 0.01f0 ? 0f0 : clamp((h1 - 4.5f0)/(temaht - 4.5f0), 0f0, 1f0)
                dadj = delmax*relh*relh - 2f0*delmax*relh + 0.65f0
                d1pf = h1 > 4.5f0 ? BC_RG_HHT1[sp]*fpow(h1 - 4.5f0, BC_RG_HHT2[sp]) + dadj : BC_RG_DIAM[sp] + dadj
                d2pf = BC_RG_HHT1[sp]*fpow(h2 - 4.5f0, BC_RG_HHT2[sp]) + dadj
                xrdgro = active_multiplier(s.control, :regd, sp, yr_now)
                dgj = (d2pf - d1pf) * xrdgro; dgj < 0f0 && (dgj = 0f0)
                d2 = d + dgj
                c1 = bc_tree_ccf(sp, d1pf, pr); c2 = bc_tree_ccf(sp, d2, pr)
                rdnext[j+1] += Float32(ky) * (c2 - c1) / 10f0 * decay
                banext[j+1] += (0.005454154f0*d2*d2 - 0.005454154f0*d*d) * pr * decay
            end
        end
    end
    # DO 30: HTGR1 + ZZRAN, XWT=0, size cap, DBH/DG (LESTB: DBH=DG=DK), DGBND.
    @inbounds for i in order
        sp = Int(t.species[i]); d = t.dbh[i]
        (!v2 && ip[sp] < 1) && continue
        xmx = v2 ? BC_RG_V2_XMAX[sp] : BC_RG_XMAX[sp]
        d >= xmx && continue
        h = t.height[i]
        local htg::Float32
        if lskiph
            htg = 0f0
        else
            htgr1 = wk3[i] - h; htgr1 < 0f0 && (htgr1 = 0f0)
            zzran = 0f0
            if dgsd >= 1f0
                while true; zzran = bachlo(s.rng, 0f0, 1f0); (zzran <= 1f0 && zzran >= -1.5f0) && break; end
            end
            htg = v2 ? htgr1 * fexp(zzran * BC_RG_V2_HSIGMA) : max(0f0, htgr1 + zzran * BC_STCOEF[ip[sp]].SD)
            cap = s.control.sp_size_cap[sp, 4]
            (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        end
        t.ht_growth[i] = htg
        if d < 3f0
            hk = h + htg
            if hk < 4.5f0
                t.dbh[i] = 0.1f0 + BC_RG_DIAM[sp] * 0.01f0 + hk * 0.001f0
                t.diam_growth[i] = 0f0
            else
                relh = abs(temaht - 4.5f0) < 0.01f0 ? 0f0 : clamp((h - 4.5f0)/(temaht - 4.5f0), 0f0, 1f0)
                dadj = delmax*relh*relh - 2f0*delmax*relh + 0.65f0
                dk = BC_RG_HHT1[sp]*fpow(hk - 4.5f0, BC_RG_HHT2[sp]) + dadj
                dk < BC_RG_DIAM[sp] && (dk = BC_RG_DIAM[sp])
                dk += hk * 0.001f0
                t.dbh[i] = dk; t.diam_growth[i] = dk
                (t.dbh[i] + t.diam_growth[i]) < BC_RG_DIAM[sp] && (t.diam_growth[i] = BC_RG_DIAM[sp] - t.dbh[i])
            end
            # DGBND (dgbnd.f): DBH+DG past SIZCAP(1) with SIZCAP(3)<1.5 ⇒ DG = SIZCAP(1)−DBH, ≥0.01
            if t.dbh[i] + t.diam_growth[i] > s.control.sp_size_cap[sp, 1] && s.control.sp_size_cap[sp, 3] < 1.5f0
                t.diam_growth[i] = max(s.control.sp_size_cap[sp, 1] - t.dbh[i], 0.01f0)
            end
        end
    end
    # esgent.f:50-70 — HTG·WK4, HT += HTG, WK4<1 DBH rescale, HHTMAX clamp (storage order).
    hhtmax = _BC_ES_HHTMAX
    @inbounds for i in (nstart + 1):t.n
        sp = Int(t.species[i])
        htemp = t.height[i] + t.ht_growth[i]
        wk4 = t.htimlt[i]
        t.ht_growth[i] = t.ht_growth[i] * wk4
        t.height[i] = t.height[i] + t.ht_growth[i]
        if wk4 < 1f0
            if t.height[i] < 4.5f0
                t.dbh[i] = 0.1f0 + 0.001f0 * t.height[i]
                t.diam_growth[i] = 0f0
            else
                t.dbh[i] = t.dbh[i] * (t.height[i] / htemp)
                t.diam_growth[i] = t.dbh[i] * (t.height[i] / htemp)
            end
        end
        t.height[i] > hhtmax[sp] && (t.height[i] = hhtmax[sp])
    end
    return s
end
