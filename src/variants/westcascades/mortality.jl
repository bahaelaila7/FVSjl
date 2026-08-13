# =============================================================================
# mortality.jl (westcascades) — WC periodic mortality (vwc/morts.f). Chunk 7.
#
# WC mortality is the ORGANON logistic annual-RIP model (NOT the shared Reineke/Hamilton MORTS driver):
# per species a MORTMAP-selected equation (CASE 1-6) gives a 5-yr logit → annual survival^(1/5) → annual
# RIP, times a low-crown CRADJ; sub-3" trees use the Gould-Harrington small-tree RIP; redwood(17)/Oregon-
# white-oak(WO) have their own forms. WKI = P·(1−(1−RIP)^FINT). Then the density self-thin: scale every
# tree's kill by the smallest integer PASS that brings post-mortality SDI < SDIMAX AND BA < 550 ft² (≤100
# passes). XSITE1 = King-converted DF SI (5.21486+0.66486·SITEAR[DF]); XSITE2 = SITEAR[WH].
#
# On wct01 (SDI 184 ≪ SDIMAX 815) the density iteration is inert at PASS=1 ⇒ TPA decline is the base RIP.
# =============================================================================

const WC_MORTMAP = Int[2,2,2,2,2,2,2,4,4,2,4,4,4,4,4,1,6,3,3,3,4,4,4,4,4,4,4,5,4,4,4,4,4,4,4,4,4,1,1]
const WC_MCLASS  = Int[1,2,2,2,2,2,3,2,3,2,4,4,3,3,4,3,1,1,1,1,1,4,3,4,3,5,5,5,5,4,4,5,1,1,3,3,5,1,5]
const WC_MORT_MVALUES = Float32[1.0,1.5,2.25,3.375,5.062]
const WC_MORT_ALPHA   = Float32[-4.4384,0.0053,-0.6001]
const WC_MORT_BETA = Float32[0.247354,0.217481,0.179705,0.205647,0.216823,0.216823,0.282203,0.216823,0.281542,0.170425,0.168227,0.216823,0.216823,0.216823,0.236925,0.163506,0.216823,0.18294,0.17269,0.302866,0.216823,0.216823,0.216823,0.216823,0.216823,0.216823,0.216823,0.216823,0.216823,0.216823,0.216823,0.216823,0.216823,0.216823,0.216823,0.216823,0.216823,0.216823,0.216823]
const WC_BM0 = Float32[-7.60159,-7.60159,-7.60159,-7.60159,-7.60159,-7.60159,-7.60159,-1.922689902,-1.922689902,-7.60159,-1.050000682,-1.050000682,-1.050000682,-1.050000682,-1.050000682,-4.13142,-4.13142,-0.761609,-0.761609,-0.761609,-2.976822456,-2.0,-2.0,-2.0,-4.317549852,-2.0,-2.0,-0.0,-1.050000682,-1.050000682,-1.050000682,-1.050000682,-4.072781265,-3.020345211,-3.020345211,-3.020345211,-2.0,-4.13142,-4.13412]
const WC_BM1 = Float32[-0.200523,-0.200523,-0.200523,-0.200523,-0.200523,-0.200523,-0.200523,-0.13608199,-0.13608199,-0.200523,-0.194363402,-0.194363402,-0.194363402,-0.194363402,-0.194363402,-1.13736,-1.13736,-0.529366,-0.529366,-0.529366,0.0,-0.5,-0.5,-0.5,-0.057696253,-0.5,-0.5,-0.0,-0.194363402,-0.194363402,-0.194363402,-0.194363402,-0.176433475,0.0,0.0,0.0,-0.5,-1.13736,-1.13736]
const WC_BM2 = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.002479863,0.002479863,0.0,0.0038031,0.0038031,0.0038031,0.0038031,0.0038031,0.0,0.0,0.0,0.0,0.0,0.0,0.015,0.015,0.015,0.0,0.015,0.015,0.0,0.0038031,0.0038031,0.0038031,0.0038031,0.0,0.0,0.0,0.0,0.015,0.0,0.0]
const WC_BM3 = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,-3.178123293,-3.178123293,0.0,-3.557300286,-3.557300286,-3.557300286,-3.557300286,-3.557300286,-0.823305,-0.823305,-4.74019,-4.74019,-4.74019,-6.223250962,-3.0,-3.0,-3.0,0.0,-3.0,-3.0,0.0,-3.557300286,-3.557300286,-3.557300286,-3.557300286,-1.729453975,-8.467882343,-8.467882343,-8.467882343,-3.0,-0.823305,-0.823305]
const WC_BM4 = Float32[0.0441333,0.0441333,0.0441333,0.0441333,0.0441333,0.0441333,0.0441333,0.0,0.0,0.0441333,0.003971638,0.003971638,0.003971638,0.003971638,0.003971638,0.0307749,0.0307749,0.0119587,0.0119587,0.0119587,0.0,0.015,0.015,0.015,0.004861355,0.015,0.015,0.0,0.003971638,0.003971638,0.003971638,0.003971638,0.0,0.013966388,0.013966388,0.013966388,0.015,0.0307749,0.0307749]
const WC_BM5 = Float32[0.00063849,0.00063849,0.00063849,0.00063849,0.00063849,0.00063849,0.00063849,0.004684133,0.004684133,0.00063849,0.005573601,0.005573601,0.005573601,0.005573601,0.005573601,0.00991005,0.00991005,0.00756365,0.00756365,0.00756365,0.0,0.01,0.01,0.01,0.00998129,0.01,0.01,0.0,0.005573601,0.005573601,0.005573601,0.005573601,0.012525642,0.009461545,0.009461545,0.009461545,0.01,0.00991005,0.00991005]

# vwc/morts.f per-tree annual RIP. Returns the ANNUAL mortality probability for one tree.
@inline function wc_mort_rip(sp::Int, d::Float32, cr::Float32, bal::Float32, ptbal::Float32,
                             ht::Float32, avh::Float32, ba::Float32, xsite1::Float32, xsite2::Float32)::Float32
    # --- large-tree ORGANON RIP (MORTMAP CASE) ---
    cradj = cr <= 0.17f0 ? 1.0f0 - exp(-(25.0f0 * cr)^2) : 1.0f0
    b0 = WC_BM0[sp]; b1 = WC_BM1[sp]; b2 = WC_BM2[sp]; b3 = WC_BM3[sp]; b4 = WC_BM4[sp]; b5 = WC_BM5[sp]
    local rip::Float32
    mm = WC_MORTMAP[sp]
    if mm == 1
        r = b0 + b1 * sqrt(d) + b3 * cr^0.25f0 + b4 * (xsite1 + 4.5f0) + b5 * bal
        r = 1.0f0 / (1.0f0 + exp(-r)); r = (1.0f0 - r)^0.2f0; rip = 1.0f0 - r * cradj
    elseif mm == 2
        r = b0 + b1 * d + b4 * (xsite1 + 4.5f0) + b5 * (bal / d)
        r = 1.0f0 / (1.0f0 + exp(-r)); r = (1.0f0 - r)^0.2f0; rip = 1.0f0 - r * cradj
    elseif mm == 3
        r = b0 + b1 * d + b2 * d * d + b3 * cr + b4 * (xsite2 + 4.5f0) + b5 * bal
        r = 1.0f0 / (1.0f0 + exp(-r)); r = (1.0f0 - r)^0.2f0; rip = 1.0f0 - r * cradj
    elseif mm == 4
        r = b0 + b1 * d + b2 * d * d + b3 * cr + b4 * (xsite1 + 4.5f0) + b5 * bal
        r = 1.0f0 / (1.0f0 + exp(-r)); r = (1.0f0 - r)^0.2f0; rip = 1.0f0 - r * cradj
    elseif mm == 5                                     # WO — Gould-Harrington (no ^0.2, no CRADJ)
        relht = avh > 0f0 ? ht / avh : 0f0; relht > 1.5f0 && (relht = 1.5f0)
        r = -6.6707f0 + 0.5105f0 * log(5f0 + ba) - 1.3183f0 * relht
        rip = 1.0f0 - 1.0f0 / (1.0f0 + exp(r))
    else                                               # mm == 6 — redwood (Castle 2021)
        r = 2.901447f0 + 0.578694f0 * d - 0.001793f0 * ptbal
        rip = 1.0f0 / (1.0f0 + exp(r))
    end
    # --- small-tree Gould-Harrington override (D<3, not redwood) ---
    if d < 3.0f0 && sp != 17
        relht = avh > 0f0 ? ht / avh : 0f0; relht > 1.5f0 && (relht = 1.5f0)
        hbh = ht >= 4.5f0 ? 4.5f0 : ht
        dbha = d + WC_MORT_BETA[sp] * hbh
        avalue = WC_MORT_MVALUES[WC_MCLASS[sp]]
        r = ptbal * avalue / sqrt(dbha + 1.0f0)
        r = WC_MORT_ALPHA[1] + WC_MORT_ALPHA[2] * r + WC_MORT_ALPHA[3] * relht
        rip = 1.0f0 - 1.0f0 / (1.0f0 + exp(r))
    end
    return rip
end

function mortality!(s::StandState, ::WestCascades; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    sd = s.coef.species
    ba = p.basal_area; avh = p.avg_height
    xsite1 = 5.21486f0 + 0.66486f0 * p.sp_site_index[16]   # WC: Curtis→King DF SI
    xsite2 = p.sp_site_index[19]                            # WH SI
    dbhstage = s.control.dbh_stage                          # WC LZEIDE=.FALSE. ⇒ Reineke/Stage min DBH
    sdimax = stand_sdimax(s)
    # --- DO 20: grown-QMD (DQ10) inputs + per-tree CIOBDS/G for the SDI self-thin iteration ---
    ciobds = @view s.scratch.mort_efftr[1:n]; g1 = @view s.scratch.mort_temwk2[1:n]   # reuse mort work buffers
    @inbounds for i in 1:n
        d = t.dbh[i]; sp = Int(t.species[i])
        bark = wc_bratio(sd, sp, d); g = t.diam_growth[i] / bark
        g1[i] = g; ciobds[i] = 2.0f0 * d * g + g * g
    end
    # AVED (average DBH weighted by tpa, over ALL trees)
    dsum = 0f0; wprob = 0f0
    @inbounds for i in 1:n; wprob += t.tpa[i]; dsum += t.dbh[i] * t.tpa[i]; end
    aved = wprob > 0f0 ? dsum / wprob : 0.0001f0
    # --- per-tree base kill WK2 ---
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    @inbounds for i in 1:n
        pr = t.tpa[i]; pr <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]
        d <= 0.5f0 && (d = 0.5f0)                           # morts.f:305 floor
        cr = Float32(t.crown_pct[i]) * 0.01f0
        bal = (1.0f0 - Float32(t.crown_ratio[i]) / 100.0f0) * ba   # BA in larger trees (PCT percentile)
        ptbal = s.density.point_bal[i]
        rip = wc_mort_rip(sp, d, cr, bal, ptbal, t.height[i], avh, ba, xsite1, xsite2)
        rip < 0.001f0 && (rip = 0.001f0)                    # morts.f:425 floor
        wki = pr * (1.0f0 - (1.0f0 - rip)^fint)             # X=1 (no MORTMULT default)
        wki > pr && (wki = pr)
        sdimax < 5.0f0 && (wki = pr)                        # climate: site can't support trees
        killed[i] = wki
    end
    # --- density self-thin: scale kill by the smallest integer PASS s.t. SDIA<SDIMAX AND BAA<550 ---
    if sdimax >= 5.0f0
        pass = 1
        while pass <= 100
            sd2sqa = 0f0; ta = 0f0
            @inbounds for i in 1:n
                pr = t.tpa[i]; wki = killed[i] * pass; wki > pr && (wki = pr)
                d = t.dbh[i]; d < dbhstage && continue
                sd2sqa += (pr - wki) * (d * d + ciobds[i]); ta += (pr - wki)
            end
            ta <= 0f0 && break
            dq10a = sqrt(sd2sqa / ta); baa = 0.005454154f0 * dq10a * dq10a * ta
            sdia = ta * (dq10a / 10.0f0)^1.605f0
            (sdia < sdimax && baa < 550.0f0) && break
            pass += 1
        end
        pass > 100 && (pass = 100)
        if pass > 1
            @inbounds for i in 1:n
                wki = killed[i] * pass; wki > t.tpa[i] && (wki = t.tpa[i]); killed[i] = wki
            end
        end
    end
    apply_fixmort!(s, killed, n, fint)
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end
