# =============================================================================
# fire/cr_crown_biomass.jl — CR (western) FFE crown biomass FMCROWW (cr/fmcroww.f)
#
# Brown & Johnston 1976 "Debris Prediction System" crown-weight model. Returns XV(0:5) =
# (foliage, woody 0-0.25", 0.25-1", 1-3", 3-6", >6") in POUNDS, per tree.
#
# CR dispatch (fmcrow.f:143-166): SPIW = CR species (1-38); SPIE = ISPMAP(SPIW) = crown group.
#   SELECT CASE(SPIW): CASE(20,21,22,28,38) → FMCROWE (Jenkins, the generic `crown_biomass`);
#   CASE DEFAULT → FMCROWW(SPIE, …) (this file). ISPMAP feeds the SPI group index below.
#
# Algorithm (fmcroww.f): a small-tree model (D ≤ DCTHGH, weight-from-H) and a large-tree model
# (LIVEWT/DEADWT + cumulative live-crown proportions P1-P4 & dead DP1-DP3), blended by SMWGT =
# ALGSLP(D, [DCTLOW,DBHCUT,DCTHGH], [1,0.5,0]). For D > DCTHGH SMWGT=0 (pure large tree, X=1);
# for D < DCTLOW pure small tree. crt01 trees are all > DCTHGH so the large-tree path dominates.
#
# STATUS: crt01-exercised species groups (SPIE 4 grand/white fir, 13 ponderosa, 15 western white
# pine, 18 Engelmann spruce) are ported + validated per-tree vs the live FMCROWW dump
# (/workspace/.crwork/fmcroww_live_dump.txt). Other groups error loudly (doctrine #5) until ported.
# =============================================================================

# ISPMAP (fmcrow.f:100-103): CR species 1..38 → crown group SPIE. (KODFOR 203/207 ⇒ ISPMAP[13]=25.)
const _CR_ISPMAP = Int[
    1, 1, 3, 4, 4, 24, 7, 8, 9, 11,
    11, 12, 13, 14, 15, 16, 18, 18, 18, 41,
    17, 17, 22, 22, 22, 22, 22, 43, 16, 16,
    16, 16, 12, 12, 12, 13, 11, 17,
]

# CR species that use the Jenkins FMCROWE (the generic `crown_biomass`) instead of FMCROWW.
@inline _cr_uses_fmcrowe(spiw::Integer) = spiw == 20 || spiw == 21 || spiw == 22 || spiw == 28 || spiw == 38

# NC (Klamath) crown-biomass group map — nc/fmcrow.f ISPMAP (NC species 1..12 → the crown-equation group
# passed to FMCROWW). NC's fmcroww.f is byte-identical to CR's and dispatches directly on this SPI (no
# internal remap), and nc/fmcrow.f calls FMCROWW for ALL species (none use the eastern FMCROWE). So NC
# routes through `cr_crownw` with this map. Groups {3 DF/OS, 4 WF/RF, 13 PP, 15 SP} are ported; the
# hardwood/cedar groups {10,17,19,20,21} error loudly (doctrine #5) until ported (not in nct01).
const _NC_ISPMAP = Int[3, 15, 3, 4, 10, 20, 21, 17, 4, 13, 17, 19]

# SPIE groups whose FMCROWW large-tree LIVEWT branches on the height percentile HP<DOMPCT(60):
# ponderosa (13), Douglas-fir (3), western larch (8), Black-Hills PP (25). Others ignore HP.
@inline _cr_crownw_needs_hp(spie::Integer) =
    spie == 13 || spie == 3 || spie == 8 || spie == 11 || spie == 25

"""
    cr_hpct_of_height(s, h) -> Float32

The FFE height percentile HP (0-100) of a tree of height `h`, matching FVS FMCROW's
`RDPSRT(HT,desc)+PCTILE(PROB=TPA)` (fmcrow.f:121-130): HP = 100·(ΣTPA of records with height ≤ h)/ΣTPA,
i.e. the reverse-cumulative TPA from the tallest. Tallest record → 100. Computed over the live tree list.
"""
function cr_hpct_of_height(s::StandState, h::Float32)::Float32
    t = s.trees; tot = 0f0; le = 0f0
    @inbounds for i in 1:t.n
        p = t.tpa[i]; p > 0f0 || continue
        tot += p
        t.height[i] <= h && (le += p)
    end
    tot <= 0f0 ? 100f0 : (le / tot) * 100f0
end

# The FMCROWW small-tree breakpoints (fmcroww.f:140-157), by SPIE group.
@inline function _cr_crownw_breaks(spi::Int)
    (spi == 9 || spi == 12 || spi == 16) && return (2f0, 1f0, 3f0)   # DBHCUT, DCTLOW, DCTHGH
    spi == 22 && return (4f0, 3f0, 5f0)
    spi == 11 && return (1f0, 0.5f0, 3f0)
    return (1f0, 0.5f0, 2f0)                                         # DEFAULT
end

"""
    cr_crownw(spie, d, h, itrnc, ic, hp, sg) -> NTuple{6,Float32}

FMCROWW crown weight XV(0:5) in lb for CR species-group `spie`. `d` DBH(in), `h` HT(ft),
`ic` crown ratio %, `hp` height percentile (0-100), `sg` V2T (only used by Gambel oak).
Mirrors cr/fmcroww.f. Errors on species groups not yet ported.
"""
function cr_crownw(spie::Integer, d::Float32, h::Float32, itrnc::Integer, ic::Integer,
                   hp::Float32, sg::Float32)::NTuple{6,Float32}
    spi = Int(spie)
    xv = zeros(Float32, 6)          # XV(0:5) as 1-based xv[1..6]
    tcr = Float32(ic) / 100f0
    c = tcr * h                     # live crown length
    r = h > 0f0 ? (c / h) * 10f0 : 0f0
    (d == 0f0 || h == 0f0 || c == 0f0) && return Tuple(xv)

    dbhcut, dctlow, dcthgh = _cr_crownw_breaks(spi)
    # SMWGT = ALGSLP(D, [DCTLOW,DBHCUT,DCTHGH], [1,0.5,0], 3) — segmented-linear, flat outside.
    smwgt = _cr_algslp3(d, (dctlow, dbhcut, dcthgh), (1f0, 0.5f0, 0f0))

    # DBRx-path groups (bristlecone 9, pinyon 12, juniper 16, Gambel oak 22): DIRECT small-tree XV + metric
    # branch biomass DFOL/DBR1-5 → ×2.2046 assembly (fmcroww.f:319-360 small, 643-1030 large, 1246-1262),
    # NOT the LIVEWT/P1-P4 conifer path.
    (spi == 9 || spi == 12 || spi == 16 || spi == 22) &&
        return _cr_crownw_dbrx(spi, d, h, sg, smwgt, dctlow, dcthgh)

    dd = d
    if d <= dcthgh
        # ---- small-tree model (fmcroww.f:174-364) ----
        totwt = _cr_crownw_small_totwt(spi, h)
        f = _cr_crownw_small_prop(spi)          # (p0,p1,p2) proportions of TOTWT into XV(0),(1),(2)
        xv[1] = f[1] * totwt; xv[2] = f[2] * totwt; xv[3] = f[3] * totwt
        for j in 1:6; xv[j] = max(0f0, xv[j] * smwgt); end
        d < dctlow && return Tuple(xv)          # pure small tree
    end

    # ---- large-tree model (fmcroww.f:299-1264) ----
    dd = dd > 40f0 ? 40f0 : dd
    livewt, deadwt, p1, p2, p3, p4, dp1, dp2, dp3 = _cr_crownw_large(spi, dd, r, c, hp)
    xx = 1f0 - smwgt
    # clamp (fmcroww.f:1164-1187)
    livewt = max(0f0, livewt); deadwt = max(0f0, deadwt)
    p1 = clamp(p1, 0f0, 1f0); p2 = clamp(p2, 0f0, 1f0); p3 = clamp(p3, 0f0, 1f0); p4 = clamp(p4, 0f0, 1f0)
    dp1 = clamp(dp1, 0f0, 1f0); dp2 = clamp(dp2, 0f0, 1f0); dp3 = clamp(dp3, 0f0, 1f0)
    p2 = max(p2, p1); p3 = max(p3, p2); p4 = max(p4, p3)
    dp2 = max(dp2, dp1); dp3 = max(dp3, dp2)
    # assembly (CASE DEFAULT, fmcroww.f:1202-1208)
    xv[1] += (livewt * p1) * xx
    xv[2] += (livewt * (p2 - p1) + deadwt * dp1) * xx
    xv[3] += (livewt * (p3 - p2) + deadwt * (dp2 - dp1)) * xx
    xv[4] += (livewt * (p4 - p3) + deadwt * (dp3 - dp2)) * xx
    xv[5] += (livewt * (1f0 - p4) + deadwt * (1f0 - dp3)) * xx
    for j in 1:6; xv[j] = max(0f0, xv[j]); end
    return Tuple(xv)
end

# DBRx-path crown weight (pinyon/juniper/bristlecone/oak): direct small-tree XV + metric branch biomass.
# Small-tree XV (fmcroww.f:319-360) is in lb (already), ×SMWGT; large-tree DFOL/DBR1-5 (kg) ×2.2046×(1-SMWGT).
function _cr_crownw_dbrx(spi::Int, d::Float32, h::Float32, sg::Float32,
                         smwgt::Float32, dctlow::Float32, dcthgh::Float32)::NTuple{6,Float32}
    xv = zeros(Float32, 6)
    if d <= dcthgh
        # direct small-tree XV (fmcroww.f:319-352), in lb
        if spi == 9                        # bristlecone: X=D/2
            x = d / 2f0
            xv[1] = 8.874f0 * x; xv[2] = 3.428f0 * x; xv[3] = 4.429f0 * x; xv[4] = 0.590f0 * x
        elseif spi == 12                   # pinyon: X=D/2
            x = d / 2f0
            xv[1] = 3.177f0 * x; xv[2] = 0.977f0 * x; xv[3] = 1.084f0 * x; xv[4] = 0.079f0 * x
        elseif spi == 16                   # juniper: X=D/2
            x = d / 2f0
            xv[1] = 0.256f0 * x; xv[2] = 0.126f0 * x; xv[3] = 1.388f0 * x; xv[4] = 0.454f0 * x
        elseif spi == 22                   # Gambel oak: X=D/4
            x = d / 4f0
            xv[1] = 3.2119f0 * x; xv[2] = 1.606f0 * x; xv[3] = 15.8115f0 * x; xv[4] = 7.7877f0 * x
        else
            error("cr_crownw: DBRx small-tree for SPIE group $spi not ported")
        end
        for j in 1:6; xv[j] = max(0f0, xv[j] * smwgt); end
        d < dctlow && return Tuple(xv)
    end
    dd = d > 40f0 ? 40f0 : d
    dfol, dbr1, dbr2, dbr3, dbr4, dbr5 = _cr_dbrx_large(spi, dd, h, sg)
    xx = 1f0 - smwgt
    cf = spi == 22 ? 1f0 : 2.2046f0        # oak DBRx already in lb; pinyon/juniper metric kg→lb
    xv[1] += max(0f0, dfol * cf) * xx
    xv[2] += max(0f0, dbr1 * cf) * xx
    xv[3] += max(0f0, dbr2 * cf) * xx
    xv[4] += max(0f0, dbr3 * cf) * xx
    xv[5] += max(0f0, dbr4 * cf) * xx
    xv[6] += max(0f0, dbr5 * cf) * xx
    for j in 1:6; xv[j] = max(0f0, xv[j]); end
    return Tuple(xv)
end

# metric branch biomass DFOL/DBR1-5 (kg, except oak in lb) for the DBRx groups (fmcroww.f:643-1030).
function _cr_dbrx_large(spi::Int, d::Float32, h::Float32, sg::Float32)
    if spi == 12 || spi == 9               # pinyon / bristlecone — same eqs, different DRC (fmcroww.f:653-654)
        drc = spi == 9 ? 2.54f0 * (d + 1.9410f0) / 1.0222f0 : d * 2.54f0
        l = log10(drc)
        dfol = 10f0^(-0.946f0 + 1.565f0 * l)
        base1 = 10f0^(-1.613f0 + 2.088f0 * l)
        dbr1 = base1 * 0.33f0
        dbr2 = base1 * 0.67f0
        dbr3 = 10f0^(-2.971f0 + 3.007f0 * l) * 0.25f0
        dbr4 = 0f0; dbr5 = 0f0
        dbr1 += 10f0^(-1.873f0 + 1.675f0 * l)                 # old twigs → <1/4"
        x = dbr2 + dbr3                                        # dead branches redistributed over 0.25-1.5"
        temp = 10f0^(-5.400f0 + 4.470f0 * l)
        if x > 1f-6
            dbr2 += temp * dbr2 / x; dbr3 += temp * dbr3 / x
        else
            dbr2 = temp
        end
        return (dfol, dbr1, dbr2, dbr3, dbr4, dbr5)
    elseif spi == 16                       # western juniper (fmcroww.f:871-891): DRC = D·2.54
        drc = d * 2.54f0
        l = log10(drc)
        base = 10f0^(-1.737f0 + 1.382f0 * l)
        dfol = base * 0.67f0
        dbr1 = base * 0.33f0
        dbr2 = 10f0^(-1.476f0 + 1.787f0 * l)
        dbr3 = 10f0^(-1.356f0 + 1.782f0 * l) * 0.25f0
        dbr4 = 0f0; dbr5 = 0f0
        x = dbr2 + dbr3
        temp = 10f0^(-3.543f0 + 2.774f0 * l)
        if x > 1f-6
            dbr2 += temp * dbr2 / x; dbr3 += temp * dbr3 / x
        else
            dbr2 = temp
        end
        return (dfol, dbr1, dbr2, dbr3, dbr4, dbr5)
    elseif spi == 22                       # Gambel oak (Chojnacky 1992, fmcroww.f:995-1030): DRC=D, in lb
        drc = d
        x = drc * drc * h / 1000f0
        x0 = 7.1046f0
        v = x <= x0 ? -0.0534f0 + 2.3077f0 * x + 0.0467f0 * x * x :
            -0.0534f0 + 2.3077f0 * x + 0.0467f0 * 3f0 * (x0 * x0 - x0 * x0 * x0 / x)
        v <= 0.01f0 && (v = 0.01f0)
        v *= sg * 2000f0                                       # SG = raw V2T (lb/ft³) → wt of plant+branches >1.5"
        lv = log10(v)
        dfol = 10f0^(-0.5655f0 + 0.8382f0 * lv - 0.0094f0 * h)
        dbr3 = 10f0^(0.3036f0 + 0.7752f0 * lv - 0.0049f0 * h)
        dbr1 = dfol * 0.5f0
        dbr2 = max(0f0, (dbr3 - dbr1) * 0.67f0)
        dbr3 = max(0f0, (dbr3 - dbr1) * 0.33f0)
        return (dfol, dbr1, dbr2, dbr3, 0f0, 0f0)
    end
    error("cr_crownw: DBRx large-tree for SPIE group $spi not ported")
end

# ALGSLP for 3 breakpoints (algslp.f): flat below x[1] / above x[3], linear between.
@inline function _cr_algslp3(x::Float32, xv::NTuple{3,Float32}, yv::NTuple{3,Float32})::Float32
    x <= xv[1] && return yv[1]
    x >= xv[3] && return yv[3]
    x <= xv[2] && return yv[1] + (yv[2] - yv[1]) * (x - xv[1]) / (xv[2] - xv[1])
    return yv[2] + (yv[3] - yv[2]) * (x - xv[2]) / (xv[3] - xv[2])
end

# small-tree TOTWT (fmcroww.f:174-255)
@inline function _cr_crownw_small_totwt(spi::Int, h::Float32)::Float32
    spi == 1 && return exp(-3.335f0 + 2.303f0 * log(h))      # subalpine/corkbark fir
    spi == 3 && return exp(-4.212f0 + 2.7168f0 * log(h))     # Douglas-fir
    spi == 4 && return 0.4284f0 * h                          # grand fir
    spi == 7 && return 0.04833f0 * h * h                      # western redcedar (CASE 7,19,20)
    spi == 8 && return 0.1128f0 * h + 0.00813f0 * h * h      # western larch
    spi == 11 && return 0.03111f0 * h * h                     # lodgepole pine
    spi == 13 && return 0.3451f0 * h                         # ponderosa
    spi == 14 && return 0.070f0 + 0.02446f0 * h * h          # whitebark pine
    spi == 15 && return 0.3292f0 * h                         # western white pine
    spi == 17 && return 0.81135f0 * h / 5.1213f0             # tanoak / CA black oak
    spi == 18 && return exp(-3.932f0 + 2.571f0 * log(h))     # Engelmann spruce
    spi == 24 && return exp(-5.126f0 + 2.563f0 * log(h))     # mountain hemlock (CASE 6,24)
    error("cr_crownw: small-tree TOTWT for SPIE group $spi not ported")
end

# small-tree XV(0),(1),(2) proportions of TOTWT (fmcroww.f:257-360)
@inline function _cr_crownw_small_prop(spi::Int)::NTuple{3,Float32}
    (spi == 1 || spi == 4 || spi == 7 || spi == 18 || spi == 24) &&
        return (0.62f0, 0.26f0, 0.12f0)                                            # CASE (1,4,6,7,18,24)
    spi == 13 && return (0.57f0, 0.14f0, 0.29f0)                                   # CASE (13,25)
    (spi == 3 || spi == 11 || spi == 14 || spi == 15) &&
        return (0.52f0, 0.27f0, 0.21f0)                                            # CASE (3,11,14,15)
    spi == 8 && return (0.40f0, 0.42f0, 0.18f0)                                    # western larch
    spi == 17 && return (0.38f0, 0.32f0, 0.30f0)                                   # CASE (17,21) tanoak/oak
    error("cr_crownw: small-tree proportions for SPIE group $spi not ported")
end

# large-tree LIVEWT/DEADWT/P1-P4/DP1-3 (fmcroww.f:381-1147). DOMPCT=60.
function _cr_crownw_large(spi::Int, d::Float32, r::Float32, c::Float32, hp::Float32)
    dompct = 60f0
    dp1 = 0f0; dp2 = 0f0; dp3 = 1f0     # fmcroww.f:108-110 init (DP3=1)
    p4 = 1f0
    if spi == 1                          # subalpine / corkbark fir
        livewt = 0.1862f0 * d * d * r + 1.066f0
        deadwt = d <= 16f0 ? exp(4.0365f0 * log(d) - 6.5431f0) : 0.31f0 * livewt
        p1 = 0.5966f0 * exp(-0.04247f0 * d)
        p2 = 0.8643f0 * exp(-0.03733f0 * d)
        p3 = d <= 2.9f0 ? 1f0 : 1.0221f0 - 0.01083f0 * d
        dp1 = d < 1.5f0 ? 1f0 : 1.2105f0 * d^(-0.565f0)
        dp2 = 1f0
        return (livewt, deadwt, p1, p2, p3, p4, dp1, dp2, dp3)
    elseif spi == 3                      # Douglas-fir
        if hp < dompct
            livewt = exp(0.1508f0 + 1.8621f0 * log(d))
            deadwt = exp(-1.928f0 + 2.353f0 * log(d))
        else
            livewt = d < 17f0 ? exp(1.1368f0 + 1.5819f0 * log(d)) : 1.0237f0 * d * d - 20.74f0
            deadwt = 0.01094f0 * d * d * d
        end
        if d > 36f0
            p1 = 0.227f0; p2 = 0.315f0; p3 = 0.465f0
        else
            p1 = 0.484f0 * exp(-0.02102f0 * d)
            p2 = 0.7289f0 * exp(-0.02332f0 * d)
            p3 = d <= 2.9f0 ? 1f0 : 1.0342f0 - 0.01584f0 * d
        end
        p4 = d <= 14f0 ? 1f0 : 1.0221f0 - 0.001821f0 * d
        dp1 = d < 1.8f0 ? 1f0 : 0.08355f0 + (1.5893f0 / d)
        dp2 = d < 9f0 ? 1f0 : 1.5673f0 * exp(-0.05232f0 * d)
        return (livewt, deadwt, p1, p2, p3, p4, dp1, dp2, dp3)
    elseif spi == 4                      # grand fir
        livewt = exp(1.3094f0 + 1.6076f0 * log(d))
        deadwt = d <= 18f0 ? exp(3.5638f0 * log(d) - 5.3154f0) : 0.38f0 * livewt
        if d > 36f0
            p1 = 0.286f0; p2 = 0.378f0; p3 = 0.488f0
        else
            p1 = 1f0 / (1.5916f0 + 0.05294f0 * d)
            p2 = 1f0 / (1.1495f0 + 0.04165f0 * d)
            p3 = d <= 2.9f0 ? 1f0 : 1.0267f0 - 0.01495f0 * d
        end
        dp1 = d < 3f0 ? 1f0 : (d > 27f0 ? 0.01f0 : 1.4336f0 * exp(-0.1816f0 * d))
        dp2 = d < 8f0 ? 1f0 : 1.2623f0 * exp(-0.0347f0 * d)
        return (livewt, deadwt, p1, p2, p3, p4, dp1, dp2, dp3)
    elseif spi == 11                     # lodgepole pine
        livewt = 0.02238f0 * d * d * d + 0.1233f0 * d * d * r - 2f0
        deadwt = d <= 10f0 ? (0.026f0 * d - 0.025f0) * livewt : 0.235f0 * livewt
        if hp < dompct
            if d <= 7.5f0
                livewt *= 0.5f0; deadwt *= 0.5f0
            else
                livewt *= 0.6f0; deadwt *= 0.6f0
            end
        end
        p1 = 0.4933f0 - 0.01167f0 * d
        p2 = 0.7767f0 - 0.01464f0 * d
        p3 = d <= 3.9f0 ? 1f0 : 1.0494f0 - 0.01402f0 * d
        if d > 20f0
            dp1 = 0.139f0; dp2 = 0.226f0
        else
            dp1 = d < 1.5f0 ? 1f0 : 1.3527f0 * d^(-0.7585f0)
            dp2 = d < 9f0 ? 1f0 : 2.7979f0 * exp(-0.1257f0 * d)
        end
        return (livewt, deadwt, p1, p2, p3, p4, dp1, dp2, dp3)
    elseif spi == 13                     # ponderosa pine
        if hp < dompct
            livewt = exp(-0.7572f0 + 2.216f0 * log(d))
            deadwt = exp(-2.5176f0 + 2.51f0 * log(d))
            p1 = 0.6501f0 * exp(-0.1544f0 * d)
            p2 = 0.8435f0 * exp(-0.1665f0 * d)
            p3 = 1.0865f0 * exp(-0.0833f0 * d)
            p4 = 1f0
        else
            livewt = exp(2.2812f0 * log(d) + 1.5098f0 * log(r) - 3.0957f0)
            deadwt = exp(2.8376f0 * log(d) - 3.7398f0)
            p1 = 0.5578f0 * exp(-0.04754f0 * d)
            p2 = d >= 31f0 ? p1 + 0.01f0 : 0.6254f0 * exp(-0.05114f0 * d)
            if d <= 1f0
                p3 = 1f0; p4 = 1f0
            else
                p3 = 0.985f0 * exp(-0.03102f0 * d)
                p4 = d <= 6.5f0 ? 1f0 : 1.083f0 - 0.01306f0 * d
            end
        end
        if d > 30f0
            dp1 = 0.004f0; dp2 = 0.06f0
        else
            dp1 = (1.4114f0 / d) - 0.04345f0
            dp2 = 1.0621f0 - 0.03342f0 * d
        end
        return (livewt, deadwt, p1, p2, p3, p4, dp1, dp2, dp3)
    elseif spi == 15                     # western white pine
        livewt = 0.0947f0 * d * d * r
        deadwt = exp(2.6076f0 * log(d) - 4.397f0)
        p1 = 0.5497f0 * exp(-0.0345f0 * d)
        p2 = 0.9138f0 - 0.0978f0 * sqrt(d)
        p3 = d <= 3.9f0 ? 1f0 : 1.0564f0 * exp(-0.0181f0 * d)
        dp1 = 1.0077f0 * d^(-0.4556f0)
        dp2 = d < 7f0 ? 1f0 : 1.0291f0 - 0.004964f0 * d
        return (livewt, deadwt, p1, p2, p3, p4, dp1, dp2, dp3)
    elseif spi == 18                     # Engelmann spruce
        livewt = exp(1.0404f0 + 1.7096f0 * log(d))
        deadwt = exp(3.6172f0 * log(d) - 6.686f0)
        if d < 40f0
            p1 = 0.5738f0 * exp(-0.0325f0 * d)
            p2 = 0.8519f0 * exp(-0.02811f0 * d)
            p3 = d <= 2.9f0 ? 1f0 : 1.03781f0 - 0.01537f0 * d
        else
            p1 = 0.158f0; p2 = 0.277f0; p3 = 0.423f0
        end
        dp1 = d < 1.8f0 ? 1f0 : 1.4657f0 * d^(-0.6454f0)
        dp2 = d < 10f0 ? 1f0 : 1f0 / (0.847f0 + 0.01678f0 * d)
        return (livewt, deadwt, p1, p2, p3, p4, dp1, dp2, dp3)
    elseif spi == 7                      # western redcedar (CASE 7,20)
        livewt = exp(1.7273f0 * log(d * r) - 2.8086f0)
        deadwt = 0.01063f0 * d * d * d
        p1 = 0.6174f0 * exp(-0.02326f0 * d)
        p2 = 0.7562f0 * exp(-0.02411f0 * d)
        p3 = d <= 2.9f0 ? 1f0 : 1.0602f0 * exp(-0.02226f0 * d)
        dp1 = d < 1.5f0 ? 1f0 : -0.01578f0 + (1.4673f0 / d)
        dp2 = d < 8f0 ? 1f0 : 1.4534f0 * exp(-0.05395f0 * d)
        return (livewt, deadwt, p1, p2, p3, p4, dp1, dp2, dp3)
    elseif spi == 8                      # western larch (no dead crown)
        livewt = exp(0.4373f0 + 1.6786f0 * log(d))
        if hp < dompct
            livewt *= d <= 7.5f0 ? 0.5f0 : 0.6f0
        end
        p1 = 0.3468f0 * exp(-0.04343f0 * d)
        p2 = 0.745f0 * exp(-0.03622f0 * d)
        p3 = d <= 2.9f0 ? 1f0 : 1.05448f0 * exp(-0.0213f0 * d)
        p4 = d <= 11f0 ? 1f0 : 0.9223f0 + 0.7197f0 / d
        deadwt = 0f0; dp1 = 0f0; dp2 = 0f0
        return (livewt, deadwt, p1, p2, p3, p4, dp1, dp2, dp3)
    elseif spi == 14                     # whitebark pine
        livewt = 0.06056f0 * d * d * d + 0.05477f0 * d * d * r + 0.646f0
        deadwt = 0.001713f0 * d * d * c + 0.33f0
        p1 = d > 20f0 ? 0.242f0 : 0.5120f0 * exp(-0.03737f0 * d)
        p2 = d > 20f0 ? 0.268f0 : 0.8644f0 * exp(-0.05854f0 * d)
        p3 = d <= 3.9f0 ? 1f0 : (d > 20f0 ? 0.670f0 : 1.0733f0 * exp(-0.02376f0 * d))
        dp1 = d < 1.4f0 ? 1f0 : 0.268f0 + (1.1733f0 / d)
        dp2 = 1f0
        return (livewt, deadwt, p1, p2, p3, p4, dp1, dp2, dp3)
    elseif spi == 17                     # tanoak / CA black oak (Snell & Little 1983)
        livewt = exp(-0.3169f0 + 2.2774f0 * log(d))
        deadwt = exp(-2.4895f0 + 2.0374f0 * log(d))
        p1 = 1f0 / (1.7936f0 + 0.5952f0 * d^0.7239f0)
        p2 = 1f0 / (0.9940f0 + 0.4229f0 * d^0.6520f0)
        p3 = d < 1.5f0 ? 1f0 : 1f0 / (0.8759f0 + 0.0927f0 * d^0.7843f0)
        dp1 = -0.1424f0 + (0.7684f0 * d^0.25f0) - (0.4730f0 * log(d))
        dp2 = d < 4.1f0 ? 1f0 : exp(-2.810f0 + (4.379f0 * d^0.25f0) - (1.691f0 * d^0.5f0))
        dp3 = d < 7.9f0 ? 1f0 : 1.027f0 - (0.003439f0 * d)
        return (livewt, deadwt, p1, p2, p3, p4, dp1, dp2, dp3)
    elseif spi == 24                     # mountain hemlock (Gholz 1979 weights, hemlock proportions)
        ldm = log(d * 2.54f0)
        livewt = (exp(-3.8169f0 + 1.9756f0 * ldm) + exp(-5.2581f0 + 2.6045f0 * ldm)) * 2.2046f0
        deadwt = exp(-9.9449f0 + 3.2845f0 * ldm) * 2.2046f0
        if d <= 40f0
            p1 = 0.5474f0 * exp(-0.03697f0 * d)
            p2 = 0.8352f0 * exp(-0.03802f0 * d)
            p3 = d <= 2.9f0 ? 1f0 : 1.0781f0 * exp(-0.02735f0 * d)
        else
            p1 = 0.125f0; p2 = 0.183f0; p3 = 0.361f0
        end
        dp1 = d < 4f0 ? 1f0 : (d > 28f0 ? 0.005f0 : 1.9608f0 * exp(-0.2064f0 * d))
        dp2 = d < 12f0 ? 1f0 : 1f0 / (0.2772f0 + 0.06141f0 * d)
        return (livewt, deadwt, p1, p2, p3, p4, dp1, dp2, dp3)
    end
    error("cr_crownw: large-tree model for SPIE group $spi not ported")
end
