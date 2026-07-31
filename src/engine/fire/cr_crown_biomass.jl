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

# SPIE groups whose FMCROWW large-tree LIVEWT branches on the height percentile HP<DOMPCT(60):
# ponderosa (13), Douglas-fir (3), western larch (8), Black-Hills PP (25). Others ignore HP.
@inline _cr_crownw_needs_hp(spie::Integer) = spie == 13 || spie == 3 || spie == 8 || spie == 25

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

# ALGSLP for 3 breakpoints (algslp.f): flat below x[1] / above x[3], linear between.
@inline function _cr_algslp3(x::Float32, xv::NTuple{3,Float32}, yv::NTuple{3,Float32})::Float32
    x <= xv[1] && return yv[1]
    x >= xv[3] && return yv[3]
    x <= xv[2] && return yv[1] + (yv[2] - yv[1]) * (x - xv[1]) / (xv[2] - xv[1])
    return yv[2] + (yv[3] - yv[2]) * (x - xv[2]) / (xv[3] - xv[2])
end

# small-tree TOTWT (fmcroww.f:174-255)
@inline function _cr_crownw_small_totwt(spi::Int, h::Float32)::Float32
    spi == 4 && return 0.4284f0 * h                          # grand fir
    spi == 13 && return 0.3451f0 * h                         # ponderosa
    spi == 15 && return 0.3292f0 * h                         # western white pine
    spi == 18 && return exp(-3.932f0 + 2.571f0 * log(h))     # Engelmann spruce
    error("cr_crownw: small-tree TOTWT for SPIE group $spi not ported")
end

# small-tree XV(0),(1),(2) proportions of TOTWT (fmcroww.f:257-360)
@inline function _cr_crownw_small_prop(spi::Int)::NTuple{3,Float32}
    (spi == 4 || spi == 18) && return (0.62f0, 0.26f0, 0.12f0)   # CASE (1,4,6,7,18,24)
    spi == 13 && return (0.57f0, 0.14f0, 0.29f0)                 # CASE (13,25)
    spi == 15 && return (0.52f0, 0.27f0, 0.21f0)                 # CASE (3,11,14,15)
    error("cr_crownw: small-tree proportions for SPIE group $spi not ported")
end

# large-tree LIVEWT/DEADWT/P1-P4/DP1-3 (fmcroww.f:381-1147). DOMPCT=60.
function _cr_crownw_large(spi::Int, d::Float32, r::Float32, c::Float32, hp::Float32)
    dompct = 60f0
    dp1 = 0f0; dp2 = 0f0; dp3 = 1f0     # fmcroww.f:108-110 init (DP3=1)
    p4 = 1f0
    if spi == 4                          # grand fir
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
    end
    error("cr_crownw: large-tree model for SPIE group $spi not ported")
end
