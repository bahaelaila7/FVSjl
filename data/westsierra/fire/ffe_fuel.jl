# WS (WestSierra) FFE initial surface-fuel loading — ws/fmcba.f DATA FULIVE/FULIVI/FUINIE/FUINII.
# WS is a CALIFORNIA/R5 variant with the same TOP-TWO-cover-species structure as NC/CA (ws/fmcba.f:519-597):
# the live + initial-dead fuel pools are interpolated from the top two cover-type species (COVCA(1..2),
# weighted by their share of the top-2 basal area COVCAWT), NOT from the single dominant COVTYP the interior
# western variants (CR/IE/EM/CI/TT/UT/BM) use. Tables indexed by WS species 1..43; columns (herb, shrub) for
# the live tables, and the 11 FFE size classes (<.25, .25-1, 1-3, 3-6, 6-12, 12-20, 20-35, 35-50, >50, Litter,
# Duff) for the dead tables. Established (60% cover) vs initiating (10% cover); interpolated by PERCOV over
# XCOV=[10,60]. Values ported verbatim from ws/fmcba.f (from J. Brown / Ottmar 2000). Reuses `_cr_algslp2`.
#
# WS species (ISPC 1..43): 1 SP 2 DF 3 WF 4 GS 5 IC 6 JP 7 RF 8 PP 9 LP 10 WB 11 WP 12 PM 13 SF 14 KP 15 FP
# 16 CP 17 LM 18 MP 19 GP 20 WE 21 GB 22 BD 23 RW 24 MH 25 WJ 26 UJ 27 CJ 28 LO 29 CY 30 BL 31 BO 32 VO 33 IO
# 34 TO 35 GC 36 AS 37 CL 38 MA 39 DG 40 BM 41 MC 42 OS 43 OH

const _WS_FULIVE = Float32[    # established stands (60% cover)  (herb, shrub)  — ws/fmcba.f:143-187
    0.20 0.10;  #  1 sugar pine (use lodgepole pine NI)
    0.20 0.20;  #  2 Douglas-fir
    0.15 0.10;  #  3 white fir
    0.20 0.20;  #  4 giant sequoia (use Douglas-fir)
    0.20 0.20;  #  5 incense cedar (use Douglas-fir)
    0.20 0.10;  #  6 Jeffrey pine (use sugar pine)
    0.15 0.10;  #  7 red fir (use white fir)
    0.20 0.25;  #  8 ponderosa pine
    0.20 0.10;  #  9 lodgepole pine
    0.20 0.10;  # 10 whitebark pine
    0.15 0.10;  # 11 western white pine
    0.04 0.05;  # 12 singleleaf pinyon
    0.15 0.10;  # 13 Pacific silver fir (use white fir)
    0.20 0.10;  # 14 knobcone pine
    0.20 0.10;  # 15 foxtail pine
    0.20 0.10;  # 16 Coulter pine
    0.20 0.10;  # 17 limber pine
    0.20 0.20;  # 18 Monterey pine
    0.23 0.22;  # 19 gray/CA foothill pine
    0.20 0.10;  # 20 washoe pine
    0.04 0.05;  # 21 GB bristlecone pine
    0.20 0.20;  # 22 bigcone Douglas-fir (use Douglas-fir)
    0.20 0.20;  # 23 redwood (use giant sequoia)
    0.15 0.20;  # 24 mountain hemlock
    0.14 0.35;  # 25 western juniper
    0.14 0.35;  # 26 Utah juniper
    0.14 0.35;  # 27 California juniper
    0.23 0.22;  # 28 California live oak
    0.25 0.25;  # 29 canyon live oak
    0.23 0.22;  # 30 blue oak
    0.23 0.22;  # 31 California black oak
    0.23 0.22;  # 32 California white/valley oak
    0.23 0.22;  # 33 interior live oak
    0.25 0.25;  # 34 tanoak
    0.25 0.25;  # 35 giant chinkapin
    0.25 0.25;  # 36 quaking aspen
    0.20 0.20;  # 37 California-laurel
    0.20 0.20;  # 38 Pacific madrone
    0.20 0.20;  # 39 Pacific dogwood
    0.20 0.20;  # 40 bigleaf maple
    0.25 0.25;  # 41 curlleaf mtn-mahogany
    0.20 0.10;  # 42 other softwoods (use lodgepole pine NI)
    0.23 0.22;  # 43 other hardwoods (use California black oak)
]

const _WS_FULIVI = Float32[    # initiating stands (10% cover)  (herb, shrub)  — ws/fmcba.f:192-236
    0.40 1.00;  #  1 sugar pine (use lodgepole pine NI)
    0.40 2.00;  #  2 Douglas-fir
    0.30 2.00;  #  3 white fir
    0.40 2.00;  #  4 giant sequoia (use Douglas-fir)
    0.40 2.00;  #  5 incense cedar (use Douglas-fir)
    0.40 1.00;  #  6 Jeffrey pine (use sugar pine)
    0.30 2.00;  #  7 red fir (use white fir)
    0.25 1.00;  #  8 ponderosa pine
    0.40 1.00;  #  9 lodgepole pine
    0.40 1.00;  # 10 whitebark pine
    0.30 2.00;  # 11 western white pine
    0.13 1.63;  # 12 singleleaf pinyon
    0.30 2.00;  # 13 Pacific silver fir (use white fir)
    0.40 1.00;  # 14 knobcone pine
    0.40 1.00;  # 15 foxtail pine
    0.40 1.00;  # 16 Coulter pine
    0.40 1.00;  # 17 limber pine
    0.40 2.00;  # 18 Monterey pine
    0.55 0.35;  # 19 gray/CA foothill pine
    0.40 1.00;  # 20 washoe pine
    0.13 1.63;  # 21 GB bristlecone pine
    0.40 2.00;  # 22 bigcone Douglas-fir (use Douglas-fir)
    0.40 2.00;  # 23 redwood (use giant sequoia)
    0.30 2.00;  # 24 mountain hemlock
    0.10 2.06;  # 25 western juniper
    0.10 2.06;  # 26 Utah juniper
    0.10 2.06;  # 27 California juniper
    0.55 0.35;  # 28 California live oak
    0.18 2.00;  # 29 canyon live oak
    0.55 0.35;  # 30 blue oak
    0.55 0.35;  # 31 California black oak
    0.55 0.35;  # 32 California white/valley oak
    0.55 0.35;  # 33 interior live oak
    0.18 2.00;  # 34 tanoak
    0.18 2.00;  # 35 giant chinkapin
    0.18 1.32;  # 36 quaking aspen
    0.40 2.00;  # 37 California-laurel
    0.40 2.00;  # 38 Pacific madrone
    0.40 2.00;  # 39 Pacific dogwood
    0.40 2.00;  # 40 bigleaf maple
    0.18 1.32;  # 41 curlleaf mtn-mahogany
    0.40 1.00;  # 42 other softwoods (use lodgepole pine NI)
    0.55 0.35;  # 43 other hardwoods (use California black oak)
]

const _WS_FUINIE = Float32[    # established (60% cover), 11 size classes  — ws/fmcba.f:242-285
    0.9 0.9 1.2 7.0 8.0  0.0 0.0 0.0 0.0 0.6 15.0;  #  1 sugar pine
    0.9 0.9 1.6 3.5 3.5  0.0 0.0 0.0 0.0 0.6 10.0;  #  2 Douglas-fir
    0.7 0.7 3.0 7.0 7.0  0.0 0.0 0.0 0.0 0.6 25.0;  #  3 white fir
    0.9 0.9 1.6 3.5 3.5  0.0 0.0 0.0 0.0 0.6 10.0;  #  4 giant sequoia (use DF)
    0.9 0.9 1.6 3.5 3.5  0.0 0.0 0.0 0.0 0.6 10.0;  #  5 incense cedar (use DF)
    0.9 0.9 1.2 7.0 8.0  0.0 0.0 0.0 0.0 0.6 15.0;  #  6 Jeffrey pine (use sugar pine)
    0.7 0.7 3.0 7.0 7.0  0.0 0.0 0.0 0.0 0.6 25.0;  #  7 red fir (use white fir)
    0.9 0.9 1.2 7.0 8.0  0.0 0.0 0.0 0.0 0.6 15.0;  #  8 ponderosa pine (use Jeffrey pine)
    0.7 0.7 1.6 2.5 2.5  0.0 0.0 0.0 0.0 1.4  5.0;  #  9 lodgepole pine
    0.7 0.7 1.6 2.5 2.5  0.0 0.0 0.0 0.0 1.4  5.0;  # 10 whitebark pine (use LP)
    1.0 1.0 1.6 10.0 10.0 10.0 0.0 0.0 0.0 0.8 30.0;  # 11 western white pine
    0.2 0.8 2.3 1.4 3.0  0.0 0.0 0.0 0.0 0.5  0.0;  # 12 singleleaf pinyon
    0.7 0.7 3.0 7.0 7.0  0.0 0.0 0.0 0.0 0.6 25.0;  # 13 Pacific silver fir (use WF)
    0.7 0.7 1.6 2.5 2.5  0.0 0.0 0.0 0.0 1.4  5.0;  # 14 knobcone pine
    0.7 0.7 1.6 2.5 2.5  0.0 0.0 0.0 0.0 1.4  5.0;  # 15 foxtail pine
    0.7 0.7 1.6 2.5 2.5  0.0 0.0 0.0 0.0 1.4  5.0;  # 16 Coulter pine
    0.7 0.7 1.6 2.5 2.5  0.0 0.0 0.0 0.0 1.4  5.0;  # 17 limber pine
    0.9 0.9 1.6 3.5 3.5  0.0 0.0 0.0 0.0 0.6 10.0;  # 18 Monterey pine
    0.3 0.7 1.4 0.2 0.1  0.0 0.0 0.0 0.0 3.9  0.0;  # 19 gray/CA foothill pine
    0.7 0.7 1.6 2.5 2.5  0.0 0.0 0.0 0.0 1.4  5.0;  # 20 washoe pine
    0.2 0.8 2.3 1.4 3.0  0.0 0.0 0.0 0.0 0.5  0.0;  # 21 GB bristlecone pine
    0.9 0.9 1.6 3.5 3.5  0.0 0.0 0.0 0.0 0.6 10.0;  # 22 bigcone Douglas-fir (use DF)
    0.9 0.9 1.6 3.5 3.5  0.0 0.0 0.0 0.0 0.6 10.0;  # 23 redwood (use giant sequoia)
    1.1 1.1 2.2 10.0 10.0 0.0 0.0 0.0 0.0 0.6 30.0;  # 24 mountain hemlock
    0.1 0.2 0.4 0.5 0.8  1.0 0.0 0.0 0.0 0.1  0.0;  # 25 western juniper
    0.1 0.2 0.4 0.5 0.8  1.0 0.0 0.0 0.0 0.1  0.0;  # 26 Utah juniper
    0.1 0.2 0.4 0.5 0.8  1.0 0.0 0.0 0.0 0.1  0.0;  # 27 California juniper
    0.1 0.1 0.0 0.0 0.0  0.0 0.0 0.0 0.0 2.9  0.0;  # 28 California live oak
    0.2 0.6 2.4 3.6 5.6  0.0 0.0 0.0 0.0 1.4 16.8;  # 29 canyon live oak
    0.3 0.7 1.4 0.2 0.1  0.0 0.0 0.0 0.0 3.9  0.0;  # 30 blue oak
    0.3 0.7 1.4 0.2 0.1  0.0 0.0 0.0 0.0 3.9  0.0;  # 31 California black oak
    0.3 0.7 1.4 0.2 0.1  0.0 0.0 0.0 0.0 3.9  0.0;  # 32 CA white/valley oak
    0.3 0.7 1.4 0.2 0.1  0.0 0.0 0.0 0.0 3.9  0.0;  # 33 interior live oak
    0.2 0.6 2.4 3.6 5.6  0.0 0.0 0.0 0.0 1.4 16.8;  # 34 tanoak (aspen — Ottmar)
    0.2 0.6 2.4 3.6 5.6  0.0 0.0 0.0 0.0 1.4 16.8;  # 35 giant chinkapin
    0.2 0.6 2.4 3.6 5.6  0.0 0.0 0.0 0.0 1.4 16.8;  # 36 quaking aspen
    0.9 0.9 1.6 3.5 3.5  0.0 0.0 0.0 0.0 0.6 10.0;  # 37 California-laurel
    0.9 0.9 1.6 3.5 3.5  0.0 0.0 0.0 0.0 0.6 10.0;  # 38 Pacific madrone
    0.9 0.9 1.6 3.5 3.5  0.0 0.0 0.0 0.0 0.6 10.0;  # 39 Pacific dogwood
    0.9 0.9 1.6 3.5 3.5  0.0 0.0 0.0 0.0 0.6 10.0;  # 40 bigleaf maple
    0.2 0.6 2.4 3.6 5.6  0.0 0.0 0.0 0.0 1.4 16.8;  # 41 curlleaf mtn-mahogany (CR AS)
    0.9 0.9 1.2 7.0 8.0  0.0 0.0 0.0 0.0 0.6 15.0;  # 42 other softwoods (LP NI)
    0.3 0.7 1.4 0.2 0.1  0.0 0.0 0.0 0.0 3.9  0.0;  # 43 other hardwoods (Gambel oak — Ottmar)
]

const _WS_FUINII = Float32[    # initiating (10% cover), 11 size classes  — ws/fmcba.f:291-334
    0.6 0.7 0.8 2.8 3.2  0.0 0.0 0.0 0.0 0.3  7.0;  #  1 sugar pine
    0.5 0.5 1.0 1.4 1.4  0.0 0.0 0.0 0.0 0.3  5.0;  #  2 Douglas-fir
    0.5 0.5 2.0 2.8 2.8  0.0 0.0 0.0 0.0 0.3 12.0;  #  3 white fir
    0.5 0.5 1.0 1.4 1.4  0.0 0.0 0.0 0.0 0.3  5.0;  #  4 giant sequoia (use DF)
    0.5 0.5 1.0 1.4 1.4  0.0 0.0 0.0 0.0 0.3  5.0;  #  5 incense cedar (use DF)
    0.6 0.7 0.8 2.8 3.2  0.0 0.0 0.0 0.0 0.3  7.0;  #  6 Jeffrey pine (use sugar pine)
    0.5 0.5 2.0 2.8 2.8  0.0 0.0 0.0 0.0 0.3 12.0;  #  7 red fir (use white fir)
    0.6 0.7 0.8 2.8 3.2  0.0 0.0 0.0 0.0 0.3  7.0;  #  8 ponderosa pine (use Jeffrey pine)
    0.1 0.1 0.2 0.5 0.5  0.0 0.0 0.0 0.0 0.5  0.8;  #  9 lodgepole pine
    0.1 0.1 0.2 0.5 0.5  0.0 0.0 0.0 0.0 0.5  0.8;  # 10 whitebark pine (use LP)
    0.6 0.6 0.8 6.0 6.0  6.0 0.0 0.0 0.0 0.4 12.0;  # 11 western white pine
    0.0 0.1 0.0 0.0 0.0  0.0 0.0 0.0 0.0 0.3  0.0;  # 12 singleleaf pinyon
    0.5 0.5 2.0 2.8 2.8  0.0 0.0 0.0 0.0 0.3 12.0;  # 13 Pacific silver fir (use WF)
    0.1 0.1 0.2 0.5 0.5  0.0 0.0 0.0 0.0 0.5  0.8;  # 14 knobcone pine
    0.1 0.1 0.2 0.5 0.5  0.0 0.0 0.0 0.0 0.5  0.8;  # 15 foxtail pine
    0.1 0.1 0.2 0.5 0.5  0.0 0.0 0.0 0.0 0.5  0.8;  # 16 Coulter pine
    0.1 0.1 0.2 0.5 0.5  0.0 0.0 0.0 0.0 0.5  0.8;  # 17 limber pine
    0.5 0.5 1.0 1.4 1.4  0.0 0.0 0.0 0.0 0.3  5.0;  # 18 Monterey pine
    0.1 0.1 0.0 0.0 0.0  0.0 0.0 0.0 0.0 2.9  0.0;  # 19 gray/CA foothill pine
    0.1 0.1 0.2 0.5 0.5  0.0 0.0 0.0 0.0 0.5  0.8;  # 20 washoe pine
    0.0 0.1 0.0 0.0 0.0  0.0 0.0 0.0 0.0 0.3  0.0;  # 21 GB bristlecone pine
    0.5 0.5 1.0 1.4 1.4  0.0 0.0 0.0 0.0 0.3  5.0;  # 22 bigcone Douglas-fir (use DF)
    0.5 0.5 1.0 1.4 1.4  0.0 0.0 0.0 0.0 0.3  5.0;  # 23 redwood (use giant sequoia)
    0.7 0.7 1.6 4.0 4.0  0.0 0.0 0.0 0.0 0.3 12.0;  # 24 mountain hemlock
    0.2 0.4 0.2 0.0 0.0  0.0 0.0 0.0 0.0 0.2  0.0;  # 25 western juniper
    0.2 0.4 0.2 0.0 0.0  0.0 0.0 0.0 0.0 0.2  0.0;  # 26 Utah juniper
    0.2 0.4 0.2 0.0 0.0  0.0 0.0 0.0 0.0 0.2  0.0;  # 27 California juniper
    0.1 0.1 0.0 0.0 0.0  0.0 0.0 0.0 0.0 2.9  0.0;  # 28 California live oak
    0.1 0.4 5.0 2.2 2.3  0.0 0.0 0.0 0.0 0.8  5.6;  # 29 canyon live oak
    0.1 0.1 0.0 0.0 0.0  0.0 0.0 0.0 0.0 2.9  0.0;  # 30 blue oak
    0.1 0.1 0.0 0.0 0.0  0.0 0.0 0.0 0.0 2.9  0.0;  # 31 California black oak
    0.1 0.1 0.0 0.0 0.0  0.0 0.0 0.0 0.0 2.9  0.0;  # 32 CA white/valley oak
    0.1 0.1 0.0 0.0 0.0  0.0 0.0 0.0 0.0 2.9  0.0;  # 33 interior live oak
    0.1 0.4 5.0 2.2 2.3  0.0 0.0 0.0 0.0 0.8  5.6;  # 34 tanoak (aspen — Ottmar)
    0.1 0.4 5.0 2.2 2.3  0.0 0.0 0.0 0.0 0.8  5.6;  # 35 giant chinkapin
    0.1 0.4 5.0 2.2 2.3  0.0 0.0 0.0 0.0 0.8  5.6;  # 36 quaking aspen
    0.5 0.5 1.0 1.4 1.4  0.0 0.0 0.0 0.0 0.3  5.0;  # 37 California-laurel
    0.5 0.5 1.0 1.4 1.4  0.0 0.0 0.0 0.0 0.3  5.0;  # 38 Pacific madrone
    0.5 0.5 1.0 1.4 1.4  0.0 0.0 0.0 0.0 0.3  5.0;  # 39 Pacific dogwood
    0.5 0.5 1.0 1.4 1.4  0.0 0.0 0.0 0.0 0.3  5.0;  # 40 bigleaf maple
    0.1 0.4 5.0 2.2 2.3  0.0 0.0 0.0 0.0 0.8  5.6;  # 41 curlleaf mtn-mahogany (CR AS)
    0.6 0.7 0.8 2.8 3.2  0.0 0.0 0.0 0.0 0.3  7.0;  # 42 other softwoods (LP NI)
    0.1 0.1 0.0 0.0 0.0  0.0 0.0 0.0 0.0 2.9  0.0;  # 43 other hardwoods (Gambel oak — Ottmar)
]

# WS live herb/shrub fuel (ws/fmcba.f:519-533): FLIVE(I) = Σ_{j=1,2} ALGSLP(PERCOV, [10,60],
# [FULIVI(I,COVCA(j))·COVCAWT(j), FULIVE(I,COVCA(j))·COVCAWT(j)]) — the TOP-2 cover-type interpolation.
@inline function ws_live_fuel_loading(covca::NTuple{2,Int}, covcawt::NTuple{2,Float32}, percov::Float32)::NTuple{2,Float32}
    herb = 0f0; shrub = 0f0
    @inbounds for j in 1:2
        c = covca[j]
        (1 <= c <= 43) || continue
        wt = covcawt[j]
        herb  += _cr_algslp2(percov, 10f0, 60f0, _WS_FULIVI[c, 1] * wt, _WS_FULIVE[c, 1] * wt)
        shrub += _cr_algslp2(percov, 10f0, 60f0, _WS_FULIVI[c, 2] * wt, _WS_FULIVE[c, 2] * wt)
    end
    return (herb, shrub)
end

# WS initial dead surface fuel (ws/fmcba.f:587-597): STFUEL(ISZ,2) = Σ_{j=1,2} ALGSLP(PERCOV, [10,60],
# [FUINII(ISZ,COVCA(j))·COVCAWT(j), FUINIE(ISZ,COVCA(j))·COVCAWT(j)]) — the TOP-2 cover-type interpolation.
function ws_dead_fuel_loading(covca::NTuple{2,Int}, covcawt::NTuple{2,Float32}, percov::Float32)::Vector{Float32}
    out = zeros(Float32, 11)
    @inbounds for isz in 1:11
        for j in 1:2
            c = covca[j]
            (1 <= c <= 43) || continue
            wt = covcawt[j]
            out[isz] += _cr_algslp2(percov, 10f0, 60f0, _WS_FUINII[c, isz] * wt, _WS_FUINIE[c, isz] * wt)
        end
    end
    return out
end

# =============================================================================
# ws_r5crwd — WS crown width (ft), the R5 California crown-width routine (ws/r5crwd.f). WS ALWAYS uses this
# (ws/cwcalc.f:379-381 `IF(VARACD.EQ.'WS') CALL R5CRWD; GO TO 9000`) — the WSMAP/Crookston-R6 equations in
# cwcalc.f are declared but never reached for WS. Needed by FMCBA's PERCOV (percent canopy cover); the generic
# `crown_width` returns 0.5 for every WS species, collapsing PERCOV≈0 and mis-selecting the initiating-stand
# fuel loads (measured: PERCOV 0.27 → wrong FUINII path; R5CRWD → PERCOV 41.99 = oracle). MAPWS maps WS species
# 1..43 → the coefficient index INDX (1..35). Warbington/Levitan for D≥SPLINE (eqn type 1 linear / 2 power /
# 3 quadratic), Dixon (DX1+DX2·D) for 4.5'≤H below the spline DBH, and SM·H below 4.5 ft. Verified bit-exact
# vs FVSws_g16 (RF 6.71+0.421·14.47=12.80=oracle; WF 5.82+0.591·20.27=17.80=oracle). Ported verbatim.
const _WS_R5_MAPWS = Int[
     2,  1,  4, 27,  7, 22,  5, 21, 23, 24,
     3, 27,  4, 27, 27, 32, 27, 27, 33, 27,
    27, 29, 27,  6, 25, 25, 25, 34,  9, 26,
     8, 35, 28, 10, 20, 18, 30, 11, 28, 28,
    28, 27, 28]
const _WS_R5_WB1 = Float32[
    6.81, -1.476, -0.997, 5.82, 6.71, 4.72, 7.11, 10.0,
    5.0, 10.0, 1.0, 6.19, 6.50, 4.57, 4.2, 4.00,
    8.00, 0.50, 3.08, 2.98, 2.24, 1.52, 1.91, 2.37,
    4.31, 4.49, 6.0, 2.0, 27.030, 12.733, 9.0684, 3.9347,
    3.8273, 5.3732, 4.5628]
const _WS_R5_WB2 = Float32[
    0.732, 1.01, 0.92, 0.591, 0.421, 0.608, 0.470, 1.20,
    1.69, 1.05, 1.43, 1.01, 1.80, 1.41, 1.42, 1.60,
    1.53, 1.62, 1.92, 1.55, 0.763, 0.891, 0.784, 0.736,
    0.628, 0.688, 0.6, 1.5, 0.8612, 2.249, 0.4702, 0.7086,
    0.7624, 0.7707, 0.6925]
const _WS_R5_WB3 = Float32[fill(0f0, 19); -0.014f0; fill(0f0, 15)]   # only INDX 20 nonzero
const _WS_R5_DX1 = Float32[
    3.62, 3.5, 3.5, 3.26, 3.5, 3.5, 3.5, 2.5,
    2.5, 2.23, 3.11, 3.5, 3.5, 3.5, 3.5, 3.5,
    2.5, 2.5, 2.5, 2.15, 3.77, 3.5, 3.5, 3.5,
    3.5, 2.5, 3.5, 2.5, 3.5, 2.5, 3.5, 3.5,
    3.5, 2.5, 2.5]
const _WS_R5_DX2 = Float32[
    1.370, 0.338, 0.329, 1.103, 1.063, 0.852, 1.192, 2.700,
    2.190, 1.630, 1.008, 1.548, 2.400, 1.624, 1.560, 1.700,
    2.630, 1.220, 2.036, 1.646, 0.7756, 0.5754, 0.6492, 0.8496,
    1.6684, 2.2175, 1.1, 1.4, 5.5672, 4.2956, 3.1654, 1.7618,
    1.9108, 3.2150, 2.2816]
# IEQN: 19*1, 3, 6*2, 2*1, 2*1, 5*2  (r5crwd.f)
const _WS_R5_IEQN = Int[fill(1, 19); 3; fill(2, 6); 1; 1; 1; 1; fill(2, 5)]
const _WS_R5_SPLINE = Float32[
    5.0, 7.4, 7.6, 5.0, 5.0, 5.0, 5.0, 5.0,
    5.0, 13.4, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
    5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
    5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
    5.0, 5.0, 5.0]
const _WS_R5_SM = Float32[
    0.7778, 0.7778, 0.7778, 0.7778, 0.7778, 0.7778, 0.7778, 0.5556,
    0.5556, 0.5556, 0.5556, 0.7778, 0.7778, 0.7778, 0.7778, 0.7778,
    0.5556, 0.5556, 0.5556, 0.5556, 0.7778, 0.7778, 0.7778, 0.7778,
    0.7778, 0.5556, 0.7778, 0.5556, 0.7778, 0.5556, 0.7778, 0.7778,
    0.7778, 0.5556, 0.5556]

"""
    ws_r5crwd(sp, d, h) -> Float32

WS crown width (ft) for WS species `sp` (1..43), DBH `d` (in), total height `h` (ft) — ws/r5crwd.f.
D≥SPLINE(INDX): Warbington/Levitan (linear/power/quadratic per IEQN); 4.5'≤H<spline-DBH: Dixon DX1+DX2·D;
H<4.5': SM·H straight line. INDX 99 (unmapped variant) ⇒ 0.
"""
@inline function ws_r5crwd(sp::Int, d::Float32, h::Float32)::Float32
    (1 <= sp <= 43) || return 0f0
    indx = _WS_R5_MAPWS[sp]
    itype = _WS_R5_IEQN[indx]; spdbh = _WS_R5_SPLINE[indx]
    if d >= spdbh
        if     itype == 1; return _WS_R5_WB1[indx] + _WS_R5_WB2[indx] * d
        elseif itype == 2; return _WS_R5_WB1[indx] * fpow(d, _WS_R5_WB2[indx])
        elseif itype == 3; return _WS_R5_WB1[indx] + _WS_R5_WB2[indx] * d + _WS_R5_WB3[indx] * d * d
        else               return 0f0
        end
    elseif h >= 4.5f0
        return _WS_R5_DX1[indx] + _WS_R5_DX2[indx] * d
    else
        return _WS_R5_SM[indx] * h
    end
end
