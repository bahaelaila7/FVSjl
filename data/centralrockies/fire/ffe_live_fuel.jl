# CR FFE live surface-fuel loading (cr fire model, FVScr_buildDir/fmcba.f DATA FULIVE/FULIVI).
# WESTERN structure (≠ eastern forest-type tables): live herb/shrub fuel (tons/acre) keyed by SPECIES,
# with separate ESTABLISHED (FULIVE) and INITIATING-stand-<10%-cover (FULIVI) tables. Selected via the
# stand's seral COVER TYPE (dominant species; COVINI2/COVINI3 for bare stands). Values from J. Brown /
# Ottmar 2000. Index = CR species 1..38. Columns = (herb, shrub).
#
# Ported verbatim from fmcba.f:76-159. The cover-type dispatch (fmcba.f flow) + FUINIT dead fuel + fuel-model
# assignment are the remaining FFE-chunk pieces; this is the FULIV data layer.

const _CR_FULIVE = Float32[    # established stands  (herb, shrub)
    0.15 0.20;  # 1 subalpine fir
    0.15 0.20;  # 2 corkbark fir
    0.20 0.20;  # 3 Douglas-fir
    0.15 0.10;  # 4 grand fir
    0.15 0.10;  # 5 white fir
    0.15 0.20;  # 6 mountain hemlock
    0.20 0.20;  # 7 western redcedar
    0.20 0.20;  # 8 western larch
    0.04 0.05;  # 9 bristlecone pine
    0.20 0.10;  # 10 limber pine
    0.20 0.10;  # 11 lodgepole pine
    0.04 0.05;  # 12 pinyon pine
    0.20 0.25;  # 13 ponderosa pine
    0.20 0.10;  # 14 whitebark pine
    0.15 0.10;  # 15 southwestern white pine
    0.04 0.05;  # 16 Utah juniper
    0.15 0.20;  # 17 blue spruce
    0.15 0.20;  # 18 Engelmann spruce
    0.15 0.20;  # 19 white spruce
    0.25 0.25;  # 20 quaking aspen
    0.25 0.25;  # 21 narrowleaf cottonwood
    0.25 0.25;  # 22 plains cottonwood
    0.23 0.22;  # 23 Gambel oak
    0.23 0.22;  # 24 Arizona white oak
    0.23 0.22;  # 25 emory oak
    0.23 0.22;  # 26 bur oak
    0.23 0.22;  # 27 silverleaf oak
    0.25 0.25;  # 28 paper birch
    0.04 0.05;  # 29 alligator juniper
    0.04 0.05;  # 30 Rocky Mountain juniper
    0.04 0.05;  # 31 oneseed juniper
    0.04 0.05;  # 32 Eastern redcedar
    0.04 0.05;  # 33 singleleaf pinyon
    0.04 0.05;  # 34 border pinyon
    0.04 0.05;  # 35 Arizona pinyon
    0.20 0.25;  # 36 Chihuahua pine
    0.20 0.10;  # 37 other softwoods
    0.25 0.25;  # 38 other hardwoods
]

const _CR_FULIVI = Float32[    # initiating stands (10% cover)  (herb, shrub)
    0.30 2.00;  # 1 subalpine fir
    0.30 2.00;  # 2 corkbark fir
    0.40 2.00;  # 3 Douglas-fir
    0.30 2.00;  # 4 grand fir
    0.30 2.00;  # 5 white fir
    0.30 2.00;  # 6 mountain hemlock
    0.40 2.00;  # 7 western redcedar
    0.40 2.00;  # 8 western larch
    0.13 1.63;  # 9 bristlecone pine
    0.40 1.00;  # 10 limber pine
    0.40 1.00;  # 11 lodgepole pine
    0.13 1.63;  # 12 pinyon pine
    0.25 0.10;  # 13 ponderosa pine
    0.40 1.00;  # 14 whitebark pine
    0.30 2.00;  # 15 southwestern white pine
    0.13 1.63;  # 16 Utah juniper
    0.30 2.00;  # 17 blue spruce
    0.30 2.00;  # 18 Engelmann spruce
    0.30 2.00;  # 19 white spruce
    0.18 1.32;  # 20 quaking aspen
    0.18 1.32;  # 21 narrowleaf cottonwood
    0.18 1.32;  # 22 plains cottonwood
    0.55 0.35;  # 23 Gambel oak
    0.55 0.35;  # 24 Arizona white oak
    0.55 0.35;  # 25 emory oak
    0.55 0.35;  # 26 bur oak
    0.55 0.35;  # 27 silverleaf oak
    0.18 1.32;  # 28 paper birch
    0.13 1.63;  # 29 alligator juniper
    0.13 1.63;  # 30 Rocky Mountain juniper
    0.13 1.63;  # 31 oneseed juniper
    0.13 1.63;  # 32 Eastern redcedar
    0.13 1.63;  # 33 singleleaf pinyon
    0.13 1.63;  # 34 border pinyon
    0.13 1.63;  # 35 Arizona pinyon
    0.25 0.10;  # 36 Chihuahua pine
    0.40 1.00;  # 37 other softwoods
    0.18 1.32;  # 38 other hardwoods
]

# ALGSLP 2-point linear interpolation with end-clamping (FVS ALGSLP for n=2).
@inline function _cr_algslp2(x::Float32, x1::Float32, x2::Float32, y1::Float32, y2::Float32)::Float32
    x <= x1 && return y1
    x >= x2 && return y2
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
end

# CR live herb/shrub fuel (fmcba.f:443-449): interpolate between INITIATING (10% cover, FULIVI) and
# ESTABLISHED (60% cover, FULIVE) by the stand's PERCOV, keyed by the dominant-species cover type COVTYP.
@inline function cr_live_fuel_loading(covtyp::Int, percov::Float32)::NTuple{2,Float32}
    (covtyp < 1 || covtyp > 38) && (covtyp = 11)     # LP default (fmcba.f:432)
    herb  = _cr_algslp2(percov, 10f0, 60f0, _CR_FULIVI[covtyp, 1], _CR_FULIVE[covtyp, 1])
    shrub = _cr_algslp2(percov, 10f0, 60f0, _CR_FULIVI[covtyp, 2], _CR_FULIVE[covtyp, 2])
    return (herb, shrub)
end

# CR FFE dead surface-fuel loading (fmcba.f:165-249 DATA FUINIE/FUINII). Per species × 11 size classes
# (<.25, .25-1, 1-3, 3-6, 6-12, 12-20, 20-35, 35-50, >50, Litter, Duff). Established (60% cover) vs
# initiating (10%); interpolated by PERCOV (fmcba.f:465-472). From J. Brown / Ottmar 2000.
const _CR_FUINIE = Float32[    # established (60% cover)
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;  # 1 subalpine fir
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;  # 2 corkbark fir
    0.9 0.9 1.6  3.5  3.5  0.0 0.0 0.0 0.0 0.6 10.0;  # 3 Douglas-fir
    0.7 0.7 3.0  7.0  7.0  0.0 0.0 0.0 0.0 0.6 25.0;  # 4 grand fir
    0.7 0.7 3.0  7.0  7.0  0.0 0.0 0.0 0.0 0.6 25.0;  # 5 white fir
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;  # 6 mountain hemlock
    1.6 1.6 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;  # 7 western redcedar
    0.9 0.9 1.6  3.5  3.5  0.0 0.0 0.0 0.0 0.6 10.0;  # 8 western larch
    0.2 0.8 2.3  1.4  3.0  0.0 0.0 0.0 0.0 0.5  0.0;  # 9 bristlecone pine
    0.9 0.9 1.2  7.0  8.0  0.0 0.0 0.0 0.0 0.6 15.0;  # 10 limber pine
    0.9 0.9 1.2  7.0  8.0  0.0 0.0 0.0 0.0 0.6 15.0;  # 11 lodgepole pine
    0.2 0.8 2.3  1.4  3.0  0.0 0.0 0.0 0.0 0.5  0.0;  # 12 pinyon pine
    0.7 0.7 1.6  2.5  2.5  0.0 0.0 0.0 0.0 1.4  5.0;  # 13 ponderosa pine
    0.9 0.9 1.2  7.0  8.0  0.0 0.0 0.0 0.0 0.6 15.0;  # 14 whitebark pine
    1.0 1.0 1.6 10.0 10.0 10.0 0.0 0.0 0.0 0.8 30.0;  # 15 SW white pine
    0.2 0.8 2.3  1.4  3.0  0.0 0.0 0.0 0.0 0.5  0.0;  # 16 Utah juniper
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;  # 17 blue spruce
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;  # 18 Engelmann spruce
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;  # 19 white spruce
    0.2 0.6 2.4  3.6  5.6  0.0 0.0 0.0 0.0 1.4 16.8;  # 20 quaking aspen
    0.2 0.6 2.4  3.6  5.6  0.0 0.0 0.0 0.0 1.4 16.8;  # 21 narrowleaf cottonwood
    0.2 0.6 2.4  3.6  5.6  0.0 0.0 0.0 0.0 1.4 16.8;  # 22 plains cottonwood
    0.3 0.7 1.4  0.2  0.1  0.0 0.0 0.0 0.0 3.9  0.0;  # 23 Gambel oak
    0.3 0.7 1.4  0.2  0.1  0.0 0.0 0.0 0.0 3.9  0.0;  # 24 Arizona white oak
    0.3 0.7 1.4  0.2  0.1  0.0 0.0 0.0 0.0 3.9  0.0;  # 25 emory oak
    0.3 0.7 1.4  0.2  0.1  0.0 0.0 0.0 0.0 3.9  0.0;  # 26 bur oak
    0.3 0.7 1.4  0.2  0.1  0.0 0.0 0.0 0.0 3.9  0.0;  # 27 silverleaf oak
    0.2 0.6 2.4  3.6  5.6  0.0 0.0 0.0 0.0 1.4 16.8;  # 28 paper birch
    0.2 0.8 2.3  1.4  3.0  0.0 0.0 0.0 0.0 0.5  0.0;  # 29 alligator juniper
    0.2 0.8 2.3  1.4  3.0  0.0 0.0 0.0 0.0 0.5  0.0;  # 30 Rocky Mtn juniper
    0.2 0.8 2.3  1.4  3.0  0.0 0.0 0.0 0.0 0.5  0.0;  # 31 oneseed juniper
    0.2 0.8 2.3  1.4  3.0  0.0 0.0 0.0 0.0 0.5  0.0;  # 32 Eastern redcedar
    0.2 0.8 2.3  1.4  3.0  0.0 0.0 0.0 0.0 0.5  0.0;  # 33 singleleaf pinyon
    0.2 0.8 2.3  1.4  3.0  0.0 0.0 0.0 0.0 0.5  0.0;  # 34 border pinyon
    0.2 0.8 2.3  1.4  3.0  0.0 0.0 0.0 0.0 0.5  0.0;  # 35 Arizona pinyon
    0.7 0.7 1.6  2.5  2.5  0.0 0.0 0.0 0.0 1.4  5.0;  # 36 Chihuahua pine
    0.9 0.9 1.2  7.0  8.0  0.0 0.0 0.0 0.0 0.6 15.0;  # 37 other softwoods
    0.2 0.6 2.4  3.6  5.6  0.0 0.0 0.0 0.0 1.4 16.8;  # 38 other hardwoods
]
const _CR_FUINII = Float32[    # initiating (10% cover)
    0.7 0.7 1.6 4.0 4.0 0.0 0.0 0.0 0.0 0.3 12.0;  # 1 subalpine fir
    0.7 0.7 1.6 4.0 4.0 0.0 0.0 0.0 0.0 0.3 12.0;  # 2 corkbark fir
    0.5 0.5 1.0 1.4 1.4 0.0 0.0 0.0 0.0 0.3  5.0;  # 3 Douglas-fir
    0.5 0.5 2.0 2.8 2.8 0.0 0.0 0.0 0.0 0.3 12.0;  # 4 grand fir
    0.5 0.5 2.0 2.8 2.8 0.0 0.0 0.0 0.0 0.3 12.0;  # 5 white fir
    0.7 0.7 1.6 4.0 4.0 0.0 0.0 0.0 0.0 0.3 12.0;  # 6 mountain hemlock
    1.6 1.6 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 7 western redcedar
    0.5 0.5 1.0 1.4 1.4 0.0 0.0 0.0 0.0 0.3  5.0;  # 8 western larch
    0.0 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3  0.0;  # 9 bristlecone pine
    0.6 0.7 0.8 2.8 3.2 0.0 0.0 0.0 0.0 0.3  7.0;  # 10 limber pine
    0.6 0.7 0.8 2.8 3.2 0.0 0.0 0.0 0.0 0.3  7.0;  # 11 lodgepole pine
    0.0 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3  0.0;  # 12 pinyon pine
    0.1 0.1 0.2 0.5 0.5 0.0 0.0 0.0 0.0 0.5  0.8;  # 13 ponderosa pine
    0.6 0.7 0.8 2.8 3.2 0.0 0.0 0.0 0.0 0.3  7.0;  # 14 whitebark pine
    0.6 0.6 0.8 6.0 6.0 6.0 0.0 0.0 0.0 0.4 12.0;  # 15 SW white pine
    0.0 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3  0.0;  # 16 Utah juniper
    0.7 0.7 1.6 4.0 4.0 0.0 0.0 0.0 0.0 0.3 12.0;  # 17 blue spruce
    0.7 0.7 1.6 4.0 4.0 0.0 0.0 0.0 0.0 0.3 12.0;  # 18 Engelmann spruce
    0.7 0.7 1.6 4.0 4.0 0.0 0.0 0.0 0.0 0.3 12.0;  # 19 white spruce
    0.1 0.4 5.0 2.2 2.3 0.0 0.0 0.0 0.0 0.8  5.6;  # 20 quaking aspen
    0.1 0.4 5.0 2.2 2.3 0.0 0.0 0.0 0.0 0.8  5.6;  # 21 narrowleaf cottonwood
    0.1 0.4 5.0 2.2 2.3 0.0 0.0 0.0 0.0 0.8  5.6;  # 22 plains cottonwood
    0.1 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.9  0.0;  # 23 Gambel oak
    0.1 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.9  0.0;  # 24 Arizona white oak
    0.1 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.9  0.0;  # 25 emory oak
    0.1 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.9  0.0;  # 26 bur oak
    0.1 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.9  0.0;  # 27 silverleaf oak
    0.1 0.4 5.0 2.2 2.3 0.0 0.0 0.0 0.0 0.8  5.6;  # 28 paper birch
    0.0 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3  0.0;  # 29 alligator juniper
    0.0 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3  0.0;  # 30 Rocky Mtn juniper
    0.0 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3  0.0;  # 31 oneseed juniper
    0.0 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3  0.0;  # 32 Eastern redcedar
    0.0 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3  0.0;  # 33 singleleaf pinyon
    0.0 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3  0.0;  # 34 border pinyon
    0.0 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3  0.0;  # 35 Arizona pinyon
    0.1 0.1 0.2 0.5 0.5 0.0 0.0 0.0 0.0 0.5  0.8;  # 36 Chihuahua pine
    0.6 0.7 0.8 2.8 3.2 0.0 0.0 0.0 0.0 0.3  7.0;  # 37 other softwoods
    0.1 0.4 5.0 2.2 2.3 0.0 0.0 0.0 0.0 0.8  5.6;  # 38 other hardwoods
]

# CR dead-fuel loading (fmcba.f:465-472): interpolate FUINII↔FUINIE by PERCOV per size class → 11-vector (hard).
function cr_dead_fuel_loading(covtyp::Int, percov::Float32)::Vector{Float32}
    (covtyp < 1 || covtyp > 38) && (covtyp = 11)
    out = Vector{Float32}(undef, 11)
    @inbounds for isz in 1:11
        out[isz] = _cr_algslp2(percov, 10f0, 60f0, _CR_FUINII[covtyp, isz], _CR_FUINIE[covtyp, isz])
    end
    return out
end
