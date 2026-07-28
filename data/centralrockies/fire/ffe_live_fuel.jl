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
