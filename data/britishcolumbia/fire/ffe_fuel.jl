# BC FFE surface-fuel loading (canada/fire/bc/fmcba.f DATA FULIVE/FULIVI/FUINIE/FUINII/COVINI). BC's FMCBA is the IE
# FMCBA (ie/fmcba.f) with 15-species tables: live herb/shrub (tons/ac) + the 11-class initial dead load keyed by the
# cover type (the species with most BA), ESTABLISHED (60% cover) and INITIATING (10% cover) rows, ALGSLP-interpolated
# on PERCOV (XCOV=[10,60]). COVTYP==0 (no trees) ⇒ COVINI(ITYPE) — the IE table verbatim (fmcba.f:179-181).
# Index = BC species 1..15. Ported verbatim.

const _BC_FULIVE = Float32[
    0.15 0.1;   # 1 WP
    0.2 0.2;    # 2 WL
    0.2 0.2;    # 3 FD
    0.15 0.1;   # 4 BG
    0.2 0.2;    # 5 HW
    0.2 0.2;    # 6 CW
    0.2 0.1;    # 7 PL
    0.15 0.2;   # 8 SE
    0.15 0.2;   # 9 BL
    0.2 0.25;   # 10 PY
    0.25 0.25;  # 11 EP (birch)
    0.25 0.25;  # 12 AT (aspen)
    0.25 0.25;  # 13 AC (cottonwood)
    0.2 0.2;    # 14 OC (as FD)
    0.25 0.25;  # 15 OH (as birch)
]

const _BC_FULIVI = Float32[
    0.3 2;      # 1
    0.4 2;      # 2
    0.4 2;      # 3
    0.3 2;      # 4
    0.4 2;      # 5
    0.4 2;      # 6
    0.4 1;      # 7
    0.3 2;      # 8
    0.3 2;      # 9
    0.25 0.1;   # 10
    0.18 1.32;  # 11
    0.18 1.32;  # 12
    0.18 1.32;  # 13
    0.4 2;      # 14
    0.18 1.32;  # 15
]

const _BC_FUINIE = Float32[
    1 1 1.6 10 10 10 0 0 0 0.8 30;          # 1
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;     # 2
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;     # 3
    0.7 0.7 3 7 7 0 0 0 0 0.6 25;           # 4
    2.2 2.2 5.2 15 20 15 0 0 0 1 35;        # 5
    2.2 2.2 5.2 15 20 15 0 0 0 1 35;        # 6
    0.9 0.9 1.2 7 8 0 0 0 0 0.6 15;         # 7
    1.1 1.1 2.2 10 10 0 0 0 0 0.6 30;       # 8
    1.1 1.1 2.2 10 10 0 0 0 0 0.6 30;       # 9
    0.7 0.7 1.6 2.5 2.5 0 0 0 0 1.4 5;      # 10
    0.2 0.6 2.4 3.6 5.6 0 0 0 0 1.4 16.8;   # 11
    0.2 0.6 2.4 3.6 5.6 0 0 0 0 1.4 16.8;   # 12
    0.2 0.6 2.4 3.6 5.6 0 0 0 0 1.4 16.8;   # 13
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;     # 14
    0.2 0.6 2.4 3.6 5.6 0 0 0 0 1.4 16.8;   # 15
]

const _BC_FUINII = Float32[
    0.6 0.6 0.8 6 6 6 0 0 0 0.4 12;         # 1
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;        # 2
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;        # 3
    0.5 0.5 2 2.8 2.8 0 0 0 0 0.3 12;       # 4
    1.6 1.6 3.6 6 8 6 0 0 0 0.5 12;         # 5
    1.6 1.6 3.6 6 8 6 0 0 0 0.5 12;         # 6
    0.6 0.7 0.8 2.8 3.2 0 0 0 0 0.3 7;      # 7
    0.7 0.7 1.6 4 4 0 0 0 0 0.3 12;         # 8
    0.7 0.7 1.6 4 4 0 0 0 0 0.3 12;         # 9
    0.1 0.1 0.2 0.5 0.5 0 0 0 0 0.5 0.8;    # 10
    0.1 0.4 5 2.2 2.3 0 0 0 0 0.8 5.6;      # 11
    0.1 0.4 5 2.2 2.3 0 0 0 0 0.8 5.6;      # 12
    0.1 0.4 5 2.2 2.3 0 0 0 0 0.8 5.6;      # 13
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;        # 14
    0.1 0.4 5 2.2 2.3 0 0 0 0 0.8 5.6;      # 15
]

# BC live herb/shrub (bc/fmcba.f:265-271): ALGSLP(PERCOV, [10,60], FULIVI↔FULIVE).
@inline function bc_live_fuel_loading(covtyp::Int, percov::Float32)::NTuple{2,Float32}
    herb  = _cr_algslp2(percov, 10f0, 60f0, _BC_FULIVI[covtyp, 1], _BC_FULIVE[covtyp, 1])
    shrub = _cr_algslp2(percov, 10f0, 60f0, _BC_FULIVI[covtyp, 2], _BC_FULIVE[covtyp, 2])
    return (herb, shrub)
end

# BC initial dead fuel (bc/fmcba.f:289-294): FUINII↔FUINIE by PERCOV per size class → 11-vector (hard).
function bc_dead_fuel_loading(covtyp::Int, percov::Float32)::Vector{Float32}
    out = Vector{Float32}(undef, 11)
    @inbounds for isz in 1:11
        out[isz] = _cr_algslp2(percov, 10f0, 60f0, _BC_FUINII[covtyp, isz], _BC_FUINIE[covtyp, isz])
    end
    return out
end
