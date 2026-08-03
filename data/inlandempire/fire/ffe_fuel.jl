# IE FFE surface-fuel loading (ie/fmcba.f DATA FULIVE/FULIVI/FUINIE/FUINII/COVINI). WESTERN structure:
# live herb/shrub (tons/ac) + dead 11-size-class loading keyed by dominant-species COVER TYPE, with
# separate ESTABLISHED (E, 60% cover) and INITIATING (I, 10% cover) tables, PERCOV-interpolated (XCOV=[10,60]).
# COVTYP==0 (bare) ⇒ COVINI(ITYPE) seral cover species. Index = IE species 1..23. Ported verbatim.

const _IE_FULIVE = Float32[
    0.15 0.1;  # 1
    0.2 0.2;  # 2
    0.2 0.2;  # 3
    0.15 0.1;  # 4
    0.2 0.2;  # 5
    0.2 0.2;  # 6
    0.2 0.1;  # 7
    0.15 0.2;  # 8
    0.15 0.2;  # 9
    0.2 0.25;  # 10
    0.15 0.2;  # 11
    0.15 0.2;  # 12
    0.2 0.25;  # 13
    0.15 0.2;  # 14
    0.2 0.25;  # 15
    0.2 0.25;  # 16
    0.2 0.2;  # 17
    0.25 0.25;  # 18
    0.25 0.25;  # 19
    0.2 0.2;  # 20
    0.25 0.25;  # 21
    0.25 0.25;  # 22
    0.15 0.2;  # 23
]

const _IE_FULIVI = Float32[
    0.3 2;  # 1
    0.4 2;  # 2
    0.4 2;  # 3
    0.3 2;  # 4
    0.4 2;  # 5
    0.4 2;  # 6
    0.4 1;  # 7
    0.3 2;  # 8
    0.3 2;  # 9
    0.25 0.1;  # 10
    0.3 2;  # 11
    0.3 2;  # 12
    0.25 0.1;  # 13
    0.3 2;  # 14
    0.25 0.1;  # 15
    0.25 0.1;  # 16
    0.4 2;  # 17
    0.18 1.32;  # 18
    0.18 1.32;  # 19
    0.4 2;  # 20
    0.18 1.32;  # 21
    0.18 1.32;  # 22
    0.3 2;  # 23
]

const _IE_FUINIE = Float32[
    1 1 1.6 10 10 10 0 0 0 0.8 30;  # 1
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;  # 2
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;  # 3
    0.7 0.7 3 7 7 0 0 0 0 0.6 25;  # 4
    2.2 2.2 5.2 15 20 15 0 0 0 1 35;  # 5
    2.2 2.2 5.2 15 20 15 0 0 0 1 35;  # 6
    0.9 0.9 1.2 7 8 0 0 0 0 0.6 15;  # 7
    1.1 1.1 2.2 10 10 0 0 0 0 0.6 30;  # 8
    1.1 1.1 2.2 10 10 0 0 0 0 0.6 30;  # 9
    0.7 0.7 1.6 2.5 2.5 0 0 0 0 1.4 5;  # 10
    1.1 1.1 2.2 10 10 0 0 0 0 0.6 30;  # 11
    1.1 1.1 2.2 10 10 0 0 0 0 0.6 30;  # 12
    0.7 0.7 1.6 2.5 2.5 0 0 0 0 1.4 5;  # 13
    1.1 1.1 2.2 10 10 0 0 0 0 0.6 30;  # 14
    0.7 0.7 1.6 2.5 2.5 0 0 0 0 1.4 5;  # 15
    0.7 0.7 1.6 2.5 2.5 0 0 0 0 1.4 5;  # 16
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;  # 17
    0.2 0.6 2.4 3.6 5.6 0 0 0 0 1.4 16.8;  # 18
    0.2 0.6 2.4 3.6 5.6 0 0 0 0 1.4 16.8;  # 19
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;  # 20
    0.2 0.6 2.4 3.6 5.6 0 0 0 0 1.4 16.8;  # 21
    0.2 0.6 2.4 3.6 5.6 0 0 0 0 1.4 16.8;  # 22
    1.1 1.1 2.2 10 10 0 0 0 0 0.6 30;  # 23
]

const _IE_FUINII = Float32[
    0.6 0.6 0.8 6 6 6 0 0 0 0.4 12;  # 1
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;  # 2
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;  # 3
    0.5 0.5 2 2.8 2.8 0 0 0 0 0.3 12;  # 4
    1.6 1.6 3.6 6 8 6 0 0 0 0.5 12;  # 5
    1.6 1.6 3.6 6 8 6 0 0 0 0.5 12;  # 6
    0.6 0.7 0.8 2.8 3.2 0 0 0 0 0.3 7;  # 7
    0.7 0.7 1.6 4 4 0 0 0 0 0.3 12;  # 8
    0.7 0.7 1.6 4 4 0 0 0 0 0.3 12;  # 9
    0.1 0.1 0.2 0.5 0.5 0 0 0 0 0.5 0.8;  # 10
    0.7 0.7 1.6 4 4 0 0 0 0 0.3 12;  # 11
    0.7 0.7 1.6 4 4 0 0 0 0 0.3 12;  # 12
    0.1 0.1 0.2 0.5 0.5 0 0 0 0 0.5 0.8;  # 13
    0.7 0.7 1.6 4 4 0 0 0 0 0.3 12;  # 14
    0.1 0.1 0.2 0.5 0.5 0 0 0 0 0.5 0.8;  # 15
    0.1 0.1 0.2 0.5 0.5 0 0 0 0 0.5 0.8;  # 16
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;  # 17
    0.1 0.4 5 2.2 2.3 0 0 0 0 0.8 5.6;  # 18
    0.1 0.4 5 2.2 2.3 0 0 0 0 0.8 5.6;  # 19
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;  # 20
    0.1 0.4 5 2.2 2.3 0 0 0 0 0.8 5.6;  # 21
    0.1 0.4 5 2.2 2.3 0 0 0 0 0.8 5.6;  # 22
    0.7 0.7 1.6 4 4 0 0 0 0 0.3 12;  # 23
]

const _IE_COVINI = Int[10, 10, 3, 3, 5, 4, 3, 3, 3, 8, 8, 4, 4, 6, 6, 6, 5, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9]  # habitat ITYPE(1..30) → seral cover species

# IE live herb/shrub (fmcba.f:283-289): ALGSLP(PERCOV, [10,60], FULIVI↔FULIVE) — reuses the generic 2-pt interp.
@inline function ie_live_fuel_loading(covtyp::Int, percov::Float32)::NTuple{2,Float32}
    (covtyp < 1 || covtyp > 23) && (covtyp = 3)      # Douglas-fir default
    herb  = _cr_algslp2(percov, 10f0, 60f0, _IE_FULIVI[covtyp, 1], _IE_FULIVE[covtyp, 1])
    shrub = _cr_algslp2(percov, 10f0, 60f0, _IE_FULIVI[covtyp, 2], _IE_FULIVE[covtyp, 2])
    return (herb, shrub)
end

# IE dead-fuel loading (fmcba.f:307-312): interpolate FUINII↔FUINIE by PERCOV per size class → 11-vector (hard).
function ie_dead_fuel_loading(covtyp::Int, percov::Float32)::Vector{Float32}
    (covtyp < 1 || covtyp > 23) && (covtyp = 3)
    out = Vector{Float32}(undef, 11)
    @inbounds for isz in 1:11
        out[isz] = _cr_algslp2(percov, 10f0, 60f0, _IE_FUINII[covtyp, isz], _IE_FUINIE[covtyp, isz])
    end
    return out
end

# IE bare-stand seral cover species (fmcba.f:279 COVTYP=COVINI(ITYPE)).
@inline ie_covini_default(itype::Int)::Int = (1 <= itype <= 30) ? _IE_COVINI[itype] : 3
