# =============================================================================
# sprout.jl — stump-sprout regeneration sub-routines (NSPREC / SPRTHT / ESSPRT)
#
# Ported from: bin/FVSsn_buildDir/essprt.f (SN-variant SELECT CASE blocks).
# These are the three pure helpers consumed by ESUCKR (esuckr.f) when a harvested
# stump regenerates as stump/root sprouts:
#
#   nsprec_sn  — number of sprouts produced per stump   (species + stump DBH)
#   sprtht_sn  — sprout height at a given sprout age     (species + site index)
#   essprt_sn  — per-record survival multiplier on the carried sprout TPA (PREM)
#
# Each is bit-faithful to the Fortran SN block. The ESSPRT per-species coefficient
# blob lives in data/southern/sprout_essprt.csv (loaded as the per-species columns
# essprt_kind / essprt_p1 / essprt_p2 / essprt_fsp). NSPREC and SPRTHT are tiny
# piecewise rules, so they stay inline rather than in a CSV.
#
# Status: pure functions, unit-tested in isolation. Wired into the ESUCKR
# generation loop in Chunk C (until then they have no .sum effect).
# =============================================================================

"""
    nsprec_sn(ispc, dstmp) -> Int

Number of sprouts produced by one cut stump (`NSPREC`, essprt.f:1102-1120, SN).
`dstmp` is the stump diameter (in). Most species yield a single sprout; species
5 yields one only below 7 in, and the oak/sweetgum group {33,61,80,82} ramps
1→3 over the 5–10 in stump-DBH range.
"""
function nsprec_sn(ispc::Integer, dstmp::Float32)::Int
    if ispc == 5
        return dstmp < 7f0 ? 1 : 0
    elseif ispc == 33 || ispc == 61 || ispc == 80 || ispc == 82
        if dstmp < 5f0
            return 1
        elseif dstmp <= 10f0                       # 5.0 ≤ DSTMP ≤ 10.0
            return Int(nint(-1f0 + 0.4f0 * dstmp)) # NINT (ties away from zero)
        else
            return 3
        end
    else
        return 1
    end
end

"SPRTHT SN species set that uses the `(0.1 + SI/50)·age` curve (essprt.f:1389)."
@inline _sprtht_sn_curve(ispc::Integer) =
    ispc == 5 || ispc == 15 || ispc == 16 ||
    (18 <= ispc <= 57) || (59 <= ispc <= 87)

"""
    sprtht_sn(ispc, si, iag) -> Float32

Sprout height (ft) at sprout age `iag` for site index `si` (`SPRTHT`,
essprt.f:1387-1393, SN). The sprouting hardwoods use `(0.1 + SI/50)·age`;
everything else falls back to the original NI regen rule `0.5 + 0.5·age`.
"""
@inline function sprtht_sn(ispc::Integer, si::Float32, iag::Real)::Float32
    a = Float32(iag)
    return _sprtht_sn_curve(ispc) ? (0.1f0 + si / 50f0) * a : 0.5f0 + 0.5f0 * a
end

"Special-establishment forests (R8/R9 NFs) that trigger the ESSPRT overrides."
@inline _es_special_forest(isefor::Integer) =
    isefor == 809 || isefor == 810 || isefor == 905 || isefor == 908

"""
    essprt_sn(coef, ispc, prem, dstmp, isefor) -> Float32

Apply the per-record sprout-survival multiplier to `prem` (carried sprout TPA)
for species `ispc` and stump diameter `dstmp` (`ESSPRT`, essprt.f:514-590, SN).

Most species use either a constant multiplier or a logistic in stump DBH,
`1/(1 + exp(-(a + b·DSTMP)))`; both forms (plus a per-species flag) are read
from `sprout_essprt.csv`. Five species (64/66/70/75/77) carry a distinct
special-forest variant (forests 809/810/905/908) handled explicitly here; for
all other forests their CSV row already holds the common-forest (ELSE) form.
"""
function essprt_sn(coef::SpeciesCoefficients, ispc::Integer, prem::Float32,
                   dstmp::Float32, isefor::Integer)::Float32
    if coef_col(coef, :essprt_fsp)[ispc] == 1f0 && _es_special_forest(isefor)
        d = dstmp
        m = if ispc == 64 || ispc == 66 || ispc == 75
                (57.3f0 - 0.0032f0 * d^3) / 100f0          # essprt.f:547/554/571
            elseif ispc == 70
                1f0 / (1f0 + exp(-(2.3656f0 - 0.2781f0 * (d / 0.7801f0))))  # :561
            else # ispc == 77
                1f0 / (1f0 + exp(-(-2.8058f0 + 22.6839f0 *
                                    (1f0 / ((d / 0.7788f0) - 0.4403f0)))))  # :578
            end
        return prem * Float32(m)
    end
    kind = coef_col(coef, :essprt_kind)[ispc]
    p1 = coef_col(coef, :essprt_p1)[ispc]
    if kind == 0f0
        return prem * p1                                    # constant multiplier
    end
    p2 = coef_col(coef, :essprt_p2)[ispc]
    return prem * (1f0 / (1f0 + exp(-(p1 + p2 * dstmp))))   # logistic in DSTMP
end
