# =============================================================================
# volume.jl (southeastalaska) — AK Region-10 volume. Chunk 8 — NOT YET PORTED.
#
# AK volume runs the standard VOLEQDEF(VAR='AK', IREGN=10) → VEQNNC → National Volume Estimator
# Library (NVEL) path plus AK-local taper (ak/logs.f cubic/board via cubrds.f). That NVEL crosswalk
# is a substantial separate chunk (same class as the other western variants' volume chunks). Until it
# lands, this fills the per-tree volume fields with 0 so the growth/mortality pipeline runs end-to-end
# and the .sum TopHt/QMD/TPA/BA columns are validatable; the cuft/bdft columns are NOT yet meaningful.
# =============================================================================

function compute_volumes_ak!(s::StandState)
    t = s.trees
    @inbounds for i in 1:(t.n + t.ndead)
        t.cuft_vol[i] = 0f0
        t.merch_cuft_vol[i] = 0f0
        t.saw_cuft_vol[i] = 0f0
        t.bdft_vol[i] = 0f0
    end
    return s
end
