# =============================================================================
# small_tree_growth.jl (olympic) — OP small-tree (REGENT) growth hook.
#
# op/regent.f grows trees below breast height / below the small-tree DBH threshold (height-first
# small-tree model). For the S248112 reference stand the inventory is all large trees (HT>4.5, the
# ORGANON/Wykoff large-tree paths), so REGENT overrides nothing here. This hook is a FAITHFUL no-op
# GUARDED on the measured absence of sub-4.5' live records: if any appear it warns (once) so the gap
# is visible rather than silently mis-grown. Port op/regent.f before running a small-tree OP stand.
# =============================================================================

const _OP_SMALLTREE_WARNED = Ref(false)

function small_tree_growth!(s::StandState, stash, ::Olympic; fint::Float32 = 5.0f0)
    t = s.trees
    nsmall = 0
    @inbounds for i in 1:t.n
        (t.tpa[i] > 0f0 && t.height[i] <= 4.5f0) && (nsmall += 1)
    end
    if nsmall > 0 && !_OP_SMALLTREE_WARNED[]
        @warn "OP small_tree_growth!(::Olympic) is a no-op but the stand has $nsmall sub-4.5' live " *
              "record(s); op/regent.f small-tree growth is UNPORTED — multi-cycle may drift." maxlog=1
        _OP_SMALLTREE_WARNED[] = true
    end
    return s
end
