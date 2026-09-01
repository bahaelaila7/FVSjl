# =============================================================================
# PPE (Parallel Processing Extension) — MXHRVP multistand harvest scheduling.
# =============================================================================
# Faithful port of the PPBASE "multistand treatment scheduling" subsystem
# (recovered source scratchpad/ppe/recovered/ppbase/hv*.f, from FVS rev bc6e2377^;
# every `HV`-prefixed name = multistand policy simulation). A landscape harvest-flow
# allocator: once per master cycle, for each user-defined harvest *policy* it selects
# which stands to cut so a landscape TARGET is met, cutting stands in PRIORITY order.
#
# This file ports the DETERMINISTIC allocation kernels; each is validated bit-exact
# against a gfortran-16 driver-golden over the pristine Fortran (the wwpb/bmatct recipe:
# scratchpad/ppe/mxhrvp/driver_hvsel.f linked over the historical FVSppe objects).
# The `MXHRVP=10` policy cap and the status-code legend are from PPHVCM.F77.
#
# Status codes (IHVSTA, per stand × policy — PPHVCM.F77):
#   0  = not a candidate for this policy
#  ±1  = candidate, not (yet) selected  (sign = would-be clearcut; +/- only live
#        between HVHRV1 and inside HVSEL)
#   2  = not needed to reach the target (its thinning yield, if any, still counts)
#  ±3  = selected (fully); all signs are positive after HVSEL
#   4  = partially selected — replicate the stand and split its sampling weight by
#        HVPART (only one such stand per policy)
# =============================================================================

"""
    hvsel!(status, priority, yield_sel, yield_notsel, target; lprtct=false) -> (status, hvpart)

Port of `HVSEL` (hvsel.f) for one harvest policy — the greedy, priority-ranked,
target-constrained stand selection, with the constraint options (max-contiguous-clearcut
`LHVMXC`, coordinated units `LHVUNT`, policy hierarchy `LHIER`, external selection
`IHVEXT`) OFF. In place, it rewrites `status` (the policy's `IHVSTA` column) from the
candidate coding (±1) to the final selection coding (2 / 3 / 4) and returns `hvpart`
(the partial-cut fraction of the one status-4 stand, 0 if none).

Arguments (all per-stand vectors of length `NOSTND`, `Float32` for bit-exactness):
- `status`     — `IHVSTA(:,policy)`: ±1 candidate / 0 non-candidate on entry.
- `priority`   — `HVPRI(:,policy)`: the ranking key (higher = cut first).
- `yield_sel`  — `HVYLDS(:,policy)`: resource a stand contributes IF selected.
- `yield_notsel` — `HVTHIN(:,policy)`: resource IF NOT selected (its scheduled thinning).
- `target`     — `TRGETS(policy)`: the landscape resource target for the master cycle.
- `lprtct`     — `LPRTCT`: allow partially cutting the marginal stand (the EXACT option).

Mirrors hvsel.f:79-276 (single active policy, options off): sort by priority descending
(`_rdpsrt!` = RDPSRT), sum the thinning yield, short-circuit if the target is already met,
else walk stands in priority order accumulating harvest yield until the target is reached.
"""
function hvsel!(status::AbstractVector{<:Integer},
                priority::AbstractVector{Float32},
                yield_sel::AbstractVector{Float32},
                yield_notsel::AbstractVector{Float32},
                target::Float32;
                lprtct::Bool = false)
    n = length(status)
    hvpart = 0.0f0
    n == 0 && return status, hvpart

    # HVSEL:95 — sort the stand priorities into descending order (RDPSRT, LODSRT=.TRUE.).
    # _rdpsrt! reads `priority` and fills `isnsrt` with the descending order; `priority`
    # is untouched. (LODSRT=.TRUE. ⇒ lseq=true initializes the index.)
    isnsrt = Vector{Int32}(undef, n)
    _rdpsrt!(priority, isnsrt; lseq = true)

    # HVSEL:100-105 — sum the yield due to thinning (status-1 candidates only). Float32.
    thnyld = 0.0f0
    @inbounds for ist in 1:n
        if status[ist] == 1
            thnyld += yield_notsel[ist]
        end
    end

    # HVSEL:110-113 — if the target is already reached (or ≤0), nothing is CUT: every
    # status-±1 candidate drops to 2 (its thinning still counts downstream).
    if thnyld >= target || target <= 0.0f0
        @inbounds for ist in 1:n
            if abs(status[ist]) == 1
                status[ist] = 2
            end
        end
    else
        # HVSEL:118-269 — walk the stands in descending-priority order, accumulating
        # harvest yield until the target is met.
        hrvyld = 0.0f0
        lmore = true
        @inbounds for ii in 1:n
            ist = Int(isnsrt[ii])
            # HVSEL:126 — skip stands whose status is not ±1.
            abs(status[ist]) != 1 && continue
            # (LHIER hierarchy / LHVUNT coordinated-unit branches omitted — options off.)
            if lmore
                # HVSEL:185-188 — if the selected-yield is ~0, don't cut it: mark 2.
                if yield_sel[ist] < 0.000001f0
                    status[ist] = 2
                    continue
                end
                # (LHVMXC max-contiguous-clearcut branch omitted — option off.)
                # HVSEL:209 — proportion of THIS stand needed to meet the target.
                pneed = (target - (thnyld + hrvyld)) / yield_sel[ist]
                if pneed >= 0.999f0
                    # HVSEL:217-224 — need (at least) the whole stand: select it (±3),
                    # add its harvest yield, remove its thinning yield from the pool.
                    status[ist] = 3 * status[ist]
                    hrvyld += yield_sel[ist]
                    thnyld -= yield_notsel[ist]
                else
                    # HVSEL:238-256 — only part of this (marginal) stand is needed.
                    if lprtct
                        if pneed > 0.01f0
                            status[ist] = 4                       # partial cut + thin
                            hrvyld += yield_sel[ist] * pneed
                            hvpart = pneed
                        else
                            status[ist] = 2                       # too small — thin it
                        end
                    else
                        if pneed > 0.5f0
                            status[ist] = 3 * status[ist]         # round up: select
                            hrvyld += yield_sel[ist]
                            thnyld -= yield_notsel[ist]
                        else
                            status[ist] = 2                       # round down: don't
                        end
                    end
                    lmore = false                                 # target reached
                end
            else
                # HVSEL:262 — past the target: this stand is not needed.
                status[ist] = 2
            end
        end
    end

    # HVSEL:274-276 — the sign bit's job is done; make all status codes positive.
    @inbounds for ist in 1:n
        status[ist] = abs(status[ist])
    end
    return status, hvpart
end
