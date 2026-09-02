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
                lprtct::Bool = false,
                lhvmxc::Bool = false, hvmxcc::Float32 = 0.0f0,
                areas::Union{Nothing,AbstractVector{Float32}} = nothing,
                border = nothing,
                ihvext::Bool = false)
    n = length(status)
    hvpart = 0.0f0
    # hvin.f:153-157 — external selection (IHVEXT=1) forces the EXACT partial-cut off.
    ihvext && (lprtct = false)
    n == 0 && return status, hvpart

    # HVSEL:95 — sort the stand priorities into descending order (RDPSRT, LODSRT=.TRUE.).
    # _rdpsrt! reads `priority` and fills `isnsrt` with the descending order; `priority`
    # is untouched. (LODSRT=.TRUE. ⇒ lseq=true initializes the index.)
    isnsrt = Vector{Int32}(undef, n)
    _rdpsrt!(priority, isnsrt; lseq = true)

    # HVSEL:100-105 — sum the yield due to thinning (status-1 candidates only). Float32.
    thnyld = 0.0f0
    hrvyld = 0.0f0
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
        lmore = true
        @inbounds for ii in 1:n
            ist = Int(isnsrt[ii])
            # HVSEL:126 — skip stands whose status is not ±1.
            abs(status[ist]) != 1 && continue
            # (LHIER hierarchy / LHVUNT coordinated-unit branches omitted — options off.)
            # HVSEL:167-175 — external selection (IHVEXT=1): never select a stand whose
            # externally-supplied priority is ≤ 0 (mark 2, skip to the next stand).
            if ihvext && priority[ist] <= 0.0f0
                status[ist] = 2
                continue
            end
            if lmore
                # HVSEL:185-188 — if the selected-yield is ~0, don't cut it: mark 2.
                if yield_sel[ist] < 0.000001f0
                    status[ist] = 2
                    continue
                end
                # HVSEL:189-205 — max-contiguous-clearcut constraint. If selecting this
                # (would-be-clearcut, status -1) stand would create more than HVMXCC
                # contiguous clear-cut acres with the already-selected clearcut stands,
                # do NOT select it: mark 2 and skip. Non-clearcut candidates (+1) return
                # contig 0 from hvcntg and are never vetoed.
                if lhvmxc
                    contig = hvcntg(ii, isnsrt, status, areas, border)
                    if contig > hvmxcc
                        status[ist] = 2
                        continue
                    end
                end
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
    # `hrvyld`/`thnyld` are the RUNNING report values (hvsel.f:352-354): the selected-
    # resource total and the residual non-selected supply (the Float32 subtraction residual
    # matters — the oracle prints e.g. 1.5e-5, not 0). Returned after (status, hvpart) so the
    # 2-value `st, hvpart = hvsel!(...)` call form stays valid.
    return status, hvpart, hrvyld, thnyld
end

"""
    spcntg!(isdwk1, areas, border) -> (contig::Float32, nneig::Int)

Port of `SPCNTG` (spcntg.f) — the contiguous clear-cut acreage of a clear-cut stand set.
`isdwk1::Vector{Int}` is stand indices, `isdwk1[1]` the SUBJECT stand; it is REORDERED in
place so the stands forming a contiguous region with the subject are at the top (`1:nneig`).
`areas[s]` = stand `s`'s area (acres); `border(i,j)::Float32` = the shared-border length of
stands `i,j` (≥ 0 if they share a border, < 0 if not — SPNBBD's `X.GE.0.0` test). Returns
the total contiguous area and `nneig`. Float32 area accumulation, bit-exact vs the Fortran
(driver-golden scratchpad/ppe/mxhrvp/driver_spcntg.f). Grows the region breadth-first from
the subject: each contiguous stand `i` pulls in any not-yet-added stand `j` sharing a border
(swapping it up into the contiguous block), until no frontier stand has an outside neighbor.
"""
function spcntg!(isdwk1::Vector{Int}, areas::AbstractVector{Float32}, border)
    nlst = length(isdwk1)
    nneig = 1
    contig = areas[isdwk1[1]]                       # SPLAAR(isdwk1(1))
    nlst == 1 && return contig, nneig
    @inbounds for ii in 1:nlst                      # DO 30 II=1,NLST
        inei = 0
        nneig >= nlst && break                      # GOTO 40
        i = isdwk1[ii]
        jj = nneig + 1                              # DO 20 JJ=NNEIG+1,NLST (start captured once)
        while jj <= nlst
            j = isdwk1[jj]
            if border(i, j) >= 0.0f0                # SPNBBD(i,j): shared border (IRC=0)
                inei = 1
                nneig += 1
                if jj > nneig
                    isdwk1[jj] = isdwk1[nneig]
                    isdwk1[nneig] = j
                end
                contig += areas[j]                  # SPLAAR(j) summed
            end
            jj += 1
        end
        inei == 0 && break                          # stand i had no outside neighbor → done
    end
    return contig, nneig
end

"""
    hvcntg(ii, isnsrt, status, areas, border) -> contig::Float32

Port of `HVCNTG` (hvcntg.f) — contiguous clear-cut acres created if the stand at
priority position `ii` (`isnsrt[ii]`) is selected. Returns 0 unless that subject stand is a
would-be-clearcut candidate (`status == -1`); otherwise builds the clear-cut set — the
subject plus every already-selected would-be-clearcut stand (`status == -3`) at higher
priority (positions `1:ii-1`) — and runs `spcntg!`.
"""
function hvcntg(ii::Integer, isnsrt::AbstractVector, status::AbstractVector{<:Integer},
                areas, border)
    ist = Int(isnsrt[ii])
    status[ist] != -1 && return 0.0f0
    isdwk1 = Int[ist]
    @inbounds for k in 1:(ii-1)
        s = Int(isnsrt[k])
        status[s] == -3 && push!(isdwk1, s)
    end
    contig, _ = spcntg!(isdwk1, areas, border)
    return contig
end

"""
    hvccut(tstv2_1, tstv2_7, dctpa, dctop) -> Int

Port of `HVCCUT` (hvccut.f) — the clearcut-condition test for the current stand.
Returns `-1` if the stand meets the DEFCCUT clearcut definition (after-thin trees/acre
`TSTV2(1) < DCTPA` **and** after-thin top height `TSTV2(7) < DCTOP`), else `+1`.
`HVHRV1` multiplies the stand's `IHVSTA` by this to set the ± sign bit that `hvsel!`
reads (sign = "would be in clearcut condition if selected"). `DCTPA`/`DCTOP` are the
`DEFCCUT` thresholds (default 0 ⇒ never a clearcut). Pure Float32 comparison — no
arithmetic, so bit-exact by construction (hvccut.f:36-37).
"""
hvccut(tstv2_1::Float32, tstv2_7::Float32, dctpa::Float32, dctop::Float32)::Int =
    (tstv2_1 < dctpa && tstv2_7 < dctop) ? -1 : 1

"""
    eval_policy_expr(expr, harvest; cycle=1, year=0, state=nothing) -> Float32

Compile-and-evaluate a MXHRVP policy expression (TARGET / PRIORITY / CREDIT) — the
`ALGCMP`+`ALGEVL` path (hvin.f/hvaloc.f) — by reusing the already-validated event-monitor
parser (`parse_event_condition`) and evaluator (`eval_event`). `harvest::HarvestVars`
supplies the PPE landscape/policy variables (AVB*/TOTALWT/OLDTARG/SELECTED); `state`
(a `StandState`, optional) supplies per-stand variables (BBA/AGE/…) for PRIORITY/CREDIT.
The evaluator's Float32 arithmetic is bit-exact to ALGEVL (validated by the event-monitor
tests); this only widens the variable set. In hot paths parse the expression once with
`parse_event_condition` and reuse the AST across stands via `eval_event(ast, ctx)`.
"""
function eval_policy_expr(expr::AbstractString, harvest::HarvestVars;
                          cycle::Integer = 1, year::Integer = 0, state = nothing)::Float32
    ast = parse_event_condition(expr)
    return eval_event(ast, EventCtx(Int(cycle), Int(year), state, harvest))
end

# --- label-set membership (base FVS lbmemr.f/lb1mem.f) ----------------------
# A stand's label set (SLSET) is a comma-space-delimited string of labels, e.g.
# "ALL, STAND1, GROUP2". A multistand policy's MSPLABEL (up to a '.') makes a stand a
# candidate iff it is a member of that stand's label set (hvaloc.f:143 CALL LBMEMR).

# lb1mem.f — length of the first comma-delimited member of setb[ip1:ip2] (1-based bytes),
# 0 if none. INDEX(set(ip1:ip2),',') → member is up to the comma; else the whole segment.
@inline function _lb1mem_len(setb, ip1::Int, ip2::Int)::Int
    (ip1 <= 0 || ip2 <= 0 || ip2 - ip1 < 0) && return 0
    @inbounds for k in ip1:ip2
        setb[k] == UInt8(',') && return k - ip1        # (lc-1): chars before the comma
    end
    return ip2 - ip1 + 1                                # no comma ⇒ whole segment
end

# lbmemr.f — true if the `lenmem`-byte member is one of the comma-space-delimited members
# of the `lenset`-byte set. Walks members (LB1MEM), comparing exactly on equal length; the
# +2 step past each member skips its trailing ", " (the delimiter is 2 chars).
function _lbmemr(memb, lenmem::Int, setb, lenset::Int)::Bool
    ip = 1
    @inbounds while ip <= lenset
        lenwrk = _lb1mem_len(setb, ip, lenset)
        if lenmem == lenwrk
            eq = true
            for k in 1:lenmem
                if memb[k] != setb[ip + k - 1]; eq = false; break; end
            end
            eq && return true
        end
        ip = ip + lenwrk + 2
    end
    return false
end

"""
    lbmemr(mem, set) -> Bool

Port of `LBMEMR` (base lbmemr.f) — true if label `mem` is a member of the comma-space
-delimited label `set` (e.g. `lbmemr("STAND1", "ALL, STAND1, GROUP2") == true`). Trailing
blanks are trimmed (LEN_TRIM); comparison is exact and length-gated. This is the MXHRVP
candidacy test (which stands a harvest policy's MSPLABEL applies to).
"""
function lbmemr(mem::AbstractString, set::AbstractString)::Bool
    m = String(rstrip(mem)); st = String(rstrip(set))
    return _lbmemr(codeunits(m), ncodeunits(m), codeunits(st), ncodeunits(st))
end

"""
    lbunin(set1, set2) -> (union::String, kode::Int)

Port of `LBUNIN` (base lbunin.f) — the label-set union: start with `set1`, append each
member of `set2` not already present (`lbmemr` de-dup), comma-space delimited, capped at
250 chars (`kode=1` if truncated). Empty/undefined operands short-circuit to the other set.
This is the primitive `LBDSET` uses to build a stand's default label set as the union over
its activity-group labels (the union-over-activity-groups default). Byte-exact / ASCII.
"""
function lbunin(set1::AbstractString, set2::AbstractString)
    s1 = String(rstrip(set1)); s2 = String(rstrip(set2))
    len1 = ncodeunits(s1); len2 = ncodeunits(s2)
    len1 <= 0 && return (len2 <= 0 ? "" : s2, 0)
    len2 <= 0 && return (s1, 0)
    union = s1; lnunin = len1; kode = 0
    s2b = codeunits(s2)
    ip = 1
    while ip <= len2
        lenwrk = _lb1mem_len(s2b, ip, len2)
        member = lenwrk > 0 ? s2[ip:ip+lenwrk-1] : ""
        if !lbmemr(member, union)
            if lnunin + lenwrk + 2 <= 250
                union = union * ", " * member
                lnunin = lnunin + lenwrk + 2
            else
                kode = 1
            end
        end
        ip = ip + lenwrk + 2
    end
    return (union, kode)
end

# --- the MXHRVP master-cycle coordinator (hvaloc.f) ------------------------
# MSPLABEL up to the '.' — the candidate-pool key (the part after '.' ties coordinated
# policies to one pool; hvin.f LNHPLB(.,2)).
_mslabel_base(lab::AbstractString) =
    (i = findfirst('.', lab); i === nothing ? String(strip(lab)) : String(strip(lab[1:prevind(lab, i)])))

"""
    ppe_read_rdaccess(path) -> Dict{String,Float32}

Port of the IHVEXT read-back — `SPRDRD` + `SPRDIS` (sprdrd.f/sprdis.f) reading the external
selector's `PPE_FFERdAccess.txt` into a stand-id → priority map (which HVALOC assigns to
`HVPRI`). Format `(A26,T30,F10.0)`: the stand id is columns 1-26 (trimmed), the priority
value columns 30-39. A `-999` in the id column terminates; a blank id line is skipped; a
blank value field maps to `-99999.0` (SPRDRD's missing-value sentinel). Bit-exact vs the
gfortran-16 driver-golden scratchpad/ppe/mxhrvp/driver_sprd.f (the real SPRDRD/SPRDIS).
"""
function ppe_read_rdaccess(path::AbstractString)
    d = Dict{String,Float32}()
    for raw in eachline(path)
        line = rpad(raw, 39)                       # pad so the fixed columns exist
        id26 = line[1:26]
        occursin("-999", id26) && break            # '-999' terminates (sprdrd.f:20)
        id = strip(id26)
        isempty(id) && continue                    # blank id skipped
        valfield = strip(line[30:39])              # T30,F10.0
        d[String(id)] = isempty(valfield) ? -99999.0f0 : parse(Float32, valfield)
    end
    return d
end

"""
    ppe_run_landscape_harvest!(stands; variant, labels, mslabel, target_expr, priority_expr,
                               credit_expr, master_years, period=5, lprtct=false) -> Vector

The MXHRVP master-cycle coordinator (hvaloc.f) composed from the validated kernels. Per
master-cycle year it selects which member stands to cut so the landscape `target_expr` is
met, in `priority_expr` order — faithful to hvaloc.f with the constraint options OFF
(LHVMXC max-contiguous-clearcut, LHVUNT coordinated units, LHIER hierarchy, IHVEXT external
selection — deferred/documented).

Each member `stands[i]` is a keyfile with area weight `.area`; `labels[i]` is its stand
label set (SPLABEL). For each stand it projects once (the already-validated engine via
`write_sum_file`) and, at each master-cycle year, runs the HVTHN1/HVHRV1 trial-cut — the
stand's own scheduled cuts under `SELECTED=0` (→ HVTHIN, yield-if-not-selected) and
`SELECTED=1` (→ HVYLDS, yield-if-selected) — then evaluates `credit_expr` and `priority_expr`
(`eval_policy_expr`, the ALGEVL path). Candidacy is `lbmemr(mslabel, labels[i])`. Then per
year: accumulate the area-weighted landscape stats, evaluate `target_expr`, and call `hvsel!`.

Returns one NamedTuple per master year: `(year, target, stand_ids, priority, credit,
selected, selected_resource, nonselected_supply, pct_of_target)` — the HVSEL selection table.
"""
function ppe_run_landscape_harvest!(stands::AbstractVector{PPEStand}; variant,
        labels::AbstractVector{<:AbstractString}, mslabel::AbstractString,
        target_expr::AbstractString, priority_expr::AbstractString, credit_expr::AbstractString,
        master_years::AbstractVector{<:Integer}, period::Integer = 5, lprtct::Bool = false,
        rdaccess::Union{Nothing,AbstractString} = nothing,
        faithful::Bool = true)
    # IHVEXT=1 external harvest selection (hvaloc.f/hvsel.f): when a road-access file is
    # supplied, the per-stand PRIORITY is OVERRIDDEN by the external selector's values
    # (SPRDRD/SPRDIS read of PPE_FFERdAccess.txt → HVPRI) and HVSEL runs in external mode
    # (never selects a stand whose supplied priority is ≤ 0). The read-back file is the
    # USER-approved staged read: stage the file the external selector would emit.
    rdmap = rdaccess === nothing ? nothing : ppe_read_rdaccess(rdaccess)
    n = length(stands)
    length(labels) == n || error("ppe_run_landscape_harvest!: labels/stands length mismatch")
    yrs = Int.(master_years); yrset = Set(yrs)
    pri  = [Dict{Int,Float32}() for _ in 1:n]     # PRIORITY per master year
    ythn = [Dict{Int,Float32}() for _ in 1:n]     # HVTHIN (SELECTED=0)
    yhrv = [Dict{Int,Float32}() for _ in 1:n]     # HVYLDS (SELECTED=1)
    ids  = Vector{String}(undef, n)
    for (i, ps) in enumerate(stands)
        got = false
        for s in each_stand(ps.keyfile; variant = variant, faithful = faithful)
            notre!(s); setup_growth!(s); compute_volumes!(s)
            ids[i] = String(strip(s.plot.stand_id)); got = true
            # HVTHN1/HVHRV1: the stand's own scheduled cuts run under SELECTED=0 / SELECTED=1,
            # then CREDIT is evaluated. A stand with no SELECTED-gated cut ⇒ CUTS is a no-op and
            # both reduce to CREDIT(state); a SELECTED-gated cut makes HVYLDS ≠ HVTHIN.
            hook = (st, yr, _pl, cy) -> begin
                y = Int(yr)
                if y in yrset
                    pri[i][y]  = eval_policy_expr(priority_expr, HarvestVars(); state = st, cycle = Int(cy), year = y)
                    ythn[i][y] = eval_policy_expr(credit_expr, HarvestVars(; selected = 0f0); state = st, cycle = Int(cy), year = y)
                    yhrv[i][y] = eval_policy_expr(credit_expr, HarvestVars(; selected = 1f0); state = st, cycle = Int(cy), year = y)
                end
            end
            write_sum_file(IOBuffer(), s; period = Int(period), stand_id = ids[i], mgmt_id = "NONE",
                           variant = variant_code(s.variant), date = "x", time = "y", cycle_hook = hook)
            break
        end
        got || error("ppe_run_landscape_harvest!: keyfile $(ps.keyfile) produced no stand")
    end
    base = _mslabel_base(mslabel)
    out = NamedTuple[]
    for y in yrs
        cand = Int[i for i in 1:n if lbmemr(base, labels[i]) && haskey(yhrv[i], y)]
        isempty(cand) && continue
        prio = Float32[pri[c][y]  for c in cand]
        ysel = Float32[yhrv[c][y] for c in cand]
        ynot = Float32[ythn[c][y] for c in cand]
        # IHVEXT=1: override the computed priorities with the external selector's values
        # (by stand id), keyed exactly as SPRDIS looks them up.
        if rdmap !== nothing
            for (k, c) in enumerate(cand)
                prio[k] = get(rdmap, strip(ids[c]), 0.0f0)
            end
        end
        status = fill(1, length(cand))
        # area-weighted landscape stats (TRGSTS/PTSTV1) for a stats-based TARGET; a constant
        # TARGET (msp.key) ignores them. AVBBA is the weighted mean priority proxy here.
        totalwt = sum(Float32(stands[c].area) for c in cand)
        avbba = totalwt > 0f0 ?
            sum(Float32(stands[c].area) * prio[k] for (k, c) in enumerate(cand)) / totalwt : 0f0
        hv = HarvestVars(; avbba = avbba, totalwt = totalwt)
        target = eval_policy_expr(target_expr, hv; year = y)
        st2, hvpart, hrvyld, thnyld = hvsel!(status, prio, ysel, ynot, target;
                                             lprtct = lprtct, ihvext = rdmap !== nothing)
        sel = Bool[st2[k] != 2 for k in eachindex(cand)]          # status 3/4 = selected, 2 = not
        partial = Int[st2[k] == 4 ? k : 0 for k in eachindex(cand)]  # the one status-4 stand (0 = none)
        pct = target > 0f0 ? (hrvyld + thnyld) / target * 100f0 : 0f0
        push!(out, (year = y, target = target, stand_ids = String[ids[c] for c in cand],
                    priority = prio, credit_sel = ysel, credit_notsel = ynot,
                    status = copy(st2), selected = sel, hvpart = hvpart,
                    selected_resource = hrvyld, nonselected_supply = thnyld, pct_of_target = pct))
    end
    return out
end
