# =============================================================================
# cuts.jl — scheduled thinning / harvest (CUTS)
#
# Ported from: base/cuts.f + base/cutstk.f (CLSSTK).
#
# Applies the management activities parsed into `control.schedule` (THIN* keywords)
# at the start of the cycle whose calendar year matches, before growth — matching
# the Fortran GRINCR order (CUTS runs before DGDRIV). Milestone scope: THINDBH
# (icflag 8), the residual-stocking DBH-class thin. Other methods land next.
# =============================================================================

const _BA_PER_TREE = 0.005454154f0   # cuts.f basal-area-per-tree factor

"""
    cuts!(state) -> state

Run any thinning/harvest scheduled for the current cycle's year. Reduces
`trees.tpa` (PROB) in place. Call at the top of `grow_cycle!`, before growth.
"""
function cuts!(s::StandState; fint::Float32 = 5f0)
    sched = s.control.schedule
    isempty(sched) && return false
    # Year of the current cycle (matches summary_row): inventory year + cycle·period.
    # (cycle_year only stores the inventory year; later years are derived.)
    yr = Int32(Int(s.control.cycle_year[1]) + Int(s.control.cycle) * round(Int, fint))
    cut = false
    @inbounds for act in sched
        act.year == yr || continue
        act.icflag == Int32(8) && (cut |= _thindbh!(s, act))   # THINDBH; more methods later
    end
    return cut
end

# CLSSTK (cutstk.f): TPA (jtyp=1) or basal-area (jtyp=2) stocking over the trees in
# the DBH / species class, using the pre-thin TPA copy `wk4`.
@inline function _clsstk(t::TreeList, wk4, n, jtyp, ispcut, dl, du)
    cstock = 0f0
    @inbounds for i in 1:n
        _cut_eligible(t, i, ispcut, dl, du) || continue
        cstock += jtyp == 2 ? wk4[i] * t.dbh[i]^2 * _BA_PER_TREE : wk4[i]
    end
    return cstock
end

# WK2 eligibility for THINDBH/lspecl (cuts.f:611-628): DBH in [dl,du) and species
# included. (Group/LEAVESP handling is added when those keywords are ported.)
@inline function _cut_eligible(t::TreeList, i, ispcut, dl, du)
    d = t.dbh[i]
    (d < dl || d >= du) && return false
    return ispcut == 0 || ispcut == t.species[i]
end

# THINDBH (cuts.f label_325→355→550→1100): thin the DBH/species class to a residual
# TPA (ctpa) and/or basal area (cba). cuteff = remove/cstock is applied to each
# eligible record's TPA in record order until the removal budget is spent (the last
# record is partial); the residual replaces PROB.
function _thindbh!(s::StandState, act::ScheduledActivity)
    t = s.trees; n = t.n
    n == 0 && return
    dbhlo, dbhhi, _eff, spf, ctpa_in, cba_in = act.params
    ispcut = floor(Int32, spf)
    valmin = dbhlo
    valmax = dbhhi < dbhlo ? 999f0 : dbhhi
    ctpa = max(0f0, ctpa_in); cba = max(0f0, cba_in)
    lbarea = cba > 0f0
    jtyp = lbarea ? 2 : 1

    wk4 = Float32[t.tpa[i] for i in 1:n]               # pre-thin PROB copy
    cstock = _clsstk(t, wk4, n, jtyp, ispcut, valmin, valmax)
    cstock <= 0f0 && return false
    remove = cstock - (ctpa + cba)
    remove <= 0f0 && return false
    cuteff = min(1f0, remove / cstock)

    totcut = 0f0
    @inbounds for i in 1:n
        _cut_eligible(t, i, ispcut, valmin, valmax) || continue
        d = t.dbh[i]
        prem = wk4[i] * cuteff
        prem > wk4[i] && (prem = wk4[i])
        prem <= 0f0 && continue
        cut_v = lbarea ? prem * d * d * _BA_PER_TREE : prem
        xleft = remove - (totcut + cut_v)
        if xleft < 0f0                                  # last record: partial cut
            prem = ((xleft + cut_v) / cut_v) * prem
            cut_v = remove - totcut
        end
        totcut += cut_v
        wk4[i] -= prem
    end
    @inbounds for i in 1:n
        t.tpa[i] = wk4[i]                               # residual replaces PROB
    end
    return true
end
