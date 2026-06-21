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

"Per-acre removed totals from a thin (zeros when nothing is cut)."
const _NO_REMOVAL = (tpa = 0f0, cuft = 0f0, mcuft = 0f0, scuft = 0f0, bdft = 0f0)

"""
    cuts!(state; fint) -> removed

Run any thinning/harvest scheduled for the current cycle's year. Reduces
`trees.tpa` (PROB) in place and returns the period's removed totals (TPA + the four
volumes, summed over the cut). Call at the top of `grow_cycle!`, before growth.
"""
function cuts!(s::StandState; fint::Float32 = 5f0)
    sched = s.control.schedule
    isempty(sched) && return _NO_REMOVAL
    # Year of the current cycle (matches summary_row): inventory year + cycle·period.
    # (cycle_year only stores the inventory year; later years are derived.)
    yr = Int32(Int(s.control.cycle_year[1]) + Int(s.control.cycle) * round(Int, fint))
    yr in s.control.years_cut && return _NO_REMOVAL   # idempotent: already cut this year
    rem = _NO_REMOVAL
    applied = false
    @inbounds for act in sched
        act.year == yr || continue
        applied = true
        ic = act.icflag
        r = ic == Int32(8) ? _thindbh!(s, act) :                       # proportional DBH-class
            (ic in (Int32(3), Int32(4), Int32(5), Int32(6))) ? _thin_sorted!(s, act) :  # BTA/ATA/BBA/ABA
            _NO_REMOVAL
        rem = (tpa = rem.tpa + r.tpa, cuft = rem.cuft + r.cuft,
               mcuft = rem.mcuft + r.mcuft, scuft = rem.scuft + r.scuft,
               bdft = rem.bdft + r.bdft)
    end
    applied && push!(s.control.years_cut, yr)
    return rem
end

# CLSSTK (cutstk.f): TPA (jtyp=1) or basal-area (jtyp=2) stocking over the trees in
# the DBH/HT/species class, using the pre-thin TPA copy `wk4`.
@inline function _clsstk(t::TreeList, wk4, n, jtyp, ispcut, dl, du, hl=0f0, hu=999f0)
    cstock = 0f0
    @inbounds for i in 1:n
        _cut_eligible(t, i, ispcut, dl, du, hl, hu) || continue
        cstock += jtyp == 2 ? wk4[i] * t.dbh[i]^2 * _BA_PER_TREE : wk4[i]
    end
    return cstock
end

# WK2 eligibility (cuts.f:611-648): DBH in [dl,du) and HT in [hl,hu) and species
# included. (Group/LEAVESP handling is added when those keywords are ported.)
@inline function _cut_eligible(t::TreeList, i, ispcut, dl, du, hl=0f0, hu=999f0)
    d = t.dbh[i]
    (d < dl || d >= du) && return false
    h = t.height[i]
    (h < hl || h >= hu) && return false
    return ispcut == 0 || ispcut == t.species[i]
end

# THINDBH (cuts.f label_325→355→550→1100): thin the DBH/species class to a residual
# TPA (ctpa) and/or basal area (cba). cuteff = remove/cstock is applied to each
# eligible record's TPA in record order until the removal budget is spent (the last
# record is partial); the residual replaces PROB.
function _thindbh!(s::StandState, act::ScheduledActivity)
    t = s.trees; n = t.n
    n == 0 && return _NO_REMOVAL
    dbhlo, dbhhi, _eff, spf, ctpa_in, cba_in = act.params
    ispcut = floor(Int32, spf)
    valmin = dbhlo
    valmax = dbhhi < dbhlo ? 999f0 : dbhhi
    ctpa = max(0f0, ctpa_in); cba = max(0f0, cba_in)
    lbarea = cba > 0f0
    jtyp = lbarea ? 2 : 1

    wk4 = Float32[t.tpa[i] for i in 1:n]               # pre-thin PROB copy
    cstock = _clsstk(t, wk4, n, jtyp, ispcut, valmin, valmax)
    cstock <= 0f0 && return _NO_REMOVAL
    remove = cstock - (ctpa + cba)
    remove <= 0f0 && return _NO_REMOVAL
    cuteff = min(1f0, remove / cstock)

    # removed totals (per-acre, gross): Σ prem·{1, CFV, MCFV, SCFV, BFV}
    rtpa = 0f0; rcuft = 0f0; rmcuft = 0f0; rscuft = 0f0; rbdft = 0f0
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
        rtpa += prem; rcuft += prem * t.cuft_vol[i]
        rmcuft += prem * t.merch_cuft_vol[i]; rscuft += prem * t.saw_cuft_vol[i]
        rbdft += prem * t.bdft_vol[i]
    end
    @inbounds for i in 1:n
        t.tpa[i] = wk4[i]                               # residual replaces PROB
    end
    return (tpa = rtpa, cuft = rcuft, mcuft = rmcuft, scuft = rscuft, bdft = rbdft)
end

# THINBTA/ATA/BBA/ABA (cuts.f label_200/225/250/275): thin the DBH/HT class from
# below (BTA/BBA) or above (ATA/ABA) to a residual TPA (BTA/ATA) or BA (BBA/ABA).
# Trees are ranked by size (RDPSRT, descending priority: −DBH from below ⇒ smallest
# first; +DBH from above ⇒ largest first) and whole records are removed (×cuteff)
# until the removal budget is spent (last record partial).
function _thin_sorted!(s::StandState, act::ScheduledActivity)
    t = s.trees; n = t.n
    n == 0 && return _NO_REMOVAL
    ic = act.icflag
    lbelow = ic == Int32(3) || ic == Int32(5)        # BTA, BBA from below
    lbarea = ic == Int32(5) || ic == Int32(6)        # BBA, ABA to basal area
    jtyp = lbarea ? 2 : 1
    target, cuteff_p, dbhlo, dbhhi_p, htlo, hthi_p = act.params
    cuteff = cuteff_p > 0f0 ? cuteff_p : 1f0
    target = max(0f0, target)
    valmax = dbhhi_p < dbhlo ? 9999f0 : dbhhi_p
    hthi   = hthi_p   < htlo  ? 9999f0 : hthi_p

    wk4 = Float32[t.tpa[i] for i in 1:n]
    cstock = _clsstk(t, wk4, n, jtyp, 0, dbhlo, valmax, htlo, hthi)
    cstock <= 0f0 && return _NO_REMOVAL
    remove = cstock - target
    remove <= 0f0 && return _NO_REMOVAL

    # priority key = ±DBH (size); eligible records ranked, ineligible pushed last.
    elig = falses(n); key = Vector{Float32}(undef, n)
    @inbounds for i in 1:n
        e = _cut_eligible(t, i, 0, dbhlo, valmax, htlo, hthi)
        elig[i] = e
        key[i] = e ? (lbelow ? -t.dbh[i] : t.dbh[i]) : -Inf32
    end
    order = sortperm(key; rev = true)               # descending (RDPSRT order)

    rtpa = 0f0; rcuft = 0f0; rmcuft = 0f0; rscuft = 0f0; rbdft = 0f0
    totcut = 0f0
    @inbounds for it in order
        elig[it] || continue
        d = t.dbh[it]
        prem = wk4[it] * cuteff
        prem > wk4[it] && (prem = wk4[it])
        prem <= 0f0 && continue
        cut_v = lbarea ? prem * d * d * _BA_PER_TREE : prem
        xleft = remove - (totcut + cut_v)
        if xleft < 0f0
            prem = ((xleft + cut_v) / cut_v) * prem
            cut_v = remove - totcut
        end
        totcut += cut_v
        wk4[it] -= prem
        rtpa += prem; rcuft += prem * t.cuft_vol[it]
        rmcuft += prem * t.merch_cuft_vol[it]; rscuft += prem * t.saw_cuft_vol[it]
        rbdft += prem * t.bdft_vol[it]
        totcut >= remove && break
    end
    @inbounds for i in 1:n
        t.tpa[i] = wk4[i]
    end
    return (tpa = rtpa, cuft = rcuft, mcuft = rmcuft, scuft = rscuft, bdft = rbdft)
end
