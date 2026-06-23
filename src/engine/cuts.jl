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

# RDPSRT (rdpsrt.f): real descending indirect quicksort — fill `index[1..n]` with a
# permutation so that key[index[1]] >= key[index[2]] >= … A faithful port of the
# Fortran partition (NOT Julia's `sortperm`) because the TIE-BREAK among equal keys
# must match FVS: snt01 has identical input trees that triple into exact-DBH-tied
# records, so the from-below cut's tie order decides WHICH lineage is removed — and
# the surviving lineages set the post-thin DGSCOR traversal order (and thus the RNG).
function _rdpsrt!(key::AbstractVector{Float32}, index::AbstractVector{Int32})
    n = length(key)
    @inbounds for i in 1:n; index[i] = Int32(i); end
    n < 2 && return index
    ipush = zeros(Int, 64)
    itop = 0; il = 1; iu = n
    indil = indiu = indip = indkl = indku = 0
    ip = kl = ku = jl = ju = 0
    t = 0f0
    @inbounds while true
        @label l30
        if iu <= il; @goto l40; end
        indil = Int(index[il]); indiu = Int(index[iu])
        if iu > il + 1; @goto l50; end
        if key[indil] >= key[indiu]; @goto l40; end
        index[il] = Int32(indiu); index[iu] = Int32(indil)
        @label l40
        itop == 0 && break
        il = ipush[itop-1]; iu = ipush[itop]; itop -= 2
        @goto l30
        @label l50
        ip = (il + iu) ÷ 2
        indip = Int(index[ip]); t = key[indip]
        index[ip] = Int32(indil)
        kl = il; ku = iu
        @label l60
        kl += 1
        if kl > ku; @goto l90; end
        indkl = Int(index[kl])
        if key[indkl] >= t; @goto l60; end
        @label l70
        indku = Int(index[ku])
        if ku < kl; @goto l100; end
        if key[indku] > t; @goto l80; end
        ku -= 1
        @goto l70
        @label l80
        index[kl] = Int32(indku); index[ku] = Int32(indkl); ku -= 1
        @goto l60
        @label l90
        indku = Int(index[ku])
        @label l100
        index[il] = Int32(indku); index[ku] = Int32(indip)
        if ku <= ip; @goto l110; end
        jl = il; ju = ku - 1; il = ku + 1
        @goto l120
        @label l110
        jl = ku + 1; ju = iu; iu = ku - 1
        @label l120
        itop += 2
        ipush[itop-1] = jl; ipush[itop] = ju
    end
    return index
end

"Per-acre removed totals from a thin (zeros when nothing is cut)."
const _NO_REMOVAL = (tpa = 0f0, cuft = 0f0, mcuft = 0f0, scuft = 0f0, bdft = 0f0)

"""
    cuts!(state; fint) -> removed

Run any thinning/harvest scheduled for the current cycle's year. Reduces
`trees.tpa` (PROB) in place and returns the period's removed totals (TPA + the four
volumes, summed over the cut). Call at the top of `grow_cycle!`, before growth.
"""
# ESTUMP (estump.f): log one removed tree of a *sprouting* species for stump
# sprouting, in removal order, gated on LSPRUT. `prem` is the removed TPA (less
# standing snags, which FVSjl's basic cut path does not model ⇒ SSNG=0). The record
# carries the stump DBH (= the cut tree's DBH), the point, and the sprout age ISHAG
# (= IFINT, the cycle length). Only the 72 SN sprouting species (ISPSPE, is_sprouting
# CSV flag) are logged — ESTUMP returns immediately for any other species. Consumed
# by `esuckr!`; inert (write-only) until that loop is wired into the cycle (Chunk C2).
@inline function _log_cut!(s::StandState, t, i::Integer, prem::Float32)
    prem > 0f0 || return
    sp = Int(t.species[i])
    # ECON: value this removed tree (DBH-class cost/revenue) at the removal point, before
    # the tree list is compacted (eccalc.f/echarv.f). Accumulated for the cycle's harvest.
    if s.econ !== nothing && s.econ.active
        s.econ.cycle_cost += harvest_value(s.econ.hrv_cost, sp, t.dbh[i], prem, t.cuft_vol[i], t.bdft_vol[i])
        s.econ.cycle_rev  += harvest_value(s.econ.hrv_rev,  sp, t.dbh[i], prem, t.cuft_vol[i], t.bdft_vol[i])
    end
    # ESTUMP cut log (sprouting species only, when sprouting is on)
    (s.control.lsprut && coef_col(s.coef, :is_sprouting)[sp] == 1f0) || return
    push!(s.control.cut_log,
          (species = Int32(sp), dstmp = t.dbh[i], prem = prem,
           plot = Int32(t.plot_id[i]), ishag = round(Int32, s.plot.cycle_length)))
    return
end

function cuts!(s::StandState; fint::Float32 = 5f0)
    s.control.lsprut && (s.plot.cycle_length = fint)  # IFINT (FINT) — sprout age for ESTUMP/SPRTHT
    sched = s.control.schedule
    conds = s.control.conditionals
    (isempty(sched) && isempty(conds)) && return _NO_REMOVAL
    # Year of the current cycle (matches summary_row): inventory year + cycle·period.
    # (cycle_year only stores the inventory year; later years are derived.)
    yr = Int32(Int(s.control.cycle_year[1]) + Int(s.control.cycle) * round(Int, fint))
    yr in s.control.years_cut && return _NO_REMOVAL   # idempotent: already cut this year
    # Fresh ESTUMP cut log for this cycle (for ESUCKR). Cleared HERE — after the
    # idempotency guard — so the second (summary) cuts! call of the cycle, which
    # early-returns above, does not wipe the log the first call populated.
    s.control.lsprut && empty!(s.control.cut_log)
    # Effective activities this cycle: the dated ones (year==yr) plus any IF/THEN block
    # whose algebraic condition is true this cycle (EVMON), with their year set to yr.
    # THINAUTO (icflag 1) is a RECURRING auto-thin: once scheduled it re-evaluates the
    # AUTMAX stocking gate every cycle from its start year (cuts.f auto path), not once.
    # A date < 1000 is a cycle number (FVS 1-based), not a calendar year: it fires in the
    # cycle whose 1-based index equals the date (OPNEW/OPFIND date convention).
    fvscyc = Int(s.control.cycle) + 1
    acts = ScheduledActivity[a for a in sched
                             if a.year == yr || (a.icflag == Int32(1) && yr >= a.year) ||
                                (0 < Int(a.year) < 1000 && Int(a.year) == fvscyc)]
    if !isempty(conds)
        ctx = EventCtx(Int(s.control.cycle) + 1, Int(yr), s)   # FVS CYCLE is 1-based
        # COMPUTE: evaluate event-monitor user variables this cycle BEFORE the IF conditions
        # read them (defs are in declaration order, so a later one may use an earlier one).
        @inbounds for (cd, nm, ast) in s.control.compute_defs
            yr >= cd && (s.control.compute_vars[nm] = eval_event(ast, ctx))
        end
        for c in conds
            eval_event(c.cond, ctx) != 0f0 || continue
            for a in c.acts
                push!(acts, ScheduledActivity(yr, a.icflag, a.params))
            end
        end
    end
    isempty(acts) && return _NO_REMOVAL
    # PASS 1 — cut MODIFIERS for this year (set state the methods read), before any
    # method runs (cuts.f processes SPECPREF/MINHARV/… then the thin in the cycle).
    cc = s.control
    @inbounds for act in acts
        act.icflag == Int32(201) && _apply_specpref!(s, act)
        act.icflag == Int32(206) && _apply_spleave!(s, act)   # SPLEAVE: per-species leave flag
        if act.icflag == Int32(202)        # TCONDMLT (cuts.f:1424): TCWT·IMC + SPCLWT·ISPECL weights
            cc.total_wt = act.params[1]; cc.special_wt = act.params[2]
        end
        if act.icflag == Int32(200)        # MINHARV (cuts.f:400): set the harvest-minimum thresholds
            cc.ba_min = act.params[1]; cc.tcf_min = act.params[2]; cc.cf_min = act.params[3]
            cc.scf_min = act.params[4]; cc.bf_min = act.params[5]
        end
    end
    # MINHARV gate is live whenever any threshold is set (persists across cycles once set).
    minharv_on = cc.ba_min > 0f0 || cc.tcf_min > 0f0 || cc.cf_min > 0f0 || cc.scf_min > 0f0 || cc.bf_min > 0f0
    tpa_snap = minharv_on ? copy(@view s.trees.tpa[1:s.trees.n]) : Float32[]
    # SETPTHIN (icflag 248) prescription this cycle → (point, metric) read by THINPT.
    # (same-cycle prescription; cross-cycle persistence would need control state.)
    pt_point = Int32(0); pt_metric = Int32(0); pt_set = false
    @inbounds for act in acts
        if act.icflag == Int32(248)
            pt_point  = Int32(round(act.params[1]))
            pt_metric = Int32(round(act.params[2]))
            pt_set = true
        end
    end
    # PASS 2 — cut METHODS.
    rem = _NO_REMOVAL
    applied = false
    @inbounds for act in acts
        ic = act.icflag
        # only CUTS methods here; establishment (427/430/431), the SPECPREF (201) and
        # SETPTHIN (248) modifiers are consumed elsewhere.
        ic in (Int32(3), Int32(4), Int32(5), Int32(6), Int32(7), Int32(8), Int32(10), Int32(12), Int32(14), Int32(1), Int32(11), Int32(17), Int32(15)) || continue
        applied = true
        r = (ic == Int32(8) || ic == Int32(12)) ? _thindbh!(s, act) : # DBH-class / HT-class residual
            ic == Int32(7)  ? _thinprsc!(s, act) :                     # prescription (cut-code marked)
            ic == Int32(10) ? _thin_sdi!(s, act) :                     # THINSDI (Zeide target SDI)
            ic == Int32(14) ? _thin_rden!(s, act) :                    # THINRDEN (Curtis RD)
            ic == Int32(1)  ? _thin_auto!(s, act) :                    # THINAUTO (auto to normal stocking)
            ic == Int32(11) ? _thin_cc!(s, act) :                      # THINCC (crown cover)
            ic == Int32(17) ? _thin_qfa!(s, act) :                     # THINQFA (Q-factor distribution)
            ic == Int32(15) ? (pt_set ? _thin_pt!(s, act, pt_point, pt_metric) : _NO_REMOVAL) :  # THINPT (point)
            (ic in (Int32(3), Int32(4), Int32(5), Int32(6))) ? _thin_sorted!(s, act) :  # BTA/ATA/BBA/ABA
            _NO_REMOVAL
        rem = (tpa = rem.tpa + r.tpa, cuft = rem.cuft + r.cuft,
               mcuft = rem.mcuft + r.mcuft, scuft = rem.scuft + r.scuft,
               bdft = rem.bdft + r.bdft)
    end
    # MINHARV (cuts.f:1556): if the cycle's total removal falls below ANY harvest minimum
    # (BA / total / merch / sawlog cubic / board feet), the whole cut is CANCELED — restore the
    # pre-thin TPA and report no removal. Default thresholds are 0, so the gate is a no-op then.
    if applied && minharv_on
        t = s.trees
        ba_rem = 0f0
        @inbounds for i in 1:length(tpa_snap)
            ba_rem += (tpa_snap[i] - t.tpa[i]) * t.dbh[i]^2 * _BA_PER_TREE
        end
        if !(ba_rem >= cc.ba_min && rem.cuft >= cc.tcf_min && rem.mcuft >= cc.cf_min &&
             rem.scuft >= cc.scf_min && rem.bdft >= cc.bf_min)
            @inbounds for i in 1:length(tpa_snap); t.tpa[i] = tpa_snap[i]; end
            return _NO_REMOVAL
        end
    end
    if applied
        push!(s.control.years_cut, yr)
        rem.tpa > 0f0 && tredel_compact!(s.trees)   # TREDEL: swap-from-end (oracle's exact post-thin layout)
    end
    return rem
end

# SPECPREF (cuts.f label_1200): set the per-species cut preference IORDER (added to
# the RDPSRT priority, so a higher value ⇒ that species is removed first). params[1]
# = species (>0 single, 0 all, <0 group via SPGROUP — group path pending), params[2]
# = preference value.
function _apply_specpref!(s::StandState, act::ScheduledActivity)
    isp = floor(Int, act.params[1]); pref = floor(Int32, act.params[2])
    pref_v = s.control.cut_pref
    if isp == 0
        fill!(pref_v, pref)
    elseif isp > 0
        1 <= isp <= MAXSP && (pref_v[isp] = pref)
    elseif -isp <= length(s.control.sp_groups)        # SPGROUP group −isp (cuts.f label_1200)
        @inbounds for sp in s.control.sp_groups[-isp]
            1 <= sp <= MAXSP && (pref_v[sp] = pref)
        end
    end
    return
end

# SPLEAVE (cuts.f:1466): set the per-species LEAVESP flag — a "left" species is excluded from
# every thin. params[1] = species (<0 SPGROUP group / 0 reset ALL to false / >0 single), params[2]
# = flag (>0 ⇒ leave, else don't leave).
function _apply_spleave!(s::StandState, act::ScheduledActivity)
    lv = s.control.leave_species
    isp = round(Int, act.params[1]); flag = act.params[2] >= 0.5f0
    if isp < 0
        g = -isp
        (1 <= g <= length(s.control.sp_groups)) || return
        @inbounds for sp in s.control.sp_groups[g]; 1 <= sp <= length(lv) && (lv[sp] = flag); end
    elseif isp == 0
        fill!(lv, false)
    else
        1 <= isp <= length(lv) && (lv[isp] = flag)
    end
    return
end

# Empty species-group table — the default for cuts that don't filter by a SPGROUP group
# (a const sentinel, so the common path allocates nothing).
const _NOGROUPS = Vector{Int32}[]

# Does species `sp` match a cut's species field `ispcut`? 0 = all, >0 = that species,
# <0 = every member of SPGROUP group −ispcut (cuts.f resolves group names to −index).
@inline function _sp_in_cut(sp_groups::Vector{Vector{Int32}}, ispcut, sp)
    ispcut == 0 && return true
    ispcut > 0 && return ispcut == sp
    g = -ispcut
    return g <= length(sp_groups) && (Int32(sp) in sp_groups[g])
end

# CLSSTK (cutstk.f): TPA (jtyp=1) or basal-area (jtyp=2) stocking over the trees in
# the DBH/HT/species class, using the pre-thin TPA copy `wk4`.
@inline function _clsstk(s::StandState, wk4, n, jtyp, ispcut, dl, du, hl=0f0, hu=999f0,
                         sp_groups::Vector{Vector{Int32}}=_NOGROUPS)
    t = s.trees; cstock = 0f0
    @inbounds for i in 1:n
        _cut_eligible(s, i, ispcut, dl, du, hl, hu, sp_groups) || continue
        cstock += jtyp == 2 ? wk4[i] * t.dbh[i]^2 * _BA_PER_TREE : wk4[i]
    end
    return cstock
end

# WK2 eligibility (cuts.f:611-648): DBH in [dl,du) and HT in [hl,hu) and species included
# (ispcut: 0 = all, >0 = one species, <0 = SPGROUP group, via `sp_groups`). A species flagged
# by SPLEAVE (leave_species) is never eligible — it's left out of both stocking and removal.
@inline function _cut_eligible(s::StandState, i, ispcut, dl, du, hl=0f0, hu=999f0,
                               sp_groups::Vector{Vector{Int32}}=_NOGROUPS)
    t = s.trees
    sp = t.species[i]
    @inbounds s.control.leave_species[sp] && return false   # SPLEAVE: leave this species
    d = t.dbh[i]
    (d < dl || d >= du) && return false
    h = t.height[i]
    (h < hl || h >= hu) && return false
    return _sp_in_cut(sp_groups, ispcut, sp)
end

# THINDBH/THINHT (cuts.f label_325→355→550→1100): thin the species class to a residual
# TPA (ctpa) and/or basal area (cba). The class is bounded by DBH (THINDBH, ICFLAG 8) or
# by HEIGHT (THINHT, ICFLAG 12) — VALMIN/VALMAX hold DBH or HT accordingly (cuts.f:637).
# cuteff = remove/cstock is applied to each eligible record's TPA in record order until
# the removal budget is spent (the last record is partial); the residual replaces PROB.
function _thindbh!(s::StandState, act::ScheduledActivity)
    t = s.trees; n = t.n
    n == 0 && return _NO_REMOVAL
    p1, p2, _eff, spf, ctpa_in, cba_in = act.params
    ispcut = floor(Int32, spf)
    if act.icflag == Int32(12)                         # THINHT: class on height
        dlo = 0f0; dhi = 999f0; hlo = p1; hhi = p2 < p1 ? 9999f0 : p2
    else                                               # THINDBH: class on DBH
        dlo = p1; dhi = p2 < p1 ? 9999f0 : p2; hlo = 0f0; hhi = 999f0
    end
    ctpa = max(0f0, ctpa_in); cba = max(0f0, cba_in)
    lbarea = cba > 0f0
    jtyp = lbarea ? 2 : 1

    grps = s.control.sp_groups                         # SPGROUP table (for ispcut<0)
    wk4 = Float32[t.tpa[i] for i in 1:n]               # pre-thin PROB copy
    cstock = _clsstk(s, wk4, n, jtyp, ispcut, dlo, dhi, hlo, hhi, grps)
    cstock <= 0f0 && return _NO_REMOVAL
    remove = cstock - (ctpa + cba)
    remove <= 0f0 && return _NO_REMOVAL
    cuteff = min(1f0, remove / cstock)

    # removed totals (per-acre, gross): Σ prem·{1, CFV, MCFV, SCFV, BFV}
    rtpa = 0f0; rcuft = 0f0; rmcuft = 0f0; rscuft = 0f0; rbdft = 0f0
    totcut = 0f0
    @inbounds for i in 1:n
        _cut_eligible(s, i, ispcut, dlo, dhi, hlo, hhi, grps) || continue
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
        _log_cut!(s, t, i, prem)             # ESTUMP
        rtpa += prem; rcuft += prem * t.cuft_vol[i]
        rmcuft += prem * t.merch_cuft_vol[i]; rscuft += prem * t.saw_cuft_vol[i]
        rbdft += prem * t.bdft_vol[i]
    end
    @inbounds for i in 1:n
        t.tpa[i] = wk4[i]                               # residual replaces PROB
    end
    return (tpa = rtpa, cuft = rcuft, mcuft = rmcuft, scuft = rscuft, bdft = rbdft)
end

# THINPRSC (cuts.f label_300 → label_550 lspecl path, ICFLAG 7): "prescription" thin.
# Despite the name it is not a per-DBH-class table — it removes the records the user
# PRE-MARKED with a cut code in the tree data (KUTKOD≥2; cuts.f:616-617 for the single-
# prescription case nps==1), at cutting efficiency `cuteff` (params[1]). All species,
# all DBH (label_300 leaves valmin=0/valmax=9999). Each marked record loses cuteff·tpa.
# (Multiple same-year prescriptions, nps>1, use the kutnow/KUTKOD-match path — deferred;
# snt01 uses one THINPRSC so nps==1.)
function _thinprsc!(s::StandState, act::ScheduledActivity)
    t = s.trees; n = t.n
    n == 0 && return _NO_REMOVAL
    cuteff = act.params[1] > 0f0 ? act.params[1] : s.control.cut_eff   # blank ⇒ EFF (CUTEFF default)
    cuteff <= 0f0 && return _NO_REMOVAL
    cuteff > 1f0 && (cuteff = 1f0)
    rtpa = 0f0; rcuft = 0f0; rmcuft = 0f0; rscuft = 0f0; rbdft = 0f0
    @inbounds for i in 1:n
        t.cut_code[i] < Int32(2) && continue        # only pre-marked (KUTKOD≥2) records
        prem = t.tpa[i] * cuteff
        prem <= 0f0 && continue
        t.tpa[i] -= prem
        _log_cut!(s, t, i, prem)             # ESTUMP
        rtpa  += prem;                  rcuft  += prem * t.cuft_vol[i]
        rmcuft += prem * t.merch_cuft_vol[i]; rscuft += prem * t.saw_cuft_vol[i]
        rbdft += prem * t.bdft_vol[i]
    end
    return (tpa = rtpa, cuft = rcuft, mcuft = rmcuft, scuft = rscuft, bdft = rbdft)
end

# Static (non-size) part of the RDPSRT cut-priority weight (cuts.f:1074):
#   IORDER[sp] (SPECPREF) + TCWT·IMC (TCONDMLT condition weight) + SPCLWT·ISPECL (special status).
# IMC is the tree's mortality/condition code clamped to 1..3 for live trees (intree.f:621). A higher
# weight ⇒ removed earlier. Defaults (TCWT=SPCLWT=0) leave just the species preference, so the
# common path is unchanged.
@inline function _cut_pref_wt(s::StandState, i::Integer)
    c = s.control; t = s.trees
    imc = Int(t.mort_code[i]); imc = imc <= 0 ? 1 : imc > 3 ? 3 : imc
    return Float32(c.cut_pref[t.species[i]]) + c.total_wt * Float32(imc) +
           c.special_wt * Float32(t.special[i])
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
    # A blank upper bound (parsed as 0) means "no upper limit" — guard with ≤ so the
    # common sparse form (e.g. `THINBTA <yr> <tpa>`, blank dbhhi/hthi) selects all
    # trees, not an empty [0,0) class. (Only an explicit dbhhi>dbhlo restricts.)
    valmax = dbhhi_p <= dbhlo ? 9999f0 : dbhhi_p
    hthi   = hthi_p  <= htlo  ? 9999f0 : hthi_p

    wk4 = Float32[t.tpa[i] for i in 1:n]
    cstock = _clsstk(s, wk4, n, jtyp, 0, dbhlo, valmax, htlo, hthi)
    cstock <= 0f0 && return _NO_REMOVAL
    remove = cstock - target
    remove <= 0f0 && return _NO_REMOVAL

    # priority key = ±DBH (size) + IORDER[sp] species preference (SPECPREF). cuts.f:635
    # computes WK2 = xsz + IORDER + weights for ALL records (the condition/density
    # weights default 0 until TCONDMLT/point keywords land) and ranks them with RDPSRT;
    # eligibility is applied in the removal loop, NOT by excluding from the sort — this
    # preserves the exact RDPSRT tie-break among equal keys (critical: see _rdpsrt!).
    elig = falses(n); key = Vector{Float32}(undef, n)
    @inbounds for i in 1:n
        elig[i] = _cut_eligible(s, i, 0, dbhlo, valmax, htlo, hthi)
        key[i] = (lbelow ? -t.dbh[i] : t.dbh[i]) + _cut_pref_wt(s, i)
    end
    order = Vector{Int32}(undef, n)
    _rdpsrt!(key, order)                             # descending, FVS tie-break

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
        _log_cut!(s, t, it, prem)            # ESTUMP (cuts.f:1713), in removal order
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

# THINSDI (cuts.f label_400, ICFLAG 10): thin to a residual Stand Density Index.
# Ported from FVSjulia (Oracle A) cuts.jl label_400 + sdical.jl SDICLS — NOT re-derived.
# SN sets LZEIDE=.TRUE. (sn/grinit.f:129), so SDICLS returns the ZEIDE *summation* SDI
#   SDI(class) = Σ_class tpa·(D/10)^1.605   (D ≥ DBHZEIDE, =0 for SN)
# — NOT the Stage linearised form (which equals Reineke N·(QMD/10)^1.605). For a stand
# with many small tripled stems the two differ materially (here 261 vs 293), and SDICLS
# uses Zeide. Removal (cuts.f:888): CUTEFF = REMOVE/SDIC (capped 1), REMOVE = SDI−target.
# Sparse form (icut=0) = LSPECL "throughout": apply CUTEFF proportionally to every
# eligible tree (residual = SDI·(1−CUTEFF) = target exactly, since the Zeide weight is
# fixed per tree). icut=1/2 = from below/above: RDPSRT by ∓DBH and remove whole records
# (cut_v = prem·(D/10)^1.605) until the SDI budget is met, last record partial.
@inline function _sdi_zeide(s::StandState, wk4, n, ispcut, dlo, dhi,
                            sp_groups::Vector{Vector{Int32}}=_NOGROUPS)
    t = s.trees; acc = 0f0
    @inbounds for i in 1:n
        d = t.dbh[i]; d <= 0f0 && continue          # DBHZEIDE = 0 (SN)
        _cut_eligible(s, i, ispcut, dlo, dhi, 0f0, 999f0, sp_groups) || continue
        acc += wk4[i] * (d / 10f0)^1.605f0
    end
    return acc
end

function _thin_sdi!(s::StandState, act::ScheduledActivity)
    t = s.trees; n = t.n
    n == 0 && return _NO_REMOVAL
    target, cuteff_p, ispcut_f, dbhlo, dbhhi_p, icut_f = act.params
    target = max(0f0, target)
    ispcut = Int(round(ispcut_f))
    icut   = Int(round(icut_f))
    valmin = dbhlo
    valmax = dbhhi_p <= dbhlo ? 999f0 : dbhhi_p

    grps = s.control.sp_groups                       # SPGROUP table (for ispcut<0)
    wk4  = Float32[t.tpa[i] for i in 1:n]
    sdic = _sdi_zeide(s, wk4, n, ispcut, valmin, valmax, grps)
    sdic <= target && return _NO_REMOVAL
    remove = sdic - target
    cuteff = remove / sdic
    cuteff > 1f0 && (cuteff = 1f0)

    rtpa = rcuft = rmcuft = rscuft = rbdft = 0f0
    @inline function _remove!(i, prem)
        prem <= 0f0 && return
        wk4[i] -= prem
        _log_cut!(s, t, i, prem)             # ESTUMP
        rtpa  += prem;                          rcuft  += prem * t.cuft_vol[i]
        rmcuft += prem * t.merch_cuft_vol[i];   rscuft += prem * t.saw_cuft_vol[i]
        rbdft += prem * t.bdft_vol[i]
        return
    end

    if icut == 0
        # LSPECL "throughout": proportional removal of CUTEFF from every eligible tree.
        @inbounds for i in 1:n
            d = t.dbh[i]; d <= 0f0 && continue
            _cut_eligible(s, i, ispcut, valmin, valmax, 0f0, 999f0, grps) || continue
            prem = wk4[i] * cuteff
            prem > wk4[i] && (prem = wk4[i])
            _remove!(i, prem)
        end
    else
        # from below (icut==1) / above: rank by ∓DBH, remove whole records until the SDI
        # budget is spent; the last record is partial (prem = remaining / (D/10)^1.605).
        lbelow = icut == 1
        ce = cuteff_p > 0f0 ? cuteff_p : 1f0
        key = Vector{Float32}(undef, n)
        @inbounds for i in 1:n
            key[i] = (lbelow ? -t.dbh[i] : t.dbh[i]) + _cut_pref_wt(s, i)
        end
        order = Vector{Int32}(undef, n); _rdpsrt!(key, order)
        totcut = 0f0
        @inbounds for it in order
            d = t.dbh[it]; d <= 0f0 && continue
            _cut_eligible(s, it, ispcut, valmin, valmax, 0f0, 999f0, grps) || continue
            totcut >= remove && break
            w  = (d / 10f0)^1.605f0
            prem = wk4[it] * ce
            prem > wk4[it] && (prem = wk4[it])
            cut_v = prem * w
            if totcut + cut_v > remove                      # last (partial) record
                prem = (remove - totcut) / w
                prem > wk4[it] && (prem = wk4[it])
                cut_v = prem * w
            end
            totcut += cut_v
            _remove!(it, prem)
        end
    end

    @inbounds for i in 1:n
        t.tpa[i] = wk4[i]
    end
    return (tpa = rtpa, cuft = rcuft, mcuft = rmcuft, scuft = rscuft, bdft = rbdft)
end

# THINRDEN (cuts.f label_400, ICFLAG 14): thin to a residual Curtis relative density.
# Ported from Oracle A cuts.jl label_400 + sdical.jl RDCLS. Curtis RD linearises about
# the stand QMD (q = Σ D²·tpa / Σ tpa over DBH≥DBHSTAGE=0):
#   tpafac = (0.25π/24²)·q^0.75 ;  diamfac = (0.75π/24²)·q^(−0.25)
#   RD(class) = Σ_class tpa·(tpafac + diamfac·D²)
# REMOVE = RD − target; CUTEFF = REMOVE/RD; sparse (icut=0) = proportional throughout;
# icut=1/2 = from below/above (sorted), per-tree weight = tpafac + diamfac·D².
@inline function _rd_curtis(s::StandState, wk4, n, ispcut, dlo, dhi,
                            sp_groups::Vector{Vector{Int32}}=_NOGROUPS)
    t = s.trees; clsd2 = 0f0; clstpa = 0f0
    @inbounds for i in 1:n
        d = t.dbh[i]; d <= 0f0 && continue
        clsd2 += d * d * wk4[i]; clstpa += wk4[i]
    end
    clstpa <= 0f0 && return (0f0, 0f0, 0f0)
    q = clsd2 / clstpa
    c2 = 24f0^2
    tpafac  = (0.25f0 * 3.14159f0 / c2) * q^0.75f0
    diamfac = (0.75f0 * 3.14159f0 / c2) * q^(0.75f0 - 1f0)
    crd = 0f0
    @inbounds for i in 1:n
        d = t.dbh[i]; d <= 0f0 && continue
        _cut_eligible(s, i, ispcut, dlo, dhi, 0f0, 999f0, sp_groups) || continue
        crd += (tpafac + diamfac * d * d) * wk4[i]
    end
    return (crd, tpafac, diamfac)
end

function _thin_rden!(s::StandState, act::ScheduledActivity)
    t = s.trees; n = t.n
    n == 0 && return _NO_REMOVAL
    target, cuteff_p, ispcut_f, dbhlo, dbhhi_p, icut_f = act.params
    target = max(0f0, target)
    ispcut = Int(round(ispcut_f))
    icut   = Int(round(icut_f))
    valmin = dbhlo
    valmax = dbhhi_p <= dbhlo ? 999f0 : dbhhi_p

    grps = s.control.sp_groups                       # SPGROUP table (for ispcut<0)
    wk4 = Float32[t.tpa[i] for i in 1:n]
    crd, tpafac, diamfac = _rd_curtis(s, wk4, n, ispcut, valmin, valmax, grps)
    crd <= target && return _NO_REMOVAL
    remove = crd - target
    cuteff = remove / crd
    cuteff > 1f0 && (cuteff = 1f0)

    rtpa = rcuft = rmcuft = rscuft = rbdft = 0f0
    @inline function _rm!(i, prem)
        prem <= 0f0 && return
        wk4[i] -= prem
        _log_cut!(s, t, i, prem)             # ESTUMP
        rtpa  += prem;                          rcuft  += prem * t.cuft_vol[i]
        rmcuft += prem * t.merch_cuft_vol[i];   rscuft += prem * t.saw_cuft_vol[i]
        rbdft += prem * t.bdft_vol[i]
        return
    end

    if icut == 0
        @inbounds for i in 1:n
            d = t.dbh[i]; d <= 0f0 && continue
            _cut_eligible(s, i, ispcut, valmin, valmax, 0f0, 999f0, grps) || continue
            prem = wk4[i] * cuteff
            prem > wk4[i] && (prem = wk4[i])
            _rm!(i, prem)
        end
    else
        lbelow = icut == 1
        ce = cuteff_p > 0f0 ? cuteff_p : 1f0
        key = Vector{Float32}(undef, n)
        @inbounds for i in 1:n
            key[i] = (lbelow ? -t.dbh[i] : t.dbh[i]) + _cut_pref_wt(s, i)
        end
        order = Vector{Int32}(undef, n); _rdpsrt!(key, order)
        totcut = 0f0
        @inbounds for it in order
            d = t.dbh[it]; d <= 0f0 && continue
            _cut_eligible(s, it, ispcut, valmin, valmax, 0f0, 999f0, grps) || continue
            totcut >= remove && break
            w = tpafac + diamfac * d * d
            w <= 0f0 && continue
            prem = wk4[it] * ce
            prem > wk4[it] && (prem = wk4[it])
            cut_v = prem * w
            if totcut + cut_v > remove
                prem = (remove - totcut) / w
                prem > wk4[it] && (prem = wk4[it])
                cut_v = prem * w
            end
            totcut += cut_v
            _rm!(it, prem)
        end
    end

    @inbounds for i in 1:n
        t.tpa[i] = wk4[i]
    end
    return (tpa = rtpa, cuft = rcuft, mcuft = rmcuft, scuft = rscuft, bdft = rbdft)
end

# THINAUTO (cuts.f label_150, ICFLAG 1): automatic thinning to a normal-stocking target.
# Ported from Oracle A cutstk.jl AUTSTK + cuts.f:582-597. AUTSTK gives FULSTK = normal
# full stocking (stems/acre) from the BA-weighted mean species SDImax (SDIDEF) and the
# stand QMD (RMSQD): FULSTK = 1 / (0.02483133/tmpmax · QMD^1.605). The thin fires only
# when STOCK ≥ AUTMAX%·FULSTK (SN default 60), and removes FROM BELOW down to a residual
# of AUTMIN%·FULSTK (SN default 45) TPA — i.e. exactly a THINBTA-below to that residual,
# so it delegates to the validated _thin_sorted! (icflag 3).
@inline function _autstk(t::TreeList, wk4, n, sp_sdi_def)
    totba = 0f0; temba = 0f0; sdsq = 0f0; stpa = 0f0
    @inbounds for i in 1:n
        d = t.dbh[i]
        tba = _BA_PER_TREE * d * d * wk4[i]
        temba += tba * sp_sdi_def[t.species[i]]
        totba += tba
        sdsq += d * d * wk4[i]; stpa += wk4[i]
    end
    (totba <= 1f0 || temba <= 1f0) && return 0f0
    tmpmax = temba / totba
    rmsqd  = stpa > 0f0 ? sqrt(sdsq / stpa) : 0f0      # current stand QMD (RMSQD)
    q = rmsqd > 2f0 ? rmsqd : 2f0
    return 1f0 / (0.02483133f0 / tmpmax * q^1.605f0)
end

function _thin_auto!(s::StandState, act::ScheduledActivity)
    t = s.trees; n = t.n
    n == 0 && return _NO_REMOVAL
    p = act.params
    autmin = p[1] > 0f0 ? p[1] : 45f0          # SN grinit defaults when fields blank
    autmax = p[2] > 0f0 ? p[2] : 60f0
    eff    = p[3] > 0f0 ? p[3] : s.control.cut_eff   # blank ⇒ EFF (CUTEFF default)
    wk4 = Float32[t.tpa[i] for i in 1:n]
    fulstk = _autstk(t, wk4, n, s.plot.sp_sdi_def)
    fulstk <= 0f0 && return _NO_REMOVAL
    stock = 0f0
    @inbounds for i in 1:n; stock += wk4[i]; end
    stock < (autmax / 100f0) * fulstk && return _NO_REMOVAL    # cuts.f:588 gate
    rstock = (autmin / 100f0) * fulstk
    # from-below thin to rstock TPA — reuse the validated THINBTA-below path (icflag 3).
    return _thin_sorted!(s, ScheduledActivity(act.year, Int32(3),
                                              (rstock, eff, 0f0, 0f0, 0f0, 0f0)))
end

# THINCC (cuts.f label_400, ICFLAG 11): thin to a residual canopy cover %.
# Ported from Oracle A cuts.jl label_400 + sdical.jl CCCLS + cwidth.jl CWIDTH. The class
# metric is crown cover AREA: CCCLS = Σ_class tpa·CW², where CW is each tree's FOREST-grown
# crown width (CWIDTH/CWCALC with the tree's actual crown ratio — not the open-grown
# iwho=1 form CCF uses). The keyword target is a cover PERCENT → equivalent crown area
# (cuts.f:816, CCC=1): CSDI = −(43560·ln(1−CC/100)/0.785398). REMOVE = area − CSDI;
# CUTEFF = REMOVE/area; sparse (icut=0) = proportional throughout; icut=1/2 sorted.
function _thin_cc!(s::StandState, act::ScheduledActivity)
    t = s.trees; n = t.n
    n == 0 && return _NO_REMOVAL
    cc_target, cuteff_p, ispcut_f, dbhlo, dbhhi_p, icut_f = act.params
    cc_target >= 100f0 && return _NO_REMOVAL           # cuts.f:823 cancel if ≥100%
    ispcut = Int(round(ispcut_f))
    icut   = Int(round(icut_f))
    valmin = dbhlo
    valmax = dbhhi_p <= dbhlo ? 999f0 : dbhhi_p
    csdi = cc_target <= 0f0 ? 0f0 :
           -(43560f0 * log(1f0 - cc_target / 100f0) / 0.785398f0)

    # per-tree forest-grown crown width (CRWDTH array, cwidth.f)
    p = s.plot
    cw = Vector{Float32}(undef, n)
    @inbounds for i in 1:n
        sp2 = s.species.class_codes[t.species[i], 1][1:2]
        cw[i] = crown_width(s.coef, sp2, t.dbh[i], t.height[i],
                            Float32(t.crown_pct[i]), 0, p.latitude, p.longitude, p.elevation)
    end

    grps = s.control.sp_groups                         # SPGROUP table (for ispcut<0)
    wk4 = Float32[t.tpa[i] for i in 1:n]
    area = 0f0
    @inbounds for i in 1:n
        d = t.dbh[i]; d <= 0f0 && continue
        _cut_eligible(s, i, ispcut, valmin, valmax, 0f0, 999f0, grps) || continue
        area += cw[i] * cw[i] * wk4[i]
    end
    area <= csdi && return _NO_REMOVAL
    remove = area - csdi
    cuteff = remove / area
    cuteff > 1f0 && (cuteff = 1f0)

    rtpa = rcuft = rmcuft = rscuft = rbdft = 0f0
    @inline function _rmc!(i, prem)
        prem <= 0f0 && return
        wk4[i] -= prem
        _log_cut!(s, t, i, prem)             # ESTUMP
        rtpa  += prem;                          rcuft  += prem * t.cuft_vol[i]
        rmcuft += prem * t.merch_cuft_vol[i];   rscuft += prem * t.saw_cuft_vol[i]
        rbdft += prem * t.bdft_vol[i]
        return
    end

    if icut == 0
        @inbounds for i in 1:n
            d = t.dbh[i]; d <= 0f0 && continue
            _cut_eligible(s, i, ispcut, valmin, valmax, 0f0, 999f0, grps) || continue
            prem = wk4[i] * cuteff
            prem > wk4[i] && (prem = wk4[i])
            _rmc!(i, prem)
        end
    else
        lbelow = icut == 1
        ce = cuteff_p > 0f0 ? cuteff_p : 1f0
        key = Vector{Float32}(undef, n)
        @inbounds for i in 1:n
            key[i] = (lbelow ? -t.dbh[i] : t.dbh[i]) + _cut_pref_wt(s, i)
        end
        order = Vector{Int32}(undef, n); _rdpsrt!(key, order)
        totcut = 0f0
        @inbounds for it in order
            d = t.dbh[it]; d <= 0f0 && continue
            _cut_eligible(s, it, ispcut, valmin, valmax, 0f0, 999f0, grps) || continue
            totcut >= remove && break
            w = cw[it] * cw[it]
            w <= 0f0 && continue
            prem = wk4[it] * ce
            prem > wk4[it] && (prem = wk4[it])
            cut_v = prem * w
            if totcut + cut_v > remove
                prem = (remove - totcut) / w
                prem > wk4[it] && (prem = wk4[it])
                cut_v = prem * w
            end
            totcut += cut_v
            _rmc!(it, prem)
        end
    end

    @inbounds for i in 1:n
        t.tpa[i] = wk4[i]
    end
    return (tpa = rtpa, cuft = rcuft, mcuft = rmcuft, scuft = rscuft, bdft = rbdft)
end

# THINQFA (cuts.f label_350, ICFLAG 17): Q-factor diameter-distribution thin. Ported from
# Oracle A cutqfa.jl (CUTQFA + CYCQFA) + cuts.f:671-759. Divides [valmin,valmax] into
# diameter classes of width DIACW; builds a target negative-exponential (Q-factor)
# distribution that meets the stand-level target (TARQFA in TPA/BA/SDI per QFATAR),
# yielding a residual target per class (clstar2); then thins each class that has excess
# down to its residual via the validated per-class residual path (_thindbh! for TPA/BA,
# _thin_sdi! for SDI). aux carries QFATAR (0=BA, 1=TPA, 2=SDI).
function _thin_qfa!(s::StandState, act::ScheduledActivity)
    t = s.trees; n = t.n
    n == 0 && return _NO_REMOVAL
    valmin, valmax_p, spec_f, qfac, diacw, tarqfa_p = act.params
    valmax = valmax_p < valmin ? 9999f0 : valmax_p
    ispcut = Int(round(spec_f))
    qfatar = Int(round(act.aux))                      # 0=BA, 1=TPA, 2=SDI
    ctar0  = max(0f0, tarqfa_p)
    (qfac <= 0f0 || diacw <= 0f0) && return _NO_REMOVAL
    ndcls = Int(floor((valmax - valmin) / diacw))
    (ndcls < 1 || ndcls > 30) && return _NO_REMOVAL

    # diameter-class midpoints (top class centred at valmax−DIACW/2, descending)
    dcls = Vector{Float32}(undef, ndcls)
    dcls[ndcls] = valmax - diacw / 2f0
    for i in ndcls-1:-1:1; dcls[i] = dcls[i+1] - diacw; end

    # per-class current TPA (tpacls1) and per-tree metric (clstar1: BA/TPA/SDI per tree)
    grps = s.control.sp_groups                          # SPGROUP table (for ispcut<0)
    wk4 = Float32[t.tpa[i] for i in 1:n]
    tpacls1 = zeros(Float32, ndcls); clstar1 = zeros(Float32, ndcls)
    for i in 1:ndcls
        dlo = dcls[i] - diacw / 2f0; dhi = dcls[i] + diacw / 2f0
        ctpa = _clsstk(s, wk4, n, 1, ispcut, dlo, dhi, 0f0, 999f0, grps)
        tpacls1[i] = ctpa
        clstar1[i] = if ctpa <= 0f0
            0f0
        elseif qfatar <= 0      # BA per tree
            _clsstk(s, wk4, n, 2, ispcut, dlo, dhi, 0f0, 999f0, grps) / ctpa
        elseif qfatar <= 1      # TPA per tree = 1
            1f0
        else                    # Zeide SDI per tree
            _sdi_zeide(s, wk4, n, ispcut, dlo, dhi, grps) / ctpa
        end
    end
    suminv = 0f0
    for i in 1:ndcls; suminv += clstar1[i] * tpacls1[i]; end
    suminv < ctar0 && return _NO_REMOVAL              # not enough inventory to meet target

    # Q-factor convergence (cutqfa.jl:111-164) → clstar2 = residual target per class
    tpacls4 = copy(tpacls1)
    tpacls2 = zeros(Float32, ndcls); tpacls3 = zeros(Float32, ndcls)
    clstar2 = zeros(Float32, ndcls)
    dinom = 0f0
    for i in 1:ndcls; dinom += clstar1[i] / qfac^(i-1); end
    dinom <= 0f0 && return _NO_REMOVAL
    tpa1 = ctar0 / dinom
    ctar = ctar0
    for _iter in 1:200
        for i in 1:ndcls
            tpacls2[i] = tpa1 / qfac^(i-1)
            tpacls3[i] = tpacls4[i] - tpacls2[i]
        end
        sumtar = 0f0
        for i in 1:ndcls
            sumtar += (tpacls3[i] <= 0f0 ? tpacls4[i] : tpacls2[i]) * clstar1[i]
        end
        dinom2 = 0f0
        for i in 1:ndcls
            if tpacls3[i] > 0f0
                dinom2 += clstar1[i] / qfac^(i-1); tpacls4[i] = tpacls3[i]
            else
                tpacls4[i] = 0f0
            end
        end
        sumv = 0f0
        for i in 1:ndcls
            clstar2[i] = tpacls3[i] > 0f0 ? clstar2[i] + tpacls2[i] * clstar1[i] :
                                            tpacls1[i] * clstar1[i]
            sumv += clstar2[i]
        end
        abs(ctar0 - sumv) < 0.1f0 && break
        dinom2 > 0f0 && (tpa1 = (ctar - sumtar) / dinom2)
        ctar -= sumtar
        abs(tpa1) <= 0.1f0 && break
    end

    # thin each excess class down to its residual target (per-class residual path)
    rtpa = rcuft = rmcuft = rscuft = rbdft = 0f0
    for i in 1:ndcls
        (tpacls1[i] * clstar1[i] - clstar2[i]) > 1f-5 || continue
        dlo = dcls[i] - diacw / 2f0; dhi = dcls[i] + diacw / 2f0
        r = if qfatar <= 1
            ctpa_t = qfatar <= 0 ? 0f0 : clstar2[i]
            cba_t  = qfatar <= 0 ? clstar2[i] : 0f0
            _thindbh!(s, ScheduledActivity(act.year, Int32(8),
                                           (dlo, dhi, 0f0, Float32(ispcut), ctpa_t, cba_t)))
        else
            _thin_sdi!(s, ScheduledActivity(act.year, Int32(10),
                                            (clstar2[i], 0f0, Float32(ispcut), dlo, dhi, 1f0)))
        end
        rtpa += r.tpa; rcuft += r.cuft; rmcuft += r.mcuft; rscuft += r.scuft; rbdft += r.bdft
    end
    return (tpa = rtpa, cuft = rcuft, mcuft = rmcuft, scuft = rscuft, bdft = rbdft)
end

# THINPT (cuts.f label_475→400, ICFLAG 15): point-level thin. SETPTHIN sets the point
# (pt_point; 0 = all points / LPTALL) and the residual metric (ithnpa: 1=TPA, 2=BA, 3=SDI
# Zeide, 4=CC crown cover, 5=Curtis RD). THINPT gives the residual + DBH class + direction
# (0 throughout / 1 below / 2 above). Each targeted point is thinned independently to the
# residual in its metric; the per-point metric scales the point's per-acre stocking by
# (PI−NONSTK) (the jpnum path in CLSSTK/SDICLS/…). Ported from Oracle A cuts.jl label_475.
function _thin_pt!(s::StandState, act::ScheduledActivity, pt_point::Int32, ithnpa::Int32)
    t = s.trees; n = t.n
    n == 0 && return _NO_REMOVAL
    (ithnpa < 1 || ithnpa > 5) && return _NO_REMOVAL
    residual, _eff, spec_f, dbhlo, dbhhi_p, dir_f = act.params
    residual = max(0f0, residual)
    ispcut = Int(round(spec_f))
    dir    = Int(round(dir_f))
    valmin = dbhlo; valmax = dbhhi_p <= dbhlo ? 999f0 : dbhhi_p
    grps = s.control.sp_groups                          # SPGROUP table (for ispcut<0)
    scale  = Float32(s.plot.points_inv) - Float32(s.plot.nonstockable)
    scale <= 0f0 && (scale = 1f0)
    pmax = max(1, Int(s.plot.points_inv))
    points = pt_point == 0 ? (1:pmax) : (Int(pt_point):Int(pt_point))

    # per-tree weight w[i] for the metric (point-independent for TPA/BA/SDI/CC; RD is
    # point-dependent and filled inside the loop).
    w = zeros(Float32, n)
    @inbounds for i in 1:n
        d = t.dbh[i]; d <= 0f0 && continue
        if ithnpa == 1
            w[i] = 1f0
        elseif ithnpa == 2
            w[i] = d * d * _BA_PER_TREE
        elseif ithnpa == 3
            w[i] = (d / 10f0)^1.605f0
        elseif ithnpa == 4
            sp2 = s.species.class_codes[t.species[i], 1][1:2]
            cw = crown_width(s.coef, sp2, d, t.height[i], Float32(t.crown_pct[i]), 0,
                             s.plot.latitude, s.plot.longitude, s.plot.elevation)
            w[i] = cw * cw
        end
    end

    wk4 = Float32[t.tpa[i] for i in 1:n]
    rtpa = rcuft = rmcuft = rscuft = rbdft = 0f0
    @inline function _rmp!(i, prem)
        prem <= 0f0 && return
        wk4[i] -= prem
        _log_cut!(s, t, i, prem)             # ESTUMP
        rtpa  += prem;                          rcuft  += prem * t.cuft_vol[i]
        rmcuft += prem * t.merch_cuft_vol[i];   rscuft += prem * t.saw_cuft_vol[i]
        rbdft += prem * t.bdft_vol[i]
        return
    end

    for jp in points
        # Curtis RD weight is point-specific (linearised about the point's QMD)
        if ithnpa == 5
            cd2 = 0f0; ct = 0f0
            @inbounds for i in 1:n
                (t.plot_id[i] == jp && t.dbh[i] > 0f0) || continue
                cd2 += t.dbh[i]^2 * wk4[i]; ct += wk4[i]
            end
            ct <= 0f0 && continue
            q = cd2 / ct; c2 = 24f0^2
            a_rd = (0.25f0 * 3.14159f0 / c2) * q^0.75f0
            b_rd = (0.75f0 * 3.14159f0 / c2) * q^(0.75f0 - 1f0)
            @inbounds for i in 1:n
                t.plot_id[i] == jp && t.dbh[i] > 0f0 && (w[i] = a_rd + b_rd * t.dbh[i]^2)
            end
        end
        # point metric (scaled to the point's per-acre stocking)
        metric = 0f0
        @inbounds for i in 1:n
            t.dbh[i] <= 0f0 && continue
            (t.plot_id[i] == jp) || continue
            _cut_eligible(s, i, ispcut, valmin, valmax, 0f0, 999f0, grps) || continue
            metric += w[i] * wk4[i]
        end
        metric *= scale
        metric <= residual && continue
        cuteff = (metric - residual) / metric

        if dir == 0
            @inbounds for i in 1:n
                t.dbh[i] <= 0f0 && continue
                (t.plot_id[i] == jp) || continue
                _cut_eligible(s, i, ispcut, valmin, valmax, 0f0, 999f0, grps) || continue
                prem = wk4[i] * cuteff
                prem > wk4[i] && (prem = wk4[i])
                _rmp!(i, prem)
            end
        else                                            # from below (1) / above (2)
            lbelow = dir == 1
            remove = metric - residual                  # in metric units (already ×scale)
            idx = Int32[i for i in 1:n if t.plot_id[i] == jp && t.dbh[i] > 0f0 &&
                        _cut_eligible(s, i, ispcut, valmin, valmax, 0f0, 999f0, grps)]
            isempty(idx) && continue
            key = Float32[(lbelow ? -t.dbh[i] : t.dbh[i]) + _cut_pref_wt(s, i) for i in idx]
            ord = Vector{Int32}(undef, length(idx)); _rdpsrt!(key, ord)
            totcut = 0f0
            @inbounds for k in ord
                i = idx[k]; totcut >= remove && break
                cv = w[i] * scale
                cv <= 0f0 && continue
                prem = wk4[i]
                if totcut + prem * cv > remove
                    prem = (remove - totcut) / cv
                    prem > wk4[i] && (prem = wk4[i])
                end
                totcut += prem * cv
                _rmp!(i, prem)
            end
        end
    end

    @inbounds for i in 1:n
        t.tpa[i] = wk4[i]
    end
    return (tpa = rtpa, cuft = rcuft, mcuft = rmcuft, scuft = rscuft, bdft = rbdft)
end
