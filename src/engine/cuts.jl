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
function cuts!(s::StandState; fint::Float32 = 5f0)
    sched = s.control.schedule
    conds = s.control.conditionals
    (isempty(sched) && isempty(conds)) && return _NO_REMOVAL
    # Year of the current cycle (matches summary_row): inventory year + cycle·period.
    # (cycle_year only stores the inventory year; later years are derived.)
    yr = Int32(Int(s.control.cycle_year[1]) + Int(s.control.cycle) * round(Int, fint))
    yr in s.control.years_cut && return _NO_REMOVAL   # idempotent: already cut this year
    # Effective activities this cycle: the dated ones (year==yr) plus any IF/THEN block
    # whose algebraic condition is true this cycle (EVMON), with their year set to yr.
    # THINAUTO (icflag 1) is a RECURRING auto-thin: once scheduled it re-evaluates the
    # AUTMAX stocking gate every cycle from its start year (cuts.f auto path), not once.
    acts = ScheduledActivity[a for a in sched
                             if a.year == yr || (a.icflag == Int32(1) && yr >= a.year)]
    if !isempty(conds)
        ctx = EventCtx(Int(s.control.cycle) + 1, Int(yr), s)   # FVS CYCLE is 1-based
        for c in conds
            eval_event(c.cond, ctx) != 0f0 || continue
            for a in c.acts
                push!(acts, ScheduledActivity(yr, a.icflag, a.params))
            end
        end
    end
    isempty(acts) && return _NO_REMOVAL
    # PASS 1 — cut MODIFIERS for this year (set state the methods read), before any
    # method runs (cuts.f processes SPECPREF/SPLEAVE/… then the thin in the cycle).
    @inbounds for act in acts
        act.icflag == Int32(201) && _apply_specpref!(s, act)
    end
    # PASS 2 — cut METHODS.
    rem = _NO_REMOVAL
    applied = false
    @inbounds for act in acts
        ic = act.icflag
        # only CUTS methods here; establishment (427/430/431) + the SPECPREF modifier
        # (201, applied above) are consumed elsewhere (ESNUTR).
        ic in (Int32(3), Int32(4), Int32(5), Int32(6), Int32(7), Int32(8), Int32(10), Int32(12), Int32(14), Int32(1), Int32(11), Int32(17)) || continue
        applied = true
        r = (ic == Int32(8) || ic == Int32(12)) ? _thindbh!(s, act) : # DBH-class / HT-class residual
            ic == Int32(7)  ? _thinprsc!(s, act) :                     # prescription (cut-code marked)
            ic == Int32(10) ? _thin_sdi!(s, act) :                     # THINSDI (Zeide target SDI)
            ic == Int32(14) ? _thin_rden!(s, act) :                    # THINRDEN (Curtis RD)
            ic == Int32(1)  ? _thin_auto!(s, act) :                    # THINAUTO (auto to normal stocking)
            ic == Int32(11) ? _thin_cc!(s, act) :                      # THINCC (crown cover)
            ic == Int32(17) ? _thin_qfa!(s, act) :                     # THINQFA (Q-factor distribution)
            (ic in (Int32(3), Int32(4), Int32(5), Int32(6))) ? _thin_sorted!(s, act) :  # BTA/ATA/BBA/ABA
            _NO_REMOVAL
        rem = (tpa = rem.tpa + r.tpa, cuft = rem.cuft + r.cuft,
               mcuft = rem.mcuft + r.mcuft, scuft = rem.scuft + r.scuft,
               bdft = rem.bdft + r.bdft)
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
    end
    # isp < 0 (species-group) needs SPGROUP/ISPGRP — handled when that lands.
    return
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

    wk4 = Float32[t.tpa[i] for i in 1:n]               # pre-thin PROB copy
    cstock = _clsstk(t, wk4, n, jtyp, ispcut, dlo, dhi, hlo, hhi)
    cstock <= 0f0 && return _NO_REMOVAL
    remove = cstock - (ctpa + cba)
    remove <= 0f0 && return _NO_REMOVAL
    cuteff = min(1f0, remove / cstock)

    # removed totals (per-acre, gross): Σ prem·{1, CFV, MCFV, SCFV, BFV}
    rtpa = 0f0; rcuft = 0f0; rmcuft = 0f0; rscuft = 0f0; rbdft = 0f0
    totcut = 0f0
    @inbounds for i in 1:n
        _cut_eligible(t, i, ispcut, dlo, dhi, hlo, hhi) || continue
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
    cuteff = act.params[1]
    cuteff <= 0f0 && return _NO_REMOVAL
    cuteff > 1f0 && (cuteff = 1f0)
    rtpa = 0f0; rcuft = 0f0; rmcuft = 0f0; rscuft = 0f0; rbdft = 0f0
    @inbounds for i in 1:n
        t.cut_code[i] < Int32(2) && continue        # only pre-marked (KUTKOD≥2) records
        prem = t.tpa[i] * cuteff
        prem <= 0f0 && continue
        t.tpa[i] -= prem
        rtpa  += prem;                  rcuft  += prem * t.cuft_vol[i]
        rmcuft += prem * t.merch_cuft_vol[i]; rscuft += prem * t.saw_cuft_vol[i]
        rbdft += prem * t.bdft_vol[i]
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
    # A blank upper bound (parsed as 0) means "no upper limit" — guard with ≤ so the
    # common sparse form (e.g. `THINBTA <yr> <tpa>`, blank dbhhi/hthi) selects all
    # trees, not an empty [0,0) class. (Only an explicit dbhhi>dbhlo restricts.)
    valmax = dbhhi_p <= dbhlo ? 9999f0 : dbhhi_p
    hthi   = hthi_p  <= htlo  ? 9999f0 : hthi_p

    wk4 = Float32[t.tpa[i] for i in 1:n]
    cstock = _clsstk(t, wk4, n, jtyp, 0, dbhlo, valmax, htlo, hthi)
    cstock <= 0f0 && return _NO_REMOVAL
    remove = cstock - target
    remove <= 0f0 && return _NO_REMOVAL

    # priority key = ±DBH (size) + IORDER[sp] species preference (SPECPREF). cuts.f:635
    # computes WK2 = xsz + IORDER + weights for ALL records (the condition/density
    # weights default 0 until TCONDMLT/point keywords land) and ranks them with RDPSRT;
    # eligibility is applied in the removal loop, NOT by excluding from the sort — this
    # preserves the exact RDPSRT tie-break among equal keys (critical: see _rdpsrt!).
    pref = s.control.cut_pref
    elig = falses(n); key = Vector{Float32}(undef, n)
    @inbounds for i in 1:n
        elig[i] = _cut_eligible(t, i, 0, dbhlo, valmax, htlo, hthi)
        key[i] = (lbelow ? -t.dbh[i] : t.dbh[i]) + Float32(pref[t.species[i]])
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
@inline function _sdi_zeide(t::TreeList, wk4, n, ispcut, dlo, dhi)
    s = 0f0
    @inbounds for i in 1:n
        d = t.dbh[i]; d <= 0f0 && continue          # DBHZEIDE = 0 (SN)
        _cut_eligible(t, i, ispcut, dlo, dhi) || continue
        s += wk4[i] * (d / 10f0)^1.605f0
    end
    return s
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

    wk4  = Float32[t.tpa[i] for i in 1:n]
    sdic = _sdi_zeide(t, wk4, n, ispcut, valmin, valmax)
    sdic <= target && return _NO_REMOVAL
    remove = sdic - target
    cuteff = remove / sdic
    cuteff > 1f0 && (cuteff = 1f0)

    rtpa = rcuft = rmcuft = rscuft = rbdft = 0f0
    @inline function _remove!(i, prem)
        prem <= 0f0 && return
        wk4[i] -= prem
        rtpa  += prem;                          rcuft  += prem * t.cuft_vol[i]
        rmcuft += prem * t.merch_cuft_vol[i];   rscuft += prem * t.saw_cuft_vol[i]
        rbdft += prem * t.bdft_vol[i]
        return
    end

    if icut == 0
        # LSPECL "throughout": proportional removal of CUTEFF from every eligible tree.
        @inbounds for i in 1:n
            d = t.dbh[i]; d <= 0f0 && continue
            _cut_eligible(t, i, ispcut, valmin, valmax) || continue
            prem = wk4[i] * cuteff
            prem > wk4[i] && (prem = wk4[i])
            _remove!(i, prem)
        end
    else
        # from below (icut==1) / above: rank by ∓DBH, remove whole records until the SDI
        # budget is spent; the last record is partial (prem = remaining / (D/10)^1.605).
        lbelow = icut == 1
        ce = cuteff_p > 0f0 ? cuteff_p : 1f0
        pref = s.control.cut_pref
        key = Vector{Float32}(undef, n)
        @inbounds for i in 1:n
            key[i] = (lbelow ? -t.dbh[i] : t.dbh[i]) + Float32(pref[t.species[i]])
        end
        order = Vector{Int32}(undef, n); _rdpsrt!(key, order)
        totcut = 0f0
        @inbounds for it in order
            d = t.dbh[it]; d <= 0f0 && continue
            _cut_eligible(t, it, ispcut, valmin, valmax) || continue
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
@inline function _rd_curtis(t::TreeList, wk4, n, ispcut, dlo, dhi)
    clsd2 = 0f0; clstpa = 0f0
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
        _cut_eligible(t, i, ispcut, dlo, dhi) || continue
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

    wk4 = Float32[t.tpa[i] for i in 1:n]
    crd, tpafac, diamfac = _rd_curtis(t, wk4, n, ispcut, valmin, valmax)
    crd <= target && return _NO_REMOVAL
    remove = crd - target
    cuteff = remove / crd
    cuteff > 1f0 && (cuteff = 1f0)

    rtpa = rcuft = rmcuft = rscuft = rbdft = 0f0
    @inline function _rm!(i, prem)
        prem <= 0f0 && return
        wk4[i] -= prem
        rtpa  += prem;                          rcuft  += prem * t.cuft_vol[i]
        rmcuft += prem * t.merch_cuft_vol[i];   rscuft += prem * t.saw_cuft_vol[i]
        rbdft += prem * t.bdft_vol[i]
        return
    end

    if icut == 0
        @inbounds for i in 1:n
            d = t.dbh[i]; d <= 0f0 && continue
            _cut_eligible(t, i, ispcut, valmin, valmax) || continue
            prem = wk4[i] * cuteff
            prem > wk4[i] && (prem = wk4[i])
            _rm!(i, prem)
        end
    else
        lbelow = icut == 1
        ce = cuteff_p > 0f0 ? cuteff_p : 1f0
        pref = s.control.cut_pref
        key = Vector{Float32}(undef, n)
        @inbounds for i in 1:n
            key[i] = (lbelow ? -t.dbh[i] : t.dbh[i]) + Float32(pref[t.species[i]])
        end
        order = Vector{Int32}(undef, n); _rdpsrt!(key, order)
        totcut = 0f0
        @inbounds for it in order
            d = t.dbh[it]; d <= 0f0 && continue
            _cut_eligible(t, it, ispcut, valmin, valmax) || continue
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
    eff    = p[3] > 0f0 ? p[3] : 1f0
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

    wk4 = Float32[t.tpa[i] for i in 1:n]
    area = 0f0
    @inbounds for i in 1:n
        d = t.dbh[i]; d <= 0f0 && continue
        _cut_eligible(t, i, ispcut, valmin, valmax) || continue
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
        rtpa  += prem;                          rcuft  += prem * t.cuft_vol[i]
        rmcuft += prem * t.merch_cuft_vol[i];   rscuft += prem * t.saw_cuft_vol[i]
        rbdft += prem * t.bdft_vol[i]
        return
    end

    if icut == 0
        @inbounds for i in 1:n
            d = t.dbh[i]; d <= 0f0 && continue
            _cut_eligible(t, i, ispcut, valmin, valmax) || continue
            prem = wk4[i] * cuteff
            prem > wk4[i] && (prem = wk4[i])
            _rmc!(i, prem)
        end
    else
        lbelow = icut == 1
        ce = cuteff_p > 0f0 ? cuteff_p : 1f0
        pref = s.control.cut_pref
        key = Vector{Float32}(undef, n)
        @inbounds for i in 1:n
            key[i] = (lbelow ? -t.dbh[i] : t.dbh[i]) + Float32(pref[t.species[i]])
        end
        order = Vector{Int32}(undef, n); _rdpsrt!(key, order)
        totcut = 0f0
        @inbounds for it in order
            d = t.dbh[it]; d <= 0f0 && continue
            _cut_eligible(t, it, ispcut, valmin, valmax) || continue
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
    wk4 = Float32[t.tpa[i] for i in 1:n]
    tpacls1 = zeros(Float32, ndcls); clstar1 = zeros(Float32, ndcls)
    for i in 1:ndcls
        dlo = dcls[i] - diacw / 2f0; dhi = dcls[i] + diacw / 2f0
        ctpa = _clsstk(t, wk4, n, 1, ispcut, dlo, dhi)
        tpacls1[i] = ctpa
        clstar1[i] = if ctpa <= 0f0
            0f0
        elseif qfatar <= 0      # BA per tree
            _clsstk(t, wk4, n, 2, ispcut, dlo, dhi) / ctpa
        elseif qfatar <= 1      # TPA per tree = 1
            1f0
        else                    # Zeide SDI per tree
            _sdi_zeide(t, wk4, n, ispcut, dlo, dhi) / ctpa
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
