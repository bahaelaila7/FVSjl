# =============================================================================
# econ_calc.jl — faithful ECON per-cycle engine: ECIN keyword storage, ECSETP, ECSTATUS,
# ECHARV, ECCALC (valueHarvest / calcEcon / computeIRR / calcAppreAmt / calcAppreSev /
# sevHrvCosts / sevHrvRevenues / calcAnnCostRevSEV) and the DBSECSUM FVS_EconSummary rows.
#
# Ported line-for-line from bin/FVS{v}_buildDir/{ecin,ecsetp,ecstatus,echarv,eccalc,dbsecsum}.f.
# State = EconCalc (core/state.jl: /ECNCOM/ + /ECONSAVE/). Float32 throughout (Fortran REAL);
# REAL**INTEGER uses the libgcc __powisf2 square-and-multiply (`_ec_powi`), REAL**REAL glibc powf.
#
# FVS quirks reproduced deliberately (faithful, not "fixed"):
#  * valueHarvest (eccalc.f:207-332) values a harvest ONLY when PCTSPEC set pctMinUnits>0 — the
#    commercial-harvest branch is nested inside `if (pctMinUnits>0)`; with no PCTSPEC no harvest
#    cost/revenue is ever accrued (live FVS confirmed: Undiscounted_Cost = annual cost only).
#  * VAR_PCT costs in valueHarvest read the varHrv* arrays indexed by the PCTVRCST counter
#    (eccalc.f:227-241), while sevHrvCosts reads varPct* (with varHrvUnits(kw)).
#  * HRVVRCST writes its supplemental DURATIONS into fixHrvDur(fixHrvCnt,:) (ecin.f:203-205), so
#    varHrvDur stays 0 (no appreciation) and the last HRVFXCST's durations are overwritten/zeroed.
# =============================================================================

const ECON_TPA_U = 1; const ECON_BF_U = 2; const ECON_FT3_U = 3; const ECON_BFLOG_U = 4; const ECON_FT3LOG_U = 5
const ECON_PER_ACRE = 6; const ECON_TPA_1000 = 7
const ECON_PRETEND_ACT = Int32(2605); const ECON_SPEC_COST_ACT = Int32(2607)
const ECON_SPEC_REV_ACT = Int32(2608); const ECON_START_ACT = Int32(2609)
const ECON_FIX_PCT = 1; const ECON_FIX_HRV = 2; const ECON_VAR_PCT = 3; const ECON_VAR_HRV = 4
const ECON_NEAR_ZERO_C = 0.01f0

"libgcc __powisf2: Float32 x^m by square-and-multiply (gfortran REAL**INTEGER)."
@inline function _ec_powi(x::Float32, m::Integer)
    n = m < 0 ? -Int(m) : Int(m)
    y = isodd(n) ? x : 1f0
    n >>= 1
    while n != 0
        x = x * x
        isodd(n) && (y = y * x)
        n >>= 1
    end
    return m < 0 ? 1f0 / y : y
end

# ---------------------------------------------------------------------------------------------
# ECIN keyword storage (called from kw_econ!)
# ---------------------------------------------------------------------------------------------

"ratesAndDurations (ecin.f:708): read the supplemental rate/duration record when `field` is '&'."
function _ec_rates_durations(field::AbstractString, kr)
    sc = EconSched()
    strip(field) == "&" || return sc
    rec = read_raw_line!(kr)
    rec = rpad(rec, 80)
    vals = Float32[]
    j = 1
    for _ in 1:(2 * ECON_MAX_RATES)
        f = strip(rec[j:j+4]); j += 5
        v = isempty(f) ? 0.0 : tryparse(Float64, f)
        (v === nothing || v < -100) && return EconSched()      # read error ⇒ ALL rate changes ignored
        push!(vals, Float32(v))
    end
    n = 0
    for i in 1:2:(2 * ECON_MAX_RATES)
        d = vals[i + 1]
        if d >= 1f0
            n += 1
            sc.rates[n] = vals[i]
            sc.durs[n] = Int32(trunc(d))
        end
    end
    return sc
end

_ec_calc!(s) = (ec = s.econ; ec.calc === nothing && (ec.calc = EconCalc(nspecies(s.variant))); ec.calc)

"""Register an ECON Event-Monitor activity (addEvent, ecin.f:781). Normal path: date = max(minDate, field 1).
Inside an IF/THEN block (OPMODE LMODE, ecin.f:785-789) field 1 is a WAIT time with no minimum — captured as a
template on the conditional (`c.if_capture`) and dated IY(ICYC)+wait when the condition fires (econ_evmon!)."""
function _ec_add_event!(c::EconCalc, code::Int32, mindate::Int, date::Float32, params::Vector{Float32})
    if c.if_capture !== nothing
        push!(c.if_capture, EconEvent(code, Int32(trunc(Int, date)), Int32(0), params, Int32(0)))
        return
    end
    push!(c.events, EconEvent(code, Int32(max(mindate, trunc(Int, date))), Int32(length(c.events) + 1), params, Int32(0)))
end

"The ECON event activity codes an IF/THEN block may carry (ecin.f addEvent callers)."
const ECON_EVENT_ACTS = (ECON_PRETEND_ACT, ECON_SPEC_COST_ACT, ECON_SPEC_REV_ACT, ECON_START_ACT)

"""
    econ_evmon!(s)

EVMON phase 1 (grincr.f:260, before ECSTATUS(…,0)) for the ECON event templates carried by IF/THEN blocks:
when a block's condition is true this cycle, each template is registered dated IY(ICYC)+wait (OPSCHD,
opadd.f `IDATE(IMGL)=IDT+IDATE(II)`). Same condition evaluation (and per-cycle re-firing) as the thinning
activities of IF blocks (cuts!); idempotent within a cycle.
"""
function econ_evmon!(s::StandState)
    c = s.econ.calc
    conds = s.control.conditionals
    isempty(conds) && return
    any(cd -> any(a -> a.icflag in ECON_EVENT_ACTS, cd.acts), conds) || return
    yr = Int(current_cycle_year(s)); fvscyc = Int(s.control.cycle) + 1
    ctx = EventCtx(fvscyc, yr, s)
    for (cd, nm, ast) in s.control.compute_defs
        _compute_due(Int(cd), s, yr, fvscyc) && (s.control.compute_vars[nm] = eval_event(ast, ctx))
    end
    for (ci, cd) in enumerate(conds)
        any(a -> a.icflag in ECON_EVENT_ACTS, cd.acts) || continue
        (ci, yr) in c.ev_fired && continue
        eval_event(cd.cond, ctx) != 0f0 || continue
        push!(c.ev_fired, (ci, yr))
        for a in cd.acts
            a.icflag in ECON_EVENT_ACTS || continue
            np = a.icflag == ECON_START_ACT ? 3 : 1
            push!(c.events, EconEvent(a.icflag, Int32(yr + Int(a.year)), Int32(length(c.events) + 1),
                                      Float32[a.params[k] for k in 1:np], Int32(0)))
        end
    end
    return
end

"""
    econ_keyword!(s, name, r, kr) -> Bool

Faithful ECIN (ecin.f) storage of one ECON-block keyword into `s.econ.calc`. Returns true if the
keyword was recognized (including ignored-with-warning cases).
"""
function econ_keyword!(s::StandState, k::AbstractString, r, kr)
    c = _ec_calc!(s)
    f(i) = i <= length(r.values) ? Float32(r.values[i]) : 0f0
    nb(i) = i <= length(r.present) && r.present[i]
    cf(i) = i <= length(r.fields) ? r.fields[i] : ""
    okunit(kw, u) = kw == "PLANTCST" ? (u == ECON_PER_ACRE || u == ECON_TPA_1000) :
                    kw in ("HRVVRCST", "PCTVRCST") ? (u in (ECON_TPA_U, ECON_BF_U, ECON_FT3_U)) :
                    kw == "HRVRVN" ? (u in (ECON_TPA_U, ECON_BF_U, ECON_BFLOG_U, ECON_FT3_U, ECON_FT3LOG_U)) :
                    kw == "PCTSPEC" ? (u in (ECON_BF_U, ECON_FT3_U)) : false
    if k == "ANNUCST" || k == "ANNURVN"
        v = k == "ANNUCST" ? c.ann_cost : c.ann_rev
        (length(v) == ECON_MAX_KEYWORDS || f(1) <= 0f0) && return true
        push!(v, EconAmtKw(f(1), _ec_rates_durations(cf(3), kr)))
    elseif k == "BURNCST" || k == "MECHCST"
        f(1) <= 0f0 && return true
        kw = EconAmtKw(f(1), _ec_rates_durations(cf(2), kr))
        k == "BURNCST" ? (c.burn = kw) : (c.mech = kw)
    elseif k == "HRVFXCST" || k == "PCTFXCST"
        v = k == "HRVFXCST" ? c.fix_hrv : c.fix_pct
        (length(v) == ECON_MAX_KEYWORDS || f(1) <= 0f0) && return true
        push!(v, EconAmtKw(f(1), _ec_rates_durations(cf(3), kr)))
    elseif k == "HRVVRCST"
        length(c.var_hrv) == ECON_MAX_KEYWORDS && return true
        (!nb(1) || f(1) < 0f0) && return true
        u = trunc(Int, f(2)); okunit(k, u) || return true
        sc = _ec_rates_durations(cf(5), kr)
        # ecin.f:203-205: rates → varHrvRate(varHrvCnt,:), DURATIONS → fixHrvDur(fixHrvCnt,:) (always
        # re-zeroed by ratesAndDurations). varHrvDur stays 0 ⇒ variable harvest costs never appreciate.
        push!(c.var_hrv, EconVarKw(f(1), Int32(u), f(3), f(4) > 0f0 ? f(4) : 999f0,
                                   EconSched(copy(sc.rates), zeros(Int32, ECON_MAX_RATES))))
        isempty(c.fix_hrv) || (c.fix_hrv[end].sched.durs = copy(sc.durs))
    elseif k == "PCTVRCST"
        length(c.var_pct) == ECON_MAX_KEYWORDS && return true
        f(1) <= 0f0 && return true
        u = trunc(Int, f(2)); okunit(k, u) || return true
        push!(c.var_pct, EconVarKw(f(1), Int32(u), f(3), f(4) > 0f0 ? f(4) : 999f0,
                                   _ec_rates_durations(cf(4), kr)))
    elseif k == "PCTSPEC"
        u = trunc(Int, f(2)); okunit(k, u) || return true
        c.pct_min_units = Int32(u)
        f(1) > 0f0 && (c.pct_min_dbh = f(1))
        f(3) > 0f0 && (c.pct_min_volume = f(3))
    elseif k == "PLANTCST"
        length(c.plant) == ECON_MAX_PLANT_COSTS && return true
        if !isempty(c.plant)
            c.plant[1].units == trunc(Int, f(2)) && return true           # identical units already set
        elseif f(1) <= 0f0
            return true
        elseif !okunit(k, trunc(Int, f(2)))
            return true
        end
        push!(c.plant, EconPlantKw(f(1), Int32(trunc(Int, f(2))), _ec_rates_durations(cf(3), kr)))
    elseif k == "HRVRVN"
        (!nb(1) || f(1) < 0f0) && return true
        u = trunc(Int, f(2)); okunit(k, u) || return true
        spfld = strip(cf(4))
        # SPDECD(4,…): blank/0/ALL ⇒ 0, numeric ⇒ index/group; an unresolvable alpha code ⇒ −999 (keyword ignored)
        spid = if isempty(spfld) || uppercase(spfld) == "ALL"
            0
        elseif tryparse(Float64, spfld) !== nothing
            species_selector(s, spfld)
        else
            (idx = first(resolve_species(spfld, s.variant, s.species, s.coef)); idx > 0 ? idx : -999)
        end
        spid == -999 && return true
        sc = _ec_rates_durations(cf(5), kr)
        price = f(1)
        (u == ECON_FT3_U || u == ECON_FT3LOG_U) && (price = f(1) / 100f0)
        (u == ECON_BF_U || u == ECON_BFLOG_U) && (price = f(1) / 1000f0)
        dia = f(3)
        nsp = size(c.rev, 1)
        dup(i) = any(r -> r.dia == dia, c.rev[i, u])
        assign!(i) = push!(c.rev[i, u], EconRevKw(price, dia, EconSched(copy(sc.rates), copy(sc.durs))))
        if spid == 0
            for i in 1:nsp
                length(c.rev[i, u]) >= ECON_MAX_KEYWORDS && continue
                c.has_rev_amt[i, u] && continue
                dup(i) && continue
                assign!(i)
            end
        elseif spid < 0
            g = -spid
            grp = 1 <= g <= length(s.control.sp_groups) ? s.control.sp_groups[g] : Int32[]
            for i in grp
                length(c.rev[i, u]) >= ECON_MAX_KEYWORDS && continue
                dup(i) && continue
                c.has_rev_amt[i, u] = true
                assign!(i)
            end
        else
            i = spid
            length(c.rev[i, u]) >= ECON_MAX_KEYWORDS && return true
            dup(i) && return true
            c.has_rev_amt[i, u] = true
            assign!(i)
        end
    elseif k == "NOTABLE"
        n = trunc(Int, f(1))
        n == 1 && (c.no_output_tables = true; c.no_log_stock_table = true)
        n == 2 && (c.no_log_stock_table = true)
    elseif k == "PRETEND"
        f(1) < 0f0 && return true
        _ec_add_event!(c, ECON_PRETEND_ACT, 1, f(1), Float32[nb(2) ? f(2) : 999f0])
    elseif k == "SPECCST" || k == "SPECRVN"
        (f(1) < 0f0 || f(3) <= 0f0) && return true
        _ec_add_event!(c, k == "SPECCST" ? ECON_SPEC_COST_ACT : ECON_SPEC_REV_ACT, 0, f(1), Float32[f(3)])
    elseif k == "STRTECON"
        f(1) < 0f0 && return true
        f4 = (f(2) <= 0f0 || f(3) > 0f0) ? 0f0 : f(4)
        _ec_add_event!(c, ECON_START_ACT, 1, f(1), Float32[f(2), f(3), f4])
        c.econ_start_year = Int32(9999)
    elseif k == "LBSCFV"
        # ecin.f:317-357: pounds per cubic foot by species (Tons column of FVS_EconHarvestValue). SPDECD(2,…).
        f(1) <= 0f0 && return true
        spfld = strip(cf(2))
        spid = if isempty(spfld) || uppercase(spfld) == "ALL"
            0
        elseif tryparse(Float64, spfld) !== nothing
            species_selector(s, spfld)
        else
            (idx = first(resolve_species(spfld, s.variant, s.species, s.coef)); idx > 0 ? idx : -999)
        end
        spid == -999 && return true
        if spid == 0
            for i in eachindex(c.lbs_ft3); c.lbs_ft3[i] <= 0f0 && (c.lbs_ft3[i] = f(1)); end
        elseif spid < 0
            g = -spid
            for i in (1 <= g <= length(s.control.sp_groups) ? s.control.sp_groups[g] : Int32[])
                c.lbs_ft3[i] > 0f0 && continue
                c.lbs_ft3[i] = f(1)
            end
        else
            c.lbs_ft3[spid] > 0f0 || (c.lbs_ft3[spid] = f(1))
        end
    else
        return false
    end
    return true
end

# ---------------------------------------------------------------------------------------------
# Activity table helpers (OPFIND/OPGET/OPDONE, OPGET2/OPDON2, OPSTUS/OPGET3 semantics)
# ---------------------------------------------------------------------------------------------

"OPEXPN date resolution: a date in 1..MAXCYC is a cycle number ⇒ IY(date)."
_ec_date(s, d::Integer) = (1 <= d <= MAXCYC) ? Int(cycle_year_at(s.control, d - 1)) : Int(d)

"OPCYCL: the cycle index (1-based ICYC) an activity dated `yr` is assigned to (first cycle with yr < IY(ICYC+1))."
function _ec_cycle_of(s, yr::Int)
    n = Int(s.control.ncycle_eff > 0 ? s.control.ncycle_eff : s.control.ncycle)
    for i in 1:n
        yr < Int(cycle_year_at(s.control, i)) && return i
    end
    return 0
end

_ec_sorted(evs) = sort(evs; by = e -> (e.date, e.seq))

"OPFIND(code) for the current cycle: not-done activities assigned to cycle `icyc`, date/seq ordered."
function _ec_opfind(s, c::EconCalc, code::Int32, icyc::Int)
    [e for e in _ec_sorted([x for x in c.events if x.code == code])
     if e.done == 0 && _ec_cycle_of(s, _ec_date(s, e.date)) == icyc]
end

"Generic OPSTUS/OPGET3 view over a list of (date, done, params) activities of one code."
struct _EcAct
    date::Int
    done::Int
    params::Vector{Float32}
end

# ---------------------------------------------------------------------------------------------
# ECSETP (fvs.f:148) — before cycling
# ---------------------------------------------------------------------------------------------
function econ_setp!(s::StandState)
    ec = s.econ
    (ec === nothing || !ec.active) && return
    c = _ec_calc!(s)
    if c.econ_start_year == -9999                                        # no valid STRTECON keyword
        push!(c.events, EconEvent(ECON_START_ACT, Int32(cycle_year_at(s.control, 0)),
                                  Int32(length(c.events) + 1), Float32[0, 0, 0], Int32(0)))
        c.econ_start_year = Int32(9999)
    end
    for i in axes(c.rev, 1), j in 1:ECON_MAX_REV_UNITS                   # RDPSRT descending-diameter index
        v = c.rev[i, j]
        c.rev_idx[i, j] = isempty(v) ? Int[] : sortperm([r.dia for r in v]; rev = true)
    end
    return
end

# ---------------------------------------------------------------------------------------------
# ECSTATUS (grincr.f:273 before CUTS, :370 after CUTS)
# ---------------------------------------------------------------------------------------------
function econ_status!(s::StandState, icyc::Int, before_cuts::Int)
    ec = s.econ
    (ec === nothing || !ec.active || ec.calc === nothing) && return
    c = ec.calc
    iy = Int(cycle_year_at(s.control, icyc - 1))
    if c.econ_start_year > iy
        evs = _ec_opfind(s, c, ECON_START_ACT, icyc)
        if !isempty(evs)
            idt = 0; prm = Float32[0, 0, 0]
            for e in evs
                idt = _ec_date(s, e.date); prm = e.params
                e.done = Int32(idt == 0 ? 1 : idt)
            end
            c.econ_start_year = Int32(idt)
            c.pretend_start_year = c.econ_start_year
            c.discount_pct = prm[1]
            c.sev_input = prm[2]
            c.do_sev = prm[3] > 0f0
            ec.discount_rate = prm[1] / 100f0
        end
    end
    for e in _ec_opfind(s, c, ECON_PRETEND_ACT, icyc)
        idt = _ec_date(s, e.date)
        (idt + before_cuts > c.pretend_end_year) && (c.pretend_start_year = Int32(idt + before_cuts))
        c.pretend_end_year = Int32(round(Int, e.params[1], RoundNearestTiesAway) + idt + before_cuts - 1)
        d = idt + before_cuts
        e.done = Int32(d == 0 ? 1 : d)
    end
    if before_cuts == 0 && c.pretend_start_year >= c.econ_start_year
        c.is_pretend_active = c.pretend_start_year <= iy && c.pretend_end_year >= iy
    end
    return
end

"""
    econ_cycle_start!(s)

ECSETP (once, before the first cycle) + ECSTATUS(…,0) for the current cycle. Must run BEFORE the cycle's
CUTS; `write_sum_file` applies CUTS ahead of `grow_cycle!`, so both call this. Idempotent within a cycle
(consumed STRTECON/PRETEND events are marked done; ECSETP registers the default start only once).
"""
function econ_cycle_start!(s::StandState)
    (s.econ === nothing || !s.econ.active) && return
    s.control.cycle == Int32(0) && econ_setp!(s)
    s.econ.calc === nothing || econ_evmon!(s)             # EVMON phase 1 (grincr.f:260) precedes ECSTATUS(…,0)
    econ_status!(s, Int(s.control.cycle) + 1, 0)
    return
end

# ---------------------------------------------------------------------------------------------
# ECHARV (cuts.f:1673) — per removed tree
# ---------------------------------------------------------------------------------------------
"getDiaGrp (echarv.f:146): index of the largest class ≤ value (descending-sorted index), or -1."
function _ec_dia_grp(c::EconCalc, sp::Int, u::Int, value::Float32)
    idx = c.rev_idx[sp, u]; v = c.rev[sp, u]
    isempty(idx) && return -1
    value >= v[idx[1]].dia && return idx[1]
    for i in 2:length(idx)
        value >= v[idx[i]].dia && return idx[i]
    end
    return -1
end

function econ_harv!(s::StandState, sp::Int, tree::Int, dbh::Float32, bf::Float32, ft3::Float32,
                    prem2::Float32, gross_space::Float32)
    ec = s.econ
    (ec === nothing || !ec.active || ec.calc === nothing) && return
    c = ec.calc
    iy = Int(current_cycle_year(s))
    (c.econ_start_year > iy || prem2 <= 0f0) && return
    harv = prem2 / gross_space
    c.dbh_sq = c.dbh_sq + (dbh * dbh * harv)
    c.harvest[ECON_TPA_U] += harv
    c.harvest[ECON_BF_U]  += bf * harv
    c.harvest[ECON_FT3_U] += ft3 * harv
    for (i, v) in enumerate(c.var_pct)
        if dbh >= v.lo && dbh < v.hi
            c.pct_bf[i] += bf * harv; c.pct_ft3[i] += ft3 * harv; c.pct_tpa[i] += harv
        end
    end
    for (i, v) in enumerate(c.var_hrv)
        if dbh >= v.lo && dbh < v.hi
            c.hrv_cost_bf[i] += bf * harv; c.hrv_cost_ft3[i] += ft3 * harv; c.hrv_cost_tpa[i] += harv
        end
    end
    if bf > 0f0
        logs = get(c.tree_logs_bf, tree, nothing)
        logs === nothing && (logs = get(ec.tree_log_bf, tree, nothing))
        if logs !== nothing && !isempty(_ec_log_iter(logs)) && !isempty(c.rev[sp, ECON_BFLOG_U])
            _ec_accum_logs!(c, sp, ECON_BFLOG_U, logs, bf, harv)
        elseif !isempty(c.rev[sp, ECON_BF_U])
            k = _ec_dia_grp(c, sp, ECON_BF_U, dbh)
            k > 0 && (c.rev_volume[sp, ECON_BF_U, k] += bf * harv)
        end
    end
    if ft3 > 0f0
        logs = get(c.tree_logs_ft3, tree, nothing)
        logs === nothing && (logs = get(ec.tree_log_ft3, tree, nothing))
        if logs !== nothing && !isempty(_ec_log_iter(logs)) && !isempty(c.rev[sp, ECON_FT3LOG_U])
            _ec_accum_logs!(c, sp, ECON_FT3LOG_U, logs, ft3, harv)
        elseif !isempty(c.rev[sp, ECON_FT3_U])
            k = _ec_dia_grp(c, sp, ECON_FT3_U, dbh)
            k > 0 && (c.rev_volume[sp, ECON_FT3_U, k] += ft3 * harv)
        end
    end
    if !isempty(c.rev[sp, ECON_TPA_U])
        k = _ec_dia_grp(c, sp, ECON_TPA_U, dbh)
        k > 0 && (c.rev_volume[sp, ECON_TPA_U, k] += harv)
    end
    return
end

_ec_cut_on(s) = s.econ !== nothing && s.econ.active && s.econ.calc !== nothing

"CUTS entry: snapshot PROB, zero the per-record PREM/yard pools, IND2 = identity (cuts.f:1600 DO 1700 prep)."
function econ_cuts_begin!(s::StandState)
    _ec_cut_on(s) || return
    c = s.econ.calc; n = s.trees.n
    c.cut_prob0 = Float32[s.trees.tpa[i] for i in 1:n]
    c.cut_prem = zeros(Float32, n); c.cut_dsng = zeros(Float32, n); c.cut_ssng = zeros(Float32, n)
    c.cut_ind2 = Int32.(1:n)
    return
end

"Each cut METHOD starts with IND2 = identity (LSPECL methods keep it; sorted methods overwrite via econ_cut_order!)."
econ_cut_method!(s::StandState) = (_ec_cut_on(s) && (s.econ.calc.cut_ind2 = Int32.(1:s.trees.n)); nothing)

"A sorted cut method's full-list RDPSRT priority order = FVS IND2 (cuts.f:1135)."
econ_cut_order!(s::StandState, order::AbstractVector) = (_ec_cut_on(s) && (s.econ.calc.cut_ind2 = Int32.(order)); nothing)

"Accumulate one method's PREM for record `i` + the YARDLOSS DSNG/SSNG pools (cuts.f:1382-1386)."
function econ_cut_accum!(s::StandState, i::Int, prem::Float32)
    _ec_cut_on(s) || return
    c = s.econ.calc
    i <= length(c.cut_prem) || return
    c.cut_prem[i] += prem
    pl = s.control.yardloss_prlost
    if pl > 0f0
        loss = prem * pl
        c.cut_dsng[i] = c.cut_dsng[i] + loss * s.control.yardloss_prdsng
        c.cut_ssng[i] = c.cut_ssng[i] + loss * (1f0 - s.control.yardloss_prdsng)
    end
    return
end

"""
    econ_cuts_replay!(s)

The CUTS DO-1700 loop's ECHARV calls (cuts.f:1600-1673): for I = 1..ITRN, IT = IND2(I); PREM = PROB−WK4
(skip < 0.00001); PREM2 = max(0, PREM−DSNG−SSNG); if PREM2 > 0.00001 call ECHARV(BFV(IT), DBH(IT),
CFVOLI, …) where CFVOLI = SCFV(I) (CS/LS/NE/SN) or MCFV(I) — the record at LOOP POSITION I, as FVS indexes it.
Run after the MINHARV gate (a canceled thin never reaches DO 1700) and before TREDEL compaction.
"""
function econ_cuts_replay!(s::StandState)
    _ec_cut_on(s) || return
    c = s.econ.calc; t = s.trees
    n = length(c.cut_prem)
    n == 0 && return
    east = s.variant isa Southern || s.variant isa Northeast || s.variant isa CentralStates ||
           s.variant isa LakeStates
    g = s.plot.gross_space; g <= 0f0 && (g = 1f0)
    for I in 1:n
        it = Int(c.cut_ind2[I])
        c.cut_prem[it] > 0f0 || continue
        prem = c.cut_prob0[it] - t.tpa[it]                  # PREM = PROB − WK4
        prem < 0.00001f0 && continue
        prem2 = max(0f0, prem - c.cut_dsng[it] - c.cut_ssng[it])
        prem2 > 0.00001f0 || continue
        cf = east ? t.saw_cuft_vol[I] : t.merch_cuft_vol[I]  # MCFV(I)/SCFV(I) — loop position I
        econ_harv!(s, Int(t.species[it]), it, t.dbh[it], t.bdft_vol[it], cf, prem2, g)
    end
    empty!(c.cut_prem)
    return
end

"ECVOL capture target (vols.f:339-426): the EconCalc when ECON is active with a unit-4/5 revenue record, else nothing."
function econ_log_capture(s::StandState)
    _ec_cut_on(s) || return nothing
    c = s.econ.calc
    any(!isempty, @view c.rev[:, ECON_BFLOG_U]) || any(!isempty, @view c.rev[:, ECON_FT3LOG_U]) || return nothing
    return c
end

"ECVOL (ecvol.f): record tree i's per-log (class DIB, gross vol) when BFV>0 (board) / MCFV>0 (cubic)."
function econ_log_store!(c, i::Int, bf::Float32, mcf::Float32, lbf, lft)
    c === nothing && return
    bf > 0f0 && lbf !== nothing ? (c.tree_logs_bf[i] = lbf) : delete!(c.tree_logs_bf, i)
    mcf > 0f0 && lft !== nothing ? (c.tree_logs_ft3[i] = lft) : delete!(c.tree_logs_ft3, i)
    return
end

# echarv.f:70-91 (board) / :102-123 (cubic): per-log volume bucketed by log-end DIB class. `logs` is the
# per-log list (dib => gross vol) stashed by compute_volumes!; defProp = Σlog / per-tree net volume.
function _ec_accum_logs!(c::EconCalc, sp::Int, u::Int, logs, treevol_net::Float32, harv::Float32)
    tv = 0f0
    for (_, v) in _ec_log_iter(logs); tv += v; end
    defprop = tv / treevol_net
    for (dib, v) in _ec_log_iter(logs)
        k = _ec_dia_grp(c, sp, u, Float32(dib))
        k > 0 && (c.rev_volume[sp, u, k] += (v * harv) / defprop)
    end
end
_ec_log_iter(logs::AbstractDict) = sort!(collect(logs); by = first)
# ECHARV walks logs `do while (logBfVol(treeId,logId) > 0.0)` — stop at the first zero-volume log.
function _ec_log_iter(logs::AbstractVector)
    k = findfirst(x -> !(x[2] > 0f0), logs)
    return k === nothing ? logs : logs[1:k-1]
end

# ---------------------------------------------------------------------------------------------
# ECCALC (gradd.f:239) — once per cycle, after establishment
# ---------------------------------------------------------------------------------------------

"calcAppreAmt (eccalc.f:1029) → (value, done)."
function _ec_appre(amt::Float32, sc::EconSched, t::Integer)
    val = amt; done = true
    sc.durs[1] == 0 && return (val, done)
    dur = 0
    for i in 1:ECON_MAX_RATES
        sc.durs[i] <= 0 && break
        dur += Int(sc.durs[i])
    end
    dur <= 0 && return (val, done)
    t < dur && (done = false)
    dt = Int(t)
    i = 1
    while i <= ECON_MAX_RATES && sc.durs[i] > 0
        if dt > sc.durs[i]
            val = val * _ec_powi(1f0 + sc.rates[i] / 100f0, sc.durs[i])
            dt -= Int(sc.durs[i])
        else
            val = val * _ec_powi(1f0 + sc.rates[i] / 100f0, dt)
            break
        end
        i += 1
    end
    return (val, done)
end

"computePV (eccalc.f:991)."
@inline _ec_pv(amt::Float32, t::Integer, rate::Float32) = amt <= 0f0 ? 0f0 : amt / _ec_powi(1f0 + rate, t)

"computePNV (eccalc.f:1006) over years 1..endTime."
function _ec_pnv(c::EconCalc, rate::Float32, endtime::Int)
    dc = 0f0; dr = 0f0
    for i in 1:endtime
        dc = dc + _ec_pv(c.undisc_cost[i], i - 1, rate)
        dr = dr + _ec_pv(c.undisc_rev[i], i, rate)
    end
    return dr - dc
end

"calcAppreSev (eccalc.f:1071)."
function _ec_appre_sev(c::EconCalc, amt::Float32, price::Float32, sc::EconSched, evnt::Int, endtime::Int)
    npv = 0f0
    at = evnt + Int(c.start_year) - Int(c.econ_start_year)
    dt = evnt
    local undisc::Float32
    if sc.durs[1] != 0
        done = false
        while !done
            v, done = _ec_appre(price, sc, at)
            u = amt * v
            npv = npv + u / _ec_powi(1f0 + c.rate, dt)
            at += endtime; dt += endtime
        end
        v, _ = _ec_appre(price, sc, at)
        undisc = amt * v
    else
        undisc = amt * price
    end
    factor = _ec_powi(1f0 + c.rate, endtime)
    return npv + (undisc * factor) / (factor - 1f0) / _ec_powi(1f0 + c.rate, dt)
end

function _ec_sev_hrv_costs(c::EconCalc, endtime::Int)
    tot = 0f0
    for (tm, typ, kwf, amt) in c.hrv_cst
        kw = Int(kwf); evnt = Int(tm); cost = 0f0
        units = kw <= length(c.var_hrv) ? Int(c.var_hrv[kw].units) : 0     # varHrvUnits(kw)
        if typ == ECON_FIX_PCT
            k = c.fix_pct[kw]; cost = _ec_appre_sev(c, amt, k.amt, k.sched, evnt, endtime)
        elseif typ == ECON_FIX_HRV
            k = c.fix_hrv[kw]; cost = _ec_appre_sev(c, amt, k.amt, k.sched, evnt, endtime)
        elseif typ == ECON_VAR_PCT || typ == ECON_VAR_HRV
            k = typ == ECON_VAR_PCT ? c.var_pct[kw] : c.var_hrv[kw]
            a = units == ECON_TPA_U ? amt : units == ECON_BF_U ? amt / 1000f0 :
                units == ECON_FT3_U ? amt / 100f0 : NaN32
            isnan(a) || (cost = _ec_appre_sev(c, a, k.amt, k.sched, evnt, endtime))
        end
        tot = tot + cost
    end
    return tot
end

function _ec_sev_hrv_revenues(c::EconCalc, endtime::Int)
    tot = 0f0
    for (tm, spf, uf, kf, amt) in c.hrv_rvn
        r = c.rev[Int(spf), Int(uf)][Int(kf)]
        tot = tot + _ec_appre_sev(c, amt, r.price, r.sched, Int(tm), endtime)
    end
    return tot
end

function _ec_init_saved!(c::EconCalc)
    c.cost_disc = 0f0; c.cost_undisc = 0f0
    empty!(c.hrv_cst); empty!(c.hrv_rvn)
    c.burn_cnt = 0; c.mech_cnt = 0; c.spec_cst_cnt = 0; c.spec_rvn_cnt = 0
    c.rev_disc = 0f0; c.rev_undisc = 0f0
    fill!(c.undisc_cost, 0f0); fill!(c.undisc_rev, 0f0)
end

function _ec_reset_cycle!(c::EconCalc)
    c.dbh_sq = 0f0; fill!(c.harvest, 0f0)
    fill!(c.hrv_cost_bf, 0f0); fill!(c.hrv_cost_ft3, 0f0); fill!(c.hrv_cost_tpa, 0f0)
    fill!(c.pct_bf, 0f0); fill!(c.pct_ft3, 0f0); fill!(c.pct_tpa, 0f0)
    fill!(c.rev_volume, 0f0)
end

"calcAnnCostRevSEV (eccalc.f:1178)."
function _ec_ann_sev!(c::EconCalc)
    c.sev_ann_cst = 0f0; c.sev_ann_rvn = 0f0
    for k in c.ann_cost; c.sev_ann_cst = c.sev_ann_cst + _ec_appre_sev(c, 1f0, k.amt, k.sched, 1, 1); end
    for k in c.ann_rev;  c.sev_ann_rvn = c.sev_ann_rvn + _ec_appre_sev(c, 1f0, k.amt, k.sched, 1, 1); end
end

"The ESTAB activities of `code` (PLANT 430 / BURN 491 / MECH 493) as (date, done-year, params)."
function _ec_estab_acts(s::StandState, code::Int32)
    out = _EcAct[]
    for (idx, a) in enumerate(s.control.schedule)
        a.icflag == code || continue
        yr = _ec_date(s, Int(a.year))
        # IACT(,4): PLANT is done when establish! booked it; BURNPREP/MECHPREP carry the status ESTAB's ESETPR /
        # NTALLY>1 cancel recorded (year done, −1 deleted, absent = pending) — see engine/establishment.jl.
        done = code == Int32(430) ? (Int32(yr) in s.estab.years_done ? yr : 0) :
               Int(get(s.estab.prep_status, idx, Int32(0)))
        prm = code == Int32(430) ? Float32[a.params[1], a.params[2]] : Float32[a.params[2]]
        push!(out, _EcAct(yr, done, prm))
    end
    return sort!(out; by = x -> x.date)
end

"OPSTUS: k-th activity (any status) dated in [y1,y2]; OPGET3: k-th DONE one."
_ec_opstus(acts, y1, y2, k) = (w = [a for a in acts if y1 <= a.date <= y2]; k <= length(w) ? w[k] : nothing)
_ec_opget3(acts, y1, y2, k) = (w = [a for a in acts if y1 <= a.date <= y2 && a.done != 0]; k <= length(w) ? w[k] : nothing)

"""
    econ_calc!(s, icyc)

ECCALC (eccalc.f) for cycle `icyc` (1-based): accrues the cycle's annual/harvest/special/establishment
cash flows, updates the discounted accumulators and PNV/IRR/B-C/RRR/SEV, and appends the FVS_EconSummary
row (DBSECSUM) to `s.econ.calc.rows`.
"""
function econ_calc!(s::StandState, icyc::Int)
    ec = s.econ
    (ec === nothing || !ec.active || ec.calc === nothing) && return
    c = ec.calc
    iy_c = Int(cycle_year_at(s.control, icyc - 1)); iy_n = Int(cycle_year_at(s.control, icyc))
    c.econ_start_year >= iy_n && return
    st = (harv_cst = 0f0, harv_rvn = 0f0, pct_cst = 0f0, sev_sum = 0f0, is_pct = false)
    L = Ref(st)
    pretend = c.is_pretend_active
    local begin_year::Int, end_year::Int, begin_time::Int, end_time::Int
    function process_start!()
        evs = _ec_opfind(s, c, ECON_START_ACT, icyc)
        isempty(evs) && return
        idt = 0; prm = Float32[0, 0, 0]
        for e in evs
            idt = _ec_date(s, e.date); prm = e.params
            e.done = Int32(idt == 0 ? 1 : idt)
        end
        if idt > begin_year
            begin_time = begin_year - Int(c.start_year) + 1
            end_time = idt - Int(c.start_year)
            end_year = idt - 1
            value_harvest!(); calc_econ!()
            _ec_reset_cycle!(c)
            L[] = (harv_cst = 0f0, harv_rvn = 0f0, pct_cst = 0f0, sev_sum = 0f0, is_pct = false)
        end
        _ec_init_saved!(c)
        c.discount_pct = prm[1]; ec.discount_rate = prm[1] / 100f0
        c.rate = c.discount_pct / 100f0
        c.sev_input = prm[2]
        c.do_sev = prm[3] > 0f0
        c.start_year = Int32(idt)
        begin_year = idt
        end_year = iy_n - 1
        if c.do_sev
            end_time = 1                                                 # calcAnnCostRevSEV sets endTime=1
            _ec_ann_sev!(c)
        end
    end
    function value_harvest!()
        cst_cnt = 0; rvn_cnt = 0
        l = L[]; harv_cst = l.harv_cst; harv_rvn = l.harv_rvn; pct_cst = l.pct_cst
        sev_sum = l.sev_sum; is_pct = l.is_pct
        if c.harvest[ECON_TPA_U] > 0f0
            t = begin_year - Int(c.econ_start_year) + 1
            evnt = begin_time
            if c.pct_min_units > 0
                if (c.harvest[c.pct_min_units] < c.pct_min_volume) ||
                   (sqrt(c.dbh_sq / c.harvest[ECON_TPA_U]) < c.pct_min_dbh)       # PCT
                    is_pct = true
                    for (i, k) in enumerate(c.fix_pct)
                        cost, _ = _ec_appre(k.amt, k.sched, t)
                        pct_cst = pct_cst + cost
                        c.do_sev && (push!(c.hrv_cst, (Float32(evnt), Float32(ECON_FIX_PCT), Float32(i), 1f0)); cst_cnt += 1)
                    end
                    for i in 1:length(c.var_pct)
                        # eccalc.f:227-241: VAR_PCT reads varHrvUnits(i)/varHrvAmt(i)/varHrvRate/Dur(i) (FVS quirk)
                        vh = i <= length(c.var_hrv) ? c.var_hrv[i] :
                             EconVarKw(0f0, Int32(0), 0f0, 999f0, EconSched())
                        amt = 0f0; units = 1f0
                        vh.units == ECON_TPA_U && (amt = c.pct_tpa[i]; units = 1f0)
                        vh.units == ECON_BF_U  && (amt = c.pct_bf[i];  units = 1000f0)
                        vh.units == ECON_FT3_U && (amt = c.pct_ft3[i]; units = 100f0)
                        cv, _ = _ec_appre(vh.amt, vh.sched, t)
                        cost = cv / units
                        if amt > 0f0 && cost > 0f0
                            pct_cst = pct_cst + amt * cost
                            c.do_sev && (push!(c.hrv_cst, (Float32(evnt), Float32(ECON_VAR_PCT), Float32(i), amt)); cst_cnt += 1)
                        end
                    end
                    c.undisc_cost[begin_time] = c.undisc_cost[begin_time] + pct_cst
                else                                                              # commercial harvest
                    for (i, k) in enumerate(c.fix_hrv)
                        cost, _ = _ec_appre(k.amt, k.sched, t)
                        harv_cst = harv_cst + cost
                        c.do_sev && (push!(c.hrv_cst, (Float32(evnt), Float32(ECON_FIX_HRV), Float32(i), 1f0)); cst_cnt += 1)
                    end
                    for (i, k) in enumerate(c.var_hrv)
                        amt = 0f0; cost = 0f0
                        cv, _ = _ec_appre(k.amt, k.sched, t)
                        if k.units == ECON_TPA_U
                            amt = c.hrv_cost_tpa[i]; cost = cv
                        elseif k.units == ECON_BF_U
                            amt = c.hrv_cost_bf[i]; cost = cv / 1000f0
                        elseif k.units == ECON_FT3_U
                            amt = c.hrv_cost_ft3[i]; cost = cv / 100f0
                        end
                        if amt > 0f0
                            harv_cst = harv_cst + amt * cost
                            c.do_sev && (push!(c.hrv_cst, (Float32(evnt), Float32(ECON_VAR_HRV), Float32(i), amt)); cst_cnt += 1)
                        end
                    end
                    for i in axes(c.rev, 1), j in 1:ECON_MAX_REV_UNITS, k in 1:length(c.rev[i, j])
                        rv = c.rev_volume[i, j, k]
                        if rv > 0f0
                            r = c.rev[i, j][k]
                            price, _ = _ec_appre(r.price, r.sched, t)
                            harv_rvn = harv_rvn + rv * price
                            c.do_sev && (push!(c.hrv_rvn, (Float32(evnt), Float32(i), Float32(j), Float32(k), rv)); rvn_cnt += 1)
                        end
                    end
                    c.undisc_cost[begin_time] = c.undisc_cost[begin_time] + harv_cst
                    c.undisc_rev[begin_time] = c.undisc_rev[begin_time] + harv_rvn
                end
            end
        end
        if c.do_sev && (!isempty(c.hrv_rvn) || !isempty(c.hrv_cst))
            sev_sum = sev_sum + _ec_sev_hrv_revenues(c, end_time) - _ec_sev_hrv_costs(c, end_time)
            if pretend
                for _ in 1:cst_cnt; pop!(c.hrv_cst); end
                for _ in 1:rvn_cnt; pop!(c.hrv_rvn); end
            end
        end
        L[] = (harv_cst = harv_cst, harv_rvn = harv_rvn, pct_cst = pct_cst, sev_sum = sev_sum, is_pct = is_pct)
    end
    function calc_econ!()
        l = L[]; sev_sum = l.sev_sum
        c.do_sev && (sev_sum = sev_sum - c.sev_ann_cst + c.sev_ann_rvn)
        for k in c.ann_cost, j in begin_time:end_time
            t = begin_year - Int(c.econ_start_year) + j - 1
            cost, _ = _ec_appre(k.amt, k.sched, t)
            c.undisc_cost[j] = c.undisc_cost[j] + cost
        end
        for k in c.ann_rev, j in begin_time:end_time
            t = begin_year - Int(c.econ_start_year) + j
            price, _ = _ec_appre(k.amt, k.sched, t)
            c.undisc_rev[j] = c.undisc_rev[j] + price
        end
        factor_sev(amt, evnt) = (fct = _ec_powi(1f0 + c.rate, end_time);
                                 ((amt * fct) / (fct - 1f0)) / _ec_powi(1f0 + c.rate, evnt))
        for (code, isc) in ((ECON_SPEC_COST_ACT, true), (ECON_SPEC_REV_ACT, false))
            # OPGET2/OPDON2 loop: every not-done activity dated in [beginAnalYear, endAnalYear]
            for e in _ec_sorted([x for x in c.events if x.code == code])
                d = _ec_date(s, e.date)
                (e.done == 0 && begin_year <= d <= end_year) || continue
                evnt = d - Int(c.start_year) + 1
                if isc
                    c.undisc_cost[evnt] = c.undisc_cost[evnt] + e.params[1]
                else
                    c.undisc_rev[evnt] = c.undisc_rev[evnt] + e.params[1]
                end
                e.done = Int32(d == 0 ? 1 : d)
                c.do_sev && (isc ? (c.spec_cst_cnt += 1) : (c.spec_rvn_cnt += 1))
            end
            cnt = isc ? c.spec_cst_cnt : c.spec_rvn_cnt
            if cnt > 0 && c.do_sev
                acts = [_EcAct(_ec_date(s, x.date), Int(x.done), x.params)
                        for x in _ec_sorted([x for x in c.events if x.code == code])]
                k = 0
                while true
                    k += 1
                    a = _ec_opstus(acts, Int(c.start_year), end_year, k)
                    a === nothing && break
                    if a.done > 0
                        g = _ec_opget3(acts, Int(c.start_year), end_year, k)
                        g === nothing && continue
                        evnt = a.done - Int(c.start_year) + 1
                        v = factor_sev(g.params[1], evnt)
                        sev_sum = isc ? sev_sum - v : sev_sum + v
                    end
                end
            end
        end
        # establishment costs: PLANT (430), MECH (493), BURN (491)
        if !isempty(c.plant)
            acts = _ec_estab_acts(s, Int32(430))
            prev = -1; k = 0
            while true
                k += 1
                a = _ec_opstus(acts, begin_year, end_year, k); a === nothing && break
                a.done > 0 || continue
                g = _ec_opget3(acts, begin_year, end_year, k); g === nothing && continue
                t = a.done - Int(c.econ_start_year) + 1
                evnt = a.done - Int(c.start_year) + 1
                for p in c.plant
                    amt = 0f0
                    if p.units == ECON_PER_ACRE
                        prev != a.done && (prev = a.done; amt = 1f0)
                    elseif p.units == ECON_TPA_1000
                        amt = g.params[2] / 1000f0
                    end
                    if amt > 0f0
                        cost, _ = _ec_appre(p.amt, p.sched, t)
                        c.undisc_cost[evnt] = c.undisc_cost[evnt] + amt * cost
                    end
                end
            end
            if c.do_sev
                prev = -1; k = 0
                while true
                    k += 1
                    a = _ec_opstus(acts, Int(c.start_year), end_year, k); a === nothing && break
                    a.done > 0 || continue
                    g = _ec_opget3(acts, Int(c.start_year), end_year, k); g === nothing && continue
                    evnt = a.done - Int(c.start_year) + 1
                    for p in c.plant
                        amt = 0f0
                        if p.units == ECON_PER_ACRE
                            prev != a.done && (prev = a.done; amt = 1f0)
                        elseif p.units == ECON_TPA_1000
                            amt = g.params[2] / 1000f0
                        end
                        amt > 0f0 && (sev_sum = sev_sum - _ec_appre_sev(c, amt, p.amt, p.sched, evnt, end_time))
                    end
                end
            end
        end
        for (code, kw, iscnt) in ((Int32(493), c.mech, :mech), (Int32(491), c.burn, :burn))
            if kw.amt > 0f0
                acts = _ec_estab_acts(s, code); k = 0
                while true
                    k += 1
                    a = _ec_opstus(acts, begin_year, end_year, k); a === nothing && break
                    a.done > 0 || continue
                    g = _ec_opget3(acts, begin_year, end_year, k); g === nothing && continue
                    t = a.done - Int(c.econ_start_year) + 1
                    evnt = a.done - Int(c.start_year) + 1
                    amt = g.params[1] / 100f0
                    cost, _ = _ec_appre(kw.amt, kw.sched, t)
                    c.undisc_cost[evnt] = c.undisc_cost[evnt] + amt * cost
                    c.do_sev && (iscnt === :mech ? (c.mech_cnt += 1) : (c.burn_cnt += 1))
                end
            end
            cnt = iscnt === :mech ? c.mech_cnt : c.burn_cnt
            if cnt > 0 && c.do_sev
                acts = _ec_estab_acts(s, code); k = 0
                while true
                    k += 1
                    a = _ec_opstus(acts, Int(c.start_year), end_year, k); a === nothing && break
                    a.done > 0 || continue
                    g = _ec_opget3(acts, Int(c.start_year), end_year, k); g === nothing && continue
                    evnt = a.done - Int(c.start_year) + 1
                    sev_sum = sev_sum - _ec_appre_sev(c, g.params[1] / 100f0, kw.amt, kw.sched, evnt, end_time)
                end
            end
        end
        for i in begin_time:end_time
            c.cost_undisc = c.cost_undisc + c.undisc_cost[i]
            c.rev_undisc = c.rev_undisc + c.undisc_rev[i]
            c.cost_disc = c.cost_disc + _ec_pv(c.undisc_cost[i], i - 1, c.rate)
            c.rev_disc = c.rev_disc + _ec_pv(c.undisc_rev[i], i, c.rate)
        end
        pnv = c.rev_disc - c.cost_disc
        forest = nothing; reprod = nothing; sev = nothing; given = nothing
        if c.sev_input > 0f0
            dsev = _ec_pv(c.sev_input, end_time, c.rate)
            forest = pnv + dsev
            reprod = pnv + dsev - c.sev_input
            given = c.sev_input
        elseif c.do_sev
            sev = sev_sum
        end
        irr_ok, irr = _ec_irr(c, c.rate, pnv, end_time)
        irr = irr * 100f0
        bc = nothing; rrr = nothing
        if c.cost_disc > ECON_NEAR_ZERO_C
            bc = c.rev_disc / c.cost_disc
            if c.rev_disc > ECON_NEAR_ZERO_C
                rrr = 100f0 * (FMath.fpow(c.rev_disc / c.cost_disc, 1f0 / Float32(end_time)) * (1f0 + c.rate) - 1f0)
            end
        end
        ft3tot = 0f0; bftot = 0f0
        if c.harvest[ECON_TPA_U] > 0f0 && !l.is_pct
            for i in axes(c.rev, 1), j in 1:ECON_MAX_REV_UNITS, k in 1:length(c.rev[i, j])
                (j == ECON_FT3_U || j == ECON_FT3LOG_U) && (ft3tot = ft3tot + c.rev_volume[i, j, k])
                (j == ECON_BF_U || j == ECON_BFLOG_U) && (bftot = bftot + c.rev_volume[i, j, k])
            end
        end
        # Harvest Volume/Value report → FVS_EconHarvestValue (eccalc.f:737-855; DBSECHARV_insert needs IDBSECON=2).
        if c.dbs_econ == 2 && c.harvest[ECON_TPA_U] > 0f0 && !l.is_pct
            t = begin_year - Int(c.econ_start_year) + 1
            nintf(x) = round(Int, x, RoundNearestTiesAway)
            for i in axes(c.rev, 1), j in 1:ECON_MAX_REV_UNITS
                idx = c.rev_idx[i, j]; v = c.rev[i, j]
                for k in length(idx):-1:1
                    kk = idx[k]; rv = c.rev_volume[i, j, kk]
                    rv <= 0f0 && continue
                    diamax = k == 1 ? 999.9f0 : v[idx[k - 1]].dia
                    mindia = maxdia = mindbh = maxdbh = -1f0
                    if j < ECON_BFLOG_U
                        maxdbh = diamax; mindbh = v[kk].dia
                    else
                        maxdia = diamax; mindia = v[kk].dia
                    end
                    price, _ = _ec_appre(v[kk].price, v[kk].sched, t)
                    tpacut = tpaval = tons = ft3v = ft3val = bfv = bfval = -1
                    amt = 0f0
                    if j == ECON_TPA_U
                        tpacut = nintf(rv); amt = rv * price; tpaval = nintf(amt)
                    elseif j == ECON_FT3_U || j == ECON_FT3LOG_U
                        if c.lbs_ft3[i] > 0f0
                            amt = rv * (c.lbs_ft3[i] / 2000f0); tons = nintf(amt)
                        end
                        ft3v = nintf(rv); amt = rv * price; ft3val = nintf(amt)
                    else
                        bfv = nintf(rv); amt = rv * price; bfval = nintf(amt)
                    end
                    push!(c.hv_rows, (year = begin_year, sp = i, min_dia = mindia, max_dia = maxdia,
                                      min_dbh = mindbh, max_dbh = maxdbh, tpa_cut = tpacut, tpa_value = tpaval,
                                      tons = tons, ft3_vol = ft3v, ft3_value = ft3val, bf_vol = bfv,
                                      bf_value = bfval, total = nintf(amt)))
                end
            end
        end
        if c.dbs_econ > 0 || !c.no_output_tables
            push!(c.rows, (year = begin_year, period = end_time, pretend = pretend ? "YES" : "NO ",
                           cost_undisc = c.cost_undisc, rev_undisc = c.rev_undisc,
                           cost_disc = c.cost_disc, rev_disc = c.rev_disc, pnv = pnv,
                           irr = irr_ok ? irr : nothing, bc = bc, rrr = rrr, sev = sev,
                           forest = forest, reprod = reprod,
                           ft3 = round(Int, ft3tot, RoundNearestTiesAway), bf = round(Int, bftot, RoundNearestTiesAway),
                           rate = c.discount_pct, given = given))
        end
        L[] = (harv_cst = l.harv_cst, harv_rvn = l.harv_rvn, pct_cst = l.pct_cst, sev_sum = sev_sum, is_pct = l.is_pct)
    end
    L[] = (harv_cst = 0f0, harv_rvn = 0f0, pct_cst = 0f0, sev_sum = 0f0, is_pct = false)
    if c.is_first_econ
        c.is_first_econ = false
        _ec_init_saved!(c)
        c.start_year = c.econ_start_year
        begin_year = Int(c.econ_start_year)
        end_year = iy_n - 1
        c.rate = c.discount_pct / 100f0
        end_time = 0
        if c.do_sev
            end_time = 1
            _ec_ann_sev!(c)
        end
        process_start!()
    else
        begin_year = iy_c
        end_year = iy_n - 1
        end_time = end_year - Int(c.start_year) + 1
        process_start!()
    end
    begin_time = begin_year - Int(c.start_year) + 1
    end_time = end_year - Int(c.start_year) + 1
    value_harvest!()
    calc_econ!()
    l = L[]
    if pretend && c.harvest[ECON_TPA_U] > 0f0
        c.undisc_cost[begin_time] = c.undisc_cost[begin_time] - (l.harv_cst + l.pct_cst)
        c.undisc_rev[begin_time] = c.undisc_rev[begin_time] - l.harv_rvn
        c.cost_undisc = c.cost_undisc - (l.harv_cst + l.pct_cst)
        c.rev_undisc = c.rev_undisc - l.harv_rvn
        c.cost_disc = c.cost_disc - _ec_pv(l.harv_cst + l.pct_cst, begin_time - 1, c.rate)
        c.rev_disc = c.rev_disc - _ec_pv(l.harv_rvn, begin_time, c.rate)
    end
    _ec_reset_cycle!(c)
    return
end

"computeIRR (eccalc.f:916) → (irrCalculated, irr as a decimal)."
function _ec_irr(c::EconCalc, rate::Float32, pnv::Float32, endtime::Int)
    (c.cost_disc < ECON_NEAR_ZERO_C || c.rev_disc < ECON_NEAR_ZERO_C) && return (false, 0f0)
    abs(pnv) < ECON_NEAR_ZERO_C && return (true, rate)
    tol = 0.00001f0; inc = 0.0001f0
    ra = rate; na = pnv
    if pnv >= 0f0
        while na > 0f0 && ra <= 0.51f0
            ra = ra + inc; na = _ec_pnv(c, ra, endtime)
        end
        (na > 0f0 && ra > 0.51f0) && return (true, ra)
        rb = ra - inc; rc = (ra + rb) / 2f0; nc = _ec_pnv(c, rc, endtime)
        while abs(nc) > ECON_NEAR_ZERO_C && abs(ra - rb) > tol
            nc > 0f0 ? (rb = rc) : (ra = rc)
            rc = (ra + rb) / 2f0; nc = _ec_pnv(c, rc, endtime)
        end
        return (true, rc)
    else
        while na < 0f0 && ra >= 0f0
            ra = ra - inc; na = _ec_pnv(c, ra, endtime)
        end
        (na < 0f0 && ra < 0f0) && return (true, ra)
        rb = ra + inc; rc = (ra + rb) / 2f0; nc = _ec_pnv(c, rc, endtime)
        while abs(nc) > ECON_NEAR_ZERO_C && abs(ra - rb) > tol
            nc < 0f0 ? (rb = rc) : (ra = rc)
            rc = (ra + rb) / 2f0; nc = _ec_pnv(c, rc, endtime)
        end
        return (true, rc)
    end
end
