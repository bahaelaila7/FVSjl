# =============================================================================
# simulate.jl — the projection cycle loop (FVS / TREGRO orchestration)
#
# Ported from: base/fvs.f (cycle loop) + base/tregro.f (per-cycle driver).
#
# After initialization, FVS projects the stand cycle by cycle: recompute stand
# density, grow diameters and heights, apply mortality and regeneration, then
# recompute statistics. This is the skeleton; mortality (C4), regeneration (C4)
# and volume/output (C5) are wired in as those chunks land.
# =============================================================================

"""
    setup_growth!(state)

One-time growth setup after initialization (the FVS LSTART pass): compute the
site-dependent diameter-growth constants (DGCONS) and run the diameter-growth
calibration against the input measured growth (COR). Needs density set first.
"""
function setup_growth!(s::StandState)
    build_cycle_schedule!(s)             # CYCLEAT/TIMEINT → cycle-boundary year array (IY)
    dub_missing_heights!(s)              # CRATET — dub HT=0 / resolve broken-top NORMHT
    apply_growth_input_types!(s)         # GROWTH IDG/IHTG=1/3 — past DBH/HT field ⇒ increment
    setup_volume_equations!(s)           # VOLEQDEF — per-species NVEL equation ids
    isempty(s.control.sp_bf_vol_eq) && (s.control.sp_bf_vol_eq = copy(s.species.vol_eq)) # VEQNNB = default
    apply_voleqnum_overrides!(s)         # VOLEQNUM — user overrides of those equation ids (cubic only)
    compute_forest_type!(s)              # FORTYP — needed by dgf!'s forest-type term
    compute_density!(s)
    sdi_max_check!(s)                     # SDICHK — reset species SDImax if over-dense
    dgcons!(s)
    # LSTART calibration uses the input measured DG. SCALE = YR/FINT (dgdriv.f:325): YR = 5 (the SN
    # base measurement period), FINT = the GROWTH keyword's DIAMETER measurement period (default 5 →
    # SCALE = 1). A FINT≠5 means the input DG increment spans FINT years, so the measured DDS is
    # rescaled to the 5-yr basis the model predicts. (The projection-cycle length is a SEPARATE FINT,
    # threaded into diameter_growth! by TIMEINT; this scale is the input-measurement period only.)
    dfint = s.control.growth_fint
    calibrate_diameter_growth!(s; scale = dfint > 0f0 ? 5f0 / dfint : 1f0)
    return s
end

"""
    build_cycle_schedule!(s)

Build the cycle-boundary year array (`control.cycle_year`, the FVS IY) from NUMCYCLE +
TIMEINT + CYCLEAT, mirroring base/fvs.f:106-135. Each cycle's length is its TIMEINT
per-cycle override (`cycle_lengths[k]`) or the uniform `control.year` (default 5); the
lengths are cumulated onto the inventory year (`cycle_year[1]`, INVYEAR) to give the
calendar year at each boundary. CYCLEAT then inserts each requested year as a NEW boundary
strictly inside the run (never extending the end or moving the start), bumping the effective
cycle count `ncycle_eff` (capped at MAXCYC). Idempotent: recomputes purely from the
immutable inputs (cycle_year[1], ncycle, cycle_lengths, cycleat_years), so the per-stand
double-call from setup is safe. For uniform cycles `cycle_year[k+1] == cycle_year[1] + k·per`
exactly, so routing the year derivations through this array is bit-exact for snt01.
"""
function build_cycle_schedule!(s::StandState)
    c = s.control
    ncyc = Int(c.ncycle); ncyc < 1 && (ncyc = 1); ncyc > MAXCYC && (ncyc = MAXCYC)
    per = round(Int, c.year); per < 1 && (per = 5)
    iy = c.cycle_year
    @inbounds for k in 2:MAXCY1                         # cumulate per-cycle lengths → boundary years
        len = c.cycle_lengths[k] > 0 ? Int(c.cycle_lengths[k]) : per
        iy[k] = iy[k-1] + Int32(len)
    end
    for yr in c.cycleat_years                           # CYCLEAT insert (fvs.f:116-135)
        (yr <= iy[1] || yr >= iy[ncyc+1]) && continue   # don't extend end / move start
        @inbounds for i in 1:ncyc
            if iy[i] < yr < iy[i+1]
                ncyc += 1; ncyc > MAXCYC && (ncyc = MAXCYC)
                for k in ncyc+1:-1:i+2; iy[k] = iy[k-1]; end
                iy[i+1] = Int32(yr)
                break
            end
        end
    end
    c.ncycle_eff = Int32(ncyc)
    return s
end

"""
Calendar year at the start of cycle `cyc` (0-based) from the IY schedule (build_cycle_schedule!).
Falls back to the uniform derivation `cycle_year[1] + cyc·per` when the schedule has not been
built yet (a boundary of 0), so direct callers that skip `setup_growth!` still get the right year.
"""
function cycle_year_at(c::Control, cyc::Integer)
    y = Int(c.cycle_year[cyc + 1])
    y > 0 && return y
    per = round(Int, c.year); per < 1 && (per = 5)
    return Int(c.cycle_year[1]) + Int(cyc) * per
end
"Calendar year at the start of the current cycle (`control.cycle`)."
current_cycle_year(s::StandState) = cycle_year_at(s.control, Int(s.control.cycle))
"Length in years of cycle `cyc` (0-based) = next boundary − this boundary."
cycle_period_at(c::Control, cyc::Integer) = Int(c.cycle_year[cyc + 2] - c.cycle_year[cyc + 1])

"""
    compute_density!(state)

Recompute the per-cycle stand density quantities the growth models read:
basal area, average dominant height (AVH), and per-point basal area (PTBAA).
"""
function compute_density!(s::StandState)
    s.plot.basal_area = stand_ba(s)
    s.plot.avg_height = stand_top_height(s)
    point_basal_area!(s)
    stand_pct!(s)                      # PCT = stand BA percentile (for DGF competition)
    return s
end

"Number of early cycles that use deterministic record tripling (FVS ICL4)."
const TRIPLE_CYCLE_LIMIT = 2

"""
    fertilizer_growth!(s; fint)

FERTILIZE / FFERT (ffert.f): the 200-lb-N fertilizer response — a multiplicative boost to each
tree's squared-diameter change (DDS) and height growth for up to 10 years after application,
applied for the `iflen` of the cycle's years that fall in that window and scaled by the application
efficacy. Carries over across cycles via `ifert_date`/`ifert_eff`. No-op until a FERTILIZE keyword
fires (so default stands are unchanged). SN is outside the model's calibrated DF/GF range — Fortran
warns but still applies the (species-agnostic) factor, so we match.
"""
function fertilizer_growth!(s::StandState; fint::Float32 = 5f0)
    c = s.control
    (isempty(c.fertilize_events) && c.ifert_date < 0) && return s
    yr = cycle_year_at(c, Int(c.cycle))   # IY schedule (TIMEINT/CYCLEAT-aware)
    @inbounds for ev in c.fertilize_events           # a fert scheduled this cycle becomes active (OPDONE)
        Int(ev.year) == yr && (c.ifert_date = Int32(yr); c.ifert_eff = ev.params[1])
    end
    c.ifert_date < 0 && return s
    ifint  = round(Int, fint)
    ifstrt = yr - Int(c.ifert_date)
    ifstrt > 10 && return s                          # > 10 yr since application ⇒ effect gone
    iflen  = min(ifstrt + ifint, 10) - ifstrt        # years of fertilizer effect within this cycle
    iflen <= 0 && return s
    feff = c.ifert_eff
    t = s.trees; ba = s.plot.basal_area
    ba_a = s.calib.bark_a; ba_b = s.calib.bark_b
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        dib = d * bark_ratio(ba_a, ba_b, t.species[i], d)
        bal = (1f0 - t.crown_ratio[i] / 100f0) * ba   # basal area in larger trees (PCT in crown_ratio)
        rdds = exp(0.1108f0 * log(d) + 0.003004f0 * bal / log(d + 1f0))
        rdds > 2.6f0 && (rdds = 2.6f0)
        dg  = t.diam_growth[i]
        dds = 2f0 * dib * dg + dg * dg                # squared-diameter change this cycle
        ddsit = (dds / 5f0) * (rdds * iflen * feff + ifint - iflen)
        t.diam_growth[i] = sqrt(dib * dib + ddsit * 5f0 / fint) - dib
        htgit = (t.ht_growth[i] / 5f0) * (1.1626f0 * iflen * feff + ifint - iflen)
        t.ht_growth[i]   = htgit * 5f0 / fint
    end
    return s
end

"""
    _maybe_burn!(s, fint) -> fire_mort

Run a scheduled SIMFIRE if this cycle's year matches `fire.fire_year` (FMMAIN). Operates
on the current (post-MORTS, post-TRIPLE) records at cycle-start dimensions, and returns
the periodic mortality VOLUME of the fire-killed TPA (each record's lost TPA × its cycle-
start cubic volume), for the caller to add to OMORT. No-op (returns 0) otherwise.
"""
function _maybe_burn!(s::StandState, fint::Float32)::Float32
    (s.fire === nothing || !s.fire.active || s.fire.fire_year == 0) && return 0f0
    yr = current_cycle_year(s)   # IY schedule (TIMEINT/CYCLEAT-aware)
    yr == Int(s.fire.fire_year) || return 0f0
    t = s.trees
    pre_tpa = Float32[t.tpa[i] for i in 1:t.n]
    pre_cfv = Float32[t.cuft_vol[i] for i in 1:t.n]
    fmburn!(s; atemp = s.fire.atemp, wind = s.fire.swind, fmois = Int(s.fire.fmois),
            psburn = s.fire.psburn, mortcode = Int(s.fire.mortcode),
            burnseas = Int(s.fire.burnseas), flmult = s.fire.flmult, crburn = s.fire.crburn,
            year = yr)
    fm = 0f0
    @inbounds for i in 1:length(pre_tpa)
        fm += (pre_tpa[i] - t.tpa[i]) * pre_cfv[i]
    end
    s.fire.fire_year = Int32(0)            # one-shot
    compute_density!(s)
    return fm
end

"""
    grow_cycle!(state; fint=5f0) -> (; accretion, mortality)

Advance the stand by one growth cycle: recompute density, grow diameters/heights
(with record tripling), apply mortality, update dimensions, and recompute volumes.
Returns the period's per-acre cubic accretion and mortality (OACC/OMORT,
vols.f:190 / update.f:60): accretion = Σ(newCFV−oldCFV)·survivingTPA / fint,
mortality = Σ killedTPA·oldCFV / fint, both ÷ gross area. Requires the cycle-start
volumes to be present in `trees.cuft_vol` (run `compute_volumes!` once at setup).
"""
function grow_cycle!(s::StandState; fint::Float32 = 5f0)
    compute_density!(s)
    apply_setsite!(s)                                      # SETSITE (act 120): mid-run site change (RCON), before growth
    compressed = apply_compress!(s)                        # COMPRESS (act 250): cluster records → NCLAS (suppresses tripling)
    # ECON: zero the cycle's harvest accumulators; cuts!/_log_cut! values each removed tree.
    econ_on = s.econ !== nothing && s.econ.active
    econ_on && (s.econ.cycle_cost = 0f0; s.econ.cycle_rev = 0f0)
    rem = cuts!(s; fint = fint)                             # CUTS — thin (accrues econ per cut tree)
    rem.tpa > 0f0 && compute_density!(s)                    # recompute post-thin density
    if econ_on
        yr = current_cycle_year(s)   # IY schedule (TIMEINT/CYCLEAT-aware)
        s.econ.base_year < 0 && (s.econ.base_year = Int32(yr))
        (s.econ.cycle_cost > 0f0 || s.econ.cycle_rev > 0f0) &&
            push!(s.econ.harvests, (Float32(yr), s.econ.cycle_cost, s.econ.cycle_rev))
    end
    # FFE: a scheduled SIMFIRE burns this cycle (FMMAIN). The fire kill is periodic
    # MORTALITY (booked at the cycle-start volume, like OMORT) — added to `mort` below.
    fire_mort = _maybe_burn!(s, fint)
    # FFE: age the standing snag cohorts (falldown + decay) over the cycle.
    if s.fire !== nothing && s.fire.active && !isempty(s.fire.snags.sp)
        update_snags!(s, round(Int, fint))
    end
    apply_volume_overrides!(s; fint = fint)  # VOLUME/BFVOLUME merch-standard overrides (volkey.f)
    t = s.trees
    nlive = t.n                              # ORIGINAL live records (pre-tripling)
    # Cycle-start volume + TPA of the originals, for the period accounting.
    old_cfv = Float32[t.cuft_vol[i] for i in 1:nlive]
    old_tpa = Float32[t.tpa[i]      for i in 1:nlive]
    # Tripling is active only for the first ICL4 cycles (s.control.icl4; default 2, set to 0
    # by NOTRIPLE / to n by NUMTRIP); afterwards growth is the stochastic serial-correlation path.
    trip = !compressed && Int(s.control.cycle) < Int(s.control.icl4)   # COMPRESS suppresses tripling (NOTRIP)
    stash = diameter_growth!(s, s.variant; tripling = trip, sfint = fint)  # DGs only; no records yet
    height_growth!(s, s.variant; scale = fint / 5f0)         # HTG scaled to the cycle length
    small_tree_growth!(s, stash; fint = fint)  # REGENT overrides DG/HTG for dbh < 3"
    apply_fix_scalers!(s, stash, :fixdg, fint)   # FIXDG/FIXHTG: one-shot DG/HTG scalers,
    apply_fix_scalers!(s, stash, :fixhtg, fint)  # after all growth, before MORTS (grincr.f:451)
    mortality!(s, s.variant; fint = fint)  # MORTS on the ORIGINAL records (FVS order)
    g = s.plot.gross_space
    # Mortality volume (OMORT): the fire kill (booked above) plus this cycle's MORTS deaths,
    # accounted on the originals before tripling.
    mort = fire_mort
    @inbounds for i in 1:nlive
        mort += (old_tpa[i] - t.tpa[i]) * old_cfv[i]
    end
    triple_records!(s, stash)              # TRIPLE after mortality (splits surviving TPA)
    fertilizer_growth!(s; fint = fint)     # FFERT fertilizer DG/HTG boost (grincr.f:564, after TRIPLE)
    htgstp!(s; fint = fint)                # HTGSTOP/TOPKILL top damage (gradd.f:158, before UPDATE)
    # Per-record cycle-start CFV (tripled records inherit the originals' cycle-0 vol).
    n = t.n
    old_cfv2 = Float32[t.cuft_vol[i] for i in 1:n]
    sd = s.coef.species
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    @inbounds for i in 1:n
        # DG is the INSIDE-bark increment; outside-bark DBH grows by DG/bark, with
        # bark evaluated at the pre-growth DBH (update.f:115 / update.jl:75).
        bark = bark_ratio(bark_a, bark_b, t.species[i], t.dbh[i])
        t.dbh[i]    += t.diam_growth[i] / bark
        t.height[i] += t.ht_growth[i]
        # Broken-top trees: the full (NORMHT) height grows by the same increment as the
        # standing height (update.f:33-35), so the topkill cubic volume keeps growing.
        t.norm_ht[i] > 0 &&
            (t.norm_ht[i] = floor(Int32, Float32(t.norm_ht[i]) + t.ht_growth[i] * 100f0 + 0.5f0))
    end
    compute_volumes!(s)                     # end-of-period volumes
    accr = 0f0
    @inbounds for i in 1:n
        d = t.cuft_vol[i] - old_cfv2[i]     # OACC over the tripled set; FVS clamps
        d > 0f0 && (accr += d * t.tpa[i])   # negative growth to 0 (vols.f: CFV>tcf ⇒ WK5=0)
    end
    comcup!(t)                              # COMCUP (grincr.f:318, end of GRINCR): drop
                                            # PROB≤1e-5 records before GRADD/next cycle
    # GRADD order (gradd.f): UPDATE → DENSE → ESNUTR → DENSE → CROWN → VOLS. Establish
    # scheduled regen AFTER growth+mortality (fresh, full TPA this period) but BEFORE
    # CROWN, so the new trees' crown ratio (ICR) is computed this cycle (not carried
    # bogus into next cycle's DGF/mortality).
    esuckr!(s; fint = fint)                 # ESNUTR — stump/root sprouts (LSPRUT; before ESTAB)
    establish!(s; fint = fint)              # ESNUTR — adds regen (ICR=0), recomputes density
    crown_ratio_update!(s; fint = fint)     # CROWN — crown ratio for ALL trees incl. new regen
    # NOTE: newly-established trees get NO volume in their birth cycle. The oracle's
    # VOLS in the establishment cycle runs before the records are inserted, so a planted
    # stand reports CFV=0 at cyc1 (verified: bare_plant 1997 cuft=0) and the regen first
    # gets volume from the next cycle's VOLS (next grow_cycle!'s compute_volumes!). The
    # crown pass above DOES set the new trees' ICR this cycle (DGF/mortality read it next).
    s.control.cycle += Int32(1)
    return (; accretion = accr / fint / g, mortality = mort / fint / g)
end

"""
    run_keyfile(keypath; variant=Southern(), faithful=true, period=5) -> String

Full multi-stand run: project EVERY stand in `keypath` (the FVS `main.f` stand loop)
and return the concatenated `.sum` text — one `-999` header + per-cycle rows per
stand. Each stand is set up (`notre!`/`setup_growth!`/`compute_volumes!`) and projected
by `write_sum_file`, which runs the per-cycle loop including scheduled management
(CUTS/ESTAB/fire). Stands are independent (`each_stand` gives each a fresh state with
the tree format carried across), so this is also the unit of thread-parallelism.
"""
function run_keyfile(keypath::AbstractString; variant::AbstractVariant = Southern(),
                     faithful::Bool = true, period::Integer = 5,
                     date::AbstractString = "", time::AbstractString = "")
    out = IOBuffer()
    case = 0
    for s in each_stand(keypath; variant = variant, faithful = faithful)
        notre!(s)
        setup_growth!(s)
        compute_volumes!(s)
        sid = strip(s.plot.stand_id)
        mid = strip(s.plot.mgmt_id); mid = isempty(mid) ? "NONE" : String(mid)
        # DBS output (DATABASE block): collect this stand's summary rows and/or per-cycle tree
        # snapshots and append them to the DSNOUT database, in addition to the text `.sum`.
        has_db = !isempty(s.control.dbs_out_file)
        sum_on = s.control.dbs_summary && has_db
        tl_on = s.control.dbs_treelist && has_db
        cp_on = s.control.dbs_compute && has_db && !isempty(s.control.compute_defs)
        cl_on = s.control.dbs_cutlist && has_db
        rows = sum_on ? SummaryRow[] : nothing
        tl_cycles = tl_on ? Tuple[] : nothing
        cp_rows = cp_on ? Tuple[] : nothing
        cl_cycles = cl_on ? Tuple[] : nothing
        # FFE Stand Carbon Report (CARBREPT): collect a row per cycle from this same simulation.
        carb_rows = (s.control.carbon_report_on && s.fire !== nothing && s.fire.active) ? Tuple[] : nothing
        hook = tl_on ? (st, yr, pl) -> push!(tl_cycles, treelist_snapshot(st, yr, pl)) : nothing
        write_sum_file(out, s; period = Int(period), stand_id = String(sid),
                       mgmt_id = mid, date = date, time = time,
                       collect_rows = rows, cycle_hook = hook, compute_collect = cp_rows,
                       cutlist_collect = cl_cycles, carbon_collect = carb_rows)
        carb_rows === nothing ||
            write_carbon_report_block(out, carb_rows; stand_id = String(sid), mgmt_id = mid)
        if has_db
            case += 1
            caseid = string(sid, "-", case)
            # FVS_Cases registry + FVS_InvReference reference dump accompany any DBS output.
            kwfile = isempty(keypath) ? "" : first(splitext(basename(keypath)))
            write_dbs_cases!(s.control.dbs_out_file, caseid, String(sid);
                             mgmt_id = mid, variant = variant_code(s.variant),
                             keyword_file = kwfile, sampling_wt = s.plot.sample_weight,
                             run_datetime = strip(string(date, " ", time)))
            write_dbs_invref!(s.control.dbs_out_file, caseid, String(sid), s)
            sum_on && write_dbs_summary!(s.control.dbs_out_file, caseid, String(sid), rows;
                                         mgmt_id = mid, variant = variant_code(s.variant))
            tl_on && write_dbs_treelist!(s.control.dbs_out_file, caseid, String(sid), tl_cycles)
            cl_on && write_dbs_cutlist!(s.control.dbs_out_file, caseid, String(sid), cl_cycles)
            if carb_rows !== nothing
                write_dbs_carbon!(s.control.dbs_out_file, caseid, String(sid), carb_rows)
                write_dbs_fuels!(s.control.dbs_out_file, caseid, String(sid), carb_rows)
                write_dbs_snagsum!(s.control.dbs_out_file, caseid, String(sid), carb_rows)
                write_dbs_dwd_vol!(s.control.dbs_out_file, caseid, String(sid), carb_rows)
                write_dbs_dwd_cov!(s.control.dbs_out_file, caseid, String(sid), carb_rows)
            end
            # Fire-EVENT DBS tables: one row per SIMFIRE event (captured by fmburn!), independent of CARBREPT
            if s.fire !== nothing && s.fire.active && !isempty(s.fire.burn_reports)
                br = s.fire.burn_reports
                write_dbs_burnreport!(s.control.dbs_out_file, caseid, String(sid), br)
                write_dbs_mortality!(s.control.dbs_out_file, caseid, String(sid), br)
                write_dbs_consumption!(s.control.dbs_out_file, caseid, String(sid), br)
            end
            if cp_on
                var_names = String[nm for (_, nm, _) in s.control.compute_defs]
                write_dbs_compute!(s.control.dbs_out_file, caseid, String(sid), var_names, cp_rows)
            end
        end
    end
    return String(take!(out))
end
