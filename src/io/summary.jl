# =============================================================================
# summary.jl — the `.sum` stand-summary table writer (SUMOUT)
#
# Ported from: base/sumout.jl (SUMOUT). Writes the machine-readable summary file:
# a `-999` header line per stand, then one fixed-format row per simulation period.
#
# Column order of the 29-field row (sumout.f / IOSUM, with the before/after-
# treatment split):
#   year age TPA  BA SDI CCF TopHt QMD                         (start of period)
#   Tcuft Mcuft Scuft Bdft                                     (start volumes)
#   remTPA remTcuft remMcuft remScuft remBdft                  (removals)
#   atBA atSDI atCCF atTopHt atQMD                             (after treatment)
#   period accretion mortality MAI  fortype sizecls stockcls   (growth + class)
# =============================================================================

const _SUM_ROW_FMT = Printf.Format(
    "%4d%4d%6d%4d%5d%4d%4d%5.1f" *
    "%6d%6d%6d%6d%6d%6d%6d%6d%6d" *
    "%4d%5d%4d%4d%5.1f  %6d%5d%6d  %6.1f %3d %1d%1d\n")

"""
    SummaryRow

One `.sum` period row. Mirrors the SUMOUT/IOSUM columns; the `at_*`
(after-treatment) fields equal the start values when there is no cut in the
period, and the removal fields are 0.
"""
Base.@kwdef mutable struct SummaryRow
    year::Int; age::Int; tpa::Int
    ba::Int; sdi::Int; ccf::Int; topht::Int; qmd::Float64
    cuft::Int; mcuft::Int; scuft::Int; bdft::Int
    rem_tpa::Int = 0; rem_cuft::Int = 0; rem_mcuft::Int = 0; rem_scuft::Int = 0; rem_bdft::Int = 0
    at_ba::Int; at_sdi::Int; at_ccf::Int; at_topht::Int; at_qmd::Float64
    period::Int = 0; accretion::Int = 0; mortality::Int = 0; mai::Float64 = 0.0
    fortype::Int; sizecls::Int; stockcls::Int
end

"Write one SUMOUT-format period row."
function write_sum_row(io::IO, r::SummaryRow)
    Printf.format(io, _SUM_ROW_FMT,
        r.year, r.age, r.tpa, r.ba, r.sdi, r.ccf, r.topht, r.qmd,
        r.cuft, r.mcuft, r.scuft, r.bdft,
        r.rem_tpa, r.rem_cuft, r.rem_mcuft, r.rem_scuft, r.rem_bdft,
        r.at_ba, r.at_sdi, r.at_ccf, r.at_topht, r.at_qmd,
        r.period, r.accretion, r.mortality, r.mai,
        r.fortype, r.sizecls, r.stockcls)
    return io
end

# =============================================================================
# CSV form of the .sum — the modern, named-column, machine-readable summary (the
# output analog of the input .tre→.csv / .key→.yaml modernization). Same data as the
# fixed-column .sum (columns documented in docs/FORMATS.md §4), flattened across stands
# with a leading StandID/MgmtID so it loads straight into pandas/R/a spreadsheet.
# =============================================================================
"Named-column header for the summary CSV (StandID/MgmtID/Title + the SUMOUT columns)."
const SUM_CSV_HEADER = [
    "StandID", "MgmtID", "Title", "Year", "Age", "Tpa", "BA", "SDI", "CCF", "TopHt", "QMD",
    "TCuFt", "MCuFt", "SCuFt", "BdFt", "RTpa", "RTCuFt", "RMCuFt", "RSCuFt", "RBdFt",
    "ATBA", "ATSDI", "ATCCF", "ATTopHt", "ATQMD", "PrdLen", "Accret", "Mort", "MAI",
    "ForType", "SizeCls", "StkCls"]

_sumf1(x) = @sprintf("%.1f", Float64(x))   # the float columns (QMD/ATQMD/MAI) print like the .sum
# Minimal CSV quoting: wrap in double-quotes (and double any embedded quote) iff the field
# carries a comma / quote / newline. The Title (the STDIDENT description) can contain commas.
_csvq(s) = (t = String(s); occursin(r"[\",\n]", t) ? '"' * replace(t, '"' => "\"\"") * '"' : t)

"Render one `SummaryRow` as a CSV row (SUM_CSV_HEADER order) for stand `sid`/`mid`/`title`."
function sum_csv_row(sid::AbstractString, mid::AbstractString, title::AbstractString, r::SummaryRow)
    join((_csvq(sid), _csvq(mid), _csvq(title),
          r.year, r.age, r.tpa, r.ba, r.sdi, r.ccf, r.topht, _sumf1(r.qmd),
          r.cuft, r.mcuft, r.scuft, r.bdft, r.rem_tpa, r.rem_cuft, r.rem_mcuft, r.rem_scuft, r.rem_bdft,
          r.at_ba, r.at_sdi, r.at_ccf, r.at_topht, _sumf1(r.at_qmd),
          r.period, r.accretion, r.mortality, _sumf1(r.mai), r.fortype, r.sizecls, r.stockcls), ',')
end

"""
    write_sum_csv(io, stands)

Write the summary CSV: the header then, for each `(stand_id, mgmt_id, title, rows)` in `stands`,
one line per `SummaryRow`. `stands` is a vector of such 4-tuples (one per projected stand). The
`Title` (the STDIDENT description that the fixed-column .sum drops) is carried as its own column.
"""
function write_sum_csv(io::IO, stands)
    println(io, join(SUM_CSV_HEADER, ','))
    for (sid, mid, title, rows) in stands
        for r in rows
            println(io, sum_csv_row(String(sid), String(mid), String(title), r))
        end
    end
    return io
end

"""
    write_sum_header(io, nperiods, stand_id, mgmt_id, sample_wt, variant, date, time, nplots)

Write the per-stand `-999` header line that precedes the period rows.
"""
# Fortran E15.7 edit descriptor: value = 0.DDDDDDD·10^p (0.1≤|m|<1), 7 mantissa digits, 2-digit exponent,
# right-justified in width 15. FVS's .sum header sample-weight uses this — NOT C %15.7E (which prints
# D.DDDDDDD·10^(p-1), 8 sig digits): 11.0 → Fortran "0.1100000E+02" vs C "1.1000000E+01". (io-serialization #2)
function _fortran_e15_7(x::Real)
    x = Float64(x)
    x == 0.0 && return "  0.0000000E+00"
    neg = x < 0; a = abs(x)
    p = floor(Int, log10(a)) + 1
    mi = round(Int, a / 10.0^p * 1.0e7)
    mi >= 10_000_000 && (mi = 1_000_000; p += 1)
    lpad(@sprintf("%s0.%07dE%s%02d", neg ? "-" : "", mi, p < 0 ? "-" : "+", abs(p)), 15)
end

function write_sum_header(io::IO, nperiods::Integer, stand_id::AbstractString,
                          mgmt_id::AbstractString, sample_wt::Real, variant::AbstractString,
                          date::AbstractString, time::AbstractString, nplots::Integer)
    @printf(io, "-999%5d %-26s %-4s%s %-2s %-10s %-8s %-10s %-11s%3d\n",
        nperiods, stand_id, mgmt_id, _fortran_e15_7(Float32(sample_wt)), variant, date, time,
        "          ", "           ", nplots)
    return io
end

"""
    write_sum_file(io, state; period=5, stand_id="", mgmt_id="NONE",
                   sample_wt=nothing, variant="SN", date="", time="", header=true)

Run the full projection and write the complete `.sum` table: the `-999` header
(when `header`) plus one period row per cycle. Each row's start-of-period stats
come from `summary_row` at that cycle; the growth columns (accretion/mortality)
come from `grow_cycle!` advancing to the next period. The final cycle has no
growth period. Cumulative removed merch volume feeds MAI. Requires `state` set up
through `setup_growth!` + `compute_forest_type!` + `compute_volumes!`.
"""
function write_sum_file(io::IO, s::StandState; period::Int = 5,
                        stand_id::AbstractString = "", mgmt_id::AbstractString = "NONE",
                        sample_wt = nothing, variant::AbstractString = "SN",
                        date::AbstractString = "", time::AbstractString = "", header::Bool = true,
                        collect_rows::Union{Nothing,Vector} = nothing, cycle_hook = nothing,
                        compute_collect::Union{Nothing,Vector} = nothing,
                        cutlist_collect::Union{Nothing,Vector} = nothing,
                        carbon_collect::Union{Nothing,Vector} = nothing,
                        potfire_collect::Union{Nothing,Vector} = nothing,
                        hrvcarbon_collect::Union{Nothing,Vector} = nothing)
    build_cycle_schedule!(s)                 # ensure the IY boundary-year array is current (idempotent)
    ncyc = Int(s.control.ncycle_eff)         # rows = ncyc + 1 (inventory + each cycle, post-CYCLEAT)
    ncyc < 1 && (ncyc = Int(s.control.ncycle))
    # FFE Stand Carbon Report (CARBREPT): carbon_on gates only the REPORT-row collection. The per-cycle
    # fuel DYNAMICS (below) run for any FFE-active stand.
    carbon_on = carbon_collect !== nothing && s.fire !== nothing && s.fire.active
    # FFE surface-fuel dynamics (decay + snag falldown + crown-lift → down wood) run EVERY cycle for any
    # FFE-active stand, exactly as FVS fmmain.f does — NOT only when the Carbon report is on. The down-
    # wood pools (FireState.cwd) feed FMCFMD's (SMALL,LARGE) fuel-model selection, so a SIMFIRE stand's
    # fire behavior depends on them having evolved since inventory. The prior carbon-only gate froze cwd
    # at the inventory value (fire_fuel9 2005 sm=7.02 == 1990, vs FVS's accumulated 9.19), which selected
    # fuel model 5 over 12 and gave byram 2905 vs FVS's 4194 (~6× low on the FM5 component).
    ffe_on = s.fire !== nothing && s.fire.active
    if ffe_on
        ffe_seed_input_snags!(s)             # inventory snags from the input dead records (FMSADD ITYP=3)
        fill!(s.fire.crown_lift_annual, 0f0)
        snapshot_ffe_oldcrown!(s)            # FMOLDC at inventory: gives the 1st cycle's crown-lift a valid
                                             # OLDCRW (else the 1st cycle's fine down-wood is lost; DDW gap)
    end
    # SNAGINIT (act 2522, FMSADD) is a SCHEDULED activity that FVS runs in FMMAIN DURING cycle 1, AFTER the
    # inventory carbon report — NOT at inventory like the input dead-tree snags. So the user snags first appear
    # in cycle-1's report (live net01 SnagDet is empty at the inventory year; the SNAGINIT cohort shows up the
    # next cycle). Defer the add to the start of the first growing cycle.
    snaginit_pending = ffe_on
    g = s.plot.gross_space
    # .sum header sample weight = SAMWT (s.plot.sample_weight), NOT gross_space. (gross_space
    # is the non-stockable expansion used for per-acre normalization below; SAMWT is the
    # stand's sampling weight, which FVS prints in the -999 header — e.g. 11, not 1.1.)
    sw = sample_wt === nothing ? s.plot.sample_weight : sample_wt
    if header
        write_sum_header(io, ncyc + 1, stand_id, mgmt_id, sw, variant, date, time, Int(s.plot.pi))
    end
    cum_rem_merch = 0f0
    prev_increment = 0f0   # removed-merch added in the most recent growing cycle (for the MAI final-row quirk)
    di(x) = trunc(Int, x + 0.5)
    for c in 0:ncyc
        compute_forest_type!(s)
        last = c == ncyc
        per = last ? 0 : cycle_period_at(s.control, c)   # THIS cycle's length (varies w/ TIMEINT/CYCLEAT)
        # main columns reflect the start-of-cycle (pre-thin) stand.
        # MAI terminal-row quirk (evtstv.f:414 + disply.f:392): intermediate rows accumulate
        # removed merch with a one-cycle lag, but the FINAL row's MAI is loaded from the
        # un-incremented TOTREM — i.e. it EXCLUDES the last growing cycle's removal.
        r = summary_row(s; period = per,
                        total_removed_merch = cum_rem_merch - (last ? prev_increment : 0f0))
        # per-cycle hook (DBS TreeList): the start-of-cycle (pre-thin) tree list at year r.year
        cycle_hook === nothing || cycle_hook(s, r.year, per)
        # FFE Stand Carbon Report row (FMCRBOUT, fmmain.f:206) — sampled at the FVS phase: AFTER FMBURN
        # (fire kill + snag booking + consumption) but BEFORE UPDATE grows the stand. For a non-fire cycle
        # that phase equals the cycle-top, pre-growth stand (sampled here). For a SIMFIRE cycle the row
        # must be POST-fire, so defer it to grow_cycle!'s `carbon_hook` (fired right after the fire, on the
        # post-fire cycle-start-size stand) — else the fire's snag/AGL/Released effects surface one row late.
        # Carbon-Released-from-Fire: 0 unless a SIMFIRE burned in r.year (fmburn! records it in burn_reports);
        # convert tons-C/ac → the report units (same factor as stand_carbon_report's pools, carbon.jl).
        _carb_push(st) = begin
            rel = 0f0
            if st.fire !== nothing
                @inbounds for br in st.fire.burn_reports
                    br.year == Int(r.year) && (rel = br.released)
                end
            end
            uf = st.control.carbon_units == 1 ? 0.90718474f0 / 0.40468564f0 :
                 st.control.carbon_units == 2 ? 0.90718474f0 : 1f0
            push!(carbon_collect, (r.year, stand_carbon_report(st), ffe_fuel_loadings(st),
                                   snag_summary(st), ffe_down_wood(st), rel * uf))
        end
        # A SIMFIRE cycle: the fire (inside grow_cycle!'s mortality_and_fire!) must consume + snag the
        # START-of-cycle fuels, so this cycle's pre-grow ffe_fuel_update! is WITHHELD and its period handed
        # to grow_cycle! (run post-fire, FVS FMBURN→annual-loop order). The carbon row is likewise deferred
        # to the post-fire `carbon_hook`. (`fire_cycle` adds the carbon-report gate.)
        fire_this_cycle = !last && _fire_due(s) && per > 1   # OPCYCL: cycle range contains fire_year
        fire_cycle = carbon_on && fire_this_cycle
        if carbon_on && !fire_cycle
            compute_density!(s); fmcba!(s)                    # refresh cover type + live fuels (FLIVE)
            _carb_push(s)
        end
        # FVS_PotFire: the potential-fire behavior under fixed severe/moderate weather (FMPOFL), per cycle
        if potfire_collect !== nothing && s.fire !== nothing && s.fire.active
            compute_density!(s)
            pfr = potential_fire_report(s)
            pfr !== nothing && push!(potfire_collect, (r.year, pfr))
        end
        if !last
            # Add the deferred SNAGINIT snags at the start of the first growing cycle (FMMAIN, post-inventory-
            # report). They are then present for this cycle's snag falldown + any SIMFIRE, and first surface in
            # the NEXT cycle's carbon/snag report — matching live FVS.
            if snaginit_pending
                ffe_add_snaginit!(s); snaginit_pending = false
            end
            # DBS FVS_Compute: snapshot the active COMPUTE variables at this (growing) cycle's
            # start — only the growing cycles get a row (the event monitor runs during growth).
            compute_collect === nothing ||
                push!(compute_collect, (r.year, snapshot_compute!(s, r.year, c)))
            # apply this cycle's scheduled thin (CUTS) BEFORE growth; report the removed
            # + after-treatment columns on THIS row (matching the Fortran .sum). cuts! is
            # idempotent, so grow_cycle!'s own cuts! call below is then a no-op.
            # FVS_CutList: arm the per-record cut sink for this (real) thin, then stash + disarm.
            cutlist_collect === nothing || (s.control.cutlist_capture = Any[])
            rem = cuts!(s; fint = Float32(per))
            if cutlist_collect !== nothing
                push!(cutlist_collect, (r.year, per, s.control.cutlist_capture))
                s.control.cutlist_capture = nothing
            end
            if rem.tpa > 0f0
                compute_density!(s)
                r.rem_tpa  = di(rem.tpa / g);  r.rem_cuft  = di(rem.cuft / g)
                r.rem_mcuft = di(rem.mcuft / g); r.rem_scuft = di(rem.scuft / g)
                r.rem_bdft = di(rem.bdft / g)
                r.at_ba = di(stand_ba(s) / g);  r.at_sdi = di(stand_sdi(s) / g)
                r.at_ccf = di(stand_ccf(s) / g); r.at_topht = di(stand_top_height(s))
                r.at_qmd = round(stand_qmd(s); digits = 1)
                prev_increment = rem.mcuft / g
                cum_rem_merch += prev_increment
            else
                prev_increment = 0f0   # this growing cycle had no removal (final-row MAI subtracts 0)
            end
            # FVS_Hrv_Carbon: collect AFTER this cycle's cut so year r.year's harvest is booked (KYR=1),
            # matching FMCHRVOUT. Habitat default 1 (SC); SE (2) applies to specific national forest codes.
            hrvcarbon_collect === nothing || s.fire === nothing || !s.fire.active ||
                push!(hrvcarbon_collect, (r.year, harvested_carbon_report(s, r.year, 1)))
            # FFE annual fuel loop BEFORE growth (report→fuel→grow). FVS interleaves the annual fuel steps
            # with FMBURN, so a fire fires on the cycle-start + fire-year's single annual step — NOT the
            # period-end fuel. When a SIMFIRE burns this cycle, split the loop: advance 1 year, stash the
            # (SMALL,LARGE) the fire burns on, then advance the rest. Non-fire cycles run the full loop once.
            if ffe_on
                # FMMAIN runs FMBURN (the fire, fmmain.f:170) BEFORE the annual fuel loop (FMSNAG/FMCWD/
                # FMCADD, fmmain.f:228), so the fire samples the START-OF-CYCLE down wood. Stash (SMALL,LARGE)
                # for the fire basis here. For a fire cycle the annual loop is DEFERRED into grow_cycle!
                # (run post-fire); non-fire cycles advance the pools the full period now (report→fuel→grow).
                # A fire in the FIRST FFE cycle (s10_fire SIMFIRE@cycle1) burns before any prior ffe_fuel_update!
                # loaded the dead-fuel pools, so init them now (FVS FMCBA precedes the first FMBURN). Otherwise
                # the stash would read zero cwd ⇒ low-fuel model ⇒ under-kill (jl 119 vs live 57 TPA). For
                # cycle≥2 fires the pools are already loaded (fuels_init), so fire_carbon stays bit-exact.
                fire_this_cycle && !s.fire.fuels_init && (compute_density!(s); fmcba!(s))
                fire_this_cycle && (s.fire.fire_smlg = _small_large_fuel(s.fire))
                fire_this_cycle || ffe_fuel_update!(s, per)
            end
            chook = fire_cycle ? (st -> (compute_density!(st); fmcba!(st); _carb_push(st))) : nothing
            gr = grow_cycle!(s; fint = Float32(per), carbon_hook = chook,
                             fuel_period = fire_this_cycle ? per : nothing)   # advances cycle
            r.accretion = trunc(Int, gr.accretion + 0.5)
            r.mortality = trunc(Int, gr.mortality + 0.5)
            if ffe_on                                   # crown-lift from THIS growth (FMSDIT) + FMOLDC snapshot
                compute_crown_lift!(s, per); snapshot_ffe_oldcrown!(s)
            end
        elseif hrvcarbon_collect !== nothing && s.fire !== nothing && s.fire.active
            push!(hrvcarbon_collect, (r.year, harvested_carbon_report(s, r.year, 1)))  # final cycle (no cut block)
        end
        write_sum_row(io, r)
        collect_rows === nothing || push!(collect_rows, r)
    end
    return io
end

"""
    summary_row(state) -> SummaryRow

Build the start-of-period `.sum` row from the current (already-grown-to-cycle)
stand state: per-acre TPA/BA/SDI/CCF/top-height/QMD, the four stand volumes, and
the forest-type / size / stocking classes. Removal, after-treatment and growth
(accretion/mortality/MAI) fields are filled by the cycle driver. The integer
columns use FVS's truncate-after-+0.5 rounding (`_dtrunc`)."""
function summary_row(s::StandState; period::Int = 0, total_removed_merch::Real = 0,
                     accretion::Real = 0, mortality::Real = 0)
    g = s.plot.gross_space
    dt(x) = trunc(Int, x + 0.5f0)
    tpa  = dt(stand_tpa(s) / g)
    ba   = dt(stand_ba(s) / g)
    sdi  = dt(stand_sdi(s) / g)
    ccf  = dt(stand_ccf(s) / g)
    toph = dt(stand_top_height(s))
    qmd  = round(stand_qmd(s); digits = 1)
    t = s.trees
    # STRICTLY SEQUENTIAL Float32 accumulation (ACC += VOL[i]·PROB[i], i=1..n) to match FVS's DISPLY DO-loop
    # order — Julia's `sum(generator)` may use PAIRWISE reduction, which reorders the Float32 adds and flips
    # the rendered integer by 1 on knife-edge rows (the non-associative tree-SUM residual).
    function vtot(f)
        # `f` is a runtime Symbol ⇒ getfield(t, f) infers as Any, boxing `fld[i]` on every add (measured
        # ~50 KB/cycle + a type-instability). All vtot fields are Vector{Float32}, so assert it: concrete
        # `fld` ⇒ allocation-free, type-stable, and the sequential Float32 accumulation order is unchanged.
        fld = getfield(t, f)::Vector{Float32}; acc = 0f0
        @inbounds for i in 1:t.n
            acc += fld[i] * t.tpa[i]
        end
        return dt(acc / g)
    end
    # Year/age come from the cycle-boundary schedule (IY, build_cycle_schedule!): the calendar
    # year at this cycle's start, and the age advanced by the elapsed years from the inventory.
    # For uniform cycles this is exactly cycle_year[1] + cyc·per (bit-exact); non-uniform TIMEINT
    # and CYCLEAT-inserted boundaries make the steps vary.
    cyc = Int(s.control.cycle)
    yr = cycle_year_at(s.control, cyc)
    age = Int(s.plot.stand_age) + (yr - Int(s.control.cycle_year[1]))
    # RESETAGE (resage.f): after the reset year (run after DISPLY, so the reset row itself keeps
    # the old age) the stand age is rebased — age(Y) = age_reset_age + (Y − reset_year).
    ry = Int(s.control.age_reset_year)
    (ry >= 0 && yr > ry) && (age = Int(s.control.age_reset_age) + (yr - ry))
    mcuft = vtot(:merch_cuft_vol)
    # MAI (BCYMAI, disply.f:383): (merch cuft + cumulative removed merch) / age.
    # `total_removed_merch` carries cross-cycle removals (0 at the inventory).
    # Computed in Float32 to match FVS REAL*4 (the %.1f rounding differs from Float64).
    # After a RESETAGE that rebased the age to ZERO, FVS shuts off MAI (disply.f:391-394 BCYMAI=0 when
    # MAIFLG≠0; evtstv.f:396 sets it when ZERO=age−period==0, i.e. the age was reset to 0, and persists it).
    # A RESETAGE to a NON-zero age keeps MAI on (e.g. s17_managed resets to 40 → MAI stays 62.5). Non-RESETAGE
    # runs have ry<0 ⇒ untouched (bit-exact); bare-ground (NEWSTD) has no RESETAGE ⇒ its own MAI path unchanged.
    mai = (ry >= 0 && yr > ry && Int(s.control.age_reset_age) == 0) ? 0f0 :
          (age > 0 ? Float32(mcuft + total_removed_merch) / Float32(age) : 0f0)
    SummaryRow(
        year = yr, age = age, tpa = tpa,
        ba = ba, sdi = sdi, ccf = ccf, topht = toph, qmd = qmd,
        cuft = vtot(:cuft_vol), mcuft = mcuft,
        scuft = vtot(:saw_cuft_vol), bdft = vtot(:bdft_vol),
        at_ba = ba, at_sdi = sdi, at_ccf = ccf, at_topht = toph, at_qmd = qmd,
        period = period, mai = mai,
        accretion = trunc(Int, accretion + 0.5), mortality = trunc(Int, mortality + 0.5),
        fortype = Int(s.plot.forest_type), sizecls = Int(s.plot.size_class),
        stockcls = Int(s.plot.stocking_class))
end
