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

"""
    write_sum_header(io, nperiods, stand_id, mgmt_id, sample_wt, variant, date, time, nplots)

Write the per-stand `-999` header line that precedes the period rows.
"""
function write_sum_header(io::IO, nperiods::Integer, stand_id::AbstractString,
                          mgmt_id::AbstractString, sample_wt::Real, variant::AbstractString,
                          date::AbstractString, time::AbstractString, nplots::Integer)
    @printf(io, "-999%5d %-26s %-4s%15.7E %-2s %-10s %-8s %-10s %-11s%3d\n",
        nperiods, stand_id, mgmt_id, Float32(sample_wt), variant, date, time,
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
                        carbon_collect::Union{Nothing,Vector} = nothing)
    build_cycle_schedule!(s)                 # ensure the IY boundary-year array is current (idempotent)
    ncyc = Int(s.control.ncycle_eff)         # rows = ncyc + 1 (inventory + each cycle, post-CYCLEAT)
    ncyc < 1 && (ncyc = Int(s.control.ncycle))
    # FFE Stand Carbon Report (CARBREPT): drive the per-cycle fuel dynamics on THIS same simulation and
    # collect a report row each cycle. Gated on an active FireState so non-FFE / fire-only stands are
    # untouched (the carbon report's fuel accumulation must not perturb a SIMFIRE stand's fuels).
    carbon_on = carbon_collect !== nothing && s.fire !== nothing && s.fire.active
    if carbon_on
        ffe_seed_input_snags!(s)             # inventory snags from the input dead records (FMSADD ITYP=3)
        fill!(s.fire.crown_lift_annual, 0f0)
    end
    g = s.plot.gross_space
    sw = sample_wt === nothing ? g : sample_wt
    if header
        write_sum_header(io, ncyc + 1, stand_id, mgmt_id, sw, variant, date, time, Int(s.plot.pi))
    end
    cum_rem_merch = 0f0
    di(x) = trunc(Int, x + 0.5)
    for c in 0:ncyc
        compute_forest_type!(s)
        last = c == ncyc
        per = last ? 0 : cycle_period_at(s.control, c)   # THIS cycle's length (varies w/ TIMEINT/CYCLEAT)
        # main columns reflect the start-of-cycle (pre-thin) stand
        r = summary_row(s; period = per, total_removed_merch = cum_rem_merch)
        # per-cycle hook (DBS TreeList): the start-of-cycle (pre-thin) tree list at year r.year
        cycle_hook === nothing || cycle_hook(s, r.year, per)
        # FFE Stand Carbon Report row (FMCRBOUT) — emitted BEFORE this cycle's fuel loop + growth, on the
        # post-growth/pre-fuel stand (matching FVS: FMCRBOUT runs before the annual FMSNAG/FMCWD/FMCADD loop).
        if carbon_on
            compute_density!(s); fmcba!(s)                    # refresh cover type + live fuels (FLIVE)
            push!(carbon_collect, (r.year, stand_carbon_report(s), ffe_fuel_loadings(s), snag_summary(s)))
        end
        if !last
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
                cum_rem_merch += rem.mcuft / g
            end
            carbon_on && ffe_fuel_update!(s, per)      # FFE annual fuel loop BEFORE growth (report→fuel→grow)
            gr = grow_cycle!(s; fint = Float32(per))   # advances cycle, returns period accr/mort
            r.accretion = trunc(Int, gr.accretion + 0.5)
            r.mortality = trunc(Int, gr.mortality + 0.5)
            if carbon_on                                # crown-lift from THIS growth (FMSDIT) + FMOLDC snapshot
                compute_crown_lift!(s, per); snapshot_ffe_oldcrown!(s)
            end
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
    vtot(f) = dt(sum((getfield(t, f)[i] * t.tpa[i] for i in 1:t.n); init = 0f0) / g)  # init for bare/empty stands
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
    mai = age > 0 ? Float32(mcuft + total_removed_merch) / Float32(age) : 0f0
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
