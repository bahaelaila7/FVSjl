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
Base.@kwdef struct SummaryRow
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
    summary_row(state) -> SummaryRow

Build the start-of-period `.sum` row from the current (already-grown-to-cycle)
stand state: per-acre TPA/BA/SDI/CCF/top-height/QMD, the four stand volumes, and
the forest-type / size / stocking classes. Removal, after-treatment and growth
(accretion/mortality/MAI) fields are filled by the cycle driver. The integer
columns use FVS's truncate-after-+0.5 rounding (`_dtrunc`)."""
function summary_row(s::StandState; period::Int = 0)
    g = s.plot.gross_space
    dt(x) = trunc(Int, x + 0.5f0)
    tpa  = dt(stand_tpa(s) / g)
    ba   = dt(stand_ba(s) / g)
    sdi  = dt(stand_sdi(s) / g)
    ccf  = dt(stand_ccf(s) / g)
    toph = dt(stand_top_height(s))
    qmd  = round(stand_qmd(s); digits = 1)
    t = s.trees
    vtot(f) = dt(sum(getfield(t, f)[i] * t.tpa[i] for i in 1:t.n) / g)
    ci = Int(s.control.cycle) + 1                       # year at start of this cycle (IY)
    yr = ci <= length(s.control.cycle_year) ? Int(s.control.cycle_year[ci]) : 0
    SummaryRow(
        year = yr, age = Int(s.plot.stand_age), tpa = tpa,
        ba = ba, sdi = sdi, ccf = ccf, topht = toph, qmd = qmd,
        cuft = vtot(:cuft_vol), mcuft = vtot(:merch_cuft_vol),
        scuft = vtot(:saw_cuft_vol), bdft = vtot(:bdft_vol),
        at_ba = ba, at_sdi = sdi, at_ccf = ccf, at_topht = toph, at_qmd = qmd,
        period = period,
        fortype = Int(s.plot.forest_type), sizecls = Int(s.plot.size_class),
        stockcls = Int(s.plot.stocking_class))
end
