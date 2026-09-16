"""
    Aggregate

Aggregate an AOI's plots into distributions for the charts.

Each 30 m pixel is 900 m² = 0.222395 ac. A plot imputed to `n` pixels inside the
AOI represents `n × 0.222395` acres; a tree's per-acre expansion (`TPA_UNADJ`)
scaled by those acres is its expanded count in the AOI.

Cycle-0 species × DBH-class **count** and **basal area** are exact from the tree
list. AOI totals of biomass / volume / carbon come from the VAT (already FVS/FIA
computed per acre), area-weighted. Per-tree biomass/volume/carbon by species×size
(and forward projection) arrive with the FVSjl adapter.
"""
module Aggregate

using ..TreeList
using ..VAT

export aggregate_aoi, ACRES_PER_PIXEL, dbh_class_lo

const ACRES_PER_PIXEL = 900.0 / 4046.8564224          # 0.222395 ac
const DBH_CAP = 40                                     # top open class "40+"

ba_per_tree(dia) = 0.005454154 * dia * dia            # ft² for one tree at DBH=dia in

"Lower bound of the 2-inch DBH class holding `dia` (inches)."
dbh_class_lo(dia) = min(2 * floor(Int, dia / 2), DBH_CAP)

struct SpeciesDbhCell
    spcd::Int16
    dbh_lo::Int
    count::Float64      # expanded trees in AOI
    ba::Float64         # ft² total in AOI
end

struct AoiAggregate
    nplots::Int
    npixels::Int
    acres::Float64
    # per-(species,dbh) cells
    cells::Vector{SpeciesDbhCell}
    species_name::Dict{Int16,String}      # spcd -> common name
    species_symbol::Dict{Int16,String}
    # AOI totals (expanded)
    ba_total::Float64                     # ft² over AOI
    tpa_expanded::Float64                 # total live trees in AOI
    carbon_l::Float64                     # tons (live aboveground)
    carbon_d::Float64
    carbon_dwn::Float64
    drybio_l::Float64                     # tons aboveground dry live biomass
    volcfnet_l::Float64                   # ft³ live net volume
end

"""
    aggregate_aoi(counts, tt, st) -> AoiAggregate

`counts` :: Dict{tm_id => pixel_count} from `Raster.aoi_tally`.
"""
function aggregate_aoi(counts::AbstractDict, tt::TreeList.TreeTable, st::VAT.StandTable)
    acc = Dict{Tuple{Int16,Int},NTuple{2,Float64}}()   # (spcd,dbh_lo) -> (count,ba)
    names = Dict{Int16,String}()
    syms = Dict{Int16,String}()
    nplots = 0; npixels = 0
    acres_tot = 0.0
    ba_total = 0.0; tpa_exp = 0.0
    cL = 0.0; cD = 0.0; cDwn = 0.0; bioL = 0.0; volL = 0.0

    pairs = collect(counts)                       # (tm_id, pixels)
    # cheap serial pass: VAT-precomputed stand totals + area/counts (no I/O)
    for (tm, px) in pairs
        acres = px * ACRES_PER_PIXEL
        npixels += px; nplots += 1; acres_tot += acres
        rec = VAT.stand(st, tm)
        rec === nothing && continue
        isnan(rec.carbon_l)   || (cL   += rec.carbon_l   * acres)
        isnan(rec.carbon_d)   || (cD   += rec.carbon_d   * acres)
        isnan(rec.carbon_dwn) || (cDwn += rec.carbon_dwn * acres)
        isnan(rec.drybio_l)   || (bioL += rec.drybio_l   * acres)
        isnan(rec.volcfnet_l) || (volL += rec.volcfnet_l * acres)
    end

    # parallel pass: the per-plot tree-list reads (CSV seek+parse) fanned across tasks,
    # each with a local accumulator, then merged — the species×DBH tally + BA/TPA totals.
    nchunks = clamp(Threads.nthreads(), 1, max(1, length(pairs)))
    chunks = [pairs[i:nchunks:end] for i in 1:nchunks]
    tasks = map(chunks) do chunk
        Threads.@spawn begin
            la = Dict{Tuple{Int16,Int},NTuple{2,Float64}}()
            lba = 0.0; ltpa = 0.0
            ln = Dict{Int16,String}(); ls = Dict{Int16,String}()
            for (tm, px) in chunk
                acres = px * ACRES_PER_PIXEL
                for tr in TreeList.treelist(tt, tm)
                    tr.statuscd == 1 || continue
                    (isnan(tr.dia) || tr.dia <= 0) && continue
                    w = tr.tpa_unadj * acres
                    ba = ba_per_tree(tr.dia) * w
                    k = (tr.spcd, dbh_class_lo(tr.dia))
                    c, b = get(la, k, (0.0, 0.0)); la[k] = (c + w, b + ba)
                    lba += ba; ltpa += w
                    get!(ln, tr.spcd, tr.common_name); get!(ls, tr.spcd, tr.species_symbol)
                end
            end
            (la, lba, ltpa, ln, ls)
        end
    end
    for t in tasks
        la, lba, ltpa, ln, ls = fetch(t)
        for (k, v) in la; c, b = get(acc, k, (0.0, 0.0)); acc[k] = (c + v[1], b + v[2]); end
        ba_total += lba; tpa_exp += ltpa
        merge!(names, ln); merge!(syms, ls)
    end

    cells = [SpeciesDbhCell(k[1], k[2], v[1], v[2]) for (k, v) in acc]
    sort!(cells; by = c -> (-(sum_species_ba(acc, c.spcd)), c.spcd, c.dbh_lo))
    AoiAggregate(nplots, npixels, acres_tot, cells, names, syms,
                 ba_total, tpa_exp, cL, cD, cDwn, bioL, volL)
end

sum_species_ba(acc, spcd) = sum(v[2] for (k, v) in acc if k[1] == spcd; init = 0.0)

end # module
