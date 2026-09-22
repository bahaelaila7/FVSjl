"""
    VAT

Loader for the TreeMap 2022 CONUS raster value attribute table
(`TreeMap2022_CONUS.tif.vat.dbf`, 65,043 plots × 27 fields).

The VAT is both the tm_id → PLT_CN crosswalk AND a table of precomputed
cycle-0 stand state (BALIVE, QMD, TPA, volume, biomass, carbon, canopy,
stand height, forest type). So the "current conditions" map layer and its
aggregates need no FVS run at all — only the forward projection does.

Keyed by TM_ID (== the raster pixel value).
"""
module VAT

using DBFTables

export StandRecord, StandTable, load_vat, stand

"Precomputed cycle-0 stand state for one TreeMap plot (from the raster VAT)."
struct StandRecord
    tm_id::Int32
    plt_cn::Int64
    count::Int32          # pixels imputed to this plot across CONUS
    fortypcd::Int16
    fortyp_name::String
    balive::Float64       # ft2/acre live basal area
    qmd::Float64          # in quadratic mean diameter
    sdisum::Float64       # stand density index
    tpa_live::Float64     # trees/acre
    tpa_dead::Float64
    volcfnet_l::Float64   # ft3/acre live net volume
    volcfnet_d::Float64
    volbfnet_l::Float64   # bd ft/acre live sawlog
    drybio_l::Float64     # tons/acre aboveground dry live biomass
    drybio_d::Float64
    carbon_l::Float64     # tons/acre live aboveground carbon
    carbon_d::Float64     # tons/acre standing dead carbon
    carbon_dwn::Float64   # tons/acre down dead carbon
    standht::Float64      # ft height of dominant trees
    canopypct::Float64    # percent live canopy cover
end

struct StandTable
    by_tm::Dict{Int32,StandRecord}
end

stand(st::StandTable, tm_id::Integer) = get(st.by_tm, Int32(tm_id), nothing)

_f(x) = x === missing ? NaN : Float64(x)
_i(x) = x === missing ? Int32(-1) : round(Int32, Float64(x))

"""
    load_vat(dbf_path) -> StandTable

Read the whole VAT into a `TM_ID => StandRecord` dictionary.
"""
function load_vat(dbf_path::AbstractString)::StandTable
    tbl = DBFTables.Table(dbf_path)
    by_tm = Dict{Int32,StandRecord}()
    sizehint!(by_tm, 70000)
    for row in tbl
        tm = _i(row.TM_ID)
        rec = StandRecord(
            tm,
            round(Int64, _f(row.PLT_CN)),
            _i(row.Count),
            _i(row.FORTYPCD),
            row.ForTypName === missing ? "" : String(strip(row.ForTypName)),
            _f(row.BALIVE),
            _f(row.QMD),
            _f(row.SDIsum),
            _f(row.TPA_LIVE),
            _f(row.TPA_DEAD),
            _f(row.VOLCFNET_L),
            _f(row.VOLCFNET_D),
            _f(row.VOLBFNET_L),
            _f(row.DRYBIO_L),
            _f(row.DRYBIO_D),
            _f(row.CARBON_L),
            _f(row.CARBON_D),
            _f(row.CARBON_DWN),
            _f(row.STANDHT),
            _f(row.CANOPYPCT),
        )
        by_tm[tm] = rec
    end
    StandTable(by_tm)
end

end # module
