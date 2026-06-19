# =============================================================================
# csv_trees.jl — tree data as CSV (the readable native form of a .tre file)
#
# A `.tre` file is fixed-column Fortran text; its modern equivalent here is a CSV
# with a named header. Both produce the same `Vector{TreeRecord}` (the schema the
# engine consumes), so the simulation never cares which one it came from.
#
# Dep-free: tree fields are numbers + a short species code with no embedded commas,
# so a hand-rolled reader/writer is enough (no CSV.jl needed — keeps the core lean
# and trim-friendly).
# =============================================================================

"CSV column order for a TreeRecord (stable; matches `TreeRecord` field semantics)."
const TREE_CSV_HEADER = [
    "plot", "id", "tpa", "history", "species", "dbh", "diam_growth", "height",
    "top_height", "ht_growth", "crown_pct",
    "damage1", "damage2", "damage3", "damage4", "damage5", "damage6",
    "mort_code", "cut_code", "pest1", "pest2", "pest3", "pest4", "pest5", "birth_age",
]

# Format a Float32 compactly but losslessly enough for round-tripping.
_csv_num(x::Float32) = x == round(x) ? string(Int(x)) : string(x)
_csv_num(x::Int32)   = string(x)

"Render one TreeRecord as a CSV row (25 fields, in TREE_CSV_HEADER order)."
function tree_to_csv_row(r::TreeRecord)
    join((
        _csv_num(r.plot), _csv_num(r.id), _csv_num(r.tpa), _csv_num(r.history),
        strip(r.species_code), _csv_num(r.dbh), _csv_num(r.diam_growth),
        _csv_num(r.height), _csv_num(r.top_height), _csv_num(r.ht_growth),
        _csv_num(r.crown_pct),
        (_csv_num(d) for d in r.damage)..., _csv_num(r.mort_code), _csv_num(r.cut_code),
        (_csv_num(p) for p in r.pest_vars)..., _csv_num(r.birth_age),
    ), ',')
end

"""
    write_trees_csv(records, path)

Write tree records to a CSV file with a named header.
"""
function write_trees_csv(records::AbstractVector{TreeRecord}, path::AbstractString)
    open(path, "w") do io
        println(io, join(TREE_CSV_HEADER, ','))
        for r in records
            println(io, tree_to_csv_row(r))
        end
    end
    return path
end

# Parse a CSV cell to Int32 / Float32 (used positionally per header).
_csv_i(s) = (v = tryparse(Int32, strip(s)); v === nothing ? Int32(0) : v)
_csv_f(s) = (v = tryparse(Float32, strip(s)); v === nothing ? 0.0f0 : v)

"""
    read_trees_csv(path) -> Vector{TreeRecord}

Read tree records from a CSV produced by `write_trees_csv` (header required).
"""
function read_trees_csv(path::AbstractString)
    recs = TreeRecord[]
    lines = readlines(path)
    isempty(lines) && return recs
    header = strip.(split(lines[1], ','))
    @assert header == TREE_CSV_HEADER "unexpected tree CSV header: $header"
    for ln in @view lines[2:end]
        isempty(strip(ln)) && continue
        c = split(ln, ',')
        push!(recs, TreeRecord(
            _csv_i(c[1]), _csv_i(c[2]), _csv_f(c[3]), _csv_i(c[4]),
            rpad(strip(c[5]), 8), _csv_f(c[6]), _csv_f(c[7]), _csv_f(c[8]),
            _csv_f(c[9]), _csv_f(c[10]), _csv_i(c[11]),
            (_csv_i(c[12]), _csv_i(c[13]), _csv_i(c[14]), _csv_i(c[15]), _csv_i(c[16]), _csv_i(c[17])),
            _csv_i(c[18]), _csv_i(c[19]),
            (_csv_i(c[20]), _csv_i(c[21]), _csv_i(c[22]), _csv_i(c[23]), _csv_i(c[24])),
            _csv_f(c[25]),
        ))
    end
    return recs
end
