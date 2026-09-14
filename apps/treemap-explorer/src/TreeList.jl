"""
    TreeList

Random-access reader for the TreeMap 2022 CONUS tree table
(`TreeMap2022_CONUS_Tree_Table.csv`, ~2.39M rows, sorted ascending by TM_ID).

A prebuilt offset index (`derived/treetable_index.bin`) maps each TM_ID to the
byte offset of its first row and its row count, so one plot's tree list is a
single seek + sequential read — no full-file scan, no in-memory 238 MB table.

The tree table already carries everything FVS TREEINIT needs, keyed by
TM_ID / PLT_CN: SPCD, DIA, HT, ACTUALHT, CR, TPA_UNADJ. No separate FIA tree DB.
"""
module TreeList

export TreeRecord, TreeTable, treelist, plt_cn

const HDR = ["TM_ID","PLT_CN","STATUSCD","TPA_UNADJ","SPCD","COMMON_NAME",
             "SCIENTIFIC_NAME","SPECIES_SYMBOL","DIA","HT","ACTUALHT","CR",
             "SUBP","TREE","AGENTCD"]

"One tree record from the TreeMap tree table (FIA TREE-derived)."
struct TreeRecord
    tm_id::Int32
    plt_cn::Int64
    statuscd::Int8        # 1 live, 2 dead
    tpa_unadj::Float64    # trees/acre this record represents
    spcd::Int16           # FIA species code
    common_name::String
    species_symbol::String
    dia::Float64          # in
    ht::Float64           # ft
    actualht::Float64     # ft
    cr::Float64           # compacted crown ratio, percent
    subp::Int16
    tree::Int16
    agentcd::Int16        # cause of death (dead only), -1 if NA
end

struct TreeTable
    csv_path::String
    tm_id::Vector{Int32}       # sorted
    offset::Vector{Int64}
    nrows::Vector{Int32}
    pos::Dict{Int32,Int}       # tm_id -> index into the arrays
end

"""
    TreeTable(csv_path, index_path)

Load the offset index built by `scripts/build_index.py` and bind it to the CSV.
"""
function TreeTable(csv_path::AbstractString, index_path::AbstractString)
    io = open(index_path, "r")
    n = read(io, UInt32) |> Int
    tm = Vector{Int32}(undef, n)
    off = Vector{Int64}(undef, n)
    nr = Vector{Int32}(undef, n)
    for i in 1:n
        tm[i]  = reinterpret(Int32, read(io, UInt32))
        off[i] = reinterpret(Int64, read(io, UInt64))
        nr[i]  = reinterpret(Int32, read(io, UInt32))
    end
    close(io)
    pos = Dict{Int32,Int}(tm[i] => i for i in 1:n)
    TreeTable(String(csv_path), tm, off, nr, pos)
end

# --- minimal CSV field splitter honoring double-quoted fields (commas inside quotes) ---
function _split_csv(line::AbstractString)
    fields = String[]
    buf = IOBuffer()
    inq = false
    i = firstindex(line)
    n = lastindex(line)
    while i <= n
        c = line[i]
        if inq
            if c == '"'
                # peek for escaped quote ""
                if i < n && line[nextind(line,i)] == '"'
                    write(buf, '"'); i = nextind(line, i)
                else
                    inq = false
                end
            else
                write(buf, c)
            end
        else
            if c == '"'
                inq = true
            elseif c == ','
                push!(fields, String(take!(buf)))
            else
                write(buf, c)
            end
        end
        i = nextind(line, i)
    end
    push!(fields, String(take!(buf)))
    fields
end

_pint(s) = (s == "NA" || isempty(s)) ? -1 : parse(Int, s)
_pflt(s) = (s == "NA" || isempty(s)) ? NaN  : parse(Float64, s)

function _parse_row(line::AbstractString)::TreeRecord
    f = _split_csv(line)
    TreeRecord(
        parse(Int32, f[1]),
        parse(Int64, f[2]),
        Int8(_pint(f[3])),
        _pflt(f[4]),
        Int16(_pint(f[5])),
        f[6],
        f[8],
        _pflt(f[9]),
        _pflt(f[10]),
        _pflt(f[11]),
        _pflt(f[12]),
        Int16(_pint(f[13])),
        Int16(_pint(f[14])),
        Int16(_pint(f[15])),
    )
end

"Return the tree list for one TM_ID (empty if the plot is not present)."
function treelist(tt::TreeTable, tm_id::Integer)::Vector{TreeRecord}
    i = get(tt.pos, Int32(tm_id), 0)
    i == 0 && return TreeRecord[]
    out = Vector{TreeRecord}(undef, tt.nrows[i])
    open(tt.csv_path, "r") do io
        seek(io, tt.offset[i])
        for k in 1:tt.nrows[i]
            out[k] = _parse_row(readline(io))
        end
    end
    out
end

end # module
