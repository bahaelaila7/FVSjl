# =============================================================================
# treedata.jl — reading tree-data (.tre) records
#
# Ported from: base/intree.f  (the FORTRAN-format machinery + record parse).
#
# Tree records are fixed-column text whose layout is given by a Fortran FORMAT
# string (TREFMT). We (1) parse that FORMAT string into a list of column fields,
# then (2) slice each data line by those fields. Both steps are pure functions —
# no globals, no I/O — so they unit-test trivially against Oracle A.
#
# The default Southern layout (sn/blkdat.f) is:
#   (I4,T1,I7,F6.0,I1,A3,F4.1,F3.1,2F3.0,F4.1,I1,3(I2,I2),2I1,I2,2I3,2I1,F3.0)
# Note the `T1` tab: I4 reads cols 1-4, then the column pointer resets to 1 and I7
# reads cols 1-7 — i.e. the id field overlaps the plot field. This is intentional
# in FVS and the parser reproduces it.
# =============================================================================

"Default Southern-variant tree-record FORMAT string (sn/blkdat.f)."
const DEFAULT_TREE_FORMAT =
    "(I4,T1,I7,F6.0,I1,A3,F4.1,F3.1,2F3.0,F4.1,I1,3(I2,I2),2I1,I2,2I3,2I1,F3.0)"

"One column field extracted from a Fortran FORMAT string."
struct FormatField
    col1::Int        # first column (1-based, inclusive)
    col2::Int        # last column (inclusive)
    kind::Symbol     # :int | :float | :string
    decimals::Int    # implied decimal places (F descriptor) when text has no '.'
end

# Split a FORMAT body on top-level commas (not inside parens).
function _split_format(s::AbstractString)
    toks = String[]; depth = 0; start = 1
    for (i, c) in enumerate(s)
        if c == '('; depth += 1
        elseif c == ')'; depth -= 1
        elseif c == ',' && depth == 0
            push!(toks, s[start:i-1]); start = i + 1
        end
    end
    push!(toks, s[start:end])
    return toks
end

function _expand_format!(fields::Vector{FormatField}, s::AbstractString, col::Base.RefValue{Int})
    for tok in _split_format(s)
        tok = strip(tok)
        isempty(tok) && continue
        up = uppercase(tok)
        if (m = match(r"^(\d+)\((.+)\)$", tok)) !== nothing      # repeated group  n(...)
            for _ in 1:parse(Int, m.captures[1])
                _expand_format!(fields, m.captures[2], col)
            end
        elseif (m = match(r"^(\d+)([IFA])(.*)$", up)) !== nothing # repeated descriptor nXw
            for _ in 1:parse(Int, m.captures[1])
                _single_format!(fields, m.captures[2] * m.captures[3], col)
            end
        else
            _single_format!(fields, up, col)
        end
    end
end

function _single_format!(fields::Vector{FormatField}, tok::AbstractString, col::Base.RefValue{Int})
    if startswith(tok, "T")            # tab — set absolute column
        col[] = parse(Int, tok[2:end])
    elseif startswith(tok, "I")
        w = parse(Int, tok[2:end])
        push!(fields, FormatField(col[], col[] + w - 1, :int, 0)); col[] += w
    elseif startswith(tok, "F")
        m = match(r"^F(\d+)\.(\d+)$", tok)
        if m !== nothing
            w = parse(Int, m.captures[1]); d = parse(Int, m.captures[2])
            push!(fields, FormatField(col[], col[] + w - 1, :float, d)); col[] += w
        end
    elseif startswith(tok, "A")
        w = parse(Int, tok[2:end])
        push!(fields, FormatField(col[], col[] + w - 1, :string, 0)); col[] += w
    elseif startswith(tok, "X")        # horizontal skip
        col[] += length(tok) > 1 ? parse(Int, tok[2:end]) : 1
    end
end

"""
    parse_tree_format(fmt) -> Vector{FormatField}

Parse a Fortran FORMAT string into ordered column fields. Pure.
"""
function parse_tree_format(fmt::AbstractString)
    s = strip(fmt)
    if !isempty(s) && s[1] == '('       # strip the outer parens
        depth = 0
        for (i, c) in enumerate(s)
            c == '(' && (depth += 1)
            c == ')' && (depth -= 1)
            if depth == 0
                s = strip(s[2:i-1]); break
            end
        end
    end
    fields = FormatField[]
    _expand_format!(fields, s, Ref(1))
    return fields
end

# --- field extractors (Fortran implied-decimal semantics for F) ---------------
@inline function _field_int(rec::AbstractString, f::FormatField)::Int32
    f.col1 > length(rec) && return Int32(0)
    s = strip(rec[f.col1:min(f.col2, length(rec))])
    isempty(s) && return Int32(0)
    v = tryparse(Int32, s); v === nothing ? Int32(0) : v
end

@inline function _field_float(rec::AbstractString, f::FormatField)::Float32
    f.col1 > length(rec) && return 0.0f0
    s = strip(rec[f.col1:min(f.col2, length(rec))])
    isempty(s) && return 0.0f0
    v = tryparse(Float32, s)
    v === nothing && return 0.0f0
    # explicit decimal point overrides the format's implied decimals
    (occursin('.', s) || f.decimals == 0) ? v : v / Float32(10.0^f.decimals)
end

@inline function _field_str(rec::AbstractString, f::FormatField, width::Int)::String
    f.col1 > length(rec) && return " "^width
    rpad(strip(rec[f.col1:min(f.col2, length(rec))]), width)
end

"""
    TreeRecord

One parsed tree-data record, with readable field names (Fortran name in comment).
`species_code` is the raw 8-char species text (alpha/FIA/PLANTS) — it is resolved
to a species index later, in INTREE, against the variant's species tables.
"""
struct TreeRecord
    plot::Int32              # point/plot number              (ITREI)
    id::Int32                # tree id                        (IDTREE)
    tpa::Float32             # trees per acre                 (PROB)
    history::Int32           # history/status code            (ITH)
    species_code::String     # raw species text (8 wide)      (CSPI)
    dbh::Float32             # diameter breast height         (DBH)
    diam_growth::Float32     # measured diameter growth       (DG)
    height::Float32          # total height                   (HT)
    top_height::Float32      # height to top (broken/dead)    (THT)
    ht_growth::Float32       # measured height growth         (HTG)
    crown_pct::Int32         # crown ratio percent            (ICR)
    damage::NTuple{6,Int32}  # 3 damage-agent/severity pairs  (IDAMCD)
    mort_code::Int32         # mortality code                 (IMC1)
    cut_code::Int32          # cut/removal code               (KUTKOD)
    pest_vars::NTuple{5,Int32} # pest extension variables     (IPVARS)
    birth_age::Float32       # age at birth, if supplied      (ABIRTH)
end

"""
    parse_tree_record(fields, line) -> TreeRecord | nothing

Slice one tree-data `line` using pre-parsed `fields`. Returns `nothing` if the
format yields too few fields. Field order matches intree.f lines 185-188.
"""
function parse_tree_record(fields::Vector{FormatField}, line::AbstractString)
    nf = length(fields)
    nf == 0 && return nothing
    maxcol = maximum(f.col2 for f in fields)
    rec = rpad(line, max(maxcol, 80))
    fi(n) = _field_int(rec, fields[n])
    ff(n) = _field_float(rec, fields[n])
    return TreeRecord(
        fi(1), fi(2), ff(3), fi(4), _field_str(rec, fields[5], 8),
        ff(6), ff(7), ff(8), ff(9), ff(10), fi(11),
        (fi(12), fi(13), fi(14), fi(15), fi(16), fi(17)),
        fi(18), fi(19),
        (fi(20), fi(21), fi(22), fi(23), fi(24)),
        nf >= 25 ? ff(25) : 0.0f0,
    )
end

"Parse every line of a `.tre` file using `fmt` (defaults to the SN layout)."
function read_tree_file(path::AbstractString; fmt::AbstractString = DEFAULT_TREE_FORMAT)
    fields = parse_tree_format(fmt)
    recs = TreeRecord[]
    for line in eachline(path)
        isempty(strip(line)) && continue
        r = parse_tree_record(fields, line)
        r === nothing || push!(recs, r)
    end
    return recs
end
