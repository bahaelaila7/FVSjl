# =============================================================================
# keyword.jl — reading keyword (.key) records
#
# Ported from: base/keywd.f  (KEYRDR + UPKEY + PARMS handling).
#
# A keyword record is fixed-column text: columns 1-8 are the keyword name, then up
# to 12 ten-column parameter fields (cols 11-130). A field that looks numeric is
# also parsed to a Float32. The reader skips `!`/`*` comment lines, blank lines,
# and `COMMENT ... END` blocks, and recognizes `STOP` and end-of-file.
#
# This is *just* the lexer — what each keyword DOES (dispatch, supplemental data
# lines, keyword blocks) is the init layer's job (C2). The original KEYRDR also
# echoed records to the listing file; that presentation concern is handled by the
# listing writer, not here, so this stays a pure reader.
# =============================================================================

"Number of 10-column parameter fields in a keyword record (cols 11-130)."
const N_KEY_FIELDS = 12

"""
    KeywordStatus

Outcome of reading one record:
  `KW_OK`   — a normal keyword record
  `KW_EOF`  — end of file
  `KW_STOP` — a `STOP` record
  `KW_PARMS`— a record with a PARMS continuation (see `parms_field`)
"""
@enum KeywordStatus KW_OK KW_EOF KW_STOP KW_PARMS

"""
    KeywordRecord

One parsed keyword record. `fields` are the raw 10-char field texts (Fortran KARD),
`values` their numeric parse (ARRAY; 0 when non-numeric), `present` whether each
field is non-blank (LNOTBK). `parms_field` is the 1-based field where a PARMS
continuation begins (only when `status == KW_PARMS`).
"""
struct KeywordRecord
    name::String
    raw::String
    fields::Vector{String}
    values::Vector{Float32}
    present::Vector{Bool}
    nfields::Int
    status::KeywordStatus
    parms_field::Int
end

"Uppercase a keyword name (UPKEY)."
upkey(s::AbstractString) = uppercase(s)

"""
    KeywordReader(io)

Stateful reader over an open keyword file. `record_count` mirrors Fortran IRECNT
(used in error messages / record numbering). `heading_done` mirrors the inverse of
KEYRDR's `lflag` (the "heading printed yet" flag) — it gates a quirk where blank
lines are skipped before the first keyword but returned as blank-name records after.
"""
mutable struct KeywordReader
    io::IO
    record_count::Int
    heading_done::Bool
end
KeywordReader(io::IO) = KeywordReader(io, 0, false)

@inline _looks_numeric(field::AbstractString) =
    all(c -> isspace(c) || occursin(c, " .+-eE0123456789"), field)

# Locate a PARMS continuation in cols 11-73; returns the number of leading fields
# (nf) and whether PARMS was found. Mirrors keywd.f:134-161.
function _scan_parms(rec::AbstractString)
    k = 11
    while k <= 73
        seg = rec[k:min(73, length(rec))]
        ip = findfirst(c -> c == 'P' || c == 'p', seg)
        ip === nothing && break
        pos = k + ip - 1
        if pos + 4 <= length(rec) && uppercase(rec[pos:pos+4]) == "PARMS"
            return div(pos - 11, 10), true       # nf, found
        end
        k += ip
    end
    return N_KEY_FIELDS, false
end

"""
    read_keyword!(reader) -> KeywordRecord

Read the next keyword record. Mirrors keywd.f:KEYRDR control flow exactly so the
keyword stream (and IRECNT record numbering) is identical to Fortran:
  * `!` comment lines are always skipped;
  * blank lines (any) are skipped *before* the first keyword (heading phase);
  * after that, `*` comments and space-padded blanks are skipped, but a truly
    empty (length-0) line is returned as a blank-name record;
  * `COMMENT … END` blocks are consumed; `STOP` and EOF are reported.
"""
function read_keyword!(r::KeywordReader)
    while true
        eof(r.io) && return _record("", "", N_KEY_FIELDS, false, 0, KW_EOF, 0)
        record = readline(r.io)
        r.record_count += 1

        # ! comment — always skipped (keywd.f:56)
        !isempty(record) && record[1] == '!' && continue
        # blank line before the heading is printed (keywd.f:61)
        !r.heading_done && strip(record) == "" && continue

        head8 = upkey(rpad(length(record) >= 8 ? record[1:8] : record, 8))
        head8[1:4] == "STOP" && return _record("STOP", record, N_KEY_FIELDS, false, 0, KW_STOP, 0)

        # first non-comment, non-heading-blank record prints the heading (keywd.f:76)
        r.heading_done = true

        # * comment or space-padded blank (but NOT a length-0 line) — skipped (keywd.f:89)
        !isempty(record) && (record[1] == '*' || strip(record) == "") && continue

        if head8[1:7] == "COMMENT"                       # consume through END (keywd.f:103)
            while !eof(r.io)
                c = readline(r.io); r.record_count += 1
                uppercase(rpad(length(c) >= 4 ? c[1:4] : c, 4)) == "END " && break
            end
            continue
        end

        return _decode_keyword(record)                   # length-0 → blank-name record
    end
end

# Decode a non-comment record into name + fields (keywd.f:133-208).
function _decode_keyword(record::AbstractString)
    rec = rpad(record, 130)
    nf, found_parms = _scan_parms(rec)

    name = upkey(rpad(rec[1:8], 8))
    fields  = fill(" "^10, N_KEY_FIELDS)
    values  = zeros(Float32, N_KEY_FIELDS)
    present = falses(N_KEY_FIELDS)

    j = 1
    for fi in 1:nf
        j += 10
        field = rpad(rec[j:min(j + 9, length(rec))], 10)
        fields[fi] = field
        if _looks_numeric(field)
            v = tryparse(Float32, strip(field))
            v === nothing || (values[fi] = v)
        end
        present[fi] = strip(field) != ""
    end

    status = KW_OK
    parms_field = 0
    if found_parms && nf < N_KEY_FIELDS
        status = KW_PARMS
        parms_field = nf + 1
        for fi in (nf + 1):7                              # supplemental fields (keywd.f:201)
            j2 = fi * 10 + 1
            fields[fi] = rpad(rec[j2:min(j2 + 9, length(rec))], 10)
        end
    end

    return KeywordRecord(name, record, fields, values, present, nf, status, parms_field)
end

_record(name, raw, nf, _present, _vals, status, pf) = KeywordRecord(
    name, raw, fill(" "^10, N_KEY_FIELDS), zeros(Float32, N_KEY_FIELDS),
    falses(N_KEY_FIELDS), nf, status, pf,
)
