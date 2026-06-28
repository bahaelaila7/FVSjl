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
function KeywordReader(io::IO)
    # Normalize line endings to LF so `readline` works for every style FVS accepts: LF (Unix),
    # CRLF (DOS, e.g. snt01.key), AND CR-only (old Mac, e.g. net01.key — 0 newlines, all CR). Without
    # this a CR-only keyfile reads as one giant line and the whole parse desyncs.
    content = replace(read(io, String), "\r\n" => "\n", "\r" => "\n")
    KeywordReader(IOBuffer(content), 0, false)
end

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
        # keyrdr.f:55-59: STOP is matched on the FULL 8-char field (TMP=RECORD(1:8), TMP.EQ.'STOP'), NOT a
        # 4-char prefix (that's END, keyrdr.f:92 TMP(1:4)). So "STOPPED"/"STOPxxx" must NOT trigger STOP.
        head8 == "STOP    " && return _record("STOP", record, N_KEY_FIELDS, false, 0, KW_STOP, 0)

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
    # Non-ASCII free-text lines (an em-dash in a STDIDENT id, a Unicode comment, …) carry no
    # fixed-column keyword fields; the byte-indexed column slicing below would land mid-character
    # and throw. Carry them verbatim as a raw record (name = first ≤8 chars, char-safe — only used
    # for the plain/non-plain test, so `.raw` is re-emitted verbatim). Keyword CARDS are ASCII.
    if !isascii(record)
        nm = upkey(rpad(String(first(record, 8)), 8))
        return KeywordRecord(nm, record, fill(" "^10, N_KEY_FIELDS), zeros(Float32, N_KEY_FIELDS),
                             falses(N_KEY_FIELDS), N_KEY_FIELDS, KW_OK, 0)
    end
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

# Index of the last present field (so trailing blank fields are dropped on write).
_last_present_field(rec::KeywordRecord) = (i = findlast(rec.present); i === nothing ? 0 : i)

"""
    write_keyfile(records, path) -> path

Write keyword records back to a legacy fixed-column `.key` file: the keyword name in
cols 1-8 and each present field left-justified in its 10-column slot (field `i` at cols
`10i+1`..`10i+10`), trailing blank fields dropped. This reproduces a card that the KEYRDR
reader (`_decode_keyword`) parses back to the same name / values / presence — a *semantic*
round-trip (the original's exact column padding is not meaningful and is not preserved,
matching the YAML form). Lets a modern YAML keyword file be converted back for legacy FVS.
"""
# A "plain" keyword record reconstructs faithfully from name + 10-col fields. Free-form
# supplemental lines (a TREEFMT FORMAT string, inline tree data) do not — they carry
# punctuation in the keyword columns — so those are round-tripped by their raw text.
_is_plain_keyword(rec::KeywordRecord) =
    occursin(r"^[A-Za-z][A-Za-z0-9]*$", strip(rec.name)) || isempty(strip(rec.name))

# Render keyword records to fixed-column `.key` text (the inverse of KEYRDR decode).
function _render_keyfile(io::IO, records::AbstractVector{KeywordRecord})
    for rec in records
        if !_is_plain_keyword(rec) && !isempty(rec.raw)
            println(io, rstrip(rec.raw))                # free-form line: emit verbatim
            continue
        end
        name = strip(rec.name)
        np = _last_present_field(rec)
        line = rpad(name, 10)                           # cols 1-8 name (+ blank 9-10)
        for i in 1:np
            line *= rpad(strip(rec.fields[i]), 10)      # field i at cols 10i+1 .. 10i+10
        end
        println(io, rstrip(line))
    end
end

"Render keyword records to a `.key`-format string (so a record list can feed `KeywordReader`)."
keyfile_string(records::AbstractVector{KeywordRecord}) =
    (io = IOBuffer(); _render_keyfile(io, records); String(take!(io)))

function write_keyfile(records::AbstractVector{KeywordRecord}, path::AbstractString)
    open(path, "w") do io
        _render_keyfile(io, records)
    end
    return path
end
