# =============================================================================
# yaml_keywords.jl — keywords as YAML (the readable native form of a .key file)
#
# A `.key` file is fixed-column Fortran text (8-col keyword + 12 ten-col fields).
# Its modern equivalent is an ordered YAML sequence — order is preserved because
# FVS keyword order is significant (activities schedule, later keywords override).
# Both produce the same `Vector{KeywordRecord}` the engine dispatches on.
#
# We store each record's keyword name and its non-empty trailing parameters as
# readable strings. The dispatch-relevant content (name, numeric values, presence
# flags, stripped field text) round-trips losslessly; the exact 10-column padding
# of the original card is intentionally not preserved (it carries no meaning the
# keyword handlers use — see PORTING.md). Reading uses YAML.jl; writing is emitted
# by hand for clean, stable formatting.
# =============================================================================

using YAML

# number of trailing params to keep = index of the last present field
_last_present(rec::KeywordRecord) = (i = findlast(rec.present); i === nothing ? 0 : i)

"""
    write_keywords_yaml(records, path)

Write keyword records as an ordered YAML `keywords:` sequence.
"""
function write_keywords_yaml(records::AbstractVector{KeywordRecord}, path::AbstractString)
    open(path, "w") do io
        println(io, "# FVS keywords — order is significant. `params` are the keyword's")
        println(io, "# fields (left to right); omitted = no parameters.")
        println(io, "keywords:")
        for rec in records
            # Free-form supplemental lines (a TREEFMT FORMAT string, inline tree data) are
            # not keyword+params; carry them verbatim so the form stays lossless.
            if !_is_plain_keyword(rec) && !isempty(rec.raw)
                println(io, "  - raw: ", _yq(rstrip(rec.raw)))
                continue
            end
            name = strip(rec.name)
            np = _last_present(rec)
            print(io, "  - keyword: ", _yq(name))
            if np > 0
                params = (strip(rec.fields[i]) for i in 1:np)
                print(io, "\n    params: [", join((_yq(p) for p in params), ", "), "]")
            end
            if rec.status == KW_PARMS
                print(io, "\n    parms_field: ", rec.parms_field)
            end
            println(io)
        end
    end
    return path
end

# YAML-quote a string (always quote → unambiguous for codes like "60." or "9999").
_yq(s) = '"' * replace(String(s), "\\" => "\\\\", "\"" => "\\\"") * '"'

"Reconstruct a KeywordRecord from a YAML entry (name + optional params, or a raw line)."
function _record_from_yaml(entry::AbstractDict)
    if haskey(entry, "raw")                          # free-form line: decode verbatim
        return _decode_keyword(rpad(string(entry["raw"]), 130))
    end
    name = uppercase(rpad(strip(string(get(entry, "keyword", ""))), 8))
    raw_params = get(entry, "params", Any[])
    params = String[strip(string(p)) for p in raw_params]

    fields  = fill(" "^10, N_KEY_FIELDS)
    values  = zeros(Float32, N_KEY_FIELDS)
    present = falses(N_KEY_FIELDS)
    for (i, p) in enumerate(params)
        i > N_KEY_FIELDS && break
        fields[i] = rpad(p, 10)
        present[i] = !isempty(p)
        if _looks_numeric(rpad(p, 10))
            v = tryparse(Float32, p)
            v === nothing || (values[i] = v)
        end
    end

    pf = get(entry, "parms_field", 0)
    status = strip(name) == "STOP" ? KW_STOP : (pf == 0 ? KW_OK : KW_PARMS)
    return KeywordRecord(name, "", fields, values, present, N_KEY_FIELDS, status, Int(pf))
end

"""
    read_keywords_yaml(path) -> Vector{KeywordRecord}

Read keyword records from a YAML file produced by `write_keywords_yaml`.
"""
function read_keywords_yaml(path::AbstractString)
    doc = YAML.load_file(path)
    entries = get(doc, "keywords", Any[])
    return KeywordRecord[_record_from_yaml(e) for e in entries]
end
