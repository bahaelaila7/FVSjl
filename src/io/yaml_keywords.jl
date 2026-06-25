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

# =============================================================================
# STRUCTURED keyword form: each list entry is a single-key map `{KEYWORD: {named
# params}}`, parameters are NAMED, and values keep their YAML type (numbers stay
# numbers). The schema below maps each keyword's param names to its 1-based field
# position (FVS fields are positional, so e.g. STDINFO's stand_origin is field 9).
# A keyword not in the schema still works via the positional/`raw` forms.
# =============================================================================
const _KW_SCHEMA = Dict{String,Vector{Pair{String,Int}}}(
    "DESIGN"   => ["basal_area_factor"=>1, "fixed_plot_area_inverse"=>2, "break_dbh"=>3,
                   "number_of_plots"=>4, "nonstockable_code"=>5, "sample_weight"=>6,
                   "stockable_proportion"=>7],
    "STDINFO"  => ["forest_code"=>1, "habitat"=>2, "stand_age"=>3, "aspect"=>4,
                   "slope"=>5, "elevation"=>6, "stand_origin"=>9],
    "SITECODE" => ["site_species"=>1, "site_index"=>2],
    "INVYEAR"  => ["year"=>1],
    "NUMCYCLE" => ["cycles"=>1],
    "THINBBA"  => ["year"=>1, "residual_basal_area"=>2, "cut_efficiency"=>3,
                   "dbh_min"=>4, "dbh_max"=>5, "species"=>6, "plot"=>7],
    "THINABA"  => ["year"=>1, "residual_basal_area"=>2, "cut_efficiency"=>3,
                   "dbh_min"=>4, "dbh_max"=>5, "species"=>6, "plot"=>7],
    "THINSDI"  => ["year"=>1, "residual_sdi"=>2, "cut_efficiency"=>3,
                   "dbh_min"=>4, "dbh_max"=>5, "species"=>6, "plot"=>7],
    "THINDBH"  => ["year"=>1, "dbh_min"=>2, "dbh_max"=>3, "cut_efficiency"=>4,
                   "residual_tpa"=>6, "species"=>7],
    "THINBTA"  => ["year"=>1, "residual_tpa"=>2, "cut_efficiency"=>3,
                   "dbh_min"=>4, "dbh_max"=>5, "species"=>6, "plot"=>7],
    "THINATA"  => ["year"=>1, "residual_tpa"=>2, "cut_efficiency"=>3,
                   "dbh_min"=>4, "dbh_max"=>5, "species"=>6, "plot"=>7],
    "THINCC"   => ["year"=>1, "residual_ccf"=>2, "cut_efficiency"=>3,
                   "dbh_min"=>4, "dbh_max"=>5, "species"=>6, "plot"=>7],
    # growth / mortality modifiers (date, species, value, dbh window)
    "BAIMULT"  => ["year"=>1, "species"=>2, "multiplier"=>3, "dbh_min"=>4, "dbh_max"=>5],
    "HTGMULT"  => ["year"=>1, "species"=>2, "multiplier"=>3, "dbh_min"=>4, "dbh_max"=>5],
    "CRNMULT"  => ["year"=>1, "species"=>2, "multiplier"=>3],
    "REGDMULT" => ["year"=>1, "species"=>2, "multiplier"=>3],
    "REGHMULT" => ["year"=>1, "species"=>2, "multiplier"=>3],
    "MORTMULT" => ["year"=>1, "species"=>2, "multiplier"=>3, "dbh_min"=>4, "dbh_max"=>5],
    "FIXMORT"  => ["year"=>1, "species"=>2, "rate"=>3, "dbh_min"=>4, "dbh_max"=>5, "option"=>6],
    "FIXDG"    => ["year"=>1, "species"=>2, "value"=>3, "dbh_min"=>4, "dbh_max"=>5],
    "FIXHTG"   => ["year"=>1, "species"=>2, "value"=>3, "dbh_min"=>4, "dbh_max"=>5],
    # density / site
    "BAMAX"    => ["bamax"=>1],
    "SDIMAX"   => ["species"=>1, "sdimax"=>2, "pct_lo"=>5, "pct_hi"=>6],
    "SDICALC"  => ["method"=>1],
    "FERTILIZ" => ["year"=>1, "nitrogen"=>2],
    "SETSITE"  => ["year"=>1, "habitat"=>2, "bamax"=>3, "species"=>4, "site_index"=>5,
                   "si_flag"=>6, "sdimax"=>7],
    # calibration / RNG
    "RANNSEED" => ["seed"=>1],
    "NOCALIB"  => ["species"=>1],
    "DGSTDEV"  => ["value"=>1],
    # establishment cards
    "PLANT"    => ["year"=>1, "species"=>2, "tpa"=>3, "survival_pct"=>4, "age"=>5,
                   "height"=>6, "shade"=>7],
    "NATURAL"  => ["year"=>1, "species"=>2, "tpa"=>3, "survival_pct"=>4, "age"=>5,
                   "height"=>6, "shade"=>7],
    "ESTAB"    => ["disturbance_date"=>1],
    # volume
    "VOLEQNUM" => ["species"=>1, "equation"=>2],
    # cycle / setup extras
    "TIMEINT"  => ["cycle"=>1, "length"=>2],
    "CYCLEAT"  => ["year"=>1],
    "MANAGED"  => ["year"=>1],
    "TFIXAREA" => ["area"=>1],
    "NUMTRIP"  => ["count"=>1],
)
# Keywords whose single named value is carried on the NEXT (free-form) line in a
# .key file, so the structured block expands to TWO records: the keyword + a raw line.
const _KW_TRAILING = Dict("STDIDENT" => "id", "TREEFMT" => "format")

# Build a positional KeywordRecord from a name + a field-index→text map.
function _kw_record_from_fields(name::AbstractString, fld::AbstractDict{Int,String})
    nm = uppercase(rpad(strip(name), 8))
    fields = fill(" "^10, N_KEY_FIELDS); values = zeros(Float32, N_KEY_FIELDS)
    present = falses(N_KEY_FIELDS)
    for (i, txt) in fld
        (i < 1 || i > N_KEY_FIELDS) && continue
        fields[i] = rpad(txt, 10); present[i] = !isempty(strip(txt))
        if _looks_numeric(rpad(txt, 10))
            v = tryparse(Float32, strip(txt)); v === nothing || (values[i] = v)
        end
    end
    status = strip(nm) == "STOP" ? KW_STOP : KW_OK
    return KeywordRecord(nm, "", fields, values, present, N_KEY_FIELDS, status, 0)
end

# A free-form line record (rendered verbatim by `_render_keyfile`) — for the lines
# that trail STDIDENT/TREEFMT.
_raw_record(text) = _decode_keyword(rpad(string(text), 130))

# A YAML scalar → the 10-col field text. Integers print without a trailing ".0".
_yaml_field_text(v) = v isa Integer ? string(v) :
                      v isa AbstractFloat ? (v == floor(v) ? string(Int(v)) : string(v)) :
                      string(v)

# Is this entry the STRUCTURED form? A single-key map whose key is not one of the
# positional/raw reserved keys, with a map (or empty) value.
function _is_structured(entry)
    entry isa AbstractDict || return false
    ks = collect(keys(entry)); length(ks) == 1 || return false
    k = string(ks[1])
    return !(k in ("keyword", "raw", "params", "parms_field"))
end

# Expand ONE structured entry `{KEYWORD: {named params}}` into 1–2 KeywordRecords.
function _records_from_structured(entry::AbstractDict)
    name = string(first(keys(entry)))
    body = first(values(entry))
    body === nothing && (body = Dict{Any,Any}())
    out = KeywordRecord[]
    # keyword carrying its value on the next line(s) (STDIDENT id / TREEFMT format)
    if haskey(_KW_TRAILING, uppercase(strip(name)))
        key = _KW_TRAILING[uppercase(strip(name))]
        push!(out, _kw_record_from_fields(name, Dict{Int,String}()))
        val = string(body isa AbstractDict ? get(body, key, "") : "")
        if uppercase(strip(name)) == "TREEFMT"
            # kw_treefmt! reads EXACTLY two ≤80-col lines and concatenates them, so a
            # FORMAT longer than 80 chars must be split. Break at a comma ≤72 in.
            cut = length(val) <= 72 ? length(val) :
                  (something(findlast(==(','), val[1:72]), 72))
            push!(out, _raw_record(val[1:cut]))
            push!(out, _raw_record(cut < length(val) ? val[cut+1:end] : ""))
        else
            push!(out, _raw_record(val))
        end
        return out
    end
    fld = Dict{Int,String}()
    if body isa AbstractDict
        schema = get(_KW_SCHEMA, uppercase(strip(name)), Pair{String,Int}[])
        pos = Dict(p.first => p.second for p in schema)
        for (pname, pval) in body
            idx = get(pos, string(pname), 0)
            idx > 0 && (fld[idx] = _yaml_field_text(pval))
        end
    end
    push!(out, _kw_record_from_fields(name, fld))
    return out
end

"""
    read_keywords_yaml(path) -> Vector{KeywordRecord}

Read keyword records from a YAML keyword file. Supports both the positional form
(`keyword:`/`params:`) and the **structured** form (`{KEYWORD: {named params}}`),
and `raw:` free-form lines; the two forms may be mixed.
"""
function read_keywords_yaml(path::AbstractString)
    doc = YAML.load_file(path)
    entries = get(doc, "keywords", Any[])
    recs = KeywordRecord[]
    for e in entries
        if _is_structured(e)
            append!(recs, _records_from_structured(e))
        else
            push!(recs, _record_from_yaml(e))
        end
    end
    return recs
end
