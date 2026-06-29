# =============================================================================
# yaml_keywords.jl — keywords as YAML (the readable native form of a .key file)
#
# A `.key` file is fixed-column Fortran text (8-col keyword + 12 ten-col fields).
# Its modern equivalent is an ORDER-AWARE, HIERARCHICAL YAML document. The whole
# keyword stream stays an ordered SEQUENCE because FVS keyword order is significant
# (activities schedule in input order, SPGROUP must precede a THIN that names it,
# COMPUTE variables must be defined before use, same-cycle activities run in input
# order, later keywords override earlier). The hierarchy only *groups* that ordered
# stream into named, readable sections — it never reorders it.
#
# Two on-disk forms are read (the engine consumes the same `Vector{KeywordRecord}`
# from either, and from a legacy `.key`):
#
#   stand:                      # HIERARCHICAL form (default emitted form). Ordered
#     - <section>:              # list of section blocks; each block is an ordered
#         - <KEYWORD>: {...}    # list of keyword entries. Blocks appear in original
#         - ...                 # order, so the flattened stream is byte-for-byte the
#     - <section>:              # original keyword order — grouping is order-PRESERVING.
#         - ...
#
#   keywords:                   # FLAT form (still accepted for back-compat). A single
#     - <KEYWORD>: {...}        # ordered list of keyword entries.
#     - ...
#
# Each keyword entry keeps NAMED, typed parameters (numbers stay numbers); a keyword
# not in the schema falls back to a positional `params:` list or a verbatim `raw:`
# line. The dispatch-relevant content (name, values, presence, stripped field text)
# round-trips losslessly; the exact 10-column padding of the original card is not
# preserved (it carries no meaning the handlers use — see PORTING.md). Reading uses
# YAML.jl; writing is emitted by hand for clean, stable, ordered formatting.
#
# ORDER-SIGNIFICANT sections (see `_SECTION_ORDER_NOTE`): every section is an ordered
# list, but `treatments` (thinning/harvest/salvage by cycle) and `event_monitor`
# (COMPUTE/IF def-before-use) are the ones whose *relative* order changes results,
# so they are emitted as a single contiguous, clearly-labelled ordered block.
# =============================================================================

using YAML

# number of trailing params to keep = index of the last present field
_last_present(rec::KeywordRecord) = (i = findlast(rec.present); i === nothing ? 0 : i)

# =============================================================================
# Section grouping. Each keyword belongs to ONE named section (purpose-grouped for
# readability). Grouping is ORDER-PRESERVING: when emitting, consecutive records of
# the same section are gathered into one block, and blocks appear in their original
# order, so flattening the hierarchy reproduces the exact keyword sequence — no
# reordering, ever. A keyword not listed here falls into `other` (still ordered).
# =============================================================================
const _KW_SECTION = Dict{String,String}(
    # stand identification & run setup
    "STDIDENT"=>"setup", "STDINFO"=>"setup", "DESIGN"=>"setup", "SITECODE"=>"setup",
    "INVYEAR"=>"setup", "NUMCYCLE"=>"setup", "TIMEINT"=>"setup", "CYCLEAT"=>"setup",
    "MANAGED"=>"setup", "RESETAGE"=>"setup", "TFIXAREA"=>"setup", "TREEFMT"=>"setup",
    "TREEDATA"=>"setup", "NOTREES"=>"setup", "NOAUTOES"=>"setup", "PROCESS"=>"setup",
    "STOP"=>"setup", "DESIGN"=>"setup",
    # density limits
    "SDIMAX"=>"density", "SDICALC"=>"density", "BAMAX"=>"density",
    # growth calibration & modifiers
    "GROWTH"=>"growth", "NOCALIB"=>"growth", "READCORD"=>"growth", "REUSCORD"=>"growth",
    "READCORH"=>"growth", "REUSCORH"=>"growth", "READCORR"=>"growth", "REUSCORR"=>"growth",
    "BAIMULT"=>"growth", "HTGMULT"=>"growth", "CRNMULT"=>"growth", "REGDMULT"=>"growth",
    "REGHMULT"=>"growth", "DGSTDEV"=>"growth", "SERLCORR"=>"growth", "RANNSEED"=>"growth",
    "FIXDG"=>"growth", "FIXHTG"=>"growth", "FIXMORT"=>"growth", "MORTMULT"=>"growth",
    "TREESZCP"=>"growth", "SIZCAP"=>"growth",
    # thinning / harvest (ORDER-SIGNIFICANT — same-cycle activities run in input order)
    "THINBBA"=>"treatments", "THINABA"=>"treatments", "THINBTA"=>"treatments",
    "THINATA"=>"treatments", "THINSDI"=>"treatments", "THINCC"=>"treatments",
    "THINHT"=>"treatments", "THINQFA"=>"treatments", "THINRDEN"=>"treatments",
    "THINDBH"=>"treatments", "THINPT"=>"treatments", "SETPTHIN"=>"treatments",
    "THINPRSC"=>"treatments", "THINAUTO"=>"treatments", "SPECPREF"=>"treatments",
    "LEAVESP"=>"treatments", "SPLEAVE"=>"treatments", "CUTEFF"=>"treatments",
    "MINHARV"=>"treatments", "SALVAGE"=>"treatments", "YARDLOSS"=>"treatments",
    # establishment & regeneration
    "ESTAB"=>"regeneration", "PLANT"=>"regeneration", "NATURAL"=>"regeneration",
    "SPROUT"=>"regeneration", "NOSPROUT"=>"regeneration",
    # volume & merchandising
    "VOLUME"=>"volume", "BFVOLUME"=>"volume", "VOLEQNUM"=>"volume", "MCDEFECT"=>"volume",
    "BFDEFECT"=>"volume", "MCFDLN"=>"volume", "BFFDLN"=>"volume",
    # species groups (ORDER-SIGNIFICANT — SPGROUP must precede a THIN that names it)
    "SPGROUP"=>"species_groups",
    # site & treatments
    "SETSITE"=>"site", "FERTILIZ"=>"site",
    # output & database
    "ECHOSUM"=>"output", "DATABASE"=>"output", "DSNOUT"=>"output", "SUMMARY"=>"output",
    "TREELIDB"=>"output", "CUTLIST"=>"output", "STRCLASS"=>"output", "COMPUTDB"=>"output",
    # event monitor (ORDER-SIGNIFICANT — COMPUTE vars must be defined before use)
    "COMPUTE"=>"event_monitor", "IF"=>"event_monitor",
    # compression & tripling
    "COMPRESS"=>"control", "NOTRIPLE"=>"control", "NUMTRIP"=>"control",
)
# Canonical section ordering, used only as a stable tie-break for documentation; the
# emitter preserves the ORIGINAL keyword order regardless of this list.
const _SECTION_LABEL = Dict(
    "setup"=>"stand setup & inventory", "density"=>"density limits",
    "growth"=>"growth calibration & modifiers", "treatments"=>"treatments (ORDER MATTERS)",
    "regeneration"=>"establishment & regeneration", "volume"=>"volume & merchandising",
    "species_groups"=>"species groups (define before use)", "site"=>"site & fertilization",
    "output"=>"output & database", "event_monitor"=>"event monitor (define before use)",
    "control"=>"compression & tripling", "other"=>"other",
)

# Which section a record belongs to. Free-form / unnamed lines (raw continuation
# lines, an IF condition, a COMPUTE body line, SPGROUP members) inherit the section
# of the record they follow so the block stays contiguous (handled by the emitter).
_section_of(name::AbstractString) = get(_KW_SECTION, uppercase(strip(name)), "other")

"""
    write_keywords_yaml(records, path; flat=false)

Write keyword records as an ORDER-AWARE, HIERARCHICAL YAML document. The default
`stand:` form groups the (still fully ordered) keyword stream into named sections;
`flat=true` emits the legacy single `keywords:` list. Both round-trip losslessly.
"""
function write_keywords_yaml(records::AbstractVector{KeywordRecord}, path::AbstractString;
                             flat::Bool=false)
    open(path, "w") do io
        flat ? _write_flat_yaml(io, records) : _write_hierarchical_yaml(io, records)
    end
    return path
end

# ---- entry emitter: one keyword record → its YAML block (named params if possible) ----
# `indent` is the leading whitespace before the `- `; emits to `io`.
function _emit_entry(io::IO, rec::KeywordRecord, indent::AbstractString)
    # Free-form supplemental lines (TREEFMT FORMAT, an IF condition, a COMPUTE body
    # line, SPGROUP member list, inline tree data) carry verbatim so the form is lossless.
    if !_is_plain_keyword(rec) && !isempty(rec.raw)
        println(io, indent, "- raw: ", _yq(rstrip(rec.raw)))
        return
    end
    name = strip(rec.name)
    np = _last_present(rec)
    schema = get(_KW_SCHEMA, uppercase(name), Pair{String,Int}[])
    # Use the NAMED structured form when (a) the keyword has a schema, (b) no PARMS
    # continuation, and (c) every present field maps to a named slot (otherwise we'd
    # silently drop a positional field). Else fall back to the positional form.
    named_ok = !isempty(schema) && rec.status != KW_PARMS &&
               all(any(p.second == i for p in schema) for i in 1:np if rec.present[i])
    if named_ok && np > 0
        println(io, indent, "- ", name, ":")
        for (pname, pos) in schema
            pos <= np && rec.present[pos] || continue
            println(io, indent, "    ", pname, ": ", _yaml_scalar(strip(rec.fields[pos])))
        end
        return
    end
    # positional / no-schema form
    print(io, indent, "- keyword: ", _yq(name))
    if np > 0
        params = (strip(rec.fields[i]) for i in 1:np)
        print(io, "\n", indent, "  params: [", join((_yq(p) for p in params), ", "), "]")
    end
    if rec.status == KW_PARMS
        print(io, "\n", indent, "  parms_field: ", rec.parms_field)
    end
    println(io)
end

# Render a stripped field string as a YAML scalar. A number stays UNQUOTED (readable,
# natural type) ONLY when it survives the YAML→number→text round-trip byte-for-byte —
# i.e. `_yaml_field_text(YAML.parse(s)) == s`. Forms that lose their exact text under
# that trip (e.g. "60." → 60.0 → "60.0", "1.00" → 1.0) are QUOTED so the field text is
# preserved exactly; that keeps the conversion lossless while most numbers read clean.
function _yaml_scalar(s::AbstractString)
    isempty(s) && return "\"\""
    if tryparse(Int, s) !== nothing
        return s                               # integer text always round-trips
    end
    if _looks_numeric(rpad(s, 1)) && tryparse(Float64, s) !== nothing
        f = parse(Float64, s)
        return _yaml_field_text(f) == s ? s : _yq(s)   # unquoted iff text-stable
    end
    return _yq(s)                              # genuine string code (e.g. "231Dd")
end

# ---- FLAT form (legacy, back-compat) ----
function _write_flat_yaml(io::IO, records::AbstractVector{KeywordRecord})
    println(io, "# FVS keywords — order is significant. Flat ordered list.")
    println(io, "keywords:")
    for rec in records
        _emit_entry(io, rec, "  ")
    end
end

# ---- HIERARCHICAL form (default) ----
# Group consecutive records into section blocks (ORDER-PRESERVING) and emit
# `stand:` as an ordered list of `{section: [entries]}` blocks. A free-form/unnamed
# record (no section) extends the CURRENT block so multi-line keywords (SPGROUP +
# members, IF + condition, COMPUTE + body, TREEFMT + format) stay contiguous.
function _write_hierarchical_yaml(io::IO, records::AbstractVector{KeywordRecord})
    println(io, "# FVS keywords — ORDER-AWARE hierarchical form. The keyword stream is")
    println(io, "# an ordered sequence; sections only GROUP it (flattening top-to-bottom")
    println(io, "# reproduces the exact order). `treatments`, `species_groups` and")
    println(io, "# `event_monitor` are order-significant (define-before-use, same-cycle")
    println(io, "# activities run in input order). See docs/KEYWORDS.md.")
    println(io, "stand:")
    cur = ""           # current section name
    started = false
    for rec in records
        sec = (_is_plain_keyword(rec) && !isempty(strip(rec.name))) ?
              _section_of(rec.name) : cur     # unnamed/raw lines stay in the current block
        if !started || sec != cur
            println(io, "  - ", sec, ":")
            cur = sec; started = true
        end
        _emit_entry(io, rec, "      ")
    end
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
    "THINHT"   => ["year"=>1, "ht_min"=>2, "ht_max"=>3, "cut_efficiency"=>4,
                   "residual_tpa"=>6, "species"=>7],
    "THINRDEN" => ["year"=>1, "residual_relsdi"=>2, "cut_efficiency"=>3,
                   "dbh_min"=>4, "dbh_max"=>5, "species"=>6, "plot"=>7],
    "THINAUTO" => ["year"=>1, "cut_efficiency"=>2],
    "LEAVESP"  => ["species"=>1],
    "SPLEAVE"  => ["species"=>1],
    "SPECPREF" => ["species"=>1, "preference"=>2],
    "CUTEFF"   => ["proportion"=>1],
    "MINHARV"  => ["min_volume"=>1],
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
    "SERLCORR" => ["phi"=>1, "theta"=>2],
    "RESETAGE" => ["year"=>1, "age"=>2],
    # establishment cards
    "PLANT"    => ["year"=>1, "species"=>2, "tpa"=>3, "survival_pct"=>4, "age"=>5,
                   "height"=>6, "shade"=>7],
    "NATURAL"  => ["year"=>1, "species"=>2, "tpa"=>3, "survival_pct"=>4, "age"=>5,
                   "height"=>6, "shade"=>7],
    "ESTAB"    => ["disturbance_date"=>1],
    # volume / defect (merch standards; species 0=all, <0=group)
    "VOLEQNUM" => ["species"=>1, "equation"=>2],
    "VOLUME"   => ["year"=>1, "species"=>2, "dbh_min"=>3, "top_diam"=>4, "stump"=>5,
                   "form_class"=>6, "method"=>7, "scf_min_dbh"=>8, "scf_top_dib"=>9, "scf_stump"=>10],
    "BFVOLUME" => ["year"=>1, "species"=>2, "bf_min_dbh"=>3, "bf_top_dib"=>4, "bf_stump"=>5],
    "MCDEFECT" => ["year"=>1, "species"=>2, "defect_5"=>3, "defect_10"=>4, "defect_15"=>5,
                   "defect_20"=>6, "defect_25"=>7],
    "BFDEFECT" => ["year"=>1, "species"=>2, "defect_5"=>3, "defect_10"=>4, "defect_15"=>5,
                   "defect_20"=>6, "defect_25"=>7],
    # species groups — the group NAME is field 1 of the SPGROUP card; the member
    # species are carried as a `raw:` line that follows (handled as a free-form line).
    "SPGROUP"  => ["group"=>1],
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

# A free-form line record (rendered verbatim by `_render_keyfile`) — for the lines that trail
# STDIDENT/TREEFMT, inline tree data, SPGROUP members. A MULTI-WORD line (a STDIDENT id+title,
# a FORMAT string, a tree-data row) must round-trip byte-for-byte, NOT be re-parsed into 10-col
# fields: `_decode_keyword` would set name=cols-1-8 — if that is a bare token (e.g. "SCENARIO"
# from "SCENARIO1  control"), `_is_plain_keyword` reads it as a keyword card and the re-render
# drops col 9-10 + reflows the rest (mangling it). So force a verbatim record whose name keeps
# a space ⇒ non-plain ⇒ emitted via `.raw`. A bare single token has no space and round-trips as
# a plain card whose render IS the token, so keep using `_decode_keyword` for it (also handles
# non-ASCII without byte-crashing).
function _raw_record(text)
    s = rstrip(string(text))
    occursin(' ', s) || return _decode_keyword(rpad(s, 130))
    return KeywordRecord(s, s, fill(" "^10, N_KEY_FIELDS), zeros(Float32, N_KEY_FIELDS),
                         falses(N_KEY_FIELDS), N_KEY_FIELDS, KW_OK, 0)
end

# A YAML scalar → the 10-col field text. The YAML number TYPE carries the original
# field's form: an integer scalar (`60`) → "60"; a float scalar (`60.0`) → "60.0"
# (Julia's `string(60.0)` keeps the trailing ".0"). This makes the named structured
# form preserve the exact field text — e.g. ".key" `60.0` round-trips as `60.0`, not
# `60` — so the dispatch-relevant text round-trips losslessly while staying readable.
_yaml_field_text(v::Integer) = string(v)
_yaml_field_text(v::AbstractFloat) = string(v)
_yaml_field_text(v) = string(v)

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

# Turn one keyword/structured/raw entry into its KeywordRecords, appending in order.
function _append_entry!(recs::Vector{KeywordRecord}, e)
    if _is_structured(e)
        append!(recs, _records_from_structured(e))
    else
        push!(recs, _record_from_yaml(e))
    end
end

# Flatten the hierarchical `stand:` form (an ordered list of `{section: [entries]}`
# blocks) into the ordered keyword stream. Section *names* are ignored on read — they
# are only a grouping label; the list order IS the keyword order, so this is exact.
function _read_hierarchical(stand)
    recs = KeywordRecord[]
    blocks = stand isa AbstractVector ? stand :          # list-of-blocks (canonical)
             stand isa AbstractDict   ? [Dict(k=>v) for (k,v) in stand] : Any[]
    for block in blocks
        block isa AbstractDict || continue
        for (_sec, entries) in block
            entries === nothing && continue
            for e in (entries isa AbstractVector ? entries : Any[entries])
                _append_entry!(recs, e)
            end
        end
    end
    return recs
end

"""
    read_keywords_yaml(path) -> Vector{KeywordRecord}

Read keyword records from a YAML keyword file. Accepts the ORDER-AWARE hierarchical
form (`stand:` → ordered section blocks) AND the flat form (`keywords:` list); within
either, each entry may be the **structured** `{KEYWORD: {named params}}` map, the
positional `keyword:`/`params:` map, or a verbatim `raw:` line. The list order is the
keyword order in every case, so reading is a lossless inverse of the writer.
"""
function read_keywords_yaml(path::AbstractString)
    doc = YAML.load_file(path)
    _is_stand_doc(doc) && return read_stand_yaml_doc(doc)    # SEMANTIC form (format: fvs-stand/*)
    if doc isa AbstractDict && haskey(doc, "stand")          # hierarchical form
        return _read_hierarchical(doc["stand"])
    end
    entries = doc isa AbstractDict ? get(doc, "keywords", Any[]) : Any[]
    recs = KeywordRecord[]
    for e in entries
        _append_entry!(recs, e)
    end
    return recs
end
