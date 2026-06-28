# =============================================================================
# yaml_stand.jl — the SEMANTIC stand format (`format: fvs-stand/v1`)
#
# A second, idiomatic YAML form that describes a stand by INTENT rather than by
# mirroring the ordered `.key` keyword stream. Where `yaml_keywords.jl` keeps the
# keyword sequence (just grouped), this form is a declarative map — `invyr`,
# `numcycle`, `stdinfo`, `treatments`, `treelist` … — and the reader "unravels" it
# into FVS keyword records in the canonical emission order (so the engine and a
# `.key` round-trip both get a correctly-ordered stream).
#
# Discriminator: the top-level `format:` key. `fvs-stand/*` ⇒ this reader;
# anything else (or its absence) ⇒ the keyword-stream reader in yaml_keywords.jl.
#
# Keys mirror FVS keyword/parameter names (invyr, numcycle, thinsdi …) so an FVS
# user maps them at a glance; full field→FVS-card mapping is in docs/FORMATS.md.
#
# Coverage: the common keywords are modeled below; ANYTHING ELSE rides along in a
# `raw_keywords:` list (the keyword-stream entry form) emitted just before the tree
# list — so the format is usable and lossless from day one and extends over time.
#
# Single stand:           Multiple stands:
#   format: fvs-stand/v1     format: fvs-stand/v1
#   stand:                   stands:
#     invyr: 1990              - { invyr: 1990, ... }
#     ...                      - { invyr: 1995, ... }
# =============================================================================

# Boolean flag keywords: `key: true` ⇒ emit the bare card; `false`/absent ⇒ nothing.
const _SEM_FLAG = ["noautoes"=>"NOAUTOES", "screen"=>"SCREEN", "notriple"=>"NOTRIPLE",
                   "echosum"=>"ECHOSUM"]

# Single-scalar keywords: `key: <v>` ⇒ emit the card with <v> in the given field.
const _SEM_SCALAR = ["invyr"=>("INVYEAR",1), "numcycle"=>("NUMCYCLE",1),
                     "bamax"=>("BAMAX",1), "sdicalc"=>("SDICALC",1),
                     "managed"=>("MANAGED",1), "numtrip"=>("NUMTRIP",1),
                     "rannseed"=>("RANNSEED",1), "tfixarea"=>("TFIXAREA",1)]

# Map keywords: `key: {param: v, …}` ⇒ one card; params map to FVS field positions.
const _SEM_MAP = ["design"   => ("DESIGN", ["baf"=>1,"fixed_plot"=>2,"break_dbh"=>3,
                                            "plots"=>4,"nsc"=>5,"sample_weight"=>6,"stockable"=>7]),
                  "stdinfo"  => ("STDINFO", ["forest"=>1,"habitat"=>2,"age"=>3,"aspect"=>4,
                                             "slope"=>5,"elev"=>6,"origin"=>9]),
                  "sitecode" => ("SITECODE", ["species"=>1,"index"=>2])]

# Treatment keywords (each `treatments:` list entry is `{<key>: {params}}`). FVS-term
# param names; the cut activities schedule in the list's author order (same-year order
# matters, so the list is preserved verbatim).
const _SEM_TREAT = Dict(
    "thinbba" =>("THINBBA", ["year"=>1,"ba"=>2,"eff"=>3,"dmin"=>4,"dmax"=>5,"species"=>6,"plot"=>7]),
    "thinaba" =>("THINABA", ["year"=>1,"ba"=>2,"eff"=>3,"dmin"=>4,"dmax"=>5,"species"=>6,"plot"=>7]),
    "thinbta" =>("THINBTA", ["year"=>1,"tpa"=>2,"eff"=>3,"dmin"=>4,"dmax"=>5,"species"=>6,"plot"=>7]),
    "thinata" =>("THINATA", ["year"=>1,"tpa"=>2,"eff"=>3,"dmin"=>4,"dmax"=>5,"species"=>6,"plot"=>7]),
    "thinsdi" =>("THINSDI", ["year"=>1,"sdi"=>2,"eff"=>3,"dmin"=>4,"dmax"=>5,"species"=>6,"plot"=>7]),
    "thincc"  =>("THINCC",  ["year"=>1,"ccf"=>2,"eff"=>3,"dmin"=>4,"dmax"=>5,"species"=>6,"plot"=>7]),
    "thinrden"=>("THINRDEN",["year"=>1,"rsdi"=>2,"eff"=>3,"dmin"=>4,"dmax"=>5,"species"=>6,"plot"=>7]),
    "thindbh" =>("THINDBH", ["year"=>1,"dmin"=>2,"dmax"=>3,"eff"=>4,"tpa"=>6,"species"=>7]),
    "thinht"  =>("THINHT",  ["year"=>1,"hmin"=>2,"hmax"=>3,"eff"=>4,"tpa"=>6,"species"=>7]),
    "thinauto"=>("THINAUTO",["year"=>1,"eff"=>2]),
    "salvage" =>("SALVAGE", ["year"=>1,"dmin"=>2,"dmax"=>3,"eff"=>4,"species"=>5]),
)

# Regeneration list keywords (`regeneration:` entries `{plant|natural: {params}}`).
const _SEM_REGEN = Dict(
    "plant"  =>("PLANT",  ["year"=>1,"species"=>2,"tpa"=>3,"survival"=>4,"age"=>5,"height"=>6,"shade"=>7]),
    "natural"=>("NATURAL",["year"=>1,"species"=>2,"tpa"=>3,"survival"=>4,"age"=>5,"height"=>6,"shade"=>7]),
)

# ---- field helpers -------------------------------------------------------------
# Build a card from {pos => text}; reuse the keyword-stream record builder.
_sem_card(name, fld::AbstractDict{Int,String}) = _kw_record_from_fields(name, fld)

# Map a `{param=>value}` body through a `[param=>pos]` schema → `{pos => fieldtext}`.
function _sem_fields(body, schema)
    fld = Dict{Int,String}()
    body isa AbstractDict || return fld
    pos = Dict(string(p.first) => p.second for p in schema)
    for (k, v) in body
        i = get(pos, string(k), 0)
        i > 0 && v !== nothing && (fld[i] = _yaml_field_text(v))
    end
    return fld
end

# Append the cards for one list-of-`{kw: {params}}` (treatments / regeneration).
function _emit_tagged_list!(recs, list, table)
    list === nothing && return
    for entry in (list isa AbstractVector ? list : Any[list])
        entry isa AbstractDict || continue
        for (tag, body) in entry
            spec = get(table, lowercase(string(tag)), nothing)
            if spec === nothing
                error("fvs-stand: unknown entry `$tag` (no semantic mapping; put it in raw_keywords)")
            end
            name, schema = spec
            push!(recs, _sem_card(name, _sem_fields(body, schema)))
        end
    end
end

# ---- one stand → ordered KeywordRecords ----------------------------------------
# Emit in canonical FVS order so the flattened stream satisfies every hard ordering
# constraint (STDIDENT first; species groups before the treatments that name them;
# COMPUTE before IF; tree list + PROCESS last). raw_keywords ride just before the
# tree list (good for extension blocks: ECON, FFE/FMIN, COMPUTE — all pre-PROCESS).
function _stand_records!(recs::Vector{KeywordRecord}, st::AbstractDict)
    g(k) = get(st, k, nothing)
    # 1. STDIDENT (+ its id line)
    if g("stdident") !== nothing
        push!(recs, _sem_card("STDIDENT", Dict{Int,String}()))
        push!(recs, _raw_record(string(g("stdident"))))
    end
    # 2. control flags
    for (k, kw) in _SEM_FLAG
        k == "echosum" && continue                    # ECHOSUM emitted with output (12) below
        g(k) === true && push!(recs, _sem_card(kw, Dict{Int,String}()))
    end
    # 3-5. design / stand info / site
    for (k, spec) in _SEM_MAP
        b = g(k); b === nothing && continue
        name, schema = spec
        push!(recs, _sem_card(name, _sem_fields(b, schema)))
    end
    _emit_tagged_list!(recs, g("setsite"),
        Dict("setsite"=>("SETSITE",["year"=>1,"habitat"=>2,"bamax"=>3,"species"=>4,
                                    "index"=>5,"flag"=>6,"sdimax"=>7])))
    # 6-7. inventory year & cycling
    for (k, spec) in _SEM_SCALAR
        v = g(k); v === nothing && continue
        name, fpos = spec
        push!(recs, _sem_card(name, Dict(fpos => _yaml_field_text(v))))
    end
    _emit_tagged_list!(recs, g("timeint"), Dict("timeint"=>("TIMEINT",["cycle"=>1,"length"=>2])))
    if g("cycleat") !== nothing
        for y in (g("cycleat") isa AbstractVector ? g("cycleat") : Any[g("cycleat")])
            push!(recs, _sem_card("CYCLEAT", Dict(1 => _yaml_field_text(y))))
        end
    end
    # 8. density
    _emit_tagged_list!(recs, g("sdimax"),
        Dict("sdimax"=>("SDIMAX",["species"=>1,"sdimax"=>2,"pct_lo"=>5,"pct_hi"=>6])))
    # 9. species groups (define BEFORE treatments)
    if g("spgroup") !== nothing
        for grp in (g("spgroup") isa AbstractVector ? g("spgroup") : Any[g("spgroup")])
            grp isa AbstractDict || continue
            push!(recs, _sem_card("SPGROUP", Dict(1 => string(get(grp, "name", "")))))
            members = get(grp, "species", nothing)
            members !== nothing &&
                push!(recs, _raw_record(members isa AbstractVector ?
                      join((string(m) for m in members), " ") : string(members)))
        end
    end
    # 10. growth modifiers (date/species/value/window cards)
    _emit_tagged_list!(recs, g("growth"), Dict(
        "baimult"=>("BAIMULT",["year"=>1,"species"=>2,"mult"=>3,"dmin"=>4,"dmax"=>5]),
        "htgmult"=>("HTGMULT",["year"=>1,"species"=>2,"mult"=>3,"dmin"=>4,"dmax"=>5]),
        "mortmult"=>("MORTMULT",["year"=>1,"species"=>2,"mult"=>3,"dmin"=>4,"dmax"=>5]),
        "fixmort"=>("FIXMORT",["year"=>1,"species"=>2,"rate"=>3,"dmin"=>4,"dmax"=>5,"option"=>6]),
        "fixdg"=>("FIXDG",["year"=>1,"species"=>2,"value"=>3,"dmin"=>4,"dmax"=>5])))
    # 12. output
    g("echosum") === true && push!(recs, _sem_card("ECHOSUM", Dict{Int,String}()))
    g("summary") === true && push!(recs, _sem_card("SUMMARY", Dict{Int,String}()))
    # 13. raw_keywords escape hatch (keyword-stream entry form; emitted verbatim)
    if g("raw_keywords") !== nothing
        for e in (g("raw_keywords") isa AbstractVector ? g("raw_keywords") : Any[g("raw_keywords")])
            _append_entry!(recs, e)
        end
    end
    # 14. treatments (author order)
    _emit_tagged_list!(recs, g("treatments"), _SEM_TREAT)
    # 15. regeneration
    _emit_tagged_list!(recs, g("regeneration"), _SEM_REGEN)
    g("sprout") === true && push!(recs, _sem_card("SPROUT", Dict{Int,String}()))
    g("estab")  !== nothing && push!(recs, _sem_card("ESTAB",
        Dict(1 => _yaml_field_text(g("estab") isa AbstractDict ? get(g("estab"),"date","") : g("estab")))))
    # 16. event monitor (COMPUTE before IF — list order preserved)
    if g("event_monitor") !== nothing
        for e in (g("event_monitor") isa AbstractVector ? g("event_monitor") : Any[g("event_monitor")])
            _append_entry!(recs, e)
        end
    end
    # 17. tree list: TREEFMT (the .tre column layout, optional) + TREEDATA. The inventory is
    # NOT inline — it is the companion file resolved by the keyfile's base name (`<key>.csv`
    # preferred, else `<key>.tre`); every stand reads that same file (see read_stand_yaml_doc's
    # REWIND). `treelist.format` only sets the .tre layout (ignored when a .csv companion exists);
    # `treelist:` may be `{}` to just emit TREEDATA.
    tl = g("treelist")
    if tl !== nothing
        fmt = tl isa AbstractDict ? get(tl, "format", nothing) : nothing
        if fmt !== nothing
            push!(recs, _sem_card("TREEFMT", Dict{Int,String}()))
            s = string(fmt)                          # kw_treefmt! reads two ≤80-col lines
            cut = length(s) <= 72 ? length(s) : something(findlast(==(','), s[1:72]), 72)
            push!(recs, _raw_record(s[1:cut]))
            cut < length(s) && push!(recs, _raw_record(s[cut+1:end]))
        end
        push!(recs, _sem_card("TREEDATA", Dict{Int,String}()))
    end
    # 18. PROCESS this stand
    push!(recs, _sem_card("PROCESS", Dict{Int,String}()))
    return recs
end

"""
    read_stand_yaml_doc(doc) -> Vector{KeywordRecord}

Unravel a parsed `format: fvs-stand/*` document (`stand:` map, or `stands:` list of
maps) into the canonical ordered keyword-record stream. A terminating `STOP` is
appended after the last stand.

Multiple stands describe **N scenarios on the same stand** — every stand reads the SAME
companion tree file (`<keyfile>.csv`/`.tre`, resolved by base name). Stock FVS reads that
file sequentially, so without intervention only the first stand would get trees; a
`REWIND 2` is therefore emitted before every stand after the first to re-read the
tree-data unit (snt01's pattern). FVSjl re-reads implicitly (REWIND is a no-op for it),
so the same stream is correct for both — the converted `.key` runs identically in stock FVS.
"""
function read_stand_yaml_doc(doc::AbstractDict)
    recs = KeywordRecord[]
    stands = haskey(doc, "stands") ? doc["stands"] :
             haskey(doc, "stand")  ? Any[doc["stand"]] : Any[]
    first = true
    for st in (stands isa AbstractVector ? stands : Any[stands])
        st isa AbstractDict || continue
        first || push!(recs, _sem_card("REWIND", Dict(1 => "2")))   # re-read the shared tree file
        first = false
        _stand_records!(recs, st)
    end
    push!(recs, _sem_card("STOP", Dict{Int,String}()))
    return recs
end

"""
    yaml_variant_code(path) -> Union{String,Nothing}

Peek a keyword YAML file's top-level `variant:` field (e.g. `variant: SN`) — works for
both YAML flavors. `nothing` when the file has no `variant:` (or isn't a mapping). The
entry points map this through `variant_from_code` so a YAML stand file selects its own
model; a `.key` has no variant (it's the binary choice), so this is YAML-only.
"""
function yaml_variant_code(path::AbstractString)
    doc = YAML.load_file(path)
    doc isa AbstractDict || return nothing
    v = get(doc, "variant", nothing)
    v === nothing ? nothing : string(v)
end

"Is this a semantic (`format: fvs-stand/*`) document?"
function _is_stand_doc(doc)
    doc isa AbstractDict || return false
    f = get(doc, "format", nothing)
    f !== nothing && startswith(lowercase(string(f)), "fvs-stand") && return true
    # tolerate a missing discriminator: `stand:` as a MAP (the keyword-stream form
    # uses `stand:` as a LIST) or a `stands:` list of maps.
    (haskey(doc, "stand") && doc["stand"] isa AbstractDict) && return true
    haskey(doc, "stands") && return true
    return false
end
