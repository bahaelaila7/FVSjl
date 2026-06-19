# =============================================================================
# input.jl — format-agnostic input loading
#
# The engine consumes `TreeRecord` / `KeywordRecord`, never a file format. These
# entry points pick the right front-end by file extension, so legacy `.tre`/`.key`
# and the modern `.csv`/`.yaml` are fully interchangeable — legacy files keep
# working unchanged (drop-in #0), while CSV/YAML are the readable native form.
# =============================================================================

"""
    read_tree_records(path; fmt=DEFAULT_TREE_FORMAT) -> Vector{TreeRecord}

Load tree records from `.tre` (fixed-column, using `fmt`) or `.csv` (named header).
"""
function read_tree_records(path::AbstractString; fmt::AbstractString = DEFAULT_TREE_FORMAT)
    ext = lowercase(splitext(path)[2])
    ext == ".csv" ? read_trees_csv(path) : read_tree_file(path; fmt = fmt)
end

"""
    convert_tre_to_csv(tre, csv; fmt=DEFAULT_TREE_FORMAT) -> csv

Translate a legacy `.tre` file into the readable CSV form (lossless w.r.t. the
parsed tree records).
"""
function convert_tre_to_csv(tre::AbstractString, csv::AbstractString;
                            fmt::AbstractString = DEFAULT_TREE_FORMAT)
    write_trees_csv(read_tree_file(tre; fmt = fmt), csv)
end

"""
    read_keyfile_records(path) -> Vector{KeywordRecord}

Read every record of a legacy `.key` file via the KEYRDR lexer, up to EOF or STOP
(STOP is included as a terminating record).
"""
function read_keyfile_records(path::AbstractString)
    recs = KeywordRecord[]
    open(path) do io
        r = KeywordReader(io)
        while true
            rec = read_keyword!(r)
            rec.status == KW_EOF && break
            push!(recs, rec)
            rec.status == KW_STOP && break
        end
    end
    return recs
end

"""
    read_keyword_records(path) -> Vector{KeywordRecord}

Load keyword records from a legacy `.key` file or a modern `.yaml` file.
"""
function read_keyword_records(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    (ext == ".yaml" || ext == ".yml") ? read_keywords_yaml(path) : read_keyfile_records(path)
end

"""
    convert_key_to_yaml(key, yaml) -> yaml

Translate a legacy `.key` file into the readable, order-preserving YAML form.
"""
function convert_key_to_yaml(key::AbstractString, yaml::AbstractString)
    write_keywords_yaml(read_keyfile_records(key), yaml)
end
