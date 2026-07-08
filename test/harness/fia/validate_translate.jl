# Dev harness: prove translate_fia_* reproduces the DB's FVS-ready rows BIT-EXACT on the
# derivable fields (the external-reference fields — SITE_INDEX/LOCATION/design/DG — are
# out of scope; see src/io/fia_translate.jl). Needs the full FIA DB (READ-ONLY).
using FVSjl, SQLite, DBInterface

const DB = get(ENV, "FIA_DB", "/workspace/SQLite_FIADB_ENTIRE.db")

row_dict(r) = Dict{String,Any}(uppercase(String(n)) => getproperty(r, n) for n in propertynames(r))
one(db, sql, a...) = (rs = collect(DBInterface.execute(db, sql, a...)); isempty(rs) ? nothing : row_dict(rs[1]))

# tolerant equality for rendered numeric/text values
same(a, b) = (a === missing || a === nothing) ? (b === missing || b === nothing) :
             (b === missing || b === nothing) ? false :
             (a isa Number && b isa Number) ? Float64(a) == Float64(b) : strip(string(a)) == strip(string(b))

function validate(stand_cns)
    db = SQLite.DB("file:$(DB)?mode=ro&immutable=1")
    tree_fields = ["SPECIES","DIAMETER","HT","CRRATIO","CRCLASS","HISTORY","CULL","PLOT_ID","TREE_ID"]
    stand_fields = ["INV_YEAR","AGE","ASPECT","SLOPE","ELEVFT","SITE_SPECIES",
                    "SITE_INDEX_BASE_AG","FOREST_TYPE_FIA","PHYSIO_REGION","STATE","COUNTY"]
    thit = Dict(f=>0 for f in tree_fields); ttot = Dict(f=>0 for f in tree_fields)
    shit = Dict(f=>0 for f in stand_fields); stot = Dict(f=>0 for f in stand_fields)
    for cn in stand_cns
        # STAND: raw COND (CN=STAND_CN) + PLOT (CN via FVS-ready PLOT_CN) → compare to FVS-ready
        fr = one(db, "SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = ?", (cn,)); fr === nothing && continue
        craw = one(db, "SELECT * FROM COND WHERE CN = ?", (cn,))
        plot_cn = FVSjl._raw(fr, "PLOT_CN", nothing)
        praw = plot_cn === nothing ? Dict{String,Any}() : (one(db, "SELECT * FROM PLOT WHERE CN = ?", (plot_cn,)))
        (craw === nothing || praw === nothing) && continue
        tr = FVSjl.translate_fia_stand(craw, praw)
        for f in stand_fields
            haskey(fr, f) || continue
            stot[f] += 1; same(get(tr, f, missing), fr[f]) && (shit[f] += 1)
        end
        # TREES: each FVS-ready tree ⋈ raw TREE on TREE_CN
        for r in DBInterface.execute(db, "SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = ?", (cn,))
            frt = row_dict(r); tcn = FVSjl._raw(frt, "TREE_CN", nothing); tcn === nothing && continue
            traw = one(db, "SELECT * FROM TREE WHERE CN = ?", (tcn,)); traw === nothing && continue
            traw["TREE_COUNT"] = FVSjl._raw(frt, "TREE_COUNT", 1.0)   # design-derived; supplied
            tt = FVSjl.translate_fia_tree(traw)
            for f in tree_fields
                haskey(frt, f) || continue
                ttot[f] += 1; same(get(tt, f, missing), frt[f]) && (thit[f] += 1)
            end
        end
    end
    println("=== TREE fields: translated vs FVS-ready (bit-exact hit / total) ===")
    for f in tree_fields; println("  ", rpad(f,12), thit[f], " / ", ttot[f]); end
    println("=== STAND direct fields (bit-exact hit / total) ===")
    for f in stand_fields; println("  ", rpad(f,20), shit[f], " / ", stot[f]); end
end

validate(ARGS)
