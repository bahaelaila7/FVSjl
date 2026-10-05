# make_fixtures.jl — build the committed tiered-suite fixture for one variant from the LIVE oracle.
#
# Usage: julia --project=. test/harness/tiered/make_fixtures.jl <VARIANT> <K> [outroot]
#   Picks K stratified FIA stands (extract_sample.jl's deterministic even stride over ECOREGION/LOCATION), writes a
#   self-contained sub-DB, one keyfile per stand × regime (tiered_common.REGIMES), runs the live oracle (BIN from
#   ledger_fia.jl — FIA-DB path = DBSTREESIN, unaffected by the g16 -std=legacy TREEFMT read bug) and stores the
#   goldens: <cn>_<regime>.live.sum, <cn>_<regime>.<Table>.csv, <cn>_<regime>.tables, plus PROVENANCE.toml.
#   A stand whose live oracle produces no .sum under `none` is replaced by the next stratified candidate
#   (deterministic), and recorded in PROVENANCE as excluded.
#   outroot defaults to test/fixtures/tiered. The slow tier passes a temp dir to diff regenerated goldens.
include(joinpath(@__DIR__, "tiered_common.jl"))
include(joinpath(@__DIR__, "..", "fia", "extract_sample.jl"))

# Variants with NO FIA population of their own draw their stands from another variant's FIA stands (deterministic, so
# the slow tier regenerates the same sample): (source VARIANT, extra SQL predicate on the stand row `s`).
#   KT: FVS_STANDINIT_COND.VARIANT is never 'KT' (KT = Kootenai/Kaniksu/Tally Lake, superseded by IE). Its home
#       forests — 110 Flathead (Tally Lake RD), 113 Idaho Panhandle (Kaniksu), 114 Kootenai — are in KT's forkod.f
#       JFOR table and their FIA stands are assigned VARIANT='IE'.
#   OC: ORGANON Southwest Oregon overlay on CA (OC links ca/forkod.f; JFOR 505..518, 610/611, 710-712). FIA never
#       assigns VARIANT='OC'; its stands are the CA-assigned FIA stands of the two Oregon national forests in that
#       table — 610 Rogue River, 611 Siskiyou (the BLM 710-712 locations have no FIA rows).
#   OP: ORGANON SMC (western Washington) overlay on PN (OP links pn/forkod.f; JFOR 609,612,800,708,709,712). FIA never
#       assigns VARIANT='OP'; its stands are the PN-assigned FIA stands of 609 Olympic NF, its namesake forest.
const SAMPLE_SOURCE = Dict("KT" => ("IE", "AND s.LOCATION IN (110,113,114)"),
                           "OC" => ("CA", "AND s.LOCATION IN (610,611)"),
                           "OP" => ("PN", "AND s.LOCATION IN (609)"))

# BC (British Columbia, metric) has no FIA population at all — FIA is US-only. Its stands come from FVS's own BC test
# databases (tests/FVSbc: FVS-BC.YSM-SkyRanch.db = 15 young-stand-monitoring / SkyRanch spacing-trial stands, IDFdk4;
# FVS-BC.Stand.Structure.Data.db = 2 Fd plantations, IDFxh), copied into FIA-named tables (FVS_STANDINIT_COND /
# FVS_TREEINIT_COND; FVS reads the columns by name, the BC column set is the FIA one). The Stand.Structure trees carry
# the stand id with a ".3200" suffix (its StandStructure.key selects them separately); they are re-keyed to their
# stand so the one-SQL tiered keyfile reads them. Pool order: (Ecoregion, Stand_ID).
const LOCAL_POOL = Dict("BC" => ["/workspace/ForestVegetationSimulator/tests/FVSbc/FVS-BC.YSM-SkyRanch.db",
                                 "/workspace/ForestVegetationSimulator/tests/FVSbc/FVS-BC.Stand.Structure.Data.db"])

# Variants whose K-stand pick is ROUND-ROBIN over the inventory tree-count classes (extract_sample.jl _tree_class:
# 0 | 1-9 | 10-49 | 50+ records). The plain stride over the strata-ordered pool picks whatever classes the stride lands
# on (CS came out 0 bare / 5 of 10 with 1-9 trees), so a fixture can miss the establishment path or the dense mature
# path entirely. Opt-in per variant, so the variants built before this (whose goldens the slow tier regenerates)
# keep their stand list.
const CLASS_BALANCED = Set(["CS", "LS", "NE", "OC", "OP", "ON"])

"`order` re-sequenced round-robin over the tree-count classes (each class keeps its relative order)."
function class_balanced_order(order::Vector, pdb::AbstractString)
    h = SQLite.DB(pdb); nt = Dict{String,Int}()
    for r in DBInterface.execute(h, "SELECT STAND_CN, COUNT(*) FROM FVS_TREEINIT_COND GROUP BY STAND_CN")
        r[1] === missing || (nt[string(r[1])] = Int(r[2]))
    end
    SQLite.close(h)
    byc = [String[] for _ in 1:4]
    for cn in order; push!(byc[_tree_class(get(nt, String(cn), 0))], String(cn)); end
    out = String[]
    while any(!isempty, byc)
        for c in 1:4; isempty(byc[c]) || push!(out, popfirst!(byc[c])); end
    end
    out
end

# ON (Ontario, metric) has no FIA population and FVS ships ONE Ontario stand (tests/FVSon FVSDataHardwood.db, LD3001).
# Its stands are drawn from the Lake States FIA population it borders — northern Minnesota, LOCATION 903 Chippewa /
# 909 Superior (both in ON's forkod.f JFOR table: ON is the LS-derived variant and keeps the LS forest codes), with the
# FIA fixed 1-acre design (BASAL_AREA_FACTOR 0, INV_PLOT_SIZE 1) — and CONVERTED TO METRIC, because the metric DB
# reader (metric/dbsqlite) takes DIAMETER/DG in cm, heights in m, SITE_INDEX in m and TREE_COUNT per hectare of the
# (1/INV_PLOT_SIZE) plot. FIA species codes are read through ON's FIAJSP table (blkdat.f). ELEVFT is left in feet (the
# reader converts it). Both engines read the identical converted sub-DB, so the instrument is oracle-comparable.
const METRIC_POOL = Dict("ON" => ("LS", "AND s.BASAL_AREA_FACTOR = 0 AND s.INV_PLOT_SIZE = 1 AND s.STATE = 27 AND s.LOCATION IN (903,909)"))
const _IN2CM = 2.54; const _FT2M = 0.3048; const _PERAC2PERHA = 2.4710538

function convert_pool_metric!(pdb::AbstractString)
    h = SQLite.DB(pdb)
    DBInterface.execute(h, """UPDATE FVS_TREEINIT_COND SET DIAMETER = DIAMETER * $_IN2CM, DG = DG * $_IN2CM,
        HT = HT * $_FT2M, HTG = HTG * $_FT2M, HTTOPK = HTTOPK * $_FT2M, HTBOLE = HTBOLE * $_FT2M, HTSAW = HTSAW * $_FT2M,
        HT_TO_CROWN_BASE = HT_TO_CROWN_BASE * $_FT2M, TREE_COUNT = TREE_COUNT * $_PERAC2PERHA""")
    DBInterface.execute(h, "UPDATE FVS_STANDINIT_COND SET SITE_INDEX = SITE_INDEX * $_FT2M, VARIANT = 'ON'")
    SQLite.close(h)
end

function build_local_pool(srcs::Vector{String}, pdb::AbstractString)
    isfile(pdb) && rm(pdb)
    h = SQLite.DB(pdb)
    for (k, src) in enumerate(srcs)
        DBInterface.execute(h, "ATTACH DATABASE '$src' AS src")
        if k == 1
            DBInterface.execute(h, "CREATE TABLE FVS_STANDINIT_COND AS SELECT * FROM src.FVS_StandInit")
            DBInterface.execute(h, "CREATE TABLE FVS_TREEINIT_COND AS SELECT * FROM src.FVS_TreeInit")
        else
            DBInterface.execute(h, "INSERT INTO FVS_STANDINIT_COND SELECT * FROM src.FVS_StandInit")
            DBInterface.execute(h, "INSERT INTO FVS_TREEINIT_COND SELECT * FROM src.FVS_TreeInit")
        end
        DBInterface.execute(h, "DETACH DATABASE src")
    end
    for col in ("Stand_CN", "Stand_ID")
        DBInterface.execute(h, """UPDATE FVS_TREEINIT_COND SET $col = (SELECT s.Stand_ID FROM FVS_STANDINIT_COND s
            WHERE FVS_TREEINIT_COND.$col LIKE s.Stand_ID || '.%') WHERE $col NOT IN (SELECT Stand_ID FROM FVS_STANDINIT_COND)""")
    end
    ids = [string(r[1]) for r in DBInterface.execute(h,
           "SELECT Stand_ID FROM FVS_STANDINIT_COND ORDER BY Ecoregion, Stand_ID")]
    SQLite.close(h)
    ids
end

function subset_local_db(pdb::AbstractString, cns::Vector{String}, out::AbstractString)
    isfile(out) && rm(out)
    h = SQLite.DB(out)
    DBInterface.execute(h, "ATTACH DATABASE '$pdb' AS p")
    lst = join(["'" * c * "'" for c in cns], ",")
    DBInterface.execute(h, "CREATE TABLE FVS_STANDINIT_COND AS SELECT * FROM p.FVS_STANDINIT_COND WHERE Stand_CN IN ($lst)")
    DBInterface.execute(h, "CREATE TABLE FVS_TREEINIT_COND AS SELECT * FROM p.FVS_TREEINIT_COND WHERE Stand_CN IN ($lst)")
    DBInterface.execute(h, "DETACH DATABASE p")
    SQLite.close(h)
end

# A regime whose extension the oracle build STUBS (ex*.f in FVS<v>_buildDir, e.g. exrd.f in CA/AK, exclim.f in AK)
# makes FVS print "FVS11 ERROR:  REQUESTED EXTENSION IS NOT PART OF THIS PROGRAM" and ignore the block — the case
# would test nothing. make_fixtures probes each regime on the first chosen stand and drops such regimes (recorded
# in PROVENANCE `skipped_regimes`; the fast tier runs PROVENANCE `regimes` only).
stubbed_extension(outtext::AbstractString) = occursin("FVS11 ERROR", outtext)

function _sha256(path)
    try strip(first(split(read(`sha256sum $path`, String)))) catch; "unavailable" end
end
function _git(dir, args...)
    try strip(read(Cmd(`git -C $dir $args`), String)) catch; "unavailable" end
end

function run_oracle(bin, dir)
    for f in ("s.sum", "s.out", "out.db"); p = joinpath(dir, f); isfile(p) && rm(p); end
    crashed = false
    try
        p = cd(dir) do
            run(pipeline(ignorestatus(`$bin --keywordfile=s.key`); stdout = devnull, stderr = devnull))
        end
        crashed = (p.termsignal != 0) || (p.exitcode > 128)
    catch
        crashed = true
    end
    sp = joinpath(dir, "s.sum")
    (isfile(sp) ? read(sp, String) : "", crashed)
end

function invyears(db)
    d = Dict{String,Int}()
    h = SQLite.DB(db)
    for r in DBInterface.execute(h, "SELECT STAND_CN, INV_YEAR FROM FVS_STANDINIT_COND")
        (r[1] === missing || r[2] === missing) && continue
        d[String(r[1])] = Int(r[2])
    end
    SQLite.close(h); d
end

function make_fixtures(v::AbstractString, K::Int, outroot::AbstractString)
    bin = BIN[v]
    fx = joinpath(outroot, lowercase(v)); mkpath(fx)
    tmp = mktempdir()
    # candidate pool: 3K stratified stands; walk the K-sample first, substituting in stride order on failure
    src, wh = get(SAMPLE_SOURCE, uppercase(v), (uppercase(v), ""))
    local_pool = haskey(LOCAL_POOL, uppercase(v))
    metric_pool = haskey(METRIC_POOL, uppercase(v))
    metric_pool && ((src, wh) = METRIC_POOL[uppercase(v)])
    pdb = joinpath(tmp, "pool.db")
    if local_pool
        pool = build_local_pool(LOCAL_POOL[uppercase(v)], pdb)
    else
        pool_path = joinpath(tmp, "pool.txt"); extract(src, 3K, pool_path; where_extra = wh)
        pool = [split(strip(l), '\t')[1] for l in eachline(pool_path) if !isempty(strip(l))]
        build_subdb(pool, pdb)
        metric_pool && convert_pool_metric!(pdb)
    end
    order = vcat(pool[2:3:end], pool[1:3:end], pool[3:3:end])   # K evenly spread first, then the rest
    uppercase(v) in CLASS_BALANCED && (order = class_balanced_order(order, pdb))
    iy = invyears(pdb)
    chosen = String[]; excluded = String[]
    probe = mktempdir(); cp(pdb, joinpath(probe, "stands.db"))
    for cn in order
        length(chosen) >= K && break
        write(joinpath(probe, "s.key"), tiered_keytext(cn, "none", get(iy, cn, 0) + 10))
        txt, crashed = run_oracle(bin, probe)
        if isempty(sum_rows(txt; metric = metric_sum(v))) || crashed
            push!(excluded, cn); continue
        end
        push!(chosen, cn)
    end
    length(chosen) == K || @warn "$v: only $(length(chosen)) runnable stands (wanted $K)"
    # the committed sub-DB holds exactly the chosen stands
    sdb = joinpath(fx, "stands.db"); (local_pool || metric_pool) ? subset_local_db(pdb, chosen, sdb) : build_subdb(chosen, sdb)
    let h = SQLite.DB(sdb); DBInterface.execute(h, "VACUUM"); SQLite.close(h); end
    run_dir = mktempdir(); cp(sdb, joinpath(run_dir, "stands.db"))
    open(joinpath(fx, "stands.txt"), "w") do io; foreach(c -> println(io, c), chosen); end
    # stub probe: the first chosen stand under each regime
    regimes = String[]; skipped = String[]; probe_err = String[]
    for r in REGIMES
        write(joinpath(run_dir, "s.key"), tiered_keytext(chosen[1], r, get(iy, chosen[1], 0) + 10))
        run_oracle(bin, run_dir)
        op = joinpath(run_dir, "s.out"); ot = isfile(op) ? read(op, String) : ""
        codes = sort(unique([m.match for m in eachmatch(r"FVS\d\d ERROR", ot)]))
        isempty(codes) || push!(probe_err, "$r: $(join(codes, "/"))")
        if stubbed_extension(ot)
            push!(skipped, r); continue
        end
        push!(regimes, r)
    end
    nfail = 0
    for cn in chosen, r in regimes
        stem = "$(cn)_$(r)"
        key = tiered_keytext(cn, r, get(iy, cn, 0) + 10)
        write(joinpath(fx, stem * ".key"), key)
        write(joinpath(run_dir, "s.key"), key)
        txt, crashed = run_oracle(bin, run_dir)
        rows = sum_rows(txt; metric = metric_sum(v))
        if crashed || isempty(rows)
            nfail += 1
            write(joinpath(fx, stem * ".live.sum"), "# LIVE_NO_OUTPUT crashed=$crashed\n")
            write(joinpath(fx, stem * ".tables"), "")
            continue
        end
        write(joinpath(fx, stem * ".live.sum"), join(rows, "\n") * "\n")
        outdb = joinpath(run_dir, "out.db")
        tbls = db_tables(outdb)
        write(joinpath(fx, stem * ".tables"), join(tbls, "\n") * (isempty(tbls) ? "" : "\n"))
        for t in tbls
            hdr, trows = db_table_rows(outdb, t)
            write_csv(joinpath(fx, "$(stem).$(t).csv"), hdr, trows)
        end
    end
    here = normpath(joinpath(@__DIR__, "..", "..", ".."))
    fvs = "/workspace/ForestVegetationSimulator"
    bscript = "/workspace/." * lowercase(v) * "work/build_g16.sh"
    # BC (build_dbfix.sh) / OC, OP, ON (build_tiered.sh): private relinks of the stock objects, each with its own script
    privb = joinpath(dirname(bin), isfile(joinpath(dirname(bin), "build_dbfix.sh")) ? "build_dbfix.sh" : "build_tiered.sh")
    isfile(privb) && (bscript = privb)
    open(joinpath(fx, "PROVENANCE.toml"), "w") do io
        println(io, "variant = \"$v\"")
        println(io, "k = $K")
        println(io, "oracle = \"$bin\"")
        println(io, "oracle_sha256 = \"$(_sha256(bin))\"")
        println(io, "oracle_build_script = \"$(isfile(bscript) ? bscript : "n/a (prebuilt/relinked)")\"")
        println(io, "oracle_main_std_legacy = \"$(bscript == privb ? (occursin("g16obj", read(privb, String)) ? "yes — ON g16 objects (-std=legacy main kept, ORACLE_SOURCE_AUDIT §3; FIA-DB path unaffected)" : "no (stock buildDir main.o)") : isfile(bscript) ? (occursin("NOLEGACY_MAIN", read(bscript, String)) ? "no (NOLEGACY_MAIN)" : "yes — main.f compiled -std=legacy (TREEFMT back-tab bug; FIA-DB path unaffected)") : "unknown")\"")
        println(io, "fvs_source_commit = \"$(_git(fvs, "rev-parse", "HEAD"))\"")
        println(io, "fvs_source_dirty = \"$(isempty(_git(fvs, "status", "--porcelain", "--untracked-files=no")) ? "no" : "yes")\"")
        println(io, "generator_commit = \"$(_git(here, "rev-parse", "HEAD"))\"")
        println(io, "generated = \"$(Base.Libc.strftime("%Y-%m-%d %H:%M:%S", Base.time()))\"")
        println(io, "regimes = [", join(["\"$r\"" for r in regimes], ", "), "]")
        println(io, "skipped_regimes = [", join(["\"$r\"" for r in skipped], ", "), "]   # FVS11: extension stubbed in this oracle build")
        println(io, "probe_fvs_errors = [", join(["\"$e\"" for e in probe_err], ", "), "]   # FVSnn ERROR codes in the stub-probe .out (first stand)")
        println(io, "sample_source = \"", local_pool ? "LOCAL " * join(basename.(LOCAL_POOL[uppercase(v)]), " + ") : (metric_pool ? "METRIC-CONVERTED " : "") * "VARIANT=$(src)$(isempty(wh) ? "" : " " * wh)", "\"")
        println(io, "stands = [", join(["\"$c\"" for c in chosen], ", "), "]")
        println(io, "excluded_no_live_output = [", join(["\"$c\"" for c in excluded], ", "), "]")
        println(io, "live_no_output_cases = $nfail")
    end
    println("$v: $(length(chosen)) stands × $(length(regimes)) regimes → $fx  (skipped (stubbed): $(join(skipped, ",")), live no-output cases: $nfail, excluded: $(length(excluded)))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("usage: make_fixtures.jl <VARIANT> <K> [outroot]")
    make_fixtures(ARGS[1], parse(Int, ARGS[2]), length(ARGS) >= 3 ? ARGS[3] : TIERED_ROOT)
end
