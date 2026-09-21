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
    pool_path = joinpath(tmp, "pool.txt"); extract(v, 3K, pool_path)
    pool = [split(strip(l), '\t')[1] for l in eachline(pool_path) if !isempty(strip(l))]
    order = vcat(pool[2:3:end], pool[1:3:end], pool[3:3:end])   # K evenly spread first, then the rest
    pdb = joinpath(tmp, "pool.db"); build_subdb(pool, pdb)
    iy = invyears(pdb)
    chosen = String[]; excluded = String[]
    probe = mktempdir(); cp(pdb, joinpath(probe, "stands.db"))
    for cn in order
        length(chosen) >= K && break
        write(joinpath(probe, "s.key"), tiered_keytext(cn, "none", get(iy, cn, 0) + 10))
        txt, crashed = run_oracle(bin, probe)
        if isempty(sum_rows(txt)) || crashed
            push!(excluded, cn); continue
        end
        push!(chosen, cn)
    end
    length(chosen) == K || @warn "$v: only $(length(chosen)) runnable stands (wanted $K)"
    # the committed sub-DB holds exactly the chosen stands
    sdb = joinpath(fx, "stands.db"); build_subdb(chosen, sdb)
    let h = SQLite.DB(sdb); DBInterface.execute(h, "VACUUM"); SQLite.close(h); end
    run_dir = mktempdir(); cp(sdb, joinpath(run_dir, "stands.db"))
    open(joinpath(fx, "stands.txt"), "w") do io; foreach(c -> println(io, c), chosen); end
    nfail = 0
    for cn in chosen, r in REGIMES
        stem = "$(cn)_$(r)"
        key = tiered_keytext(cn, r, get(iy, cn, 0) + 10)
        write(joinpath(fx, stem * ".key"), key)
        write(joinpath(run_dir, "s.key"), key)
        txt, crashed = run_oracle(bin, run_dir)
        rows = sum_rows(txt)
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
    open(joinpath(fx, "PROVENANCE.toml"), "w") do io
        println(io, "variant = \"$v\"")
        println(io, "k = $K")
        println(io, "oracle = \"$bin\"")
        println(io, "oracle_sha256 = \"$(_sha256(bin))\"")
        println(io, "oracle_build_script = \"$(isfile(bscript) ? bscript : "n/a (prebuilt/relinked)")\"")
        println(io, "oracle_main_std_legacy = \"$(isfile(bscript) ? (occursin("NOLEGACY_MAIN", read(bscript, String)) ? "no (NOLEGACY_MAIN)" : "yes — main.f compiled -std=legacy (TREEFMT back-tab bug; FIA-DB path unaffected)") : "unknown")\"")
        println(io, "fvs_source_commit = \"$(_git(fvs, "rev-parse", "HEAD"))\"")
        println(io, "fvs_source_dirty = \"$(isempty(_git(fvs, "status", "--porcelain", "--untracked-files=no")) ? "no" : "yes")\"")
        println(io, "generator_commit = \"$(_git(here, "rev-parse", "HEAD"))\"")
        println(io, "generated = \"$(Base.Libc.strftime("%Y-%m-%d %H:%M:%S", Base.time()))\"")
        println(io, "regimes = [", join(["\"$r\"" for r in REGIMES], ", "), "]")
        println(io, "stands = [", join(["\"$c\"" for c in chosen], ", "), "]")
        println(io, "excluded_no_live_output = [", join(["\"$c\"" for c in excluded], ", "), "]")
        println(io, "live_no_output_cases = $nfail")
    end
    println("$v: $(length(chosen)) stands × $(length(REGIMES)) regimes → $fx  (live no-output cases: $nfail, excluded: $(length(excluded)))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("usage: make_fixtures.jl <VARIANT> <K> [outroot]")
    make_fixtures(ARGS[1], parse(Int, ARGS[2]), length(ARGS) >= 3 ? ARGS[3] : TIERED_ROOT)
end
