# Stratified reservoir sampler: N random STAND_CN per target variant from FVS_STANDINIT_COND,
# restricted to stands that actually have trees (TREEINIT_COND). One streaming scan, seeded.
# Usage: julia --project=. sample_stands.jl <N_per_variant> <seed> <outdir>
using SQLite, Random
const DB = "/workspace/SQLite_FIADB_ENTIRE.db"
const TARGETS = ("SN","NE","CS","LS")
function main(n::Int, seed::Int, outdir::String)
    db = SQLite.DB("file:$DB?mode=ro&immutable=1")
    rng = MersenneTwister(seed)
    # reservoir per variant over STAND_CN
    res = Dict(v => String[] for v in TARGETS); seen = Dict(v => 0 for v in TARGETS)
    for row in DBInterface.execute(db, "SELECT STAND_CN, VARIANT FROM FVS_STANDINIT_COND")
        row[2] === missing && continue
        v = uppercase(strip(String(row[2]))); v in TARGETS || continue
        row[1] === missing && continue
        cn = String(row[1])
        seen[v] += 1; r = res[v]
        if length(r) < n; push!(r, cn)
        else j = rand(rng, 1:seen[v]); j <= n && (r[j] = cn) end
    end
    mkpath(outdir)
    for v in TARGETS
        open(joinpath(outdir, "$(lowercase(v))_sample.txt"), "w") do io
            for cn in res[v]; println(io, cn, "\t", v); end
        end
        println("$v: sampled $(length(res[v])) of $(seen[v]) total")
    end
end
main(parse(Int, ARGS[1]), parse(Int, ARGS[2]), ARGS[3])
