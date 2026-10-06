# build_ti10_fixture.jl — goldens for test/unit/test_ocop_timeint10.jl from live_ti10.sh runs.
# usage: julia --project=<FVSjl> build_ti10_fixture.jl OUT V:STAND:REGIME ...   (OUT = the live_ti10.sh output root)
# writes test/fixtures/tiered_ti10/<v>/<stand>_<regime>_ti10.{key,live.sum,tables,<Table>.csv} (+ stands.db copied from the tiered fixture).
include(joinpath(@__DIR__, "..", "tiered", "tiered_common.jl"))
out = ARGS[1]; root = joinpath(@__DIR__, "..", "..", "fixtures")
for spec in ARGS[2:end]
    v, st, rg = split(spec, ':'); src = joinpath(out, "$(v)_$(st)_$(rg)", "live"); fx = joinpath(root, "tiered_ti10", lowercase(v))
    mkpath(fx); cp(joinpath(root, "tiered", lowercase(v), "stands.db"), joinpath(fx, "stands.db"); force = true)
    stem = "$(st)_$(rg)_ti10"
    write(joinpath(fx, stem * ".key"), replace(read(joinpath(src, "s.key"), String), "out_live.db" => "out.db"))
    cp(joinpath(src, "s.sum"), joinpath(fx, stem * ".live.sum"); force = true)
    t = db_tables(joinpath(src, "out_live.db")); write(joinpath(fx, stem * ".tables"), join(t, "\n") * "\n")
    for x in t; h, r = db_table_rows(joinpath(src, "out_live.db"), x); write_csv(joinpath(fx, "$stem.$x.csv"), h, r); end
    println(stem, " ", t)
end
