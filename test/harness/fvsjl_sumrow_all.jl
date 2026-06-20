# fvsjl_sumrow_all.jl — run FVSjl's cycle-0 .sum row for every scenario .key and
# print "name | <row>" (or "name | ERROR <msg>"). Loads FVSjl once and loops, so
# it's the fast FVSjl leg of the C10 quality-control sweep. Compared against the
# oracle's .sum cycle-0 row by sweep_compare.jl.
using FVSjl

dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "scenarios")
keys = sort(filter(f -> endswith(f, ".key"), readdir(dir)))

for k in keys
    name = replace(k, ".key" => "")
    path = joinpath(dir, k)
    try
        s, _ = initialize(path)
        FVSjl.notre!(s)
        FVSjl.setup_growth!(s)
        FVSjl.compute_forest_type!(s)
        FVSjl.compute_volumes!(s)
        r = FVSjl.summary_row(s; period = 5)
        io = IOBuffer(); FVSjl.write_sum_row(io, r)
        println(name, " | ", rstrip(String(take!(io))))
    catch e
        println(name, " | ERROR ", sprint(showerror, e))
    end
end
