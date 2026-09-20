# classify.jl — per (variant,regime,file,col) precision-class counts of a report.jl TSV (see classify_fns.jl).
# Usage: julia test/harness/tiered/classify.jl <report.tsv>
include(joinpath(@__DIR__, "classify_fns.jl"))

function main(path)
    agg = Dict{Tuple{String,String,String,String},Dict{String,Int}}()
    for (i, l) in enumerate(eachline(path))
        i == 1 && continue
        f = split(l, '\t'); length(f) < 8 && continue
        k = (f[1], f[3], f[4], f[5]); c = cls(f[7], f[8])
        d = get!(agg, k, Dict{String,Int}()); d[c] = get(d, c, 0) + 1
    end
    for k in sort(collect(keys(agg)))
        println(join(k, "\t"), "\t", join(["$c=$n" for (c, n) in sort(collect(agg[k]))], " "))
    end
end
main(ARGS[1])
