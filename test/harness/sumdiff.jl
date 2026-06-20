# sumdiff.jl — tolerance-aware .sum comparison for the oracle gate.
#
# Two FVS .sum files agree if every numeric field matches within a small tolerance.
# FVS volumes/biomass go through Float32 transcendentals (log/exp/pow), so a strict
# byte diff over-flags ±1 ulp noise (e.g. board-foot volume 1219 vs 1220) that is not
# a real model difference. We allow abs diff ≤ 1 OR relative diff ≤ 0.1% per field.
#
# Usage: julia sumdiff.jl <a.sum> <b.sum>   → prints "MATCH" / "DIFF ..." ; exit 0/1.

function rows(path)
    out = Vector{Vector{Float64}}()
    for ln in eachline(path)
        occursin("-999", ln) && continue
        toks = split(strip(ln))
        isempty(toks) && continue
        nums = Float64[]
        ok = true
        for t in toks
            v = tryparse(Float64, t)
            v === nothing && (ok = false; break)
            push!(nums, v)
        end
        ok && push!(out, nums)
    end
    return out
end

a = rows(ARGS[1]); b = rows(ARGS[2])
atol = 1.0; rtol = 0.001
if length(a) != length(b)
    println("DIFF row count: $(length(a)) vs $(length(b))"); exit(1)
end
for (r, (ra, rb)) in enumerate(zip(a, b))
    length(ra) == length(rb) || (println("DIFF row $r width"); exit(1))
    for (c, (x, y)) in enumerate(zip(ra, rb))
        d = abs(x - y)
        if d > atol && d > rtol * max(abs(x), abs(y))
            println("DIFF row $r col $c: $x vs $y"); exit(1)
        end
    end
end
println("MATCH")
