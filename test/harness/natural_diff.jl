# natural_diff.jl — per-scenario congruence of FVSjl vs Oracle A over the FULL .sum
# (all cycles), reporting the FIRST divergence in the key stand columns. For the
# natural-process sweep (thinning pinned). Usage: julia natural_diff.jl <minedir> <oradir>
function rows(path)
    out = Vector{Vector{Float64}}()
    isfile(path) || return out
    for ln in eachline(path)
        occursin("-999", ln) && continue
        toks = split(strip(ln)); isempty(toks) && continue
        v = tryparse.(Float64, toks)
        all(!isnothing, v) && push!(out, Float64.(v))
    end
    return out
end

function run_sweep(minedir, oradir, atol, rtol)
cols = get(ENV, "DYNONLY", "") == "1" ?
       ((3, "TPA"), (4, "BA"), (5, "SDI"), (7, "TopHt"), (8, "QMD")) :
       ((3, "TPA"), (4, "BA"), (7, "TopHt"), (8, "QMD"), (11, "scuft"), (12, "bdft"))
names = sort([replace(f, ".sum" => "") for f in readdir(minedir) if endswith(f, ".sum")])

npass = 0; ndiff = 0; nerr = 0
ok(x, y) = abs(x - y) <= atol || abs(x - y) <= rtol * max(abs(x), abs(y))
for nm in names
    a = rows(joinpath(minedir, nm * ".sum")); b = rows(joinpath(oradir, nm * ".sum"))
    if isempty(a) || isempty(b)
        nerr += 1; println("ERR  ", nm, " (mine ", length(a), " ora ", length(b), " rows)"); continue
    end
    worst = ""
    for r in 1:min(length(a), length(b))
        for (ci, cn) in cols
            ci <= length(a[r]) && ci <= length(b[r]) || continue
            if !ok(a[r][ci], b[r][ci])
                worst = "cyc$(r-1) $cn " * string(a[r][ci]) * " vs " * string(b[r][ci]); break
            end
        end
        isempty(worst) || break
    end
    if isempty(worst) && length(a) == length(b)
        npass += 1
    else
        ndiff += 1
        println("DIFF ", nm, ": ", isempty(worst) ? "rowcount $(length(a)) vs $(length(b))" : worst)
    end
end
println("\n=== ", npass, " MATCH, ", ndiff, " DIFF, ", nerr, " ERR  of ", length(names), " (atol=$atol rtol=$rtol) ===")
end

run_sweep(ARGS[1], ARGS[2],
          length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 1.0,
          length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 0.01)
