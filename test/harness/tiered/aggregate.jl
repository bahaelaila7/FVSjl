# aggregate.jl — AGGREGATE STATISTICS for the fast tier. Print-only: it never gates anything (the per-cell
# allowlist in test_tiered.jl remains the gate). The per-cell view answers "is anything new broken?"; these
# aggregates answer "how close are we, where is the mass, and is it moving?" — the actionable view.
#
# Denominators come from the committed goldens themselves (rows x columns of every DBS golden + every printed
# field of every .sum row), i.e. exactly the cells the fast tier compares.

_agg_num(x) = tryparse(Float64, strip(x))

"Compared-cell denominators per (variant, regime) and per (variant, file), from the committed goldens."
function golden_cell_counts(vs::Vector{String})
    per_vr = Dict{Tuple{String,String},Int}(); per_vf = Dict{Tuple{String,String},Int}()
    for v in vs
        dir = fixture_dir(v)
        isdir(dir) || continue
        for f in readdir(dir)
            path = joinpath(dir, f)
            if endswith(f, ".csv")
                stem = f[1:end-4]; parts = split(stem, '.')
                length(parts) >= 2 || continue
                regime = String(last(split(parts[1], '_'))); file = String(parts[2])
                ls = readlines(path); length(ls) >= 2 || continue
                cells = (length(ls) - 1) * (count(==(','), ls[1]) + 1)
                per_vr[(v, regime)] = get(per_vr, (v, regime), 0) + cells
                per_vf[(v, file)]   = get(per_vf, (v, file), 0) + cells
            elseif endswith(f, ".live.sum")
                regime = String(last(split(f[1:end-9], '_')))
                cells = 0
                for l in readlines(path)
                    occursin("-999", l) && continue
                    cells += length(split(l))
                end
                per_vr[(v, regime)] = get(per_vr, (v, regime), 0) + cells
                per_vf[(v, "sum")]  = get(per_vf, (v, "sum"), 0) + cells
            end
        end
    end
    (per_vr, per_vf)
end

"|Δ| of a mismatch when both sides parse as numbers, else `nothing` (text/presence difference)."
function _absdiff(m)
    g = _agg_num(m.gold); j = _agg_num(m.got)
    (g === nothing || j === nothing) ? nothing : abs(g - j)
end

"""
    print_aggregate_stats(all_ms, vs)

Print the aggregate view: per-variant match rate, magnitude distribution, per-variant×regime rates, and the
columns carrying the most mismatching cells (with the worst |Δ| and an example). Print-only.
"""
function print_aggregate_stats(all_ms, vs::Vector{String})
    per_vr, per_vf = golden_cell_counts(vs)
    tot_v = Dict{String,Int}(); for ((v, _), n) in per_vr; tot_v[v] = get(tot_v, v, 0) + n; end
    mm_v = Dict{String,Int}(); mm_vr = Dict{Tuple{String,String},Int}(); mm_vfc = Dict{Tuple{String,String,String},Int}()
    bands = Dict{String,Dict{String,Int}}(); stands = Dict{String,Set{String}}()
    worst = Dict{Tuple{String,String,String},Tuple{Float64,Any}}()
    for m in all_ms
        m.file in ("TALLY", "FIXTURE") && continue
        v = m.variant
        mm_v[v] = get(mm_v, v, 0) + 1
        mm_vr[(v, m.regime)] = get(mm_vr, (v, m.regime), 0) + 1
        k = (v, m.file, m.col); mm_vfc[k] = get(mm_vfc, k, 0) + 1
        push!(get!(stands, v, Set{String}()), m.stand)
        d = _absdiff(m); b = get!(bands, v, Dict{String,Int}())
        lbl = d === nothing ? "text" : d <= 1e-4 ? "<=1e-4" : d <= 1.0 ? "<=1" : ">1"
        b[lbl] = get(b, lbl, 0) + 1
        if d !== nothing
            w = get(worst, k, (-1.0, nothing)); d > w[1] && (worst[k] = (d, m))
        end
    end
    println("\n== aggregate statistics (print-only; the per-cell allowlist above is the gate) ==")
    println(rpad("variant", 9), lpad("compared", 11), lpad("mismatch", 10), lpad("match%", 9),
            lpad("stands", 8), "   magnitude of mismatches")
    for v in sort(vs)
        tot = get(tot_v, v, 0); mm = get(mm_v, v, 0); tot == 0 && continue
        b = get(bands, v, Dict{String,Int}())
        mag = join(["$k=$(get(b,k,0))" for k in ("text", "<=1e-4", "<=1", ">1") if get(b, k, 0) > 0], " ")
        println(rpad(v, 9), lpad(tot, 11), lpad(mm, 10),
                lpad(string(round(100 * (tot - mm) / tot; digits = 2)), 9),
                lpad(length(get(stands, v, Set{String}())), 8), "   ", mag)
    end
    tt = sum(values(tot_v)); tm = sum(values(mm_v))
    tt > 0 && println(rpad("ALL", 9), lpad(tt, 11), lpad(tm, 10), lpad(string(round(100 * (tt - tm) / tt; digits = 2)), 9))
    println("\n-- match% per variant/regime --")
    for k in sort(collect(keys(per_vr)))
        n = per_vr[k]; n == 0 && continue
        mm = get(mm_vr, k, 0)
        println("  ", rpad("$(k[1])/$(k[2])", 22), lpad(n, 10), lpad(mm, 9),
                lpad(string(round(100 * (n - mm) / n; digits = 2)), 9), "%")
    end
    println("\n-- columns carrying the most mismatching cells --")
    for (k, n) in first(sort(collect(mm_vfc); by = x -> -x[2]), 15)
        w = get(worst, k, (0.0, nothing)); m = w[2]
        ex = m === nothing ? "" : "  e.g. $(m.stand)@$(m.year) gold=$(m.gold) jl=$(m.got)"
        println("  ", rpad("$(k[1])/$(k[2])/$(k[3])", 46), lpad(n, 8),
                m === nothing ? "" : "  max|Δ|=$(round(w[1]; sigdigits = 4))", ex)
    end
    return nothing
end
