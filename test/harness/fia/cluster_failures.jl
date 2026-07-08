# cluster_failures.jl — triage a sweep's failing-stand list by CAUSE proxy so systematic bug CLASSES pop out
# (instead of deep-diving individual stands). Groups failures by dominant species + a structure bucket
# (stand avg DBH → seedling/sapling/pole/sawtimber), against the same sub-DB the sweep ran on.
#
# Usage: julia --project=. test/harness/fia/cluster_failures.jl <failout.txt> <subdb.db> [allstands.txt]
#   <failout.txt>  = FIA_FAILOUT list (STAND_CN<TAB>VARIANT) from the sweep
#   <subdb.db>     = the indexed sub-DB used for the sweep
#   [allstands.txt]= optional full stand list (.stands) → also reports each species' FAIL/TOTAL rate

import SQLite, DBInterface

struct_bucket(d) = d < 1 ? "seed<1" : d < 5 ? "sap1-5" : d < 9 ? "pole5-9" : d < 15 ? "saw9-15" : "lg15+"

function domspecies(db, cn)
    r = first(DBInterface.execute(db, "SELECT SPECIES sp, COUNT(*) n FROM FVS_TREEINIT_COND WHERE STAND_CN='$cn' GROUP BY SPECIES ORDER BY n DESC LIMIT 1"))
    a = first(DBInterface.execute(db, "SELECT ROUND(AVG(DIAMETER),1) d, COUNT(*) c FROM FVS_TREEINIT_COND WHERE STAND_CN='$cn'"))
    return (sp = Int(r.sp), frac = r.n / a.c, avgdbh = a.d === missing ? 0.0 : Float64(a.d))
end

function main(failout, subdb, allstands)
    db = SQLite.DB(subdb)
    fails = [split(strip(l), '\t')[1] for l in eachline(failout) if !isempty(strip(l))]
    # cluster fails by (dominant species, structure bucket)
    byspec = Dict{Int,Int}(); bystruct = Dict{String,Int}(); byboth = Dict{Tuple{Int,String},Int}()
    for cn in fails
        d = domspecies(db, cn)
        byspec[d.sp] = get(byspec, d.sp, 0) + 1
        b = struct_bucket(d.avgdbh); bystruct[b] = get(bystruct, b, 0) + 1
        byboth[(d.sp, b)] = get(byboth, (d.sp, b), 0) + 1
    end
    println("=== $(length(fails)) failures — by dominant FIA species (top 15) ===")
    for (sp, n) in sort(collect(byspec), by = x -> -x[2])[1:min(15, end)]
        println("  FIA $sp : $n failures")
    end
    println("=== by structure (stand avg DBH) ===")
    for (b, n) in sort(collect(bystruct), by = x -> -x[2]); println("  $b : $n"); end
    # if the full stand list is given, compute per-species FAIL RATE (fails / total for that dom species)
    if allstands !== nothing && isfile(allstands)
        tot = Dict{Int,Int}()
        for l in eachline(allstands)
            isempty(strip(l)) && continue
            cn = split(strip(l), '\t')[1]
            sp = domspecies(db, cn).sp; tot[sp] = get(tot, sp, 0) + 1
        end
        println("=== fail RATE by dominant species (fails/total, species with ≥3 fails) ===")
        for (sp, n) in sort(collect(byspec), by = x -> -x[2])
            n < 3 && continue
            t = get(tot, sp, n)
            println("  FIA $sp : $n / $t = $(round(100n/t, digits=1))%")
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("usage: cluster_failures.jl <failout.txt> <subdb.db> [allstands.txt]")
    main(ARGS[1], ARGS[2], length(ARGS) >= 3 ? ARGS[3] : nothing)
end
