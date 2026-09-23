# extract_sample.jl — Pillar-1: build a REPRODUCIBLE, stratified per-variant FIA sample from the FVS-ready
# tables (FVS_STANDINIT_COND, VARIANT column). Deterministic (no RNG): stands are ordered by the stratum
# key then STAND_CN, and every K-th is taken so the sample spreads evenly across strata.
#
# Strata: ECOREGION (ecological unit — 100% populated; drives the SN dgf EUT term + species/geography),
# LOCATION (national forest), and the stand's INVENTORY TREE-COUNT CLASS. Sampling across ECOREGION
# guarantees the sample exercises many EUT DG coefficients + species mixes — exactly the axis that
# surfaced the eco_unit bug.
#
# ⚠ The tree-count class is NOT cosmetic. Many western variants' FIA populations are dominated by
# stands with ZERO FVS_TREEINIT_COND records (nonstocked / seedling-only conditions that FVS fills by
# AUTOES regeneration). Without this stratum an EM draw came out 37/40 BARE, and the 12-stand EM
# TIERED FIXTURE came out 12/12 bare — so both instruments measured only the establishment path and
# were blind to the growth/crown/mortality models they were built to certify. Every sample now spans
# all four classes present in the population, so a change to either path is visible.
#
# Usage:  julia --project=. test/harness/fia/extract_sample.jl <VARIANT> <N> [out.txt]
#   VARIANT ∈ {SN,NE,CS,LS,...} (the FVS_STANDINIT_COND.VARIANT value)
#   writes <N> lines "STAND_CN<TAB>VARIANT" to out.txt (default test/harness/fia/<variant>_sample.txt),
#   plus a strata summary to stderr. Read-only on the DB.
#   TREE_CLASS_STRATA=0 restores the pre-2026-09-23 ECOREGION/LOCATION-only stride (for reproducing an
#   older sample); the class boundaries are fixed so a given (variant, N) is still deterministic.

import SQLite, DBInterface

const DB = "/workspace/SQLite_FIADB_ENTIRE.db"

# Inventory tree-count class. Fixed boundaries so a (variant, N) request stays reproducible.
#   1 = 0 records (bare — FVS regenerates it), 2 = 1-9, 3 = 10-49, 4 = 50+
_tree_class(nt::Integer) = nt == 0 ? 1 : nt < 10 ? 2 : nt < 50 ? 3 : 4

# Even stride over a strata-ordered list: takes k items spread across the whole list.
_stride_take(rows, k) =
    isempty(rows) || k <= 0 ? eltype(rows)[] :
    rows[unique(clamp.(round.(Int, (0.5:1.0:min(k, length(rows))) .* (length(rows) / min(k, length(rows)))),
                       1, length(rows)))]

function extract(variant::AbstractString, n::Int, out::AbstractString)
    db = SQLite.DB(DB)
    # Pull all candidate stands for the variant with their inventory tree count, ordered by
    # (ECOREGION, LOCATION, STAND_CN) — deterministic. LEFT JOIN so the ZERO-tree stands are kept
    # (they are a legitimate, and in several variants dominant, part of the population).
    rows = NamedTuple[]
    q = """SELECT s.STAND_CN cn, s.ECOREGION eco, s.LOCATION loc, COUNT(t.STAND_CN) nt
           FROM FVS_STANDINIT_COND s LEFT JOIN FVS_TREEINIT_COND t ON t.STAND_CN = s.STAND_CN
           WHERE s.VARIANT = '$(variant)' AND s.STAND_CN IS NOT NULL
           GROUP BY s.STAND_CN
           ORDER BY s.ECOREGION, s.LOCATION, s.STAND_CN"""
    for r in DBInterface.execute(db, q)
        push!(rows, (cn = string(r.cn), eco = string(something(r.eco, "")),
                     loc = r.loc === missing ? 0 : Int(r.loc), nt = Int(r.nt)))
    end
    total = length(rows)
    total == 0 && error("no stands for VARIANT=$variant")
    n = min(n, total)
    if get(ENV, "TREE_CLASS_STRATA", "1") == "0"
        sample = _stride_take(rows, n)                       # legacy ECOREGION/LOCATION-only stride
    else
        # BALANCED (not proportional) allocation across the tree-count classes present in the population.
        # These samples are CODE-PATH COVERAGE instruments, not population estimates: EM's FIA population
        # is ~88% bare, so a proportional draw still spends 35 of 40 stands on the establishment path and
        # leaves the growth/crown/mortality models nearly unexercised. Equal shares (remainder to the
        # largest classes, shortfall from a small class redistributed) keep every path visible.
        # ⚠ Because of this, a pass-rate over one of these samples is NOT a population pass-rate.
        byclass = Dict{Int,Vector{eltype(rows)}}()
        for r in rows; push!(get!(byclass, _tree_class(r.nt), eltype(rows)[]), r); end
        present = sort(collect(keys(byclass)))
        quota = Dict(c => 0 for c in present)
        left = n
        # Round-robin from the smallest class up, so a class with few stands is filled before the big
        # ones absorb the remainder; stop handing a class more than it has.
        while left > 0 && any(quota[c] < length(byclass[c]) for c in present)
            for c in sort(present; by = c -> length(byclass[c]))
                left <= 0 && break
                quota[c] < length(byclass[c]) || continue
                quota[c] += 1; left -= 1
            end
        end
        sample = eltype(rows)[]
        for c in present
            append!(sample, _stride_take(byclass[c], min(quota[c], length(byclass[c]))))
        end
        sort!(sample; by = r -> (r.eco, r.loc, r.cn))          # restore the deterministic strata order
    end
    open(out, "w") do io
        for s in sample
            println(io, s.cn, '\t', variant)
        end
    end
    necos = length(unique(getfield.(sample, :eco)))
    nlocs = length(unique(getfield.(sample, :loc)))
    cls = [count(r -> _tree_class(r.nt) == c, sample) for c in 1:4]
    println(stderr, "VARIANT=$variant  population=$total  sampled=$(length(sample))  " *
                    "distinct ECOREGION=$necos  distinct LOCATION=$nlocs  " *
                    "tree-count classes [0 | 1-9 | 10-49 | 50+] = $cls  → $out")
    return length(sample)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("usage: extract_sample.jl <VARIANT> <N> [out.txt]")
    v = ARGS[1]; n = parse(Int, ARGS[2])
    out = length(ARGS) >= 3 ? ARGS[3] : joinpath(@__DIR__, lowercase(v) * "_sample.txt")
    extract(v, n, out)
end
