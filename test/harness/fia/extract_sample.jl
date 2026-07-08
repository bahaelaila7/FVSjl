# extract_sample.jl — Pillar-1: build a REPRODUCIBLE, stratified per-variant FIA sample from the FVS-ready
# tables (FVS_STANDINIT_COND, VARIANT column). Deterministic (no RNG): stands are ordered by the stratum
# key then STAND_CN, and every K-th is taken so the sample spreads evenly across strata.
#
# Strata: ECOREGION (ecological unit — 100% populated; drives the SN dgf EUT term + species/geography) and
# LOCATION (national forest). Sampling across ECOREGION guarantees the sample exercises many EUT DG
# coefficients + species mixes — exactly the axis that surfaced the eco_unit bug.
#
# Usage:  julia --project=. test/harness/fia/extract_sample.jl <VARIANT> <N> [out.txt]
#   VARIANT ∈ {SN,NE,CS,LS,...} (the FVS_STANDINIT_COND.VARIANT value)
#   writes <N> lines "STAND_CN<TAB>VARIANT" to out.txt (default test/harness/fia/<variant>_sample.txt),
#   plus a strata summary to stderr. Read-only on the DB.

import SQLite, DBInterface

const DB = "/workspace/SQLite_FIADB_ENTIRE.db"

function extract(variant::AbstractString, n::Int, out::AbstractString)
    db = SQLite.DB(DB)
    # Pull all candidate stands for the variant, ordered by (ECOREGION, LOCATION, STAND_CN) — deterministic.
    rows = NamedTuple[]
    q = """SELECT STAND_CN, ECOREGION, LOCATION FROM FVS_STANDINIT_COND
           WHERE VARIANT = '$(variant)' AND STAND_CN IS NOT NULL
           ORDER BY ECOREGION, LOCATION, STAND_CN"""
    for r in DBInterface.execute(db, q)
        push!(rows, (cn = string(r.STAND_CN), eco = string(something(r.ECOREGION, "")),
                     loc = r.LOCATION === missing ? 0 : Int(r.LOCATION)))
    end
    total = length(rows)
    total == 0 && error("no stands for VARIANT=$variant")
    n = min(n, total)
    # Even stride across the strata-ordered list → spreads the sample across ECOREGION/LOCATION.
    stride = total / n
    idx = unique(clamp.(round.(Int, (0.5:1.0:n) .* stride), 1, total))
    sample = rows[idx]
    open(out, "w") do io
        for s in sample
            println(io, s.cn, '\t', variant)
        end
    end
    necos = length(unique(getfield.(sample, :eco)))
    nlocs = length(unique(getfield.(sample, :loc)))
    println(stderr, "VARIANT=$variant  population=$total  sampled=$(length(sample))  " *
                    "distinct ECOREGION=$necos  distinct LOCATION=$nlocs  → $out")
    return length(sample)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("usage: extract_sample.jl <VARIANT> <N> [out.txt]")
    v = ARGS[1]; n = parse(Int, ARGS[2])
    out = length(ARGS) >= 3 ? ARGS[3] : joinpath(@__DIR__, lowercase(v) * "_sample.txt")
    extract(v, n, out)
end
