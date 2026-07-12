# annotate_manifest.jl — Pillar-1: annotate the existing per-variant sample manifests with their STRATA
# (forest type, stand-structure/age class, site class, geography) and emit a coverage report. Read-only on
# the FVS-ready DB; does NOT touch the CN\tVARIANT `_sample.txt` files the sweep/validate harnesses consume.
#
# Produces (per variant, next to the sample): <v>_sample_strata.csv  (STAND_CN,forest_type,age_class,
# site_class,ecoregion) and prints a distinct-strata coverage summary — the documented "plot IDs + strata"
# manifest the Pillar-1 done-state calls for.
#
# Usage: julia --project=. test/harness/fia/annotate_manifest.jl [SN NE CS LS]

import SQLite, DBInterface

const DB = "/workspace/SQLite_FIADB_ENTIRE.db"
const HERE = @__DIR__

# stand-structure class from stand age (seedling→sawtimber); FIA-conventional breaks
age_class(a) = a === missing || a == 0 ? "unknown" :
               a < 20 ? "1_seedsap" : a < 40 ? "2_pole" : a < 80 ? "3_smallsaw" :
               a < 120 ? "4_largesaw" : "5_oldgrowth"
# site class from site index (feet); coarse low/med/high/vhigh
site_class(s) = s === missing || s == 0 ? "unknown" :
                s < 40 ? "lo" : s < 60 ? "med" : s < 80 ? "hi" : "vhi"

function annotate(variant)
    smp = joinpath(HERE, lowercase(variant) * "_sample.txt")
    isfile(smp) || (println(stderr, "no manifest for $variant"); return nothing)
    cns = String[]
    for ln in eachline(smp)
        f = split(strip(ln), '\t'); isempty(f[1]) || push!(cns, String(f[1]))
    end
    isempty(cns) && return nothing
    db = SQLite.DB(DB)
    inlist = join(("'" * c * "'" for c in cns), ",")
    q = """SELECT STAND_CN, FOREST_TYPE, FOREST_TYPE_FIA, AGE, SITE_INDEX, ECOREGION
           FROM FVS_STANDINIT_COND WHERE STAND_CN IN ($inlist)"""
    rows = NamedTuple[]
    for r in DBInterface.execute(db, q)
        ft = r.FOREST_TYPE_FIA !== missing ? string(r.FOREST_TYPE_FIA) :
             (r.FOREST_TYPE !== missing ? string(r.FOREST_TYPE) : "unknown")
        push!(rows, (cn = string(r.STAND_CN), ft = ft,
                     ac = age_class(r.AGE === missing ? missing : Int(round(Float64(r.AGE)))),
                     sc = site_class(r.SITE_INDEX === missing ? missing : Float64(r.SITE_INDEX)),
                     eco = r.ECOREGION === missing ? "unknown" : string(r.ECOREGION)))
    end
    out = joinpath(HERE, lowercase(variant) * "_sample_strata.csv")
    open(out, "w") do io
        println(io, "STAND_CN,forest_type,age_class,site_class,ecoregion")
        for x in rows; println(io, x.cn, ',', x.ft, ',', x.ac, ',', x.sc, ',', x.eco); end
    end
    nft = length(unique(getfield.(rows, :ft)))
    nac = length(unique(getfield.(rows, :ac)))
    nsc = length(unique(getfield.(rows, :sc)))
    neco = length(unique(getfield.(rows, :eco)))
    println("$variant: n=$(length(rows))  forest_types=$nft  age_classes=$nac  site_classes=$nsc  ecoregions=$neco  → $(basename(out))")
    println("   age spread: ", join(sort(collect(Set(getfield.(rows, :ac)))), " "))
    println("   site spread: ", join(sort(collect(Set(getfield.(rows, :sc)))), " "))
    return (variant=variant, n=length(rows), nft=nft, nac=nac, nsc=nsc, neco=neco)
end

vs = length(ARGS) > 0 ? ARGS : ["SN", "NE", "CS", "LS"]
println("=== Pillar-1 strata coverage of the existing sample manifests ===")
for v in vs; annotate(v); end
