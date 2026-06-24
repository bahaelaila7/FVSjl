# test_dbs_treelist.jl — C6 DBS: the FVS_TreeList SQLite table (dbstrls.f).
#
# DATABASE TREELIDB makes FVSjl write a per-cycle, per-tree FVS_TreeList table. The exact
# per-tree row set can differ from Fortran by the tripling/COMCUP record PARTITION (same total,
# different #records), so this validates DATA INTEGRITY: each cycle's Σ(TPA) and the volume-
# weighted Σ(TCuFt·TPA) must reconstruct the stand TPA / cubic volume from the (Fortran-bit-
# exact) text `.sum`. Confirms the per-tree TPA (PROB/GROSPC) + per-tree volumes are correct.

using Test, FVSjl
using SQLite

const _TL_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")

@testset "DBS FVS_TreeList table — per-tree aggregates vs .sum" begin
    key = joinpath(_TL_DIR, "dbs_treelist.key")
    if !isfile(key)
        @test_skip "dbs_treelist scenario not available"
    else
        dbpath = tempname() * ".db"; isfile(dbpath) && rm(dbpath)
        tmpkey = joinpath(_TL_DIR, "_tl_run.key")
        cp(joinpath(_TL_DIR, "dbs_treelist.tre"), joinpath(_TL_DIR, "_tl_run.tre"); force=true)
        open(tmpkey, "w") do io
            for l in eachline(key); println(io, replace(l, "__DSNOUT__" => dbpath)); end
        end
        try
            sumtxt = FVSjl.run_keyfile(tmpkey; faithful = true)
            @test isfile(dbpath)
            # per-cycle TreeList aggregates: Σ TPA and Σ(TCuFt·TPA)
            db = SQLite.DB(dbpath); agg = Dict{Int,Vector{Float64}}()
            for r in SQLite.DBInterface.execute(db, "SELECT Year,TPA,TCuFt FROM FVS_TreeList")
                v = get!(agg, r.Year, [0.0, 0.0])
                v[1] += Float64(r.TPA); v[2] += Float64(r.TCuFt) * Float64(r.TPA)
            end
            SQLite.close(db)
            @test !isempty(agg)

            # the text .sum gives the validated stand TPA (col 3) + total cuft (col 9) per year
            for l in split(sumtxt, "\n")
                f = split(l)
                length(f) >= 9 || continue
                y = tryparse(Int, f[1]); (y === nothing && continue)
                haskey(agg, y) || continue
                @test abs(agg[y][1] - parse(Int, f[3])) <= 1          # Σ TPA ≈ stand TPA
                @test abs(agg[y][2] - parse(Int, f[9])) <= 3          # Σ(TCuFt·TPA) ≈ total cuft
            end
        finally
            rm(tmpkey; force=true); rm(joinpath(_TL_DIR, "_tl_run.tre"); force=true)
            rm(dbpath; force=true)
        end
    end
end
