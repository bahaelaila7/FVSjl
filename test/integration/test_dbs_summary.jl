# test_dbs_summary.jl — C6 DBS: the FVS_Summary SQLite table (dbssumry.f) vs live Fortran.
#
# The DATABASE block (DSNOUT file + SUMMARY) makes FVSjl write its per-cycle summary into a
# SQLite FVS_Summary table — the same data as the text `.sum`, into a database (the "modern
# IO / same SQLite outputs" goal). This runs the scenario (rewriting the placeholder DSNOUT to
# a temp .db), reads FVS_Summary back, and compares to the Fortran FVSOut.db dump baseline.

using Test, FVSjl
using SQLite

const _DB_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")

function _read_fvs_summary(dbpath)
    db = SQLite.DB(dbpath); rows = NTuple{9,Any}[]
    for r in SQLite.DBInterface.execute(db,
            "SELECT Year,Tpa,BA,SDI,CCF,TopHt,QMD,MCuFt,BdFt FROM FVS_Summary ORDER BY Year")
        push!(rows, (r.Year, r.Tpa, r.BA, r.SDI, r.CCF, r.TopHt,
                     round(Float64(r.QMD), digits=1), r.MCuFt, r.BdFt))
    end
    SQLite.close(db)
    return rows
end

@testset "DBS FVS_Summary table vs Fortran" begin
    key = joinpath(_DB_DIR, "dbs_summary.key")
    sav = joinpath(_DB_DIR, "dbs_summary.csv.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "dbs_summary scenario not available"
    else
        # point DSNOUT at a temp .db + keep the .tre basename
        dbpath = tempname() * ".db"; isfile(dbpath) && rm(dbpath)
        tmpkey = joinpath(_DB_DIR, "_dbs_run.key")
        cp(joinpath(_DB_DIR, "dbs_summary.tre"), joinpath(_DB_DIR, "_dbs_run.tre"); force=true)
        open(tmpkey, "w") do io
            for l in eachline(key); println(io, replace(l, "__DSNOUT__" => dbpath)); end
        end
        try
            FVSjl.run_keyfile(tmpkey; faithful = true)
            @test isfile(dbpath)                              # the SQLite db was written
            jl = _read_fvs_summary(dbpath)

            # Fortran FVS_Summary baseline (Year,Tpa,BA,SDI,CCF,TopHt,QMD,MCuFt,BdFt)
            ft = NTuple{9,Any}[]
            for l in Iterators.drop(eachline(sav), 1)
                c = split(strip(l), ',')
                isempty(c[1]) && continue
                push!(ft, (parse(Int,c[1]), parse(Int,c[2]), parse(Int,c[3]), parse(Int,c[4]),
                           parse(Int,c[5]), parse(Int,c[6]), parse(Float64,c[7]),
                           parse(Int,c[8]), parse(Int,c[9])))
            end

            @test length(jl) == length(ft)
            for (j, f) in zip(jl, ft)
                @test j[1] == f[1]                            # Year
                for k in 2:6; @test j[k] == f[k]; end         # Tpa,BA,SDI,CCF,TopHt
                @test abs(j[7] - f[7]) <= 0.05                # QMD
                @test j[8] == f[8]                            # MCuFt — BIT-EXACT (measured Δ=0; was over-cautious ≤2)
                @test j[9] == f[9]                            # BdFt  — BIT-EXACT (measured Δ=0; was over-cautious ≤5, closed by the BFTOPK fix)
            end
        finally
            rm(tmpkey; force=true); rm(joinpath(_DB_DIR, "_dbs_run.tre"); force=true)
            rm(dbpath; force=true)
        end
    end
end
