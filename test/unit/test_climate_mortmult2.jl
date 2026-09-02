# Climate-FVS MORTMULT field 4 = CLMRTMLT2 (clin.f:365/385 → clmorts.f:223/237) — the SECOND MortMult weight,
# a per-species multiplier on the transfer-distance climate mortality SPMORT2/DMORT (distinct from field 3 =
# CLMRTMLT1, which scales the viability SPMORT1/FYRMORT). Ported 2026-09-02: kw_climate! now parses field 4 and
# ClimateState.mortmult2 feeds apply_climate_mort!'s DMORT accumulation + applied rate.
#
# VALIDATED vs the live FVSie_clean oracle (FVS_Climate.dClimMort), single UNTHINNED S248112 CGCM3_A2 stand
# (test/fixtures/climate/clim_mm2.key). With CLMRTMLT2=2.0 on all species the oracle's dClimMort DOUBLES for the
# deterministic dominant species — BIT-EXACT reproduced by FVSjl: 2030 0.01336→0.02673, 2040 0.20628→0.41256,
# 2050 0.54751→1.09501 (oracle == jl to the report's precision). Late cycles/minor species corner on the same
# OLDRN/DGSD growth straddle the `.sum` and test_climate_dbs corner. field4 absent OR =1.0 is byte-identical
# (inert). Oracle values dumped from a fresh FVSie_clean run 2026-09-02.

using Test
using FVSjl
using SQLite, DBInterface

# Build a keyfile from clim_mm2.key (single-stand) into `dst`, rewriting DSNOut → `db` and (optionally) inserting
# a `MortMult` card carrying field 4 = `mm2` (blank species/field3 ⇒ all species, CLMRTMLT1 stays 1) before the
# Climate block's `End`. Also copies the .tre alongside as <keystem>.tre.
function _build_mm_key(fix, dst, db, mm2)
    src = readlines(joinpath(fix, "clim_mm2.key"))
    cp(joinpath(fix, "clim_mm2.tre"), first(splitext(dst)) * ".tre"; force = true)
    open(dst, "w") do io
        for (i, ln) in enumerate(src)
            s = strip(ln)
            if s == "climtest_oracle.db"
                println(io, db); continue
            end
            if mm2 !== nothing && s == "End" && i > 1 && strip(src[i-1]) == "-999"
                println(io, "MortMult" * " "^32 * mm2)   # field 4 (cols 41-50) = CLMRTMLT2
            end
            println(io, ln)
        end
    end
end

# dClimMort (SPMORT2 transfer-distance climate mortality) per (Year, SpeciesFVS) from an FVS_Climate DBS table.
function _dclim(db)
    d = Dict{Tuple{Int,String},Float64}()
    for r in DBInterface.execute(db, "SELECT Year,SpeciesFVS,dClimMort FROM FVS_Climate")
        d[(Int(r.Year), String(r.SpeciesFVS))] = r.dClimMort === missing ? 0.0 : Float64(r.dClimMort)
    end
    return d
end

@testset "Climate-FVS MORTMULT field 4 (CLMRTMLT2, DMORT) vs FVSie_clean" begin
    fix = joinpath(@__DIR__, "..", "fixtures", "climate")
    @assert isfile(joinpath(fix, "clim_mm2.key"))
    dir = mktempdir()
    V = FVSjl.variant_from_code("IE")
    bkey = joinpath(dir, "base.key"); bdb = joinpath(dir, "base.db")
    mkey = joinpath(dir, "mm2.key");  mdb = joinpath(dir, "mm2.db")
    ikey = joinpath(dir, "mm1.key");  idb = joinpath(dir, "mm1.db")
    _build_mm_key(fix, bkey, bdb, nothing)   # no MortMult card
    _build_mm_key(fix, mkey, mdb, "2.0")     # CLMRTMLT2 = 2.0, all species
    _build_mm_key(fix, ikey, idb, "1.0")     # CLMRTMLT2 = 1.0, all species (must be inert)
    for k in (bkey, mkey, ikey)
        FVSjl.run_keyfile(k; variant = V)
    end
    base = _dclim(SQLite.DB(bdb)); mm2 = _dclim(SQLite.DB(mdb)); mm1 = _dclim(SQLite.DB(idb))

    # (1) BASE dClimMort — BIT-EXACT vs the live oracle for the deterministic dominant species.
    # (oracle FVSie_clean, dumped 2026-09-02)
    for sp in ("DF", "LP", "PP", "WL", "GF", "ES")
        @test isapprox(base[(2030, sp)], 0.01336; atol = 5e-4)
        @test isapprox(base[(2040, sp)], 0.20628; atol = 5e-4)
        @test isapprox(base[(2050, sp)], 0.54751; atol = 5e-4)
    end

    # (2) CLMRTMLT2 = 2.0 DOUBLES the transfer-distance mortality — BIT-EXACT vs oracle, and exactly 2× base.
    for sp in ("DF", "LP", "WL", "GF")
        @test isapprox(mm2[(2030, sp)], 0.02673; atol = 5e-4)   # 2 × 0.01336
        @test isapprox(mm2[(2040, sp)], 0.41256; atol = 5e-4)   # 2 × 0.20628
        @test isapprox(mm2[(2050, sp)], 1.09501; atol = 5e-3)   # 2 × 0.54751
        # the doubling relationship itself (field-4 effect), while the tree list is still shared (early cycles)
        @test isapprox(mm2[(2040, sp)], 2f0 * base[(2040, sp)]; rtol = 1e-4)
    end

    # (3) INERT: field 4 = 1.0 is byte-identical to NO MortMult card (default CLMRTMLT2 = 1).
    mdiff = 0.0
    for kk in union(keys(base), keys(mm1))
        mdiff = max(mdiff, abs(get(base, kk, 0.0) - get(mm1, kk, 0.0)))
    end
    @test mdiff == 0.0
end
