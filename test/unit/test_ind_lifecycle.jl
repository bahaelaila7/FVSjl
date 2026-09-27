# test_ind_lifecycle.jl — cycle-0 AVHT40/DENSE walk CRATET's IND, then gradd.f:186's fresh RDPSRT, in EVERY variant.
#
# Every variant's cratet.f seeds IND from IND1 and calls RDPSRT(ITRN,DBH,IND,.FALSE.) before the LBKDEN DENSE, and
# re-sorts with RDPSRT(.TRUE.) only when dead records were deleted (cr 139-146/250, ec 195-202/306, ...). jl walked that
# IND only for BM/CI/UT/TT/CA/SO/WS and used an empirical double sort elsewhere, which picks a different tree at the
# 40-TPA cutoff when DBHs tie. Fixture: two FIA stands (sub-DB of SQLite_FIADB_ENTIRE), live .sum from FVScr_clean and
# FVSec_g16:
#   CR 12166988010690: 2007 TopHt was 29 vs live 30; the cycle-1 DG then read the wrong PCT ⇒ 1/7 rows exact, now 7/7.
#   EC 22404660010497: 2035/2045 rows diverged (5/7 exact), now 7/7.
using FVSjl, Test

const _INDLC = joinpath(@__DIR__, "..", "fixtures", "ind_lifecycle")

function _indlc_sum_rows(path::AbstractString)
    rows = Dict{Int,Vector{Float64}}()
    for ln in eachline(path)
        f = split(ln)
        (length(f) > 12 && all(isdigit, f[1]) && length(f[1]) == 4) || continue
        rows[parse(Int, f[1])] = parse.(Float64, f[3:12])  # TPA BA SDI CCF TopHt QMD TCuFt MCuFt SCuFt BdFt
    end
    return rows
end

@testset "cycle-0 CRATET IND lifecycle in every variant (vs live)" begin
    for (v, key) in ((FVSjl.CentralRockies(), "cr_12166988010690"), (FVSjl.EastCascades(), "ec_22404660010497"))
        dir = mktempdir()
        cp(joinpath(_INDLC, key * ".key"), joinpath(dir, key * ".key"))
        cp(joinpath(_INDLC, "stands.db"), joinpath(dir, "stands.db"))
        p = tempname()
        write(p, cd(() -> FVSjl.run_keyfile(key * ".key"; variant = v, output = :sum), dir))
        jl = _indlc_sum_rows(p)
        live = _indlc_sum_rows(joinpath(_INDLC, key * ".live.sum"))
        @test !isempty(live) && keys(jl) == keys(live)
        for y in sort(collect(keys(live)))
            @test (key, y, get(jl, y, nothing)) == (key, y, live[y])      # every printed column, exact
        end
    end
end
