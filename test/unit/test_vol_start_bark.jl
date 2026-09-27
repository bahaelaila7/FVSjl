# test_vol_start_bark.jl — projected-cycle volumes use the START-of-cycle bark ratio (vols.f:132,150-151), vs live.
#
# Every western vols.f does D=DBH(I); BARK=BRATIO(ISPC,D,H); IF(.NOT.LSTART) D=D+DG(I)/BARK — the bark that sets the
# merch/board tops (MTOPS=TOPD·BARK, fvsvol.f:172), the FW2 DBTBH and CFTOPK is the pre-growth one. jl recomputed it
# from the grown DBH in UT/CI/CR/EM/IE/KT/EC/PN (BM/SO/WS/NC/CA/OP already used the stashed t.vol_bark).
# Fixture: two FIA stands (sub-DB of SQLite_FIADB_ENTIRE) with live .sum from FVScr_clean / FVSut_g16:
#   CR 316922874489998 (pure PP, 300FW2W122): 2024 TCuFt/MCuFt/BdFt were 663/512/2368 vs live 658/490/2299.
#   UT 2875463010690 (402MATW122 PP): 2004 MCuFt 502 vs live 504 (PP D9.90, D_start 8.30: merch 5.1 vs 5.4).
#   (UT keyed to 2 cycles: master UT growth diverges from 2024 — TopHt 35 vs 34 — a separate, non-volume residual.)
using FVSjl, Test

const _NVX = joinpath(@__DIR__, "..", "fixtures", "nvel_volume")

function _vsb_sum_rows(path::AbstractString)
    rows = Dict{Int,Vector{Float64}}()
    for ln in eachline(path)
        f = split(ln)
        (length(f) > 12 && all(isdigit, f[1]) && length(f[1]) == 4) || continue
        rows[parse(Int, f[1])] = parse.(Float64, f[3:12])  # TPA BA SDI CCF TopHt QMD TCuFt MCuFt SCuFt BdFt
    end
    return rows
end
_sum_rows_text(txt::AbstractString) = (p = tempname(); write(p, txt); _vsb_sum_rows(p))

@testset "projected-cycle volume uses the start-of-cycle bark (vs live)" begin
    for (v, key) in ((FVSjl.CentralRockies(), "cr_316922874489998"), (FVSjl.Utah(), "ut_2875463010690"))
        dir = mktempdir()
        cp(joinpath(_NVX, key * ".key"), joinpath(dir, key * ".key"))
        cp(joinpath(_NVX, "stands.db"), joinpath(dir, "stands.db"))
        jl = _sum_rows_text(cd(() -> FVSjl.run_keyfile(key * ".key"; variant = v, output = :sum), dir))
        live = _vsb_sum_rows(joinpath(_NVX, key * ".live.sum"))
        @test !isempty(live) && keys(jl) == keys(live)
        for y in sort(collect(keys(live)))
            @test (key, y, jl[y]) == (key, y, live[y])      # every printed column, exact
        end
    end
end
