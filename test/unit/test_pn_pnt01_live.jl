# PN pnt01 (S248112, 5 stands: control, THINDBH, shelterwood+ECON, FFE SIMFIRE+SALVAGE, NOTREES+PLANT) vs live
# FVSpn_g16. pnt01 used to stop in establishment (`KeyError :estab_min_ht`). FVSpn compiles WC's estab/esgent/
# regent/smhgdg/dgbnd/cratet/morts sources, so PN runs the shared WC paths with its own data (smhgdg.f:260 and
# morts.f:337 raw DF SITEAR, pn/htdbh.f, misintpn.f). Checked: every .sum row's non-volume columns (TPA..QMD) equal
# live; the volume columns are PN's volume path (separate work). Stand 4 from 2013 (post SIMFIRE+SALVAGE) is open.
using Test, FVSjl

function _pn_sum_rows(lines)
    rows = Dict{Tuple{Int,Int},Vector{SubString{String}}}(); sid = 0
    for l in lines
        startswith(l, "-999") && (sid += 1; continue)
        f = split(l)
        (sid > 0 && length(f) > 10 && length(f[1]) == 4 && all(isdigit, f[1])) || continue
        rows[(sid, parse(Int, f[1]))] = f
    end
    return rows
end

@testset "PN pnt01 non-volume .sum == live FVSpn" begin
    dir = joinpath(@__DIR__, "..", "fixtures", "pacificnorthwest", "pnt01live")
    gold = _pn_sum_rows(readlines(joinpath(dir, "pnt01.live.sum")))
    out = mktempdir()
    cp(joinpath(dir, "pnt01.key"), joinpath(out, "pnt01.key")); cp(joinpath(dir, "pnt01.tre"), joinpath(out, "pnt01.tre"))
    jl = cd(out) do
        _pn_sum_rows(split(FVSjl.run_keyfile("pnt01.key"; variant = FVSjl.PacificNorthwest(), output = :sum), '\n'))
    end
    @test length(gold) == 56
    @test sort(collect(keys(jl))) == sort(collect(keys(gold)))
    for k in sort(collect(keys(gold)))
        g = gold[k][3:8]; j = get(jl, k, fill(SubString(""), 12))[3:8]   # TPA BA SDI CCF TopHt QMD
        if k[1] == 4 && 2013 <= k[2] <= 2033
            # FFE SIMFIRE+FLAMEADJ+SALVAGE stand after the 2003 fire: with the R6 snag/fuel layer (ffe-westside) the fire
            # and FVS_Mortality match live, and TPA is within 1 (2013 70 vs 71); 2043 is exact again.
            @test_broken j == g
        else
            @test j == g
        end
    end
end
