# test_west_ffe_crash_live.jl — western FFE (salvage/SIMFIRE) paths that crashed, vs the LIVE .sum rows.
#
# Keys and goldens: the tiered-west fixtures (test/fixtures/tiered/<v>/<stand>_<regime>.key / .live.sum, stands.db).
# Each case must run, and every summary row must equal the live oracle's.

using Test, FVSjl
const _FF = FVSjl
const _FF_TIER = joinpath(@__DIR__, "..", "fixtures", "tiered")
const _FF_VAR = Dict("ec" => _FF.EastCascades(), "nc" => _FF.Klamath(), "ca" => _FF.CentralCalifornia(),
                     "ws" => _FF.WestSierra())

function _ff_run(v, stem)
    dir = mktempdir()
    cp(joinpath(_FF_TIER, v, "stands.db"), joinpath(dir, "stands.db"))
    key = join([l == "out.db" ? joinpath(dir, "out.db") : l == "stands.db" ? joinpath(dir, "stands.db") : l
                for l in split(read(joinpath(_FF_TIER, v, stem * ".key"), String), '\n')], '\n')
    write(joinpath(dir, "k.key"), key)
    txt = ""; err = ""
    try
        txt = _FF.run_keyfile(joinpath(dir, "k.key"); variant = _FF_VAR[v])
    catch e
        err = sprint(showerror, e)
    end
    return txt, err
end

# (variant, stand_regime, what crashed)
const _FF_CASES = [
    ("ec", "30193987010497_salvage",  "fmcrowe.f DOBF = 4/BRATIO: EC had no bark in the FFE crown-biomass chain (KeyError :bark_intercept)"),
    ("ec", "504392203126144_simfire", "same"),
    ("nc", "449523860489998_salvage", "fmcroww.f SPIE group 21 (madrone) small/large-tree crown weight not ported; then nc_cwcalc 81802"),
    ("nc", "23721711010900_simfire",  "fmcroww.f SPIE group 21 large-tree; nc_cwcalc 81802/63102 (MA/TO)"),
    ("nc", "23660512010900_simfire",  "FMCBA/FMCFMD crown width = CRWDTH: R5CRWD on the R5 forests (was R6 model 2: fire kill 11 TPA short)"),
    ("ca", "23999387010900_simfire",  "FMCBA live fuel: CA fell into the eastern FULIV path (\"FFE fuel tables not ported\")"),
    ("ws", "23771657010900_salvage",  "same, WS"),
]

@testset "western FFE crash cases: run, and .sum rows == live" begin
    for (v, stem, why) in _FF_CASES
        txt, err = _ff_run(v, stem)
        @test (v, stem, err) == (v, stem, "")
        isempty(err) || continue
        # summary rows only: jl's .sum text also carries the FFE report tables (their year rows are short)
        jl = [split(l) for l in split(txt, '\n') if occursin(r"^\d{4} ", l) && length(split(l)) >= 20]
        lv = [split(l) for l in readlines(joinpath(_FF_TIER, v, stem * ".live.sum")) if occursin(r"^\d{4} ", l)]
        @test (v, stem, length(jl)) == (v, stem, length(lv))
        for (a, b) in zip(jl, lv)
            @test (v, stem, a) == (v, stem, b)
        end
    end
end
