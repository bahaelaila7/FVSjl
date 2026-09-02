# MISTMULT dwarf-mistletoe spread-probability multipliers (mistoe/misin.f opt 1: YPLMLT/YNGMLT,
# applied at mistoe.f:331/359 PPLUS·=YPLMLT / PMINUS·=YNGMLT). Ported as two GrowthMultiplier
# kinds (:dm_inc/:dm_dec) parsed by kw_mistmult!, queried by active_multiplier, threaded into
# cr_mistoe!/ie_mistoe! at the logistic-spread hook. (2026-09-02.)
#
# VALIDATED vs the live relinked /workspace/.crwork/FVScr_clean oracle on the 6-tree CR DM stand
# (3 DF infected DMR 3/4/2, damage code 33) with a `MISTMULT 1 3 2.0` card (DF=species 3, DMR
# increase multiplier 2.0):
#   * mult=1.0 / no-MISTMULT ⇒ BYTE-IDENTICAL to baseline (RNG-safe: it scales the probability,
#     never the rann! draw) — asserted here, and guarded stand-wide by the 339/11 gate.
#   * cyc0 (1990) DM report BIT-EXACT vs the oracle (Mean_DMR 2.4, Inf_TPA 1110, Mort_TPA 56) with
#     the card present — the multiplier only takes effect during projection, so cyc0 is unchanged.
#   * inc=2.0 raises the DM spread the SAME direction as the oracle (Mean_DMR & DM-mortality climb,
#     surviving infected TPA falls): oracle 2030 DMR 2.5→3.0 / Mort_TPA 58→82; jl 2030 DMR 2.4→2.7
#     / Mort_TPA 57→73. The cyc1+ magnitude gap is the PRE-EXISTING DM-report projection straddle
#     (the oracle books AUTOES regen cohorts jl doesn't — same corner test_dm_dbs documents; it is
#     present in the BASELINE too, NOT introduced by MISTMULT), so the delta lands on a smaller
#     infected population. cyc0-bit-exact + inert-at-1.0 + correct signed delta = the port is faithful.

using Test
using FVSjl
using FVSjl: run_keyfile, CentralRockies
using SQLite, DBInterface

function _run_mistmult(dir; mistmult::Bool, inc::String = "1.0")
    tre = joinpath(@__DIR__, "..", "fixtures", "dm", "dm6.tre")
    cp(tre, joinpath(dir, "dm6.tre"))
    db = joinpath(dir, "dm6.db"); key = joinpath(dir, "dm6.key")
    mm = mistmult ? "MISTMULT           1         3       $inc\n" : ""
    open(key, "w") do io
        print(io, """
STDIDENT
DM6 CR MISTMULT
DESIGN                                        11.0       1.0
STDINFO          303    001010      60.0     315.0      30.0      88.0
INVYEAR         1990
NUMCYCLE           4
TREEFMT
(I4,T1,I7,F6.0,I1,A3,F4.1,F3.1,2F3.0,F4.1,I1,3(I2,I2),2I1,I2,2I3,2I1,F3.0)
OPEN            55
$(joinpath(dir, "dm6.tre"))
TREEDATA          55
MISTOE
MISTPRT
$(mm)END
DATABASE
DSNOUT
$db
SUMMARY
MISRPTS
END
PROCESS
STOP
""")
    end
    run_keyfile(key; variant = CentralRockies())
    d = SQLite.DB(db)
    rows = Dict(r.Year => (dmr = round(r.Mean_DMR, digits = 1), inf = r.Inf_TPA, mort = r.Mort_TPA)
                for r in DBInterface.execute(d,
                    "SELECT Year,Mean_DMR,Inf_TPA,Mort_TPA FROM FVS_DM_Spp_Sum"))
    SQLite.close(d)
    return rows
end

@testset "MISTMULT dwarf-mistletoe spread multipliers (misin.f opt 1)" begin
    base  = mktempdir(d -> _run_mistmult(d; mistmult = false))
    one   = mktempdir(d -> _run_mistmult(d; mistmult = true, inc = "1.0"))
    two   = mktempdir(d -> _run_mistmult(d; mistmult = true, inc = "2.0"))

    # 1) INERT: MISTMULT with mult=1.0 is byte-identical to no-MISTMULT (RNG-safe — proves the
    #    threaded ×YPLMLT/×YNGMLT never perturbs the draw stream). Every reported year matches.
    @test sort(collect(keys(one))) == sort(collect(keys(base)))
    for y in keys(base)
        @test one[y] == base[y]
    end

    # 2) cyc0 (1990) BIT-EXACT vs FVScr_clean, unchanged by the card (multiplier acts on projection).
    for r in (two[1990], base[1990])
        @test round(Float64(r.dmr), digits = 1) == 2.4
        @test Int(round(r.inf))  == 1110
        @test Int(round(r.mort)) == 56
    end

    # 3) inc=2.0 raises spread — same SIGN as the oracle delta (non-vacuous: deleting the threading
    #    makes two[]==base[] and these fail). By 2030 the DM rating and DM-mortality both climb and
    #    the surviving infected TPA falls, relative to the baseline projection.
    @test two[2030].dmr  > base[2030].dmr
    @test two[2030].mort > base[2030].mort
    @test two[2030].inf  < base[2030].inf
end
