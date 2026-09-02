# MISTPINF forced initial dwarf-mistletoe infection (mistoe/misin.f opt 10 → misinf.f MISINF,
# activity 2006). Ported as kw_mistpinf! (parse + validate + schedule a ScheduledActivity in
# s.control.mistpinf) and dm_misinf! (applied each cycle after the spread model, mistoe.f:517),
# with the MISRAN JRAN LCG (misran.f) for the random method. (2026-09-02.)
#
# VALIDATED vs the live relinked /workspace/.crwork/FVScr_clean oracle on the 6-tree CR DM stand
# (dm6.tre: DF trees 1-3 infected DMR 3/4/2, DF trees 4-5 + PP tree 6 uninfected) with a
# `MISTPINF 1 3 <prop> <level> <method>` card (DF=species 3, date=cycle 1):
#   * NO MISTPINF card ⇒ BYTE-IDENTICAL to baseline (dm_misinf! returns on an empty schedule —
#     RNG-safe: the JRAN LCG is a SEPARATE stream from rann!, advanced only by a due card).
#   * LEVEL=1, prop=1.0 (or any prop that saturates the uninfectable trees), ALL methods (0/1/2):
#     the introduction cycle (2000) is BIT-EXACT vs the oracle — Mean_DMR 2.581, Inf_TPA 1414,
#     Mort_TPA 55 (baseline 1990 Inf_TPA 1110 raised to 1414: the two uninfected DF trees are
#     infected to DMR 1). Method is irrelevant at LEVEL=1 because every infected tree gets DMR 1,
#     so the visit ORDER and the tripled-record COUNT don't change the aggregate.
#   * CORNERED (LEVEL>1 / partial proportion): FVS runs MISINF over the TRIPLED treelist (ITRN
#     includes the tripling sub-records — the oracle assigns DMR 1,2,3 round-robin across the 3
#     sub-records of one original tree), whereas FVSjl carries ONE DMR per CENTRAL record (the
#     tripling sub-records are materialized later, at UPDATE). So the LEVEL>1 round-robin DMR
#     distribution and the exact random tree SELECTION for a partial proportion straddle the
#     oracle — the SAME central-record DMR / tripling-materialization gap the DM spread report
#     documents (test_dm_dbs). The DIRECTION is correct (higher LEVEL ⇒ higher Mean_DMR: oracle
#     LEVEL=3 raises 2000 Mean_DMR 2.581→~2.63; jl raises it the same way).

using Test
using FVSjl
using FVSjl: run_keyfile, CentralRockies
using SQLite, DBInterface

const _DM_FIX = joinpath(@__DIR__, "..", "fixtures", "dm", "dm6.tre")

function _run_mistpinf(dir; card::String = "")
    cp(_DM_FIX, joinpath(dir, "dm6.tre"))
    db = joinpath(dir, "dm6.db"); key = joinpath(dir, "dm6.key")
    open(key, "w") do io
        print(io, """
STDIDENT
DM6 CR MISTPINF
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
$(card)END
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
    rows = Dict{Int,NamedTuple}()
    for r in DBInterface.execute(d,
            "SELECT Year,Mean_DMR,Inf_TPA,Mort_TPA FROM FVS_DM_Spp_Sum ORDER BY Year")
        rows[Int(r.Year)] = (dmr = round(Float64(r.Mean_DMR), digits = 3),
                             inf = Int(round(r.Inf_TPA)), mort = Int(round(r.Mort_TPA)))
    end
    SQLite.close(d)
    return rows
end

@testset "MISTPINF forced initial dwarf-mistletoe infection (misin.f opt 10)" begin
    base = mktempdir(d -> _run_mistpinf(d; card = ""))
    # LEVEL=1, prop=1.0, methods 0 (random), 1 (hi→lo), 2 (lo→hi)
    m0 = mktempdir(d -> _run_mistpinf(d; card = "MISTPINF           1         3       1.0       1.0       0.0\n"))
    m1 = mktempdir(d -> _run_mistpinf(d; card = "MISTPINF           1         3       1.0       1.0       1.0\n"))
    m2 = mktempdir(d -> _run_mistpinf(d; card = "MISTPINF           1         3       1.0       1.0       2.0\n"))

    # 1) cyc0 (1990) is UNCHANGED by the card (the infection is introduced during cycle 1 → reported
    #    at 2000): the initial DMR from the damage codes is bit-exact vs the oracle.
    @test base[1990] == (dmr = 2.359, inf = 1110, mort = 56)
    for m in (m0, m1, m2)
        @test m[1990] == base[1990]
    end

    # 2) LEVEL=1 introduction cycle (2000) is BIT-EXACT vs FVScr_clean, and METHOD-INDEPENDENT
    #    (every infected tree ⇒ DMR 1, so the visit order / tripled-record count is irrelevant).
    #    Oracle 2000: Mean_DMR 2.581, Inf_TPA 1414, Mort_TPA 55.
    @test m0[2000] == (dmr = 2.581, inf = 1414, mort = 55)
    @test m1[2000] == m0[2000]
    @test m2[2000] == m0[2000]

    # 3) The forced infection RAISES the infected TPA over the no-card baseline (the two previously
    #    uninfected DF trees are now infected) — direction is correct.
    @test m0[2000].inf > base[2000].inf
    @test m0[2000].inf == 1414

    # 4) An out-of-range card is REJECTED at parse (misin.f validation) ⇒ NOT scheduled ⇒ the run is
    #    byte-identical to the no-card baseline. Bad species (99 > MAXSP 38), bad level (7 > 6),
    #    bad proportion (1.5 > 1), bad method (3 > 2).
    for bad in ("MISTPINF           1        99       1.0       1.0       0.0\n",
                "MISTPINF           1         3       1.0       7.0       0.0\n",
                "MISTPINF           1         3       1.5       1.0       0.0\n",
                "MISTPINF           1         3       1.0       1.0       3.0\n")
        r = mktempdir(d -> _run_mistpinf(d; card = bad))
        @test r == base
    end

    # 5) CORNERED (documented): LEVEL=3 raises the introduction-cycle Mean_DMR above LEVEL=1 (the
    #    round-robin assigns DMR>1) — DIRECTION matches the oracle (2000 Mean_DMR climbs from 2.581),
    #    but the exact value straddles the oracle because FVS round-robins over the tripled treelist
    #    while jl carries one DMR per central record (see the header note).
    l3 = mktempdir(d -> _run_mistpinf(d; card = "MISTPINF           1         3       1.0       3.0       2.0\n"))
    @test l3[2000].dmr > m0[2000].dmr        # higher LEVEL ⇒ higher Mean_DMR (both jl and oracle)
    @test l3[2000].inf == 1414               # same infected population, just a higher rating
end
