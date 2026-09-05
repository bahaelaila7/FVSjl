# coverage_breadth.jl — measure the breadth of FVSjl source code a VARIANT's run touches, to order the
# enriched regime-matrix sweep "widest-breadth-first" (USER 2026-09-05). Run UNDER coverage:
#   julia --project=. --code-coverage=user test/harness/fia/coverage_breadth.jl <VARIANT> <subdb> [regimes]
# It runs FVSjl.run_keyfile on every stand in <subdb> across the given regimes (default: the existing
# regime set none,plant,thinbba,salvage,simfire — which exercises base growth/mortality/volume + ESTAB
# (plant) + FFE (simfire) + management (thinbba/salvage)). Coverage counts are written by Julia to
# src/**/*.jl.<pid>.cov; the driver (coverage_breadth.sh) aggregates DISTINCT covered lines per variant.
using FVSjl, SQLite, DBInterface

const VARIANT = ARGS[1]
const SUBDB   = ARGS[2]
const REGIMES = length(ARGS) >= 3 ? split(ARGS[3], ",") : ["none","plant","thinbba","salvage","simfire"]

kwrec(kw, f...) = rpad(kw,10) * join(lpad(string(x),10) for x in f)
regime_block(r) =
    r == "simfire" ? "FMIn\n" * kwrec("SIMFIRE","2.0","10.00","1","50.0") * "\nEnd" :
    r == "thinbba" ? kwrec("THINBBA","2.0","40.0") :
    r == "salvage" ? kwrec("SALVAGE","2.0","0.0","999.0","0.9") :
    r == "plant"   ? "ESTAB\n" * kwrec("PLANT","2.0","3","400") * "\nEnd" : ""

keytext(cn, db, regime) = """
STDIDENT
$cn
DATABASE
DSNin
$db
StandSQL
SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
TreeSQL
SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
END
NUMCYCLE         5.0
$(regime_block(regime))
ECHOSUM
PROCESS
STOP
"""

var = FVSjl.variant_from_code(VARIANT)
db  = SQLite.DB(SUBDB)
cns = String[]
const NSTANDS = parse(Int, get(ENV, "COV_NSTANDS", "3"))
for r in DBInterface.execute(db, "SELECT DISTINCT STAND_CN FROM FVS_STANDINIT_COND WHERE STAND_CN IS NOT NULL LIMIT $NSTANDS")
    push!(cns, String(r.STAND_CN))
end
close(db)

nrun = Ref(0); nok = Ref(0)
dir = mktempdir()
key = joinpath(dir, "cov.key")
for cn in cns, rg in REGIMES
    write(key, keytext(cn, SUBDB, rg))
    nrun[] += 1
    try
        FVSjl.run_keyfile(key; variant=var); nok[] += 1
    catch e
        @warn "run failed" variant=VARIANT cn=cn regime=rg err=sprint(showerror, e)
    end
end
rm(dir; recursive=true, force=true)
println(stderr, "coverage_breadth: variant=$VARIANT stands=$(length(cns)) regimes=$(length(REGIMES)) runs=$(nrun[]) ok=$(nok[])")
