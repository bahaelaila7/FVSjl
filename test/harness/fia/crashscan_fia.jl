# crashscan_fia.jl — FAST jl-only robustness crash-hunt over real FIA plots at LARGE scale.
#
# The full differential sweeps (validate_fia.jl / manage_fia.jl) run live FVS per stand, so they
# are oracle-bound (~1 stand/s under FFE). A CRASH-HUNT does not need the oracle: it only asks
# "does FVSjl throw / segfault on this real plot?". Dropping the live subprocess + .sum diff makes
# this ~10× faster, so it can cover tens of thousands of stands per variant — the "exhaust all
# FVS-ready stands" dimension. It does NOT check fidelity (that's the differential harnesses); it
# only catches exceptions, OOB, and non-finite .sum output. Optional regime injects the fire path
# (the historically crash-prone FFE branch — slices 27 & 33 both lived there).
#
# Usage: FIA_DB=<subdb> julia --project=. test/harness/fia/crashscan_fia.jl <standlist> <SN|NE|CS|LS> [regime]
#   regime ∈ {none (default), simfire, thinbba, salvage, plant} — matches manage_fia's regime blocks.

using FVSjl
const DB = get(ENV, "FIA_DB", "/workspace/SQLite_FIADB_ENTIRE.db")
const VAR = Dict("SN"=>FVSjl.Southern(),"NE"=>FVSjl.Northeast(),"CS"=>FVSjl.CentralStates(),"LS"=>FVSjl.LakeStates())

kwrec(kw, fields...) = rpad(kw, 10) * join(lpad(string(f), 10) for f in fields)
regime_block(r) =
    r == "simfire" ? "FMIn\n" * kwrec("SIMFIRE", "2.0", "10.00", "1", "50.0") * "\nEnd" :
    r == "thinbba" ? kwrec("THINBBA", "2.0", "40.0") :
    r == "salvage" ? kwrec("SALVAGE", "2.0", "0.0", "999.0", "0.9") :
    r == "plant"   ? "ESTAB\n" * kwrec("PLANT", "2.0", "3", "400") * "\nEnd" :
                     ""   # none

keytext(cn, regime) = """
STDIDENT
$cn
DATABASE
DSNin
$DB
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

# a .sum with a non-finite / absurd value is as bad as a throw — flag it too
function sane_sum(text)
    got = false
    for ln in split(text, '\n')
        f = split(strip(ln)); length(f) < 8 && continue
        y = tryparse(Int, f[1]); (y===nothing || y<1000 || y>3000) && continue
        for i in 3:8
            v = tryparse(Float64, f[i]); v === nothing && continue
            (isfinite(v) && v >= 0 && v < 1e9) || return false
            got = true
        end
    end
    got
end

function main(listfile, v, regime)
    var = VAR[v]
    stands = [split(strip(l), '\t')[1] for l in eachline(listfile) if !isempty(strip(l))]
    n=0; ok=0; empty=0; crash=0; insane=0
    crashes = Tuple{String,String}[]
    dir = mktempdir()
    for cn in stands
        n += 1
        n % 500 == 0 && (print(stderr, "[$n/$(length(stands))] crash=$crash insane=$insane\r"); flush(stderr))
        kf = joinpath(dir, "s.key"); write(kf, keytext(cn, regime))
        out = try
            FVSjl.run_keyfile(kf; variant=var)
        catch e
            crash += 1
            length(crashes) < 40 && push!(crashes, (cn, first(sprint(showerror, e), 160)))
            continue
        end
        if isempty(out); empty += 1
        elseif !sane_sum(out); insane += 1; length(crashes) < 40 && push!(crashes, (cn, "INSANE_SUM"))
        else; ok += 1; end
    end
    println(stderr)
    println("\n===== CRASHSCAN ($regime): $v =====")
    println("stands=$n  ok=$ok  empty(no-sum)=$empty  CRASH(throw)=$crash  INSANE(sum)=$insane")
    if !isempty(crashes)
        println("--- first $(length(crashes)) crash/insane stands ---")
        for (cn, msg) in crashes; println("  $cn  $msg"); end
    else
        println("NO CRASHES — jl ran every stand to a finite .sum.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("usage: crashscan_fia.jl <standlist> <variant> [regime]")
    main(ARGS[1], ARGS[2], length(ARGS) >= 3 ? ARGS[3] : "none")
end
