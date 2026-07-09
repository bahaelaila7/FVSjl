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

# A crash-hunt success = jl produced a .sum with at least one projection DATA ROW. FVS's .sum is FIXED-WIDTH
# and adjacent columns can ABUT (e.g. "1995-999" = year 1995 + a -999 code with no separating space), so a
# whitespace-split of the whole row is unreliable — but the LEADING 4-digit year is always present and
# unambiguous. We only need "did jl run to a real projection row?"; value-level fidelity is the differential
# harness's job (validate_fia / manage_fia), not this crash scan. Non-finite/absurd values would surface there.
function has_data_row(text)
    for ln in split(text, '\n')
        s = strip(ln)
        length(s) >= 4 || continue
        y = tryparse(Int, s[1:4])
        (y !== nothing && 1000 <= y <= 3000) && return true
    end
    false
end

function main(listfile, v, regime)
    var = VAR[v]
    stands = [split(strip(l), '\t')[1] for l in eachline(listfile) if !isempty(strip(l))]
    n=0; ok=0; empty=0; crash=0
    crashes = Tuple{String,String}[]
    dir = mktempdir()
    for cn in stands
        n += 1
        n % 500 == 0 && (print(stderr, "[$n/$(length(stands))] crash=$crash\r"); flush(stderr))
        kf = joinpath(dir, "s.key"); write(kf, keytext(cn, regime))
        out = try
            FVSjl.run_keyfile(kf; variant=var)
        catch e
            crash += 1
            msg = first(sprint(showerror, e), 200)
            println(stderr, "\nCRASH $cn :: $msg"); flush(stderr)   # log immediately (survives interruption)
            length(crashes) < 40 && push!(crashes, (cn, msg))
            continue
        end
        has_data_row(out) ? (ok += 1) : (empty += 1)   # empty = no projection row (e.g. nonstocked; not a crash)
    end
    println(stderr)
    println("\n===== CRASHSCAN ($regime): $v =====")
    println("stands=$n  ok(projected)=$ok  empty(no-data-row)=$empty  CRASH(throw/OOB)=$crash")
    if !isempty(crashes)
        println("--- first $(length(crashes)) CRASH stands ---")
        for (cn, msg) in crashes; println("  $cn  $msg"); end
    else
        println("NO CRASHES — jl ran every stand without throwing.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("usage: crashscan_fia.jl <standlist> <variant> [regime]")
    main(ARGS[1], ARGS[2], length(ARGS) >= 3 ? ARGS[3] : "none")
end
