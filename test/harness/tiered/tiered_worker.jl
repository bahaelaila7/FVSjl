# tiered_worker.jl — crash-isolated case runner. A jl SIGSEGV (e.g. an @inbounds out-of-range write) kills the
# process, so the fast tier / report / snapshot run cases in this worker subprocess and the parent restarts it past a
# crashing case (see run_isolated in tiered_runner.jl).
# Usage: julia --project=. tiered_worker.jl <mode compare|snapshot> <VARIANT> <cases.tsv> <out.tsv> <progress.log>
#   cases.tsv: "<cn>\t<regime>" per line. Appends result lines to out.tsv and "START/DONE <cn> <regime>" to progress.
#   out lines: compare → "M\t<file>\t<col>\t<year>\t<gold>\t<jl>\t<cn>\t<regime>" and "S\t<cn>\t<regime>\t<tpa>\t<ba>\t<tcuft>"
#              snapshot → "H\t<cn>\t<regime>\t<kind>\t<hash>"
mode, v, casesf, outf, progf = ARGS
include(joinpath(@__DIR__, mode == "snapshot" ? "snapshot.jl" : "tiered_runner.jl"))
_clean(s) = replace(string(s), r"[\t\r\n]+" => " ")

cases = [Tuple(split(l, '\t')) for l in eachline(casesf) if !isempty(strip(l))]
const _OUTLOCK = ReentrantLock()
emit(io, lines) = lock(_OUTLOCK) do
    foreach(l -> println(io, l), lines); flush(io)
end
open(outf, "a") do out
    open(progf, "a") do prog
        Threads.@threads for i in eachindex(cases)
            cn, r = String(cases[i][1]), String(cases[i][2])
            emit(prog, ["START $cn $r"])
            lines = String[]
            if mode == "snapshot"
                for (k, h) in snapshot_case(v, cn, r); push!(lines, join(("H", cn, r, k, h), '\t')); end
            else
                txt, db, crashed, err = run_case(v, cn, r)
                if crashed
                    push!(lines, join(("M", "CRASH", "*", "*", "", _clean(first(err, 300)), cn, r), '\t'))
                else
                    for m in compare_case(v, cn, r, txt, db)
                        push!(lines, join(("M", m.file, m.col, m.year, _clean(m.gold), _clean(m.got), cn, r), '\t'))
                    end
                    sg = case_signs(v, cn, r, txt)
                    sg === nothing || push!(lines, join(("S", cn, r, sg["TPA"], sg["BA"], sg["TCuFt"]), '\t'))
                end
            end
            emit(out, lines)
            emit(prog, ["DONE $cn $r"])
        end
    end
end
