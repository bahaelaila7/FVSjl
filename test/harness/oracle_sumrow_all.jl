# oracle_sumrow_all.jl — FVSjulia (Oracle A) cycle-0 .sum row for every scenario.
# Runs FVSjulia.main per key (fresh .sum), prints "name | yr age TPA BA SDI CCF
# TopHt QMD Tcuft Mcuft Scuft Bdft ... fortype size stock". The C10 QC compares
# this to FVSjl's fvsjl_sumrow_all.jl output (cols 1-12 + classes) to surface
# where the FVSjl re-derivation diverges across species/fortypes.
import Pkg
using FVSjulia

dir = length(ARGS) >= 1 ? ARGS[1] : "/workspace/FVSjl/test/harness/scenarios"
keys = sort(filter(f -> endswith(f, ".key"), readdir(dir)))

for k in keys
    name = replace(k, ".key" => "")
    path = joinpath(dir, k)
    sumf = joinpath(dir, name * ".sum")
    isfile(sumf) && rm(sumf)
    try
        redirect_stdout(devnull) do
            FVSjulia.main(["--keywordfile=" * path])
        end
        if isfile(sumf)
            rows = String[]
            for ln in eachline(sumf)
                if !occursin("-999", ln) && !isempty(strip(ln))
                    push!(rows, strip(ln)); length(rows) >= 2 && break
                end
            end
            length(rows) >= 1 && println(name, " | ", rows[1])
            length(rows) >= 2 && println(name, "@1 | ", rows[2])
        else
            println(name, " | NOSUM")
        end
    catch e
        println(name, " | ERROR ", sprint(showerror, e))
    end
end
