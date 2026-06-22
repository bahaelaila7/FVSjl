# fvsjl_fullsum_all.jl — write FVSjl's FULL multi-cycle .sum for every (single-stand)
# scenario key into <outdir>, for the natural-process congruence sweep vs Oracle A.
# Usage: julia --project fvsjl_fullsum_all.jl <scenario-dir> <outdir> [glob]
using FVSjl

dir    = ARGS[1]
outdir = ARGS[2]; mkpath(outdir)
pat    = length(ARGS) >= 3 ? ARGS[3] : ""
keys   = sort(filter(f -> endswith(f, ".key") && (isempty(pat) || occursin(pat, f)), readdir(dir)))

for k in keys
    name = replace(k, ".key" => "")
    out  = joinpath(outdir, name * ".sum")
    try
        s, _ = initialize(joinpath(dir, k))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
        io = IOBuffer(); FVSjl.write_sum_file(io, s)
        open(out, "w") do f; print(f, String(take!(io))); end
    catch e
        open(out, "w") do f; println(f, "-999 ERROR ", sprint(showerror, e)); end
    end
end
println("wrote ", length(keys), " FVSjl .sum files to ", outdir)
