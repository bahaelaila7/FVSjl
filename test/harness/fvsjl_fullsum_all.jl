# fvsjl_fullsum_all.jl — write FVSjl's FULL multi-cycle .sum for every scenario key
# into <outdir>, for the 3-way congruence sweep vs Oracle A / live Fortran.
# Uses run_keyfile (loops each_stand) so MULTI-STAND keys emit every stand — single-
# stand keys still emit one stand, but a 5-stand key like snt01 now produces all 5
# stand blocks (matching Fortran's 55-row .sum) instead of only stand 1.
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
        txt = FVSjl.run_keyfile(joinpath(dir, k); faithful = true)
        open(out, "w") do f; print(f, txt); end
    catch e
        open(out, "w") do f; println(f, "-999 ERROR ", sprint(showerror, e)); end
    end
end
println("wrote ", length(keys), " FVSjl .sum files (multi-stand aware) to ", outdir)
