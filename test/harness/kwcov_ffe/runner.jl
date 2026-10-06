# runner.jl — jl-vs-live comparison for the keyword-isolation fixture directories (test/fixtures/kwcov_ffe/<v>/,
# test/fixtures/kwcov_on/): run jl on <stem>.key (absolute out.db path in a scratch dir) and diff every .sum field and
# every DBS table cell against the live goldens with tiered_runner's compare_case (print precision, no tolerance).
include(joinpath(@__DIR__, "..", "tiered", "tiered_runner.jl"))

const KCASE_ROOTS = Dict("ffe" => normpath(joinpath(@__DIR__, "..", "..", "fixtures", "kwcov_ffe")),
                         "on"  => normpath(joinpath(@__DIR__, "..", "..", "fixtures", "kwcov_on")),
                         "ext" => normpath(joinpath(@__DIR__, "..", "..", "fixtures", "kwcov_ext")))

"Stems (without .key) under <root>/<v>/ that have a live golden."
kcase_stems(root, v) = sort([f[1:end-4] for f in readdir(joinpath(root, v)) if endswith(f, ".key") &&
                            isfile(joinpath(root, v, f[1:end-4] * ".live.sum"))])

"Run jl on one case; returns (mismatches::Vector{Mismatch}, crashed::Bool, err)."
function run_kcase(root::AbstractString, v::AbstractString, stem::AbstractString)
    fx = joinpath(root, v); dir = mktempdir()
    cp(joinpath(fx, stem * ".key"), joinpath(dir, stem * ".key"))
    for e in (".tre", ".trl"); isfile(joinpath(fx, stem * e)) && cp(joinpath(fx, stem * e), joinpath(dir, stem * e)); end
    outdb = joinpath(dir, "out.db")
    indb = joinpath(dir, "stands.db")      # DB-read cases (ON): the shared sub-DB of the fixture directory
    isfile(joinpath(fx, "stands.db")) && cp(joinpath(fx, "stands.db"), indb)
    key = join([l == "out.db" ? outdb : l == "stands.db" ? indb : l for l in split(read(joinpath(dir, stem * ".key"), String), '\n')], '\n')
    write(joinpath(dir, stem * ".key"), key)
    txt = ""; crashed = false; err = ""
    try
        txt = FVSjl.run_keyfile(joinpath(dir, stem * ".key"); variant = VAR[uppercase(v)])
    catch e
        crashed = true; err = sprint(showerror, e)
    end
    ms = crashed ? [Mismatch(v, stem, "kcase", "CRASH", "crash", "*", "no crash", first(err, 200))] :
         compare_case(v, stem, "kcase", txt, isfile(outdb) ? outdb : ""; fx = fx, stem = stem)
    rm(dir; recursive = true, force = true)
    (ms, crashed, err)
end
