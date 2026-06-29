#!/usr/bin/env julia
# =============================================================================
# fvsjl-run.jl — run an FVS stand and write its summary.
#
#   julia --project bin/fvsjl-run.jl <key.{key,yaml}> [--variant SN|NE] [--output sum|csv] [-o out]
#
# Output goes to stdout (or a file with -o). Two run-level choices, because a legacy
# `.key`/`.tre` carries neither:
#   --variant SN|NE   which model to run as. A YAML's `variant:` is used when omitted; a
#                     `.key` defaults to SN — pass --variant to run it as NE. (flag overrides)
#   --output  sum|csv the summary format. `.sum` (legacy fixed-column) is the default; `csv`
#                     is the modern named-column form. A YAML's `output_format:` is used when
#                     omitted; an explicit --output overrides it. Default sum.
#
# The companion tree file (`<stem>.csv`/`.tre`) is found by base name (see docs/FORMATS.md).
# =============================================================================
using FVSjl

function main(args)
    pos = String[]; variant = nothing; output = nothing; outfile = nothing
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("--variant", "-v")
            variant = variant_from_code(args[i+1]); i += 2
        elseif a in ("--output", "-f")
            output = Symbol(lowercase(args[i+1])); i += 2
        elseif a == "-o"
            outfile = args[i+1]; i += 2
        elseif startswith(a, "--variant=")
            variant = variant_from_code(split(a, "=", limit = 2)[2]); i += 1
        elseif startswith(a, "--output=")
            output = Symbol(lowercase(split(a, "=", limit = 2)[2])); i += 1
        else
            push!(pos, a); i += 1
        end
    end
    if isempty(pos)
        println(stderr, "usage: fvsjl-run <key.{key,yaml}> [--variant SN|NE] [--output sum|csv] [-o outfile]")
        return 1
    end
    key = pos[1]
    isfile(key) || (println(stderr, "error: no such file: $key"); return 1)
    local txt
    try
        txt = run_keyfile(key; variant = variant, output = output)
    catch e
        println(stderr, "error: ", sprint(showerror, e)); return 1
    end
    if outfile === nothing
        print(txt)
    else
        write(outfile, txt); println("wrote ", outfile)
    end
    return 0
end

exit(main(ARGS))
