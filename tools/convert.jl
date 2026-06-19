#!/usr/bin/env julia
# convert.jl — translate legacy FVS inputs to the readable native formats.
#
#   julia --project tools/convert.jl stand.key [stand.yaml]
#   julia --project tools/convert.jl stand.tre [stand.csv]
#
# .key → .yaml (ordered keyword sequence), .tre → .csv (named tree columns).
# Output path defaults to the input with the new extension. The legacy files keep
# working directly (the engine auto-detects by extension), so this is for
# migration / human inspection, not a required step.

using FVSjl

function main(args)
    if isempty(args)
        println("usage: convert.jl <input.key|input.tre> [output.yaml|output.csv]")
        return 1
    end
    inp = args[1]
    ext = lowercase(splitext(inp)[2])
    if ext == ".key"
        out = length(args) >= 2 ? args[2] : first(splitext(inp)) * ".yaml"
        convert_key_to_yaml(inp, out)
        println("wrote $out")
    elseif ext == ".tre"
        out = length(args) >= 2 ? args[2] : first(splitext(inp)) * ".csv"
        convert_tre_to_csv(inp, out)
        println("wrote $out")
    else
        println("unrecognized input extension '$ext' (expected .key or .tre)")
        return 1
    end
    return 0
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    exit(main(ARGS))
end
