#!/usr/bin/env julia
# =============================================================================
# fvsjl-translate.jl — convert FVS stand inputs between the legacy fixed-column
# forms (.key / .tre) and the modern readable forms (.yaml / .csv), both ways.
#
#   julia --project bin/fvsjl-translate.jl <src> <dst> [tree-format]
#
# Direction is inferred from the extensions:
#   .key  → .yaml     keywords to readable YAML        (and .yaml → .key back)
#   .tre  → .csv      tree records to named-column CSV  (and .csv → .tre back)
#
# The engine reads either form directly (run_keyfile / read_tree_records), so this
# tool is for modernizing a legacy stand or producing a legacy file for stock FVS.
# `tree-format` overrides the .tre FORMAT string (default: the SN layout).
# =============================================================================
using FVSjl

function main(args)
    if length(args) < 2
        println(stderr, "usage: fvsjl-translate <src.{key,yaml,tre,csv}> <dst.{yaml,key,csv,tre}> [tree-format]")
        return 1
    end
    src, dst = args[1], args[2]
    isfile(src) || (println(stderr, "error: no such file: $src"); return 1)
    tree_fmt = length(args) >= 3 ? args[3] : FVSjl.DEFAULT_TREE_FORMAT
    try
        translate_io(src, dst; tree_fmt = tree_fmt)
    catch e
        println(stderr, "error: ", sprint(showerror, e))
        return 1
    end
    println("wrote ", dst)
    return 0
end

exit(main(ARGS))
