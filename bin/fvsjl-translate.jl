#!/usr/bin/env julia
# =============================================================================
# fvsjl-translate.jl — convert FVS stand inputs between the legacy fixed-column
# forms (.key / .tre) and the modern readable forms (.yaml / .csv), both ways.
#
#   julia --project bin/fvsjl-translate.jl <src> <dst> [tree-format] [--flat]
#
# Direction is inferred from the extensions:
#   .key  → .yaml     keywords to ORDER-AWARE hierarchical YAML (and .yaml → .key)
#   .tre  → .csv      tree records to named-column CSV          (and .csv → .tre)
#
# `.key → .yaml` emits the grouped `stand:` form by default (sections only group the
# still-ordered keyword stream — see docs/KEYWORDS.md); pass `--flat` for the legacy
# single `keywords:` list. The engine reads either form directly (run_keyfile /
# read_tree_records), so this tool is for modernizing a legacy stand or producing a
# legacy file for stock FVS. `tree-format` overrides the .tre FORMAT (default: SN).
# =============================================================================
using FVSjl

function main(args)
    flat = "--flat" in args
    pos  = filter(a -> !startswith(a, "--"), args)
    if length(pos) < 2
        println(stderr, "usage: fvsjl-translate <src.{key,yaml,tre,csv}> <dst.{yaml,key,csv,tre}> [tree-format] [--flat]")
        return 1
    end
    src, dst = pos[1], pos[2]
    isfile(src) || (println(stderr, "error: no such file: $src"); return 1)
    tree_fmt = length(pos) >= 3 ? pos[3] : FVSjl.DEFAULT_TREE_FORMAT
    try
        # `--flat` only affects the .key → .yaml direction; for that case route through
        # the writer with the flag, otherwise use the generic extension-driven translator.
        if flat && lowercase(splitext(src)[2]) == ".key" &&
           lowercase(splitext(dst)[2]) in (".yaml", ".yml")
            FVSjl.write_keywords_yaml(FVSjl.read_keyfile_records(src), dst; flat = true)
        else
            translate_io(src, dst; tree_fmt = tree_fmt)
        end
    catch e
        println(stderr, "error: ", sprint(showerror, e))
        return 1
    end
    println("wrote ", dst)
    return 0
end

exit(main(ARGS))
