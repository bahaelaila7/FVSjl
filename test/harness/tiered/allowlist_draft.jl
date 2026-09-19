# allowlist_draft.jl — draft KNOWN_RESIDUALS.toml entries from MEASURED fast-tier mismatches (report.jl TSV).
# Groups by (variant, regime, file, col); each entry is OPEN with an issue text = a rule-assigned cause tag + the
# measured class counts (classify.jl) + a sample, and max_cells = the measured count (spread ⇒ fail). The tag rules
# are explicit below — review every NEW tag before committing (they are findings, reported to the coordinator).
# Usage: julia test/harness/tiered/allowlist_draft.jl <report.tsv> [owner_note] >> KNOWN_RESIDUALS.toml
include(joinpath(@__DIR__, "classify_fns.jl"))

function tag(variant, file, col, classes::Dict{String,Int}, maxabs::Float64, sample)
    owner_drift = variant == "BM" ? "fork bm-base-resid3" : variant == "IE" ? "fork ie-ulp3" :
                  "the $(variant) campaign (not yet dug)"
    dom = isempty(classes) ? "" : first(sort(collect(classes); by = x -> -x[2]))[1]
    driftonly = all(k -> k in ("ulp1", "ulp2_16", "ulpbig", "tiny", "small"), keys(classes))
    if file in ("FVS_TreeList", "FVS_CutList", "FVS_ATRTList") && col in ("Ht2TDCF", "Ht2TDBF")
        return "list-table merch-top heights not filled by jl (written 0) — owner: fork list-ht2td"
    elseif file == "FVS_InvReference"
        return "NEW: jl FVS_InvReference merch-spec columns (CFMinDBH/CFTopDia/BFMinDBH/…) written as 0 — DBS writer gap"
    elseif file == "sum" && col == "EXTRA_CONTENT"
        return "NEW: jl run_keyfile's .sum text includes FFE report tables (CARBON REPORT); FVS writes those to .out, its .sum has only -999 + summary rows"
    elseif col == "PRESENCE"
        return "NEW: DBS table presence differs from live ($(sample)) — table emission/gating gap"
    elseif file == "CRASH"
        return "jl CRASH — $(sample)"
    elseif file == "TALLY"
        return "one-directional population bias ($(sample)) — owner: $(owner_drift)"
    elseif file == "FVS_Summary" && col in ("QMD", "ATQMD") && dom == "ulpbig"
        return "NEW: jl FVS_Summary stores QMD/ATQMD rounded to 1 decimal (print value); FVS stores the full REAL*4"
    elseif dom == "gold_text9"
        return "NEW: FVS stores this table's REALs as 9-significant-digit list-directed text; jl stores full double"
    elseif dom == "jl_not_f32"
        return "NEW: jl writer stores a Float64-computed value where FVS stores REAL*4"
    elseif file == "sum" && maxabs <= 1.0
        return "±1 print-boundary .sum cells from per-record Float32 drift — owner: $(owner_drift)"
    elseif driftonly
        return "per-record Float32 drift below .sum print precision — owner: $(owner_drift)"
    else
        return "material divergence (model/report) — needs dig; owner: $(owner_drift)"
    end
end

function draft(path)
    groups = Dict{NTuple{4,String},Vector{Vector{SubString{String}}}}()
    for (i, l) in enumerate(eachline(path))
        i == 1 && continue
        f = split(l, '\t'); length(f) < 8 && continue
        push!(get!(groups, (f[1], f[3], f[4], f[5]), Vector{SubString{String}}[]), f)
    end
    date = Base.Libc.strftime("%Y-%m-%d", Base.time())
    for k in sort(collect(keys(groups)))
        rs = groups[k]; cl = Dict{String,Int}(); maxabs = 0.0
        for f in rs
            c = cls(f[7], f[8]); cl[c] = get(cl, c, 0) + 1
            a = tryparse(Float64, f[7]); b = tryparse(Float64, f[8])
            (a !== nothing && b !== nothing) && (maxabs = max(maxabs, abs(a - b)))
        end
        s = rs[1]; sample = "e.g. stand $(s[2]) year $(s[6]): live=$(first(s[7], 40)) jl=$(first(s[8], 60))"
        t = tag(k[1], k[3], k[4], cl, maxabs, k[4] == "PRESENCE" ? "live=$(s[7]), jl=$(s[8])" :
                k[3] in ("TALLY", "CRASH") ? s[8] : "")
        println("\n[[residual]]")
        println("variant = \"$(k[1])\"\nregime = \"$(k[2])\"\nfile = \"$(k[3])\"\ncol = \"$(k[4])\"\nstand = \"*\"\nyear = \"*\"")
        println("status = \"OPEN\"")
        println("max_cells = $(length(rs))")
        cs = join(["$c=$n" for (c, n) in sort(collect(cl))], " ")
        println("issue = \"\"\"$t | measured $date: $(length(rs)) cells [$cs] max|Δ|=$(round(maxabs, sigdigits = 4)) | $(replace(sample, '"' => '\''))\"\"\"")
    end
end
draft(ARGS[1])
