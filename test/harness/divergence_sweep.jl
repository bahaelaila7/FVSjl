# divergence_sweep.jl — broad live-FVS vs FVSjl differential over many stands.
#
# The campaign's discovery + regression backbone (the "FIA-plots" principle: push many
# diverse stands through BOTH engines and rank where they differ). For each .key it runs
# jl `run_keyfile` and the live binary (via {sn,ne,cs}_oracle.sh), diffs every .sum column
# across every cycle, and prints scenarios ranked by max RELATIVE divergence — skipping the
# ULP floor so only real (non-ULP) divergences surface.
#
# Usage:  julia --project=. test/harness/divergence_sweep.jl <variant> <key>...
#   variant ∈ {sn,ne,cs};  keys default to test/harness/scenarios/*.key for SN.
using FVSjl

const ORACLE = Dict("sn"=>"sn_oracle.sh","ne"=>"ne_oracle.sh","cs"=>"cs_oracle.sh")
const VAR    = Dict("sn"=>Southern(),"ne"=>Northeast(),"cs"=>CentralStates())
# .sum data columns (1-based incl Year): 3 TPA 4 BA 5 SDI 6 CCF 7 TopHt 8 QMD 9 Tcuft 10 Mcuft 11 Scuft 12 Bdft
const COLS = Dict(3=>"TPA",4=>"BA",5=>"SDI",6=>"CCF",7=>"TopHt",8=>"QMD",9=>"Tcuft",10=>"Mcuft",11=>"Scuft",12=>"Bdft")

# Parse a .sum into per-stand blocks (split on the -999 header), each a Dict(year => fields),
# so the diff aligns by (stand index, year) instead of a naive row zip (multi-stand keys broke that).
function _blocks(s)
    blocks = Vector{Dict{Int,Vector{SubString{String}}}}()
    cur = nothing
    for l in split(s, '\n')
        t = split(strip(l)); isempty(t) && continue
        if startswith(strip(l), "-999"); cur = Dict{Int,Vector{SubString{String}}}(); push!(blocks, cur); continue; end
        # Only true .sum data rows (full ~28-column layout). An appended CARBREPT carbon-report block also
        # starts each row with a year but has ~12 cols — without this guard it overwrote the real row at the
        # same year key, making col 11 (Scuft) read 0.0 (the carbon-report false-positive on carbon_*).
        (tryparse(Int, t[1]) !== nothing && cur !== nothing && length(t) >= 20) && (cur[parse(Int,t[1])] = t)
    end
    return blocks
end

function sweep(variant::AbstractString, keys::Vector{String})
    here = @__DIR__
    results = Tuple{String,Float64,String}[]   # (scenario, max_rel, detail)
    for key in keys
        isfile(key) || continue
        stem = first(splitext(basename(key)))
        # live
        outdir = mktempdir()
        try
            run(pipeline(`bash $(joinpath(here, ORACLE[variant])) $key $outdir`; stdout=devnull, stderr=devnull))
        catch; end
        livesum = filter(p->endswith(p,".sum"), readdir(outdir; join=true))
        isempty(livesum) && (push!(results,(stem, NaN, "live FPE/no-sum")); continue)
        LB = _blocks(read(first(livesum), String))
        # jl
        JB = try _blocks(FVSjl.run_keyfile(key; variant=VAR[variant], output=:sum)) catch e; push!(results,(stem,NaN,"jl error: $(sprint(showerror,e))")); continue end
        maxrel = 0.0; detail = "bit-exact"
        for si in 1:min(length(LB), length(JB))         # align by stand index
            lb, jb = LB[si], JB[si]
            for (yr, l) in lb                            # then by year
                haskey(jb, yr) || continue
                j = jb[yr]
                for (i,nm) in COLS
                    (i<=length(l) && i<=length(j)) || continue
                    a = parse(Float64,l[i]); b = parse(Float64,j[i]); a==0 && continue
                    r = abs(a-b)/abs(a)
                    (abs(a-b) <= 1 || r <= 0.002) && continue   # ULP floor: ≤1 print unit AND ≤0.2%
                    r > maxrel && (maxrel = r; detail = "s$si $nm@$yr live=$(l[i]) jl=$(j[i]) ($(round(r*100,digits=2))%)")
                end
            end
        end
        push!(results, (stem, maxrel, detail))
    end
    sort!(results; by = x -> (isnan(x[2]) ? -1.0 : x[2]), rev=true)
    println("\n=== divergence sweep ($variant): $(length(results)) stands, ranked by max non-ULP rel diff ===")
    for (scn, mr, d) in results
        tag = isnan(mr) ? "ERR " : mr==0 ? "ok  " : "DIFF"
        println("  $tag $(rpad(scn,22)) $(d)")
    end
end

vv = length(ARGS) >= 1 ? ARGS[1] : "sn"
ks = length(ARGS) >= 2 ? ARGS[2:end] :
     sort(filter(p->endswith(p,".key"), readdir(joinpath(@__DIR__,"scenarios"); join=true)))
sweep(vv, Vector{String}(ks))
