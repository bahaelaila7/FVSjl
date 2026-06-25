#!/usr/bin/env bash
# validate.sh — keyword-coverage fidelity gates for every scenarios/*.key:
#   (A) FVSjl key/tre  ==  FVSjl yaml/csv   (byte-identical; same engine)
#   (B) FVSjl key/tre  ~=  live FVSsn       (sumdiff: abs<=1 or rel<=0.1% ULP)
# Batched: one Julia process for all FVSjl runs, one for all sumdiffs.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$HERE/../.." && pwd)"
SCEN="$HERE/scenarios"; FTBASE=/workspace/FVSjulia/tests/fortran_baseline.sh
JL="julia --project=$ROOT"

# 1. structured YAML for every scenario (Python writer)
for key in "$SCEN"/*.key; do
  python3 "$ROOT/examples/key_to_structured_yaml.py" "$key" > "${key%.key}.yaml" 2>/dev/null
done
# 2. FVSjl: build csvs, run key+yaml, write <name>.key.sum, report yaml==key  (ONE process)
echo "== FVSjl (key vs yaml) =="
$JL "$HERE/run_fvsjl.jl" "$SCEN" 2>/dev/null | tee /tmp/kc_fvsjl.txt
# 3. FVSsn on every key
for key in "$SCEN"/*.key; do
  name=$(basename "$key" .key)
  bash "$FTBASE" "$key" "$SCEN/ft_$name" >/dev/null 2>&1
done
# 4. sumdiff each (ONE process)
echo "== FVSjl key vs FVSsn (ULP) =="
$JL -e '
function rows(p); o=Vector{Vector{Float64}}(); for ln in eachline(p); occursin("-999",ln)&&continue; t=split(strip(ln)); isempty(t)&&continue; v=Float64[]; ok=true; for x in t; n=tryparse(Float64,x); n===nothing&&(ok=false;break); push!(v,n); end; ok&&push!(o,v); end; o; end
function cmp(a,b); length(a)!=length(b) && return "rows $(length(a))/$(length(b))"; for (r,(ra,rb)) in enumerate(zip(a,b)); length(ra)!=length(rb) && return "row$r width"; for (c,(x,y)) in enumerate(zip(ra,rb)); d=abs(x-y); (d>1.0 && d>0.001*max(abs(x),abs(y))) && return "r$r c$c $x/$y"; end; end; ""; end
scen="'"$SCEN"'"
for key in sort(filter(f->endswith(f,".key"), readdir(scen;join=true)))
    name=first(splitext(basename(key)))
    js=joinpath(scen,name*".key.sum"); fs=joinpath(scen,"ft_"*name,name*".sum")
    if isfile(js) && isfile(fs)
        d=cmp(rows(js),rows(fs)); println(name, "\t", isempty(d) ? "bPASS" : "bFAIL "*d)
    else
        println(name, "\tbNOSUM")
    end
end' 2>/dev/null
