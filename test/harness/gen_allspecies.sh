#!/usr/bin/env bash
# gen_allspecies.sh — (re)generate the CS+NE all-species coverage fixtures and their
# live golden .sum files (test/fixtures/allspecies/), used by test_allspecies.jl.
#
# Every species' growth/crown-width/volume/mortality/FORTYP coefficient row is exercised
# by a realistic mixed stand, validated against the live binary (the BW basswood crown-
# width gap was found and fixed via this sweep).
#
#   CS: one mixed stand, one tree of each of the 96 species (cst01 inventory rows with
#       the species field [cols 34-36] swapped; tree IDs at cols 24-27 kept intact).
#       Live FVScs runs it directly.
#   NE: the live binary FPEs on stands of many all-unusual species at once (a FORTYP/
#       stats limitation), so the 107 non-blank species are greedily packed (new species
#       + net01-real filler trees) into stands the live binary runs. Blank-alpha species
#       (NE 71) is excluded — unaddressable by the 2-char field.
#
# Usage:  test/harness/gen_allspecies.sh        # rebuild fixtures + refresh golden sums
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
FIX="$ROOT/fixtures/allspecies"; mkdir -p "$FIX"
JL="julia --project=$ROOT/.."
CS_TRE=/workspace/ForestVegetationSimulator/tests/FVScs/cst01.tre
NE_TRE=/workspace/ForestVegetationSimulator/tests/FVSne/net01.tre

# --- CS mixed all-species stand ----------------------------------------------------------
CS_CODES="$($JL -e 'using FVSjl; c=FVSjl.coefficients(CentralStates()).code_alpha; print(join([strip(c[i]) for i in 1:FVSjl.nspecies(CentralStates())],","))')"
NE_CODES="$($JL -e 'using FVSjl; c=FVSjl.coefficients(Northeast()).code_alpha; print(join([strip(c[i]) for i in 1:FVSjl.nspecies(Northeast()) if strip(c[i])!=""],","))')"

python3 - "$CS_TRE" "$NE_TRE" "$FIX" "$CS_CODES" "$NE_CODES" /workspace/FVSjl/tmp/oracles/FVScs_new /workspace/FVSjl/tmp/oracles/FVSne_new <<'PY'
import sys, subprocess, os
cs_tre, ne_tre, FIX, cs_codes, ne_codes, CSBIN, NEBIN = sys.argv[1:8]
cs_codes=cs_codes.split(","); ne_codes=ne_codes.split(",")

def read_tre(p): return [l.rstrip("\r") for l in open(p).read().splitlines() if l.strip()]
def swap_species(base, code):  # overwrite cols 34-36 (A3), keep tree id (cols 24-27) intact
    line=list(base.ljust(80)); line[33:36]=list((code+"   ")[:3]); return "".join(line).rstrip()+"\r"
def write_tre(path, lines): open(path,"w",newline="\n").write("\n".join(lines)+"\n")

CS_KEY="""SCREEN
NOAUTOES
STATS
STDIDENT
CSALLSP  ALL-SPECIES MIX.
DESIGN                                        11.0       1.0
STDINFO          905                60.0     315.0      30.0      10.0       40.
SITECODE          19        60
INVYEAR       1990.0
NUMCYCLE        10.0
TREEFMT
(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,
T52,I2,T66,5I1,T54,7I1,T75,F3.0)
TREEDATA
ECHOSUM
PROCESS
STOP
"""
NE_KEY=CS_KEY.replace("CSALLSP","NEALLSP").replace(
  "STDINFO          905                60.0     315.0      30.0      10.0       40.",
  "STDINFO        922.0                60.0     315.0      30.0      20.0").replace(
  "SITECODE          19        60","SITECODE          13        56")

# CS: mixed stand, one tree per species (cycle cst01 rows for DBH variety)
cst=read_tre(cs_tre)
cs_lines=[swap_species(cst[i%len(cst)], code) for i,code in enumerate(cs_codes)]
write_tre(f"{FIX}/cs_allsp.tre", cs_lines); open(f"{FIX}/cs_allsp.key","w").write(CS_KEY)

# NE: greedy-pack species into stands that the live binary runs (new species + net01 filler)
net=read_tre(ne_tre); open(f"{FIX}/_ne.key","w").write(NE_KEY)
def ne_stand_runs(specieslist):
    out=[swap_species(b, specieslist[ti]) if ti<len(specieslist) else b+"\r" for ti,b in enumerate(net)]
    write_tre(f"{FIX}/_probe.tre", out); open(f"{FIX}/_probe.key","w").write(NE_KEY)
    return subprocess.run([NEBIN,"--keywordfile=_probe.key"],cwd=FIX,
                          stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL).returncode==0
remaining=list(ne_codes); gi=0
while remaining:
    lo,hi,best=1,min(len(net),len(remaining)),1
    while lo<=hi:
        mid=(lo+hi)//2
        if ne_stand_runs(remaining[:mid]): best=mid; lo=mid+1
        else: hi=mid-1
    grp=remaining[:best]; remaining=remaining[best:]
    out=[swap_species(b, grp[ti]) if ti<len(grp) else b+"\r" for ti,b in enumerate(net)]
    write_tre(f"{FIX}/ne_cov{gi}.tre", out); open(f"{FIX}/ne_cov{gi}.key","w").write(NE_KEY)
    gi+=1
for f in ("_probe.tre","_probe.key","_ne.key"):
    try: os.remove(f"{FIX}/{f}")
    except FileNotFoundError: pass
print(f"CS: 1 stand / {len(cs_codes)} species; NE: {gi} stands / {len(ne_codes)} species")
PY

# --- refresh live golden .sum files ------------------------------------------------------
( cd "$FIX" && /workspace/FVSjl/tmp/oracles/FVScs_new --keywordfile=cs_allsp.key >/dev/null 2>&1 && cp cs_allsp.sum cs_allsp.live.sum )
for k in "$FIX"/ne_cov*.key; do b=$(basename "$k" .key)
  ( cd "$FIX" && /workspace/FVSjl/tmp/oracles/FVSne_new --keywordfile="$b.key" >/dev/null 2>&1 && cp "$b.sum" "$b.live.sum" )
done
# clean live run artifacts (keep only .key/.tre/.live.sum)
find "$FIX" -maxdepth 1 -type f ! -name '*.key' ! -name '*.tre' ! -name '*.live.sum' -delete
echo "fixtures + golden sums written to $FIX"
