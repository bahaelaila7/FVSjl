#!/usr/bin/env bash
# gen_allspecies_snls.sh — (re)generate the SN + LS all-species coverage fixtures and their
# live golden .sum files (test/fixtures/allspecies/), used by test_allspecies.jl.
#
# Companion to gen_allspecies.sh (CS+NE). Every species' growth/crown-width/volume/mortality/
# FORTYP coefficient row is exercised by a realistic mixed stand, validated against the live
# binary. Like NE, the live SN/LS binaries FPE on stands packed with too many all-unusual
# species at once (a FORTYP/stats limitation), so species are greedily packed into the largest
# stands the live binary will run (new species + canonical-inventory filler for DBH variety).
# Blank-alpha placeholder species (unaddressable by the 2-char field) are excluded.
#
#   SN: 90 species, snt01 inventory rows (STDINFO 80106 231Dd, SITECODE 63), live FVSsn.
#   LS: 68 species, lst01 inventory rows (STDINFO 903.0, SITECODE 2), live FVSls.
#
# Usage:  test/harness/gen_allspecies_snls.sh     # rebuild fixtures + refresh golden sums
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
FIX="$ROOT/fixtures/allspecies"; mkdir -p "$FIX"
JL="julia --project=$ROOT/.."
SN_TRE=/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.tre
LS_TRE=/workspace/ForestVegetationSimulator/tests/FVSls/lst01.tre

SN_CODES="$($JL -e 'using FVSjl; c=FVSjl.coefficients(Southern()).code_alpha; print(join([strip(c[i]) for i in 1:FVSjl.nspecies(Southern()) if strip(c[i])!=""],","))')"
LS_CODES="$($JL -e 'using FVSjl; c=FVSjl.coefficients(LakeStates()).code_alpha; print(join([strip(c[i]) for i in 1:FVSjl.nspecies(LakeStates()) if strip(c[i])!=""],","))')"

python3 - "$SN_TRE" "$LS_TRE" "$FIX" "$SN_CODES" "$LS_CODES" /tmp/FVSsn_new /tmp/FVSls_new <<'PY'
import sys, subprocess, os
sn_tre, ls_tre, FIX, sn_codes, ls_codes, SNBIN, LSBIN = sys.argv[1:8]
sn_codes=sn_codes.split(","); ls_codes=ls_codes.split(",")

def read_tre(p): return [l.rstrip("\r") for l in open(p).read().splitlines() if l.strip()]
def swap_species(base, code):  # overwrite cols 34-36 (A3), keep tree id (cols 24-27) intact
    line=list(base.ljust(80)); line[33:36]=list((code+"   ")[:3]); return "".join(line).rstrip()+"\r"
def write_tre(path, lines): open(path,"w",newline="\n").write("\n".join(lines)+"\n")

SN_KEY="""SCREEN
NOAUTOES
STATS
STDIDENT
SNALLSP  ALL-SPECIES MIX.
DESIGN                                        11.0       1.0
STDINFO        80106   231Dd        60.0     315.0      30.0       7.0
SITECODE          63      60.
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
LS_KEY="""SCREEN
NOAUTOES
STATS
STDIDENT
LSALLSP  ALL-SPECIES MIX.
DESIGN                                        11.0       1.0
STDINFO        903.0       1.0      60.0     315.0      30.0      10.0
SITECODE           2        60
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

def greedy_pack(prefix, codes, base_tre, key_text, BIN):
    net=read_tre(base_tre)
    def stand_runs(specieslist):
        out=[swap_species(b, specieslist[ti]) if ti<len(specieslist) else b+"\r" for ti,b in enumerate(net)]
        write_tre(f"{FIX}/_probe.tre", out); open(f"{FIX}/_probe.key","w").write(key_text.replace(f"{prefix.upper()}ALLSP","_PROBE"))
        return subprocess.run([BIN,"--keywordfile=_probe.key"],cwd=FIX,
                              stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL).returncode==0
    remaining=list(codes); gi=0
    while remaining:
        lo,hi,best=1,min(len(net),len(remaining)),1
        while lo<=hi:
            mid=(lo+hi)//2
            if stand_runs(remaining[:mid]): best=mid; lo=mid+1
            else: hi=mid-1
        grp=remaining[:best]; remaining=remaining[best:]
        out=[swap_species(b, grp[ti]) if ti<len(grp) else b+"\r" for ti,b in enumerate(net)]
        write_tre(f"{FIX}/{prefix}_cov{gi}.tre", out); open(f"{FIX}/{prefix}_cov{gi}.key","w").write(key_text)
        gi+=1
    for f in ("_probe.tre","_probe.key"):
        try: os.remove(f"{FIX}/{f}")
        except FileNotFoundError: pass
    return gi

sn_stands=greedy_pack("sn", sn_codes, sn_tre, SN_KEY, SNBIN)
ls_stands=greedy_pack("ls", ls_codes, ls_tre, LS_KEY, LSBIN)
print(f"SN: {sn_stands} stands / {len(sn_codes)} species; LS: {ls_stands} stands / {len(ls_codes)} species")
PY

# --- refresh live golden .sum files ------------------------------------------------------
for k in "$FIX"/sn_cov*.key; do b=$(basename "$k" .key)
  ( cd "$FIX" && /tmp/FVSsn_new --keywordfile="$b.key" >/dev/null 2>&1 && cp "$b.sum" "$b.live.sum" )
done
for k in "$FIX"/ls_cov*.key; do b=$(basename "$k" .key)
  ( cd "$FIX" && /tmp/FVSls_new --keywordfile="$b.key" >/dev/null 2>&1 && cp "$b.sum" "$b.live.sum" )
done
find "$FIX" -maxdepth 1 -type f ! -name '*.key' ! -name '*.tre' ! -name '*.live.sum' -delete
echo "SN+LS fixtures + golden sums written to $FIX"
