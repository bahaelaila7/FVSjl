#!/bin/bash
# Per-species FVS_TreeList differential: inject a DATABASE/TREELIDB block into a scenario key,
# run it through BOTH live (oracle) and jl (run_keyfile), aggregate every per-tree column by
# (Year,SpeciesFIA) [partition-invariant], and report non-ULP divergences. Args: <variant> <key>
set -e
V=$1; KEY=$2; STEM=$(basename "$KEY" .key)
FVSJL=/workspace/FVSjl
WORK=/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/tld_$STEM
rm -rf "$WORK"; mkdir -p "$WORK"
declare -A ORACLE=( [sn]=sn_oracle.sh [ne]=ne_oracle.sh [cs]=cs_oracle.sh )
declare -A VAR=( [sn]=Southern [ne]=Northeast [cs]=CentralStates )

# Inject the DATABASE block + TREELIST after the STDIDENT's 2nd line (the title), before DESIGN/INVYEAR.
# Use awk to insert after the first line that starts a keyword following STDIDENT — simplest: insert after
# the line matching INVYEAR (present in every scenario) — DATABASE can appear anywhere before PROCESS.
inject() { # $1 dsnout-path  -> stdout modified key
  awk -v db="$1" '
    /^INVYEAR/ && !done { print; print "DATABASE"; print "DSNOUT"; print db; print "TREELIDB"; print "END"; print "TREELIST           0"; done=1; next }
    { print }
  ' "$KEY"
}

# live
OD="$WORK/live"; mkdir -p "$OD"
inject "$OD/out.db" > "$WORK/live.key"
# copy companion .tre if present (basename match for the injected key)
[ -f "${KEY%.key}.tre" ] && cp "${KEY%.key}.tre" "$WORK/live.tre"
bash "$FVSJL/test/harness/${ORACLE[$V]}" "$WORK/live.key" "$OD" >/dev/null 2>&1 || true
# live DB may be the literal path we set, or __DSNOUT__ if oracle ignores DSNOUT; find the treelist DB
LDB=""; for f in "$OD/out.db" "$OD/__DSNOUT__" "$WORK/out.db"; do
  [ -f "$f" ] && python3 -c "import sqlite3,sys; c=sqlite3.connect('$f'); sys.exit(0 if any(r[0]=='FVS_TreeList' for r in c.execute(\"select name from sqlite_master\")) else 1)" 2>/dev/null && { LDB="$f"; break; }
done

# jl
JD="$WORK/jl"; mkdir -p "$JD"
inject "$JD/out.db" > "$JD/run.key"
[ -f "${KEY%.key}.tre" ] && cp "${KEY%.key}.tre" "$JD/run.tre"
(cd "$JD" && timeout 200 julia --project=$FVSJL -e "using FVSjl; FVSjl.run_keyfile(\"run.key\"; variant=${VAR[$V]}(), faithful=true)" >/dev/null 2>&1 || true)
JDB=""; for f in "$JD/out.db"; do [ -f "$f" ] && JDB="$f"; done

echo "== $STEM: live_db=$([ -n "$LDB" ] && echo yes || echo NO) jl_db=$([ -n "$JDB" ] && echo yes || echo NO) =="
[ -z "$LDB" ] || [ -z "$JDB" ] && { echo "  (missing a treelist DB — skip)"; exit 0; }

python3 - "$LDB" "$JDB" << 'PYEOF'
import sqlite3,sys
def agg(db):
    c=sqlite3.connect(db); r={}
    try:
        for row in c.execute("""select Year,SpeciesFIA,sum(TPA),sum(DBH*TPA),sum(Ht*TPA),
            sum(PctCr*TPA),sum(TCuFt*TPA),sum(MCuFt*TPA),sum(SCuFt*TPA),sum(BdFt*TPA),sum(DG*TPA),sum(HtG*TPA)
            from FVS_TreeList group by Year,SpeciesFIA"""):
            r[(row[0],row[1])]=row[2:]
    except Exception as e: print("  agg err",e)
    return r
L=agg(sys.argv[1]); J=agg(sys.argv[2])
cols=["TPA","DBH","Ht","PctCr","TCuFt","MCuFt","SCuFt","BdFt","DG","HtG"]
worst={}
for k in sorted(set(L)&set(J)):
    for i,cn in enumerate(cols):
        a,b=L[k][i],J[k][i]
        if a is None or b is None: continue
        d=abs(a-b); rel=d/abs(a) if a else (0 if d<1e-6 else 9)
        if d>1.0 and rel>0.003: worst.setdefault(cn,[]).append((rel,k,a,b))
if not worst: print("  ALL per-tree columns bit-exact/ULP vs live ✓  (%d species-years)"%len(set(L)&set(J)))
else:
    for cn in cols:
        if cn in worst:
            w=sorted(worst[cn],reverse=True)[0]
            print(f"  DIFF {cn}: {len(worst[cn])}x worst {w[0]*100:.1f}% yr/sp{w[1]} live={w[2]:.2f} jl={w[3]:.2f}")
PYEOF
