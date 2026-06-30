import shutil
FMT1="(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,"
FMT2="T52,I2,T66,5I1,T54,7I1,T75,F3.0)"
DESIGN="DESIGN                                        11.0       1.0"
SI={"NE":"STDINFO        922.0                60.0     315.0      30.0      20.0",
    "SN":"STDINFO        80106   232BA        60.0     315.0      30.0       7.0"}
def scen(v,idx,clen,fo,thin):
    inv=1990; ncyc=8
    rid=f"SCN{v}{idx:02d}  clen{clen} f{fo or 0} {thin or 'none'}"
    L=["SCREEN","NOAUTOES","STDIDENT",rid,DESIGN,SI[v],"SITECODE          13        56",
       f"INVYEAR       {inv}.0", f"NUMCYCLE        {ncyc}.0",
       "TIMEINT"+" "*13+f"{clen:10d}"]   # explicit cycle length so fires land on boundaries
    if thin=="THINDBH": L.append("THINDBH                               4.")
    if fo: y=inv+fo*clen; L += ["FMIN", f"SIMFIRE         {y}     10.00         1      50.0", "END"]
    L += ["TREEFMT", FMT1, FMT2, "TREEDATA","ECHOSUM","PROCESS","STOP"]
    return "\r".join(L)+"\r"
scn=[(v,i,clen,fo,thin) for v in ("NE","SN") for i,(clen,fo,thin) in enumerate(
      [(c,f,t) for c in (5,10) for f in (None,1,2,3) for t in (None,"THINDBH")])]
for (v,i,clen,fo,thin) in scn:
    n=f"{v.lower()}{i:02d}"
    open(f"/tmp/sweep/scn/{n}.key","w",newline="").write(scen(v,i,clen,fo,thin))
    src="/workspace/ForestVegetationSimulator/tests/FVSne/net01.tre" if v=="NE" else "/workspace/FVSjl/test/harness/scenarios/snt01_alpha.tre"
    shutil.copy(src,f"/tmp/sweep/scn/{n}.tre")
print(len(scn),"scenarios")
