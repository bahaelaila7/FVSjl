import re, csv, os
src=open("/workspace/ForestVegetationSimulator/canada/bc/htgf.f").read().split("\n")
def blocks(name):
    out={}; i=0
    while i<len(src):
        m=re.match(rf"\s*DATA {name}\((\d+)\)\s*/",src[i])
        if m:
            k=int(m.group(1)); buf=[]; i+=1; started=False
            while i<len(src):
                body=re.sub(r"^\s*[>&]","",src[i]); buf.append(body)
                if re.search(r"\)\s*/\s*$",body) and started: break
                if "LTHG_STR" in body: started=True
                i+=1
            out[k]="\n".join(buf)
        i+=1
    return out
def strip_comment(s):
    r=[];q=False
    for ch in s:
        if ch=="'":q=not q
        if ch=="!" and not q:break
        r.append(ch)
    return "".join(r)
def parse(txt):
    arrs=re.findall(r"\(/(.*?)/\)", txt, re.S)
    spp=[int(x) for x in arrs[0].split(",") if x.strip()!=""]
    zones=[z.strip() for z in re.findall(r"'([^']*)'", arrs[1]) if z.strip()!=""]
    tail=txt.split("/)",2)[-1]
    scal=[]
    for ln in tail.split("\n"):
        c=strip_comment(ln).strip().rstrip(",").strip()
        for tok in c.split(","):
            tok=tok.strip()
            if re.match(r"^[-+]?[0-9.][0-9.eE+-]*$",tok): scal.append(float(tok))
    return spp,zones,scal
FIELDS=["SI","B3","C","CN","CNSI","DG","DG2","DBH1","DBH2","HT2"]
bl=blocks("LTHG"); print("LTHG blocks:",len(bl))
rows=[]
for k in sorted(bl):
    spp,zones,scal=parse(bl[k]); rows.append((k,spp,zones,scal))
    if k<=2: print(f" LTHG({k}) sp={spp[:3]} nzones={len(zones)} nscal={len(scal)} scal={scal}")
print("scal-count set:", set(len(r[3]) for r in rows))
outdir="/workspace/FVSjl/data/britishcolumbia"
with open(outdir+"/htg_lthg.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["idx","spp","zones"]+FIELDS)
    for k,spp,zones,scal in rows:
        w.writerow([k, ";".join(map(str,spp)), ";".join(zones)]+scal[:len(FIELDS)])
print("wrote htg_lthg.csv:", len(rows), "rows")
