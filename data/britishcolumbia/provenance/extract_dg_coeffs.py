import re, csv
src=open("/workspace/ForestVegetationSimulator/canada/bc/dgf.f").read().split("\n")

def blocks(name):
    # collect DATA <name>(k) / ... / blocks by brace/paren tracking
    out={}
    i=0
    while i<len(src):
        m=re.match(rf"\s*DATA {name}\((\d+)\)\s*/",src[i])
        if m:
            k=int(m.group(1)); buf=[]; i+=1
            depth=0; started=False
            while i<len(src):
                line=src[i]
                # strip leading continuation '>' and trailing comments
                body=re.sub(r"^\s*[>&]","",line)
                buf.append(body)
                if "/" in body and buf and re.search(r"\)\s*/\s*$",body) and started:
                    break
                if "MD_STR" in body or "SS_STR" in body: started=True
                i+=1
            out[k]="\n".join(buf)
        i+=1
    return out

def strip_comment(s):
    # remove ! comments (not inside quotes)
    res=[];q=False
    for ch in s:
        if ch=="'":q=not q
        if ch=="!" and not q:break
        res.append(ch)
    return "".join(res)

def parse_md(txt, nscal):
    # find the two (/ ... /) arrays
    arrs=re.findall(r"\(/(.*?)/\)", txt, re.S)
    spp=[x.strip() for x in arrs[0].split(",")]
    spp=[int(x) for x in spp if x.strip()!=""]
    zones=re.findall(r"'([^']*)'", arrs[1])
    zones=[z.strip() for z in zones if z.strip()!=""]
    # scalars: text after the 2nd (/.../) close
    tail=txt.split("/)",2)[-1]  # after 2nd array
    scal=[]
    for ln in tail.split("\n"):
        c=strip_comment(ln).strip().rstrip(",").strip()
        if c in ("",")","MD_STR (","SS_STR (") or c.startswith(")"): 
            continue
        for tok in c.split(","):
            tok=tok.strip()
            if re.match(r"^[-+]?[0-9.][0-9.eE+-]*$",tok):
                scal.append(float(tok))
    return spp, zones, scal

ZN_FIELDS=["OBSERV","CON","CASP","SASP","EL","EL2","CCFA","LD","DSQ","DBAL1","DBAL2","CR","BAL","SIGMAR"]
zn=blocks("ZNKONST")
print("ZNKONST blocks:",len(zn))
rows=[]
for k in sorted(zn):
    spp,zones,scal=parse_md(zn[k],len(ZN_FIELDS))
    rows.append((k,";".join(str(x) for x in spp), ";".join(zones), scal))
    if k<=2: print(f" ZNK({k}) sp={spp[:3]} zones={zones} nscal={len(scal)} scal={scal}")
# write CSV
with open("/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/znkonst.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["idx","spp","zones"]+ZN_FIELDS)
    for k,sp,z,scal in rows:
        w.writerow([k,sp,z]+scal[:len(ZN_FIELDS)])
print("wrote znkonst.csv; scal-count check:", set(len(r[3]) for r in rows))

# --- SSKONST ---
ss=blocks("SSKONST")
print("\nSSKONST blocks:",len(ss))
srows=[]
for k in sorted(ss):
    spp,pretty,scal=parse_md(ss[k],1)
    srows.append((k,";".join(str(x) for x in spp),";".join(pretty),scal))
    if k<=2: print(f" SSK({k}) sp={spp[:1]} pretty={pretty} CON={scal}")
outdir="/workspace/FVSjl/data/britishcolumbia"
import os; os.makedirs(outdir,exist_ok=True)
with open(outdir+"/dg_znkonst.csv","w",newline="") as f:
    import csv as _c; w=_c.writer(f); w.writerow(["idx","spp","zones"]+ZN_FIELDS)
    for k,sp,z,scal in rows: w.writerow([k,sp,z]+scal[:len(ZN_FIELDS)])
with open(outdir+"/dg_sskonst.csv","w",newline="") as f:
    import csv as _c; w=_c.writer(f); w.writerow(["idx","spp","prettynames","CON"])
    for k,sp,z,scal in srows: w.writerow([k,sp,z]+(scal[:1] if scal else [0.0]))
print("CON-count check:", set(len(r[3]) for r in srows))
print("wrote", outdir+"/dg_znkonst.csv and dg_sskonst.csv")
