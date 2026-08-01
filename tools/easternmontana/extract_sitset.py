#!/usr/bin/env python3
"""Extract EM sitset site/SDI tables from em/sitset.f for site_setup!.

Tables (all clean DATA in em/sitset.f): BAMAXA(30) BA-max per ITYPE; SDICON(9) SDImax per
SDI-group; MAPSS(30) site-species per ITYPE; MAPSDI(122) IEMTYP→SDI-group; MAPSIT(30,11)
site index per ITYPE × site-species-group (col-major). Logic (sitset.f): ISISP=MAPSS(ITYPE);
per species I, SITEAR(I)=MAPSIT(ITYPE,grp(I)); ISDI=MAPSDI(IEMTYP) (else 5); SDIDEF(I)=
BAMAX>0 ? BAMAX/(0.5454154*PMSDIU/100) : SDICON(ISDI); BAMAX=BAMAXA(ITYPE) if unset.
Species→site-group grp(I): 4,5→1; 6→2; 12→3; 11,13-16,19→4; 17→5; 1,2,18→6; 10→7; 3→8;
7→9; 8→10; 9→11; else 70.0.  Emits data/easternmontana/site_*.csv.
Validated vs live FVSem emt01 (ITYPE=4): SDIDEF=696 all species (= live "SDI MAX" table).
"""
import re
def expand(body):
    body = re.sub(r"!.*", "", body)
    out = []
    for t in re.findall(r"[-\d.]+\*[-\d.]+|[-\d.]+", body):
        if "*" in t: n, v = t.split("*"); out += [float(v)]*int(float(n))
        else: out.append(float(t))
    return out
def grab(name):
    m = re.search(r"DATA\s+"+name+r"\s*/(.*?)/",
                  open("/workspace/ForestVegetationSimulator/em/sitset.f").read(), re.S)
    return expand(m.group(1))
bamaxa, sdicon, mapss = grab("BAMAXA"), grab("SDICON"), grab("MAPSS")
mapsdi, mapsit = grab("MAPSDI"), grab("MAPSIT")
assert (len(bamaxa),len(sdicon),len(mapss),len(mapsdi),len(mapsit)) == (30,9,30,122,330)
MAPSIT = lambda it,g: int(mapsit[(g-1)*30 + (it-1)])
# self-check vs live FVSem emt01 (ITYPE 4, IEMTYP 29): SDIDEF = SDICON(MAPSDI(29))
isdi = int(mapsdi[29-1]); assert int(sdicon[isdi-1]) == 696, sdicon[isdi-1]
D = "/workspace/FVSjl/data/easternmontana/"
def wr(fn, hdr, vals): open(D+fn,"w").write(hdr+"\n"+"\n".join(map(str,vals))+"\n")
wr("site_bamaxa.csv","bamax", [int(x) for x in bamaxa])
wr("site_sdicon.csv","sdimax", [int(x) for x in sdicon])
wr("site_mapss.csv","site_species", [int(x) for x in mapss])
wr("site_mapsdi.csv","sdi_group", [int(x) for x in mapsdi])
with open(D+"site_mapsit.csv","w") as f:
    f.write("itype,"+",".join(f"g{g}" for g in range(1,12))+"\n")
    for it in range(1,31):
        f.write(f"{it},"+",".join(str(MAPSIT(it,g)) for g in range(1,12))+"\n")
print("wrote site_bamaxa/sdicon/mapss/mapsdi/mapsit.csv")
print(f"self-check OK: emt01 ITYPE4 IEMTYP29 -> ISDI={isdi} SDIDEF={int(sdicon[isdi-1])} (= live 696), BAMAX={int(bamaxa[3])}, DF SI={MAPSIT(4,8)}")
