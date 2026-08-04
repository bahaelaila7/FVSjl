#!/usr/bin/env python3
# Extract BC V2-regime DG coefficient + mapping DATA blocks from canada/bc/dgf.f (task #133).
# Handles Fortran N*V repeat syntax. 2D A(m,n) column-major → n rows of m.
import re
src = open('/workspace/ForestVegetationSimulator/canada/bc/dgf.f').read()

def grab(name):
    m = re.search(r'DATA\s+'+name+r'/(.*?)/', src, re.S)
    body = '\n'.join(l for l in m.group(1).split('\n') if not l.strip().startswith(('C','c')))
    vals=[]
    for tok in re.findall(r'\d+\*-?\d+\.?\d*|-?\d+\.?\d*', body):
        if '*' in tok:
            k,v = tok.split('*'); vals += [v]*int(k)
        else:
            vals.append(tok)
    return vals

specs = [('DGHAB',6,15),('DGFOR',6,15),('DGDS',4,15),('DGCCFA',5,15),
         ('DGEL',1,15),('DGEL2',1,15),('DGSASP',1,15),('DGCASP',1,15),('DGSLOP',1,15),('DGSLSQ',1,15),
         ('MAPHAB',30,15),('MAPCCF',30,15),('MAPLOC',11,15),('MAPDSQ',11,15),('OBSERV',6,15)]
for name,m,n in specs:
    v = grab(name)
    assert len(v)==m*n, f"{name}: {len(v)} != {m*n}"
    lines=[f"{j+1},"+",".join(v[j*m:(j+1)*m]) for j in range(n)]
    open(f'data/britishcolumbia/v2dg_{name}.csv','w').write("sp,"+",".join(f"c{i+1}" for i in range(m))+"\n"+"\n".join(lines)+"\n")
print("extracted", len(specs), "blocks OK")
print("sp14 DGHAB=", grab('DGHAB')[13*6:14*6])
print("sp14 MAPHAB=", grab('MAPHAB')[13*30:14*30])
