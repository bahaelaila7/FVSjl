#!/usr/bin/env python3
# Extract BC Kozak-2002 taper coefficients (canada/bc/log.f) — the IZ=2 (FIZ=2, zones D-J) column only,
# which is the sole column BC uses (min.f IFZ=int(FIZ,2), FIZ=2 fixed in cfvol.f). Each A_k/PP is DIMENSION
# (16 Kozak-species, 3 zone-groups), Fortran column-major DATA ⇒ values 1-16 = IZ1, 17-32 = IZ2, 33-48 = IZ3.
import re, sys
src = open('/workspace/ForestVegetationSimulator/canada/bc/log.f').read()

def grab(name):
    # find 'DATA <name>/ ... /'  (values, ignoring continuation markers + comments)
    m = re.search(r'DATA\s+'+name+r'/(.*?)/', src, re.S)
    body = m.group(1)
    # strip Fortran comment lines
    body = '\n'.join(l for l in body.split('\n') if not l.strip().startswith(('C','c')))
    nums = re.findall(r'-?\d+\.\d+', body)
    vals = [float(x) for x in nums]
    assert len(vals) == 48, f"{name}: got {len(vals)}"
    return vals[16:32]   # IZ=2 middle column

arrays = {k: grab(k) for k in ['A1','A2','A3','A4','A5','A6','A7','A8','PP']}
# SPTR (min.f) — BC species (1..15) -> Kozak taper species (1..16)
sptr = [7,10,1,4,3, 2,8,5,4,9, 14,15,11,1,14]

hdr = "kozak_sp,A1,A2,A3,A4,A5,A6,A7,A8,PP"
lines = [hdr]
for i in range(16):
    row = [str(i+1)] + [f"{arrays[k][i]:.6f}" for k in ['A1','A2','A3','A4','A5','A6','A7','A8','PP']]
    lines.append(','.join(row))
open('/workspace/FVSjl/data/britishcolumbia/vol_taper_iz2.csv','w').write('\n'.join(lines)+'\n')
open('/workspace/FVSjl/data/britishcolumbia/vol_sptr.csv','w').write('bc_sp,kozak_sp\n'+'\n'.join(f"{i+1},{s}" for i,s in enumerate(sptr))+'\n')
print("wrote vol_taper_iz2.csv (16 Kozak sp x 9 coeffs, IZ=2) + vol_sptr.csv")
print("sp14 -> kozak", sptr[13], " A1..A8/PP:", [f"{arrays[k][sptr[13]-1]:.5f}" for k in ['A1','A2','A3','A4','A5','A6','A7','A8','PP']])
