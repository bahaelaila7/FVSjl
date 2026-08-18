#!/usr/bin/env python3
# Re-instrument the wpbr routines to dump REAL*4 values in HEX (Z9 edit descriptor)
# so the goldens round-trip to the exact IEEE Float32 bit pattern (F-format decimals
# are under-precise for large magnitudes → 1-ULP artifacts). Inputs AND outputs are
# dumped in hex so the Julia replay feeds bit-exact inputs.
import os
SRC='/workspace/ForestVegetationSimulator/wpbr'
OVR='/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/wpbr/ovr'
os.makedirs(OVR, exist_ok=True)
def load(n):
    with open(os.path.join(SRC,n)) as f: return f.readlines()
def save(n,l):
    with open(os.path.join(OVR,n),'w') as f: f.writelines(l)
def after(lines, anchor, block, occ=1):
    out=[]; c=0
    for ln in lines:
        out.append(ln)
        if anchor in ln:
            c+=1
            if c==occ: out.append(block)
    return out

# BRSETP: BRGD, BRHTBC in hex
b=load('brsetp.f')
blk=("         WRITE(66,7601) K,HT(K),DBH(K),ICR(K),BRGD(K),BRHTBC(K)\n"
     " 7601    FORMAT('HSETP ',I5,1X,Z8,1X,Z8,1X,I4,1X,Z8,1X,Z8)\n")
b=after(b, "BRHTBC(K)=(BRHT-(BRHT*(FLOAT(ICR(K))/100.0)))*100.0", blk)
save('brsetp.f', b)

# BRTARG: per-tree GI, TSTARG in hex (loop 30 sets GI(M)/TSTARG(M))
t=load('brtarg.f')
blk=("         WRITE(66,7701) M,GI(M),TSTARG(M)\n"
     " 7701    FORMAT('HTARG ',I5,1X,Z8,1X,Z8)\n")
t=after(t, "         TSTARG(M)=TBSUM\n", blk)
save('brtarg.f', t)

# BRIBA: BA (input) + RIDEF (output) in hex
ib=load('briba.f')
blk=("      WRITE(66,7301) ICYC,BA,RIDEF\n"
     " 7301 FORMAT('HIBA  ',I3,1X,Z8,1X,Z8)\n")
ib=after(ib, "   50 CONTINUE\n", blk)
save('briba.f', ib)

# BRECAN: HITE,SSTAR (inputs) + RITEM,TNEWC,PLI (outputs) in hex; NUMTIM int
e=load('brecan.f')
blk=("      WRITE(66,7901) ICYC,IBRN,HITE,SSTAR,RITEM,TNEWC,PLI,NUMTIM\n"
     " 7901 FORMAT('HECAN ',I3,1X,I4,1X,Z8,1X,Z8,1X,Z8,1X,Z8,1X,Z8,1X,I5)\n")
e=after(e, "      NUMTIM=INT(TNEWC)+1\n", blk)
# placement: SSTHT,CRLEN,TUP (inputs) + TOUT,PLETH (outputs) in hex
blk2=("            WRITE(66,7902) ICYC,IBRN,SSTHT,CRLEN,TUP,TOUT,PLETH\n"
      " 7902       FORMAT('HPLAC ',I3,1X,I4,1X,Z8,1X,Z8,1X,Z8,1X,Z8,1X,Z8)\n")
e=after(e, "            IF(PLETH.LT.0.0) PLETH=0.0\n", blk2)
save('brecan.f', e)

# BRCGRO: bole growth GIRAMT,GROBOL (inputs) + GIRD (output) in hex, plus gird-in
g=load('brcgro.f')
# capture GIRD before update? The bole branch: GIRD=GIRD+... We dump giramt,grobol,gird(new)
blk=("                  WRITE(66,7122) ICYC,K,NCAN,GIRAMT,GROBOL,GIRD\n"
     " 7122             FORMAT('HCGBO ',I3,1X,I5,1X,I3,1X,Z8,1X,Z8,1X,Z8)\n")
g=after(g, "                  GIRDL(NCAN,K)=GIRD\n", blk)
save('brcgro.f', g)
print("hex-instrumented:", sorted(os.listdir(OVR)))
