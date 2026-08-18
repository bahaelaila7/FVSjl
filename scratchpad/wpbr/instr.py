#!/usr/bin/env python3
# Inject dump WRITE statements into copies of the wpbr .f routines for dump-replay.
import shutil, os
SRC='/workspace/ForestVegetationSimulator/wpbr'
OVR='/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/wpbr/ovr'
os.makedirs(OVR, exist_ok=True)

def load(name):
    with open(os.path.join(SRC,name)) as f: return f.readlines()
def save(name, lines):
    with open(os.path.join(OVR,name),'w') as f: f.writelines(lines)

def inject_after(lines, anchor, block, occurrence=1):
    out=[]; cnt=0
    for ln in lines:
        out.append(ln)
        if anchor in ln:
            cnt+=1
            if cnt==occurrence:
                out.append(block)
    return out

# ---- BRSETP: dump per-tree init (BRGD, BRHTBC) at end of tree loop ----
b=load('brsetp.f')
blk=("         WRITE(66,9601) ICYC,K,IDTREE(K),ISP(K),HT(K),DBH(K),\n"
     "     &     ICR(K),BRAGE(K),BRGD(K),BRHTBC(K),ISTOTY(K)\n"
     " 9601    FORMAT('SETP ',I3,1X,I5,1X,I6,1X,I3,1X,F12.6,1X,F12.6,1X,\n"
     "     &     I4,1X,F8.2,1X,F14.7,1X,F14.7,1X,I2)\n")
# inject right after BRHTBC assignment line
b=inject_after(b, "BRHTBC(K)=(BRHT-(BRHT*(FLOAT(ICR(K))/100.0)))*100.0", blk)
save('brsetp.f', b)

# ---- BRTARG: dump per-tree GI/RI and stand DFACT/RIAF/AVG ----
t=load('brtarg.f')
# per-tree RI after RI(MM) set (inside loop 35)
blk1=("         WRITE(66,9701) ICYC,MM,GI(MM),RI(MM),TSTARG(MM),ISTOTY(MM)\n"
      " 9701    FORMAT('TARG ',I3,1X,I5,1X,F14.7,1X,F14.7,1X,F14.6,1X,I2)\n")
t=inject_after(t, "RISUM=RISUM+RI(MM)", blk1)
# stand-level after DFACT min reset (before end species loop 90)
blk2=("      WRITE(66,9702) ICYC,I4,RIDEF,RIAF(I4),DFACT(I4,1),AVGRI(I4),\n"
      "     &   AVGGI(I4),AVGSTS(I4),BRNTRECS(I4)\n"
      " 9702 FORMAT('TARGS',I3,1X,I2,1X,F14.9,1X,F12.7,1X,F12.7,1X,F14.9,\n"
      "     &   1X,F12.6,1X,F14.4,1X,I4)\n")
t=inject_after(t, "   90 CONTINUE", blk2)  # after species loop (single WP host)
save('brtarg.f', t)

# ---- BRSTAT: dump THPROB/PITCA per species ----
s=load('brstat.f')
blk=("         WRITE(66,9801) ICYC,I4,THPROB(I4),PITCA(I4),AVTCPT(I4)\n"
     " 9801    FORMAT('STAT ',I3,1X,I2,1X,F12.5,1X,F10.6,1X,F12.6)\n")
s=inject_after(s, "AVECPT(I4)=SUMEC(I4)/THPROB(I4)", blk)
save('brstat.f', s)

# ---- BRECAN: dump each BRANN draw + canker creation ----
e=load('brecan.f')
# entry dump
blk0=("      WRITE(66,9901) ICYC,IBRN,HITE,RI(IBRN),SSTAR,SSTHT,PROP,PIMX,\n"
      "     &   IBRSTAT(IBRN),ILCAN(IBRN),ITCAN(IBRN)\n"
      " 9901 FORMAT('ECAN ',I3,1X,I4,1X,F12.7,1X,F12.9,1X,F14.5,1X,F12.7,\n"
      "     &   1X,F10.7,1X,F10.7,1X,I2,1X,I3,1X,I4)\n")
e=inject_after(e, "      EXPC = 0.0\n", blk0)
# PLI/NUMTIM
blk1=("      WRITE(66,9902) ICYC,IBRN,RITEM,TNEWC,PLI,NUMTIM,CRLEN\n"
      " 9902 FORMAT('ECANP',I3,1X,I4,1X,F12.9,1X,F14.6,1X,F12.9,1X,I5,1X,F14.5)\n")
e=inject_after(e, "      NUMTIM=INT(TNEWC)+1\n", blk1)
# each infection draw
blk2=("         WRITE(66,9903) ICYC,IBRN,J,XBRAN\n"
      " 9903    FORMAT('ECAND',I3,1X,I4,1X,I4,1X,F12.9)\n")
e=inject_after(e, "         CALL BRANN(XBRAN)\n", blk2, occurrence=1)
# canker created: dump TUP/TOUT/PLETH and bole/branch (after DUP/GIRDL set)
blk3=("               WRITE(66,9904) ICYC,IBRN,ICANB,TUP,TOUT,PLETH,\n"
      "     &            DOUT(ICANB,IBRN),XBRAN\n"
      " 9904          FORMAT('ECANC',I3,1X,I4,1X,I3,1X,F12.4,1X,F12.4,1X,\n"
      "     &            F12.7,1X,F12.4,1X,F12.9)\n")
e=inject_after(e, "               ISTCAN(ICANB,IBRN)=0\n", blk3)
save('brecan.f', e)

# ---- BRCGRO: dump kills (WRITE at both kill sites; single shared FORMAT) ----
g=load('brcgro.f')
wblk=("                     WRITE(66,9111) ICYC,K,NCAN,GIRD,UP,WK2(K),\n"
      "     &                  BRPB(K),ISTCAN(NCAN,K)\n")
fblk=(" 9111 FORMAT('CGKIL',I3,1X,I5,1X,I3,1X,F10.4,1X,F12.4,\n"
      "     &   1X,F12.6,1X,F12.6,1X,I2)\n")
g3=[]
for ln in g:
    g3.append(ln)
    if 'IBRSTAT(K)=7' in ln:
        g3.append(wblk)
g3=inject_after(g3, "  500 CONTINUE", fblk)
# per-canker inactivation draw + growth: dump right after XRAN drawn
xblk=("            WRITE(66,9121) ICYC,K,NCAN,OUT,UP,GIRD,XRAN,JCSTAT\n"
      " 9121       FORMAT('CGRAN',I3,1X,I5,1X,I3,1X,F12.4,1X,F12.4,1X,\n"
      "     &         F10.4,1X,F12.9,1X,I2)\n")
g3=inject_after(g3, "            CALL BRANN(XRAN)\n", xblk)
# bole canker growth result: after GIRDL(NCAN,K)=GIRD in bole branch
bblk=("                  WRITE(66,9122) ICYC,K,NCAN,GIRAMT,GROBOL,GIRD\n"
      " 9122             FORMAT('CGBOL',I3,1X,I5,1X,I3,1X,F12.5,1X,F12.6,\n"
      "     &            1X,F10.5)\n")
g3=inject_after(g3, "                  GIRDL(NCAN,K)=GIRD\n", bblk)
save('brcgro.f', g3)

# ---- BRIBA: dump BA + RIDEF (RIMETH>=1 runs) ----
ib=load('briba.f')
blk=("      WRITE(66,9301) ICYC,BA,RIDEF\n"
     " 9301 FORMAT('IBA  ',I3,1X,F12.5,1X,F14.9)\n")
ib=inject_after(ib, "   50 CONTINUE\n", blk)
save('briba.f', ib)

print("instrumented:", os.listdir(OVR))
