import struct, numpy as np
def f(i): return struct.unpack('<f', struct.pack('<i', i))[0]
def r(x): return np.float32(x)
vals={'CRASH':1061997773,'THRESH':0,'TELIG':1106692757,'TOTSUC':1058843656,
      'SPCNUM':1088421888,'DF9KIL':1093555593,'OKILL':1093555593}
for k,v in vals.items(): print(k,'=',round(f(v),6),' 0x%08X'%(v&0xffffffff))
CRASH=r(f(1061997773)); TOTSUC=r(f(1058843656)); SPCNUM=r(f(1088421888))
avg=r(TOTSUC/SPCNUM)
for sp,winsuc,prp in [('DF',0.056,1057172953),('GF',0.139,1064514355),('LP',0.028,1048784345)]:
    p=r(r(CRASH*r(winsuc))/avg)
    print(sp,'PRPMRT jl=',round(float(p),7),'gold=',round(f(prp),7),'match=',abs(float(p)-f(prp))<1e-6)
