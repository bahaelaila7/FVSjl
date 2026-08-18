using FVSjl; const F=FVSjl
hx(x)=uppercase(string(reinterpret(UInt32,x);base=16,pad=8))
fh(h)=reinterpret(Float32,parse(UInt32,h;base=16))
# DFBMOD with OKILL (windthrow feed): MANSTART cyc2, BACHLO from first 56457 uniform + OKILL
d=F.dfb_defaults!(); d.active=true
d.okill=reinterpret(Float32,UInt32(1093555593))   # DF9KIL = 10.895883
k=F.dfb_mod(d, fh("41B78CB8"), fh("42824F2E"), 4, 2)
println("DFKILL(OKILL) jl=",hx(k)," gold=41A9B72E ", hx(k)=="41A9B72E" ? "OK" : "MISMATCH", " val=",k)
