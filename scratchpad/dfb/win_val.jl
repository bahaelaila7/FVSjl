using FVSjl; const F=FVSjl
include("/tmp/win_stand.jl")
f(i)=reinterpret(Float32,UInt32(i&0xffffffff))
hx(x)=uppercase(string(reinterpret(UInt32,x);base=16,pad=8))
n=length(WIN_STAND)
species=Int[r[2] for r in WIN_STAND]
pct=Float32[f(r[3]) for r in WIN_STAND]
ht =Float32[f(r[4]) for r in WIN_STAND]
dbh=Float32[f(r[5]) for r in WIN_STAND]
prob=Float32[f(r[6]) for r in WIN_STAND]
wk2=Float32[f(r[7]) for r in WIN_STAND]
df9kil,telig=F.dfb_win_kernel!(species,dbh,ht,prob,pct,wk2,F._DFB_IFVSSP_IE,80.0f0,20.0f0,0.8f0,0.0f0,3,23)
gold=Dict(15=>1081956374,17=>1075708295,22=>1083113028,21=>1086232353,1=>1071438053)
ok=true
for (I,gw) in sort(collect(gold))
    idx=findfirst(x->x[1]==I, WIN_STAND)
    m = wk2[idx]==f(gw)
    global ok &= m
    println("WK2[I=$I] jl=",hx(wk2[idx])," gold=",hx(f(gw))," ", m ? "OK" : "X")
end
tg=telig==f(1106692757); dg=df9kil==f(1093555593)
println("TELIG  jl=",hx(telig)," gold=",hx(f(1106692757))," ", tg ? "OK" : "X")
println("DF9KIL jl=",hx(df9kil)," gold=",hx(f(1093555593))," ", dg ? "OK" : "X")
println(ok && tg && dg ? "ALL-OK" : "FAIL")
