using FVSjl
const F = FVSjl
hex(x::Float32) = string(reinterpret(UInt32, x); base=16, pad=8) |> uppercase

dbhv = Float32[3.0,5.0,8.0,8.7,9.0,10.0,12.4,15.0,18.0,20.0,24.0,30.0]
tpav = Float32[40.0,30.0,25.0,20.0,18.0,15.0,12.0,10.0,8.0,6.0,4.0,2.0]

# BA accumulation exactly as the Fortran driver (in order, Float32)
ba = 0.0f0
for i in eachindex(dbhv)
    global ba += F.DFB_BAF * dbhv[i] * dbhv[i] * tpav[i]
end
println("BA     ", hex(ba))

r = F.dfb_er(dbhv, tpav, dbhv, tpav, ba)
println("BA9    ", hex(r.ba9))
println("BADF9  ", hex(r.badf9))
println("A45DBH ", hex(r.a45dbh))
println("PBADF4 ", hex(r.pbadf4))
println("LMIN   ", r.lmin)
p = F.dfb_prb(r.a45dbh, r.pbadf4, r.badf9, r.ba9)
println("PROTBK ", hex(p))

st = F.dfb_dbh_start(dbhv, tpav)
for i in 1:20
    println("START $i ", hex(st[i]))
end
testd = Float32[0.4,1.0,1.9,2.0,8.7,9.0,10.0,12.4,19.0,20.0,24.0,30.0,40.0,45.0]
for d in testd
    println("DFBIND ", hex(d), " ", F.dfb_ind(d))
end
