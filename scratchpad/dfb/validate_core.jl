using FVSjl
const F = FVSjl
hx(x::Float32) = uppercase(string(reinterpret(UInt32, x); base=16, pad=8))
fromhex(h) = reinterpret(Float32, parse(UInt32, h; base=16))

d = F.dfb_defaults!()
d.active = true

println("== DFBRAN stream from ORSEED+1128 = 56457 ==")
gold = ["3EE23A99","3E6A6620","3E56B516"]
draws = Float32[]
for i in 1:3
    u = F.dfb_rand!(d)
    push!(draws, u)
    println("draw $i  jl=$(hx(u))  gold=$(gold[i])  ", hx(u)==gold[i] ? "OK" : "MISMATCH")
end

# reset rng, validate BACHLO + DFBMOD DFKILL
d2 = F.dfb_defaults!(); d2.active = true
badf9 = fromhex("41B78CB8"); ba9 = fromhex("42824F2E")
dfkill = F.dfb_mod(d2, badf9, ba9, 4)
println("\n== DFBMOD DFKILL ==")
println("jl=$(hx(dfkill))  gold=412518D2  ", hx(dfkill)=="412518D2" ? "OK" : "MISMATCH")
println("dfkill value = ", dfkill)

# per-record TAMORT (DBH method), SUMDBH from all DF>=9
println("\n== DFBMRT per-record TAMORT ==")
sumdbh_gold = fromhex("43A53859")
dfrec = [  # (D, P, DCLAS, TAMORT_gold)
 ("3DDB85ED","424CABC9",1,nothing),
 ("40A12BB6","41D2B8B1",3,nothing),
 ("401B2F71","41CF0B9F",1,nothing),
 ("405E1D1E","41DB99B5",2,nothing),
 ("4147E132","40E3E549",6,"4031CDC0"),
 ("416A2885","408CA304",7,"40008AB8"),
 ("41414764","40CF9482",6,"401C9AF1"),
 ("4143886C","410157C8",6,"40456FE1"),
]
df_dbh = Float32[fromhex(r[1]) for r in dfrec]
df_tpa = Float32[fromhex(r[2]) for r in dfrec]
sumdbh = let acc=0.0f0
    for k in eachindex(df_dbh); df_dbh[k]>=9f0 && (acc += df_tpa[k]*df_dbh[k]); end
    acc
end
println("SUMDBH jl=$(hx(sumdbh)) gold=$(hx(sumdbh_gold))  ", hx(sumdbh)==hx(sumdbh_gold) ? "OK" : "MISMATCH")
for (k,r) in enumerate(dfrec)
    F.dfb_ind(df_dbh[k]) < 5 && continue
    tamort = dfkill * (df_dbh[k]*df_tpa[k]/sumdbh)
    tamort > df_tpa[k] && (tamort = df_tpa[k])
    println("rec$k D=$(r[1]) tamort jl=$(hx(tamort)) gold=$(r[4])  ", hx(tamort)==r[4] ? "OK" : "MISMATCH")
end
