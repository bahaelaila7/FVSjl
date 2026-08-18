using FVSjl; const F=FVSjl
hx(x)=uppercase(string(reinterpret(UInt32,x);base=16,pad=8))
fh(h)=reinterpret(Float32,parse(UInt32,h;base=16))
fd(i)=reinterpret(Float32,UInt32(i&0xffffffff))
# 8 DF records at cyc2 (DFBMRT order). D,P from MRT dump; incoming WK2 = windthrow-or-background;
# final WK2 (gold) from MRT dump.
# (D, P, wk2_in_decimal_or_hex, wk2_final_hex)
recs = [
  ("3DDB85ED","424CABC9", 1102856148, "41BC3FD4"),  # J=2  DCLAS1
  ("40A12BB6","41D2B8B1", 1080096000, "4060F500"),  # J=9  DCLAS3
  ("401B2F71","41CF0B9F", 1080141181, "4061A57D"),  # J=11 DCLAS1
  ("405E1D1E","41DB99B5", 1075730756, "401E5944"),  # J=12 DCLAS2
  ("4147E132","40E3E549", 0x407D5816, "40E3E549"),  # J=15 eligible windthrow WK2
  ("416A2885","408CA304", 0x401E0187, "408CA304"),  # J=17 eligible
  ("41414764","40CF9482", 1053009369, "40AD36DB"),  # J=19 not eligible (background WK2)
  ("4143886C","410157C8", 0x408EFE44, "410157C8"),  # J=22 eligible
]
n=length(recs)
df_dbh=Float32[fh(r[1]) for r in recs]
df_tpa=Float32[fh(r[2]) for r in recs]
wk2in =Float32[fd(r[3]) for r in recs]
old_tpa=copy(df_tpa)
tpa = Float32[old_tpa[i]-wk2in[i] for i in 1:n]
t=(tpa=tpa,)
dfkill=fh("41A9B72E"); badf9=fh("41B78CB8")
F.dfb_mrt!(t, old_tpa, collect(1:n), df_dbh, df_tpa, dfkill, badf9, false, 10.895883f0)
ok=true
for i in 1:n
    expect = old_tpa[i] - fh(recs[i][4])       # surviving TPA the engine stores (avoids cancellation)
    m = t.tpa[i]==expect
    global ok &= m
    println("J$i tpa jl=",hx(t.tpa[i])," expect=",hx(expect)," WK2gold=",recs[i][4]," ", m ? "OK" : "X")
end
println(ok ? "ALL-OK" : "FAIL")
