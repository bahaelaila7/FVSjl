using FVSjl, SQLite, DBInterface
using FVSjl: run_keyfile, Kootenai
v = Kootenai()
dir = mktempdir()
# same tree file as the oracle rdsum.key (S248112)
cp("/workspace/.ktwork/rdrun/rdsum.tre", joinpath(dir,"rdab.tre"))
dbp = joinpath(dir,"rdab.db"); key = joinpath(dir,"rdab.key")
write(key, """
STDIDENT
S248112  RD SUM
DATABASE
DSNOUT
$dbp
RDSUM
END
SCREEN
NOAUTOES
NOTRIPLE
STATS
DESIGN                                        11.0       1.0
STDINFO     11406001     570.0      60.0     315.0      30.0      34.0
INVYEAR       1990.0
NUMCYCLE        10.0
TREEFMT
(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,
T52,I2,T66,5I1,T54,7I1,T75,F3.0)
TREEDATA
RDIN
RRTYPE             3
RRINIT             0        10        10        20       0.1        10         3
SAREA            100
RRDOUT
BBCLEAR
END
ECHOSUM
PROCESS
STOP
""")
run_keyfile(key; variant=v)
jl = SQLite.DB(dbp); orc = SQLite.DB("/workspace/.ktwork/rdrun/rdsum_oracle.db")
jr = [NamedTuple(r) for r in DBInterface.execute(jl,"SELECT Year,Age,RD_Type,Num_Centers,RD_Area,Spread_Ft_per_Year,Stumps_per_Acre,Stumps_BA,Mort_TPA,Mort_CuFt,UnInf_TPA,Inf_TPA,Live_Merch_CuFt,Live_BA FROM FVS_RD_Sum ORDER BY Year")]
orr= [NamedTuple(r) for r in DBInterface.execute(orc,"SELECT Year,Age,RD_Type,Num_Centers,RD_Area,Spread_Ft_per_Year,Stumps_per_Acre,Stumps_BA,Mort_TPA,Mort_CuFt,UnInf_TPA,Inf_TPA,Live_Merch_CuFt,Live_BA FROM FVS_RD_Sum ORDER BY Year")]
println("jl rows: ", length(jr), "  oracle rows: ", length(orr))
for (a,b) in zip(jr,orr)
  println("YR ",a.Year," | jl  Age=",a.Age," T=",a.RD_Type," NC=",a.Num_Centers," Area=",round(a.RD_Area,digits=3)," Spr=",round(a.Spread_Ft_per_Year,digits=3)," Stmp=",round(a.Stumps_per_Acre,digits=3)," Mort=",round(a.Mort_TPA,digits=2)," Inf=",round(a.Inf_TPA,digits=2)," BA=",round(a.Live_BA,digits=2))
  println("      | orc Age=",b.Age," T=",b.RD_Type," NC=",b.Num_Centers," Area=",round(b.RD_Area,digits=3)," Spr=",round(b.Spread_Ft_per_Year,digits=3)," Stmp=",round(b.Stumps_per_Acre,digits=3)," Mort=",round(b.Mort_TPA,digits=2)," Inf=",round(b.Inf_TPA,digits=2)," BA=",round(b.Live_BA,digits=2))
end
SQLite.close(jl); SQLite.close(orc)
