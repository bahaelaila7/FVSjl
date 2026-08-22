using FVSjl, SQLite, DBInterface
v = FVSjl.InlandEmpire()
dir = mktempdir()
tre = """
   1      248112       0101   011LP 11510   0734   00111     0  0
   2      248112       0101   031DF 001     0026   00222     0  0
   9      248112       0103   011LP 09511   0603   00111     0  0
  17      248112       0106   011DF 10010   0654   00111     0  0
"""
write(joinpath(dir,"bm.tre"), tre)
dbp = joinpath(dir,"bm.db")
key = joinpath(dir,"bm.key")
write(key, """
SCREEN
NOAUTOES
NOTRIPLE
STATS
STDIDENT
S248112  WWPB PPBMMAIN
DESIGN                                        11.0       1.0
STDINFO     11406001     570.0      60.0     315.0      30.0      34.0
INVYEAR       1990.0
NUMCYCLE         5.0
TREEFMT
(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,
T52,I2,T66,5I1,T54,7I1,T75,F3.0)
TREEDATA
BMIN
DISPERSE         1.0       5.0       3.0      20.0
END
DATABASE
DSNOUT
$dbp
PPBMMAIN
END
PROCESS
STOP
""")
run_keyfile(key; variant=v)
db = SQLite.DB(dbp)
tabs=[r.name for r in DBInterface.execute(db,"SELECT name FROM sqlite_master WHERE type='table' AND name='FVS_BM_Main'")]
println("FVS_BM_Main present: ", "FVS_BM_Main" in tabs)
if "FVS_BM_Main" in tabs
  n=first(DBInterface.execute(db,"SELECT COUNT(*) c FROM FVS_BM_Main")).c
  println("rows: ", n)
  for r in DBInterface.execute(db,"SELECT Year,StandBA,TPA,TPA_BtlKld,BA_BtlKld,StandVol,VolBtlKld,PreDispBKP FROM FVS_BM_Main")
    println("  yr=",r.Year," BA=",round(r.StandBA,digits=2)," TPA=",round(r.TPA,digits=1)," TPAkld=",round(r.TPA_BtlKld,digits=2)," BAkld=",round(r.BA_BtlKld,digits=3)," Vol=",round(r.StandVol,digits=1)," Volkld=",round(r.VolBtlKld,digits=2)," preBKP=",round(r.PreDispBKP,digits=3))
  end
end
SQLite.close(db)
