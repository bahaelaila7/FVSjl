import SQLite, DBInterface
db = SQLite.DB("file:/workspace/SQLite_FIADB_ENTIRE.db?mode=ro", )
for v in ARGS
  open("/workspace/.sweep_east/pop_$(lowercase(v)).txt","w") do io
    for r in DBInterface.execute(db, "SELECT STAND_CN FROM FVS_STANDINIT_COND WHERE VARIANT='$v' AND STAND_CN IS NOT NULL ORDER BY ECOREGION, LOCATION, STAND_CN")
      println(io, string(r.STAND_CN))
    end
  end
end
