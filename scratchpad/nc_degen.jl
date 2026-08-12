using FVSjl, SQLite, DBInterface
db="/workspace/SQLite_FIADB_ENTIRE.db"; DB=SQLite.DB(db); cn="850447807290487"
# DB raw tree rows (species, DBH/DRC)
println("=== DB FVS_TREEINIT_COND (first 8 rows) ===")
tcols=[r.name for r in DBInterface.execute(DB,"PRAGMA table_info(FVS_TREEINIT_COND)")]
dcol=filter(c->occursin(r"^DIAMETER$|^DBH$|^DG$|DRC|HISTORY|TREE_COUNT|SPECIES"i,c),tcols)
println("cols of interest: ", dcol)
q="SELECT SPECIES,DIAMETER,DG,HT,TREE_COUNT FROM FVS_TREEINIT_COND WHERE STAND_CN='$cn' LIMIT 8"
for r in DBInterface.execute(DB,q); println("  sp=$(r.SPECIES) D=$(r.DIAMETER) DG=$(r.DG) HT=$(r.HT) TC=$(r.TREE_COUNT)"); end
# jl initial trees
kt="STDIDENT\n$cn\nNUMCYCLE 1.0\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nECHOSUM\nPROCESS\nSTOP\n"
d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,kt)
println("=== jl initial trees (first 8) ===")
for s in FVSjl.each_stand(kf; variant=FVSjl.Klamath())
    t=s.trees
    for i in 1:min(8,t.n); println("  sp=$(t.species[i]) D=$(round(t.dbh[i],digits=2)) H=$(round(t.height[i],digits=1)) tpa=$(round(t.tpa[i],digits=2))"); end
    println("  n=$(t.n) ΣBA=$(round(sum(0.005454f0*t.dbh[i]^2*t.tpa[i] for i in 1:t.n),digits=1))")
    break
end
