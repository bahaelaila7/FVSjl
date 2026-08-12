using FVSjl
db="/workspace/SQLite_FIADB_ENTIRE.db"; cn="850447807290487"
kt="STDIDENT\n$cn\nNUMCYCLE 2.0\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nECHOSUM\nPROCESS\nSTOP\n"
d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,kt)
for s in FVSjl.each_stand(kf; variant=FVSjl.Klamath())
    FVSjl.notre!(s); FVSjl.setup_growth!(s)   # mirror the real run's pre-growth (height estimation, etc.)
    t=s.trees
    println("AFTER setup_growth: n=$(t.n)")
    for i in 1:t.n; println("  sp=$(t.species[i]) D=$(round(t.dbh[i],digits=3)) H=$(round(t.height[i],digits=2)) tpa=$(round(t.tpa[i],digits=1)) crown=$(t.crown_pct[i])"); end
    # site index for this stand
    println("  site_idx(sp1)=$(s.plot.sp_site_index[1]) forest=$(s.plot.forest_idx)")
    break
end
