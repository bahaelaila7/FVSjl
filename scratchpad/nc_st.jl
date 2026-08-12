using FVSjl
db="/workspace/SQLite_FIADB_ENTIRE.db"; cn="850447807290487"
kt="STDIDENT\n$cn\nNUMCYCLE 2.0\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nECHOSUM\nPROCESS\nSTOP\n"
d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,kt)
FVSjl.run_keyfile(kf; variant=FVSjl.Klamath())
