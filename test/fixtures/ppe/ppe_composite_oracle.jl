# Oracle-baked PPE COMPOSITE YIELD fixture (oracle-free at test time).
#
# Source: the historical FVSppe (Parallel Processing Extension v1.0) landscape driver
# run on the 4-stand East Cascades 15-cycle test (/workspace/.ppework/run/ecpp.{key,tre,out},
# rebuilt from FVS rev bc6e2377^). Weights 1/11/11/11 (TOTAL SAMPLE WEIGHT 34).
#
# `PPE_ORACLE_STANDS` = each stand's (sample_weight, per-year before-thin yield rows) as
# read from the four per-stand "SUMMARY STATISTICS" tables (ecpp.out).
# `PPE_ORACLE_COMPOSITE` = the oracle's "COMPOSITE YIELD STATISTICS" table (ecpp.out
# ~lines 1105-1127) — the landscape aggregate FVSjl's _ppe_aggregate must reproduce.
#
# Aggregation validated (ppbase/cmadds.f + cmprt2.f): TREES/CU FT/MERCH are area-weighted
# means Σ(v·w)/Σw; ACC/MOR are weighted by area·period Σ(v·w·prd)/Σ(w·prd) (equal here
# because every period = 10); PRD = Σ(prd·w)/Σw; sample WEIGHT = Σw. Rounding = IFIX(x+.5).

# (weight, [ (year, tpa, cuft, mcuft, bdft, prd, acc, mort), ... ]) per stand.
const PPE_ORACLE_STANDS = [
    # weight 11
    (11, [
        (year=1990, tpa=536.0, cuft=1624.0, mcuft=1102.0, bdft=5567.0, prd=10.0, acc=108.0, mort=1.0),
        (year=2000, tpa=529.0, cuft=2688.0, mcuft=2167.0, bdft=11066.0, prd=10.0, acc=143.0, mort=2.0),
        (year=2010, tpa=523.0, cuft=4095.0, mcuft=3438.0, bdft=17716.0, prd=10.0, acc=177.0, mort=4.0),
        (year=2020, tpa=514.0, cuft=5829.0, mcuft=5069.0, bdft=26450.0, prd=10.0, acc=190.0, mort=31.0),
        (year=2030, tpa=468.0, cuft=7421.0, mcuft=6641.0, bdft=34923.0, prd=10.0, acc=188.0, mort=31.0),
        (year=2040, tpa=434.0, cuft=8992.0, mcuft=8150.0, bdft=43782.0, prd=10.0, acc=200.0, mort=35.0),
        (year=2050, tpa=405.0, cuft=10641.0, mcuft=9755.0, bdft=53086.0, prd=10.0, acc=201.0, mort=36.0),
        (year=2060, tpa=383.0, cuft=12298.0, mcuft=11576.0, bdft=64400.0, prd=10.0, acc=185.0, mort=85.0),
        (year=2070, tpa=338.0, cuft=13293.0, mcuft=12698.0, bdft=71692.0, prd=10.0, acc=174.0, mort=119.0),
        (year=2080, tpa=288.0, cuft=13840.0, mcuft=13339.0, bdft=76543.0, prd=10.0, acc=157.0, mort=116.0),
        (year=2090, tpa=249.0, cuft=14255.0, mcuft=13764.0, bdft=80312.0, prd=10.0, acc=157.0, mort=123.0),
        (year=2100, tpa=214.0, cuft=14588.0, mcuft=14084.0, bdft=83386.0, prd=10.0, acc=142.0, mort=109.0),
        (year=2110, tpa=188.0, cuft=14923.0, mcuft=14445.0, bdft=86895.0, prd=10.0, acc=142.0, mort=116.0),
        (year=2120, tpa=165.0, cuft=15180.0, mcuft=14720.0, bdft=90015.0, prd=10.0, acc=132.0, mort=113.0),
        (year=2130, tpa=146.0, cuft=15368.0, mcuft=14930.0, bdft=92660.0, prd=10.0, acc=138.0, mort=115.0),
        (year=2140, tpa=129.0, cuft=15598.0, mcuft=15151.0, bdft=95378.0, prd=0.0, acc=0.0, mort=0.0),
    ]),
    # weight 11
    (11, [
        (year=1990, tpa=536.0, cuft=1624.0, mcuft=1102.0, bdft=5567.0, prd=10.0, acc=108.0, mort=1.0),
        (year=2000, tpa=529.0, cuft=2688.0, mcuft=2167.0, bdft=11066.0, prd=10.0, acc=143.0, mort=2.0),
        (year=2010, tpa=523.0, cuft=4095.0, mcuft=3438.0, bdft=17716.0, prd=10.0, acc=119.0, mort=1.0),
        (year=2020, tpa=394.0, cuft=3667.0, mcuft=3150.0, bdft=16552.0, prd=10.0, acc=145.0, mort=2.0),
        (year=2030, tpa=389.0, cuft=5092.0, mcuft=4390.0, bdft=23388.0, prd=10.0, acc=156.0, mort=3.0),
        (year=2040, tpa=385.0, cuft=6624.0, mcuft=5885.0, bdft=31946.0, prd=10.0, acc=118.0, mort=3.0),
        (year=2050, tpa=286.0, cuft=4966.0, mcuft=4546.0, bdft=24687.0, prd=10.0, acc=141.0, mort=3.0),
        (year=2060, tpa=282.0, cuft=6342.0, mcuft=5905.0, bdft=32289.0, prd=10.0, acc=147.0, mort=5.0),
        (year=2070, tpa=278.0, cuft=7763.0, mcuft=7261.0, bdft=40263.0, prd=10.0, acc=91.0, mort=3.0),
        (year=2080, tpa=176.0, cuft=4384.0, mcuft=3915.0, bdft=20933.0, prd=10.0, acc=102.0, mort=4.0),
        (year=2090, tpa=173.0, cuft=5361.0, mcuft=4847.0, bdft=26411.0, prd=10.0, acc=114.0, mort=5.0),
        (year=2100, tpa=170.0, cuft=6450.0, mcuft=6013.0, bdft=33426.0, prd=10.0, acc=80.0, mort=4.0),
        (year=2110, tpa=129.0, cuft=4280.0, mcuft=4060.0, bdft=22409.0, prd=10.0, acc=91.0, mort=5.0),
        (year=2120, tpa=127.0, cuft=5139.0, mcuft=4907.0, bdft=27454.0, prd=10.0, acc=95.0, mort=6.0),
        (year=2130, tpa=124.0, cuft=6029.0, mcuft=5818.0, bdft=32827.0, prd=10.0, acc=71.0, mort=5.0),
        (year=2140, tpa=95.0, cuft=4400.0, mcuft=4213.0, bdft=23530.0, prd=0.0, acc=0.0, mort=0.0),
    ]),
    # weight 11
    (11, [
        (year=1990, tpa=536.0, cuft=1624.0, mcuft=1102.0, bdft=5567.0, prd=10.0, acc=108.0, mort=1.0),
        (year=2000, tpa=529.0, cuft=2688.0, mcuft=2167.0, bdft=11066.0, prd=10.0, acc=143.0, mort=2.0),
        (year=2010, tpa=523.0, cuft=4095.0, mcuft=3438.0, bdft=17716.0, prd=10.0, acc=177.0, mort=4.0),
        (year=2020, tpa=514.0, cuft=5829.0, mcuft=5069.0, bdft=26450.0, prd=10.0, acc=138.0, mort=5.0),
        (year=2030, tpa=140.0, cuft=4928.0, mcuft=4666.0, bdft=24536.0, prd=10.0, acc=151.0, mort=8.0),
        (year=2040, tpa=137.0, cuft=6363.0, mcuft=6084.0, bdft=32690.0, prd=10.0, acc=153.0, mort=10.0),
        (year=2050, tpa=134.0, cuft=7792.0, mcuft=7480.0, bdft=41414.0, prd=10.0, acc=61.0, mort=3.0),
        (year=2060, tpa=31.0, cuft=3733.0, mcuft=3671.0, bdft=22270.0, prd=10.0, acc=56.0, mort=4.0),
        (year=2070, tpa=31.0, cuft=4253.0, mcuft=4186.0, bdft=25944.0, prd=10.0, acc=50.0, mort=4.0),
        (year=2080, tpa=31.0, cuft=4717.0, mcuft=4674.0, bdft=29602.0, prd=10.0, acc=51.0, mort=5.0),
        (year=2090, tpa=31.0, cuft=5179.0, mcuft=5114.0, bdft=32926.0, prd=10.0, acc=48.0, mort=5.0),
        (year=2100, tpa=30.0, cuft=5603.0, mcuft=5556.0, bdft=36446.0, prd=10.0, acc=49.0, mort=6.0),
        (year=2110, tpa=30.0, cuft=6031.0, mcuft=5975.0, bdft=39599.0, prd=10.0, acc=46.0, mort=6.0),
        (year=2120, tpa=30.0, cuft=6432.0, mcuft=6365.0, bdft=42653.0, prd=10.0, acc=45.0, mort=7.0),
        (year=2130, tpa=29.0, cuft=6814.0, mcuft=6743.0, bdft=45639.0, prd=10.0, acc=46.0, mort=7.0),
        (year=2140, tpa=29.0, cuft=7204.0, mcuft=7119.0, bdft=48629.0, prd=0.0, acc=0.0, mort=0.0),
    ]),
    # weight 1
    (1, [
        (year=1990, tpa=0.0, cuft=0.0, mcuft=0.0, bdft=0.0, prd=10.0, acc=0.0, mort=0.0),
        (year=2000, tpa=800.0, cuft=0.0, mcuft=0.0, bdft=0.0, prd=10.0, acc=83.0, mort=0.0),
        (year=2010, tpa=597.0, cuft=834.0, mcuft=0.0, bdft=0.0, prd=10.0, acc=69.0, mort=23.0),
        (year=2020, tpa=329.0, cuft=1285.0, mcuft=997.0, bdft=5043.0, prd=10.0, acc=72.0, mort=30.0),
        (year=2030, tpa=221.0, cuft=1705.0, mcuft=1617.0, bdft=7051.0, prd=10.0, acc=70.0, mort=38.0),
        (year=2040, tpa=158.0, cuft=2030.0, mcuft=1935.0, bdft=8841.0, prd=10.0, acc=62.0, mort=39.0),
        (year=2050, tpa=121.0, cuft=2268.0, mcuft=2201.0, bdft=10090.0, prd=10.0, acc=55.0, mort=34.0),
        (year=2060, tpa=98.0, cuft=2471.0, mcuft=2398.0, bdft=11417.0, prd=10.0, acc=49.0, mort=32.0),
        (year=2070, tpa=83.0, cuft=2646.0, mcuft=2579.0, bdft=12649.0, prd=10.0, acc=45.0, mort=30.0),
        (year=2080, tpa=72.0, cuft=2796.0, mcuft=2750.0, bdft=13681.0, prd=10.0, acc=43.0, mort=30.0),
        (year=2090, tpa=62.0, cuft=2922.0, mcuft=2866.0, bdft=14741.0, prd=10.0, acc=41.0, mort=31.0),
        (year=2100, tpa=55.0, cuft=3024.0, mcuft=2982.0, bdft=15833.0, prd=10.0, acc=41.0, mort=29.0),
        (year=2110, tpa=48.0, cuft=3150.0, mcuft=3111.0, bdft=17057.0, prd=10.0, acc=44.0, mort=31.0),
        (year=2120, tpa=43.0, cuft=3281.0, mcuft=3208.0, bdft=17885.0, prd=10.0, acc=40.0, mort=29.0),
        (year=2130, tpa=38.0, cuft=3390.0, mcuft=3362.0, bdft=19309.0, prd=10.0, acc=33.0, mort=26.0),
        (year=2140, tpa=35.0, cuft=3458.0, mcuft=3432.0, bdft=20060.0, prd=0.0, acc=0.0, mort=0.0),
    ]),
]

# (year, trees, totcuft, mcuft, bdft, acc, mort, prd, weight) — the oracle composite.
const PPE_ORACLE_COMPOSITE = [
    (year=1990, trees=520, cuft=1576, mcuft=1070, bdft=5403, acc=105, mort=1, prd=10, weight=34),
    (year=2000, trees=537, cuft=2609, mcuft=2103, bdft=10741, acc=141, mort=2, prd=10, weight=34),
    (year=2010, trees=525, cuft=3999, mcuft=3337, bdft=17195, acc=155, mort=4, prd=10, weight=34),
    (year=2020, trees=470, cuft=4996, mcuft=4328, bdft=22618, acc=155, mort=13, prd=10, weight=34),
    (year=2030, trees=329, cuft=5693, mcuft=5126, bdft=27011, acc=162, mort=15, prd=10, weight=34),
    (year=2040, trees=314, cuft=7171, mcuft=6566, bdft=35336, acc=154, mort=17, prd=10, weight=34),
    (year=2050, trees=270, cuft=7637, mcuft=7112, bdft=38857, acc=132, mort=15, prd=10, weight=34),
    (year=2060, trees=228, cuft=7311, mcuft=6914, bdft=38823, acc=127, mort=31, prd=10, weight=34),
    (year=2070, trees=212, cuft=8266, mcuft=7887, bdft=44986, acc=103, mort=42, prd=10, weight=34),
    (year=2080, trees=162, cuft=7504, mcuft=7175, bdft=41516, acc=102, mort=41, prd=10, weight=34),
    (year=2090, trees=148, cuft=8108, mcuft=7760, bdft=45614, acc=104, mort=44, prd=10, weight=34),
    (year=2100, trees=136, cuft=8708, mcuft=8387, bdft=50049, acc=89, mort=39, prd=10, weight=34),
    (year=2110, trees=114, cuft=8257, mcuft=8012, bdft=48676, acc=92, mort=42, prd=10, weight=34),
    (year=2120, trees=105, cuft=8751, mcuft=8504, bdft=52330, acc=89, mort=42, prd=10, weight=34),
    (year=2130, trees=98, cuft=9227, mcuft=8993, bdft=55932, acc=83, mort=42, prd=10, weight=34),
    (year=2140, trees=83, cuft=8902, mcuft=8669, bdft=54793, acc=0, mort=0, prd=0, weight=34),
]
