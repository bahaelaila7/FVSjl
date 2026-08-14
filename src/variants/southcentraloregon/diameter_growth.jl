# =============================================================================
# diameter_growth.jl (southcentraloregon) — SO large-tree DDS (so/dgf.f). Chunk 3.
#
# SO's dgf.f is NOT group-compressed: every coefficient array is length MAXSP=33, indexed directly by
# ISPC (identity map). Relative to EC it carries the FULL Prognosis term set — DGMAI·RMAI and
# DGSIC·XSITE in DGCON, plus DGCCFA·RELDEN + DGMACC·RMAI·RELDEN in CONSPP, and DGBA/DGHAH/DGLBA/DGLCCF
# (sp10 only) in the ln(DDS) body. Four DDS branches:
#   • GENERAL Wykoff ln(DDS): most species (so/dgf.f:693-706)
#   • JU(11): western-juniper DF-projection (UT-variant form, so/dgf.f:568-586)
#   • AS(24): aspen DGFASP (UT-variant form, so/dgf.f:591-596) → shared `_em_dgfasp`
#   • RA(22): red alder decreasing-increment eqn (byte-identical to WC/EC/PN, so/dgf.f:603-616)
#   Warm Springs (IFOR==10) borrows EC/WC eqns for sp {1,2,6,7,18,32}; sot01 forest 601 = IFOR1 ⇒ inert,
#   scaffolded but not yet cyc-exact.
# RMAI (so/maical.f) is a stand constant (ADJMAI, capped 128) computed in `so_dgcons!` and cached in
# p.mai_adj so both DGCON and the CONSPP/dgf! body read it. MEASURED bit-exact vs FVSso_g16 dgf per-tree
# LN(DDS) on sot01 (all 33 species: general + JU + AS + RA branches).
# =============================================================================

# ── SO per-species scalar coefficients (so/dgf.f DATA), length 33, index = ISPC.
const SO_DGLD   = Float32[0.77889,1.05245,0.33160,0.89887,0.99324,0.72947,0.63249,1.32610,1.186676,0.53028,0.0,0.89887,0.91259,0.980383,0.9042530,0.213947,0.609098,0.58705,0.7224620,0.8793380,0.8895960,0.0,1.0241860,0.0,0.8895960,0.8895960,1.310111,0.8895960,0.8895960,0.8895960,0.8895960,0.33160,0.8895960]
const SO_DGCR   = Float32[3.36606,2.89738,4.43317,3.11044,1.73837,1.57139,2.35065,1.29730,2.763519,3.43377,0.0,3.11044,2.44945,1.709846,4.1231012,1.523464,1.158355,1.29360,2.1603479,1.9700520,1.7325350,0.0,0.4593870,0.0,1.7325350,1.7325350,0.271183,1.7325350,1.7325350,1.7325350,1.7325350,4.43317,1.7325350]
const SO_DGCRSQ = Float32[-1.80146,-0.95430,-2.82894,-1.13806,-0.12161,0.0,-0.63173,0.0,-0.871061,-1.84300,0.0,-1.13806,-0.72173,0.0,-2.6893401,0.0,0.0,0.0,-0.8341960,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-2.82894,0.0]
const SO_DGBA   = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.0001730,-0.000981,0.0,0.0,0.0,-0.0009810,-0.000981,0.0,-0.000981,-0.000981,-0.0009810,-0.000981,0.0,-0.000981]
const SO_DGBAL  = Float32[0.00121,0.00020,0.00479,0.00106,-0.00087,0.00483,0.00169,0.0,0.0,0.00325,0.0,0.00106,0.00380,0.0,0.0,-0.358634,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.00479,0.0]
const SO_DGDBAL = Float32[-0.00897,-0.00437,-0.01918,-0.00707,-0.00105,-0.01497,-0.00647,-0.00239,-0.003728,-0.01362,0.0,-0.00707,-0.01620,-0.000261,-0.006368,0.0,-0.004253,-0.02284,-0.004065,-0.004215,-0.001265,0.0,-0.010222,0.0,-0.001265,-0.001265,0.0,-0.001265,-0.001265,-0.001265,-0.001265,-0.01918,-0.001265]
const SO_DGMAI  = Float32[0.0,0.00375,0.00044,0.00114,0.0,0.00246,0.00778,0.0,0.0,0.00360,0.0,0.00114,0.00753,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.00044,0.0]
const SO_DGMACC = Float32[0.00001,-0.00001,0.00001,0.00001,0.00001,0.0,-0.00003,0.0,0.0,0.0,0.0,0.00001,-0.00002,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.00001,0.0]
const SO_DGPCCF = Float32[0.0,-0.00060,-0.00041,-0.00027,0.0,-0.00053,0.0,-0.00044,0.0,0.0,0.0,-0.00027,-0.00064,-0.000643,-0.0004710,0.0,-0.000568,-0.00094,0.0,0.0,0.0,0.0,-0.000757,0.0,0.0,0.0,-0.000473,0.0,0.0,0.0,0.0,-0.00041,0.0]
const SO_DGCCFA = Float32[-0.00016,0.00157,-0.00138,-0.00024,-0.00027,0.0,0.00109,0.0,0.0,0.0,0.0,-0.00024,0.00201,0.0,0.0,-0.199592,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.00138,0.0]
const SO_OBSERV = Float32[256.0,306.0,697.0,5513.0,779.0,501.0,3208.0,1185.0,2062.,7837.0,22.0,5513.0,1198.0,1210.,1467.,22.0,591.,100.0,4836.,475.,220.0,125.,78.0,22.0,220.,220.,306.,220.,220.,220.,220.,697.0,220.]
const SO_DGCASP = Float32[0.12915,0.08831,-0.17389,-0.14234,-0.29385,-0.04156,-0.39893,0.38002,-0.444594,-0.05217,0.0,-0.14234,-0.65502,-0.059062,-0.3745120,-0.609774,-0.156235,-0.06625,0.0,0.0,0.0859580,0.0,0.0,0.0,0.0859580,0.0859580,0.0,0.0859580,0.0859580,0.0859580,0.0859580,-0.17389,0.0859580]
const SO_DGSASP = Float32[-0.19278,0.56358,-0.40153,-0.17130,0.33069,-0.43281,-0.25658,-0.17911,0.139180,-0.14076,0.0,-0.17130,-0.28887,-0.128126,-0.2076590,-0.017520,0.258712,0.05534,0.0,0.0,-0.8639800,0.0,0.0,0.0,-0.86398,-0.86398,0.0,-0.86398,-0.86398,-0.86398,-0.86398,-0.40153,-0.86398]
const SO_DGSLOP = Float32[0.77922,-1.50252,-0.18923,0.02912,-0.59628,-0.90318,0.19253,-0.81780,0.0,-0.29407,0.0,0.02912,0.64936,0.240178,0.4002230,-2.057060,-0.635704,0.11931,0.4214860,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.18923,0.0]
const SO_DGSLSQ = Float32[-0.93813,1.21439,-0.49057,-0.29671,1.07549,1.36603,0.0,0.84368,0.0,0.16735,0.0,-0.29671,-0.37153,0.131356,0.0,2.11326,0.0,0.0,-0.6936100,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.49057,0.0]
const SO_DGEL   = Float32[0.00279,0.00466,0.01544,0.00362,-0.03036,-0.00637,0.01012,0.0,0.0248,-0.00331,0.0,0.00362,0.02956,-0.015087,-0.0690450,0.0,0.004379,-0.00175,-0.0400670,0.0,-0.075986,0.0,-0.012111,0.0,-0.0759860,-0.075986,0.0049,-0.075986,-0.075986,-0.075986,-0.075986,0.01544,-0.075986]
const SO_DGEL2  = Float32[-0.00001,0.00001,-0.00010,-0.00006,0.00037,0.00014,-0.00020,0.0,-0.00033429,0.00006,0.0,-0.00006,-0.00033,0.0,0.0006080,0.0,0.0,-0.000067,0.0003950,0.0,0.0011930,0.0,0.0,0.0,0.001193,0.001193,-0.00008781,0.001193,0.001193,0.001193,0.001193,-0.00010,0.001193]
const SO_DGSIC  = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.001766,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const SO_DGSITE = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.492695,0.0,0.0,0.0,0.0,0.323625,0.6849390,0.0,0.351929,0.0,0.3804160,0.2528530,0.2273070,0.0,1.965888,0.0,0.227307,0.227307,0.213526,0.227307,0.227307,0.227307,0.227307,0.0,0.227307]
const SO_DGLBA  = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.122905,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const SO_DGHAH  = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.0003580,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const SO_SLODUM = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.174404,0.0,0.0,-0.290174,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]

# DGFOR[isp, isfor] (so/dgf.f DATA DGFOR(6,33) → jl [isp,isfor]; ISPFOR∈1..6 via MAPLOC).
const SO_DGFOR = Float32[
   -0.23185   0.0       0.0        0.0       0.0      0.0;      #  1 WP
   -1.09124  -1.04630  -0.048564  -0.47130   0.0      0.0;      #  2 SP
    0.83941   0.97742   0.67207    1.17949   0.83941  0.0;      #  3 DF
   -0.16674  -0.04804   0.02336    0.47563   0.0      0.0;      #  4 WF
   -0.25368   0.0       0.0        0.0       0.0      0.0;      #  5 MH
    0.23477   0.0       0.0        0.0       0.0      0.0;      #  6 IC
   -0.03997  -0.17367   0.80024    0.40611   0.17307  0.35268;  #  7 LP
   -4.64535   0.0       0.0        0.0       0.0      0.0;      #  8 ES
   -2.073942  0.0       0.0        0.0       0.0      0.0;      #  9 SH
    1.48355   1.34282   1.46037    1.54377   0.0      0.0;      # 10 PP
    0.0       0.0       0.0        0.0       0.0      0.0;      # 11 JU
   -0.16674  -0.04804   0.02336    0.47563   0.0      0.0;      # 12 GF
   -1.68772  -0.96240  -1.38122   -1.03440  -0.76742  0.0;      # 13 AF
   -0.441408  0.0       0.0        0.0       0.0      0.0;      # 14 SF
   -1.127977  0.0       0.0        0.0       0.0      0.0;      # 15 NF
    1.911884  0.0       0.0        0.0       0.0      0.0;      # 16 WB
   -0.605649  0.0       0.0        0.0       0.0      0.0;      # 17 WL
    1.49419   0.0       0.0        0.0       0.0      0.0;      # 18 RC
   -0.147675  0.0       0.0        0.0       0.0      0.0;      # 19 WH
   -1.3100671 0.0       0.0        0.0       0.0      0.0;      # 20 PY
   -0.107648  0.0       0.0        0.0       0.0      0.0;      # 21 WA
    0.0       0.0       0.0        0.0       0.0      0.0;      # 22 RA
   -7.753469  0.0       0.0        0.0       0.0      0.0;      # 23 BM
    0.0       0.0       0.0        0.0       0.0      0.0;      # 24 AS
   -0.107648  0.0       0.0        0.0       0.0      0.0;      # 25 CW
   -0.107648  0.0       0.0        0.0       0.0      0.0;      # 26 CH
   -1.958189  0.0       0.0        0.0       0.0      0.0;      # 27 WO
   -0.107648  0.0       0.0        0.0       0.0      0.0;      # 28 WI
   -0.107648  0.0       0.0        0.0       0.0      0.0;      # 29 GC
   -0.107648  0.0       0.0        0.0       0.0      0.0;      # 30 MC
   -0.107648  0.0       0.0        0.0       0.0      0.0;      # 31 MB
    0.83941   0.97742   0.67207    1.17949   0.83941  0.0;      # 32 OS
   -0.107648  0.0       0.0        0.0       0.0      0.0]      # 33 OH

# MAPLOC[isp, ifor] (so/dgf.f DATA MAPLOC(10,33) → jl [isp,ifor]; loc class → DGFOR col, IFOR 1..10).
const SO_MAPLOC = Int[
    1 1 1 1 1 1 1 1 1 1;   #  1
    1 2 3 4 4 4 4 1 1 1;   #  2
    1 1 2 3 4 5 4 5 1 1;   #  3
    1 1 2 2 3 3 3 4 1 1;   #  4
    1 1 1 1 1 1 1 1 1 1;   #  5
    1 1 1 1 1 1 1 1 1 1;   #  6
    1 2 2 3 4 5 4 6 1 1;   #  7
    1 1 1 1 1 1 1 1 1 1;   #  8
    1 1 1 1 1 1 1 1 1 1;   #  9
    1 2 2 3 3 4 3 2 1 1;   # 10
    1 1 1 1 1 1 1 1 1 1;   # 11
    1 1 2 2 3 3 3 4 1 1;   # 12
    1 2 2 3 4 5 4 5 1 1;   # 13
    1 1 1 1 1 1 1 1 1 1;   # 14
    1 1 1 1 1 1 1 1 1 1;   # 15
    1 1 1 1 1 1 1 1 1 1;   # 16
    1 1 1 1 1 1 1 1 1 1;   # 17
    1 1 1 1 1 1 1 1 1 1;   # 18
    1 1 1 1 1 1 1 1 1 1;   # 19
    1 1 1 1 1 1 1 1 1 1;   # 20
    1 1 1 1 1 1 1 1 1 1;   # 21
    1 1 1 1 1 1 1 1 1 1;   # 22
    1 1 1 1 1 1 1 1 1 1;   # 23
    1 1 1 1 1 1 1 1 1 1;   # 24
    1 1 1 1 1 1 1 1 1 1;   # 25
    1 1 1 1 1 1 1 1 1 1;   # 26
    1 1 1 1 1 1 1 1 1 1;   # 27
    1 1 1 1 1 1 1 1 1 1;   # 28
    1 1 1 1 1 1 1 1 1 1;   # 29
    1 1 1 1 1 1 1 1 1 1;   # 30
    1 1 1 1 1 1 1 1 1 1;   # 31
    1 1 2 3 4 5 4 5 1 1;   # 32
    1 1 1 1 1 1 1 1 1 1]   # 33

# MAPDSQ[isp, ifor] (so/dgf.f DATA MAPDSQ(10,33)) — DGDS column selector; all ones except sp 4/6/10/12.
const SO_MAPDSQ = Int[
    1 1 1 1 1 1 1 1 1 1;   #  1
    1 1 1 1 1 1 1 1 1 1;   #  2
    1 1 1 1 1 1 1 1 1 1;   #  3
    1 1 1 2 3 2 3 1 1 1;   #  4
    1 1 1 1 1 1 1 1 1 1;   #  5
    1 1 1 1 2 1 2 1 1 1;   #  6
    1 1 1 1 1 1 1 1 1 1;   #  7
    1 1 1 1 1 1 1 1 1 1;   #  8
    1 1 1 1 1 1 1 1 1 1;   #  9
    1 2 3 3 3 3 3 1 1 1;   # 10
    1 1 1 1 1 1 1 1 1 1;   # 11
    1 1 1 2 3 2 3 1 1 1;   # 12
    1 1 1 1 1 1 1 1 1 1;   # 13
    1 1 1 1 1 1 1 1 1 1;   # 14
    1 1 1 1 1 1 1 1 1 1;   # 15
    1 1 1 1 1 1 1 1 1 1;   # 16
    1 1 1 1 1 1 1 1 1 1;   # 17
    1 1 1 1 1 1 1 1 1 1;   # 18
    1 1 1 1 1 1 1 1 1 1;   # 19
    1 1 1 1 1 1 1 1 1 1;   # 20
    1 1 1 1 1 1 1 1 1 1;   # 21
    1 1 1 1 1 1 1 1 1 1;   # 22
    1 1 1 1 1 1 1 1 1 1;   # 23
    1 1 1 1 1 1 1 1 1 1;   # 24
    1 1 1 1 1 1 1 1 1 1;   # 25
    1 1 1 1 1 1 1 1 1 1;   # 26
    1 1 1 1 1 1 1 1 1 1;   # 27
    1 1 1 1 1 1 1 1 1 1;   # 28
    1 1 1 1 1 1 1 1 1 1;   # 29
    1 1 1 1 1 1 1 1 1 1;   # 30
    1 1 1 1 1 1 1 1 1 1;   # 31
    1 1 1 1 1 1 1 1 1 1;   # 32
    1 1 1 1 1 1 1 1 1 1]   # 33

# DGDS[isp, idsq] (so/dgf.f DATA DGDS(4,33) → jl [isp,idsq]; column selected by MAPDSQ).
const SO_DGDS = Float32[
   -0.000091   0.0        0.0        0.0;   #  1
   -0.000295   0.0        0.0        0.0;   #  2
   -0.000066   0.0        0.0        0.0;   #  3
   -0.000414  -0.000239  -0.000325   0.0;   #  4
   -0.000252   0.0        0.0        0.0;   #  5
   -0.000251  -0.000133   0.0        0.0;   #  6
   -0.000099   0.0        0.0        0.0;   #  7
    0.0        0.0        0.0        0.0;   #  8
   -0.0004572  0.0        0.0        0.0;   #  9
   -0.000324  -0.000356  -0.000120   0.0;   # 10
    0.0        0.0        0.0        0.0;   # 11
   -0.000414  -0.000239  -0.000325   0.0;   # 12
   -0.000468   0.0        0.0        0.0;   # 13
   -0.0002189  0.0        0.0        0.0;   # 14
   -0.0003996  0.0        0.0        0.0;   # 15
   -0.0006538  0.0        0.0        0.0;   # 16
   -0.0001683  0.0        0.0        0.0;   # 17
    0.0        0.0        0.0        0.0;   # 18
   -0.0001546  0.0        0.0        0.0;   # 19
   -0.0001323  0.0        0.0        0.0;   # 20
    0.0        0.0        0.0        0.0;   # 21
    0.0        0.0        0.0        0.0;   # 22
   -0.0001737  0.0        0.0        0.0;   # 23
    0.0        0.0        0.0        0.0;   # 24
    0.0        0.0        0.0        0.0;   # 25
    0.0        0.0        0.0        0.0;   # 26
   -0.0003048  0.0        0.0        0.0;   # 27
    0.0        0.0        0.0        0.0;   # 28
    0.0        0.0        0.0        0.0;   # 29
    0.0        0.0        0.0        0.0;   # 30
    0.0        0.0        0.0        0.0;   # 31
   -0.000066   0.0        0.0        0.0;   # 32
    0.0        0.0        0.0        0.0]   # 33

# IBSERV[isic, isp] (so/dgf.f DATA IBSERV(5,33)) — ISIC-classed calibration count for sp 11/16/24.
# All species 5*1000 except WB(16)=[27,70,123,101,33] and AS(24)=[184,429,356,162,74].
const SO_IBSERV_16 = Float32[27,70,123,101,33]
const SO_IBSERV_24 = Float32[184,429,356,162,74]

# so/dgf.f species index sets.
const SO_ASPTEM_SP = (11, 16, 24)                                  # ASPECT − 45° (0.7854 rad)
const SO_ELCAP_SP  = Set{Int}([21,25,26,28,29,30,31,33])           # TEMEL capped at 30
const SO_SL0_SP    = (14, 17)                                      # slope==0 → SASP = SLODUM
const SO_ISIC_SP   = (11, 16, 24)                                  # use IBSERV(ISIC) not OBSERV

# so/maical.f DATA ISPNUM(33): SO species → the FIA-code proxy fed to ADJMAI.
const SO_MAI_ISPNUM = Int32[119,117,202,15,264,122,108,93,101,122,101,15,21,101,101,101,101,101,101,101,746,746,746,746,746,746,746,746,746,746,746,202,746]

"""index site species (ISISP). so/sitset.f default = PP(10); shared across DGCONS/MAICAL."""
@inline so_isisp(p) = (1 <= Int(p.site_species) <= 33) ? Int(p.site_species) : 10

"so/maical.f RMAI = min(ADJMAI(ISPNUM(ISISP), converted SITEAR(ISISP), 10), 128). Stand constant."
function _so_rmai(s::StandState)::Float32
    p = s.plot
    isisp = Int(p.site_species); isisp == 0 && (isisp = 3)         # maical.f: ISISP==0 ⇒ DF(3)
    sssi = p.sp_site_index[isisp]; sssi == 0f0 && (sssi = 140f0)
    isisp == 3 && (sssi = -43.78f0 + 2.16f0 * sssi)               # DF: Cochran → McArdle
    isisp == 5 && (sssi = sssi * 3.28f0)                          # MH: metric → English
    rmai = _adjmai(SO_MAI_ISPNUM[isisp], sssi, 10f0)
    rmai > 128f0 && (rmai = 128f0)
    return rmai
end

# so/dgf.f ENTRY DGCONS — per-species DGCON + ATTEN. IFOR used directly (no JFOR remap). Also caches RMAI.
function so_dgcons!(s::StandState)
    c = s.calib; p = s.plot
    ifor = Int(p.forest_idx)
    elev = p.elevation; slope = p.slope; asp = p.aspect
    sina = sin(asp); cosa = cos(asp)
    rmai = _so_rmai(s); p.mai_adj = rmai                          # so/maical.f RMAI (used by DGCON & dgf!)
    isisp = so_isisp(p)
    xsite_isisp = p.sp_site_index[isisp]                          # SITEAR(ISISP), used only in DGSIC·XSITE
    # ISIC: SITEAR(ISISP)/10 → 5 classes (so/dgf.f:790-796), for the IBSERV-classed species 11/16/24.
    lsi = Int(floor(xsite_isisp / 10f0))
    isic = lsi < 2 ? 1 : lsi >= 5 ? 5 : lsi
    @inbounds for isp in 1:33
        sitear = p.sp_site_index[isp]
        isfor = SO_MAPLOC[isp, ifor]
        asptem = isp in SO_ASPTEM_SP ? asp - 0.7854f0 : asp
        temel = elev
        (isp in SO_ELCAP_SP && temel > 30f0) && (temel = 30f0)
        if (isp in SO_SL0_SP) && slope <= 0f0
            sasp = SO_SLODUM[isp]
        else
            sasp = (SO_DGSASP[isp] * sin(asptem) + SO_DGCASP[isp] * cos(asptem) + SO_DGSLOP[isp]) * slope +
                   SO_DGSLSQ[isp] * slope * slope
        end
        dgcon = SO_DGFOR[isp, isfor] + SO_DGEL[isp] * temel + SO_DGEL2[isp] * temel * temel +
                SO_DGMAI[isp] * rmai + SO_DGSIC[isp] * xsite_isisp +
                SO_DGSITE[isp] * log(sitear) + sasp
        isp == 8 && (dgcon += 0.86756f0 * log(sitear))            # ES extra ln(SITEAR)
        c.dg_const[isp] = dgcon
        c.atten[isp] = isp == 11 ? 1000f0 :
                       isp == 16 ? SO_IBSERV_16[isic] :
                       isp == 24 ? SO_IBSERV_24[isic] : SO_OBSERV[isp]
    end
    return s
end

# so/dgf.f main body — per-tree WK2 = LN(DDS). 4 branches: JU(11), AS(24), RA(22), GENERAL.
function dgf!(s::StandState, ::SouthCentralOregon)
    p, t, c = s.plot, s.trees, s.calib
    dens = s.density; sd = s.coef.species
    wk2 = view(s.scratch.wk, 2, :)
    ba = p.basal_area; avh = p.avg_height
    rmai = p.mai_adj; relden = p.relative_density
    alccf = relden > 0f0 ? log(relden) : 0f0
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        isp = Int(t.species[i])
        pt_i = Int(t.plot_id[i])
        pccf = (1 <= pt_i <= length(dens.point_ccf)) ? dens.point_ccf[pt_i] : 0f0
        cr = Float32(t.crown_pct[i]) * 0.01f0
        pctfrac = 1f0 - t.crown_ratio[i] / 100f0
        bal = pctfrac * ba
        cor = c.dg_cor[isp]
        si = p.sp_site_index[isp]
        bark = so_bratio(sd, isp, d)
        relht = avh > 0f0 ? min(t.height[i] / avh, 1.5f0) : 0f0
        # CONSPP = DGCON + COR + 0.01·DGCCFA·RELDEN + 0.01·DGMACC·RMAI·RELDEN
        conspp = c.dg_const[isp] + cor + 0.01f0 * SO_DGCCFA[isp] * relden +
                 0.01f0 * SO_DGMACC[isp] * rmai * relden
        if isp == 11                                              # JU: western-juniper DF-projection (UT)
            dpp = d < 1f0 ? 1f0 : d
            batem = ba < 1f0 ? 1f0 : ba
            df = 0.25897f0 + 1.03129f0 * dpp - 0.0002025464f0 * batem + 0.00177f0 * si
            (df - dpp) > 1f0 && (df = dpp + 1f0)
            df < dpp && (df = dpp)
            diagr = (df - dpp) * bark
            dds = diagr <= 0f0 ? -9.21f0 : log(diagr * (2f0 * dpp * bark + diagr)) + conspp
        elseif isp == 24                                          # AS: aspen DGFASP (UT), raw crown pct
            cr_raw = Float32(t.crown_pct[i])
            rmsqd = _TT_CUR_RMSQD[] >= 0f0 ? _TT_CUR_RMSQD[] : stand_qmd(s)
            aspdg = _em_dgfasp(d, cr_raw, bark, si, rmsqd, ba)
            dds = aspdg + log(cor2_of(c, isp)) + cor
        elseif isp == 22                                          # RA: red alder decreasing-increment eqn
            const0 = 3.250531f0 - 0.003029f0 * ba
            diagr = d <= 18f0 ? const0 - 0.166496f0 * d + 0.004618f0 * d * d :
                                const0 - (const0 / 10f0) * (d - 18f0)
            diagr < 0.1f0 && (diagr = 0.1f0)
            dds = log(diagr * (2f0 * d * bark + diagr)) + log(cor2_of(c, isp)) + cor
        else                                                      # GENERAL Wykoff ln(DDS)
            dummy = isp == 14 ? -0.799079f0 : 0f0                 # SF(14) DUMMY from EC
            balx = (isp == 11 || isp == 16) ? bal / 100f0 : bal   # sp 11/16 scaled BAL
            dglccf = 0f0; dgpcf2 = 0f0
            if isp == 10                                          # PP: extra ln(RELDEN) + PCCF² terms
                dglccf = -0.20488f0
                dgpcf2 = pccf > 400f0 ? 0f0 : 0.000001f0
            end
            dgdsq = SO_DGDS[isp, SO_MAPDSQ[isp, Int(p.forest_idx)]]
            dds = conspp + SO_DGLD[isp] * log(d) + SO_DGBAL[isp] * balx +
                  cr * (SO_DGCR[isp] + cr * SO_DGCRSQ[isp]) +
                  dgdsq * d * d + SO_DGDBAL[isp] * balx / log(d + 1f0) +
                  SO_DGPCCF[isp] * pccf + dglccf * alccf + dgpcf2 * pccf * pccf +
                  SO_DGHAH[isp] * relht + SO_DGLBA[isp] * log(ba) + SO_DGBA[isp] * ba + dummy
            isp == 8 && (dds += 0.49649f0 * relht)                # ES HOAVH bonus
        end
        dds < -9.21f0 && (dds = -9.21f0)
        wk2[i] = dds
    end
    return s
end
