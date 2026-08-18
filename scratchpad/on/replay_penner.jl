# Standalone bit-exact replay of Ontario Penner large-tree diameter growth (canada/on/dgf.f).
# Reads the g16 per-tree dump (fort.771: I ISPC KSP AGS ICYC + 10 Float32-hex:
#   DIAM(in) DBHM_final(cm) SIM BAM QMDM BALM HTM BARK DIAGR DDS), replays the 10-yr
# annual iteration in Float32, and checks DBHM_final / DIAGR / DDS bit-exact vs the oracle.
# Transcendentals: glibc single-precision logf/expf via ccall (measured 0-mismatch).
using Printf

const INtoCM = 2.54f0
const CMtoIN = 0.3937f0

# --- Penner coefficient arrays (35 eqns), canada/on/dgf.f DATA blocks ---
const B0 = Float32[3.691,1.8979,5.451,1.1499,4.184, 2.1184,1.3936,1.5596,2.8199,1.0315,
 1.4697,1.4792,2.3705,3.3822,2.256, 2.0375,1.5782,1.626,0.8673,1.5838,
 2.6442,1.4321,2.2574,1.0566,1.2195, 0.,0.,2.9716,3.2426,3.388,
 18.6647,1.6447,2.5493,1.1281,1.5392]
const B1 = Float32[1.0989,0.2965,1.5444,0.,1.0302, 0.6541,0.0648,0.1237,0.,0.,
 0.1872,0.0138,0.2704,0.6565,0.6342, 0.4857,0.0787,0.0892,0.196,0.3696,
 0.3003,0.,0.3295,0.286,0.2808, 0.5888,0.6024,0.8478,1.1018,1.1824,
 7.7658,0.0835,0.6732,0.,0.038]
const B2 = Float32[0.0396,0.0045,0.0766,0.0282,0.074, 0.1046,0.0158,0.0525,0.0147,0.0803,
 0.00433,0.000774,0.,0.00127,0.026, 0.0215,0.,0.,0.0189,0.0322,
 0.,0.00166,0.,0.0168,0.0148, 0.0692,0.0644,0.0234,0.0323,0.0372,
 0.297,0.00941,0.0186,0.00176,0.]
const BBAL = Float32[0.0477,0.0305,0.0205,0.0266,0.0426, 0.0552,0.0356,0.0661,0.0217,0.0468,
 0.0225,0.000676,0.,0.0111,0.0155, 0.0185,0.00224,0.00211,0.0302,0.0325,
 0.0209,0.0486,0.015,0.00708,0.00837, 0.0295,0.07,0.0131,0.00381,0.00655,
 0.0545,0.0327,0.0045,0.0151,0.0267]
const BHT = Float32[0.,0.,0.,0.,0., 0.,0.,0.,0.,0., 0.,0.,0.,0.,0., 0.,0.,0.,0.0132,0.,
 0.,0.,0.,0.,0., 0.,0.,0.,0.,0., 0.,0.,0.,0.,0.]
const BSI = Float32[-0.0117,0.,-0.0351,-0.0836,-0.0624, -0.0379,-0.018,-0.0596,-0.093,-0.0665,
 0.,0.00166,0.,0.,0., 0.,0.,0.,0.,0.,
 0.0101,-0.025,0.,0.,0., 0.,0.,0.,0.,0.,
 0.,0.,0.,0.,-0.0211]
const BBA = Float32[0.,0.0041,0.,0.0257,0., 0.,0.,0.,0.,0., 0.0226,0.,0.00875,0.0244,0.00705,
 0.00292,0.,0.,0.,0., 0.,0.,0.,0.0181,0.0145, 0.1151,0.1273,0.,0.0215,0.0205,
 0.,0.,0.022,0.0111,0.00405]
const BDBHQ = Float32[-0.015,-0.0089,-0.0343,0.,0., 0.,0.,0.,0.,0., -0.0201,0.,0.,0.,0.,
 0.,0.,0.,0.000716,-0.012, 0.,0.,0.,0.,0., 0.,0.,0.,0.,0., 0.,0.,0.,0.,0.]
const BAGS = Float32[0.,0.,0.,0.,0., 0.,0.,0.,0.,0., 0.,0.,0.,0.,0.,
 -0.2502,0.,-0.0622,0.,-0.1575, 0.,0.,0.0612,0.,-0.199,
 0.,-1.6033,0.,0.,-0.158, 0.,0.,0.,-0.5982,0.]
const B95 = Float32[0.56,0.51667,0.52,0.46228,0.24, 0.48,0.39,0.557,0.22,0.517,
 0.56,0.422,0.32,0.2,0.55, 0.55,0.52,0.52,0.514,0.514,
 0.241,0.46,0.46,0.575,0.575, 0.48,0.48,0.614,0.57,0.57,
 0.3,0.28,0.46,0.46,0.38]

logf(x::Float32) = ccall(:logf, Float32, (Float32,), x)
expf(x::Float32) = ccall(:expf, Float32, (Float32,), x)

h2f(h) = reinterpret(Float32, parse(UInt32, h; base=16))
f2h(x::Float32) = uppercase(string(reinterpret(UInt32, x); base=16, pad=8))

# Core Penner large-tree DG: returns DBHM_final (cm) and DDS, given dumped inputs.
function penner_dds(ksp, ags, diam_in::Float32, sim, bam, qmdm, balm, htm, bark)
    dbhm = diam_in * INtoCM
    x1 = -(B0[ksp]) - (BBAL[ksp]*balm) - (BHT[ksp]*htm) - (BSI[ksp]*sim) -
         (BBA[ksp]*bam) - (BDBHQ[ksp]*qmdm) - (BAGS[ksp]*Float32(ags))
    for _ in 1:10
        dgln = x1 + (B1[ksp]*logf(dbhm)) - (B2[ksp]*dbhm)
        deld = min(max(dgln, -5f0), 5f0)
        deld = expf(deld)
        deld = min(max(deld, 0.0001f0), B95[ksp])
        dbhm = dbhm + deld
    end
    d = dbhm * CMtoIN
    diagr = (d - diam_in) * bark
    dds = diagr * (2f0*diam_in*bark + diagr)
    return dbhm, diagr, dds
end

dumpfile = length(ARGS) >= 1 ? ARGS[1] : "/workspace/.onwork/run/fort.771"
nrec = 0; okdbhm = 0; okdiagr = 0; okdds = 0
seen = Set{Tuple{Int,Int}}()
println("  I ISPC KSP AGS   DIAM(in)   DBHM_or/jl(hex)      DIAGR_or/jl        DDS_or/jl        match")
for ln in eachline(dumpfile)
    f = split(strip(ln))
    isempty(f) && continue
    I,ISPC,KSP,AGS,ICYC = parse.(Int, f[1:5])
    diam,dbhm_o,sim,bam,qmdm,balm,htm,bark,diagr_o,dds_o = h2f.(f[6:15])
    (I,ISPC) in seen && continue; push!(seen,(I,ISPC))   # dedup repeated calls
    global nrec += 1
    dbhm_j, diagr_j, dds_j = penner_dds(KSP, AGS, diam, sim, bam, qmdm, balm, htm, bark)
    mB = f2h(dbhm_j)==f2h(dbhm_o); mG = f2h(diagr_j)==f2h(diagr_o); mD = f2h(dds_j)==f2h(dds_o)
    global okdbhm += mB; global okdiagr += mG; global okdds += mD
    @printf("%3d %4d %3d %3d  %9.4f  %s/%s %s/%s %s/%s  %s\n",
      I,ISPC,KSP,AGS, diam, f2h(dbhm_o),f2h(dbhm_j), f2h(diagr_o),f2h(diagr_j),
      f2h(dds_o),f2h(dds_j), (mB&&mG&&mD) ? "OK" : "MISMATCH")
end
println()
@printf("RESULT: DBHM %d/%d  DIAGR %d/%d  DDS %d/%d bit-exact\n",
        okdbhm,nrec, okdiagr,nrec, okdds,nrec)
