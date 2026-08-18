# Validate the STAGED coefficient file + core Penner function reproduce the g16 dump bit-exact.
using Printf
include("dg_coefficients.jl")
@inline on_logf(x::Float32) = ccall(:logf, Float32, (Float32,), x)
@inline on_expf(x::Float32) = ccall(:expf, Float32, (Float32,), x)
@inline function on_penner_dds(ksp, ags, diam_in::Float32, sim, bam, qmdm, balm, htm, bark)
    dbhm = diam_in * ON_INtoCM
    x1 = -(ON_B0[ksp]) - (ON_BBAL[ksp]*balm) - (ON_BHT[ksp]*htm) - (ON_BSI[ksp]*sim) -
         (ON_BBA[ksp]*bam) - (ON_BDBHQ[ksp]*qmdm) - (ON_BAGS[ksp]*Float32(ags))
    for _ in 1:10
        dgln = x1 + (ON_B1[ksp]*on_logf(dbhm)) - (ON_B2[ksp]*dbhm)
        deld = min(max(dgln, -5f0), 5f0); deld = on_expf(deld)
        deld = min(max(deld, 0.0001f0), ON_B95[ksp]); dbhm = dbhm + deld
    end
    d = dbhm * ON_CMtoIN; diagr = (d - diam_in) * bark
    return dbhm, diagr, diagr*(2f0*diam_in*bark + diagr)
end
h2f(h) = reinterpret(Float32, parse(UInt32, h; base=16))
f2h(x::Float32) = uppercase(string(reinterpret(UInt32, x); base=16, pad=8))
n=0; ok=0; seen=Set{Tuple{Int,Int}}()
for ln in eachline(get(ARGS,1,"/workspace/.onwork/run/fort.771"))
    f=split(strip(ln)); isempty(f) && continue
    I,ISPC,KSP,AGS,ICYC=parse.(Int,f[1:5])
    diam,dbo,sim,bam,qmdm,balm,htm,bark,dgo,ddo=h2f.(f[6:15])
    (I,ISPC) in seen && continue; push!(seen,(I,ISPC)); global n+=1
    db,dg,dd=on_penner_dds(KSP,AGS,diam,sim,bam,qmdm,balm,htm,bark)
    m = f2h(db)==f2h(dbo) && f2h(dg)==f2h(dgo) && f2h(dd)==f2h(ddo)
    global ok+=m
    @printf("I=%d ISPC=%d KSP=%d AGS=%d  DDS or=%s jl=%s  %s\n",I,ISPC,KSP,AGS,f2h(ddo),f2h(dd), m ? "OK" : "MISMATCH")
end
@printf("\nSTAGED on_penner_dds: %d/%d bit-exact (DBHM+DIAGR+DDS)\n", ok, n)
