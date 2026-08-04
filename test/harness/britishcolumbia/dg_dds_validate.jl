# Validate the BC V3 DEFAULT DDS formula vs 2106 oracle data points (all sp14 -> ZNKONST(5)).
DIR="/workspace/FVSjl/data/britishcolumbia"
zn=[split(l,',') for l in split(strip(read(joinpath(DIR,"dg_znkonst.csv"),String)),'\n')]
hdr=zn[1]; col(n)=findfirst(==(n),hdr)
row5=first(r for r in zn[2:end] if r[col("idx")]=="5")
c(n)=parse(Float32, row5[col(n)])
LD,DSQ,BALc,DBAL1,DBAL2,CRc = c("LD"),c("DSQ"),c("BAL"),c("DBAL1"),c("DBAL2"),c("CR")
println("ZNK(5): LD=$LD DSQ=$DSQ BALc=$BALc DBAL1=$DBAL1 DBAL2=$DBAL2 CRc=$CRc")
CMtoIN = 0.3937f0

function dds_default(D,BAL,CR,CONSPP,BRAT)
    pre = CONSPP + LD*D + DSQ*D*D + BALc*BAL + DBAL1*(BAL/D) + DBAL2*(BAL/log(D+1f0)) + CRc*CR
    dcm = exp(max(-9.21f0, pre))
    dds = (dcm*dcm + 2f0*dcm*D*BRAT) * CMtoIN * CMtoIN
    return max(-9.21f0, log(max(0.001f0, dds)))
end

function main()
    rows=[split(l,',') for l in split(strip(read("/workspace/FVSjl/test/harness/britishcolumbia/ref_dds_all_BC_sp14.csv",String)),'\n')]
    n=0; exact=0; maxerr=0.0f0; worst=""
    for r in rows
        D=parse(Float32,r[2]); BAL=parse(Float32,r[3]); CR=parse(Float32,r[4])
        CONSPP=parse(Float32,r[5]); BRAT=parse(Float32,r[6]); ddsO=parse(Float32,r[7])
        my=dds_default(D,BAL,CR,CONSPP,BRAT)
        n+=1
        if my==ddsO; exact+=1; end
        e=abs(my-ddsO); if e>maxerr; maxerr=e; worst="D=$D BAL=$BAL my=$my or=$ddsO"; end
    end
    println("n=$n bit-exact=$exact ($(round(100*exact/n,digits=1))%)  maxabserr=$maxerr")
    println("worst: $worst")
end
main()
