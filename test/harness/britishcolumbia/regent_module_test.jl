const BC_DATADIR = "/workspace/FVSjl/data/britishcolumbia"
struct BritishColumbia end
nspecies(::BritishColumbia) = 15
include("/workspace/FVSjl/src/variants/britishcolumbia/dg_coefficients.jl")
include("/workspace/FVSjl/src/variants/britishcolumbia/regent_coefficients.jl")
function main()
    println("BC_STCOEF=", length(BC_STCOEF))
    ip = bc_resolve_stcoef(14, "ICHmw2/01", "ICH"); println("sp14 ip=$ip (expect 13)")
    zrh=[split(l,',') for l in split(strip(read("/workspace/FVSjl/test/harness/britishcolumbia/ref_regent_htg_all_BC_sp14.csv",String)),'\n')]
    n=0;ex=0;mx=0f0; asp=deg2rad(315f0); slope=0.30f0
    for r in zrh
        h1=parse(Float32,r[2]);bal=parse(Float32,r[3]);rdj=parse(Float32,r[4]);rhcon=parse(Float32,r[5]);vO=parse(Float32,r[6])
        my=bc_v3_sthg(14, ip, h1, bal, rdj, rhcon, asp, slope)
        n+=1; my==vO&&(ex+=1); e=abs(my-vO); e>mx&&(mx=e)
    end
    println("module bc_v3_sthg: $ex/$n bit-exact ($(round(100*ex/n,digits=1))%) maxerr=$mx")
end
main()
