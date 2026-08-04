const BC_DATADIR = "/workspace/FVSjl/data/britishcolumbia"
struct BritishColumbia end
nspecies(::BritishColumbia) = 15
include("/workspace/FVSjl/src/variants/britishcolumbia/dg_coefficients.jl")
include("/workspace/FVSjl/src/variants/britishcolumbia/height_growth_coefficients.jl")
function main()
    ip14 = bc_resolve_lthg(14, "ICHmw2/01", "ICH")
    println("sp14 LTHG ip = $ip14  (SI=$(ip14>0 ? BC_LTHG[ip14].SI : 0))")
    rows=[split(l,',') for l in split(strip(read("/workspace/FVSjl/test/harness/britishcolumbia/ref_htg_all_BC_sp14.csv",String)),'\n')]
    n=0; exact=0; maxe=0f0; worst=""
    for r in rows
        isp=parse(Int,r[1]); ht=parse(Float32,r[2]); dbh=parse(Float32,r[3]); dg=parse(Float32,r[4])
        vO=parse(Float32,r[5]); mlt=parse(Float32,r[6])
        ip = bc_resolve_lthg(isp, "ICHmw2/01", "ICH")
        my = bc_v3_htg(isp, ip, ht, dbh, dg; mlt=mlt)
        n+=1; my==vO && (exact+=1); e=abs(my-vO); if e>maxe; maxe=e; worst="ht=$ht dbh=$dbh dg=$dg my=$my or=$vO"; end
    end
    println("n=$n bit-exact=$exact ($(round(100*exact/n,digits=1))%) maxabserr=$maxe")
    println("worst: $worst")
end
main()
