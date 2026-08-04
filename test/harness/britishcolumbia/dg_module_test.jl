const BC_DATADIR = "/workspace/FVSjl/data/britishcolumbia"
include("/workspace/FVSjl/src/variants/britishcolumbia/dg_coefficients.jl")

function test()
    println("Loaded: ZNKONST=$(length(BC_ZNKONST)) SSKONST=$(length(BC_SSKONST))")
    # 1) matcher vs ground truth
    truth = Dict(1=>(1,1),2=>(2,2),3=>(5,3),4=>(17,21),5=>(10,5),6=>(11,7),7=>(13,11),
                 8=>(15,15),9=>(17,21),10=>(18,24),11=>(19,25),12=>(21,25),13=>(21,25),14=>(5,3),15=>(19,25))
    mok=0
    for sp in 1:15
        ip,jp = bc_resolve_ipjp(sp, "ICHmw2/01", "ICH")
        (ip,jp)==truth[sp] && (mok+=1)
    end
    println("matcher: $mok/15 ", mok==15 ? "✓" : "✗")
    # 2) DDS vs 2106 oracle points (sp14 -> ip=5, zone ICH)
    rows=[split(l,',') for l in split(strip(read("/workspace/FVSjl/test/harness/britishcolumbia/ref_dds_all_BC_sp14.csv",String)),'\n')]
    n=0; exact=0; maxe=0f0
    for r in rows
        d_cm=parse(Float32,r[2]); bal=parse(Float32,r[3]); cr=parse(Float32,r[4])
        conspp=parse(Float32,r[5]); brat=parse(Float32,r[6]); ddsO=parse(Float32,r[7])
        # NOTE ref D is already in cm; bc_v3_dds expects inches then *INtoCM. Pass d_in = d_cm/INtoCM
        my = bc_v3_dds(14, 5, "ICH", d_cm/BC_INtoCM, bal, cr, conspp, brat)
        n+=1; my==ddsO && (exact+=1); e=abs(my-ddsO); e>maxe && (maxe=e)
    end
    println("DDS: $exact/$n bit-exact ($(round(100*exact/n,digits=1))%) maxerr=$maxe")
end
test()
