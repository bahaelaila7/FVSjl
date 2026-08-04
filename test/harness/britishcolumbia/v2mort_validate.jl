# BC V2 mortality RIP-formula validation. Replays the instrumented FVSbc MORV2 dump (all_BC_essf cycle1,
# sp14) through the ported Hamilton RIP (morts.f:572-578) using the oracle's per-tree G/RELDBH/BA, and
# checks RIP bit-exact. Validates the RIP equation + REIN(IP) + B5=PMSC. (G-computation + shared RIPP/WKI
# tail validated separately / already for V3.) Constants GMULT/REIN confirmed via MORCON instrument.
using FVSjl; const F=FVSjl
const OUT = "/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/bcv2/morv2.txt"
f32(x)=parse(Float32,x)
poten1=F.BC_POT[F.BC_MORT_IPDG_44]; poten2=F.BC_POT[F.BC_MORT_IPDG2_44]
rein1=(1f0-(poten1/20f0+1f0)^(-1.605f0))/0.06821f0
rein2=(1f0-(poten2+1f0)^(-1.605f0))/0.86610f0

function rip_formula(sp,d,g,reldbh,ba)
    dd = d<=0.5f0 ? 0.5f0 : d
    ip = dd<=5f0 ? 2 : 1
    rip = 2.76253f0 + 0.222310f0*sqrt(dd) - 0.0460508f0*sqrt(ba) + 11.2007f0*g -
          0.554421f0/dd + F.BC_PMSC[sp] + 0.246301f0*reldbh + 6.07129f0*g/dd
    rip = clamp(rip,-70f0,70f0)
    rip = 1f0/(1f0+exp(rip))
    rip * (ip==1 ? rein1 : rein2)
end

function main()
    rows=[split(l) for l in readlines(OUT) if startswith(l,"MORV2")]
    n=0; ex=0; mx=0f0
    for r in rows
        sp=parse(Int,r[3]); d=f32(r[4]); g=f32(r[5]); reldbh=f32(r[6]); ba=f32(r[7]); ripO=f32(r[8])
        my=rip_formula(sp,d,g,reldbh,ba)
        n+=1; my==ripO && (ex+=1); mx=max(mx,abs(my-ripO))
    end
    println("BC V2 mortality RIP replay (all_BC_essf, sp14): $ex/$n bit-exact ($(round(100ex/n,digits=1))%)  maxerr=$mx")
end
main()
