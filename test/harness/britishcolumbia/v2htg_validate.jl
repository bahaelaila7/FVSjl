# BC V2 height-growth validation. Replays instrumented FVSbc HTV2 dump (all_BC_essf sp14 cycle1) through the
# ported exp-form HTG using oracle per-tree D/HTI/DG, checks HTG bit-exact-or-cornered. IHT=2 (ITYPE=4).
using FVSjl; const F=FVSjl
const OUT="/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/bcv2/htv2.txt"
f32(x)=parse(Float32,x)
iht=F.BC_HTG_MAPHAB[F.BC_V2_ITYPE]; hghch=F.BC_HGHC[iht]; h2cof=F.BC_HGH2[iht]; hdgcof=F.BC_HGLDD[iht]
function htg(sp,d,hti,dg)
    htcon=hghch+F.BC_HGSC[sp]
    con=htcon+h2cof*hti*hti+F.BC_HGLD[sp]*log(d)+F.BC_HTG_HGLH*log(hti)
    hexp=dg>0f0 ? exp(con+hdgcof*log(dg)) : 0f0
    (hexp+F.BC_HTG_BIAS)
end
function main()
    rows=[split(l) for l in readlines(OUT) if startswith(l,"HTV2")]
    n=0;ex=0;mx=0f0
    for r in rows
        sp=parse(Int,r[3]);d=f32(r[4]);hti=f32(r[5]);dg=f32(r[6]);htgO=f32(r[8])
        my=htg(sp,d,hti,dg); n+=1; my==htgO && (ex+=1); mx=max(mx,abs(my-htgO))
    end
    println("BC V2 HTG replay (all_BC_essf sp14): $ex/$n bit-exact ($(round(100ex/n,digits=1))%) maxerr=$mx")
end
main()
