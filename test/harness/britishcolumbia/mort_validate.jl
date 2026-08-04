# BC chunk-7 mortality validation: replay oracle-dumped mortality inputs (BCMSTAND/BCMTREE from the
# instrumented FVSbc on all_BC) through mortality!'s equation and compare RIP/RIPP/WKI bit-exactly.
using FVSjl
const OUT = "/workspace/.bcwork/mrun/all_BC.out"
f32(x) = parse(Float32, x)

function main()
    lines = readlines(OUT)
    stand = Dict{Int,NamedTuple}()
    for l in lines
        startswith(l, "BCMSTAND") || continue
        p = split(l)
        stand[parse(Int, p[2])] = (aved=f32(p[3]), bax=f32(p[4]), bacls=parse(Int, p[5]),
            bamax=f32(p[6]), rz=f32(p[7]), sdimax=f32(p[8]), dq10=f32(p[9]), ba=f32(p[10]), fint=f32(p[11]))
    end
    nR=0; exR=0; exP=0; exW=0; mxR=0f0; mxP=0f0; mxW=0f0
    for l in lines
        startswith(l, "BCMTREE") || continue
        p = split(l)
        cyc = parse(Int, p[2]); sp = parse(Int, p[3]); d = f32(p[4])
        ripO = f32(p[5]); rippO = f32(p[6]); wkiO = f32(p[7]); pr = f32(p[8])
        st = stand[cyc]
        dx = max(2.5f0, d * FVSjl.BC_INtoCM)
        mrtcls = FVSjl.BC_MORT_MRTCLS[sp]
        dbhcls = 4
        for j in 1:3
            if dx <= FVSjl.BC_DBHBRK[j]; dbhcls = j; break; end
        end
        rip = FVSjl.BC_FMRT[1, mrtcls, st.bacls, dbhcls] / 100f0
        ripp = st.ba * st.rz
        st.ba <= st.bamax && (ripp += (st.bamax - st.ba) * rip)
        ripp /= st.bamax
        ripp < rip && (ripp = rip); ripp > 1f0 && (ripp = 1f0)
        wki = pr * (1f0 - (1f0 - ripp)^st.fint)
        nR += 1
        rip  == ripO  && (exR += 1); mxR = max(mxR, abs(rip - ripO))
        ripp == rippO && (exP += 1); mxP = max(mxP, abs(ripp - rippO))
        wki  == wkiO  && (exW += 1); mxW = max(mxW, abs(wki - wkiO))
    end
    println("BC mortality replay (all_BC, all cycles, V3 FD-tabular):")
    println("  RIP : $exR/$nR bit-exact ($(round(100*exR/nR,digits=1))%)  maxerr=$mxR")
    println("  RIPP: $exP/$nR bit-exact ($(round(100*exP/nR,digits=1))%)  maxerr=$mxP")
    println("  WKI : $exW/$nR bit-exact ($(round(100*exW/nR,digits=1))%)  maxerr=$mxW")
end
main()
