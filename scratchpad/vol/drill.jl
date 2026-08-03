using FVSjl
const M = FVSjl
keyf = ARGS[1]
var = ARGS[2]
V = var=="CI" ? M.CentralIdaho() : var=="UT" ? M.Utah() : var=="EM" ? M.EasternMontana() : error(var)
for s in M.each_stand(keyf; variant=V, faithful=true)
    M.notre!(s); M.setup_growth!(s); M.compute_volumes!(s)
    t = s.trees
    veq = s.species.vol_eq
    println("STAND ", strip(s.plot.stand_id), "  n=", t.n, " ndead=", t.ndead)
    println(rpad("i",4), rpad("sp",4), rpad("eq",12), rpad("dbh",8), rpad("ht",7), rpad("tpa",9),
            rpad("tcf",9), rpad("mcf",9), rpad("bdf",9))
    tot_tcf=0.0; tot_mcf=0.0; tot_bdf=0.0
    for i in 1:(t.n)
        sp = Int(t.species[i])
        tpa = t.tpa[i]
        println(rpad(i,4), rpad(sp,4), rpad(strip(veq[sp]),12),
                rpad(round(t.dbh[i],digits=2),8), rpad(round(t.height[i],digits=1),7),
                rpad(round(tpa,digits=3),9),
                rpad(round(t.cuft_vol[i],digits=3),9), rpad(round(t.merch_cuft_vol[i],digits=3),9),
                rpad(round(t.bdft_vol[i],digits=3),9))
        tot_tcf += t.cuft_vol[i]*tpa; tot_mcf += t.merch_cuft_vol[i]*tpa; tot_bdf += t.bdft_vol[i]*tpa
    end
    println("PERACRE TCuFt=", round(tot_tcf,digits=1), " MCuFt=", round(tot_mcf,digits=1), " BdFt=", round(tot_bdf,digits=1))
    println()
end
