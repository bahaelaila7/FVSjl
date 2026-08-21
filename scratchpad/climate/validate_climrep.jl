using FVSjl
const F = FVSjl
V = F.variant_from_code("IE")
for s in F.each_stand("clim_iet.key"; variant=V)
    F.notre!(s); F.setup_growth!(s); F.compute_volumes!(s)
    (s.climate === nothing || !s.climate.active) && (println("clim inactive"); break)
    ry = Int(F.current_cycle_year(s))          # 1990 (cycle-0 start = report label)
    F.grow_cycle!(s; fint=10f0)                 # advance one 10-yr cycle (post-growth tree list)
    rep = F.climate_report(s; report_year=ry, fint=10)
    println("=== jl climate_report report_year=", ry, " (n=", length(rep), ") ===")
    for r in rep
        sym = s.climate.plant_symbols[r.sp]
        println("  ", rpad(sym,6), " viab=", rpad(round(r.viab,digits=4),7),
                " BA=", rpad(round(r.ba,digits=2),7), " TPA=", rpad(round(r.tpa,digits=2),8),
                " gm=", rpad(round(r.gmult,digits=4),7), " sm=", rpad(round(r.sitgm,digits=4),6),
                " mxd=", rpad(round(r.mxden,digits=3),6), " M1=", rpad(round(r.mort1,digits=4),7),
                " ae=", round(r.potestab,digits=2))
    end
    break
end
