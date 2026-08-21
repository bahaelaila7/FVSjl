using FVSjl
const F = FVSjl
V = F.variant_from_code("OC")
key = "ocsnag_jl.key"
io = IOBuffer()
for s in F.each_stand(key; variant=V)
    F.notre!(s); F.setup_growth!(s); F.compute_volumes!(s)
    sid = strip(s.plot.stand_id)
    carb = Tuple[]
    F.write_sum_file(io, s; period=5, stand_id=String(sid),
        mgmt_id="NONE", variant=F.variant_code(s.variant), date="", time="",
        collect_rows=F.SummaryRow[], carbon_collect=carb)
    println("=== jl snag_summary per cycle (hard slots 1,2,3,total ; soft total) ===")
    for row in carb
        yr = row[1]; ss = row[4]   # snag_summary = (hard=7tuple, soft=7tuple)
        h = ss.hard; sf = ss.soft
        println("  ", yr, "  h1=", round(h[1],digits=2), " h2=", round(h[2],digits=2),
                " h3=", round(h[3],digits=2), " hTot=", round(h[7],digits=2), " sTot=", round(sf[7],digits=2))
    end
    break
end
