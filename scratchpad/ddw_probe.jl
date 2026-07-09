using FVSjl
const F = FVSjl

function probe(key)
    println("=== ", key, " ===")
    for s in F.each_stand(key; variant=F.Southern(), faithful=true)
        F.notre!(s); F.setup_growth!(s); F.compute_volumes!(s)
        carb = Tuple[]
        io = IOBuffer()
        F.write_sum_file(io, s; period=5, stand_id="S248112", mgmt_id="NONE",
                         variant=F.variant_code(s.variant), carbon_collect=carb)
        for row in carb
            yr = row[1]
            rep = row[2]                     # stand_carbon_report NamedTuple
            dw  = row[5]                      # ffe_down_wood NamedTuple
            ddw = rep.down_wood
            rd(x) = round(x; digits=3)
            vs = dw.vol_soft; vh = dw.vol_hard
            println(rpad(string(yr),6),
                    "  DDW(C)=", round(ddw; digits=4),
                    "  softTot=", rd(vs[end]), " hardTot=", rd(vh[end]))
            println("        soft bins=", map(rd, vs))
            println("        hard bins=", map(rd, vh))
        end
        break
    end
end

# introspect the report field name once
for s in F.each_stand("test/harness/scenarios/carbon_snt.key"; variant=F.Southern(), faithful=true)
    F.notre!(s); F.setup_growth!(s); F.compute_volumes!(s)
    r = F.stand_carbon_report(s)
    println("report fields = ", fieldnames(typeof(r)))
    dw = F.ffe_down_wood(s)
    println("down_wood fields = ", fieldnames(typeof(dw)))
    break
end

probe("test/harness/scenarios/carbon_snt.key")
probe("/tmp/fia_val/csp.key")
