using FVSjl; const F=FVSjl
v=F.PacificNorthwest()
D="/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/pnrun/"
rows(kb)=filter(l->!startswith(l,"-999"), split(strip(F.run_keyfile(D*kb*".key"; variant=v, output=:sum)),'\n'))
mort(row)=parse(Float64, split(row)[25]); yr(row)=split(row)[1]
c=rows("pnctrl"); d=rows("pndfb")
println("PN FVSjl:  year  ctrlMOR  dfbMOR   (oracle cyc2: 14 -> 36)")
for i in 1:length(c)
    println("  ",yr(c[i]),"  ",lpad(mort(c[i]),6),"  ",lpad(mort(d[i]),6))
end
println("cyc1 identical (inert): ", c[1]==d[1])
println("cyc2 DFB > ctrl: ", mort(d[2]) > mort(c[2]))
