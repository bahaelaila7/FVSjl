using FVSjl; const F=FVSjl
v=F.InlandEmpire()
D="/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/mrun/"
jl(kb)=filter(l->!startswith(l,"-999"), split(strip(F.run_keyfile(D*kb*".key"; variant=v, output=:sum)),'\n'))
orc(kb)=filter(l->!startswith(l,"-999"), split(strip(read(D*kb*"_clean.sum",String)),'\n'))
mort(row)=parse(Int, split(row)[25])
year(row)=split(row)[1]
for kb in ("ransched","windthr")
    j=jl(kb); o=orc(kb)
    println("=== $kb: year  jlMOR  orcMOR ===")
    for i in 1:min(length(j),length(o))
        println("  ",year(j[i]),"  ",lpad(mort(j[i]),4),"  ",lpad(mort(o[i]),4))
    end
end
