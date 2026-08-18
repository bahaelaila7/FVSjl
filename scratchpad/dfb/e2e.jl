using FVSjl; const F=FVSjl
v=F.InlandEmpire()
D="/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/mrun/"
# ensure .tre next to each key
for kb in ("ransched","ranstart","curout_pk","curout_lv","windthr")
    cp(D*"iet01.tre", D*kb*".tre"; force=true)
end
rows(kb)=filter(l->!startswith(l,"-999"), split(strip(F.run_keyfile(D*kb*".key"; variant=v, output=:sum)),'\n'))
# control = a DFB-off run: reuse ransched.key but no outbreak? Build a plain control from ransched with DFB removed is easier: use windthr as base minus. Instead compare cycles.
mort(row)=parse(Int, split(row)[25])
for kb in ("ransched","ranstart","curout_pk","curout_lv","windthr")
    r=rows(kb)
    ms=[mort(x) for x in r]
    println(rpad(kb,10)," MOR/cyc=",ms)
end
