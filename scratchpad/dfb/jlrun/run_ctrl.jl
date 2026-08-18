using FVSjl
v = FVSjl.InlandEmpire()
key = "/tmp/claude-1000/-workspace/b4e1b3b1-495b-403e-810b-5db3604b56cc/scratchpad/dfb/jlrun/ctrl.key"
s = FVSjl.run_keyfile(key; variant=v, output=:sum)
for l in split(strip(s), '\n')
    startswith(l, "-999") && continue
    println(l)
end
