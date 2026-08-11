using FVSjl
key = "scratchpad/nct01_ctl.key"
out = FVSjl.run_keyfile(key; variant = FVSjl.Klamath(), period = 5)
print(String(take!(copy(IOBuffer(out)))) === "" ? out : out)
