using FVSjl
const V = FVSjl.EasternMontana()
cd(@__DIR__)
for k in ("host_off_jl", "host_defol_jl")
    println("=== RUN $k ===")
    try
        txt = FVSjl.run_keyfile(k*".key"; variant=V)
        open(k*".jlsum2","w") do io; write(io, txt); end
        # print data rows (strip header line 1)
        for (i,ln) in enumerate(split(txt,'\n'))
            i == 1 && continue
            isempty(strip(ln)) && continue
            println(ln)
        end
    catch e
        println("ERROR in $k: ", e)
        Base.showerror(stdout, e, catch_backtrace()); println()
    end
end
