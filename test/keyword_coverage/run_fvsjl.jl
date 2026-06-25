# run_fvsjl.jl — batch the FVSjl side of the keyword-coverage harness in ONE process.
# For every scenarios/*.key: build the modern .csv (from .tre) + .yaml, run the engine
# on BOTH the .key and the .yaml, write <name>.key.sum, and report whether yaml==key.
using FVSjl
const SNFMT = "(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,T52,I2,T66,5I1,T54,7I1,T75,F3.0)"
scen = ARGS[1]
keys = sort(filter(f -> endswith(f, ".key"), readdir(scen; join=true)))
for key in keys
    name = first(splitext(basename(key)))
    tre = joinpath(scen, name * ".tre")
    isfile(tre) && convert_tre_to_csv(tre, joinpath(scen, name * ".csv"); fmt = SNFMT)
    local ksum, ysum
    try
        ksum = run_keyfile(key)
        write(joinpath(scen, name * ".key.sum"), ksum)
    catch e
        println("$name\tKEYERR\t", sprint(showerror, e)[1:min(60,end)]); continue
    end
    yamlf = joinpath(scen, name * ".yaml")
    if isfile(yamlf)
        try
            ysum = run_keyfile(yamlf)
            println("$name\t", ksum == ysum ? "yPASS" : "yFAIL")
        catch e
            println("$name\tYAMLERR\t", sprint(showerror, e)[1:min(60,end)])
        end
    else
        println("$name\tNOYAML")
    end
end
