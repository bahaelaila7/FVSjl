using FVSjl
const KEY = "test/fixtures/canonical/snt01.key"

# warm up (compile) then measure full-run allocation, faithful path
function run_once()
    FVSjl.run_keyfile(KEY; variant=FVSjl.Southern(), faithful=true)
end
run_once()                       # compile
a1 = @allocated run_once()
run_once()
a2 = @allocated run_once()
println("run_keyfile(snt01; faithful) @allocated = ", a1, " / ", a2, " B/run")
