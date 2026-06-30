using FVSjl
strip_ts(s)=replace(s, r"\d\d-\d\d-\d{4} \d\d:\d\d:\d\d"=>"TS")
for k in sort(filter(f->endswith(f,".key"), readdir("/tmp/sweep/scn"; join=true)))
    n=basename(k)[1:end-4]
    V = startswith(n,"ne") ? FVSjl.Northeast() : FVSjl.Southern()
    y="/tmp/sweep/scn/$n.yaml"
    try
        FVSjl.translate_io(k, y)
        sk=run_keyfile(k; variant=V); sy=run_keyfile(y; variant=V)
        write("/tmp/sweep/scn/$n.jlkey.sum", sk)
        rt = strip_ts(sk)==strip_ts(sy) ? "OK" : "FAIL"
        println("$n yaml_rt=$rt")
    catch e
        println("$n ERROR: ", sprint(showerror,e)[1:min(80,end)])
    end
end
