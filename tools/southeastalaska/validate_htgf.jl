using FVSjl
const M = FVSjl
function run(path)
    B1=M.AK_HG_B1; B2=M.AK_HG_B2; B3=M.AK_HG_B3; B4=M.AK_HG_B4; B5=M.AK_HG_B5; B6=M.AK_HG_B6
    HTLO=M.AK_HTLO; HTHI=M.AK_HTHI
    spmult(sp)= (sp==14||sp==21||sp==23) ? 0.45f0 : (sp==22 ? 0.65f0 : 1.0f0)
    lines = split(read(path,String), '\n')
    worst_pot=0.0; worst_htg=0.0; n=0
    for i in 1:length(lines)
        l = lines[i]
        occursin("IN HTGF: I=", l) || continue
        g(key) = (m = match(Regex(key*"=\\s*([-0-9.]+)"), l); m===nothing ? NaN : parse(Float64, m.captures[1]))
        sp=Int(g("ISPC")); d=g("DBH"); h=g("HT"); temel=g("TEMEL"); xsite=g("XSITE"); dg10=g("DG10"); pfh=g("PFHMOD"); poto=g("POTHTG")
        htgo = NaN
        if i+1 <= length(lines) && occursin("HTGF 120F", lines[i+1])
            nums = [parse(Float64,m.match) for m in eachmatch(r"[-0-9]+\.[0-9]+", lines[i+1])]
            htgo = nums[2]
        end
        D=Float32(d)
        basehg = exp(B1[sp]+B2[sp]*D*D+B3[sp]*log(D)+B4[sp]*Float32(temel)+B5[sp]*log(Float32(xsite))+B6[sp]*log(Float32(dg10)))
        pot = 10f0*basehg*Float32(pfh)*spmult(sp)
        H=Float32(h)
        hgbnd = (H>=HTLO[sp] && H<HTHI[sp]) ? max(1f0-(H-HTLO[sp])/(HTHI[sp]-HTLO[sp]),0.1f0) : (H<HTLO[sp] ? 1f0 : 0.1f0)
        htg = pot*hgbnd; htg<=0.1f0 && (htg=0.1f0)
        dpot = abs(pot-poto)/max(abs(poto),1e-6)
        dhtg = isnan(htgo) ? 0.0 : abs(htg-htgo)/max(abs(htgo),1e-6)
        worst_pot=max(worst_pot,dpot); worst_htg=max(worst_htg,dhtg); n+=1
        (dpot>1e-4||dhtg>1e-4) && println("MISMATCH sp=$sp D=$d H=$h POT jl=$pot fvs=$poto HTG jl=$htg fvs=$htgo")
    end
    println("checked $n trees; worst POTHTG rel-err=$worst_pot ; worst HTG rel-err=$worst_htg")
end
run(ARGS[1])
