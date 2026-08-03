# Extract bm/htdbh.f forest-dependent Curtis-Arney P2/P3/P4 (MALHUR/OCHOCO/UMATIL/WALWIT, each [18,3])
# → data/bluemountains/htdbh_coeffs_bm.csv. IFOR 1-4 = Malheur/Ochoco/Umatilla/Wallowa-Whitman.
const SRC="/workspace/ForestVegetationSimulator/bm/htdbh.f"
const L=readlines(SRC)
function dv(name)
    i0=findfirst(l->occursin(Regex("^\\s+DATA\\s+"*name*"\\s*/"),l),L)
    body=String[]
    for j in i0:length(L)
        s=strip(L[j]); (startswith(s,"C")||isempty(s))&&continue
        seg = j==i0 ? String(L[j][findfirst('/',L[j])+1:end]) : (m=match(r"^\s{5}\S(.*)$",L[j]); m===nothing ? String(L[j]) : String(m.captures[1]))
        ci=findfirst('!',seg); ci!==nothing&&(seg=seg[1:ci-1])
        cl=findfirst('/',seg); if cl!==nothing; push!(body,seg[1:cl-1]); break; else push!(body,seg) end
    end
    [parse(Float64,strip(t)) for t in split(join(body," "),",") if !isempty(strip(t))]
end
forests=["MALHUR","OCHOCO","UMATIL","WALWIT"]
open("data/bluemountains/htdbh_coeffs_bm.csv","w") do io
    println(io,"forest_ifor,species_index,P2,P3,P4")
    for (fi,fn) in enumerate(forests)
        v=dv(fn); @assert length(v)==54 "$fn has $(length(v))"   # 18×3 col-major
        for sp in 1:18
            println(io,join([fi,sp,v[sp],v[18+sp],v[36+sp]],","))   # col-major: P2=v[sp],P3=v[18+sp],P4=v[36+sp]
        end
    end
end
u=dv("UMATIL")
println("wrote htdbh_coeffs_bm.csv; UMATIL DF(3): P2=",u[3]," P3=",u[21]," P4=",u[39])
