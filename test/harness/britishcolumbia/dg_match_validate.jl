# Standalone validation of the BC V3 PrettyName->IP/JP matcher vs captured oracle ground truth.
const DIR="/workspace/FVSjl/data/britishcolumbia"
readcsv(f)=[split(l,',') for l in split(strip(read(joinpath(DIR,f),String)),'\n')]

function loadblocks(fname, zonecol)
    rows=readcsv(fname); hdr=rows[1]
    idxc=findfirst(==("idx"),hdr); sppc=findfirst(==("spp"),hdr); znc=findfirst(==(zonecol),hdr)
    blocks=[]
    for r in rows[2:end]
        isempty(r[1]) && continue
        spp=[parse(Int,x) for x in split(r[sppc],';') if x!=""]
        zones=[strip(x) for x in split(r[znc],';') if strip(x)!=""]
        push!(blocks,(idx=parse(Int,r[idxc]), spp=spp, zones=zones))
    end
    blocks
end
zn=loadblocks("dg_znkonst.csv","zones")
ss=loadblocks("dg_sskonst.csv","prettynames")

# matcher: find first block (in idx order) whose spp contains I and some zone/prettyname is substring of PAT
# two passes: PAT=series, then PAT=zone*"/all "
function resolve(blocks, sp, series, zone)
    for pat in (series, zone*"/all ")
        for b in blocks
            (sp in b.spp) || continue
            for z in b.zones
                if occursin(z, pat)   # INDEX(PAT, prettyname)>0  == prettyname is substring of PAT
                    return b.idx
                end
            end
        end
    end
    return -1
end

# ground truth captured from instrumented oracle (all_BC = ICH / ICHmw2/01)
truth = Dict(1=>(1,1),2=>(2,2),3=>(5,3),4=>(17,21),5=>(10,5),6=>(11,7),7=>(13,11),
             8=>(15,15),9=>(17,21),10=>(18,24),11=>(19,25),12=>(21,25),13=>(21,25),14=>(5,3),15=>(19,25))
function run(zn,ss,truth)
    series="ICHmw2/01"; zone="ICH"
    println("sp | jl IP/JP | truth IP/JP | ok")
    allok=true
    for sp in 1:15
        ip=resolve(zn,sp,series,zone); jp=resolve(ss,sp,series,zone)
        t=truth[sp]; ok=(ip==t[1] && jp==t[2]); allok &= ok
        println("$sp | $ip/$jp | $(t[1])/$(t[2]) | $(ok ? "✓" : "✗")")
    end
    println(allok ? "\nALL 15 MATCH — matcher validated" : "\nMISMATCH")
end
run(zn,ss,truth)
