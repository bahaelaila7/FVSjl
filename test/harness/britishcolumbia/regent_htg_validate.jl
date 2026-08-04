DIR="/workspace/FVSjl/data/britishcolumbia"
readcsv(f)=[split(l,',') for l in split(strip(replace(read(joinpath(DIR,f),String),"\r"=>"")),'\n')]
rows=readcsv("regent_stcoef.csv"); h=rows[1]; ci(n)=findfirst(==(n),h)
blocks=[]
for r in rows[2:end]
    isempty(strip(r[1])) && continue
    spp=[parse(Int,x) for x in split(r[ci("spp")],';') if x!=""]
    zones=[strip(x) for x in split(r[ci("zones")],';') if strip(x)!=""]
    d=Dict(n=>parse(Float32,strip(r[ci(n)])) for n in ("CON","LNHT","HT","BAL","CCF","CASP","SASP","SLP"))
    push!(blocks,(idx=parse(Int,r[1]),spp=spp,zones=zones,c=d))
end
function resolve(sp,series,zone)
    for pat in (series, zone*"/all ")
        for b in blocks; (sp in b.spp)||continue
            for z in b.zones; occursin(z,pat)&&return b.idx; end
        end
    end; 0
end
const FTtoM=0.3048f0; const MtoFT=3.28084f0; const FT2=0.2295643f0
function main()
    ip=resolve(14,"ICHmw2/01","ICH"); c=blocks[ip].c
    slope=0.30f0; asp=deg2rad(315f0)
    conadj=c["CON"]+c["CASP"]*cos(asp)*slope+c["SASP"]*sin(asp)*slope+c["SLP"]*slope
    println("sp14 ST_COEF ip=$ip CON=$(c["CON"])")
    zrh=[split(l,',') for l in split(strip(read("/workspace/FVSjl/test/harness/britishcolumbia/ref_regent_htg_all_BC_sp14.csv",String)),'\n')]
    n=0;ex=0;mx=0f0
    for r in zrh
        h1=parse(Float32,r[2]);bal=parse(Float32,r[3]);rdj=parse(Float32,r[4]);rhcon=parse(Float32,r[5]);vO=parse(Float32,r[6])
        htm=max(2f0,h1*FTtoM)
        pre=conadj+rhcon+c["LNHT"]*log(htm)+c["HT"]*htm+c["BAL"]*(bal*FT2*0.01f0)+c["CCF"]*rdj
        my=exp(max(-88f0,pre))*MtoFT
        n+=1; my==vO&&(ex+=1); e=abs(my-vO); e>mx&&(mx=e)
    end
    println("n=$n bit-exact=$ex ($(round(100*ex/n,digits=1))%) maxabserr=$mx")
end
main()
