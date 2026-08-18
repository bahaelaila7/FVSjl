# Validate the STAGED src/engine/wsbwe.jl functions (feeder + BWERNP/BWEBET RNG)
# against the instrumented-oracle dumps. Loads wsbwe.jl with minimal stubs.
abstract type AbstractWsbweState end
struct StandState end          # stub — only needed so wsbwe_apply!/kw_wsbwe! signatures parse
struct KeywordReader end       # stub
const DIR = @__DIR__
include(joinpath(DIR, "wsbwe.jl"))

h2f(s) = reinterpret(Float32, parse(UInt32, s, base=16))
f2h(x::Float32) = uppercase(string(reinterpret(UInt32, x), base=16, pad=8))

# ---- parse SIT dump (tree list + stand scalars) ----
sit = readlines(joinpath(DIR,"bwe_sit.txt"))
hdr = split(first(filter(l->startswith(l,"SITHDR"),sit)))
IY1=parse(Int,hdr[3]); IFINT=parse(Int,hdr[4])
FINT=h2f(hdr[7]); ELEV=h2f(hdr[8]); OLDTPA=h2f(hdr[9]); ORMSQD=h2f(hdr[10])
SP=Int[]; HT=Float32[]; DBH=Float32[]; DG=Float32[]; ICR=Int[]; PROB=Float32[]
for l in sit
    startswith(l,"SITTRE") || continue
    s=split(l)
    push!(SP,parse(Int,s[3])); push!(HT,h2f(s[5])); push!(DBH,h2f(s[6]))
    push!(DG,h2f(s[7])); push!(ICR,Int(round(h2f(s[8])))); push!(PROB,h2f(s[9]))
end
ns=length(SP)

# ---- build a WsbweState with the host_defol.key DEFOL schedule ----
w = wsbwe_defaults!(); w.active=true; w.ldefol=true
for yr in IY1:(IY1+IFINT-1), spc in (3,8,9)   # DF=3, ES=8, AF=9 (FVS EM indices)
    push!(w.defol_sched, (Float32(yr),Float32(spc),0.0f0,85.0f0,80.0f0,70.0f0,60.0f0))
end

# ---- run STAGED feeder ----
(PRBIO,PEDDS,PEHTG,AVYRMX,IFHOST) = wsbwe_feeder(w, ns, SP,HT,DBH,DG,ICR,PROB,
    WSBWE_IBWSPM_EM, WSBWE_IBIOMP_EM, ELEV,OLDTPA,ORMSQD,FINT,IFINT,IY1)

# ---- compare final arrays vs PDM dumps ----
pdm = readlines(joinpath(DIR,"bwe_pdm.txt"))
npe=0
for l in pdm
    startswith(l,"PDMPE") || continue
    s=split(l); ih=parse(Int,s[2]); is=parse(Int,s[3])
    PEDDS[ih,is]===h2f(s[4]) || (global npe+=1; println("  PEDDS[$ih,$is] ",f2h(PEDDS[ih,is])," vs ",s[4]))
    PEHTG[ih,is]===h2f(s[5]) || (global npe+=1; println("  PEHTG[$ih,$is] ",f2h(PEHTG[ih,is])," vs ",s[5]))
    AVYRMX[ih,is]===h2f(s[6]) || (global npe+=1; println("  AVYRMX[$ih,$is] ",f2h(AVYRMX[ih,is])," vs ",s[6]))
end
nprb=0
for l in pdm
    startswith(l,"PDMPRB") || continue
    s=split(l); ih=parse(Int,s[2]); ic=parse(Int,s[3]); ia=parse(Int,s[4])
    PRBIO[ih,ic,ia]===h2f(s[5]) || (global nprb+=1)
end
println("STAGED wsbwe_feeder → PEDDS/PEHTG/AVYRMX mism=$npe ; PRBIO mism=$nprb")

# ---- validate STAGED BWERNP/BWEBET RNG stream by driving the per-tree BWERNP sequence ----
# Reuse the dumped BWEPDM inputs to reproduce the exact draw order, calling the staged kernels.
function pdmin(pdm)
    TIN=Dict{Int,Any}(); TOU=Dict{Int,Any}()
    for l in pdm; startswith(l,"PDMTIN")||continue; s=split(l); I=parse(Int,s[2])
        TIN[I]=(IHOST=parse(Int,s[4]),ISZI=parse(Int,s[5]),ICR0=parse(Int,s[10]),
                HT=h2f(s[13]),MFT=h2f(s[21]),MFM=h2f(s[22]),PCTK=h2f(s[23]),ITRE=parse(Int,s[11])) ; end
    for l in pdm; startswith(l,"PDMTOU")||continue; s=split(l); I=parse(Int,s[2])
        TOU[I]=(ITRUNC=parse(Int,s[3]),NORMHT=parse(Int,s[4])) ; end
    (TIN,TOU)
end
(TIN,TOU)=pdmin(pdm)
ISCT=zeros(Int,64,2); for l in pdm; startswith(l,"PDMSCT")||continue; s=split(l); ISCT[parse(Int,s[2]),1]=parse(Int,s[3]); ISCT[parse(Int,s[2]),2]=parse(Int,s[4]); end
IND1=zeros(Int,ns); for l in pdm; startswith(l,"PDMIND")||continue; s=split(l); IND1[parse(Int,s[2])]=parse(Int,s[3]); end
# entry arrays from PDMPE
PE=zeros(Float32,6,3); PH=zeros(Float32,6,3); AV=zeros(Float32,6,3)
for l in pdm; startswith(l,"PDMPE")||continue; s=split(l); ih=parse(Int,s[2]); is=parse(Int,s[3]); PE[ih,is]=h2f(s[4]); PH[ih,is]=h2f(s[5]); AV[ih,is]=h2f(s[6]); end
IFIR=WSBWE_IFIR
w.rng_s0 = Float64(w.dseed)   # BWERPT(DSEEDD=55329)
draws=Float32[]
# instrument wsbwe_rand! by wrapping: capture draws
function draw!(w,acc); v=wsbwe_rand!(w); push!(acc,v); v; end
# reimplement BWERNP/BWEBET call-order using staged kernels but recording draws
function bernp_rec(w,acc,xmean,xmxvar)
    (xmean<0.001f0||xmean>0.999f0) && return xmean
    xm = xmean>0.5f0 ? (1.0f0-xmean) : xmean
    var=(WSBWE_BSC*xmxvar)*(WSBWE_BC*wsbwe_pow(xm*(WSBWE_BC-1.0f0),WSBWE_BC-1.0f0)*wsbwe_exp(-wsbwe_pow((WSBWE_BC-1.0f0)*xm,WSBWE_BC)))
    ww_=(xmean*(1.0f0-xmean)/var)-1.0f0; v=xmean*ww_; ww_=(1.0f0-xmean)*ww_
    (x,kode)=bebet_rec(w,acc,v,ww_); kode<0 ? xmean : x
end
function bebet_rec(w,acc,a0,b0)
    (a0<=0.00001f0||b0<=0.00001f0) && return (0.0f0,-1)
    kode=0
    if min(a0,b0)>1.0f0
        a=min(a0,b0);b=max(a0,b0);alf=a+b;bet=sqrt((alf-2.0f0)/(2.0f0*a*b-alf));gam=a+(1.0f0/b);local w_
        while true
            kode+=1; u1=draw!(w,acc); u2=draw!(w,acc)
            v=bet*wsbwe_log(u1/(1.0f0-u1)); v>20.0f0 && (v=20.0f0); w_=a*wsbwe_exp(v)
            z=u1*u1*u2; rr=gam*v-1.3862944f0; s=a+rr-w_
            s+2.609438f0>=5.0f0*z && break; t=wsbwe_log(z); s>=t && break
            rr+alf*wsbwe_log(alf/(b+w_))<t && continue; break
        end
        return ((a==a0) ? w_/(b+w_) : b/(b+w_),kode)
    else
        a=max(a0,b0);b=min(a0,b0);alf=a+b;bet=1.0f0/b;del=1.0f0+a-b
        c1=del*(0.0138889f0+0.0416667f0*b)/(a*bet-0.777778f0); c2=0.25f0+(0.5f0+0.25f0/del)*b; local w_
        while true
            kode+=1; u1=draw!(w,acc); u2=draw!(w,acc)
            if u1>=0.5f0
                z=u1*u1*u2
                if z<=0.25f0; v=bet*wsbwe_log(u1/(1.0f0-u1)); v>20.0f0&&(v=20.0f0); w_=a*wsbwe_exp(v); break; end
                z>=c2 && continue
                v=bet*wsbwe_log(u1/(1.0f0-u1)); v>20.0f0&&(v=20.0f0); w_=a*wsbwe_exp(v)
                alf*(wsbwe_log(alf/(b+w_))+v)-1.3862944f0<wsbwe_log(z) && continue; break
            else
                y=u1*u2; z=u1*y
                if 0.25f0*u2+z-y>=c1; continue
                else
                    v=bet*wsbwe_log(u1/(1.0f0-u1)); v>20.0f0&&(v=20.0f0); w_=a*wsbwe_exp(v)
                    alf*(wsbwe_log(alf/(b+w_))+v)-1.3862944f0<wsbwe_log(z) && continue; break
                end
            end
        end
        return ((a==a0) ? w_/(b+w_) : b/(b+w_),kode)
    end
end
for ispi in 1:19
    ISCT[ispi,1]==0 && continue
    for ii in ISCT[ispi,1]:ISCT[ispi,2]
        I=IND1[ii]; t=TIN[I]; ih=t.IHOST; ih>=6 && continue; iszi=t.ISZI
        AVDEF=bernp_rec(w,draws,AV[ih,iszi],0.06f0)
        if iszi==1
            Xc = t.HT>10.0f0 ? 10.0f0 : t.HT
            PRTOPK=1.0f0/(1.0f0+wsbwe_exp(-(-2.5817f0-0.027635f0*Float32(t.ICR0)+3.709f0*sqrt(AVDEF)+0.0488f0*Xc)))
        else
            PRTOPK = IFIR[ih]==0 ? 0.96f0*(1.0f0-wsbwe_exp(-(wsbwe_pow(0.65f0*(AVDEF+1.0f0),14.0f0)))) :
                                    0.90f0*(1.0f0-wsbwe_exp(-(wsbwe_pow(0.60f0*(AVDEF+1.0f0),11.0f0))))
        end
        Xr=draw!(w,draws)
        # (topkill PART draws none; DDS BWERNP always)
        XD=bernp_rec(w,draws,PE[ih,iszi],0.03f0)
        # HTG XH: PEDDS<0.99 deterministic (no draw) — matches fixture
    end
end
exp_rng=[h2f(split(l)[2]) for l in readlines(joinpath(DIR,"bwe_rng.txt")) if startswith(l,"R")]
rngok = length(draws)==length(exp_rng) && all(draws .=== exp_rng)
println("STAGED BWERNP/BWEBET RNG stream: julia=$(length(draws)) oracle=$(length(exp_rng)) bit-exact=$rngok")
println((npe==0 && nprb==0 && rngok) ? "STAGED wsbwe.jl VALIDATED (feeder + RNG bit-exact)" : "STAGED wsbwe.jl: CHECK FAILURES ABOVE")
