# YEAR LOOP: BWEDR(TOTP)->BWEDEF->BWEAGE->BWEDAM  (include AFTER replay_feeder.jl)
const RECRX=Float32[0.,.2,.8,1.]; const RECRY=Float32[.6,.6,0.,0.]
const RECHST=Float32[0.5,1.0,0.5,0.5,1.0,0.5]
const STHTGR=Float32[-2.3661,-2.4757,-2.0008,-2.3661,-2.9171,0.0]
const DEFOL_HOSTS=[2,5,4]
const DEFP=Float32[85.,80.,70.,60.]

function run_cycle(SITA, IY1,IFINT,FINT)
    FOLPOT=copy(SITA.FOLPOT); FOLADJ=copy(SITA.FOLADJ)
    FNEW=copy(SITA.FNEW); FOLD1=copy(SITA.FOLD1); FOLD2=copy(SITA.FOLD2); FREM=copy(SITA.FREM)
    BWTPHA=SITA.BWTPHA; IFHOST=SITA.IFHOST
    PRBIO=ones(Float32,6,9,4)
    RDDSM1=ones(Float32,6,3); RHTGM1=ones(Float32,6,3)
    PEDDS=zeros(Float32,6,3); PEHTG=zeros(Float32,6,3)
    AVYRMX=zeros(Float32,6,3); BWMXCD=zeros(Float32,6,3)
    CUMDEF=zeros(Float32,6,3,5); APRBYR=zeros(Float32,6,3,2,5); CDEFp=zeros(Float32,6,3)
    NCUMYR=0; ICUMYR=0
    BWFINT = FINT<1.0f0 ? 1.0f0 : FINT
    IBWYR1=IY1; IBWYR2=IY1+IFINT-1
    dumps_age=String[]; dumps_dam=String[]
    for IYRCUR in IY1:(IY1+IFINT-1)
        TOTP=zeros(Float32,6,3)
        for ihost in 1:6
            IFHOST[ihost]==0 && continue
            for icrown in 1:9
                iszi=div(icrown+2,3)
                for ifage in 1:3
                    # match Fortran left-assoc grouping (bwedr.f:126-128)
                    TOTP[ihost,iszi] = (TOTP[ihost,iszi] + FOLPOT[ihost,icrown,ifage]*0.6f0) + FOLADJ[ihost,icrown,ifage]*0.4f0
                end
                TOTP[ihost,iszi] += FREM[icrown,ihost]
            end
        end
        # BWEDEF
        for ih in DEFOL_HOSTS
            IFHOST[ih]==0 && continue
            for icrown in 1:9
                FNEW[icrown,ih]  *= (1.0f0-(DEFP[1]/100.0f0))
                FOLD1[icrown,ih] *= (1.0f0-(DEFP[2]/100.0f0))
                FOLD2[icrown,ih] *= (1.0f0-(DEFP[3]/100.0f0))
                FREM[icrown,ih]  *= (1.0f0-(DEFP[4]/100.0f0))
            end
        end
        # BWEAGE
        TOTR=zeros(Float32,6,3)
        for ihost in 1:6
            IFHOST[ihost]==0 && continue
            for icrown in 1:9
                iszi=div(icrown+2,3)
                DIV=FOLADJ[ihost,icrown,1]; PRBIO[ihost,icrown,1]=0.0f0
                FNEW[icrown,ihost]<0.00001f0 && (FNEW[icrown,ihost]=0.0f0)
                DIV>0.00001f0 && (PRBIO[ihost,icrown,1]=FNEW[icrown,ihost]/DIV)
                DIV=FOLADJ[ihost,icrown,2]; PRBIO[ihost,icrown,2]=0.0f0
                FOLD1[icrown,ihost]<0.00001f0 && (FOLD1[icrown,ihost]=0.0f0)
                DIV>0.00001f0 && (PRBIO[ihost,icrown,2]=FOLD1[icrown,ihost]/DIV)
                DIV=FOLADJ[ihost,icrown,3]; PRBIO[ihost,icrown,3]=0.0f0
                FOLD2[icrown,ihost]<0.00001f0 && (FOLD2[icrown,ihost]=0.0f0)
                DIV>0.00001f0 && (PRBIO[ihost,icrown,3]=FOLD2[icrown,ihost]/DIV)
                DIV=FOLADJ[ihost,icrown,4]; PRBIO[ihost,icrown,4]=0.0f0
                FREM[icrown,ihost]<0.00001f0 && (FREM[icrown,ihost]=0.0f0)
                DIV>0.00001f0 && (PRBIO[ihost,icrown,4]=FREM[icrown,ihost]/DIV)
                if BWTPHA[ihost,iszi]>0.0f0
                    # match Fortran single-expression left-assoc grouping (bweage.f:80-82)
                    TOTR[ihost,iszi]=(((TOTR[ihost,iszi]+FNEW[icrown,ihost])+FOLD1[icrown,ihost])+FOLD2[icrown,ihost])+FREM[icrown,ihost]
                end
            end
        end
        for ihost in 1:6
            IFHOST[ihost]==0 && continue
            for icrown in 1:9
                IREM=1
                PRB=PRBIO[ihost,icrown,1]
                if PRB < 0.80f0
                    IREM=2
                    iszi=div(icrown+2,3)
                    if NCUMYR>0
                        ICMYR=NCUMYR
                        ICUMYR>0 && (ICMYR=ICUMYR)
                        CUM=0.0f0
                        ICMYR>0 && (CUM=CUMDEF[ihost,iszi,ICMYR])
                        CUM>=40.0f0 && (IREM=3)
                    end
                end
                FR=FREM[icrown,ihost]
                PR=FOLPOT[ihost,icrown,4]
                PRALL=PR
                for ifage in 1:3; PRALL+=FOLPOT[ihost,icrown,ifage]; end
                XMULT=0.0f0
                DIV=FOLPOT[ihost,icrown,3]
                DIV>0.000001f0 && (XMULT=FOLPOT[ihost,icrown,4]/DIV)
                FA=FOLD2[icrown,ihost]*XMULT
                if IREM==1
                    DIF=FR-PR
                    if DIF>=0.0f0
                        FA=PR
                        (DIF/PR>0.10f0) && (FA=FR-(DIF*0.40f0))
                        FR=FA
                    else
                        FR=FR+FA
                        FR>PR && (FR=PR)
                    end
                elseif IREM==2
                    FR=FR+FA*bweslp(PRB,RECRX,RECRY,4)*RECHST[ihost]
                    FR>PRALL && (FR=PRALL)
                else
                    FR=(FR+FA*bweslp(PRB,RECRX,RECRY,4))*0.85f0
                    FR>PRALL && (FR=PRALL)
                end
                FREM[icrown,ihost]=FR
                XMULT=0.0f0; DIV=FOLPOT[ihost,icrown,2]
                DIV>0.00001f0 && (XMULT=FOLPOT[ihost,icrown,3]/DIV)
                FOLD2[icrown,ihost]=FOLD1[icrown,ihost]*XMULT
                XMULT=0.0f0; DIV=FOLPOT[ihost,icrown,1]
                DIV>0.00001f0 && (XMULT=FOLPOT[ihost,icrown,2]/DIV)
                FOLD1[icrown,ihost]=FNEW[icrown,ihost]*XMULT
                XMULT=0.0f0; DIV=FOLPOT[ihost,icrown,2]
                DIV>0.00001f0 && (XMULT=FOLD1[icrown,ihost]/DIV)
                FNEW[icrown,ihost]=FOLPOT[ihost,icrown,1]*bweslp((PRBIO[ihost,icrown,2]*0.6f0)+(XMULT*0.4f0),
                    RELFXd[:,icrown],RELFYd[:,icrown],4)
                FAa=FOLADJ[ihost,icrown,4]
                FAa=FAa+FOLADJ[ihost,icrown,3]
                PRp=FOLPOT[ihost,icrown,4]
                FAa>PRp && (FAa=PRp)
                FRp=FREM[icrown,ihost]
                FAa<FRp && (FAa=FRp)
                FOLADJ[ihost,icrown,4]=FAa
                for ii in 1:2
                    ifage=4-ii
                    XMULT=0.0f0; DIV=FOLPOT[ihost,icrown,ifage-1]
                    DIV>0.00001f0 && (XMULT=FOLPOT[ihost,icrown,ifage]/DIV)
                    FOLADJ[ihost,icrown,ifage]=FOLADJ[ihost,icrown,ifage-1]*XMULT
                end
                FOLADJ[ihost,icrown,1]=FNEW[icrown,ihost]
            end
        end
        for ih in 1:6, ic in 1:9, ia in 1:4
            push!(dumps_age, "PRBIO $IYRCUR $ih $ic $ia "*f2h(PRBIO[ih,ic,ia]))
        end
        for ih in 1:6, ic in 1:3
            push!(dumps_age, "TOTR $IYRCUR $ih $ic "*f2h(TOTR[ih,ic]))
        end
        # BWEDAM
        AVPRBO=zeros(Float32,6,3,2); CDEFp.=0.0f0
        STREES=0.0f0
        for ihost in 1:6
            IFHOST[ihost]!=1 && continue
            for iszi in 1:3; STREES+=BWTPHA[ihost,iszi]; end
        end
        if STREES!=0.0f0
            if NCUMYR<5
                NCUMYR+=1; ICUMYR=NCUMYR
            else
                ICUMYR+=1; ICUMYR>5 && (ICUMYR=1)
            end
            for ihost in 1:6
                IFHOST[ihost]==0 && continue
                for iszi in 1:3
                    BWTPHA[ihost,iszi]<=0.0f0 && continue
                    DIV=TOTP[ihost,iszi]; CU=0.0f0
                    DIV>0.00001f0 && (CU=100.0f0-((TOTR[ihost,iszi]/DIV)*100.0f0))
                    CU<0.0f0 && (CU=0.0f0)
                    CUMDEF[ihost,iszi,ICUMYR]=CU
                    APRBYR[ihost,iszi,1,ICUMYR]=PRBIO[ihost,3*iszi-2,1]
                    APRBYR[ihost,iszi,2,ICUMYR]=0.0f0
                end
                for icrown in 1:9
                    iszi=div(icrown+2,3)
                    BWTPHA[ihost,iszi]<=0.0f0 && continue
                    APRBYR[ihost,iszi,2,ICUMYR]+=PRBIO[ihost,icrown,1]*0.3333333f0
                end
            end
            for ihost in 1:6
                IFHOST[ihost]==0 && continue
                for iszi in 1:3
                    BWTPHA[ihost,iszi]<=0.0f0 && continue
                    XMULT=1.0f0/Float32(NCUMYR)
                    for i in 1:NCUMYR
                        CDEFp[ihost,iszi]+=CUMDEF[ihost,iszi,i]
                        for itype in 1:2
                            AVPRBO[ihost,iszi,itype]+=APRBYR[ihost,iszi,itype,i]*XMULT
                        end
                    end
                end
            end
            if NCUMYR>=4
                for ihost in 1:6
                    IFHOST[ihost]==0 && continue
                    for iszi in 1:3
                        BWTPHA[ihost,iszi]<=0.0f0 && continue
                        if AVYRMX[ihost,iszi] < 1.0f0-AVPRBO[ihost,iszi,1]
                            AVYRMX[ihost,iszi]=1.0f0-AVPRBO[ihost,iszi,1]
                        end
                    end
                end
            end
            XFIFTH=0.0f0
            INRUN=IYRCUR-IBWYR1+1
            IFIFTH=INRUN%5
            IFIFTH==0 && (XFIFTH=5.0f0/BWFINT)
            (IFIFTH!=0 && IYRCUR==IBWYR2) && (XFIFTH=Float32(IFIFTH)/BWFINT)
            for ihost in 1:6
                IFHOST[ihost]==0 && continue
                for iszi in 1:3
                    BWTPHA[ihost,iszi]<=0.0f0 && continue
                    RDDS=0.083861f0*fpow(AVPRBO[ihost,iszi,2]*100.0f0, 0.4725f0+0.07f0*RDDSM1[ihost,iszi])*fpow(RDDSM1[ihost,iszi],0.3241f0)
                    RDDS>1.0f0 && (RDDS=1.0f0)
                    RDDSM1[ihost,iszi]=RDDS
                    RHTG=0.193013f0*fpow(AVPRBO[ihost,iszi,1]*100.0f0, 0.3814f0-0.0212f0*RHTGM1[ihost,iszi])*fpow(RHTGM1[ihost,iszi],0.5509f0)
                    RHTG>1.0f0 && (RHTG=1.0f0)
                    RHTGM1[ihost,iszi]=RHTG
                    PEDDS[ihost,iszi]+=RDDS/BWFINT
                    if iszi>1
                        PEHTG[ihost,iszi]+=RHTG/BWFINT
                    else
                        RHTG=1.0f0
                        AVPRBO[ihost,iszi,2]<0.98f0 && (RHTG=fexp(STHTGR[ihost]*(1.0f0-AVPRBO[ihost,iszi,1])))
                        PEHTG[ihost,iszi]+=RHTG*XFIFTH
                    end
                    if INRUN>=5 || IYRCUR==IBWYR2
                        BWMXCD[ihost,iszi]<CDEFp[ihost,iszi] && (BWMXCD[ihost,iszi]=CDEFp[ihost,iszi])
                    end
                end
            end
        end
        for ih in 1:6, is in 1:3
            push!(dumps_dam, "DAM $IYRCUR $ih $is "*f2h(AVPRBO[ih,is,1])*" "*f2h(AVPRBO[ih,is,2])*" "*
                f2h(PEDDS[ih,is])*" "*f2h(PEHTG[ih,is])*" "*f2h(AVYRMX[ih,is])*" "*f2h(CDEFp[ih,is]))
        end
    end
    return (PRBIO=PRBIO,PEDDS=PEDDS,PEHTG=PEHTG,AVYRMX=AVYRMX,dumps_age=dumps_age,dumps_dam=dumps_dam)
end

CYC = run_cycle(SITA, IY1, IFINT, FINT)

age = readlines(joinpath(DIR,"bwe_age.txt"))
ageexp = Dict{String,String}()
for l in age
    s=split(l)
    if s[1]=="PRBIO"; ageexp["PRBIO $(s[2]) $(s[3]) $(s[4]) $(s[5])"]=s[6]
    elseif s[1]=="TOTR"; ageexp["TOTR $(s[2]) $(s[3]) $(s[4])"]=s[5]; end
end
function cmp_age(dumps,ageexp)
    n=0
    for d in dumps
        s=split(d); key=join(s[1:end-1]," "); got=s[end]
        if haskey(ageexp,key) && ageexp[key]!=got
            n+=1; n<=8 && println("  AGE $key got=$got exp=$(ageexp[key])")
        end
    end
    n
end
println("BWEAGE per-year PRBIO+TOTR: mismatches=", cmp_age(CYC.dumps_age,ageexp))

dam = readlines(joinpath(DIR,"bwe_dam.txt"))
damexp=Dict{String,NTuple{6,String}}()
for l in dam
    startswith(l,"DAM ") || continue
    s=split(l); damexp["$(s[2]) $(s[3]) $(s[4])"]=(s[5],s[6],s[7],s[8],s[9],s[10])
end
function cmp_dam(dumps,damexp)
    n=0
    for d in dumps
        s=split(d); key="$(s[2]) $(s[3]) $(s[4])"; got=(s[5],s[6],s[7],s[8],s[9],s[10])
        if haskey(damexp,key) && damexp[key]!=got
            n+=1; n<=12 && println("  DAM $key got=$got exp=$(damexp[key])")
        end
    end
    n
end
println("BWEDAM per-year AVPRBO/PEDDS/PEHTG/AVYRMX/CDEF: mismatches=", cmp_dam(CYC.dumps_dam,damexp))

pdm = readlines(joinpath(DIR,"bwe_pdm.txt"))
function cmp_final(pdm,CYC)
    npe=0
    for l in pdm
        startswith(l,"PDMPE") || continue
        s=split(l); ih=parse(Int,s[2]); is=parse(Int,s[3])
        CYC.PEDDS[ih,is]===h2f(s[4]) || (npe+=1; println("  PEDDS[$ih,$is] got=",f2h(CYC.PEDDS[ih,is])," exp=",s[4]))
        CYC.PEHTG[ih,is]===h2f(s[5]) || (npe+=1; println("  PEHTG[$ih,$is] got=",f2h(CYC.PEHTG[ih,is])," exp=",s[5]))
        CYC.AVYRMX[ih,is]===h2f(s[6]) || (npe+=1; println("  AVYRMX[$ih,$is] got=",f2h(CYC.AVYRMX[ih,is])," exp=",s[6]))
    end
    nprb=0
    for l in pdm
        startswith(l,"PDMPRB") || continue
        s=split(l); ih=parse(Int,s[2]); ic=parse(Int,s[3]); ia=parse(Int,s[4])
        CYC.PRBIO[ih,ic,ia]===h2f(s[5]) || (nprb+=1; nprb<=8 && println("  PRBIO[$ih,$ic,$ia] got=",f2h(CYC.PRBIO[ih,ic,ia])," exp=",s[5]))
    end
    (npe,nprb)
end
(npe,nprb)=cmp_final(pdm,CYC)
println("FINAL PEDDS/PEHTG/AVYRMX vs PDM: mismatches=$npe")
println("FINAL PRBIO vs PDM: mismatches=$nprb")
