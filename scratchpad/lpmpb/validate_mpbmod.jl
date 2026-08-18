# MPBMOD year-loop replay (lpmpb/mpbmod.f, deterministic epidemic) — LGO path, LAGG/LREP/LPS/LDC=F.
# All params from fort.785; per-class inputs from fort.785 'C' recs. Validates final CLASS(I,IMPROB)=
# TREES(I) survivors vs fort.784 CLPOP (target: total kill, survival << 3.6e-7 cap threshold).

include("betin_pkg.jl")                       # BETIN package (extracted from validate_betin.jl)

h2d(s) = reinterpret(Float64, parse(UInt64, s, base=16))
h2f(s) = reinterpret(Float32, parse(UInt32, s, base=16))

# --- load params (fort.785) ---
P = Dict{String,Float32}(); IP = Dict{String,Int}()
G_DST = Float32[]; G_EXO = Float32[]
TREES0 = Float32[]; SURF = Float32[]; DIAM = Float32[]; PHLOEM = Float32[]
for ln in eachline("/workspace/.iework/lpmpb/run/fort.785")
    t = split(ln)
    if t[1] == "IPARM"
        IP["NG"]=parse(Int,t[2]); IP["INCRS"]=parse(Int,t[3]); IP["IB"]=parse(Int,t[4])
        IP["MPMXYR"]=parse(Int,t[5]); IP["NACLAS"]=parse(Int,t[6])
    elseif t[1] == "P";  P[t[2]] = h2f(t[3])
    elseif t[1] == "G";  push!(G_DST, h2f(t[4])); push!(G_EXO, h2f(t[5]))  # G DST <IG> <DST> <EXODUS>
    elseif t[1] == "C"
        push!(TREES0, h2f(t[3])); push!(SURF, h2f(t[4])); push!(DIAM, h2f(t[5])); push!(PHLOEM, h2f(t[6]))
    end
end
NG=IP["NG"]; INCRS=IP["INCRS"]; IB=IP["IB"]; MPMXYR=IP["MPMXYR"]; NACLAS=IP["NACLAS"]
CE=P["CE"]; EXCON=P["EXCON"]; SQFTPA=P["SQFTPA"]; STRBUG=P["STRBUG"]; STRP=P["STRP"]
SEXRAT0=P["SEXRAT"]; HS=P["HS"]; CF1=P["CF1"]; CF2=P["CF2"]; CF3=P["CF3"]
TA=P["TA"]; AMP1=P["AMP1"]; AMP2=P["AMP2"]; CRITAD=P["CRITAD"]
EFELEV=P["EFELEV"]; EFLAT=P["EFLAT"]; DST=G_DST; EXODUS=G_EXO

TPROB = sum(TREES0)                            # total LP PROB
SEXDBH = @. 0.918f0 - 0.0168f0*DIAM
EFPHLM = @. max(0.0f0, 16.67f0*PHLOEM - 0.667f0)

# EMERG (emerg.f) with persistent C threaded explicitly (the -fno-automatic static local).
mutable struct EmergState; C::Float64; end
function emerg!(st, by::Float32, t::Int, incrs::Int)::Float32
    if t == 0
        st.C = 1.0/2.0^incrs
    else
        st.C = st.C*(incrs - t + 1)/t
    end
    return by*Float32(st.C)
end

# --- MPBMOD year loop (LGO path) — wrapped in a function to avoid top-level soft-scope traps ---
function run_mpbmod(TREES0, SURF, DIAM, PHLOEM, SEXDBH, EFPHLM, TPROB;
                    NG, INCRS, IB, MPMXYR, NACLAS, CE, EXCON, SQFTPA, STRBUG, STRP,
                    SEXRAT0, HS, CF1, CF2, CF3, TA, AMP1, AMP2, CRITAD, EFELEV, EFLAT, DST, EXODUS)
TREES = copy(TREES0)
BY = STRBUG; Pp = STRP; SEXRAT = SEXRAT0
MPBYR = 0; SADLPP = 0.0f0; SSCUM = SADLPP
NINC = INCRS + 1
AGG = zeros(Float32, NACLAS, NINC)
AD  = zeros(Float32, NINC)
AMP = zeros(Float32, NINC)
BOLD = zeros(Float32, NG)

while true
    MPBYR += 1
    Q = 1.0f0 - Pp
    GENO = NG == 3 ? Float32[Pp*Pp, 2*Pp*Q, Q*Q] : Float32[Pp, Q]
    B3SUM = zeros(Float32, NG)
    fill!(BOLD, 0.0f0)
    OS = 0.0f0
    for I in 1:NACLAS; OS += SURF[I]*TREES[I]; end
    OS += SADLPP                              # LGO
    fill!(AGG, 0.0f0); fill!(AD, 0.0f0); fill!(AMP, 0.0f0)
    est = EmergState(0.0)
    # ---- EMERGENCE loop ----
    for INC in 1:NINC
        INC1 = INC - 1
        E1 = 0.0f0
        for I in 1:NACLAS; E1 += SURF[I]*TREES[I]; end
        E3 = 0.0f0; EF3 = 0.0f0; TAGG = 0.0f0
        if INC > 1
            for KK in IB:INC1
                AMP[KK] = AMP1 - AMP2*AD[KK]*(1.0f0 - SEXRAT)
                AMP[KK] < 0.0f0 && (AMP[KK] = 0.0f0)
                for I in 1:NACLAS
                    TAGG += AGG[I,KK]; E3 += SURF[I]*AGG[I,KK]; EF3 += SURF[I]*AGG[I,KK]*AMP[KK]
                end
            end
        end
        RHO1 = 0.0f0
        for I in 1:NACLAS; RHO1 += TREES[I]; end
        E2 = OS - (E1 + E3) + 0.0f0            # SNOHST=0
        EFFS = E1 + E2 + EF3
        RHO3 = TAGG/SQFTPA
        RHO2 = TPROB - TAGG - RHO1; RHO2 < 0.0f0 && (RHO2 = 0.0f0); RHO2 /= SQFTPA
        RHO1 /= SQFTPA
        BNEW = emerg!(est, BY, INC1, INCRS)
        B1INC = 0.0f0; B3INC = 0.0f0
        for IG in 1:NG
            B0 = GENO[IG]*BNEW + BOLD[IG]
            EXLOSS = B0*EXODUS[IG]
            B0 -= EXLOSS
            B1 = B0*E1/EFFS; B2 = B0*E2/EFFS; B3 = B0*EF3/EFFS
            fm(cf,rho) = (fmx = -cf*2.0f0*sqrt(rho)*DST[IG]*DST[IG]/DST[1]; fmx > -80.0f0 ? exp(fmx) : 0.0f0)
            B1 -= B1*fm(CF1,RHO1); B2 -= B2*fm(CF2,RHO2); B3 -= B3*fm(CF3,RHO3)
            BOLD[IG] = B1 + B2
            B1INC += B1; B3INC += B3; B3SUM[IG] += B3
        end
        PIODEN = Float64(B1INC/E1)                 # gfortran: REAL division then widen to DOUBLE
        XX = 1.0 - exp(-PIODEN)
        for I in 1:NACLAS
            AGG[I,INC] = 0.0f0
            DTA = Float64(TA); DSMTA = Float64(SURF[I]) - DTA + 1.0
            if XX != 0.0 && DSMTA > 0.0
                AGG[I,INC] = Float32(Float64(TREES[I])*betin(DTA, DSMTA, XX))  # REAL·DOUBLE→DOUBLE then round
            end
            TREES[I] -= AGG[I,INC]
        end
        if EF3 > 0.0f0
            for KK in IB:INC1; AD[KK] += AMP[KK]*B3INC/EF3; end
        end
    end
    # ---- PRODUCTIVITY loop ----
    BY = 0.0f0; BYT = 0.0f0
    for INC in 1:NINC
        XEG = -0.117f0*AD[INC] > -80.0f0 ? exp(-0.117f0*AD[INC]) : 0.0f0
        EGGS = 630.0f0*(1.0f0 - XEG)
        for I in 1:NACLAS
            PSURV = 1.0f0 - exp(-AD[INC]*0.04328f0)
            AD[INC] >= 2.595f0 && (PSURV = PSURV*6.812f0*exp(-1.191f0*sqrt(AD[INC])))
            YOUNG = EGGS*PSURV*EFELEV*EFLAT*EFPHLM[I]
            if AD[INC] >= CRITAD
                TRKILL = AGG[I,INC]
                SKILL = SURF[I]*TRKILL
                BY += YOUNG*SKILL*HS*SEXDBH[I]; BYT += YOUNG*SKILL*HS
            else
                TREES[I] += AGG[I,INC]; AGG[I,INC] = 0.0f0
            end
        end
    end
    BYT > 0.0f0 && (SEXRAT = BY/BYT)
    B3ALL = sum(B3SUM)
    if B3ALL != 0.0f0
        Pp = NG == 3 ? (B3SUM[1] + 0.5f0*B3SUM[2])/B3ALL : (B3SUM[1]/B3ALL)^2
    end
    println("yr $MPBYR: BY=$(round(BY,digits=3))  ΣTREES=$(round(sum(TREES),digits=4))  P=$(round(Pp,digits=4))")
    (BY >= 1.0f0 && MPBYR < MPMXYR) || break
end
return TREES
end

TREES = run_mpbmod(TREES0, SURF, DIAM, PHLOEM, SEXDBH, EFPHLM, TPROB;
    NG, INCRS, IB, MPMXYR, NACLAS, CE, EXCON, SQFTPA, STRBUG, STRP,
    SEXRAT0, HS, CF1, CF2, CF3, TA, AMP1, AMP2, CRITAD, EFELEV, EFLAT, DST, EXODUS)

# --- compare final survivors vs fort.784 CLPOP (col2 = post-MPBMOD CLASS(I,IMPROB)) ---
orc = Float32[]
for ln in eachline("/workspace/.iework/lpmpb/run/fort.784")
    t = split(ln); t[1] == "CLPOP" || continue
    push!(orc, h2f(t[4]))  # CLPOP <I> <pre=SURVIV> <post=CLASS(I,IMPROB)>; col2=t[4]=survivor
end
println("\n--- final class survivors TREES(I) vs oracle ---")
maxsurv = 0.0f0
for I in 1:NACLAS
    surv_ratio = TREES[I] / TREES0[I]
    global maxsurv = max(maxsurv, abs(surv_ratio))
    println("cls $I  jl=$(TREES[I])  orc=$(orc[I])   survival=$(surv_ratio)")
end
println("\nmax class survival ratio = $maxsurv  (cap threshold 3.6e-7)")
println(maxsurv < 3.6e-7 ? "*** TOTAL KILL — WK2(I)=PROB(I)-1e-6 for all, .sum bit-exact ***" :
        "!!! survival above cap threshold — need exact match")
