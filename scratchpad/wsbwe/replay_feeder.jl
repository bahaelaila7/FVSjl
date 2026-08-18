# Dump-replay validation of the WSBWE FEEDER chain (BWESIT->BWEBMS->BWEADV->BWEDEF->
# BWEAGE->BWEDAM) vs instrumented FVSem_wsbwe (host_defol.key, cycle 1). Float32-hex.
const DIR = @__DIR__
h2f(s) = reinterpret(Float32, parse(UInt32, s, base=16))
f2h(x::Float32) = uppercase(string(reinterpret(UInt32, x), base=16, pad=8))
const LIBM = "libm.so.6"
@inline flog(x::Float32) = ccall((:logf, LIBM), Float32, (Float32,), x)
@inline fexp(x::Float32) = ccall((:expf, LIBM), Float32, (Float32,), x)
@inline fpow(x::Float32, y::Float32) = ccall((:powf, LIBM), Float32, (Float32,Float32), x, y)
ifix(x::Float32) = trunc(Int, x)

# ---- EM block data (bwebkem.f) ----
const IBWSPM = Int[7,6,2,7,7,7,7,5,4,7,7,7,7,7,7,7,7,7,7]
const PRCRN3 = Float32[0.05,0.3,0.65,0.15,0.45,0.40,0.15,0.45,0.40]
# THEOFL(4,9): Fortran column-major DATA fills (age fastest). THEOFL[age,crown]
const THEOFL = permutedims(reshape(Float32[
 0.70,0.24,0.04,0.02, 0.45,0.30,0.20,0.05, 0.35,0.30,0.25,0.10,
 0.50,0.35,0.12,0.03, 0.30,0.25,0.20,0.25, 0.10,0.20,0.20,0.50,
 0.50,0.30,0.15,0.05, 0.30,0.25,0.20,0.25, 0.10,0.20,0.20,0.50], (4,9)))  # [crown,age] after permute? see below
# We'll index THEOFL as THEO[age,crown]; build directly:
const THEO = reshape(Float32[
 0.70,0.24,0.04,0.02, 0.45,0.30,0.20,0.05, 0.35,0.30,0.25,0.10,
 0.50,0.35,0.12,0.03, 0.30,0.25,0.20,0.25, 0.10,0.20,0.20,0.50,
 0.50,0.30,0.15,0.05, 0.30,0.25,0.20,0.25, 0.10,0.20,0.20,0.50], (4,9))  # THEO[age,crown]
const RELFXd = reshape(Float32[
 0.0,0.04,0.95,1.0, 0.0,0.04,0.95,1.0, 0.0,0.04,0.95,1.0,
 0.0,0.04,0.70,1.0, 0.0,0.12,0.80,1.0, 0.0,0.15,0.90,1.0,
 0.0,0.04,0.70,1.0, 0.0,0.12,0.80,1.0, 0.0,0.15,0.90,1.0], (4,9))  # RELFX[k,crown]
const RELFYd = reshape(Float32[
 0.05,0.05,1.0,1.0, 0.05,0.05,1.0,1.0, 0.05,0.05,1.0,1.0,
 0.05,0.05,1.0,1.0, 0.15,0.15,1.0,1.0, 0.20,0.20,1.0,1.0,
 0.05,0.05,1.0,1.0, 0.15,0.15,1.0,1.0, 0.20,0.20,1.0,1.0], (4,9))

# BWEBMS EM (bwebmsem.f)
const IBIOMP = Int[1,2,3,4,2,11,7,8,9,10,11,11,11,11,11,11,11,11,11]
const BINT2 = Float32[2.666072,1.756537,2.705866,3.115084,2.654572,3.059351,2.622505,3.300852,3.060169,2.452492,2.622505]
const BINT12 = Float32[-1.94951,-4.73762,-2.05828,-2.43200,-4.17456,-2.24876,-3.13488,-2.93508,-1.60998,-2.74410,-2.63387]
const BCL12 = Float32[1.22023,1.98479,1.25837,1.60270,2.00749,1.37600,1.62368,1.96125,1.32649,1.58171,1.35092]

# BWESLP linear interpolation
function bweslp(xx::Float32, x::Vector{Float32}, y::Vector{Float32}, n::Int)::Float32
    for i in 1:n-1
        (xx < x[i] || xx > x[i+1]) && continue
        return y[i] + ((y[i+1]-y[i])/(x[i+1]-x[i]))*(xx-x[i])
    end
    r = y[n]
    xx < x[1] && (r = y[1])
    return r
end

# BWECRC
function bwecrc(xht::Float32)
    ihtc = 1
    if xht >= 23.0f0; ihtc = 2; end
    if xht >= 46.0f0 && ihtc==2; ihtc = 3; end
    # replicate exact fortran: ihtc=1; if xht<23 goto10; ihtc=2; if xht<46 goto10; ihtc=3
    ihtc = 1
    if !(xht < 23.0f0)
        ihtc = 2
        if !(xht < 46.0f0)
            ihtc = 3
        end
    end
    icrc2 = ihtc*3; icrc1 = icrc2-2
    return (ihtc, icrc1, icrc2)
end

# ---------------- parse SIT dump ----------------
sit = readlines(joinpath(DIR,"bwe_sit.txt"))
hdr = split(first(filter(l->startswith(l,"SITHDR"),sit)))
ICYC=parse(Int,hdr[2]); IY1=parse(Int,hdr[3]); IFINT=parse(Int,hdr[4]); ITRN=parse(Int,hdr[5]); IPTINV=parse(Int,hdr[6])
FINT=h2f(hdr[7]); ELEV=h2f(hdr[8]); OLDTPA=h2f(hdr[9]); ORMSQD=h2f(hdr[10]); HOSTST_or=h2f(hdr[11])
# tree list
struct Tree; I::Int; ISP::Int; ITRE::Int; HT::Float32; DBH::Float32; DG::Float32; ICR::Int; PROB::Float32; WK4::Float32; end
trees = Tree[]
for l in sit
    startswith(l,"SITTRE") || continue
    s=split(l)
    push!(trees, Tree(parse(Int,s[2]),parse(Int,s[3]),parse(Int,s[4]),
        h2f(s[5]),h2f(s[6]),h2f(s[7]), Int(round(h2f(s[8]))), h2f(s[9]), h2f(s[10])))
end
# IY(ICYC+1) and IY(1): cycle 1 -> end = IY1+IFINT = 2000; inv = IY1 = 1990
const IYP1 = IY1 + IFINT   # IY(ICYC+1)
const IYINV = IY1          # IY(1)

# ---------------- Step 1: BWEBMS (ICVOPT=2) validate WK4 ----------------
alntpa = flog(OLDTPA)
nbad = 0
for t in trees
    D = t.DBH; ispi = IBIOMP[t.ISP]; H = t.HT
    CL = (Float32(t.ICR)*H)/100.0f0
    RD = D/ORMSQD
    DDS = (2.0f0*D*t.DG + t.DG*t.DG)/FINT
    if D < 3.5f0
        wk4 = fexp(BINT12[ispi] + BCL12[ispi]*flog(CL) - 0.12975f0*alntpa + 0.40350f0*flog(H))
        wk4 = wk4 * 1.13178f0
    else
        wk4 = fexp(BINT2[ispi] + 1.468547f0*flog(D) + 0.308847f0*flog(DDS) - 1.077047f0*flog(H) +
              0.690825f0*flog(CL) - 0.142096f0*alntpa + 0.399244f0*flog(RD))
    end
    if wk4 !== t.WK4
        nbad += 1
        nbad <= 6 && println("  WK4 mismatch tree $(t.I) sp$(t.ISP) D=$D got=",f2h(wk4)," exp=",f2h(t.WK4))
    end
end
println("BWEBMS WK4: $(length(trees)-nbad)/$(length(trees)) bit-exact")

# ---------------- Step 2: BWESIT accumulation ----------------
function bwesit_accum(trees)
FOLPOT = zeros(Float32,6,9,4)
FOLADJ = zeros(Float32,6,9,4)
FOLNH  = zeros(Float32,9)
BWTPHA = zeros(Float32,7,3)
IFHOST = zeros(Int,7)
for t in trees
    IFHOST[IBWSPM[t.ISP]] = 1
end
HOSTST = 0.0f0
for t in trees
    (ihtc,icrc1,icrc2) = bwecrc(t.HT)
    iszi = ihtc
    PROBI = t.PROB*2.47103f0
    ihost = IBWSPM[t.ISP]
    BWTPHA[ihost,iszi] += PROBI
    BIO = t.WK4*453.6f0*PROBI
    if ihost < 7
        for icrown in icrc1:icrc2, ifage in 1:4
            FOLPOT[ihost,icrown,ifage] += THEO[ifage,icrown]*PRCRN3[icrown]*BIO
        end
    else
        for icrown in icrc1:icrc2
            FOLNH[icrown] += PRCRN3[icrown]*BIO
        end
    end
end
# HOSTST (loop 267): sum BWTPHA*0.4047 for ihost != 7 present
for ihost in 1:7
    IFHOST[ihost]==0 && continue
    for iszi in 1:3
        ihost != 7 && (HOSTST += BWTPHA[ihost,iszi]*0.4047f0)
    end
end
# average FOLPOT (loop 330)
for ihost in 1:6
    IFHOST[ihost]==0 && continue
    for icrown in 1:9
        iszi = div(icrown+2,3)
        DIV = BWTPHA[ihost,iszi]
        if DIV > 0.0f0
            for ifage in 1:4
                FOLPOT[ihost,icrown,ifage] /= DIV
            end
        else
            for ifage in 1:4; FOLPOT[ihost,icrown,ifage]=0.0f0; end
        end
    end
end
# average FOLNH (loop 350)
for icrown in 1:9
    iszi = div(icrown+2,3)
    DIV = BWTPHA[7,iszi]
    FOLNH[icrown] = DIV>0.0f0 ? FOLNH[icrown]/DIV : 0.0f0
end
# FOLADJ = FOLPOT*POFPOT (POFPOT=1)
for ihost in 1:6
    IFHOST[ihost]==0 && continue
    for icrown in 1:9, ifage in 1:4
        FOLADJ[ihost,icrown,ifage] = FOLPOT[ihost,icrown,ifage]*1.0f0
    end
end
# BWEADV: larch(6) fold into nonhost then IFHOST(6)=0; then load FNEW/FOLD1/FOLD2/FREM
if IFHOST[6]==1
    for icrown in 1:9
        FOLNH[icrown] += FOLADJ[6,icrown,1]+FOLADJ[6,icrown,2]+FOLADJ[6,icrown,3]+FOLADJ[6,icrown,4]
    end
    IFHOST[6]=0
end
FNEW=zeros(Float32,9,6); FOLD1=zeros(Float32,9,6); FOLD2=zeros(Float32,9,6); FREM=zeros(Float32,9,6)
PRBIO=ones(Float32,6,9,4)   # init 1.0
for ihost in 1:6
    IFHOST[ihost]==0 && continue
    for icrown in 1:9
        FNEW[icrown,ihost]  = FOLADJ[ihost,icrown,1]
        FOLD1[icrown,ihost] = FOLADJ[ihost,icrown,2]*PRBIO[ihost,icrown,2]
        FOLD2[icrown,ihost] = FOLADJ[ihost,icrown,3]*PRBIO[ihost,icrown,3]
        FREM[icrown,ihost]  = FOLADJ[ihost,icrown,4]*PRBIO[ihost,icrown,4]
    end
end
return (FOLPOT=FOLPOT,FOLADJ=FOLADJ,FOLNH=FOLNH,BWTPHA=BWTPHA,IFHOST=IFHOST,
        HOSTST=HOSTST,FNEW=FNEW,FOLD1=FOLD1,FOLD2=FOLD2,FREM=FREM,PRBIO=PRBIO)
end

const SITA = bwesit_accum(trees)
FOLPOT=SITA.FOLPOT; FOLADJ=SITA.FOLADJ; BWTPHA=SITA.BWTPHA; IFHOST=SITA.IFHOST
HOSTST=SITA.HOSTST; FNEW=SITA.FNEW; FOLD1=SITA.FOLD1; FOLD2=SITA.FOLD2; FREM=SITA.FREM

# --- compare to SIT dump arrays ---
function cmp_sit(sit,FOLPOT,FOLADJ,BWTPHA,FNEW,FOLD1,FOLD2,FREM)
nb=0
for l in sit
    if startswith(l,"FOL ")
        s=split(l); ih=parse(Int,s[2]); ic=parse(Int,s[3]); ia=parse(Int,s[4])
        gp=FOLPOT[ih,ic,ia]; ga=FOLADJ[ih,ic,ia]
        gp===h2f(s[5]) || (nb+=1; nb<=6 && println("  FOLPOT[$ih,$ic,$ia] got=",f2h(gp)," exp=",s[5]))
        ga===h2f(s[6]) || (nb+=1; nb<=6 && println("  FOLADJ[$ih,$ic,$ia] got=",f2h(ga)," exp=",s[6]))
    elseif startswith(l,"BWTPHA")
        s=split(l); ih=parse(Int,s[2]); is=parse(Int,s[3])
        BWTPHA[ih,is]===h2f(s[4]) || (nb+=1; nb<=6 && println("  BWTPHA[$ih,$is] got=",f2h(BWTPHA[ih,is])," exp=",s[4]))
    elseif startswith(l,"FBUF")
        s=split(l); ic=parse(Int,s[2]); ih=parse(Int,s[3])
        FNEW[ic,ih]===h2f(s[4]) || (nb+=1; nb<=6 && println("  FNEW[$ic,$ih] got=",f2h(FNEW[ic,ih])," exp=",s[4]))
        FOLD1[ic,ih]===h2f(s[5]) || (nb+=1)
        FOLD2[ic,ih]===h2f(s[6]) || (nb+=1)
        FREM[ic,ih]===h2f(s[7]) || (nb+=1)
    end
end
return nb
end
let nb=cmp_sit(sit,FOLPOT,FOLADJ,BWTPHA,FNEW,FOLD1,FOLD2,FREM)
    println("BWESIT arrays (FOLPOT/FOLADJ/BWTPHA/FNEW/FOLD1/FOLD2/FREM): mismatches=$nb")
end
println("HOSTST got=",f2h(HOSTST)," exp=",f2h(HOSTST_or), "  ", HOSTST===HOSTST_or ? "OK" : "DIFF")
