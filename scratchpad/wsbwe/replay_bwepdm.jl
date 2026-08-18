# Standalone dump-replay validation of BWEPDM per-tree apply + damage RNG (BWERNP/BWEBET)
# vs instrumented FVSem_wsbwe (host_defol.key, cycle 1). Float32-hex bit-exact.
const DIR = @__DIR__
h2f(s) = reinterpret(Float32, parse(UInt32, s, base=16))
f2h(x::Float32) = uppercase(string(reinterpret(UInt32, x), base=16, pad=8))

# gfortran (glibc) single-precision transcendentals — bit-exact match to the oracle.
const LIBM = "libm.so.6"
@inline flog(x::Float32) = ccall((:logf, LIBM), Float32, (Float32,), x)
@inline fexp(x::Float32) = ccall((:expf, LIBM), Float32, (Float32,), x)
@inline fpow(x::Float32, y::Float32) = ccall((:powf, LIBM), Float32, (Float32,Float32), x, y)
@inline fsqrt(x::Float32) = sqrt(x)   # IEEE-exact

# ---------------- RNG: BWERAN + BWERNP + BWEBET ----------------
mutable struct RNG; s0::Float64; draws::Vector{Float32}; end
RNG(seed) = RNG(Float64(seed), Float32[])
@inline function bweran!(r::RNG)::Float32
    s1 = rem(16807.0 * r.s0, 2147483647.0)
    r.s0 = s1
    v = Float32(s1 / 2147483648.0)
    push!(r.draws, v)
    return v
end

# BWEBET (Cheng 1978). Returns (x, kode).
function bwebet(a0::Float32, b0::Float32, r::RNG)
    if a0 <= 0.00001f0 || b0 <= 0.00001f0
        return (0.0f0, -1)
    end
    kode = 0
    if min(a0,b0) > 1.0f0
        a = min(a0,b0); b = max(a0,b0)
        alf = a + b
        bet = sqrt((alf - 2.0f0)/(2.0f0*a*b - alf))
        gam = a + (1.0f0/b)
        local w
        while true                       # label 10
            kode += 1
            u1 = bweran!(r); u2 = bweran!(r)
            v = bet*flog(u1/(1.0f0-u1))
            if v > 20.0f0; v = 20.0f0; end
            w = a*fexp(v)
            z = u1*u1*u2
            rr = gam*v - 1.3862944f0
            s = a + rr - w
            if s + 2.609438f0 >= 5.0f0*z; break; end   # GOTO 15
            t = flog(z)
            if s >= t; break; end                       # GOTO 15
            if rr + alf*flog(alf/(b+w)) < t; continue; end # GOTO 10
            break                                        # fall to 15
        end
        x = (a == a0) ? w/(b+w) : b/(b+w)
        return (x, kode)
    else
        a = max(a0,b0); b = min(a0,b0)
        alf = a + b
        bet = 1.0f0/b
        del = 1.0f0 + a - b
        c1 = del*(0.0138889f0 + 0.0416667f0*b)/(a*bet - 0.777778f0)
        c2 = 0.25f0 + (0.5f0 + 0.25f0/del)*b
        local w
        while true                       # label 20
            kode += 1
            u1 = bweran!(r); u2 = bweran!(r)
            if u1 >= 0.5f0
                # label 30
                z = u1*u1*u2
                if z <= 0.25f0
                    v = bet*flog(u1/(1.0f0-u1))
                    if v > 20.0f0; v = 20.0f0; end
                    w = a*fexp(v)
                    break                # GOTO 60
                end
                if z >= c2; continue; end # GOTO 20
                # fall to 50
                v = bet*flog(u1/(1.0f0-u1))
                if v > 20.0f0; v = 20.0f0; end
                w = a*fexp(v)
                if alf*(flog(alf/(b+w))+v) - 1.3862944f0 < flog(z); continue; end # GOTO 20
                break                    # GOTO 60
            else
                y = u1*u2
                z = u1*y
                if 0.25f0*u2 + z - y >= c1
                    continue             # GOTO 20
                else
                    # label 50
                    v = bet*flog(u1/(1.0f0-u1))
                    if v > 20.0f0; v = 20.0f0; end
                    w = a*fexp(v)
                    if alf*(flog(alf/(b+w))+v) - 1.3862944f0 < flog(z); continue; end # GOTO 20
                    break                # GOTO 60
                end
            end
        end
        x = (a == a0) ? w/(b+w) : b/(b+w)
        return (x, kode)
    end
end

const BC = 2.7f0
const BSC = 0.9304f0
function bwernp(xmean::Float32, xmxvar::Float32, r::RNG)::Float32
    if xmean < 0.001f0 || xmean > 0.999f0
        return xmean
    end
    xm = xmean > 0.5f0 ? (1.0f0 - xmean) : xmean
    var = (BSC*xmxvar)*(BC*fpow(xm*(BC-1.0f0), BC-1.0f0)*fexp(-fpow((BC-1.0f0)*xm, BC)))
    w = (xmean*(1.0f0-xmean)/var) - 1.0f0
    v = xmean*w
    w = (1.0f0-xmean)*w
    (x, kode) = bwebet(v, w, r)
    return kode < 0 ? xmean : x
end

# gfortran real4 ** int4 (libgfortran pow_r4_i4), right-to-left binary
function pow_r4_i4(base::Float32, n::Int)::Float32
    x = base; e = n; pw = 1.0f0
    while true
        if (e & 1) != 0; pw *= x; end
        e >>= 1
        e != 0 ? (x *= x) : break
    end
    return pw
end

ifix(x::Float32) = trunc(Int, x)   # Fortran IFIX truncates toward zero

# ---------------- coefficients (bwepdm.f) ----------------
const IFIR = Int[0,1,0,0,1,1]
const B0=Float32[46.27900,57.75010,46.27900,46.27900,57.75010,57.75010]
const B1=Float32[-1.80810,-2.28210,-1.80810,-1.80810,-2.28210,-2.28210]
const B2=Float32[0.01930,0.02430,0.01930,0.01930,0.02430,0.02430]
const B3=Float32[0.00550,0.0,0.00550,0.00550,0.0,0.0]
const B4=Float32[-0.00808,-0.00793,-0.00808,-0.00808,-0.00793,-0.00793]
const B5=Float32[0.57450,0.92870,0.57450,0.57450,0.92870,0.92870]
const B6=Float32[-0.24050,-0.22330,-0.24050,-0.24050,-0.22330,-0.22330]
const B7=Float32[-0.09640,-0.13180,-0.09640,-0.09640,-0.13180,-0.13180]
const TK0=Float32[21.89880,17.76920,21.89880,21.89880,17.76920,17.76920]
const TK1=Float32[-0.00714,-0.00535,-0.00714,-0.00714,-0.00535,-0.00535]
const TK2=Float32[6.99e-6,5.22e-6,6.99e-6,6.99e-6,5.22e-6,5.22e-6]
const TK3=Float32[0.0,-0.01270,0.0,0.0,-0.01270,-0.01270]
const TK4=Float32[-5.26520,-7.94630,-5.26520,-5.26520,-7.94630,-7.94630]
const TK5=Float32[-0.13630,-0.20360,-0.13630,-0.13630,-0.20360,-0.20360]
const TK6=Float32[0.28030,0.12690,0.28030,0.28030,0.12690,0.12690]
const TK7=Float32[0.14030,0.38010,0.14030,0.14030,0.38010,0.38010]
const TK8=Float32[-0.03770,-0.04120,-0.03770,-0.03770,-0.04120,-0.04120]

# ---------------- parse dumps ----------------
sit = readlines(joinpath(DIR,"bwe_sit.txt"))
pdm = readlines(joinpath(DIR,"bwe_pdm.txt"))
rng = readlines(joinpath(DIR,"bwe_rng.txt"))

# header
hdr = split(first(filter(l->startswith(l,"PDMHDR"), pdm)))
ICYC=parse(Int,hdr[2]); IBWYR=parse(Int,hdr[3]); NOBWYR=parse(Int,hdr[4]); IPTINV=parse(Int,hdr[5]); ITRN=parse(Int,hdr[6])
FA=h2f(hdr[7]); ELEV=h2f(hdr[8])

PEDDS=zeros(Float32,6,3); PEHTG=zeros(Float32,6,3); AVYRMX=zeros(Float32,6,3)
for l in pdm; startswith(l,"PDMPE") || continue; s=split(l); ih=parse(Int,s[2]); is=parse(Int,s[3]); PEDDS[ih,is]=h2f(s[4]); PEHTG[ih,is]=h2f(s[5]); AVYRMX[ih,is]=h2f(s[6]); end
PNTBA=zeros(Float32,IPTINV); PNTHBA=zeros(Float32,IPTINV)
for l in pdm; startswith(l,"PDMPNT") || continue; s=split(l); i=parse(Int,s[2]); PNTBA[i]=h2f(s[3]); PNTHBA[i]=h2f(s[4]); end
ISCT=zeros(Int,64,2)
for l in pdm; startswith(l,"PDMSCT") || continue; s=split(l); i=parse(Int,s[2]); ISCT[i,1]=parse(Int,s[3]); ISCT[i,2]=parse(Int,s[4]); end
IND1=zeros(Int,ITRN)
for l in pdm; startswith(l,"PDMIND") || continue; s=split(l); IND1[parse(Int,s[2])]=parse(Int,s[3]); end

# per-tree inputs keyed by tree index I
TIN=Dict{Int,Any}()
for l in pdm; startswith(l,"PDMTIN") || continue; s=split(l)
    I=parse(Int,s[2])
    TIN[I]=(ISPI=parse(Int,s[3]),IHOST=parse(Int,s[4]),ISZI=parse(Int,s[5]),ICRC1=parse(Int,s[6]),
            IIT0=parse(Int,s[7]),INH0=parse(Int,s[8]),IMC0=parse(Int,s[9]),ICR0=parse(Int,s[10]),
            ITRE=parse(Int,s[11]),IMIST=parse(Int,s[12]),
            HT=h2f(s[13]),DBH=h2f(s[14]),DGI=h2f(s[15]),HTGI=h2f(s[16]),WK20=h2f(s[17]),
            PROB=h2f(s[18]),CFV=h2f(s[19]),BARK=h2f(s[20]),MFT=h2f(s[21]),MFM=h2f(s[22]),PCTK=h2f(s[23]))
end
TOU=Dict{Int,Any}()
for l in pdm; startswith(l,"PDMTOU") || continue; s=split(l)
    I=parse(Int,s[2])
    TOU[I]=(ITRUNC=parse(Int,s[3]),NORMHT=parse(Int,s[4]),IMC=parse(Int,s[5]),ICR=parse(Int,s[6]),
            CAVDEF=h2f(s[7]),CX=h2f(s[8]),CPART=h2f(s[9]),CXD=h2f(s[10]),CXH=h2f(s[11]),
            KTK=h2f(s[12]),PR=h2f(s[13]),DG=h2f(s[14]),HTG=h2f(s[15]),HT=h2f(s[16]),WK2=h2f(s[17]))
end

# ---------------- replay BWEPDM per-tree, in ISCT/IND1 order, drawing RNG ----------------
function run_replay()
r = RNG(55329.0)   # BWERPT(DSEEDD=55329) at BWEPDM entry
mism = 0; ntree = 0
LTOPK = true
MAXSP = 19
chk(name,I,got::Float32,exp::Float32) = (got===exp ? true : (println("  MISMATCH tree $I $name got=",f2h(got)," exp=",f2h(exp)); false))
chki(name,I,got::Int,exp::Int) = (got==exp ? true : (println("  MISMATCH tree $I $name got=",got," exp=",exp); false))

for ISPI in 1:MAXSP
    ISCT[ISPI,1]==0 && continue
    # IHOST from IBWSPM handled via dumped per-tree IHOST; here mimic: skip nonhost/larch
    for II in ISCT[ISPI,1]:ISCT[ISPI,2]
        I = IND1[II]
        t = TIN[I]; o = TOU[I]
        IHOST=t.IHOST; ISZI=t.ISZI; ICRC1=t.ICRC1
        IHOST>=6 && continue
        ntree += 1
        H = t.HT; BARK=t.BARK; DBH=t.DBH; DGI=t.DGI; HTGI=t.HTGI
        MFT=t.MFT; MFM=t.MFM; PCTK=t.PCTK
        ITRUNC=t.IIT0; NORMHT=t.INH0; IMC=t.IMC0; ICR=t.ICR0
        HTG = HTGI          # current htg
        cavdef=-9.0f0; cx=-9.0f0; cpart=-9.0f0; cxd=-9.0f0; cxh=-9.0f0
        FTKILL=0.0f0; PRTOPK=0.0f0
        if LTOPK
            AVDEF = bwernp(AVYRMX[IHOST,ISZI], 0.06f0, r); cavdef=AVDEF
            if ISZI==1
                X = H; if H>10.0f0; X=10.0f0; end
                PRTOPK = 1.0f0/(1.0f0+fexp(-(-2.5817f0 -0.027635f0*Float32(ICR) + 3.709f0*sqrt(AVDEF) + 0.0488f0*X)))
            else
                if IFIR[IHOST]==0
                    PRTOPK = 0.96f0*(1.0f0-fexp(-(fpow(0.65f0*(AVDEF+1.0f0), 14.0f0))))
                else
                    PRTOPK = 0.90f0*(1.0f0-fexp(-(fpow(0.60f0*(AVDEF+1.0f0), 11.0f0))))
                end
            end
            X = bweran!(r); cx=X
            if X < PRTOPK
                if MFT>0.0f0 && MFM>0.0f0
                    PART = 1.0f0/(1.0f0+fexp(TK0[IHOST]+TK1[IHOST]*ELEV+TK2[IHOST]*ELEV*ELEV+
                        TK3[IHOST]*Float32(t.IMIST)+TK4[IHOST]*PCTK+TK5[IHOST]*MFT+TK6[IHOST]*MFM+
                        TK7[IHOST]*MFM*PCTK+TK8[IHOST]*MFT*MFM))
                else
                    PART = 0.0f0
                end
                cpart=PART
                if PART>0.9f0; PART=0.9f0; end
                if PART>0.0f0
                    FTKILL = H*PART
                    TOPH = H - FTKILL
                    ITRC2 = ifix(TOPH*100.0f0+0.5f0)
                    ITRC1 = ITRUNC
                    if ITRC1>0
                        if ITRC1>ITRC2; ITRUNC=ITRC2; end
                        H = TOPH
                    else
                        D = DBH*BARK
                        if H>=25.0f0 && D>=6.0f0
                            AF = t.CFV/(0.00545415f0*D*D*H)
                            AF = 0.44244f0 - (0.99167f0/AF) - (1.43237f0*flog(AF)) + (1.68581f0*sqrt(AF)) - (0.13611f0*AF*AF)
                            DTK = FTKILL/H
                            DTK = (DTK/((AF*DTK)+(1.0f0-AF)))*D
                            if DTK>4.0f0
                                ITRUNC=ITRC2; NORMHT=ifix(H*100.0f0+0.5f0); IMC=3
                            elseif DTK>2.0f0 && IMC<2
                                IMC=2
                            end
                        end
                        H = TOPH
                        IOD = ICR
                        if IOD>=0
                            CN = (Float32(IOD)/100.0f0*H) - H + TOPH
                            NEW = ifix(CN/TOPH*100.0f0+0.5f0)
                            if NEW<5; NEW=5; end
                            ICR = -NEW
                        end
                    end
                end
                HTG = 0.0f0
            end
        end
        # DG modification
        DBHI=DBH
        DDS = DGI*(2.0f0*BARK*DBHI+DGI)
        XD = bwernp(PEDDS[IHOST,ISZI], 0.03f0, r); cxd=XD
        DDS = DDS*XD
        DG = sqrt((DBHI*BARK)^2 + DDS) - BARK*DBHI
        if DG<0.0f0; DG=0.0f0; end
        # HTG modification
        if HTG>0.0f0
            if PEDDS[IHOST,ISZI]<0.99f0
                XH = PEHTG[IHOST,ISZI]+((1.0f0-PEHTG[IHOST,ISZI])/(1.0f0-PEDDS[IHOST,ISZI])*(XD-PEDDS[IHOST,ISZI]))
            else
                XH = bwernp(PEHTG[IHOST,ISZI], 0.03f0, r)
            end
            cxh=XH
            HTG = HTG*XH
            if HTG<0.0f0; HTG=0.0f0; end
        end
        # mortality
        if ITRUNC==0 || NORMHT==0
            KTK = 0.0f0
        else
            KTK = 10.0f0 - Float32(ifix(Float32(ITRUNC)/Float32(NORMHT)*10.0f0 + 0.5f0))
        end
        if MFT>0.0f0 && MFM>0.0f0
            PR = 1.0f0/(1.0f0+fexp(B0[IHOST]+B1[IHOST]*ELEV+B2[IHOST]*ELEV*ELEV+
                B3[IHOST]*PNTBA[t.ITRE]+B4[IHOST]*PNTHBA[t.ITRE]+B5[IHOST]*MFT+
                B6[IHOST]*KTK+B7[IHOST]*MFT*MFM))
        else
            PR = 0.0f0
        end
        BASE = t.WK20/t.PROB
        PR = 1.0f0 - PR
        if NOBWYR>0
            PR = pow_r4_i4(PR,IBWYR)*((1.0f0-BASE)^FA)
        else
            PR = pow_r4_i4(PR,IBWYR)
        end
        PR = 1.0f0 - PR
        WK2 = t.WK20
        if BASE>PR
            PR = BASE
        else
            if PR>0.98f0; PR=0.98f0; end
            PR = t.PROB*PR
            WK2 = PR
        end
        # compare
        ok = true
        ok &= chk("CAVDEF",I,cavdef,o.CAVDEF)
        ok &= chk("CX",I,cx,o.CX)
        ok &= chk("CPART",I,cpart,o.CPART)
        ok &= chk("CXD",I,cxd,o.CXD)
        ok &= chk("CXH",I,cxh,o.CXH)
        ok &= chk("KTK",I,KTK,o.KTK)
        ok &= chk("PR",I,PR,o.PR)
        ok &= chk("DG",I,DG,o.DG)
        ok &= chk("HTG",I,HTG,o.HTG)
        ok &= chk("HT",I,H,o.HT)
        ok &= chk("WK2",I,WK2,o.WK2)
        ok &= chki("ITRUNC",I,ITRUNC,o.ITRUNC)
        ok &= chki("NORMHT",I,NORMHT,o.NORMHT)
        ok &= chki("IMC",I,IMC,o.IMC)
        ok &= chki("ICR",I,ICR,o.ICR)
        ok || (mism += 1)
    end
end

# compare RNG stream
exp_rng = [h2f(split(l)[2]) for l in rng if startswith(l,"R")]
rngok = length(r.draws)==length(exp_rng) && all(r.draws .=== exp_rng)
return (ntree=ntree, mism=mism, draws=r.draws, exp_rng=exp_rng, rngok=rngok)
end

res = run_replay()
ntree=res.ntree; mism=res.mism; rngok=res.rngok; exp_rng=res.exp_rng
r = (draws=res.draws,)
println("\n=== RESULTS ===")
println("trees processed: $ntree  mismatched trees: $mism")
println("RNG draws: julia=$(length(r.draws)) oracle=$(length(exp_rng))  stream bit-exact=$rngok")
if !rngok
    n = min(length(r.draws),length(exp_rng))
    for i in 1:n
        if r.draws[i] !== exp_rng[i]
            println("  first RNG diff at draw $i: jl=",f2h(r.draws[i])," or=",f2h(exp_rng[i])); break
        end
    end
end
println(mism==0 && rngok ? "BWEPDM REPLAY: PASS (bit-exact)" : "BWEPDM REPLAY: FAIL")
