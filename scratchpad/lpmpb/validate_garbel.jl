# GARBEL dump-replay: reduce 18 LP trees -> NACLAS=10 classes (lpmpb/garbel.f).
# KEYMPB=2,3,6*0,1 -> attrs A1=DBH (wt 1.0), A2=XPT/phloem (wt 4.0); IMP=1.
# Validates membership + class-PROB vs the oracle fort.780/781 dumps.

h2f(s) = reinterpret(Float32, parse(UInt32, s, base=16))

# --- parse per-tree inputs (fort.781: IN  I  DBH XPT PROB HT) ---
DBH = zeros(Float32, 18); XPT = zeros(Float32, 18); PROB = zeros(Float32, 18)
for ln in eachline("/workspace/.iework/lpmpb/run/fort.781")
    p = split(ln)
    p[1] == "IN" || continue
    i = parse(Int, p[2])
    DBH[i] = h2f(p[3]); XPT[i] = h2f(p[4]); PROB[i] = h2f(p[6])  # fields: DBH XPT WK3 PROB
end

NRECS = 18; NCLAS = 10; PN1 = 0.5f0
IPT = collect(1:NRECS)                     # packed pointer, loaded 1..N
Pw  = zeros(Float32, NRECS)

# --- GRPSUM (grpsum.f): standardized weighted score accumulation ---
function grpsum!(P, IPT, ATR, WT)
    NR = length(IPT); XN = Float32(NR)
    AVE = 0.0f0; STDV = 0.0f0
    for ii in 1:NR
        i = IPT[ii]; AVE += ATR[i]; STDV += ATR[i]^2
    end
    STDV = (STDV - AVE*AVE/XN)/XN
    STDV = STDV > 1.0f-9 ? sqrt(STDV) : 1.0f0
    AVE = AVE/XN
    for ii in 1:NR
        i = IPT[ii]; P[i] += WT*(ATR[i]-AVE)/STDV
    end
end
grpsum!(Pw, IPT, DBH, 1.0f0)               # BETTER(1)=1.0
grpsum!(Pw, IPT, XPT, 4.0f0)               # BETTER(2)=4.0

# --- RDPSRT (.FALSE.): sort IPT so P(IPT[1]) >= P(IPT[2]) >= ... (descending) ---
sort!(IPT, by = i -> -Pw[i])

# --- STEP5 method 1: NCL1 max-diff class boundaries ---
NCL1 = trunc(Int, NCLAS*PN1 + 0.5f0); NCL1 < 1 && (NCL1 = 1)   # Fortran REAL->INT truncates: 5.5->5
NCL2 = NCLAS - NCL1                                            # =5
if !(NRECS - NCLAS > 5); NCL1 = NCLAS; NCL2 = 0; end           # 8>5 -> method2 on
CLASd = zeros(Float32, NCL1)               # CLAS(J,IMP) diff storage
ISC1  = zeros(Int, max(NCLAS, NCL1) + NCL2 + 4)  # ISC(J,1)
ISC2  = zeros(Int, length(ISC1))                 # ISC(J,2)
NC1 = NCL1 - 1
ipt1 = IPT[1]
for I in 2:NRECS
    ipt2 = IPT[I]
    diffp = Pw[ipt1] - Pw[ipt2]
    Jset = NCL1
    for J in 1:NC1
        if diffp <= CLASd[J+1]; Jset = J; break; end
        CLASd[J] = CLASd[J+1]; ISC1[J] = ISC1[J+1]
    end
    CLASd[Jset] = diffp; ISC1[Jset] = I
    global ipt1 = ipt2
end
ISC1[1] = NRECS + 1                        # sentinel
sort!(view(ISC1, 1:NCL1))                  # IQRSRT ascending

# STEP6(A): class lengths, ISC(J,1)=end of class J
I1 = 1
for J in 1:NCL1
    I2 = ISC1[J] - 1
    ISC1[J] = I2; ISC2[J] = I2 - I1 + 1; global I1 = I2 + 1
end
println("jl M1CLASS ends=", [ISC1[J] for J in 1:NCL1], " (oracle [1,6,16,17,18])")

# STEP6(C): method 2 -- split the NCL2 largest classes
for K in 1:NCL2
    MAX = 0; Ibig = 0
    for J in 1:NCL1
        if ISC2[J] > MAX; MAX = ISC2[J]; Ibig = J; end
    end
    global NCL1 += 1
    MAXs = floor(Int, MAX/2.0f0 + 0.5f0)
    ISC1[NCL1] = ISC1[Ibig]
    ISC1[Ibig] = ISC1[Ibig] - MAXs
    ISC2[NCL1] = MAXs
    ISC2[Ibig] = ISC2[Ibig] - MAXs
end
sort!(view(ISC1, 1:NCL1))                  # IQRSRT resort

# STEP6(E): assign sector pointers MPISC(J,1)=start, MPISC(J,2)=end
I1 = 1
MP1 = zeros(Int, NCL1); MP2 = zeros(Int, NCL1)
for J in 1:NCL1
    I2 = ISC1[J]; MP1[J] = I1; MP2[J] = I2; global I1 = I2 + 1
end

# STEP7: class PROB = sum of PROB over members
CLSPROB = zeros(Float32, NCL1)
for I in 1:NCL1, jj in MP1[I]:MP2[I]
    CLSPROB[I] += PROB[IPT[jj]]
end

# --- compare vs oracle fort.780 dump (CLS i mp1 mp2 hex ; MEM i ipt) ---
orc_mp = Dict{Int,Tuple{Int,Int}}(); orc_hex = Dict{Int,String}(); orc_mem = Dict{Int,Vector{Int}}()
for ln in eachline("/workspace/.iework/lpmpb/run/fort.780")
    p = split(ln)
    if p[1] == "CLS"
        c = parse(Int,p[2]); orc_mp[c]=(parse(Int,p[3]),parse(Int,p[4])); orc_hex[c]=p[5]
    elseif p[1] == "MEM"
        c = parse(Int,p[2]); push!(get!(orc_mem,c,Int[]), parse(Int,p[3]))
    end
end

# --- GRCLAS class-avg DBH (CLASS,2) + SURFLP surface (surfce.f/surflp.f) ---
surflp(d::Float32) = d > 5.0f0 ? 8.835f0*d - 40.82f0 : d*0.672f0
avgDBH = zeros(Float32, NCL1); SURF = zeros(Float32, NCL1)
for c in 1:NCL1
    acc = 0.0f0
    for jj in MP1[c]:MP2[c]; i = IPT[jj]; acc += DBH[i]*PROB[i]; end
    avgDBH[c] = CLSPROB[c] > 1.0f-30 ? acc/CLSPROB[c] : 0.0f0
    SURF[c] = surflp(avgDBH[c])
end
function compare_surf()
    ok = true
    orc = Dict{Int,Tuple{String,String}}()
    for ln in eachline("/workspace/.iework/lpmpb/run/fort.783")
        p = split(ln); p[1]=="SURF" || continue
        orc[parse(Int,p[2])] = (uppercase(p[3]), uppercase(p[4]))
    end
    println("\n--- SURFCE (GRCLAS avg-DBH + SURFLP) ---")
    for c in 1:NCL1
        jd = uppercase(string(reinterpret(UInt32, avgDBH[c]), base=16, pad=8))
        js = uppercase(string(reinterpret(UInt32, SURF[c]),   base=16, pad=8))
        od, os = get(orc, c, ("",""))
        m = (jd==od && js==os); ok &= m
        println("SURF $c avgDBH jl=$jd orc=$od  surf jl=$js orc=$os  $(m ? "ok" : "XX")")
    end
    println(ok ? "*** SURFCE BIT-EXACT (avg-DBH + LP surface) ***" : "!!! SURFCE MISMATCH")
end

function compare()
    println("NACLAS jl=$NCL1 oracle=$NCLAS")
    ok = (NCL1 == NCLAS); memok = true; hexok = true
    for c in 1:NCL1
        jlmem = [IPT[jj] for jj in MP1[c]:MP2[c]]
        omem  = get(orc_mem, c, Int[])
        jlhex = uppercase(string(reinterpret(UInt32, CLSPROB[c]), base=16, pad=8))
        ohex  = uppercase(get(orc_hex, c, ""))
        m = (jlmem == omem); h = (jlhex == ohex)
        memok &= m; hexok &= h
        println("CLS $c  jl[$(MP1[c]),$(MP2[c])] mem=$jlmem  orc=$omem  $(m ? "ok" : "XX")  prob jl=$jlhex orc=$ohex $(h ? "ok" : "XX")")
    end
    println(ok && memok && hexok ? "\n*** GARBEL BIT-EXACT: NACLAS + membership + class-PROB all match ***" :
            "\n!!! MISMATCH  count=$ok mem=$memok prob=$hexok")
end
compare()
compare_surf()
