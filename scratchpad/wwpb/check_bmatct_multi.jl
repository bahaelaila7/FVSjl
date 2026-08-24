# Validate bmatct_multi! bit-exact vs driver_bmatct_multi.f golden.
push!(LOAD_PATH, "/workspace/.wt-ppe2/src")
using FVSjl
const F = FVSjl

hx(x::Float32) = uppercase(string(reinterpret(UInt32, x), base=16, pad=8))

# landscape geometry (matches driver)
xloc = Float32[0, 1000, 0, 1200]
yloc = Float32[0, 0, 1500, 1200]
area = Float32[100, 200, 150, 80]
stock = [true, true, true, true]

numer1 = Float32[1000, 300, 700, 500]
bkp0   = Float32[100, 40, 10, 60]
tfood1 = Float32[50, 80, 30, 20]

function mkstands()
    ss = F.WwpbStand[]
    for i in 1:4
        s = F.WwpbStand()
        s.numer[1] = numer1[i]; s.numer[2] = 0.0f0
        s.bkp = bkp0[i]; s.bkpips = 0.0f0
        s.tfood[1] = tfood1[i]; s.tfood[2] = 0.0f0
        s.grfstd = 1.0f0
        push!(ss, s)
    end
    ss
end

w = F.wwpb_defaults!(:ie)   # pbspec=1
usera = (1.0f0,1.0f0,1.0f0); selfa=(1.0f0,1.0f0,1.0f0)
userc = (100.0f0,100.0f0,100.0f0); urmax=(1.0f0,1.0f0,1.0f0)

results = Dict{String,Vector{String}}()

# Scenario A: OUTOFF
ls = F.WwpbLandscape(xloc,yloc,area,stock)
ss = mkstands()
F.bmatct_multi!(ls, ss, w; usera=usera,selfa=selfa,userc=userc,urmax=urmax,
                outoff=true, ipson=false)
lines = String[]
for i in 1:4
    push!(lines, "A $i BKP $(hx(ss[i].bkp)) ATTC $(hx(ls.attc[1,i]))")
    push!(lines, "A $i ATB $(hx(ls.attbyk[1,i])) BOU $(hx(ls.bkpout[1,i]))")
    push!(lines, "A $i BIN $(hx(ls.bkpin[1,i])) SLF $(hx(ls.selfbkp[1,i]))")
    push!(lines, "A $i BPS $(hx(ls.bkps[i]))")
end

# Scenario B: OW floating (ufloat=-1); fresh landscape (acdone reset)
lsB = F.WwpbLandscape(xloc,yloc,area,stock)
ssB = mkstands()
F.bmatct_multi!(lsB, ssB, w; usera=usera,selfa=selfa,userc=userc,urmax=urmax,
                outoff=false, ufloat=-1.0f0, rvod=1.0f0, stocko=1.0f0, ipson=false)
for i in 1:4
    push!(lines, "B $i BKP $(hx(ssB[i].bkp)) ATTC $(hx(lsB.attc[1,i]))")
    push!(lines, "B $i ATB $(hx(lsB.attbyk[1,i])) BOU $(hx(lsB.bkpout[1,i]))")
    push!(lines, "B $i BIN $(hx(lsB.bkpin[1,i])) SLF $(hx(lsB.selfbkp[1,i]))")
    push!(lines, "B $i BPS $(hx(lsB.bkps[i]))")
end

# Scenario C: OW fixed constants (ufloat=0)
lsC = F.WwpbLandscape(xloc,yloc,area,stock; cbao=50.0f0, crvond=1.0f0,
                      cbaho=Float32[20,-1], cbaspo=Float32[5,-1], cspo=Float32[2,-1],
                      cbkpo=Float32[10,-1])
ssC = mkstands()
F.bmatct_multi!(lsC, ssC, w; usera=usera,selfa=selfa,userc=userc,urmax=urmax,
                outoff=false, ufloat=0.0f0, rvod=1.0f0, stocko=1.0f0, ipson=false)
for i in 1:4
    push!(lines, "C $i BKP $(hx(ssC[i].bkp)) ATTC $(hx(lsC.attc[1,i]))")
    push!(lines, "C $i ATB $(hx(lsC.attbyk[1,i])) BOU $(hx(lsC.bkpout[1,i]))")
    push!(lines, "C $i BIN $(hx(lsC.bkpin[1,i])) SLF $(hx(lsC.selfbkp[1,i]))")
    push!(lines, "C $i BPS $(hx(lsC.bkps[i]))")
end

for l in lines
    println(l)
end
