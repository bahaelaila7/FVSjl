# =============================================================================
# WWPB multi-stand landscape dispersal (bmatct.f / bmdrv.f, MXSTND>1).
#
# GOLDEN: scratchpad/wwpb/driver_bmatct_multi.f relinks the PRISTINE wwpb/bmatct.f
# (ppework wwpb/src) over a 4-stand landscape with three Outside-World regimes
# (A: OW-off pure interstand redistribution; B: OW floating UFLOAT=-1; C: OW
# fixed constants). The Float32 hex below are the gfortran-16 outputs (verbatim);
# bmatct_multi! must reproduce every one BIT-EXACT. Oracle-free at test time.
# Reproduce the golden: cd scratchpad/wwpb && gfortran-16 -std=legacy -w
#   -fno-automatic -Iinc_multi driver_bmatct_multi.f
#   /workspace/.ppework/wwpb/src/bmatct.f -o drv_bmatct_multi && ./drv_bmatct_multi
# =============================================================================
using Test
using FVSjl
const F = FVSjl

hx(x::Float32) = uppercase(string(reinterpret(UInt32, x), base=16, pad=8))

# 4-stand landscape geometry (meters, meters, acres) — matches the driver.
const XLOC = Float32[0, 1000, 0, 1200]
const YLOC = Float32[0, 0, 1500, 1200]
const AREA = Float32[100, 200, 150, 80]
const STOCK = [true, true, true, true]
const NUMER1 = Float32[1000, 300, 700, 500]
const BKP0   = Float32[100, 40, 10, 60]
const TFOOD1 = Float32[50, 80, 30, 20]

function mkstands()
    ss = F.WwpbStand[]
    for i in 1:4
        s = F.WwpbStand()
        s.numer[1] = NUMER1[i]; s.numer[2] = 0.0f0
        s.bkp = BKP0[i]; s.bkpips = 0.0f0
        s.tfood[1] = TFOOD1[i]; s.tfood[2] = 0.0f0
        s.grfstd = 1.0f0
        push!(ss, s)
    end
    ss
end

# (tag => stand => (BKP, ATTC, ATB, BOU, BIN, SLF, BPS)) gfortran-16 golden hex.
const GOLD = Dict(
 ("A",1)=>("41D61A7C","00000000","3D65FEC9","00000000","00000000","42B9FA2A","3E828A22"),
 ("A",2)=>("4141D19E","00000000","3D99F86D","00000000","00000000","4202DF90","3EAC5820"),
 ("A",3)=>("409AC9F9","00000000","3D896410","00000000","00000000","4116F0CA","3EAB5F61"),
 ("A",4)=>("414BBA72","00000000","3D4E8E94","00000000","00000000","424CBBA8","3E6B0ACB"),
 ("B",1)=>("4193487E","3DA05180","3D65FEC9","422AE832","408E3A49","425507C8","3E95D9CA"),
 ("B",2)=>("40BC48BE","3D827BD7","3D99F86D","41B6A4A3","3F8AE95F","4160B44F","3EB3EB90"),
 ("B",3)=>("4062815F","3D97A452","3D896410","4093960B","403C57A6","40A2A6CE","3EAFA04D"),
 ("B",4)=>("40DF8B35","3DAEA8ED","3D4E8E94","4215C9C5","401AF36B","4199E945","3E97E996"),
 ("C",1)=>("40BDB6FE","3DA05180","3D65FEC9","42A5E0C9","3E22B7D1","417DD606","3EAF9156"),
 ("C",2)=>("3FB32EF0","3D827BD7","3D99F86D","420F71F4","3D1EEC9A","4058A8B1","3EB93D76"),
 ("C",3)=>("3F457C90","3D97A452","3D896410","4107AD35","3DD779F1","3FB791C8","3EB8777C"),
 ("C",4)=>("3FDA203C","3DAEA8ED","3D4E8E94","425BB1F6","3DB14632","408A91A7","3EB2CCD9"),
)

@testset "WWPB bmatct_multi! (multi-stand landscape dispersal)" begin
    w = F.wwpb_defaults!(:ie)  # PBSPEC=1 (MPB)
    usera = (1.0f0,1.0f0,1.0f0); selfa = (1.0f0,1.0f0,1.0f0)
    userc = (100.0f0,100.0f0,100.0f0); urmax = (1.0f0,1.0f0,1.0f0)

    checkone(tag, ls, ss) = begin
        for i in 1:4
            g = GOLD[(tag,i)]
            @test hx(ss[i].bkp)        == g[1]
            @test hx(ls.attc[1,i])     == g[2]
            @test hx(ls.attbyk[1,i])   == g[3]
            @test hx(ls.bkpout[1,i])   == g[4]
            @test hx(ls.bkpin[1,i])    == g[5]
            @test hx(ls.selfbkp[1,i])  == g[6]
            @test hx(ls.bkps[i])       == g[7]
        end
    end

    # A: Outside World OFF — pure interstand SCORE/PROP redistribution.
    lsA = F.WwpbLandscape(XLOC,YLOC,AREA,STOCK); ssA = mkstands()
    F.bmatct_multi!(lsA, ssA, w; usera=usera,selfa=selfa,userc=userc,urmax=urmax,
                    outoff=true, ipson=false)
    checkone("A", lsA, ssA)

    # B: Outside World ON, floating (UFLOAT=-1).
    lsB = F.WwpbLandscape(XLOC,YLOC,AREA,STOCK); ssB = mkstands()
    F.bmatct_multi!(lsB, ssB, w; usera=usera,selfa=selfa,userc=userc,urmax=urmax,
                    outoff=false, ufloat=-1.0f0, rvod=1.0f0, stocko=1.0f0, ipson=false)
    checkone("B", lsB, ssB)

    # C: Outside World ON, fixed constants (UFLOAT=0).
    lsC = F.WwpbLandscape(XLOC,YLOC,AREA,STOCK; cbao=50.0f0, crvond=1.0f0,
                          cbaho=Float32[20,-1], cbaspo=Float32[5,-1],
                          cspo=Float32[2,-1], cbkpo=Float32[10,-1]); ssC = mkstands()
    F.bmatct_multi!(lsC, ssC, w; usera=usera,selfa=selfa,userc=userc,urmax=urmax,
                    outoff=false, ufloat=0.0f0, rvod=1.0f0, stocko=1.0f0, ipson=false)
    checkone("C", lsC, ssC)
end

# Multi-year cascade golden (BKP per stand per year): reuse the SAME landscape so
# ACDONE persists (spatial cache computed once, reused), NUMER held fixed, BKP evolves
# through repeated BMATCT calls — the master-cycle dispersal cascade bmdrv relies on.
const GOLD_D = [  # D[year][stand]
 ["41D61A7C","4141D19E","409AC9F9","414BBA72"],
 ["4116CBBF","4074306C","4001902A","40800C19"],
 ["40607F3E","3F9D4BA8","3F540BD4","3FAE5632"],
 ["3FA8A1A1","3ECDBA40","3EAB1FBC","3EF2A7D8"],
 ["3EFC8488","3E083EC0","3E086C38","3E29EE80"]]
const GOLD_E = [
 ["4193487E","40BC48BE","4062815F","40DF8B35"],
 ["407F1DB8","3F697E60","3F730A5C","3F87119E"],
 ["3F5F3BDA","3E18FD00","3E70DF20","3E3186E0"],
 ["3E41F840","3CD48A00","3D64D2C0","3CF91780"],
 ["3D277B00","3B9C1800","3C53D100","3BB95000"]]

@testset "WWPB bmatct_multi! multi-year cascade (ACDONE persistence)" begin
    w = F.wwpb_defaults!(:ie)
    usera=(1.0f0,1.0f0,1.0f0); selfa=(1.0f0,1.0f0,1.0f0)
    userc=(100.0f0,100.0f0,100.0f0); urmax=(1.0f0,1.0f0,1.0f0)
    # D: OW-off cascade
    lsD = F.WwpbLandscape(XLOC,YLOC,AREA,STOCK); ssD = mkstands()
    for yr in 1:5
        F.bmatct_multi!(lsD, ssD, w; usera=usera,selfa=selfa,userc=userc,urmax=urmax,
                        outoff=true, ipson=false)
        for i in 1:4; @test hx(ssD[i].bkp) == GOLD_D[yr][i]; end
    end
    # E: OW-floating cascade
    lsE = F.WwpbLandscape(XLOC,YLOC,AREA,STOCK); ssE = mkstands()
    for yr in 1:5
        F.bmatct_multi!(lsE, ssE, w; usera=usera,selfa=selfa,userc=userc,urmax=urmax,
                        outoff=false, ufloat=-1.0f0, rvod=1.0f0, stocko=1.0f0, ipson=false)
        for i in 1:4; @test hx(ssE[i].bkp) == GOLD_E[yr][i]; end
    end
end

@testset "WWPB bmdrv_multi! collapses to single-stand at MXSTND=1" begin
    # At MXSTND=1 / OUTOFF=T the multi-stand BMATCT self-cancels (PROP=1), so the
    # bmdrv_multi! year loop must reproduce wwpb_outbreak_cycle! BIT-EXACT.
    w = F.wwpb_defaults!(:ie)
    coeffs = F.wwpb_init_coeffs(w.upsiz)

    # a small synthetic treelist (one host-pine record + one nonhost).
    mktrees() = (n=2,
        species=Int32[7, 3],                # LP host (sp7), DF nonhost
        dbh=Float32[12.0, 10.0],
        tpa=Float32[40.0, 30.0],
        height=Float32[60.0, 55.0],
        crown_pct=Float32[40.0, 45.0],
        ht_growth=Float32[1.0, 1.0],
        cuft_vol=Float32[15.0, 12.0])
    code = ["WP","L","DF","GF","WH","C","LP","S","AF","PP","OTH"]
    spα(sp::Int) = (1 <= sp <= length(code)) ? code[sp] : ""

    sarea = 40.0f0; seed = zeros(Float32, F.WWPB_NSCL); seed[3] = 5.0f0

    # reference: single-stand outbreak
    stRef = F.WwpbStand()
    F.wwpb_outbreak_cycle!(stRef, w, coeffs, mktrees(), spα;
                           sarea=sarea, iyr1=1, iyr2=3, seed_pbkill=copy(seed))

    # multi-stand (single stand) via bmdrv_multi!
    w2 = F.wwpb_defaults!(:ie)
    st2 = F.WwpbStand()
    F.bmsdit!(st2, mktrees(), w2, spα)
    fill!(st2.rvdsc, 1.0f0); st2.pbkill .= seed
    ls = F.WwpbLandscape(Float32[0.0], Float32[0.0], Float32[sarea], [true])
    F.bmdrv_multi!(ls, [st2], w2, coeffs; area=Float32[sarea], iyr1=1, iyr2=3,
                   outoff=true)

    @test hx(st2.bkp) == hx(stRef.bkp)
    @test all(hx.(st2.pbkill) .== hx.(stRef.pbkill))
    @test all(hx.(vec(st2.tree)) .== hx.(vec(stRef.tree)))
    @test all(hx.(vec(st2.tpbk)) .== hx.(vec(stRef.tpbk)))
end

@testset "PPE mode-2 landscape composition (ppe_landscape_mode2!)" begin
    # The chunk-4 wiring: grow-to-mortality per-stand states → landscape dispersal →
    # per-stand BMKILL handback. Verify it composes on a placed 2-stand landscape:
    # BKP is redistributed across stands (interstand coupling that mode-1 cannot do),
    # and BMKILL produces a per-record mortality.
    w = F.wwpb_defaults!(:ie)
    coeffs = F.wwpb_init_coeffs(w.upsiz)
    code = ["WP","L","DF","GF","WH","C","LP","S","AF","PP","OTH"]
    spα(sp::Int) = (1 <= sp <= length(code)) ? code[sp] : ""
    mktrees(tpa) = (n=2, species=Int32[7,3], dbh=Float32[12.0,10.0], tpa=Float32[tpa,25.0],
                    height=Float32[60.0,55.0], crown_pct=Float32[40.0,45.0],
                    ht_growth=Float32[1.0,1.0], cuft_vol=Float32[15.0,12.0])
    seedA = zeros(Float32, F.WWPB_NSCL); seedA[3] = 6.0f0

    members = [
        (trees=mktrees(50.0f0), area=100.0f0, xloc=0.0f0,    yloc=0.0f0, stock=true, sp_alpha=spα),
        (trees=mktrees(40.0f0), area=120.0f0, xloc=800.0f0,  yloc=0.0f0, stock=true, sp_alpha=spα),
    ]
    kills, ls, stands = F.ppe_landscape_mode2!(members, w, coeffs; iyr1=1, iyr2=3,
                        usera=(1.0f0,1.0f0,1.0f0), selfa=(1.0f0,1.0f0,1.0f0),
                        userc=(100.0f0,100.0f0,100.0f0), urmax=(1.0f0,1.0f0,1.0f0),
                        outoff=true, seed_pbkill=[copy(seedA), nothing])

    @test length(kills) == 2
    @test all(k -> length(k) == 2, kills)          # per-record WK2 for each stand
    @test all(k -> all(isfinite, k), kills)
    # interstand redistribution occurred: within URMAX the seeded stand's SELFBKP is a
    # PROPER fraction of its total BKP (some pressure dispersed to the neighbor), so it
    # is strictly less than the no-neighbor self-allocation (= area·BKP).
    @test ls.n == 2
    @test ls.selfbkp[1,1] > 0.0f0
    @test ls.selfbkp[1,1] < ls.area[1] * 200.0f0   # loose upper bound; coupling present
    # some mortality was produced on the seeded host stand
    @test sum(kills[1]) > 0.0f0
end
