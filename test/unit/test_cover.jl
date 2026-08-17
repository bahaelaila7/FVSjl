# test_cover.jl — COVER canopy-cover report path, dump-replay bit-exact vs FVSem_g16.
#
# COVER (covr/vcovr) is a REPORT-ONLY understory/canopy-cover extension (no tree-record
# write, no RNG), active in the shipped western binaries. These tests replay the four
# canopy kernels against per-tree / per-height-class dumps captured from an instrumented
# FVSem_g16 (single-.o swap of cvcw/cvshap/cvcbms/cvsum.f; the instrumented COVER report
# was proven byte-identical to clean), EM emt01, 2 cycles → 3 CVCNOP invocations.
#
# Golden dumps (test/unit/):
#   cover_cvcw_g16_dump.txt   CVCWTREE icyc tree isp TRECW P CRAREA ; CVCWFINAL icyc itrn CRAREA
#   cover_cvshap_g16_dump.txt CVSHAP  icyc tree isp ispi DBH HT CR CL RAD TPA ISHAPE
#   cover_cvcbms_g16_dump.txt CVCBMS  icyc tree ispi D H CL RD DDS TPA \n TRFBMS
#   cover_cvsum_g16_dump.txt  CVSUMHC icyc ip1 ithn j TXHT PCXHT PROXHT VOLXHT CFBXHT

using Test

const _COVER_HERE = @__DIR__
_cover_parse_e22(s) = Float32(parse(Float64, replace(strip(s), "E" => "e")))

# ---- CVCW: crown-area running sum ----
function _cover_load_cvcw(path)
    invs = Vector{Tuple{Int,Int,Vector{Float32},Vector{Float32},Float32}}()
    cw = Float32[]; pr = Float32[]
    for ln in eachline(path)
        f = split(ln)
        if startswith(ln, "CVCWTREE")
            push!(cw, _cover_parse_e22(f[5])); push!(pr, _cover_parse_e22(f[6]))
        elseif startswith(ln, "CVCWFINAL")
            push!(invs, (parse(Int, f[2]), parse(Int, f[3]), copy(cw), copy(pr),
                         _cover_parse_e22(f[4])))
            empty!(cw); empty!(pr)
        end
    end
    return invs
end

@testset "COVER — CVCW crown-area dump-replay (EM, bit-exact vs FVSem_g16)" begin
    invs = _cover_load_cvcw(joinpath(_COVER_HERE, "cover_cvcw_g16_dump.txt"))
    @test !isempty(invs)
    cv = FVSjl.CoverState()
    for (icyc, itrn, cw, pr, orc) in invs
        @test length(cw) == itrn
        jl = FVSjl.cover_cvcw!(cv, cw, pr, itrn)
        @test reinterpret(UInt32, jl) == reinterpret(UInt32, orc)   # hex-exact
    end
end

# ---- CVSHAP: crown-shape discriminant argmax → ISHAPE ----
@testset "COVER — CVSHAP crown-shape dump-replay (EM, 351/351 ISHAPE bit-exact)" begin
    n = 0; bad = 0
    for ln in eachline(joinpath(_COVER_HERE, "cover_cvshap_g16_dump.txt"))
        startswith(ln, "CVSHAP") || continue
        f = split(ln)
        isp  = parse(Int, f[4])
        dbh  = _cover_parse_e22(f[6]); ht = _cover_parse_e22(f[7])
        cr   = _cover_parse_e22(f[8]); rad = _cover_parse_e22(f[10])
        tpa  = _cover_parse_e22(f[11]); ish_dump = parse(Int, f[12])
        ispi = FVSjl._CV_MAP_EM[isp]
        ish  = FVSjl._cover_cvshap(ispi, dbh, ht, cr, rad, tpa)
        n += 1
        ish == ish_dump || (bad += 1)
    end
    @test n == 351
    @test bad == 0            # ISHAPE (integer argmax) is bit-exact for every tree
end

# ---- CVCBMS: COVOPT=2 foliage biomass (transcendental; 1-ULP libm tolerance) ----
@testset "COVER — CVCBMS foliage-biomass dump-replay (EM, hex / 1-ULP libm)" begin
    lines = collect(eachline(joinpath(_COVER_HERE, "cover_cvcbms_g16_dump.txt")))
    n = 0; hexbad = 0; ulpmax = 0; k = 0
    while k < length(lines)
        k += 1; ln = lines[k]
        startswith(ln, "CVCBMS") || continue
        f = split(ln)
        ispi = parse(Int, f[4]); d = _cover_parse_e22(f[5]); h = _cover_parse_e22(f[6])
        cl = _cover_parse_e22(f[7]); rd = _cover_parse_e22(f[8]); dds = _cover_parse_e22(f[9])
        tpa = _cover_parse_e22(f[10])
        k += 1; trf_dump = _cover_parse_e22(lines[k])
        trf = FVSjl._cover_cvcbms2(ispi, d, h, cl, rd, dds, tpa)
        n += 1
        if reinterpret(UInt32, trf) != reinterpret(UInt32, trf_dump)
            hexbad += 1
            ulpmax = max(ulpmax,
                abs(Int(reinterpret(UInt32, trf)) - Int(reinterpret(UInt32, trf_dump))))
        end
    end
    @test n == 351
    @test hexbad <= 14        # 337/351 hex-exact; 14 differ by exactly 1 ULP (exp/log)
    @test ulpmax <= 1         # never more than a single ULP — cannot move the integer report
end

# ---- CVSUM: per-height-class canopy geometry (frustum sums) ----
@testset "COVER — CVSUM height-class geometry dump-replay (EM, 239/240 bit-exact, 1@1-ULP)" begin
    P((ic,i)) = (ic, i)
    trecw = Dict{Tuple{Int,Int},Float32}(); prob = Dict{Tuple{Int,Int},Float32}()
    dbh = Dict{Tuple{Int,Int},Float32}(); ht = Dict{Tuple{Int,Int},Float32}()
    icrp = Dict{Tuple{Int,Int},Int}(); ishape = Dict{Tuple{Int,Int},Int}()
    trfbms = Dict{Tuple{Int,Int},Float32}(); ntree = Dict{Int,Int}()
    for ln in eachline(joinpath(_COVER_HERE, "cover_cvcw_g16_dump.txt"))
        startswith(ln, "CVCWTREE") || continue
        f = split(ln); ic = parse(Int, f[2]); i = parse(Int, f[3])
        trecw[(ic,i)] = _cover_parse_e22(f[5]); prob[(ic,i)] = _cover_parse_e22(f[6])
    end
    for ln in eachline(joinpath(_COVER_HERE, "cover_cvshap_g16_dump.txt"))
        startswith(ln, "CVSHAP") || continue
        f = split(ln); ic = parse(Int, f[2]); i = parse(Int, f[3])
        dbh[(ic,i)] = _cover_parse_e22(f[6]); ht[(ic,i)] = _cover_parse_e22(f[7])
        icrp[(ic,i)] = round(Int, _cover_parse_e22(f[8]) * 100f0)
        ishape[(ic,i)] = parse(Int, f[12]); ntree[ic] = max(get(ntree, ic, 0), i)
    end
    let lines = collect(eachline(joinpath(_COVER_HERE, "cover_cvcbms_g16_dump.txt"))), k = 0
        while k < length(lines)
            k += 1; ln = lines[k]; startswith(ln, "CVCBMS") || continue
            f = split(ln); ic = parse(Int, f[2]); i = parse(Int, f[3]); k += 1
            trfbms[(ic,i)] = _cover_parse_e22(lines[k])
        end
    end
    gold = Dict{Tuple{Int,Int},NTuple{5,Float32}}()
    for ln in eachline(joinpath(_COVER_HERE, "cover_cvsum_g16_dump.txt"))
        startswith(ln, "CVSUMHC") || continue
        f = split(ln); ic = parse(Int, f[2]); j = parse(Int, f[5])
        gold[(ic,j)] = (_cover_parse_e22(f[6]), _cover_parse_e22(f[7]), _cover_parse_e22(f[8]),
                        _cover_parse_e22(f[9]), _cover_parse_e22(f[10]))
    end

    tot = 0; bad = 0; ulpmax = 0
    for ic in sort(collect(keys(ntree)))
        txht = zeros(Float32,16); crxht = zeros(Float32,16); proxht = zeros(Float32,16)
        volxht = zeros(Float32,16); cfbxht = zeros(Float32,16); pcxht = zeros(Float32,16)
        for i in 1:ntree[ic]
            FVSjl._cover_cvsum_tree!(txht, crxht, proxht, volxht, cfbxht,
                                     ishape[(ic,i)], trecw[(ic,i)], dbh[(ic,i)], ht[(ic,i)],
                                     icrp[(ic,i)], prob[(ic,i)], trfbms[(ic,i)])
        end
        for j in 1:16
            pcxht[j] = 100f0*(1f0 - exp(-0.01f0*(crxht[j]/435.6f0)))
        end
        for j in 1:16
            g = gold[(ic,j)]
            for (a, gv) in ((txht[j],g[1]), (pcxht[j],g[2]), (proxht[j],g[3]),
                            (volxht[j],g[4]), (cfbxht[j],g[5]))
                tot += 1
                if reinterpret(UInt32,a) != reinterpret(UInt32,gv)
                    bad += 1
                    ulpmax = max(ulpmax, abs(Int(reinterpret(UInt32,a)) - Int(reinterpret(UInt32,gv))))
                end
            end
        end
    end
    @test tot == 240
    @test bad <= 1            # 239/240 hex-exact; 1 PROXHT value @1-ULP (accumulated PAREA)
    @test ulpmax <= 1
end
