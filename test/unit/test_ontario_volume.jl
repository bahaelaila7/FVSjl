# test_ontario_volume.jl — ON (Ontario) per-tree volume dump-replay BIT-EXACT vs FVSon_g16.
#
# ON's default volume method (METHC=METHB=8) routes through the canada/on volont.f ZAK/HONER
# special equations + the Mowraski net-merch cull. Golden = ontario_volont_dump.txt: an instrumented
# canada/on/vols.f (single-.o swap into FVSon_g16, /workspace/.onwork/FVSon_voldump) dumps per tree,
# for both projection cycles on ont01:
#   I ISPC ICYC  D(in) H(ft) GTV GMV CFV NMV  BARK ABIRTH DBHMIN BFMIND BFTOPD   (all Float32-hex)
# where GTV = total merch cubic ft (VN), GMV = merch sawlog cubic ft (WK1/VM), NMV = Mowraski-culled
# net merch (BFV), BARK = BRATIO, ABIRTH = dubbed tree age.
#
# We replay the ON volume kernel feeding the oracle's exact D/H/BARK/ABIRTH and assert GTV/GMV/NMV are
# bit-exact, then assert on_tree_age reproduces ABIRTH from height alone (the age Mowraski consumes).
# The transcendental chain (Zak taper integral + ZAKHT cube roots + Honer + Carmean age curves) matches
# only through glibc logf/expf/powf and gfortran's (x²)² integer-fourth-power rounding.

using Test

@testset "ON — per-tree volume (ZAK/HONER + Mowraski) dump-replay bit-exact vs FVSon_g16" begin
    _h2f(h) = reinterpret(Float32, parse(UInt32, h; base = 16))
    _u(x) = reinterpret(UInt32, Float32(x))
    # ON merch standards for ont01 (grinit/sitset, IFOR 9): TOPD = STMP-metric defaults.
    cmToIn = 0.3937f0; cmToFt = 0.0328084f0
    TOPD = 10f0 * cmToIn; STMP = 30f0 * cmToFt
    # SITEAR by species (from the wired growth state), needed for the age reproduction check.
    SI = Dict(5 => 58.486f0, 6 => 57.531f0, 26 => 52.16758f0, 28 => 52.16758f0,
              9 => 60.0f0, 11 => 60.0f0, 1 => 61.597f0, 8 => 53.36521f0)

    nrec = okG = okM = okN = okA = nA = 0
    for ln in eachline(joinpath(@__DIR__, "ontario_volont_dump.txt"))
        f = split(strip(ln)); isempty(f) && continue
        sp = parse(Int, f[2]); icyc = parse(Int, f[3])
        d, h, gtv_o, gmv_o, _cfv, nmv_o, bark, ab, dbhmin, bfmind, bftopd = _h2f.(f[4:14])
        if FVSjl.ON_IS_ZAK[sp]
            vn = FVSjl.on_zakvol(sp, d, h, h, bark, true, STMP)
            zht = FVSjl.on_zakht(sp, d, h, bark, TOPD)
            vm = zht > 0f0 ? FVSjl.on_zakvol(sp, d, h, zht, bark, false, STMP) : 0f0
        else
            vn, vm = FVSjl.on_honer(sp, d, h, TOPD, STMP)
        end
        vn < 0f0 && (vn = 0f0); vm < 0f0 && (vm = 0f0)
        nmv = (d >= bfmind && d > bftopd) ? FVSjl.on_mowraski(sp, max(vm, 0f0), ab) : 0f0
        nmv < 0f0 && (nmv = 0f0)
        nrec += 1
        okG += _u(vn) == _u(gtv_o); okM += _u(vm) == _u(gmv_o); okN += _u(nmv) == _u(nmv_o)
        @test _u(vn) == _u(gtv_o)
        @test _u(vm) == _u(gmv_o)
        @test _u(nmv) == _u(nmv_o)
        # Age reproduction from height (findag.f → htcalc.f MODE0=0). ABIRTH is dubbed once from the
        # INVENTORY height (cyc0) and frozen, so only replay it against the cyc0 records.
        if icyc == 0 && haskey(SI, sp) && ab > 0f0
            nA += 1
            age = FVSjl.on_tree_age(sp, h, SI[sp])
            okA += _u(age) == _u(ab)
            @test _u(age) == _u(ab)
        end
    end
    @test nrec == 29
    @test (okG, okM, okN) == (nrec, nrec, nrec)   # GTV / GMV / NMV all bit-exact
    @test okA == nA                                # on_tree_age reproduces ABIRTH
end
