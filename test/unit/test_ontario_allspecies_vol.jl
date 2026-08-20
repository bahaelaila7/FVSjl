# =============================================================================
# test_ontario_allspecies_vol.jl — cubic volume (volont ZAK/HONER) bit-exact across ALL 72 ON species.
#
# Companion to the DGF/htont all-species tests. ON's default cubic volume (METHC=8) routes per species
# through either the Zakrzewski (ZAK) or Honer special equations, selected by ON_IS_ZAK[sp]. The ont01
# stand exercised only 8 species; this locks the whole set. An instrumented FVSon_voldump (canada/on/
# vols.f, unit 774; .sum data rows byte-identical to clean) dumps per tree D / H / VN (total merch cuft)
# / VM (sawlog cuft) / BARK. The SHIPPED on_zakvol / on_zakht / on_honer reproduce VN and VM bit-exact
# (Float32-hex) for every distinct species.
#
# Golden = test/fixtures/ontario/ont_all_voldump.txt. Validated 2026-08-20: VN 70/70, VM 70/70 (70
# distinct ISPC — the 3 willow alpha codes WI collapse to one species on resolve). Merch standards
# TOPD = 10 cm top / STMP = 30 cm stump (IFOR 9 defaults, same STDINFO as ont01).
# =============================================================================

using Test, FVSjl
const F = FVSjl

@testset "ON — cubic volume ZAK/HONER bit-exact across all 72 species (vs FVSon_voldump)" begin
    dump = joinpath(@__DIR__, "..", "fixtures", "ontario", "ont_all_voldump.txt")
    if !isfile(dump)
        @test_skip "ont_all volume dump fixture not present"
    else
        h2f(h) = reinterpret(Float32, parse(UInt32, h; base = 16))
        u(x) = reinterpret(UInt32, Float32(x))
        cmToIn = 0.3937f0; cmToFt = 0.0328084f0
        TOPD = 10f0 * cmToIn; STMP = 30f0 * cmToFt
        seen = Set{Int}(); n = okN = okM = 0
        for ln in eachline(dump)
            startswith(strip(ln), "#") && continue
            f = split(strip(ln))
            length(f) >= 7 || continue
            _, sp = parse.(Int, f[1:2])
            d, h, vn_o, vm_o, bark = h2f.(f[3:7])
            sp in seen && continue           # one (cyc0) record per species
            push!(seen, sp)
            if F.ON_IS_ZAK[sp]
                vn = F.on_zakvol(sp, d, h, h, bark, true, STMP)
                zht = F.on_zakht(sp, d, h, bark, TOPD)
                vm = zht > 0f0 ? F.on_zakvol(sp, d, h, zht, bark, false, STMP) : 0f0
            else
                vn, vm = F.on_honer(sp, d, h, TOPD, STMP)
            end
            n += 1
            okN += u(vn) == u(vn_o)
            okM += u(vm) == u(vm_o)
        end
        @test n >= 60                # most of the 72 species (duplicate alpha codes collapse a few)
        @test okN == n               # total merch cubic bit-exact for every species
        @test okM == n               # sawlog cubic bit-exact
    end
end
