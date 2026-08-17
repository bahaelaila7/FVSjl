# test_bc_dmntrd.jl — BC NEWSPRED (canada/newmist) DMNTRD multi-cycle crown-third infection
# remap (dmntrd.f) + the tripling-offspring DMR propagation (base triple.f MISGET/MISPUT). #196.
#
# DMNTRD re-bins each tree's per-(crown-third × life-history pool) DMINF loads onto the MOVED
# crown-third breakpoints every cycle: DM infections are FIXED on branches (do not move within the
# crown), but height growth pushes fresh uninfected crown up top and crown recession removes the
# lower margin. By locating where each PREVIOUS breakpoint (ms.pbrkpt) now lies among the CURRENT
# breakpoints (ms.brkpnt), the infection density is re-apportioned (with a `Mult` density scalar for
# crown enlargement). Gated to cyc≥2 by ms.ntdn. The tests below hand-verify the remap arithmetic
# against the Fortran (dmntrd.f:84-199), the Mult scaling, and the Σpbrkpt==0 skip guard.
using Test, FVSjl
const _DMN = FVSjl

# Minimal BC stand carrying an ACTIVE mistletoe state sized to `n` trees. dm_ntrd! only reads
# s.trees.n and the mistletoe pbrkpt/brkpnt/dminf, so the rest of the stand can stay at defaults.
function _dmn_stand(n::Int)
    s = _DMN.StandState(_DMN.BritishColumbia())
    s.trees.n = n
    ms = _DMN.MistletoeState(); ms.active = true; ms.newmod = true
    ms.dmr    = zeros(Int32, n)
    ms.dminf  = zeros(Float32, n, _DMN.DM_CRTHRD, _DMN.DM_NPOOL)
    ms.brkpnt = zeros(Float32, n, _DMN.DM_BPCNT)
    ms.pbrkpt = zeros(Float32, n, _DMN.DM_BPCNT)
    s.mistletoe = ms
    return s, ms
end

@testset "BC DMNTRD crown-third infection remap (dmntrd.f)" begin
    A = _DMN.DM_ACTIVE

    @testset "no-movement invariant (pbrkpt==brkpnt ⇒ Mult=1, DMINF unchanged)" begin
        s, ms = _dmn_stand(1)
        ms.pbrkpt[1, :] .= Float32[12, 8, 4, 0]
        ms.brkpnt[1, :] .= Float32[12, 8, 4, 0]
        ms.dminf[1, 1, A] = 100f0; ms.dminf[1, 2, A] = 200f0; ms.dminf[1, 3, A] = 300f0
        _DMN.dm_ntrd!(s)
        # Each old breakpoint maps to its own new third (identity); density unchanged.
        @test ms.dminf[1, 1, A] ≈ 100f0
        @test ms.dminf[1, 2, A] ≈ 200f0
        @test ms.dminf[1, 3, A] ≈ 300f0
    end

    @testset "vertical crown shift (grow up + recede at base), Mult=1" begin
        # prev crown [12,8,4,0] (len 12); new crown shifted up 6 → [18,14,10,6] (len 12). Mult=1.
        # Hand-derivation (dmntrd.f): old bp1→new third2 split .5/.5 into thirds 2&3; old bp2→third3
        # split; old bp3 (at 4) is below the new crown bottom (6) ⇒ OUTGROWN (k=0), dropped.
        #   NewVal[1]=0 (fresh uninfected top crown)
        #   NewVal[2]=Old[1]*0.5 = 50
        #   NewVal[3]=Old[1]*0.5 + Old[2]*0.5 = 50 + 100 = 150
        s, ms = _dmn_stand(1)
        ms.pbrkpt[1, :] .= Float32[12, 8, 4, 0]
        ms.brkpnt[1, :] .= Float32[18, 14, 10, 6]
        ms.dminf[1, 1, A] = 100f0; ms.dminf[1, 2, A] = 200f0; ms.dminf[1, 3, A] = 300f0
        _DMN.dm_ntrd!(s)
        @test ms.dminf[1, 1, A] ≈ 0f0
        @test ms.dminf[1, 2, A] ≈ 50f0
        @test ms.dminf[1, 3, A] ≈ 150f0
    end

    @testset "Mult density scaling (crown doubles in length ⇒ density halved, mass conserved)" begin
        # prev [12,8,4,0] (len 12); new [24,16,8,0] (len 24, base fixed). Mult=(12-0)/(24-0)=0.5.
        #   NewVal[1]=0; NewVal[2]=Old[1]=100*0.5=50; NewVal[3]=(Old[2]+Old[3])=(200+300)*0.5=250.
        # Total after = 300 = 0.5*(100+200+300): absolute infection preserved as crown enlarges.
        s, ms = _dmn_stand(1)
        ms.pbrkpt[1, :] .= Float32[12, 8, 4, 0]
        ms.brkpnt[1, :] .= Float32[24, 16, 8, 0]
        ms.dminf[1, 1, A] = 100f0; ms.dminf[1, 2, A] = 200f0; ms.dminf[1, 3, A] = 300f0
        _DMN.dm_ntrd!(s)
        @test ms.dminf[1, 1, A] ≈ 0f0
        @test ms.dminf[1, 2, A] ≈ 50f0
        @test ms.dminf[1, 3, A] ≈ 250f0
        @test sum(ms.dminf[1, :, A]) ≈ 300f0
    end

    @testset "all life-history pools remapped, not just ACTIVE" begin
        s, ms = _dmn_stand(1)
        ms.pbrkpt[1, :] .= Float32[12, 8, 4, 0]
        ms.brkpnt[1, :] .= Float32[24, 16, 8, 0]
        for p in 1:_DMN.DM_NPOOL
            ms.dminf[1, 1, p] = 100f0; ms.dminf[1, 2, p] = 200f0; ms.dminf[1, 3, p] = 300f0
        end
        _DMN.dm_ntrd!(s)
        for p in 1:_DMN.DM_NPOOL
            @test ms.dminf[1, 2, p] ≈ 50f0
            @test ms.dminf[1, 3, p] ≈ 250f0
        end
    end

    @testset "Σpbrkpt==0 skip guard (fresh/tripled record, no remap ÷0)" begin
        # A newborn record has pbrkpt≡0 (MISPUTZ zeros PBRKPT); DMNTRD must skip it (Fortran GOTO 100
        # divide-by-zero guard), leaving its pools untouched this cycle.
        s, ms = _dmn_stand(1)
        ms.pbrkpt[1, :] .= 0f0
        ms.brkpnt[1, :] .= Float32[18, 14, 10, 6]
        ms.dminf[1, 3, A] = 42f0
        _DMN.dm_ntrd!(s)
        @test ms.dminf[1, 3, A] ≈ 42f0   # unchanged (skipped)
    end

    @testset "ms.ntdn default (first DMTREG entry gates DMNTRD off)" begin
        @test _DMN.MistletoeState().ntdn == false
    end
end

@testset "BC tripling / regen offspring DMR propagation (MISPUT)" begin
    # base triple.f:110-111 does MISGET(parent)→MISPUT(offspring): the offspring inherits the parent's
    # DMR (IMIST/DMRATE) but its DMINF pools + breakpoints are loaded from DMKLDG — a COMMON scratch
    # array that is NEVER populated ⇒ all zeros. jl's copy_tree! already carries t.dmr to offspring;
    # _dm_triple_carry! re-syncs ms.dmr for the appended records, and _dm_ensure_capacity! zero-fills
    # their pools — reproducing the exact FVS post-tripling state (parent DMR, zero pools).
    s, ms = _dmn_stand(2)
    ms.dmr = Int32[3, 0]
    ms.dminf[1, 3, _DMN.DM_ACTIVE] = 5f0            # parent-1 carries infection in its bottom third
    # simulate triple_records!: nlive=2 parents stay at 1,2; offspring appended at 3..6 with the
    # parent's t.dmr (copy_tree!). parent-1 → records 3,4 (DMR 3); parent-2 → records 5,6 (DMR 0).
    nold = length(ms.dmr)
    s.trees.n = 6
    s.trees.dmr[1] = 3; s.trees.dmr[2] = 0
    s.trees.dmr[3] = 3; s.trees.dmr[4] = 3; s.trees.dmr[5] = 0; s.trees.dmr[6] = 0
    _DMN._dm_ensure_capacity!(ms, 6)
    _DMN._dm_triple_carry!(ms, s.trees, nold)
    @test ms.dmr == Int32[3, 0, 3, 3, 0, 0]                    # offspring inherit parent DMR
    @test ms.dminf[1, 3, _DMN.DM_ACTIVE] ≈ 5f0                 # parent infection preserved
    @test all(ms.dminf[3:6, :, :] .== 0f0)                    # offspring pools ZERO (MISPUT/DMKLDG=0)
    @test all(ms.pbrkpt[3:6, :] .== 0f0)                      # offspring pbrkpt ZERO ⇒ DMNTRD skips them
    @test size(ms.dminf, 1) == 6 && size(ms.pbrkpt, 1) == 6   # capacity grew, parents preserved
end

# YSM029-271 end-to-end DMNTRD trajectory vs the live-FVS oracle (cornered). The reproducer reads
# the BC metric DB; skipped if the DB fixture is absent. jl reproduces the oracle's DM-FREE cyc0
# ingest bit-exact (TPA 2300 / SDI 311). The multi-cycle divergence (2137 jl SDI 1303 vs oracle 936)
# is the documented BC young-dense-lodgepole self-thin / OLDRN straddle (#206 class), NOT the DM:
# the oracle's OWN DM dump (YSM-SkyRanchDBdump FVS_DM_Stnd_Sum_Metric) kills only ~113 TPH total
# (~7% of mortality) with Mean_DMI held at 1-2 the whole run — a MINOR effect jl matches. DMNTRD is
# an infection-conserving REMAP, so on this NOTRIPLE sparse-seed stand it is inert on the aggregate
# .sum (A/B verified byte-identical vs DMNTRD-off). The checkpoints below guard the DMNTRD-inclusive
# engine against regression.
@testset "BC NEWSPRED YSM DMNTRD trajectory (cornered vs oracle DM dump)" begin
    key = "/workspace/.bcwork/newspred/ysm271_dm.key"
    db  = "/workspace/.bcwork/newspred/FVS-BC.YSM-SkyRanch.db"
    if !isfile(key) || !isfile(db)
        @test_skip "YSM reproducer DB fixture not available"
    else
        out = FVSjl.run_keyfile(key; variant = FVSjl.BritishColumbia())
        rows = [split(l) for l in split(out, "\n")
                if length(split(l)) >= 5 && (y = tryparse(Int, first(split(l)));
                                             y !== nothing && 2000 < y < 2200)]
        @test length(rows) == 25                              # 24 cycles + inventory
        yr(r) = parse(Int, r[1]); tpa(r) = parse(Int, r[3]); sdi(r) = parse(Int, r[5])
        c0 = rows[1]
        @test yr(c0) == 2018 && tpa(c0) == 2300 && sdi(c0) == 311   # DM-free ingest = oracle bit-exact
        cN = rows[end]
        @test yr(cN) == 2137
        @test tpa(cN) == 1075 && sdi(cN) == 1303              # deterministic DMNTRD checkpoint (cornered)
    end
end
