# test_readcor.jl — READCOR{D,H,R} / REUSCOR{D,H,R} growth-correction keywords.
#
# These read a block of MAXSP per-species correction terms (8F10.0) that modify the growth
# model CONSTANTS before the LSTART calibration (initre.f:5600/6900/7500):
#   READCORD → COR2  : DGCON += ln(COR2)   large-tree diameter growth (dgf.f:1168)
#   READCORH → HCOR2 : HTCON += ln(HCOR2)  large-tree height growth   (htgf.f:332)
#   READCORR → RCOR2 : RHCON  = RCOR2      small-tree height growth    (regent.f:585)
# REUSE* re-enables the previously-read terms without re-reading (multi-stand carry-over).
# All default off / terms = 1 ⇒ no-op (the rest of the suite, snt01 included, stays bit-exact).
# Checks: (1) the reader fills the array + sets the flag; (2) REUSE just sets the flag; (3) a
# READCORD COR2=1.3 projection matches live Fortran on every structural column (board-feet within
# the documented ±Scribner Float32 transcendental noise), and actually changes the stand.

using Test, FVSjl

const _RC_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_rc_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]
_rc_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                                y !== nothing && 1900 < y < 2100)]

@testset "READCOR / REUSCOR growth corrections" begin
    # 1/2. unit — the reader + REUSE flag, via a tiny in-memory keyword stream.
    mkstream(body) = FVSjl.KeywordReader(IOBuffer(body))
    # READCORD: a COR2 line (sp1 = 1.5, sp3 = 2.0, rest blank ⇒ 0) over 8F10.0 continuation lines.
    cor_line = string(lpad("1.5", 10), lpad("", 10), lpad("2.0", 10), repeat(lpad("", 10), 5))
    body = cor_line * "\n" * repeat(repeat(" ", 80) * "\n", 11)   # 12 lines total (90 sp / 8)
    s = FVSjl.StandState(FVSjl.Southern())
    kr = mkstream(body)
    FVSjl.kw_readcord!(s, kr)
    @test s.control.dg_cor2_on
    @test s.control.dg_cor2[1] == 1.5f0
    @test s.control.dg_cor2[3] == 2.0f0
    @test s.control.dg_cor2[2] == 0f0          # blank field ⇒ 0 (Fortran F10.0), guarded > 0 on apply
    # REUSE just re-enables (no read)
    s2 = FVSjl.StandState(FVSjl.Southern())
    FVSjl.kw_reuscord!(s2); @test s2.control.dg_cor2_on && s2.control.dg_cor2[1] == 1f0
    FVSjl.kw_reuscorh!(s2); @test s2.control.htg_cor2_on
    FVSjl.kw_reuscorr!(s2); @test s2.control.regh_cor2_on

    # 3. FORTRAN — the readcord scenario (COR2 = 1.3 all species) vs live Fortran.
    key = joinpath(_RC_DIR, "readcord.key")
    sav = joinpath(_RC_DIR, "readcord.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "readcord scenario not available"
    else
        jl = _rc_rows(FVSjl.run_keyfile(key; faithful = true))
        ft = _rc_base(sav)
        @test length(jl) == length(ft)
        for (j, f) in zip(jl, ft)
            @test j[1] == f[1]                                   # YEAR
            for col in (3, 4, 5, 6, 7, 9, 10)                    # TPA/BA/SDI/CCF/TopHt/TCuFt/MCuFt
                @test j[col] == f[col]
            end
            @test parse(Int, j[12]) == parse(Int, f[12])  # BdFt — BIT-EXACT (measured Δ=0; was over-cautious ≤2, closed by the BFTOPK fix)
        end
        # the COR2 boost actually changed the stand (vs the same key with READCORD removed)
        offkey = joinpath(_RC_DIR, "_readcord_off.key")
        cp(joinpath(_RC_DIR, "readcord.tre"), joinpath(_RC_DIR, "_readcord_off.tre"); force = true)
        open(offkey, "w") do io
            lines = collect(eachline(key)); i = 1
            while i <= length(lines)
                if strip(lines[i]) == "READCORD"; i += 13; continue; end   # drop kw + 12 data lines
                println(io, lines[i]); i += 1
            end
        end
        try
            off = _rc_rows(FVSjl.run_keyfile(offkey; faithful = true))
            @test any(jl[i][10] != off[i][10] for i in eachindex(jl))      # MCuFt differs ⇒ COR2 fired
        finally
            rm(offkey; force = true); rm(joinpath(_RC_DIR, "_readcord_off.tre"); force = true)
        end
    end
end
