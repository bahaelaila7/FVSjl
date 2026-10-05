# test_mortmsb.jl — MORTMSB / MSBMRT alternate "mature-stand breakup" mortality.
#
# MORTMSB (base keyword option 137) inflicts EXTRA mortality, concentrated in a DBH range, once a stand's
# projected QMD passes QMDMSB — simulating the break-up of overmature stands. It sets a target additional
# kill TMORE = TN − (CONST·(D10/QMDMSB)^SLPMSB·QMDMSB^−1.605)·PMSDIU and removes it across records in
# [DLOMSB, DHIMSB) from above/below/throughout (MFLMSB) at efficiency EFFMSB, via base/msbmrt.f.
#
# The integration scenario is a dense unthinned loblolly stand (dense_long) with a STEEP MORTMSB curve
# (QMDMSB=10, SLPMSB=−10): MSB fires from 2025 on, breaking up the large-tree cohort (QMD collapses from
# 10.6 to 5.7 as the [8,40") trees are killed). Validated BIT-EXACT on TPA/BA/QMD every cycle vs live FVSsn
# (mortmsb.sum.save); the volume columns carry the usual ±1–2 Float32 transcendental noise (present pre-MSB).

using Test, FVSjl

const _MSB_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_msb_rows(txt) = [split(l) for l in split(txt, "\n")
                  if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                                y !== nothing && 1980 < y < 2110)]
_msb_base(path) = [split(l) for l in eachline(path)
                   if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                                 y !== nothing && 1980 < y < 2110)]

@testset "MORTMSB alternate mortality (mature-stand breakup)" begin
    # 1. UNIT — keyword parse/validation (initre.f:13700-13714: bad field ⇒ reset all to MSB-off defaults).
    mkrec(vals) = begin
        fields = fill("", 12); values = zeros(Float32, 12); present = falses(12)
        for (i, v) in enumerate(vals)
            v === nothing && continue
            fields[i] = string(v); values[i] = Float32(v); present[i] = true
        end
        FVSjl.KeywordRecord("MORTMSB", "", fields, values, present, 1, FVSjl.KW_OK, 0)
    end
    s = FVSjl.StandState(FVSjl.Southern())
    # valid record sets all six params
    FVSjl.kw_mortmsb!(s, mkrec((10.0, -10.0, 0.9, 8.0, 40.0, 3.0)))
    @test s.control.msb_qmd   == 10f0
    @test s.control.msb_slope == -10f0
    @test s.control.msb_eff   == 0.9f0
    @test s.control.msb_dlo   == 8f0
    @test s.control.msb_dhi   == 40f0
    @test s.control.msb_flag  == Int32(3)
    # an out-of-range slope (> −1.605) is invalid ⇒ ALL params reset to defaults (MSB disabled)
    s2 = FVSjl.StandState(FVSjl.Southern())
    FVSjl.kw_mortmsb!(s2, mkrec((10.0, -0.5, 0.9, 8.0, 40.0, 3.0)))
    @test s2.control.msb_slope == 0f0 && s2.control.msb_qmd == 999f0 && s2.control.msb_flag == Int32(1)
    # an inverted DBH range (DLO ≥ DHI) is invalid ⇒ reset
    s3 = FVSjl.StandState(FVSjl.Southern())
    FVSjl.kw_mortmsb!(s3, mkrec((10.0, -10.0, 0.9, 40.0, 8.0, 3.0)))
    @test s3.control.msb_slope == 0f0 && s3.control.msb_dhi == 999f0

    # 2. FORTRAN — the dense overmature scenario vs live FVSsn (.sum.save), every cycle.
    key  = joinpath(_MSB_DIR, "mortmsb.key")
    base = _msb_base(joinpath(_MSB_DIR, "mortmsb.sum.save"))
    got  = _msb_rows(FVSjl.run_keyfile(key))
    @test length(got) == length(base) && !isempty(base)
    fired = false
    for (g, b) in zip(got, base)
        @test g[1] == b[1]                      # year
        @test g[3] == b[3]                      # TPA  — bit-exact (the MSB-killed count)
        @test g[5] == b[5]                      # BA   — bit-exact
        @test g[8] == b[8]                      # QMD  — bit-exact (breakup collapses it)
        # volume columns: col11 (SCuFt) is BIT-EXACT (measured Δ0 all cycles) → ==. col9/col10 (TCuFt/MCuFt)
        # carry a grown-cycle ±1 print flip (internal ~0.25: jl col9 2070 1554.648 vs live 1554 — ~2500× a sum-
        # order ULP, so NOT sum order). ★ LIVE-CONFIRMED (2026-07-06, FVS_TreeList DBS NOTRIPLE, matched by TreeId):
        # mortmsb's DETERMINISTIC path is BIT-EXACT at EVERY cycle 1995-2110 (0/27 DBH + crown, max|Δ|=0) — incl. the
        # stochastic dgscor! that fires under NOTRIPLE. So (unlike cst01, whose "DGSCOR" turned out to be a fixable
        # crown band-aid) mortmsb has NO hidden deterministic bug: the residual is PURELY the tripling×DGF-seed
        # interaction (the deterministic-tripling frm + SIGMAR spread amplify a sub-ULP grown-DBH Float32 accumulation
        # into the tripled records — same mechanism as treeszcp). A PERMITTED primitive (transcendental/DGF Float32
        # floor × tripling), positively confirmed. TCuFt/MCuFt (dbh²-driven) show it; SCuFt (merch-top-capped) rounds
        # bit-exact. (Corrects the prior "DGSCOR serial-correlation" label — NOTRIPLE bit-exact proves DGSCOR is fine.)
        @test parse(Float32, g[11]) == parse(Float32, b[11])   # SCuFt — BIT-EXACT
        # confirm the breakup actually fired (a year where TPA drops far more than ordinary self-thinning)
        parse(Int, b[1]) == 2025 && parse(Float32, b[3]) < 200f0 && (fired = true)
    end
    # col9/col10 (TCuFt/MCuFt) ±1 print flip = the DGSCOR diameter serial-correlation floor (see the corrected
    # verdict above): a sub-ULP grown DBH → dbh²-driven cubic Δ~0.02% → adjacent print integer. A permitted COR
    # primitive (WK3/DGSCOR), NOT sum-order and NOT the SN HTGF transcendentals (now fpow/fexp/flog-routed, inert here).
    # 2026-09-26: the ±1 flip was NOT a DGSCOR floor — the calibration DENSE's PCT (pctile.f TOT/PCTIN1 order, WK5=D·(D·P))
    # and read-order IND1 sums now match live, and col9/col10 are bit-exact every cycle.
    @test all(parse(Float32, g[9])  == parse(Float32, b[9])  for (g, b) in zip(got, base))  # col9 TCuFt
    @test all(parse(Float32, g[10]) == parse(Float32, b[10]) for (g, b) in zip(got, base))  # col10 MCuFt
    @test fired   # the test must exercise the MSBMRT path, not a degenerate no-fire stand

    # 3. FORTRAN, PER RECORD — FVS_TreeList (TREELIDB) of the same run vs live FVSsn_g16 (mortmsb_treelist.live.csv,
    # 7 years × every record, Float32 values). Pins the per-record chain the .sum totals hide: MSBMRT kills
    # (TPACLS association, morts.f:647), the top-killed volume height (vols.f:146 NORMHT/100.0), DGF's logf
    # (dgf.f:290), R9CUFT in REAL*4 (r9clark.f:961) and the PROB=0 records (vols.f skip + gradd.f:346 divide-back).
    dir = mktempdir()
    keytxt = replace(read(key, String), "ECHOSUM\n" => "ECHOSUM\nDATABASE\nDSNOUT\nmsbtl.db\nTREELIDB\nEND\n")
    write(joinpath(dir, "mortmsb_tl.key"), keytxt)
    cp(joinpath(_MSB_DIR, "mortmsb.tre"), joinpath(dir, "mortmsb_tl.tre"))
    cd(() -> FVSjl.run_keyfile("mortmsb_tl.key"), dir)
    gl = readlines(joinpath(_MSB_DIR, "mortmsb_treelist.live.csv")); hdr = split(gl[1], ',')
    db = FVSjl.SQLite.DB(joinpath(dir, "msbtl.db"))
    jl = Dict((Int(r.Year), Int(r.TreeIndex)) => Dict(String(c) => r[c] for c in propertynames(r))
              for r in FVSjl.DBInterface.execute(db, "SELECT * FROM FVS_TreeList"))
    nbad = 0; nrow = 0
    for l in gl[2:end]
        f = split(l, ','); k = (parse(Int, f[1]), parse(Int, f[2])); nrow += 1
        r = get(jl, k, nothing)
        r === nothing && (nbad += 1; continue)
        for (j, c) in enumerate(hdr[3:end])
            Float32(r[String(c)]) == parse(Float32, f[j + 2]) || (nbad += 1)
        end
    end
    @test nrow > 1500
    @test nbad == 0
end
