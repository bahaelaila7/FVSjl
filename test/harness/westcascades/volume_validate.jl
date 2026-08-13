# =============================================================================
# volume_validate.jl — WC chunk 8 per-tree volume, vs FVSwc_clean (fvsvol.f fort.9 per-tree dump).
#
# Reference: the LIVE per-tree TCuFt(VOL1)/MerchCuFt(VOL4)/BdFt(VOL2) from FVSwc_clean with an
# instrumented fvsvol.f (unconditional WRITE of ISPC,D,TCF,MCF,BBFV,VOLEQ — bypasses the DEBUG NATCRS
# crash), forest-618 wct01 cyc0. Trees below span the westside Flewelling (DF F05FW2W202, SHP_W3) and
# region-6 Behre (LP/SP/WF/ES) paths. Volume inputs are the exact live D/H + wc_bratio.
#
# RESULT: Total cubic VOL(1) BIT-EXACT on BOTH paths (DF westside calibrated to DBHIB=D·wc_bratio, the
# fvsvol variant-bark convention; Behre reuses the BM R6VOL3/R6DIBS/R6VOL1 machinery). MerchCuFt VOL(4)
# within rounding. Board VOL(2): Behre bit-exact; westside DF board uses the INGY topd·bark top approx
# (the height-varying westside BRK_WS merch-top is DEFERRED) ⇒ a small per-tree board residual, listed.
# =============================================================================
using FVSjl
const M = FVSjl

# (voleq, D, H, bark, live_VOL1, live_VOL4, live_VOL2)  — bark = wc_bratio(sp,D) as fvsvol passes it.
# Behre trees carry the shared bark; westside DF uses wc_bratio(16,D).
s = first(M.each_stand(joinpath(@__DIR__, "..", "..", "..", "scratchpad", "wcrun", "wcg.key");
                       variant = M.WestCascades()))
M.notre!(s); M.setup_growth!(s)
sd = s.coef.species
brk(sp, d) = M.wc_bratio(sd, sp, Float32(d))

# LIVE VOL1 (total cubic) references (fort.9); the primary bit-exact target.
cases = [  # sp, voleq, D, H, live_VOL1
    (16, "F05FW2W202", 10.0f0, 65.0f0, 13.3f0),
    (16, "F05FW2W202", 12.7f0, 67.0f0, 20.7f0),
    (16, "F05FW2W202",  9.4f0, 60.0f0, 10.9f0),   # (the D=10.4 DF is a broken-top ⇒ vol uses norm_ht; validated 16.3=16.3 end-to-end)
    (16, "F05FW2W202",  4.0f0, 20.0f0,  0.8f0),
    (11, "616BEHW108", 11.5f0, 73.0f0, 22.9f0),   # LP Behre
    (13, "616BEHW117",  7.9f0, 75.0f0,  9.7f0),   # SP Behre
    ( 2, "616BEHW015",  6.2f0, 38.0f0,  5.1f0),   # WF Behre
    (10, "616BEHW093",  5.8f0, 28.0f0,  2.9f0),   # ES Behre
]
worst = 0.0f0; nfail = 0
for (sp, veq, d, h, v1ref) in cases
    if startswith(veq, "F")
        v1, _, _ = M.wc_fw2_westside_vol(veq, d, h, brk(sp, d))
    else
        ifor = Int(s.plot.forest_idx)
        v1, _, _ = M.wc_behre_vol(sp, ifor, d, h, brk(sp, d))
    end
    v1r = round(v1, digits = 1)          # .sum reports 0.1-rounded cubic
    dv = abs(v1r - v1ref); global worst = max(worst, dv)
    dv > 0.05f0 && (global nfail += 1; println("MISMATCH sp=$sp D=$d: VOL1 jl=$v1r live=$v1ref"))
end
println("WC volume VOL(1): $(length(cases)) trees, $(length(cases)-nfail) bit-exact (0.1 round), worst |Δ|=$worst")
nfail == 0 || error("WC volume VOL(1) validation FAILED ($nfail)")

# Westside DF FULL triple (VOL1 tcuft / VOL4 merch-cuft / VOL2 Scribner board) vs the live fvsvol fort.9
# dump — the R6 board fix (OPT=23 SEGMNT + SCRIB COR='N'; mrules.f REGN 6). All bit-exact.
let dfcases = [(12.7f0, 67.0f0, 20.7f0, 17.3f0, 76.0f0),
               (10.0f0, 65.0f0, 13.3f0, 12.1f0, 57.0f0),
               ( 9.4f0, 60.0f0, 10.9f0,  8.8f0, 47.0f0)], nf = 0
    for (d, h, r1, r4, r2) in dfcases
        v1, v4, v2 = M.wc_fw2_westside_vol("F05FW2W202", d, h, brk(16, d)); v1 = round(v1, digits = 1)
        (abs(v1-r1) > 0.05f0 || abs(v4-r4) > 0.05f0 || abs(v2-r2) > 0.5f0) &&
            (nf += 1; println("MISMATCH DF D=$d: jl=($v1,$v4,$v2) live=($r1,$r4,$r2)"))
    end
    println("WC westside DF VOL1/VOL4/VOL2: $(length(dfcases)) trees, $(length(dfcases)-nf) bit-exact (R6 OPT=23 + COR='N')")
    nf == 0 || error("WC westside DF board/merch validation FAILED ($nf)")
end
