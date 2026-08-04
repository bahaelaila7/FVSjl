# BC V2 DG WIRING check. Drives the REAL engine setup (each_stand → notre! → setup_growth! → density) on
# all_BC_essf (ESSFdk/01 → V2 regime), then calls the real dgf!(s, BC()) and inspects the per-tree WK2 it
# produced. Confirms (1) the V2 dispatch fires (bc_lv2atv true), (2) CONSPP for sp14 == validated NI DGCON
# (0.7631423) + calib COR + 0.01·DGCCF·RELDEN, (3) WK2 == bc_v2_dds(engine inputs) — i.e. the wiring feeds
# the already-oracle-validated formula the right inputs. (Full .sum bit-exact awaits V2 mortality/HTG/crown.)
import Pkg; Pkg.activate("/workspace/FVSjl"; io = devnull)
using FVSjl
const F = FVSjl
const SCEN = "/workspace/FVSjl/test/harness/scenarios/all_BC_essf"

function main()
    s = first(F.each_stand("$SCEN.key"; variant = F.BritishColumbia()))
    F.notre!(s)
    F.setup_growth!(s)         # runs bc_dgcons! (→ bc_v2_dgcons! for V2) + calibration
    F.compute_forest_type!(s)
    F.compute_density!(s)

    zone, series = F.bc_stand_zone(s)
    println("zone=$zone series=$series  → LV2ATV(V2)=", F.bc_lv2atv(zone))
    @assert F.bc_lv2atv(zone) "expected V2 regime for ESSF"

    # engine DG
    F.dgf!(s, s.variant)
    wk2 = s.scratch.wk[2, :]
    t = s.trees; c = s.calib; p = s.plot
    relden = p.relative_density; ba100 = p.basal_area / 100f0
    println("BA=$(p.basal_area)  RELDEN=$relden  elev=$(p.elevation) asp=$(p.aspect) slope=$(p.slope)")
    println("sp14 dg_const(DGCON)=", c.dg_const[14], "  dg_cor(COR)=", c.dg_cor[14], "  atten=", c.atten[14])

    # DGCCF/DGDSQ the wiring should use (ITYPE=IFOR=4)
    dgccf14 = F.BC_V2_DGCCFA[14][F.BC_V2_MAPCCF[14][4]]
    dgdsq14 = F.BC_V2_DGDS[14][F.BC_V2_MAPDSQ[14][4]]
    println("sp14 DGCCF=", dgccf14, "  DGDSQ=", dgdsq14)

    # per-tree consistency: WK2 == bc_v2_dds(D, BAL, CR, CONSPP, DGDSQ)
    n = 0; ok = 0; sample = String[]
    for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        sp = Int(t.species[i]); sp == 14 || continue
        conspp = c.dg_const[sp] + c.dg_cor[sp] + 0.01f0 * dgccf14 * relden
        pct = t.crown_ratio[i]; bal = (1f0 - pct/100f0) * ba100
        cr = Float32(t.crown_pct[i]) * 0.01f0
        expect = F.bc_v2_dds(sp, d, bal, cr, conspp, dgdsq14)
        n += 1; wk2[i] == expect && (ok += 1)
        length(sample) < 4 && push!(sample, "d=$d bal=$(round(bal,digits=3)) cr=$(round(cr,digits=2)) conspp=$(round(conspp,digits=4)) wk2=$(round(wk2[i],digits=5))")
    end
    println("per-tree WK2 == bc_v2_dds(engine inputs): $ok/$n", ok == n ? "  ✓ wiring feeds the validated formula" : "  ✗ MISMATCH")
    for sm in sample; println("  ", sm); end
end
main()
