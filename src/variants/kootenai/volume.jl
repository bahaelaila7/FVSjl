# =============================================================================
# volume.jl (kootenai) — KT volume (chunk 8). Region-1 Flewelling FW2 (NVEL), the SAME kernel CR uses
# (cr_fw2_vol) — KT's VOLEQ "I00FW2W<FIA>" maps through _fw2_jsp's 'I'/INGY geocode branch (already ported).
# Merch standards from kt/grinit.f: cubic TOPD=4.5/DBHMIN=7 (sp7 lodgepole=6); board BFTOPD=4.5/BFMIND=7
# (sp7=6); stump=1. KT bark = bark_ratio(KT_BKRAT) (kt/bratio.f), NOT cr_bratio.
# =============================================================================
function compute_volumes_kt!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq; sd = s.coef.species
    ba_a = s.calib.bark_a; ba_b = s.calib.bark_b
    iregn = 1                                    # KT = Northern (Region 1)
    topd = 4.5f0; bftopd = 4.5f0; stump = 1.0f0
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        dbhmin = sp == 7 ? 6f0 : 7f0             # kt/grinit.f (lodgepole sp7 = 6)
        bfmind = sp == 7 ? 6f0 : 7f0
        bark = bark_ratio(ba_a, ba_b, sp, d)
        v = startswith(veq[sp], "I") ?
            cr_fw2_vol(veq[sp], d, h; bark = bark, topd = topd, bftopd = bftopd, stump = stump, iregn = iregn) :
            zeros(Float32, 15)
        tcf = max(v[1], 0f0)
        mcf = d >= dbhmin ? max(v[4] + v[7], 0f0) : 0f0
        bf  = d >= bfmind ? v[2] : 0f0
        if t.trunc[i] > 0 && tcf > 0f0 && h >= 4.5f0     # broken-top reduction (Behre taper)
            vmax = tcf
            tcf, mcf = cr_cftopk(tcf, mcf, d, h, vmax, bark, Int(t.trunc[i]), stump, topd)
            bf = cr_bftopk(bf, d, h, vmax, bark, Int(t.trunc[i]), stump, bftopd)
        end
        t.cuft_vol[i] = tcf; t.merch_cuft_vol[i] = mcf
        t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = max(bf, 0f0)
    end
    return s
end
