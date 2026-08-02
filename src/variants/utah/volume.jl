# =============================================================================
# volume.jl (utah) — UT volume (ut VEQNNC). Chunk 8. Three families, all reused:
#   400/402MATW → R4 Matney taper (r4vol_volumes/r4vol_board, same as TT)
#   407FW2W     → Flewelling FW2 (cr_fw2_vol, same as KT/EM; BS/ES only)
#   400/300DVEW → Chojnacky INT-339 D2H woodland (r4d2h_vol1, same as TT; PJ species)
# Merch standards (ut/grinit.f): TOPD=6, DBHMIN=8 (sp7 LP=7), BFTOPD=6, BFMIND=8 (sp7=7).
# =============================================================================

function compute_volumes_ut!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq
    ba_a = s.calib.bark_a; ba_b = s.calib.bark_b
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = veq[sp]; se = strip(eq); mdl = length(se) >= 7 ? se[4:6] : "   "
        bark = bark_ratio(ba_a, ba_b, sp, d)
        dbhmin = sp == 7 ? 7.0f0 : 8.0f0
        if mdl == "MAT"
            mtopp = 6.0f0 * bark                             # TOPD=6 outside-bark → inside-bark top
            tcf, mcf = r4vol_volumes(eq, d, h, mtopp, 0f0)
            t.cuft_vol[i] = max(tcf, 0f0)
            t.merch_cuft_vol[i] = d >= dbhmin ? max(mcf, 0f0) : 0f0
            bf = d >= dbhmin ? r4vol_board(eq, d, h, mtopp, 0f0) : 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = max(bf, 0f0)
        elseif mdl == "FW2"
            v = cr_fw2_vol(eq, d, h; bark = bark, topd = 6.0f0, bftopd = 6.0f0, stump = 1f0, iregn = 4)
            t.cuft_vol[i] = max(v[1], 0f0)
            t.merch_cuft_vol[i] = d >= dbhmin ? max(v[4] + v[7], 0f0) : 0f0
            t.saw_cuft_vol[i] = 0f0
            t.bdft_vol[i] = d >= dbhmin ? max(v[2], 0f0) : 0f0
        else                                                 # DVE woodland (r4d2h, region 4)
            vol1 = r4d2h_vol1(eq, d, h)
            t.cuft_vol[i] = max(vol1, 0f0)
            t.merch_cuft_vol[i] = d >= dbhmin ? max(vol1, 0f0) : 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0
        end
    end
    return s
end
