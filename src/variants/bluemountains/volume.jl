# =============================================================================
# volume.jl (bluemountains) — BM volume (bm VEQNNC). Chunk 8.
#   FW2W  → Flewelling FW2 (cr_fw2_vol, same as KT/EM/UT) — the main conifers WL/DF/GF/LP/ES/AF/PP.
#   616BEHW → region-6 Behre (minor species WP/MH/WJ/WB/LM/PY/YC/AS/CW/OS/OH — not in bmt01, DEFERRED).
# Merch (bm/grinit.f): TOPD=4.5, DBHMIN=7 (LP sp7=6), BFTOPD=4.5, BFMIND=7 (sp7=6).
# =============================================================================

function compute_volumes_bm!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq
    sd = s.coef.species
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = veq[sp]; se = strip(eq); mdl = length(se) >= 7 ? se[4:6] : "   "
        bark = bm_bratio(sd, sp, d)
        dbhmin = sp == 7 ? 6.0f0 : 7.0f0
        if mdl == "FW2"
            v = cr_fw2_vol(eq, d, h; bark = bark, topd = 4.5f0, bftopd = 4.5f0, stump = 1f0, iregn = 6)
            t.cuft_vol[i] = max(v[1], 0f0)
            t.merch_cuft_vol[i] = d >= dbhmin ? max(v[4] + v[7], 0f0) : 0f0
            t.saw_cuft_vol[i] = 0f0
            t.bdft_vol[i] = d >= dbhmin ? max(v[2], 0f0) : 0f0
        else                                                 # 616BEHW (region-6 Behre) — deferred (not in bmt01)
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0
        end
    end
    return s
end
