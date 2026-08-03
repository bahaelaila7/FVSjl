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

# bm/NVEL r6vol3.f — Behre total-cubic profile (VOLEQ 616BEH*** → ZONE 1). Smalian-integrated taper
# DR = HRATIO/(0.62·HRATIO+0.38). Inputs: DBHOB, DBTBH (=D·(1-bark)), FCLASS (form class), HTTOT.
# Used for the BM minor species (WP/MH/WJ/WB/LM/PY/YC/AS/CW/OS/OH) whose VEQNNC = 616BEHW.
function bm_r6vol3(dbhob::Float32, dbtbh::Float32, fclass::Int, httot::Float32, zone::Int)::Float32
    topd = 4.0f0; a = 0.62f0; b = 1.0f0 - a
    d17 = Float32(fclass) / 100.0f0 * dbhob
    h17 = zone == 1 ? 17.3f0 : 33.6f0
    dbhib = dbhob - dbtbh
    (dbhib <= 0f0 || dbhib > dbhob) && (dbhib = dbhob)
    httot <= h17 && return 0.00272708f0 * dbhib * dbhib * httot   # small tree (below FC height)
    D = zeros(Float32, 21)
    # butt log (stump→17.3ft): Smalian of DBHIB and D17 over H17
    vol = 0.00272708f0 * (dbhib * dbhib + d17 * d17) * h17
    htup = httot - h17
    D[1] = d17
    s = 0.0f0; hh = 0.0f0; ilast = 2; hit_top = false
    for i in 2:20
        hratio = (htup - ((i - 1) * 16.3f0)) / htup
        hratio <= 0.0f0 && (ilast = i; break)
        dr = hratio / (a * hratio + b)
        D[i] = dr * D[1]
        ilast = i
        if D[i] < topd; break; end
        vol += 0.00272708f0 * (D[i-1] * D[i-1] + D[i] * D[i]) * 16.3f0
        if D[i] == topd; s = 16.3f0; hit_top = true; break; end
    end
    if !hit_top
        dr = topd / d17
        hx = (dr * b * htup) / (1.0f0 - (a * dr))
        hh = (ilast - 2) * 16.3f0
        s = htup - hx - hh
        vol += 0.00272708f0 * (D[ilast-1] * D[ilast-1] + topd * topd) * s
    end
    htup2 = httot - (16.3f0 * (ilast - 2) + h17) - s   # top-of-tree cone
    vol += 0.00272708f0 * (topd * topd) * htup2
    return vol
end
