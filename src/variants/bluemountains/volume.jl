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
    iforst = Int(s.plot.user_forest_code) % 100          # R6 forest number (614→14 Umatilla)
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
        else                                                 # 616BEHW (region-6 Behre total cubic)
            spec = length(se) >= 10 ? String(se[8:10]) : "999"
            fclass = bm_formcl(spec, iforst, d)
            dbtbh = d * (1f0 - bark)                          # double bark thickness (fvsvol.f:153)
            dbhib = d - dbtbh
            v = if h <= 17.3f0                                # R6VOL short-tree guard (TTH≤FC_HT)
                0.00272708f0 * dbhib * dbhib * h
            else
                bm_r6vol3(d, dbtbh, fclass, h, 1)            # ZONE 1 (VOLEQ prefix 616)
            end
            t.cuft_vol[i] = max(v, 0f0)
            t.merch_cuft_vol[i] = 0f0                         # Behre supplies total cubic only
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0
        end
    end
    return s
end

# bm/NVEL formclas.f FORMCL_BM — R6 Blue Mountains form-class lookup. Binary search of the FIA
# code SPEC in FIAJSP (11 sorted codes), IFCDBH = (D−1)/10+1 clamped [1,5] (D>40.9→5), then
# FC = forest_table[ISPC, IFCDBH]. IFORST 4/7/14/16 → Malheur/Ochoco/Umatilla/Wallowa-Whitman;
# any other → Wallowa-Whitman. SPEC not in FIAJSP → FC=80. Tables stored [ifcdbh, ispc].
const BM_FCL_FIAJSP = ["   ", "017", "019", "073", "093", "108", "119", "122", "202", "264", "999"]
const BM_FCL_MALH = Int[80 76 78 78 77 80 78 78 78 75 60;
                        80 78 80 79 80 83 78 78 77 79 60;
                        80 77 80 80 82 83 79 80 77 79 60;
                        80 76 82 82 84 80 81 82 80 79 60;
                        80 76 82 77 84 80 78 83 77 78 60]
const BM_FCL_OCHO = Int[80 76 78 78 82 70 78 76 79 75 60;
                        80 78 76 78 82 75 80 78 79 78 60;
                        80 77 74 80 82 75 80 78 76 79 60;
                        80 74 74 80 82 75 82 80 76 79 60;
                        80 74 74 80 82 75 80 80 76 78 60]
const BM_FCL_UMAT = Int[80 76 74 78 77 86 78 78 77 75 60;
                        80 78 74 78 77 86 78 78 77 75 60;
                        80 77 74 78 75 86 80 80 77 75 60;
                        80 76 75 78 75 86 81 81 77 79 60;
                        80 76 75 78 75 86 81 81 77 78 60]
const BM_FCL_WLWH = Int[80 76 78 78 84 85 78 78 78 75 60;
                        80 78 79 82 84 86 78 78 77 79 60;
                        80 77 79 77 84 85 78 80 77 79 60;
                        80 76 79 75 84 85 78 82 77 79 60;
                        80 76 79 75 84 85 78 83 77 78 60]

function bm_formcl(spec::AbstractString, iforst::Int, d::Real)::Int
    ispc = findfirst(==(spec), BM_FCL_FIAJSP)
    ispc === nothing && return 80                             # not a form-class species
    ifcdbh = Int(floor((Float32(d) - 1f0) / 10f0 + 1f0))
    ifcdbh < 1 && (ifcdbh = 1)
    Float32(d) > 40.9f0 && (ifcdbh = 5)
    tbl = iforst == 4 ? BM_FCL_MALH : iforst == 7 ? BM_FCL_OCHO :
          iforst == 14 ? BM_FCL_UMAT : BM_FCL_WLWH
    return tbl[ifcdbh, ispc]
end

# bm/NVEL r6vol3.f — Behre total-cubic profile (VOLEQ 616BEH*** → ZONE 1). Smalian-integrated taper
# DR = HRATIO/(0.62·HRATIO+0.38). Inputs: DBHOB, DBTBH (=D·(1-bark)), FCLASS (form class), HTTOT.
# Used for the BM minor species (WP/MH/WJ/WB/LM/PY/YC/AS/CW/OS/OH) whose VEQNNC = 616BEHW.
function bm_r6vol3(dbhob::Float32, dbtbh::Float32, fclass::Int, httot::Float32, zone::Int)::Float32
    topd = 4.0f0
    d17 = Float32(fclass) / 100.0f0 * dbhob
    h17 = zone == 1 ? 17.3f0 : 33.6f0
    dbhib = dbhob - dbtbh
    (dbhib <= 0f0 || dbhib > dbhob) && (dbhib = dbhob)
    if dbhib < topd                                            # whole stem below top diam
        return 0.00272708f0 * dbhib * dbhib * httot
    elseif d17 < topd                                          # two-log: butt (H17) + upper cylinder
        v = 0.00272708f0 * (dbhib * dbhib + d17 * d17) * h17
        v += 0.00272708f0 * (d17 * d17) * (httot - h17)
        return v
    end
    a = 0.62f0; b = 1.0f0 - a
    D = zeros(Float32, 21)
    vol = 0.00272708f0 * (dbhib * dbhib + d17 * d17) * h17     # butt log stump→17.3ft
    htup = httot - h17
    D[1] = d17
    s = 0.0f0; goto100 = false; goto130 = false; iexit = 21
    i = 2
    while i <= 20
        hratio = (htup - ((i - 1) * 16.3f0)) / htup
        if hratio <= 0.0f0; iexit = i; goto100 = true; break; end
        dr = hratio / (a * hratio + b)
        D[i] = dr * D[1]
        if D[i] < topd; iexit = i; goto100 = true; break; end
        vol += 0.00272708f0 * (D[i-1] * D[i-1] + D[i] * D[i]) * 16.3f0
        if D[i] == topd; s = 16.3f0; iexit = i; goto130 = true; break; end
        i += 1
    end
    if !goto130                                                # label 100 (incl. natural loop end, iexit=21)
        goto100 || (iexit = 21)
        dr = topd / d17
        hx = (dr * b * htup) / (1.0f0 - (a * dr))
        hh = (iexit - 2) * 16.3f0
        s = htup - hx - hh
        vol += 0.00272708f0 * (D[iexit-1] * D[iexit-1] + topd * topd) * s
    end
    htup2 = httot - (16.3f0 * (iexit - 2) + h17) - s           # label 130: top-of-tree cone
    vol += 0.00272708f0 * (topd * topd) * htup2
    return vol
end
