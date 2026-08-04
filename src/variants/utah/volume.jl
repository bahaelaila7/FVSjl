# =============================================================================
# volume.jl (utah) — UT volume (ut VEQNNC). Chunk 8. Three families, all reused:
#   400/402MATW → R4 Matney taper (r4vol_volumes/r4vol_board, same as TT)
#   407FW2W     → Flewelling FW2 (cr_fw2_vol, same as KT/EM; BS/ES only)
#   400/300DVEW → Chojnacky INT-339 D2H woodland (r4d2h_vol1, same as TT; PJ species)
# Merch standards (ut/grinit.f): TOPD=6, DBHMIN=8 (sp7 LP=7), BFTOPD=6, BFMIND=8 (sp7=7).
# =============================================================================

# R3D2HV (volume/NVEL/r3d2hv.f): Region-3 woodland D2H. UT routes geocode-300 DVEW here (dvest.f:71).
# UT's only 300DVEW species is GO (300DVEW800 → oak). FCLASS=FRMCLS default 80 ≠ 1 ⇒ MSTEM/"multistem"
# branch; PROD='02' (fvsvol.f:184) ⇒ UNT=3 ⇒ VOL(1)=VOL(4)=GCUFT4 (entire = merch cubic). DBH<1 ⇒ 0;
# DBH≤3 ⇒ 0.1 floor (r3d2hv.f:441-442). Returns VOL(1). Bit-exact match to live GO 0.1-floor saplings.
@inline function r3d2hv_vol1(eq::AbstractString, d::Float32, h::Float32)::Float32
    d < 1f0 && return 0f0
    code = strip(eq)[8:10]
    d2h = d * d * h
    if code == "800"          # Oaks (INT-391 Juniper/Pinyon/Oak/Mesquite; VOLEQ(2:3)="00")
        d <= 3f0 && return 0.1f0
        d2ha = d2h / 1000f0
        g = d2ha <= 4f0 ? -0.028f0 + 1.9545f0 * d2ha + 0.1400f0 * d2ha * d2ha :
                           6.691f0 + 1.9545f0 * d2ha - 17.918f0 / d2ha
        return max(g, 0f0)
    end
    return 0f0
end

function compute_volumes_ut!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq; c = s.control
    utmerch = (stmp = c.sp_stump_ht, topd = c.sp_top_diam, scfstmp = c.sp_scf_stump,
               scftop = c.sp_scf_topd, bftopd = c.sp_bf_topd, bfstmp = c.sp_bf_stump)
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
        # Top-killed trees: full cubic (VMAX) uses the NORMAL height (norm_ht), then r4_topkill trims (see TT).
        hv = (t.trunc[i] > 0 && t.norm_ht[i] > 0) ? Float32(t.norm_ht[i]) / 100f0 : h
        if mdl == "MAT"
            mtopp = 6.0f0 * bark                             # TOPD=6 outside-bark → inside-bark top
            tcf, mcf = r4vol_volumes(eq, d, hv, mtopp, 0f0)
            mcf = d >= dbhmin ? max(mcf, 0f0) : 0f0
            bf = d >= dbhmin ? r4vol_board(eq, d, hv, mtopp, 0f0) : 0f0
            tcf, mcf, bf = r4_topkill(t, i, sp, d, hv, bark, max(tcf, 0f0), mcf, bf, utmerch)
            t.cuft_vol[i] = max(tcf, 0f0); t.merch_cuft_vol[i] = max(mcf, 0f0)
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = max(bf, 0f0)
        elseif mdl == "FW2"
            v = cr_fw2_vol(eq, d, hv; bark = bark, topd = 6.0f0, bftopd = 6.0f0, stump = 1f0, iregn = 4)
            tcf = max(v[1], 0f0)
            mcf = d >= dbhmin ? max(v[4] + v[7], 0f0) : 0f0
            bf  = d >= dbhmin ? max(v[2], 0f0) : 0f0
            tcf, mcf, bf = r4_topkill(t, i, sp, d, hv, bark, tcf, mcf, bf, utmerch)
            t.cuft_vol[i] = max(tcf, 0f0); t.merch_cuft_vol[i] = max(mcf, 0f0)
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = max(bf, 0f0)
        else                                                 # DVE woodland — dvest.f geocode dispatch
            # geocode '3' → Region-3 R3D2HV (GO); '4' → Region-4 R4D2H (PJ). (dvest.f:71/86)
            vol1 = se[1] == '3' ? r3d2hv_vol1(eq, d, h) : r4d2h_vol1(eq, d, h)
            t.cuft_vol[i] = max(vol1, 0f0)
            t.merch_cuft_vol[i] = d >= dbhmin ? max(vol1, 0f0) : 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0
        end
    end
    return s
end
