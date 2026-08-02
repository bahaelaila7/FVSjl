# =============================================================================
# volume.jl (teton) — TT volume (chunk 8). VOLEQ = 400MATW<code> (Region-4 Matney = r4vol) for most
# conifers + 400DVEW/401DVEW (DVE/Gevorkiantz = cr_dve_vol) for PM/RM/MC/OH/UJ. Cubic (total CF0 +
# merch CFGRS) via r4vol_volumes — validated 499/499 bit-exact. Board-foot (Scribner) = TODO.
# =============================================================================

# TT VOLEQ per species (ttt01.out; VOLCODE ≠ FIA — e.g. WB uses 108). geocode 400/401.
const TT_VOL_EQ = String[
    "400MATW108", "400MATW108", "400MATW202", "400DVEW133", "400MATW093", "400MATW746",
    "400MATW108", "400MATW093", "400MATW019", "400MATW122", "401DVEW065", "400DVEW066",
    "400MATW108", "400MATW108", "400MATW108", "400DVEW475", "400MATW108", "400DVEW998"]

function compute_volumes_tt!(s::StandState)
    t = s.trees; veq = s.species.vol_eq
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = veq[sp]; mdl = length(strip(eq)) >= 7 ? strip(eq)[4:6] : "   "
        if mdl == "MAT"
            # R4 merch cubic (tt/grinit.f): TOPD=6" outside-bark → inside-bark top = 6·bark; DBHMIN=8 (sp7=7).
            mtopp = 6.0f0 * tt_bratio(sp, d)
            dbhmin = sp == 7 ? 7.0f0 : 8.0f0
            tcf, mcf = r4vol_volumes(eq, d, h, mtopp, 0f0)   # (total CF0, merch CFGRS) — bit-exact
            t.cuft_vol[i] = max(tcf, 0f0)
            t.merch_cuft_vol[i] = d >= dbhmin ? max(mcf, 0f0) : 0f0
            bfmind = sp == 7 ? 7.0f0 : 8.0f0                 # board DBHMIN (BFMIND, tt/grinit.f)
            bf = d >= bfmind ? r4vol_board(eq, d, h, mtopp, 0f0) : 0f0   # BFGRS Scribner (M=1)
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = max(bf, 0f0)
        else
            # DVEW (PM/UJ/RM/MC/OH) — TT woodland volume. NOT CR's region-2 Chojnacky (that gave ~5× high:
            # jl UJ 2090 1252 vs live 258) ⇒ TT DVEW uses a different model/coefs (or DRC). DEFERRED — 0 for now.
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0
        end
    end
    return s
end
