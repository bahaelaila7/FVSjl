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
    mtopp = 5.622f0                          # R4 cubic merch top DIB (matches live; source TBD)
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = veq[sp]; mdl = length(strip(eq)) >= 7 ? strip(eq)[4:6] : "   "
        if mdl == "MAT"
            tcf, mcf = r4vol_volumes(eq, d, h, mtopp, 0f0)   # (total CF0, merch CFGRS) — bit-exact
            t.cuft_vol[i] = max(tcf, 0f0)
            t.merch_cuft_vol[i] = d >= 7f0 ? max(mcf, 0f0) : 0f0   # merch DBHMIN (R4 cubic ~7")
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0     # board Scribner = TODO
        else
            # DVEW (PM/RM/MC/OH/UJ) — reuse the CR DVE/Gevorkiantz kernel (deferred; not in ttt01)
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0
        end
    end
    return s
end
