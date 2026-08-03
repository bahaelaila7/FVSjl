# =============================================================================
# volume.jl (centralidaho) — CI volume (ci VEQNNC). Chunk 8. Same three families as UT:
#   400MATW → R4 Matney taper (r4vol); I15FW2W → Flewelling FW2 (cr_fw2_vol); 400DVEW → r4d2h woodland.
# VEQNNC dumped from live FVSci cit01.out. Merch (ci/grinit.f): TOPD=6, DBHMIN=8 (sp7 LP=7) — = UT.
# CI bark = ci_bratio (POWER model), so this mirrors compute_volumes_ut! but swaps the bark call.
# =============================================================================

const CI_VOL_EQ = String[
    "400MATW117", "400MATW073", "I15FW2W202", "I15FW2W017", "400MATW015",   # WP WL DF GF WH
    "400MATW081", "400MATW108", "I15FW2W093", "400MATW019", "I15FW2W122",   # RC LP ES AF PP
    "400MATW108", "400DVEW998", "400MATW746", "400DVEW064", "400DVEW475",   # WB PY AS WJ MC
    "400MATW108", "400DVEW998", "400MATW108", "400MATW108"]                 # LM CW OS OH

function compute_volumes_ci!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq; sd = s.coef.species
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = veq[sp]; se = strip(eq); mdl = length(se) >= 7 ? se[4:6] : "   "
        bark = ci_bratio(sd, sp, d)
        dbhmin = sp == 7 ? 7.0f0 : 8.0f0
        if mdl == "MAT"
            mtopp = 6.0f0 * bark
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
