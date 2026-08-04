# =============================================================================
# volume.jl (centralidaho) — CI volume (ci VEQNNC). Chunk 8. Same three families as UT:
#   400MATW → R4 Matney taper (r4vol); I15FW2W → Flewelling FW2 (cr_fw2_vol); 400DVEW → r4d2h woodland.
# VEQNNC dumped from live FVSci cit01.out. Merch (ci/grinit.f): TOPD=6, DBHMIN=8 (sp7 LP=7) — = UT.
# CI bark = ci_bratio (POWER model), so this mirrors compute_volumes_ut! but swaps the bark call.
# =============================================================================

# CI volume equations are FOREST-dependent (NVEL voleqdef.f R4_EQN, dispatched by FORNUM=forest%100). The base
# table below is the INGY-region assignment (region-4 forests FORNUM 2/12/13 = 402/412/413, incl. cit01=412):
# for the 4 mixed-conifer species GF/ES/PP/DF it uses the Flewelling I15FW2W* profile. The OTHER region-4 CI
# forests — FORNUM 6/14 = 406/414 (Salmon-Challis / Sawtooth) — use the region-4 Matney 400MATW* equation for
# those same 4 species; every other species is identical across the region-4 forests. jl previously hardcoded
# ONLY the INGY table (from cit01), so DF/GF/ES/PP volume on 406/414 stands was wrong (FW2 merch ~11% / board
# ~19% low vs Matney) — confirmed bit-exact: DBH34.1 DF live 400MATW202 gives merch 129.4/bdft 740 = live.
const CI_VOL_EQ = String[
    "400MATW117", "400MATW073", "I15FW2W202", "I15FW2W017", "400MATW015",   # WP WL DF GF WH
    "400MATW081", "400MATW108", "I15FW2W093", "400MATW019", "I15FW2W122",   # RC LP ES AF PP
    "400MATW108", "400DVEW998", "400MATW746", "400DVEW064", "400DVEW475",   # WB PY AS WJ MC
    "400MATW108", "400DVEW998", "400MATW108", "400MATW108"]                 # LM CW OS OH
# Matney override for FORNUM 6/14 (forest_idx 3=406, 6=414): GF(4) ES(8) PP(10) DF(3) → region-4 Matney.
const CI_VOL_EQ_MATNEY = let v = copy(CI_VOL_EQ)
    v[3] = "400MATW202"; v[4] = "400MATW015"; v[8] = "400MATW093"; v[10] = "400MATW122"; v
end
# forest_idx (ci_forkod! order [117,402,406,412,413,414]) → equation table. 3=406, 6=414 are the Matney forests.
@inline ci_vol_eq_table(forest_idx::Integer) = (forest_idx == 3 || forest_idx == 6) ? CI_VOL_EQ_MATNEY : CI_VOL_EQ

function compute_volumes_ci!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq; sd = s.coef.species; c = s.control
    cimerch = (stmp = c.sp_stump_ht, topd = c.sp_top_diam, scfstmp = c.sp_scf_stump,
               scftop = c.sp_scf_topd, bftopd = c.sp_bf_topd, bfstmp = c.sp_bf_stump)
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = veq[sp]; se = strip(eq); mdl = length(se) >= 7 ? se[4:6] : "   "
        bark = ci_bratio(sd, sp, d)
        dbhmin = sp == 7 ? 7.0f0 : 8.0f0
        # Top-killed trees: full cubic (VMAX) uses the NORMAL height (norm_ht), then r4_topkill trims (see TT).
        hv = (t.trunc[i] > 0 && t.norm_ht[i] > 0) ? Float32(t.norm_ht[i]) / 100f0 : h
        if mdl == "MAT"
            mtopp = 6.0f0 * bark
            tcf, mcf = r4vol_volumes(eq, d, hv, mtopp, 0f0)
            mcf = d >= dbhmin ? max(mcf, 0f0) : 0f0
            bf = d >= dbhmin ? r4vol_board(eq, d, hv, mtopp, 0f0) : 0f0
            tcf, mcf, bf = r4_topkill(t, i, sp, d, hv, bark, max(tcf, 0f0), mcf, bf, cimerch)
            t.cuft_vol[i] = max(tcf, 0f0); t.merch_cuft_vol[i] = max(mcf, 0f0)
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = max(bf, 0f0)
        elseif mdl == "FW2"
            v = cr_fw2_vol(eq, d, hv; bark = bark, topd = 6.0f0, bftopd = 6.0f0, stump = 1f0, iregn = 4)
            tcf = max(v[1], 0f0)
            mcf = d >= dbhmin ? max(v[4] + v[7], 0f0) : 0f0
            bf  = d >= dbhmin ? max(v[2], 0f0) : 0f0
            tcf, mcf, bf = r4_topkill(t, i, sp, d, hv, bark, tcf, mcf, bf, cimerch)
            t.cuft_vol[i] = max(tcf, 0f0); t.merch_cuft_vol[i] = max(mcf, 0f0)
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = max(bf, 0f0)
        else                                                 # DVE woodland (r4d2h, region 4): NO CFTOPK trim
            vol1 = r4d2h_vol1(eq, d, h)
            t.cuft_vol[i] = max(vol1, 0f0)
            t.merch_cuft_vol[i] = d >= dbhmin ? max(vol1, 0f0) : 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0
        end
    end
    return s
end
