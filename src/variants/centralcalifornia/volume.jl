# =============================================================================
# volume.jl (centralcalifornia) — CA volume (R6 NVEL). Chunk 8.
#
# CA VOLEQDEF (VAR='CA', IREGN=KODFOR/100, FORST=INTFOR; ca/sitset.f:251-290) resolves per-species NVEL
# equations, dumped from FVSca_clean cat01.out "NATIONAL VOLUME ESTIMATOR LIBRARY EQUATION NUMBERS":
#   • WF(fia15) → I00FW2W093, PP(fia122) → I00FW2W073   — INGY Flewelling (reuse cr_fw2_vol, iregn=6)
#   • DF(fia202) → F06FW2W202, WH(fia263) → F06FW2W263  — westside Flewelling SHP (reuse wc_fw2_westside_vol)
#   • everything else → Behre 616BEHW<fia>             — reuse bm_r6vol3/r6dibs/r6vol1 + CA form class (ca/formcl.f)
# Merch (ca/grinit.f): DBHMIN=7 (sp-index 11 = 6); TOPD=BFTOPD=4.5 (R6, ca/sitset.f IFOR 6-10).
# =============================================================================

# ca/formcl.f form-class tables (data/centralcalifornia/formcl_ca.csv): ROGRFC/SISKFC = [ISPC,DBHclass 1-5];
# BLM710/711/712 = [ISPC] (no DBH class). IFOR 6=Rogue,7=Siskiyou,8=Rosenburg,9=Medford,10=CoosBay; else 80.
const CA_FORMCL = let
    rogr = zeros(Int, 50, 5); sisk = zeros(Int, 50, 5)
    blm = Dict("BLM710" => zeros(Int, 50), "BLM711" => zeros(Int, 50), "BLM712" => zeros(Int, 50))
    for l in readlines(joinpath(CA_DATADIR, "formcl_ca.csv"))[2:end]
        isempty(strip(l)) && continue
        f = split(strip(l), ','); tbl = f[1]; sp = parse(Int, f[2]); dc = parse(Int, f[3]); v = parse(Int, f[4])
        if tbl == "ROGRFC"; rogr[sp, dc] = v
        elseif tbl == "SISKFC"; sisk[sp, dc] = v
        else; blm[tbl][sp] = v; end
    end
    (rogr = rogr, sisk = sisk, blm710 = blm["BLM710"], blm711 = blm["BLM711"], blm712 = blm["BLM712"])
end

# ca/formcl.f: form class FC for species `sp`, forest index `ifor`, DBH `d` (FRMCLS keyword override not modeled).
function ca_formcl(sp::Integer, ifor::Int, d::Real)::Int
    (sp < 1 || sp > 50) && return 80
    ifcdbh = Int(floor((Float32(d) - 1f0) / 10f0 + 1f0))
    ifcdbh < 1 && (ifcdbh = 1); Float32(d) > 40.9f0 && (ifcdbh = 5)
    if ifor == 6; return CA_FORMCL.rogr[sp, ifcdbh]
    elseif ifor == 7; return CA_FORMCL.sisk[sp, ifcdbh]
    elseif ifor == 8; return CA_FORMCL.blm710[sp]
    elseif ifor == 9; return CA_FORMCL.blm711[sp]
    elseif ifor == 10; return CA_FORMCL.blm712[sp]
    end
    return 80
end

# ca/sitset.f:252-278 VOLEQDEF (VAR='CA', IREGN=KODFOR/100) for the Region-6 forests (IFOR 6,7 = 610 Rogue River /
# 611 Siskiyou — identical tables), read off FVSca_g16's "NATIONAL VOLUME ESTIMATOR LIBRARY EQUATION NUMBERS"
# (STDINFO 610/611; cubic == board): westside Flewelling DF/WH (F06), INGY WF/PP (I00FW2W093/073), Behre otherwise.
const CA_R6_VOL_EQ = (
    "616BEHW000", "616BEHW081", "616BEHW242", "I00FW2W093", "616BEHW020", "616BEHW021", "F06FW2W202", "F06FW2W263",
    "616BEHW264", "616BEHW101", "616BEHW103", "616BEHW108", "616BEHW000", "616BEHW113", "616BEHW116", "616BEHW117",
    "616BEHW119", "I00FW2W073", "616BEHW000", "616BEHW000", "616BEHW064", "616BEHW000", "616BEHW000", "616BEHW231",
    "616BEHW299", "616BEHW000", "616BEHW000", "616BEHW000", "616BEHW000", "616BEHW815", "616BEHW818", "616BEHW000",
    "616BEHW000", "616BEHW312", "616BEHW000", "616BEHW351", "616BEHW361", "616BEHW431", "616BEHW492", "616BEHW000",
    "616BEHW000", "616BEHW631", "616BEHW000", "616BEHW746", "616BEHW747", "616BEHW920", "616BEHW000", "616BEHW000",
    "616BEHW998", "616BEHW211")
# … and the BLM forests (IFOR 8,9 = 710 Roseburg / 711 Medford; IFOR 10 = 712 Coos Bay): voleqdef R7_EQN ⇒
# B00BEHW<fia> / DF B01BEHW202 (710/711) or B02BEHW202 (712) ⇒ NVEL BLMVOL (volinit.f routes VOLEQ(1:1)='B').
# jl used to send these through the R6 Behre form-class path (every 710-712 tree wrong, measured vs FVSca_g16).
const CA_BLM_VOL_EQ = (
    "B00BEHW081", "B00BEHW081", "B00BEHW242", "B00BEHW015", "B00BEHW021", "B00BEHW021", "B01BEHW202", "B00BEHW263",
    "B00BEHW260", "B00BEHW119", "B00BEHW108", "B00BEHW108", "B00BEHW108", "B00BEHW108", "B00BEHW116", "B00BEHW117",
    "B00BEHW119", "B00BEHW122", "B00BEHW108", "B00BEHW108", "B00BEHW242", "B00BEHW093", "B00BEHW211", "B00BEHW231",
    "B00BEHW999", "B00BEHW800", "B00BEHW800", "B00BEHW800", "B00BEHW800", "B00BEHW800", "B00BEHW800", "B00BEHW800",
    "B00BEHW800", "B00BEHW312", "B00BEHW800", "B00BEHW351", "B00BEHW361", "B00BEHW431", "B00BEHW999", "B00BEHW312",
    "B00BEHW999", "B00BEHW631", "B00BEHW800", "B00BEHW999", "B00BEHW747", "B00BEHW999", "B00BEHW231", "B00BEHW631",
    "B00BEHW999", "B00BEHW211")

"VOLEQ for CA species `sp` (1..50) on forest index `ifor` (ca/forkod.f JFOR order 505,506,508,511,514,610,611,710,
711,712,518)."
function ca_voleq(ifor::Int, sp::Int)::String
    (ifor == 6 || ifor == 7) && return CA_R6_VOL_EQ[sp]
    8 <= ifor <= 10 && return (ifor == 10 && sp == 7) ? "B02BEHW202" : CA_BLM_VOL_EQ[sp]
    return CA_R5_VOL_EQ[sp]                           # R5 forests IFOR 1-5 and 11 (518 → R5_EQN)
end

# CA Behre per-tree volume — reuse the BM R6 machinery + CA form class (mirrors ec_behre_vol).
function ca_behre_vol(sp::Int, ifor::Int, d::Float32, h::Float32, bark::Float32)
    fclass = ca_formcl(sp, ifor, d)
    dbtbh = d * (1f0 - bark); dbhib = d - dbtbh
    vol2 = 0f0; vol4 = 0f0
    v1 = if h <= 17.3f0
        0.00272708f0 * dbhib * dbhib * h
    else
        v = bm_r6vol3(d, dbtbh, fclass, h, 1)
        mtopp = 4.5f0 * bark
        xlogs, ld1 = bm_r6dibs(d, fclass, mtopp, h)
        lv1, lv4 = bm_r6vol1(d, fclass, xlogs, ld1)
        nlog = Int(floor(xlogs)); nacc = (xlogs - nlog) > 0f0 ? nlog + 1 : nlog
        for k in 1:nacc
            vol2 += bm_anint(lv1[k]); vol4 += bm_anint(lv4[k] * 10f0) / 10f0
        end
        v
    end
    return (max(v1, 0f0), max(vol4, 0f0), max(vol2, 0f0))
end

# ca/sitset.f VOLEQDEF for the REGION-5 forests (IFOR 1-5: 505/506/508/511/514; 518 → 514): the R5_EQN FIA
# table ⇒ Wensel-Krumland R5TAP profile "500WO2W" (conifers + redwood) or the R5HARV California-hardwood DVE
# "500DVEW" — the SAME NVEL kernels NC's R5 forests use (nc_wo2w_vol / nc_r5harv_vol). Per CA species 1..50,
# read off the live FVSca_g16 "NATIONAL VOLUME ESTIMATOR LIBRARY EQUATION NUMBERS" table (STDINFO 505).
# jl used to run the R6 Behre/FW2 path on R5 stands (+10% TCuFt, +56% BdFt at cycle 0, measured).
const CA_R5_VOL_EQ = String[
    "500WO2W081", "500WO2W081", "500WO2W081", "500WO2W015", "500WO2W020", "500WO2W020", "500WO2W202", "500WO2W015",
    "500WO2W015", "500WO2W108", "500WO2W108", "500WO2W108", "500WO2W108", "500WO2W108", "500WO2W116", "500WO2W117",
    "500WO2W117", "500WO2W122", "500WO2W108", "500WO2W108", "500DVEW060", "500WO2W015", "500DVEW212", "500WO2W108",
    "500WO2W108", "500DVEW801", "500DVEW805", "500DVEW807", "500DVEW811", "500DVEW815", "500DVEW818", "500DVEW821",
    "500DVEW839", "500DVEW312", "500DVEW807", "500DVEW351", "500DVEW361", "500DVEW431", "500DVEW807", "500DVEW807",
    "500DVEW818", "500DVEW631", "500DVEW818", "500DVEW818", "500DVEW818", "500DVEW807", "500DVEW807", "500DVEW981",
    "500DVEW801", "500WO2W211"]

function compute_volumes_ca!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; sd = s.coef.species
    ifor = Int(s.plot.forest_idx)
    c = s.control
    merch = (stmp = c.sp_stump_ht, topd = c.sp_top_diam, scfstmp = c.sp_scf_stump,
             scftop = c.sp_scf_topd, bftopd = c.sp_bf_topd, bfstmp = c.sp_bf_stump)
    # vols.f:86-90 zeroes HT2TD for every record; NATCRS (fvsvol.f:337-339 / 484-487) stores the NVEL HT1PRD.
    fill!(t.merch_top_cf, 0f0); fill!(t.merch_top_bf, 0f0)
    htb = zeros(Float32, 2)
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        htb[1] = 0f0; htb[2] = 0f0
        if d < 1f0 || sp < 1 || sp > 50
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = ca_voleq(ifor, sp)
        se = strip(eq); mdl = length(se) >= 7 ? se[4:6] : "   "
        # vols.f:150 BARK=BRATIO(ISPC,D,H) is taken at the START-of-cycle DBH (before D=D+DG/BARK) — the stashed
        # vol_bark; the grown-DBH bark only at cycle 0 / for dead records. It sets the merch/board tops (TOPD·BARK).
        bark = (i <= t.n && t.vol_bark[i] > 0f0) ? t.vol_bark[i] : wc_bratio(sd[:bark1][sp], sd[:bark2][sp], Int(sd[:bark_imap][sp]), d)
        dbhmin = sp == 11 ? 6.0f0 : 7.0f0               # ca/grinit.f: sp-index 11 = 6, else 7
        hv = (t.trunc[i] > 0 && t.norm_ht[i] > 0) ? Float32(t.norm_ht[i]) / 100f0 : h
        local tcf::Float32, mcf::Float32, bf::Float32
        if se[1] == '5'                                 # Region 5 (ca/sitset.f VOLEQDEF IREGN=5; IFOR 1-5, 11)
            tcf, mcf, bf = nvel_r5_vol(eq, d, hv, bark, c.sp_top_diam[sp], c.sp_bf_topd[sp]; ht2td = htb)
        elseif se[1] == 'B'                             # BLM 710/711/712 → NVEL BLMVOL, ca/formcl.f BLM form class
            tcf, mcf, bf = _blm_natcrs(eq, ca_formcl(sp, ifor, d), d, hv, bark, c.sp_top_diam[sp], c.sp_bf_topd[sp],
                                       c.sp_bf_dbhmin[sp]; ht2td = htb)
        elseif mdl == "FW2" && (se[1] == 'F' || se[1] == 'f')
            tcf, mcf, bf = wc_fw2_westside_vol(eq, d, hv, bark; ht2td = htb)   # F06 westside SHP (DF/WH)
        elseif mdl == "FW2"
            v = cr_fw2_vol(eq, d, hv; bark = bark, topd = 4.5f0, bftopd = 4.5f0, stump = 1f0,
                           iregn = 6, board_cor = 'N', merch_opt = 23, sf_hs = true, ht2td = htb)   # INGY (WF/PP)
            tcf = max(v[1], 0f0); mcf = max(v[4] + v[7], 0f0); bf = max(v[2], 0f0)
        else                                            # 616BEHW
            tcf, mcf, bf = ca_behre_vol(sp, ifor, d, hv, bark)
            htb[1] = r6vol_ht1prd(d, ca_formcl(sp, ifor, d), 4.5f0 * bark, hv, d - d * (1f0 - bark)); htb[2] = htb[1]
        end
        # fvsvol.f:337-339 / 484-487 (BFMIND = DBHMIN, ca/grinit.f)
        d >= dbhmin && (t.merch_top_cf[i] = htb[1]; t.merch_top_bf[i] = htb[2])
        # vols.f CFTOPK/BFTOPK broken-top trim — was missing on CA (R6 DF trc49 D15.9: jl 53.8 vs live 40.2 cuft).
        tcf, mcf, bf = r4_topkill(t, i, sp, d, hv, bark, tcf, mcf, bf, merch, (6 <= ifor <= 10) ? _BM_TOPD45 : _R4_TOPD6)
        t.cuft_vol[i] = max(tcf, 0f0)
        t.merch_cuft_vol[i] = d >= dbhmin ? max(mcf, 0f0) : 0f0
        t.saw_cuft_vol[i] = 0f0
        t.bdft_vol[i] = d >= dbhmin ? max(bf, 0f0) : 0f0
    end
    return s
end
