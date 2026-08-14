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

# ca/sitset.f VOLEQDEF (VAR='CA'), resolved per FIA species (dumped from FVSca_clean cat01.out).
function _ca_r6_eqn(fia::Int)::String
    fia == 202 && return "F06FW2W202"                # DF → westside Flewelling
    fia == 263 && return "F06FW2W263"                # WH → westside Flewelling
    fia == 15  && return "I00FW2W093"                # WF → INGY (Engelmann-spruce eqn 093)
    fia == 122 && return "I00FW2W073"                # PP → INGY (western-larch eqn 073)
    return "616BEHW" * lpad(string(fia), 3, '0')     # region-6 Behre default
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

function compute_volumes_ca!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; sd = s.coef.species
    ifor = Int(s.plot.forest_idx)
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0 || sp < 1 || sp > 50
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = _ca_r6_eqn(parse(Int, strip(s.species.fia[sp])))
        se = strip(eq); mdl = length(se) >= 7 ? se[4:6] : "   "
        bark = wc_bratio(sd[:bark1][sp], sd[:bark2][sp], Int(sd[:bark_imap][sp]), d)
        dbhmin = sp == 11 ? 6.0f0 : 7.0f0               # ca/grinit.f: sp-index 11 = 6, else 7
        hv = (t.trunc[i] > 0 && t.norm_ht[i] > 0) ? Float32(t.norm_ht[i]) / 100f0 : h
        local tcf::Float32, mcf::Float32, bf::Float32
        if mdl == "FW2" && (se[1] == 'F' || se[1] == 'f')
            tcf, mcf, bf = wc_fw2_westside_vol(eq, d, hv, bark)   # F06 westside SHP (DF/WH)
        elseif mdl == "FW2"
            v = cr_fw2_vol(eq, d, hv; bark = bark, topd = 4.5f0, bftopd = 4.5f0, stump = 1f0,
                           iregn = 6, board_cor = 'N', merch_opt = 23)   # INGY (WF/PP)
            tcf = max(v[1], 0f0); mcf = max(v[4] + v[7], 0f0); bf = max(v[2], 0f0)
        else                                            # 616BEHW
            tcf, mcf, bf = ca_behre_vol(sp, ifor, d, hv, bark)
        end
        t.cuft_vol[i] = max(tcf, 0f0)
        t.merch_cuft_vol[i] = d >= dbhmin ? max(mcf, 0f0) : 0f0
        t.saw_cuft_vol[i] = 0f0
        t.bdft_vol[i] = d >= dbhmin ? max(bf, 0f0) : 0f0
    end
    return s
end
