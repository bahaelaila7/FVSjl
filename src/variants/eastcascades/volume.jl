# =============================================================================
# volume.jl (eastcascades) — EC volume (R6 NVEL). Chunk 8.
#
# EC forest-608 (Okanogan, FORNUM 8) VOLEQ (voleqdef.f R6_EQN + dumped from FVSec_clean): the FW2W species
# route to INGY Flewelling with the Inland geosubs I11/I12 (NOT the F-prefix westside SHP nor the I00 base
# WC used) — DF→I12FW2W202, LP→I12FW2W108, ES→I11FW2W093, WL/PP(fia 73/122)→I12FW2W122, GF(fia17)→
# I11FW2W017 — everything else → Behre 616BEHW<fia>. The INGY path reuses the shared cr_fw2_vol (SHP_C2,
# iregn=6); Behre reuses the BM R6 machinery + EC form class (formcl_ec.csv). Merch (ec/sitset.f westside):
# LP(sp7)=6" DBHmin, else 7"; TOPD=BFTOPD=4.5.
# =============================================================================

const EC_FORMCL = let
    A = zeros(Int, 10, 32, 5)
    for l in readlines(joinpath(EC_DATADIR, "formcl_ec.csv"))[2:end]
        f = split(strip(l), ','); (isempty(f) || isempty(f[1])) && continue
        fi = parse(Int, f[1]); sp = parse(Int, f[2]); dc = parse(Int, f[3]); v = parse(Int, f[4])
        (1 <= fi <= 10 && 1 <= sp <= 32 && 1 <= dc <= 5) && (A[fi, sp, dc] = v)
    end
    A
end
function ec_formcl(sp::Integer, ifor::Int, d::Real)::Int
    (sp < 1 || sp > 32) && return 80
    (ifor < 1 || ifor > 10) && (ifor = 2)
    ifcdbh = Int(floor((Float32(d) - 1f0) / 10f0 + 1f0))
    ifcdbh < 1 && (ifcdbh = 1); Float32(d) > 40.9f0 && (ifcdbh = 5)
    v = EC_FORMCL[ifor, sp, ifcdbh]; v == 0 ? 80 : v
end

# voleqdef.f R6_EQN, FORNUM=8 (Okanogan) + DISTNUM=0. INGY overrides by species FIA; else region-6 Behre.
function _ec_r6_eqn(fornum::Int, fia::Int)::String
    if fornum == 8 || fornum == 17
        fia == 202 && return "I12FW2W202"           # DF
        fia == 108 && return "I12FW2W108"            # LP
        fia == 93  && return "I11FW2W093"            # ES
        (fia == 122 || fia == 73) && return "I12FW2W122"   # PP / WL
        fia == 17  && return "I11FW2W017"            # GF
    end
    return "616BEHW" * lpad(string(fia), 3, '0')     # region-6 Behre default
end

# EC Behre per-tree volume — reuse the BM R6 machinery + EC form class.
function ec_behre_vol(sp::Int, ifor::Int, d::Float32, h::Float32, bark::Float32)
    fclass = ec_formcl(sp, ifor, d)
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

function compute_volumes_ec!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq; sd = s.coef.species
    ifor = Int(s.plot.forest_idx)
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0 || sp < 1 || sp > 32
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = veq[sp]; se = strip(eq); mdl = length(se) >= 7 ? se[4:6] : "   "
        bark = wc_bratio(sd, sp, d)
        dbhmin = sp == 7 ? 6.0f0 : 7.0f0                # ec/sitset.f westside: LP(sp7)=6, else 7
        hv = (t.trunc[i] > 0 && t.norm_ht[i] > 0) ? Float32(t.norm_ht[i]) / 100f0 : h
        local tcf::Float32, mcf::Float32, bf::Float32
        if mdl == "FW2" && (se[1] == 'F' || se[1] == 'f')
            tcf, mcf, bf = wc_fw2_westside_vol(eq, d, hv, bark)   # F-prefix westside SHP (not used on forest 8)
        elseif mdl == "FW2"
            v = cr_fw2_vol(eq, d, hv; bark = bark, topd = 4.5f0, bftopd = 4.5f0, stump = 1f0,
                           iregn = 6, board_cor = 'N', merch_opt = 23)   # MRULES REGN 6: COR='N', OPT=23
            tcf = max(v[1], 0f0); mcf = max(v[4] + v[7], 0f0); bf = max(v[2], 0f0)
        else                                            # 616BEHW
            tcf, mcf, bf = ec_behre_vol(sp, ifor, d, hv, bark)
        end
        t.cuft_vol[i] = max(tcf, 0f0)
        t.merch_cuft_vol[i] = d >= dbhmin ? max(mcf, 0f0) : 0f0
        t.saw_cuft_vol[i] = 0f0
        t.bdft_vol[i] = d >= dbhmin ? max(bf, 0f0) : 0f0
    end
    return s
end
