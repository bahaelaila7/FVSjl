# =============================================================================
# volume.jl (eastcascades) — EC volume (R6 NVEL). Chunk 8.
#
# EC VOLEQ = voleqdef.f R6_EQN's eastside forest table (_ec_r6_eqn): INGY Flewelling with the forest's Inland
# geosub (e.g. Okanogan DF→I12FW2W202, GF→I11FW2W017), the westside F0xFW2W for DF on Mt Hood, else region-6
# Behre 616BEHW<fia>. The INGY path reuses the shared cr_fw2_vol (SHP_C2,
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

# voleqdef.f R6_EQN — EC is an EASTSIDE variant (VAR not in the westside list), so it takes the same eastside
# forest table as BM (_bm_r6_eqn); ec/sitset.f passes DIST='  ' ⇒ DISTNUM 0. jl formerly knew only Okanogan/
# Wenatchee and sent every other forest to Behre — Mt Hood (606) runs DF on westside F05FW2W202 and SF/GF/NF/ES/
# LP/PP/WH on I11-I13 INGY (FVSec_g16 NVEL table), so its cycle-0 TCuFt was +2.4%.
_ec_r6_eqn(fornum::Int, fia::Int)::String = _bm_r6_eqn(fornum, 0, fia)

# EC Behre per-tree volume — reuse the BM R6 machinery + EC form class.
function ec_behre_vol(sp::Int, ifor::Int, d::Float32, h::Float32, bark::Float32)
    fclass = ec_formcl(sp, ifor, d)
    dbtbh = d * (1f0 - bark); dbhib = d - dbtbh
    vol2 = 0f0; vol4 = 0f0
    v1 = if h <= 17.3f0
        0.00272708f0 * (dbhib * dbhib) * h
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
    t = s.trees; veq = s.species.vol_eq; sd = s.coef.species; c = s.control
    ifor = Int(s.plot.forest_idx)
    ecmerch = (stmp = c.sp_stump_ht, topd = c.sp_top_diam, scfstmp = c.sp_scf_stump,
               scftop = c.sp_scf_topd, bftopd = c.sp_bf_topd, bfstmp = c.sp_bf_stump)
    # vols.f:86-90 zeroes HT2TD for every record; NATCRS (fvsvol.f:337-339 / 484-487) stores the NVEL HT1PRD.
    fill!(t.merch_top_cf, 0f0); fill!(t.merch_top_bf, 0f0)
    htb = zeros(Float32, 2)
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        htb[1] = 0f0; htb[2] = 0f0
        if d < 1f0 || sp < 1 || sp > 32
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = veq[sp]; se = strip(eq); mdl = length(se) >= 7 ? se[4:6] : "   "
        # vols.f:132,150-151: BARK=BRATIO(ISPC,DBH_start,H) before `D=D+DG(I)/BARK` ⇒ projected cycles use the stashed
        # start-of-cycle bark (t.vol_bark) for the merch tops / DBTBH / CFTOPK; grown-DBH bark at cycle 0 / dead records.
        bark = (i <= t.n && t.vol_bark[i] > 0f0) ? t.vol_bark[i] : wc_bratio(sd, sp, d)
        dbhmin = sp == 7 ? 6.0f0 : 7.0f0                # ec/sitset.f westside: LP(sp7)=6, else 7
        hv = (t.trunc[i] > 0 && t.norm_ht[i] > 0) ? Float32(t.norm_ht[i]) / 100f0 : h
        local tcf::Float32, mcf::Float32, bf::Float32
        if mdl == "FW2" && (se[1] == 'F' || se[1] == 'f')
            tcf, mcf, bf = wc_fw2_westside_vol(eq, d, hv, bark; ht2td = htb)   # F-prefix westside SHP (not used on forest 8)
        elseif mdl == "FW2"
            v = cr_fw2_vol(eq, d, hv; bark = bark, topd = 4.5f0, bftopd = 4.5f0, stump = 1f0,
                           iregn = 6, board_cor = 'N', merch_opt = 23, sf_hs = true, ht2td = htb)   # MRULES REGN 6
            tcf = max(v[1], 0f0); mcf = max(v[4] + v[7], 0f0); bf = max(v[2], 0f0)
        else                                            # 616BEHW
            tcf, mcf, bf = ec_behre_vol(sp, ifor, d, hv, bark)
            htb[1] = r6vol_ht1prd(d, ec_formcl(sp, ifor, d), 4.5f0 * bark, hv, d - d * (1f0 - bark)); htb[2] = htb[1]
        end
        # fvsvol.f:337-339 / 484-487 (BFMIND = DBHMIN, ec/sitset.f)
        d >= dbhmin && (t.merch_top_cf[i] = htb[1]; t.merch_top_bf[i] = htb[2])
        # Broken/killed-top reduction (vols.f:193 CFTOPK/BFTOPK) — for EVERY equation, not just INGY: the volume
        # above used the NORMAL height (norm_ht/hv, cratet.f), so trim it back to the standing break via the Behre
        # taper. EC merch TOPD=4.5 (== BM). ect01 (Okanogan): rec-6 WL / rec-22 DF broken tops, TCuFt 1617→1615.
        # Mt Hood (606) runs WL on Behre and DF on westside F05FW2W, which the INGY-only trim skipped.
        tcf, mcf, bf = r4_topkill(t, i, sp, d, hv, bark, tcf, mcf, bf, ecmerch, _BM_TOPD45)
        t.cuft_vol[i] = max(tcf, 0f0)
        t.merch_cuft_vol[i] = d >= dbhmin ? max(mcf, 0f0) : 0f0
        t.saw_cuft_vol[i] = 0f0
        t.bdft_vol[i] = d >= dbhmin ? max(bf, 0f0) : 0f0
    end
    return s
end
