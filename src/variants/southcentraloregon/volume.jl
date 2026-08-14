# =============================================================================
# volume.jl (southcentraloregon) — SO volume (R6 NVEL). Chunk 8.
#
# SO VOLEQDEF (dumped from FVSso_g16 so_c1.out "NATIONAL VOLUME ESTIMATOR LIBRARY EQUATION NUMBERS", forest 601):
#   • FW2 (all I-prefix INGY, NO westside-SHP F-prefix unlike CA): SP I13FW2W122 / DF I11FW2W202 / WF·SH·GF·SF
#     I00FW2W017 / IC I11FW2W242 / LP I11FW2W108 / PP I11FW2W122 / WL I11FW2W073  → reuse cr_fw2_vol (iregn=6).
#   • everything else → Behre 616BEHW<fia>  → reuse bm_r6vol3/r6dibs/r6vol1 + SO form class (so/formcl.f).
# Merch (so/grinit.f): DBHMIN=BFMIND=SCFMIND=9.0 ALL species (no sp-11 special); TOPD=4.5 (IFOR 1-3,10) else 6
# (so/sitset.f:167). SO bark = so_bratio (NOT wc_bratio). MEASURED vs the relinked FVSso oracle so_c1.sum.
# =============================================================================

# so/formcl.f DESHFC/FREMFC/WINEFC form-class tables [ISPC, DBHclass 1..5] (33 sp). Forest 1/10=DESH, 2=FREM, 3=WINE.
const SO_FORMCL_DESH = permutedims(Float32[
    80 78 80 80 78 65 82 82 80 76 60 80 80 82 80 80 80 69 78 56 70 76 76 77 76 76 70 77 75 77 77 80 77
    81 78 80 80 78 65 82 82 80 76 60 80 80 82 80 81 80 70 78 60 70 76 78 77 78 78 70 77 75 77 77 80 77
    81 78 78 78 78 65 82 82 78 79 60 78 78 82 78 81 80 70 78 60 70 77 78 77 78 78 70 77 75 77 77 78 77
    82 80 76 78 76 65 82 82 78 80 60 78 78 82 78 82 80 68 78 60 70 77 78 77 78 78 70 77 75 77 77 76 77
    82 80 76 76 74 65 82 82 76 80 60 76 76 82 76 82 80 68 78 60 70 76 78 76 78 78 70 76 75 76 76 76 76])
const SO_FORMCL_FREM = permutedims(Float32[
    83 70 66 68 64 64 83 77 66 70 34 76 78 82 80 83 78 69 76 56 70 76 76 77 76 76 70 77 75 77 77 66 77
    82 82 68 76 64 66 82 79 74 82 48 78 76 82 80 82 78 70 78 60 70 76 78 77 78 78 70 77 75 77 77 68 77
    80 82 68 82 66 66 80 80 80 82 48 77 74 82 78 80 80 70 80 60 70 77 78 77 78 78 70 77 75 77 77 65 77
    80 84 68 83 66 66 80 80 81 84 50 76 74 84 78 80 80 68 80 60 70 77 78 77 78 78 70 77 75 77 77 68 77
    80 78 68 80 66 66 80 81 78 80 50 76 74 84 76 80 80 68 82 60 70 76 78 77 78 78 70 77 75 77 77 68 77])
const SO_FORMCL_WINE = permutedims(Float32[
    78 72 78 76 75 62 83 80 76 73 60 76 81 82 80 81 78 69 76 56 70 76 76 77 76 76 70 77 75 77 77 78 77
    79 76 77 77 78 65 82 80 78 78 60 78 81 82 80 82 78 70 76 60 70 76 78 77 78 78 70 77 75 77 77 77 77
    80 76 75 78 79 64 80 80 78 80 60 77 81 82 78 82 80 70 76 60 70 77 78 77 78 78 70 77 75 77 77 75 77
    79 77 74 78 79 62 80 79 78 80 60 76 78 84 78 80 80 68 74 60 70 77 78 77 78 78 70 77 75 77 77 74 77
    78 77 73 76 78 62 80 78 78 80 60 76 74 84 76 80 80 68 74 60 70 76 78 77 78 78 70 77 75 77 77 73 77])

# so/formcl.f: form class FC for species `sp`, forest index `ifor`, DBH `d` (FRMCLS keyword override not modeled).
# Only IFOR≤3 or 10 use the tables; other forests default 80 (matches the .LE.3.OR.10 gate + ELSE FC=80).
@inline function so_formcl(sp::Integer, ifor::Int, d::Real)::Int
    (sp < 1 || sp > 33) && return 80
    (ifor <= 3 || ifor == 10) || return 80
    ifcdbh = Int(floor((Float32(d) - 1f0) / 10f0 + 1f0))
    ifcdbh < 1 && (ifcdbh = 1); Float32(d) > 40.9f0 && (ifcdbh = 5)
    tbl = (ifor == 1 || ifor == 10) ? SO_FORMCL_DESH : ifor == 2 ? SO_FORMCL_FREM : SO_FORMCL_WINE
    return Int(tbl[sp, ifcdbh])
end

# so/sitset.f VOLEQDEF (VAR='SO', forest 601 = IFOR1), per JSP (dumped from FVSso_g16 so_c1.out).
const SO_VOLEQ = String[
    "616BEHW119", "I13FW2W122", "I11FW2W202", "I00FW2W017", "616BEHW264",
    "I11FW2W242", "I11FW2W108", "616BEHW093", "I00FW2W017", "I11FW2W122",
    "616BEHW064", "I00FW2W017", "616BEHW019", "I00FW2W017", "616BEHW022",
    "616BEHW101", "I11FW2W073", "616BEHW242", "616BEHW263", "616BEHW231",
    "616BEHW352", "616BEHW351", "616BEHW312", "616BEHW746", "616BEHW747",
    "616BEHW768", "616BEHW815", "616BEHW920", "616BEHW431", "616BEHW475",
    "616BEHW478", "616BEHW299", "616BEHW998"]

# SO Behre per-tree volume — reuse the BM R6 machinery + SO form class (mirrors ca_behre_vol / ec_behre_vol).
function so_behre_vol(sp::Int, ifor::Int, d::Float32, h::Float32, bark::Float32)
    fclass = so_formcl(sp, ifor, d)
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

function compute_volumes_so!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; sd = s.coef.species
    ifor = Int(s.plot.forest_idx)
    topd = (ifor <= 3 || ifor == 10) ? 4.5f0 : 6.0f0     # so/sitset.f:167 TOPD/BFTOPD
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0 || sp < 1 || sp > 33
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = SO_VOLEQ[sp]; se = strip(eq); mdl = length(se) >= 7 ? se[4:6] : "   "
        bark = so_bratio(sd, sp, d)
        hv = (t.trunc[i] > 0 && t.norm_ht[i] > 0) ? Float32(t.norm_ht[i]) / 100f0 : h
        local tcf::Float32, mcf::Float32, bf::Float32
        if mdl == "FW2"                                  # all I-prefix INGY (cr_fw2_vol)
            v = cr_fw2_vol(eq, d, hv; bark = bark, topd = topd, bftopd = topd, stump = 1f0,
                           iregn = 6, board_cor = 'N', merch_opt = 23)
            tcf = max(v[1], 0f0); mcf = max(v[4] + v[7], 0f0); bf = max(v[2], 0f0)
        else                                             # 616BEHW
            tcf, mcf, bf = so_behre_vol(sp, ifor, d, hv, bark)
        end
        t.cuft_vol[i] = max(tcf, 0f0)
        t.merch_cuft_vol[i] = d >= 9f0 ? max(mcf, 0f0) : 0f0   # so/grinit.f DBHMIN=9
        t.saw_cuft_vol[i] = 0f0
        t.bdft_vol[i] = d >= 9f0 ? max(bf, 0f0) : 0f0
    end
    return s
end
