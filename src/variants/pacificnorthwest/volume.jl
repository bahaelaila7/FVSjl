# =============================================================================
# volume.jl (pacificnorthwest) — PN volume (R6 NVEL). Chunk 8.
#
# PN forest-612 (Siuslaw) VOLEQ (dumped from FVSpn_clean): DF → F00FW2W202 (westside Flewelling, geosub
# '00' = base, no regional override); WH → F03FW2W263; everything else → Behre 616BEHW<fia> (PN does NOT
# use the INGY overrides WC had for GF/NF/IC; RA → NVBM240351, not reachable on pnt01). Reuses the WC
# westside SHP_W3/W4/W5 (wc_fw2_westside_vol handles geosub '00') + BM's R6 Behre machinery + PN form class
# (formcl_pn.csv). Merch (pn/sitset.f non-BLM) = the WC westside defaults (TOPD=BFTOPD=4.5, DBHMIN=7, LP=6).
# =============================================================================

const PN_FORMCL = let
    T = Dict{Tuple{Int,Int,Int},Int}()
    for l in readlines(joinpath(PN_DATADIR, "formcl_pn.csv"))[2:end]
        f = split(strip(l), ','); (isempty(f) || isempty(f[1])) && continue
        T[(parse(Int, f[1]), parse(Int, f[2]), parse(Int, f[3]))] = parse(Int, f[4])
    end
    A = zeros(Int, 10, 39, 5)
    for ((fi, sp, dc), v) in T
        if fi <= 6; A[fi, sp, dc] = v else; for c in 1:5; A[fi, sp, c] = v; end; end
    end
    A
end
function pn_formcl(sp::Integer, ifor::Int, d::Real)::Int
    (sp < 1 || sp > 39) && return 80
    (ifor < 1 || ifor > 10) && (ifor = 2)
    ifcdbh = Int(floor((Float32(d) - 1f0) / 10f0 + 1f0))
    ifcdbh < 1 && (ifcdbh = 1); Float32(d) > 40.9f0 && (ifcdbh = 5)
    v = PN_FORMCL[ifor, sp, ifcdbh]; v == 0 ? 80 : v
end

# --- PN VOLEQDEF (pn/sitset.f:213-243 → voleqdef.f R6_EQN westside for REGN 6 / R7_EQN BLM for REGN 7; the
#     Quinault reservation 800 is volumed as Olympic 609, pn/sitset.f:224-230). The complete PN input domain is the
#     6 JFOR forests (pn/forkod.f:63) × 39 species, so this is the full decision table, read off FVSpn_g16's
#     "NATIONAL VOLUME ESTIMATOR LIBRARY EQUATION NUMBERS" block per forest (cubic == board everywhere). jl had
#     only Siuslaw (612): 609/800 put SS on F03FW2W263 and DF on F03, WH on F00; BLM 708/709/712 are B00BEHW /
#     B01-B02BEHW202 (NVEL BLMVOL). ---
const _PN_VOLEQ_BY_FOREST = Dict{Int,NTuple{39,String}}(
    609 => ("616BEHW011", "616BEHW015", "616BEHW017", "616BEHW019", "616BEHW020", "F03FW2W263", "616BEHW022", "616BEHW042", "616BEHW081", "616BEHW093", "616BEHW108", "616BEHW116", "616BEHW117", "616BEHW119", "616BEHW122", "F03FW2W202", "616BEHW211", "616BEHW242", "F00FW2W263", "616BEHW264", "616BEHW312", "616BEHW351", "616BEHW352", "616BEHW375", "616BEHW431", "616BEHW746", "616BEHW747", "616BEHW815", "616BEHW064", "616BEHW072", "616BEHW101", "616BEHW103", "616BEHW231", "616BEHW492", "616BEHW500", "616BEHW768", "616BEHW920", "616BEHW000", "616BEHW999"),
    612 => ("616BEHW011", "616BEHW015", "616BEHW017", "616BEHW019", "616BEHW020", "616BEHW098", "616BEHW022", "616BEHW042", "616BEHW081", "616BEHW093", "616BEHW108", "616BEHW116", "616BEHW117", "616BEHW119", "616BEHW122", "F00FW2W202", "616BEHW211", "616BEHW242", "F03FW2W263", "616BEHW264", "616BEHW312", "NVBM240351", "616BEHW352", "616BEHW375", "616BEHW431", "616BEHW746", "616BEHW747", "616BEHW815", "616BEHW064", "616BEHW072", "616BEHW101", "616BEHW103", "616BEHW231", "616BEHW492", "616BEHW500", "616BEHW768", "616BEHW920", "616BEHW000", "616BEHW999"),
    800 => ("616BEHW011", "616BEHW015", "616BEHW017", "616BEHW019", "616BEHW020", "F03FW2W263", "616BEHW022", "616BEHW042", "616BEHW081", "616BEHW093", "616BEHW108", "616BEHW116", "616BEHW117", "616BEHW119", "616BEHW122", "F03FW2W202", "616BEHW211", "616BEHW242", "F00FW2W263", "616BEHW264", "616BEHW312", "616BEHW351", "616BEHW352", "616BEHW375", "616BEHW431", "616BEHW746", "616BEHW747", "616BEHW815", "616BEHW064", "616BEHW072", "616BEHW101", "616BEHW103", "616BEHW231", "616BEHW492", "616BEHW500", "616BEHW768", "616BEHW920", "616BEHW000", "616BEHW999"),
    708 => ("B00BEHW011", "B00BEHW015", "B00BEHW017", "B00BEHW015", "B00BEHW021", "B00BEHW098", "B00BEHW022", "B00BEHW042", "B00BEHW081", "B00BEHW093", "B00BEHW108", "B00BEHW116", "B00BEHW117", "B00BEHW119", "B00BEHW122", "B01BEHW202", "B00BEHW211", "B00BEHW242", "B00BEHW263", "B00BEHW260", "B00BEHW312", "B00BEHW351", "B00BEHW361", "B00BEHW999", "B00BEHW431", "B00BEHW999", "B00BEHW747", "B00BEHW800", "B00BEHW242", "B00BEHW073", "B00BEHW119", "B00BEHW108", "B00BEHW231", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999"),
    709 => ("B00BEHW011", "B00BEHW015", "B00BEHW017", "B00BEHW015", "B00BEHW021", "B00BEHW098", "B00BEHW022", "B00BEHW042", "B00BEHW081", "B00BEHW093", "B00BEHW108", "B00BEHW116", "B00BEHW117", "B00BEHW119", "B00BEHW122", "B01BEHW202", "B00BEHW211", "B00BEHW242", "B00BEHW263", "B00BEHW260", "B00BEHW312", "B00BEHW351", "B00BEHW361", "B00BEHW999", "B00BEHW431", "B00BEHW999", "B00BEHW747", "B00BEHW800", "B00BEHW242", "B00BEHW073", "B00BEHW119", "B00BEHW108", "B00BEHW231", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999"),
    712 => ("B00BEHW011", "B00BEHW015", "B00BEHW017", "B00BEHW015", "B00BEHW021", "B00BEHW098", "B00BEHW022", "B00BEHW042", "B00BEHW081", "B00BEHW093", "B00BEHW108", "B00BEHW116", "B00BEHW117", "B00BEHW119", "B00BEHW122", "B02BEHW202", "B00BEHW211", "B00BEHW242", "B00BEHW263", "B00BEHW260", "B00BEHW312", "B00BEHW351", "B00BEHW361", "B00BEHW999", "B00BEHW431", "B00BEHW999", "B00BEHW747", "B00BEHW800", "B00BEHW242", "B00BEHW073", "B00BEHW119", "B00BEHW108", "B00BEHW231", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999", "B00BEHW999"),
)

"VOLEQ for PN species `sp` on forest `kodfor` (post-FORKOD JFOR code); unknown codes use the grinit 612 row."
_pn_voleq(kodfor::Int, sp::Int)::String =
    (1 <= sp <= 39) ? get(_PN_VOLEQ_BY_FOREST, kodfor, _PN_VOLEQ_BY_FOREST[612])[sp] : "           "

# PN Behre per-tree volume — reuse the BM R6 machinery + PN form class.
function pn_behre_vol(sp::Int, ifor::Int, d::Float32, h::Float32, bark::Float32; topd::Float32 = 4.5f0)
    fclass = pn_formcl(sp, ifor, d)
    dbtbh = d * (1f0 - bark); dbhib = d - dbtbh
    vol2 = 0f0; vol4 = 0f0
    v1 = if h <= 17.3f0
        0.00272708f0 * dbhib * dbhib * h
    else
        v = bm_r6vol3(d, dbtbh, fclass, h, 1)
        mtopp = topd * bark
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

function compute_volumes_pn!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    t = s.trees; veq = s.species.vol_eq; sd = s.coef.species; c = s.control
    ifor = Int(s.plot.forest_idx)                    # forkod JFOR index (1..6): form class + BLM merch specs
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; h = t.height[i]; sp = Int(t.species[i])
        if d < 1f0 || sp < 1 || sp > 39
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0
            t.saw_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0; continue
        end
        eq = veq[sp]; se = strip(eq); mdl = length(se) >= 7 ? se[4:6] : "   "
        # vols.f:132,150-151: BARK=BRATIO(ISPC,DBH_start,H) before `D=D+DG(I)/BARK` ⇒ projected cycles use the stashed
        # start-of-cycle bark (t.vol_bark) for the merch tops / DBTBH / CFTOPK; grown-DBH bark at cycle 0 / dead records.
        bark = (i <= t.n && t.vol_bark[i] > 0f0) ? t.vol_bark[i] : wc_bratio(sd, sp, d)
        # pn/sitset.f:188-209 merch specs (init_merch_standards!): BLM IFOR 4-6 TOPD=BFTOPD=5, DBHMIN=7 for all;
        # otherwise TOPD=BFTOPD=4.5, DBHMIN=BFMIND=7 (LP sp11 = 6).
        dbhmin = c.sp_dbh_min[sp]; bfmind = c.sp_bf_dbhmin[sp]
        topd = c.sp_top_diam[sp]; stmp = c.sp_stump_ht[sp]
        bftopd = c.sp_bf_topd[sp]; bfstmp = c.sp_bf_stump[sp]
        # pn/vols.f:145-146 (= wc/vols.f): a top-killed tree (H>=4.5, ITRUNC>0) is volumed at its NORMAL height …
        tkill = h >= 4.5f0 && t.trunc[i] > 0
        hv = tkill ? Float32(t.norm_ht[i]) / 100f0 : h
        local tcf::Float32, mcf::Float32, bf::Float32
        if mdl == "FW2" && (se[1] == 'F' || se[1] == 'f')
            tcf, mcf, bf = wc_fw2_westside_vol(eq, d, hv, bark; topd = topd, bftopd = bftopd, stump = stmp)
        elseif mdl == "FW2"
            v = cr_fw2_vol(eq, d, hv; bark = bark, topd = topd, bftopd = bftopd, stump = stmp, iregn = 6, board_cor = 'N', merch_opt = 23)
            tcf = max(v[1], 0f0); mcf = max(v[4] + v[7], 0f0); bf = max(v[2], 0f0)
        elseif se[1] == 'B' || se[1] == 'b'           # BLM forests (R7_EQN) → NVEL BLMVOL
            tcf, mcf, bf = _blm_natcrs(eq, pn_formcl(sp, ifor, d), d, hv, bark, topd, bftopd, bfmind)
        elseif startswith(se, "NVB")                  # RA on 612: NVBM240351 (NSVB) — the CR NVB kernel, region 6
            # fvsvol.f:85-88: a LIVE top-killed tree passes BRKHT=ITRNC/100; NSVB trims VOL(1) and HT1PRD itself
            # (pn/vols.f:191 skips CFTOPK for 'NVB'; BFTOPK at vols.f:390 still runs, below)
            brk = (tkill && i <= t.n) ? Float32(t.trunc[i]) / 100f0 : 0f0
            v = cr_nvb_vol(eq, d, hv; bark = bark, topd = topd, stump = stmp, bftopd = bftopd, iregn = 6, brkht = brk)
            tcf = max(v[1], 0f0); mcf = max(v[4], 0f0); bf = max(v[2], 0f0)
        else                                          # 616BEHW
            tcf, mcf, bf = pn_behre_vol(sp, ifor, d, hv, bark; topd = topd)
        end
        tcf = max(tcf, 0f0)
        mcf = d >= dbhmin ? max(mcf, 0f0) : 0f0      # fvsvol.f:513 MCF only for D>=DBHMIN
        bf  = d >= bfmind ? max(bf, 0f0) : 0f0       # fvsvol.f:517 BBFV=0 for D<BFMIND
        # … then trimmed to the standing broken stem by CFTOPK/BFTOPK (pn/vols.f:191-193, 390-391).
        if tkill && tcf > 0f0
            vmax = tcf        # NVB: BFMAX = the board call's TVOL(1) = the already-truncated Vtotib·Rrem (fvsvol.f:495)
            startswith(se, "NVB") ||
                ((tcf, mcf) = cr_cftopk(tcf, mcf, d, hv, vmax, bark, Int(t.trunc[i]), stmp, topd))
            bf = cr_bftopk(bf, d, hv, vmax, bark, Int(t.trunc[i]), bfstmp, bftopd)
        end
        t.cuft_vol[i] = tcf
        t.merch_cuft_vol[i] = mcf
        t.saw_cuft_vol[i] = 0f0
        t.bdft_vol[i] = bf
    end
    return s
end
