# =============================================================================
# organon_volume.jl — OC (Oregon Coast) BLM cubic volume (chunk C10a).
#
# The NVEL BLM kernel (blmvol.f + blmtap.f) lives in src/engine/blm_vol.jl (`blm_vol`, one implementation
# shared by OC/OP/NC/WC/PN); VOLEQDEF gives OC BLM Behre equations (`B00/B01 BEHW`+FIA) and volinit.f routes
# `VOLEQ(1:1)=='B'` to BLMVOL. This file holds only the OC data: OC_VOLEQ and oc_formcl (oc/formcl.f).
#
# MEASURED vs FVSoc_clean TREELIST (ocmin cyc0): per-tree total cuft (e.g. DF D12.7/H67 → 30.1,
# GF D6.2/H38 → 6.7, LP D11.5/H73 → 16.9). Board-foot (SCRIB) is C10b. See docs/OC_VARIANT_PORT_AUDIT.md.
# =============================================================================

# oracle VOLEQDEF table (VAR='OC', REGN 7, FORST 11) — per FVS species index 1..50.
const OC_VOLEQ = String[
    "B00BEHW081","B00BEHW081","B00BEHW242","B00BEHW017","B00BEHW021","B00BEHW021","B01BEHW202",
    "B00BEHW263","B00BEHW260","B00BEHW119","B00BEHW108","B00BEHW108","B00BEHW108","B00BEHW108",
    "B00BEHW116","B00BEHW117","B00BEHW119","B00BEHW122","B00BEHW108","B00BEHW108","B00BEHW242",
    "B00BEHW093","B00BEHW211","B00BEHW231","B00BEHW999","B00BEHW800","B00BEHW800","B00BEHW800",
    "B00BEHW800","B00BEHW800","B00BEHW800","B00BEHW800","B00BEHW800","B00BEHW312","B00BEHW800",
    "B00BEHW351","B00BEHW361","B00BEHW431","B00BEHW999","B00BEHW312","B00BEHW999","B00BEHW631",
    "B00BEHW800","B00BEHW999","B00BEHW747","B00BEHW999","B00BEHW231","B00BEHW631","B00BEHW999",
    "B00BEHW211"]

# oc/formcl.f BLM711 (Medford) per-species form class.
const OC_FORMCL_BLM711 = Float32[
    70.,78.,78.,91.,78.,78.,87.,91.,70.,73.,68.,68.,68.,68.,70.,84.,76.,85.,68.,68.,
    70.,74.,70.,76.,70.,80.,80.,80.,80.,80.,80.,80.,80.,84.,80.,88.,81.,83.,83.,84.,
    83.,84.,80.,72.,72.,75.,76.,84.,70.,75.]

"oc/formcl.f FORMCL (IFOR=9 Medford): per-species BLM711 form class, INT-truncated."
@inline oc_formcl(sp::Int) = trunc(Int, OC_FORMCL_BLM711[sp])

"""
    oc_tree_cuft(sp, dbh, ht) -> total cubic feet

blmvol.f BLMVOL total-cubic driver for OC. Small trees (total ht ≤ 17.8 or SMD_17 < MTOPP) use the
closed form `0.00272708·DBHIB²·TTH`; larger trees use `BLMTCUB` (Behre-taper Smalian). TTH = HT+1.5.
"""
function oc_tree_cuft(sp::Int, dbh::Float32, ht::Float32)
    dbh <= 0f0 && return 0f0
    # fvsvol MTOPP = TOPD(=4.5, sitset.f:244 IFOR 6-10)·BARK
    return blm_vol(OC_VOLEQ[sp], 4.5f0 * oc_bratio(sp, dbh), ht, dbh, oc_formcl(sp))[1]
end

"""
    oc_tree_mvol(sp, dbh, ht) -> (merch_cuft, bdft)

BLM merchantable cubic (VOL(4), Smalian on rounded log DIBs) and Scribner board-foot (VOL(2)) for OC.
Shares one 16-ft bucking (BLMMLEN→NUMLOG→SEGMNT→BLMGDIB) between the two products, since fvsvol's cubic
and board-foot VOLINIT calls use the same MTOPP=4.5·BARK, D17, STUMP=1 here. The fvsvol DBHMIN/BFMIND
gates are applied by the caller. Small trees / merch length < 8 ft ⇒ (0,0).
"""
function oc_tree_mvol(sp::Int, dbh::Float32, ht::Float32; ht_out = nothing)
    ht_out === nothing || (ht_out[] = 0f0)
    dbh <= 0f0 && return (0f0, 0f0)
    _, v2, v4 = blm_vol(OC_VOLEQ[sp], 4.5f0 * oc_bratio(sp, dbh), ht, dbh, oc_formcl(sp); bfpflg = true,
                        ht_out = ht_out)   # HT1PRD (blmvol.f:427), the same for the cubic and board calls
    return (v4, v2)
end

"""
    compute_volumes_oc!(s)

OC BLM volume (Behre taper). Fills `cuft_vol` (total cubic, `oc_tree_cuft`), `merch_cuft_vol`
(VOL(4), gated D≥DBHMIN) and `bdft_vol` (VOL(2) Scribner, gated D≥BFMIND) for every record.
The `.sum` aggregates `vol·tpa/GROSPC` in the shared reporter.
"""
function compute_volumes_oc!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    c = s.control
    t = s.trees
    # Broken-top merch standards for CFTOPK/BFTOPK: raw grinit TOPD/STMP (=4.5/1, NOT ×BARK).
    merch = (stmp = c.sp_stump_ht, topd = c.sp_top_diam, scfstmp = c.sp_scf_stump,
             scftop = c.sp_scf_topd, bftopd = c.sp_bf_topd, bfstmp = c.sp_bf_stump)
    # vols.f:104 zeroes HT2TD for every record; NATCRS (fvsvol.f:337-339 / 484-487) stores BLMVOL's HT1PRD.
    fill!(t.merch_top_cf, 0f0); fill!(t.merch_top_bf, 0f0)
    ht1 = Ref(0f0)
    @inbounds for i in 1:(t.n + t.ndead)
        sp = Int(t.species[i])
        if 1 <= sp <= 50
            d = t.dbh[i]; h = t.height[i]
            # Broken/dead-top trees (ITRUNC>0, H≥4.5): build the profile from the dubbed NORMAL
            # height, then truncate back to the break with CFTOPK/BFTOPK (vols.f:164-165, cftopk.f).
            tkill = h >= 4.5f0 && t.trunc[i] > 0
            htap = tkill ? Float32(t.norm_ht[i]) * 0.01f0 : h
            eq = r6_voleq_oc_op(s, sp)
            if eq === nothing
                # BLM forests (710/711/712 …): the NVEL BLMVOL Behre-taper path (VEQNNC B0xBEHW…).
                tcf = oc_tree_cuft(sp, d, htap)
                v4, v2 = oc_tree_mvol(sp, d, htap; ht_out = ht1)
                d >= c.sp_dbh_min[sp]   && (t.merch_top_cf[i] = ht1[])
                d >= c.sp_bf_dbhmin[sp] && (t.merch_top_bf[i] = ht1[])
            else
                # Region-6 national forests (610 Rogue River / 611 Siskiyou): organon/vols.f → NATCRS on the R6_EQN
                # VEQNNC — westside Flewelling FW2 (F06 DF/WH), INGY FW2 (I00 WF/PP) or Behre 616BEHW (ca/formcl.f form
                # class, which OC links) — the same NVEL kernels CA uses for these forests (MEASURED live FVSoc 610/611:
                # VEQNNC 616BEHW… / F06FW2W202 / I00FW2W073; the BLM-only port gave cyc-0 TCuFt 18506 vs live 12329).
                bark = (i <= t.n && t.vol_bark[i] > 0f0) ? t.vol_bark[i] : oc_bratio(sp, d)
                htb = zeros(Float32, 2)
                se = strip(eq); mdl = length(se) >= 7 ? se[4:6] : "   "
                if mdl == "FW2" && (se[1] == 'F' || se[1] == 'f')
                    tcf, v4, v2 = wc_fw2_westside_vol(eq, d, htap, bark; ht2td = htb)
                elseif mdl == "FW2"
                    v = cr_fw2_vol(eq, d, htap; bark = bark, topd = 4.5f0, bftopd = 4.5f0, stump = 1f0,
                                   iregn = 6, board_cor = 'N', merch_opt = 23, sf_hs = true, ht2td = htb)
                    tcf = max(v[1], 0f0); v4 = max(v[4] + v[7], 0f0); v2 = max(v[2], 0f0)
                else
                    ifor = Int(s.plot.forest_idx)
                    tcf, v4, v2 = ca_behre_vol(sp, ifor, d, htap, bark)
                    htb[1] = r6vol_ht1prd(d, ca_formcl(sp, ifor, d), 4.5f0 * bark, htap, d - d * (1f0 - bark)); htb[2] = htb[1]
                end
                tcf = max(tcf, 0f0); v4 = max(v4, 0f0); v2 = max(v2, 0f0)
                d >= c.sp_dbh_min[sp]   && (t.merch_top_cf[i] = htb[1])
                d >= c.sp_bf_dbhmin[sp] && (t.merch_top_bf[i] = htb[2])
            end
            mcf = d >= c.sp_dbh_min[sp]  ? v4 : 0f0
            bf  = d >= c.sp_bf_dbhmin[sp] ? v2 : 0f0
            if tkill && tcf > 0f0
                bark = t.vol_bark[i] > 0f0 ? t.vol_bark[i] : oc_bratio(sp, d)
                vmax = tcf
                tcf, mcf, _ = cftopk(merch, sp, d, htap, tcf, mcf, 0f0, vmax, bark, Int(t.trunc[i]))
                bf = bftopk(merch, sp, d, htap, bf, vmax, bark, Int(t.trunc[i]))
            end
            t.cuft_vol[i] = tcf
            t.merch_cuft_vol[i] = mcf
            t.bdft_vol[i] = bf
        else
            t.cuft_vol[i] = 0f0; t.merch_cuft_vol[i] = 0f0; t.bdft_vol[i] = 0f0
        end
        t.saw_cuft_vol[i] = 0f0
    end
    return s
end
