# =============================================================================
# volume.jl (olympic) — OP BLM Behre-taper cubic + board volume. Chunk 2.
#
# OP's VOLEQDEF equations are BLM Behre (`B00/B01 BEHW`+FIA — dumped from the live FVSop_clean
# "NATIONAL VOLUME ESTIMATOR LIBRARY EQUATION NUMBERS" table, forest 708); `volinit.f:367` routes
# `VOLEQ(1:1)=='B'` to BLMVOL. The NVEL BLM kernel (blmvol.f + blmtap.f) is BYTE-IDENTICAL across the
# buildDirs; it is the one engine implementation `blm_vol` (src/engine/blm_vol.jl). Only the OP data differ:
#   • OP_VOLEQ    — op/voleqdef.f 39-species BLM Behre equations (DF = B01BEHW202).
#   • op_formcl   — op/formcl.f per-species BLM708 (Salem) form class (forest_idx 4; no DBH-class dep).
#   • op_bratio   — op/bratio.f bark (organon_nwo.jl), NOT oc_bratio.
#   • MTOPP = TOPD·BARK with TOPD=5.0 (op/sitset.f IFOR 4,5,6 BLM CASE), vs OC's 4.5; DBHMIN=BFMIND=7.
#
# The reference-stand S248112 species (WF/ES/LP/SP/PP/DF) are all Behre; DF routes through the
# B01 (profile 1) branch, the rest through B00 (profile 10). VALIDATED BIT-EXACT per-tree (27/27
# records, TCF/MCF/BF to the F10.3 print) vs the live oracle: a WRITE was added to a temp copy of
# fvsvol.f (buildDir kept PRISTINE — marker 0), FVSop_dbg relinked, cyc0 D/HT/TCF/MCF/BF captured.
# The aggregate cyc0 .sum is TCuFt/MCuFt/BdFt = 1472/972/5003 (reached once CRATET/ORGANON height-
# dubbing supplies the two dubbed heights 66.31 DF / 62.39 LP — chunk 3). See test_op_site_and_volume.jl.
# =============================================================================

# op/voleqdef.f VOLEQDEF (VAR='OP', forest 708) per FVS species 1..39 — the authoritative table dumped
# from the live FVSop_clean NVEL equation-number listing.
const OP_VOLEQ = String[
    "B00BEHW011","B00BEHW015","B00BEHW017","B00BEHW015","B00BEHW021",   # SF WF GF AF RF
    "B00BEHW098","B00BEHW022","B00BEHW042","B00BEHW081","B00BEHW093",   # SS NF YC IC ES
    "B00BEHW108","B00BEHW116","B00BEHW117","B00BEHW119","B00BEHW122",   # LP JP SP WP PP
    "B01BEHW202","B00BEHW211","B00BEHW242","B00BEHW263","B00BEHW260",   # DF RW RC WH MH
    "B00BEHW312","B00BEHW351","B00BEHW361","B00BEHW631","B00BEHW431",   # BM RA MA TO GC
    "B00BEHW999","B00BEHW747","B00BEHW800","B00BEHW242","B00BEHW073",   # AS CW WO WJ LL
    "B00BEHW119","B00BEHW108","B00BEHW231","B00BEHW999","B00BEHW999",   # WB KP PY DG HT
    "B00BEHW999","B00BEHW999","B00BEHW999","B00BEHW999"]                # CH WI __ OT

# op/formcl.f BLM708 (Salem) per-species Girard form class (IFOR=4; single value, no DBH class).
const OP_FORMCL_BLM708 = Int[
    84,86,84,82,75, 80,84,73,73,77, 68,75,75,76,82, 80,75,76,88,72,
    84,88,70,70,75, 75,74,70,60,75, 82,82,60,70,70, 75,75,74,74]

"op/formcl.f FORMCL (IFOR=4 BLM Salem): per-species BLM708 form class."
@inline op_formcl(sp::Int) = OP_FORMCL_BLM708[sp]

"blmvol.f BLMVOL total-cubic driver for OP (mirror of oc_tree_cuft, MTOPP=TOPD·BARK with TOPD=5.0)."
function op_tree_cuft(sp::Int, dbh::Float32, ht::Float32; topd::Float32 = 5.0f0, topbark::Float32 = -1f0)
    dbh <= 0f0 && return 0f0
    # MTOPP = TOPD·BARK with BARK = BRATIO at the START-of-cycle DBH (vols.f:150-151, `BARK=BRATIO(D)`
    # BEFORE `D=D+DG/BARK`). Pass the stashed `vol_bark` via `topbark`; fall back to the grown-DBH bark
    # (cyc0 / no stash, where pre-growth==current). The 0.001 pre-vs-grown bark difference flips a Scribner
    # log class on broken-top trees (op2c tree 19: last log 16ft/dib4 vs the oracle's 14ft/dib5).
    mtopp = topd * (topbark > 0f0 ? topbark : op_bratio(sp, dbh))
    return blm_vol(OP_VOLEQ[sp], mtopp, ht, dbh, op_formcl(sp))[1]
end

"blmvol.f BLM merch-cubic VOL(4) + Scribner VOL(2) for OP (mirror of oc_tree_mvol, TOPD=5.0)."
function op_tree_mvol(sp::Int, dbh::Float32, ht::Float32; topd::Float32 = 5.0f0, topbark::Float32 = -1f0,
                      ht_out = nothing)
    ht_out === nothing || (ht_out[] = 0f0)
    dbh <= 0f0 && return (0f0, 0f0)
    mtopp = topd * (topbark > 0f0 ? topbark : op_bratio(sp, dbh))   # BRATIO(D_start) top bark (vols.f:150) — see op_tree_cuft
    _, v2, v4 = blm_vol(OP_VOLEQ[sp], mtopp, ht, dbh, op_formcl(sp); bfpflg = true, ht_out = ht_out)
    return (v4, v2)
end

"""
    compute_volumes_op!(s)

OP BLM volume (Behre taper). Fills `cuft_vol` (total cubic, `op_tree_cuft`), `merch_cuft_vol`
(VOL(4), gated D≥DBHMIN) and `bdft_vol` (VOL(2) Scribner, gated D≥BFMIND) for every record — the
exact structure of `compute_volumes_oc!` with the OP data + TOPD=5.0 (op/sitset.f BLM CASE).
"""
function compute_volumes_op!(s::StandState)
    s.control.merch_init || init_merch_standards!(s)
    c = s.control
    t = s.trees
    topd = c.sp_top_diam[1] > 0f0 ? c.sp_top_diam[1] : 5.0f0
    merch = (stmp = c.sp_stump_ht, topd = c.sp_top_diam, scfstmp = c.sp_scf_stump,
             scftop = c.sp_scf_topd, bftopd = c.sp_bf_topd, bfstmp = c.sp_bf_stump)
    # vols.f zeroes HT2TD for every record; NATCRS (fvsvol.f:337-339 / 484-487) stores BLMVOL's HT1PRD.
    fill!(t.merch_top_cf, 0f0); fill!(t.merch_top_bf, 0f0)
    ht1 = Ref(0f0)
    @inbounds for i in 1:(t.n + t.ndead)
        sp = Int(t.species[i])
        if 1 <= sp <= 39
            d = t.dbh[i]; h = t.height[i]
            tkill = h >= 4.5f0 && t.trunc[i] > 0
            htap = tkill ? Float32(t.norm_ht[i]) * 0.01f0 : h
            # Merch-top bark = BRATIO(D_start) (vols.f:150) — the stashed `vol_bark` (0 at cyc0 ⇒ grown-DBH bark).
            topbark = t.vol_bark[i] > 0f0 ? t.vol_bark[i] : op_bratio(sp, d)
            tcf = op_tree_cuft(sp, d, htap; topd = topd, topbark = topbark)
            v4, v2 = op_tree_mvol(sp, d, htap; topd = topd, topbark = topbark, ht_out = ht1)
            d >= c.sp_dbh_min[sp]   && (t.merch_top_cf[i] = ht1[])
            d >= c.sp_bf_dbhmin[sp] && (t.merch_top_bf[i] = ht1[])
            mcf = d >= c.sp_dbh_min[sp]   ? v4 : 0f0
            bf  = d >= c.sp_bf_dbhmin[sp] ? v2 : 0f0
            if tkill && tcf > 0f0
                bark = t.vol_bark[i] > 0f0 ? t.vol_bark[i] : op_bratio(sp, d)
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
