# =============================================================================
# volume.jl (olympic) — OP BLM Behre-taper cubic + board volume. Chunk 2.
#
# OP's VOLEQDEF equations are BLM Behre (`B00/B01 BEHW`+FIA — dumped from the live FVSop_clean
# "NATIONAL VOLUME ESTIMATOR LIBRARY EQUATION NUMBERS" table, forest 708); `volinit.f:367` routes
# `VOLEQ(1:1)=='B'` to BLMVOL. The NVEL BLM kernel (blmvol.f + blmtap.f) is BYTE-IDENTICAL between
# the OC and OP buildDirs, so the whole taper/bucking machinery is REUSED from OC's organon_volume.jl
# (`oc_blmtapeq`, `oc_blmtapeq_tapequ`, `oc_double_bark`, `oc_blmtap`, `oc_blmtcub`, `oc_blmmlen`,
# `oc_numlog`, `oc_segmnt!`, `oc_blmgdib!`, `oc_scrib`). Only the OP-specific data differ:
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
function op_tree_cuft(sp::Int, dbh::Float32, ht::Float32; topd::Float32 = 5.0f0)
    dbh <= 0f0 && return 0f0
    veq = OP_VOLEQ[sp]
    profile = oc_blmtapeq(veq)
    tapequ = oc_blmtapeq_tapequ(veq)
    dbhib = oc_double_bark(tapequ, dbh)
    dbhib <= 0.0001f0 && return 0f0
    mtopp = topd * op_bratio(sp, dbh)
    tth = ht + 1.5f0
    fclass = Float32(op_formcl(sp))
    if tth > 0f0
        tth <= 17.8f0 && return 0.00272708f0*(dbhib*dbhib)*tth
        smd_17 = trunc(sqrt(dbhib*dbhib - (dbhib*dbhib)*17.3f0/tth) + 0.5f0)
        smd_17 < mtopp && return 0.00272708f0*(dbhib*dbhib)*tth
    end
    d17 = round((dbh*fclass)/100.0f0, RoundNearestTiesAway)
    return oc_blmtcub(profile, dbh, tth, d17, 16.3f0)
end

"blmvol.f BLM merch-cubic VOL(4) + Scribner VOL(2) for OP (mirror of oc_tree_mvol, TOPD=5.0)."
function op_tree_mvol(sp::Int, dbh::Float32, ht::Float32; topd::Float32 = 5.0f0)
    dbh <= 0f0 && return (0f0, 0f0)
    veq = OP_VOLEQ[sp]
    profile = oc_blmtapeq(veq)
    tapequ = oc_blmtapeq_tapequ(veq)
    dbhib = oc_double_bark(tapequ, dbh)
    dbhib <= 0.0001f0 && return (0f0, 0f0)
    mtopp = topd * op_bratio(sp, dbh)
    tth = ht + 1.5f0
    if tth > 0f0
        tth <= 17.8f0 && return (0f0, 0f0)
        smd_17 = trunc(sqrt(dbhib*dbhib - (dbhib*dbhib)*17.3f0/tth) + 0.5f0)
        smd_17 < mtopp && return (0f0, 0f0)
    end
    fclass = Float32(op_formcl(sp))
    d17 = round((dbh*fclass)/100.0f0, RoundNearestTiesAway)
    stump = 1.0f0
    lmerch = oc_blmmlen(profile, tth, dbh, d17, stump, mtopp)
    lmerch < 8.0f0 && return (0f0, 0f0)
    numseg = oc_numlog(lmerch)
    loglen = zeros(Float32, 20)
    numseg = oc_segmnt!(loglen, lmerch, numseg)
    numseg == 0 && return (0f0, 0f0)
    logdia = zeros(Float32, 22)
    oc_blmgdib!(logdia, profile, mtopp, tth, dbh, dbhib, d17, stump, 0.3f0, numseg, loglen)
    v4 = 0f0
    dibl = round(dbhib, RoundNearestTiesAway)
    @inbounds for i in 1:numseg
        dibs = round(logdia[i+1], RoundNearestTiesAway)
        logv = 0.00272708f0*(dibl*dibl + dibs*dibs)*loglen[i]
        v4 += round(logv*10.0f0, RoundNearestTiesAway)/10.0f0
        dibl = dibs
    end
    v2 = 0f0
    @inbounds for i in 1:numseg
        dib = round(logdia[i+1], RoundNearestTiesAway)
        v2 += round(oc_scrib(dib, loglen[i]), RoundNearestTiesAway)
    end
    return (v4 < 0f0 ? 0f0 : v4, v2 < 0f0 ? 0f0 : v2)
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
    @inbounds for i in 1:(t.n + t.ndead)
        sp = Int(t.species[i])
        if 1 <= sp <= 39
            d = t.dbh[i]; h = t.height[i]
            tkill = h >= 4.5f0 && t.trunc[i] > 0
            htap = tkill ? Float32(t.norm_ht[i]) * 0.01f0 : h
            tcf = op_tree_cuft(sp, d, htap; topd = topd)
            v4, v2 = op_tree_mvol(sp, d, htap; topd = topd)
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
