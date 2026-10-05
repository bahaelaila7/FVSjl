# =============================================================================
# fire/west_ffe.jl — the western (non-eastern) FFE volume basis for snags and the live-carbon stem.
#
# {v}/fmsvol.f (identical in every western variant, :146-151): FMSVOL/FMSVL2 call NATCRS on the tree's DBH and
# height with XHT=−1 for an intact stem ⇒ LTKIL=.FALSE. ⇒ no CFTOPK, and return VOL2HT = MAX(X,TCF)
# (X = 0.005454154·H) — or MCF when LMERCH (the merch-carbon pool, fmcrbout.f:127 LVWEST). A snag that has lost
# height (FMSNGHT) is volumed as its death-form tree (DBHS, HTDEAD) trimmed by CFTOPK at IHT=INT(HTIH·100).
# Only CS/LS/NE/SN take the MAX(X,MCF)/SCF branch.
#
# `ffe_west_nocut` gives each western variant's NATCRS (TCF, MCF) on an arbitrary (d, h) with no top-kill — the same
# kernels as its compute_volumes (FW2 / MATW r4vol / DVE / Behre / NVB), minus the CFTOPK broken-top trim — plus the
# FMSVOL bark ratio BRATIO(JS,D,H) and whether that equation family takes the CFTOPK trim (`trim`: fmsvol.f calls
# NATCRS, which returns CTKFLG=.TRUE. (fvsvol.f:535) for every family — the CI/TT/UT DVE woodland included — and
# FMSVOL has no vols.f:191 NVB exemption; EM included (em/fmsvol.f:139-140)).
# `nothing` for a variant not yet on this layer (its own snag/live paths stay in force).
# =============================================================================

"Western variants whose FFE snag bole / live-carbon stem use `ffe_west_nocut` (fmsvol.f non-eastern branch)."
_ffe_west_vol(v) = v isa InlandEmpire || v isa Kootenai || v isa CentralIdaho || v isa Teton || v isa Utah ||
                   v isa EasternMontana || v isa CentralRockies || v isa EastCascades ||
                   v isa WestCascades || v isa PacificNorthwest || v isa SoutheastAlaska ||
                   v isa SouthCentralOregon || v isa Klamath || v isa BlueMountains ||
                   v isa WestSierra || v isa CentralCalifornia || v isa BritishColumbia

"""
FMSVOL's tiny-tree cone floor X = 0.005454154·H (VOL2HT = MAX(X,TCF), fmsvol.f:146-151) — absent in BC: canada/fire/bc/
fmsvol.f:121-124 computes X but returns VOL2HT = VN (VM when LMERCH). `max(_ffe_xfloor(v, h), vol)` = the right form.
"""
@inline _ffe_xfloor(v, h::Float32)::Float32 = v isa BritishColumbia ? 0f0 : 0.005454154f0 * h

"""
    bc_fmsvol(s, sp, d, h, xht) -> (vn, vm)

canada/fire/bc/fmsvol.f FMSVOL/FMSVL2 (as relinked in the private oracle /workspace/.bcwork/dbfix — the CFTOPK call gets
its May-2024 SCF argument, a zero): XHT ≤ −1 (an intact stem) ⇒ LTKIL=.FALSE., XHT=H; IHT = INT(XHT·100). CFVOL (the
Kozak taper, canada/bc/cfvol.f) volumes the stem truncated at HTRUNC = IHT/100 (H when IHT=0); then, with CTKFLG=LTKIL
(fmsvol.f:105), a top-killed call (XHT > −1) also runs CFTOPK (ie/cftopk.f) on that ALREADY-truncated VN with VMAX=VN,
LCONE=.FALSE., BARK=BRATIO — a second, Behre-ratio trim of the total. STMP/TOPD = canada/bc/grinit.f 30 cm / 10 cm.
"""
function bc_fmsvol(s, sp::Int, d::Float32, h::Float32, xht::Float32)
    ltkil = xht > -1f0
    ltkil || (xht = h)
    iht = unsafe_trunc(Int32, xht * 100f0)
    vn, vm = bc_tree_vol(sp, d, h, Int(iht))
    if ltkil
        bark = bc_bratio(sp)
        vn, vm = cr_cftopk(vn, vm, d, h, vn, bark, Int(iht), 30f0 * BC_CMtoFT, 10f0 * BC_CMtoIN)
    end
    return (vn, vm)
end

"""
    ffe_west_nocut(s, sp, d, h) -> (tcf, mcf, bark, trim) | nothing

NATCRS total/merch cubic for species `sp` at DBH `d`, height `h`, with NO top-kill (fmsvol.f XHT=−1), for the
variants in `_ffe_west_vol`. MCF is gated at the variant's DBHMIN like its compute_volumes; `bark` is the variant
BRATIO that FMSVOL hands CFTOPK and `trim` whether that family takes the CFTOPK trim (`ffe_west_snag_vol_at`).
"""
function ffe_west_nocut(s::StandState, sp::Int, d::Float32, h::Float32)
    v = s.variant
    _ffe_west_vol(v) || return nothing
    # canada/fire/bc/fmsvol.f FMSVL2 → CFVOL has no D<1 gate (Kozak taper at any D>0). MEASURED FVSbc_instr Fir.20
    # simfire 2020 FMDOUT: a D=0.933in H=16.1ft PL VT 0.03817 live / 0 jl ⇒ Aboveground_Total_Live 17.2257 / 17.2252.
    if v isa BritishColumbia && d > 0f0 && h > 0f0 && 1 <= sp <= length(s.species.vol_eq)
        vn, vm = bc_fmsvol(s, sp, d, h, -1f0)
        return (vn, vm, bc_bratio(sp), true)
    end
    (d < 1f0 || h <= 0f0 || sp < 1 || sp > length(s.species.vol_eq)) && return (0f0, 0f0, 1f0, false)
    eq = s.species.vol_eq[sp]; se = strip(eq); mdl = length(se) >= 7 ? se[4:6] : "   "
    if v isa BritishColumbia                                   # bc: FMSVL2(XHT=−1) — CFVOL at INT(H·100)/100, no CFTOPK
        vn, vm = bc_fmsvol(s, sp, d, h, -1f0)
        return (vn, vm, bc_bratio(sp), true)
    elseif v isa BlueMountains                                 # bm: NVEL FW2 / 616BEHW R6VOL, NATCRS CTKFLG=T (fvsvol.f)
        tcf, mcf = bm_nocut_cuft(s, sp, d, h)
        return (tcf, mcf, bm_bratio(s.coef.species, sp, d), true)
    elseif v isa SoutheastAlaska                               # ak: NVEL F32/DVE/CUR/DEM (ak_tree_vol), NATCRS CTKFLG=T
        tcf, mcf, _ = ak_tree_vol(s, sp, d, h)
        return (max(tcf, 0f0), max(mcf, 0f0), ak_bratio(sp, d), true)
    elseif v isa EasternMontana
        # em/fmsvol.f:139-140 `IF(CTKFLG .AND. LTKIL) CALL CFTOPK` — NATCRS returns CTKFLG=.TRUE. for the DVE woodland
        # equations too (fvsvol.f:531), so every EM family trims (MEASURED FVSem_g16 684750664126144 inventory AS snag
        # D6.1 HTDEAD 35 XHT 26: VMAX 2.4 → VOL2HT 2.3469481).
        tcf, mcf = em_nocut_cuft(s, sp, d, h)
        return (tcf, mcf, em_bratio(sp, d), true)
    elseif v isa CentralRockies                                # cr: NVEL DVE / NVB / FW2 (compute_volumes_cr!)
        tcf, mcf = cr_nocut_cuft(s, sp, d, h)
        return (tcf, mcf, cr_bratio(s.coef.species, sp, d, Int(s.plot.model_type)), true)
    elseif v isa Kootenai                                      # kt: FW2 only (compute_volumes_kt!)
        bark = KT_BKRAT[sp]                                    # kt/bratio.f BRATIO = BKRAT(IS)
        startswith(eq, "I") || return (0f0, 0f0, bark, false)
        w = cr_fw2_vol(eq, d, h; bark = bark, topd = 4.5f0, bftopd = 4.5f0, stump = 1f0, iregn = 1)
        return (max(w[1], 0f0), d >= (sp == 7 ? 6f0 : 7f0) ? max(w[4] + w[7], 0f0) : 0f0, bark, true)
    elseif v isa EastCascades                                  # ec: compute_volumes_ec! (FW2 westside / FW2 / 616BEHW)
        # ec/fmsvol.f = the shared western FMSVOL: NATCRS at (D,H), BRATIO(JS,D,H), no top-kill; CTKFLG=.TRUE. for
        # every family (fvsvol.f), so a height-lost snag always takes the CFTOPK trim.
        bark = wc_bratio(s.coef.species, sp, d); dbhmin = sp == 7 ? 6f0 : 7f0
        if mdl == "FW2" && (se[1] == 'F' || se[1] == 'f')
            tcf, mcf, _ = wc_fw2_westside_vol(eq, d, h, bark)
        elseif mdl == "FW2"
            w = cr_fw2_vol(eq, d, h; bark = bark, topd = 4.5f0, bftopd = 4.5f0, stump = 1f0,
                           iregn = 6, board_cor = 'N', merch_opt = 23)
            tcf = w[1]; mcf = w[4] + w[7]
        else
            tcf, mcf, _ = ec_behre_vol(sp, Int(s.plot.forest_idx), d, h, bark)
        end
        return (max(tcf, 0f0), d >= dbhmin ? max(mcf, 0f0) : 0f0, bark, true)
    elseif v isa WestCascades || v isa PacificNorthwest        # wc/pn: compute_volumes_{wc,pn}! kernels at (D,H)
        c = s.control; bark = wc_bratio(s.coef.species, sp, d); ifor = Int(s.plot.forest_idx)
        dbhmin = c.sp_dbh_min[sp]; bfmind = c.sp_bf_dbhmin[sp]
        topd = c.sp_top_diam[sp]; stmp = c.sp_stump_ht[sp]; bftopd = c.sp_bf_topd[sp]
        if mdl == "FW2" && (se[1] == 'F' || se[1] == 'f')
            tcf, mcf, _ = wc_fw2_westside_vol(eq, d, h, bark; topd = topd, bftopd = bftopd, stump = stmp)
        elseif mdl == "FW2"
            w = cr_fw2_vol(eq, d, h; bark = bark, topd = topd, bftopd = bftopd, stump = stmp, iregn = 6,
                           board_cor = 'N', merch_opt = 23)
            tcf = w[1]; mcf = w[4] + w[7]
        elseif se[1] == 'B' || se[1] == 'b'                    # BLM forests → NVEL BLMVOL
            fc = v isa WestCascades ? wc_formcl(sp, ifor, d) : pn_formcl(sp, ifor, d)
            tcf, mcf, _ = _blm_natcrs(eq, fc, d, h, bark, topd, bftopd, bfmind)
        elseif v isa WestCascades && mdl == "CUR"              # 603 red alder PROFILE + R10TAP
            tcf, mcf, _ = wc_cur_vol(d, h, bark; topd = topd, bftopd = bftopd, stump = stmp, bfmind = bfmind)
        elseif v isa PacificNorthwest && startswith(se, "NVB") # 612 red alder NSVB (intact stem: BRKHT 0)
            w = cr_nvb_vol(eq, d, h; bark = bark, topd = topd, stump = stmp, bftopd = bftopd, iregn = 6, brkht = 0f0)
            tcf = w[1]; mcf = w[4]
        else                                                   # 616BEHW
            tcf, mcf, _ = v isa WestCascades ? wc_behre_vol(sp, ifor, d, h, bark; topd = topd) :
                                               pn_behre_vol(sp, ifor, d, h, bark; topd = topd)
        end
        return (max(tcf, 0f0), d >= dbhmin ? max(mcf, 0f0) : 0f0, bark, true)
    elseif v isa SouthCentralOregon                            # so: compute_volumes_so! kernels at (D,H)
        # so/fmsvol.f, fmsnag.f, fmcwd.f, fmcbio.f are byte-identical to wc/ec's (the shared western FMSVOL layer); SO's
        # volume equations live in the driver's per-forest tables (species.vol_eq stays blank), so jl's generic snag
        # path took the eastern R8-Clark branch ⇒ 0 bole volume (MEASURED FVSso_g16 374286168489998 2015 inventory
        # Standing_Dead 2.91 live / 0.05 jl) and the live-carbon merch stem came from Jenkins (Aboveground_Merch_Live
        # 34.43 / 35.07).
        ifor = Int(s.plot.forest_idx); bark = so_bratio(s.coef.species, sp, d)
        topd = (ifor <= 3 || ifor == 10) ? 4.5f0 : 6.0f0                        # so/sitset.f:167
        if !(ifor <= 3 || ifor == 10)                                          # Region 5 (R5_EQN)
            tcf, mcf, _ = nvel_r5_vol(SO_R5_VOLEQ[sp], d, h, bark, topd, topd)
        else
            eqs = (ifor == 2 || ifor == 3) ? SO_VOLEQ_FR[sp] : SO_VOLEQ[sp]
            if eqs[4:6] == "FW2"
                w = cr_fw2_vol(eqs, d, h; bark = bark, topd = topd, bftopd = topd, stump = 1f0,
                               iregn = 6, board_cor = 'N', merch_opt = 23, sf_hs = true)
                tcf = w[1]; mcf = w[4] + w[7]
            else
                tcf, mcf, _ = so_behre_vol(sp, ifor, d, h, bark)
            end
        end
        return (max(tcf, 0f0), d >= 9f0 ? max(mcf, 0f0) : 0f0, bark, true)            # so/grinit.f DBHMIN 9
    elseif v isa Klamath                                       # nc: compute_volumes_nc! kernels at (D,H)
        # nc/fmsvol.f & co. = the shared western layer too; NC keeps its equations in NC_VOL_EQ / NC_R6_VOL_EQ /
        # NC_R7_VOL_EQ (R5 WO2W/DVE, R6 Siskiyou, R7 BLM), not in species.vol_eq.
        c = s.control; sd = s.coef.species; ifor = Int(s.plot.forest_idx)
        bark = nc_bratio(sd[:bark1][sp], sd[:bark2][sp], Int(sd[:bark_imap][sp]), d)
        if ifor == 5 || ifor == 7
            tcf, mcf, _ = nc_blmvol_vol(sp, d, h, bark; fclass = ifor == 7 ? NC_BLM712_FC[sp] : 80f0, ifor = ifor)
        elseif ifor == 4
            eq6 = NC_R6_VOL_EQ[sp]; m6 = eq6[4:6]
            if m6 == "FW2" && (eq6[1] == 'F' || eq6[1] == 'f')
                tcf, mcf, _ = wc_fw2_westside_vol(eq6, d, h, bark; topd = c.sp_top_diam[sp], bftopd = c.sp_bf_topd[sp])
            elseif m6 == "FW2"
                w = cr_fw2_vol(eq6, d, h; bark = bark, topd = c.sp_top_diam[sp], bftopd = c.sp_bf_topd[sp], stump = 1f0,
                               iregn = 6, board_cor = 'N', merch_opt = 23, sf_hs = true)
                tcf = w[1]; mcf = w[4] + w[7]
            else
                tcf, mcf, _ = nc_behre_vol(sp, d, h, bark; topd = c.sp_top_diam[sp])
            end
        else
            eqn = NC_VOL_EQ[sp]; topib = c.sp_top_diam[sp] * bark; bftib = c.sp_bf_topd[sp] * bark
            tcf, mcf, _ = eqn[4:6] == "WO2" ? nc_wo2w_vol(eqn, d, h; mtopp = topib, bftop = bftib) :
                                               nc_r5harv_vol(eqn, d, h, topib)
        end
        return (max(tcf, 0f0), d >= 9f0 ? max(mcf, 0f0) : 0f0, bark, true)            # nc/grinit.f DBHMIN 9
    elseif v isa WestSierra                                    # ws: compute_volumes_ws! kernels at (D,H)
        # ws/fmsvol.f (and fmsnag/fmcwd/fmcbio) are byte-identical to so's — the shared western layer — but WS keeps
        # its equations in WS_VOL_EQ (species.vol_eq blank), so jl's generic snag path took the eastern R8-Clark branch
        # ⇒ 0 bole (MEASURED FVSws_g16 7689156010901 2006 inventory Standing_Dead 10.61 live / 0.03 jl).
        bark = ws_bratio(s.coef.species, sp, d); weq = WS_VOL_EQ[sp]
        if weq[4:6] == "WO2"
            s5 = _nc_r5tap_sp(weq[8:10])
            (s5 == 0 || h < 5f0) && return (0f0, 0f0, bark, true)
            dibat = ht -> nc_r5tap_dib(s5, d, h, Float32(ht))
            tcf = _nint(_fw2_tcubic(dibat, h) * 10.0f0) * 1f-1
            mcf = d >= WS_VOL_DBHMIN ? nc_wo2w_merch(dibat, h; mtopp = WS_VOL_TOPD * bark, stump = 1f0, minlen = 2f0,
                                                      merchl = 8f0) : 0f0
        else                                                   # DVE California hardwood (r5harv.f)
            tcf, mcf, _ = nc_r5harv_vol(weq, d, h, WS_VOL_TOPD * bark)
        end
        return (max(tcf, 0f0), d >= WS_VOL_DBHMIN ? max(mcf, 0f0) : 0f0, bark, true)
    elseif v isa CentralCalifornia                             # ca: compute_volumes_ca! kernels at (D,H)
        # ca/fmsvol.f = the shared western FMSVOL; CA's equations come from ca_voleq (species.vol_eq blank) ⇒ the same
        # 0-bole gap as WS (MEASURED FVSca_g16 374401353489998 2015 Standing_Dead 0.77 live / 0 jl).
        sd = s.coef.species; ifor = Int(s.plot.forest_idx); c = s.control
        bark = wc_bratio(sd[:bark1][sp], sd[:bark2][sp], Int(sd[:bark_imap][sp]), d)
        ceq = ca_voleq(ifor, sp); cse = strip(ceq); cm = length(cse) >= 7 ? cse[4:6] : "   "
        if cse[1] == '5'
            tcf, mcf, _ = nvel_r5_vol(ceq, d, h, bark, c.sp_top_diam[sp], c.sp_bf_topd[sp])
        elseif cse[1] == 'B'
            tcf, mcf, _ = _blm_natcrs(ceq, ca_formcl(sp, ifor, d), d, h, bark, c.sp_top_diam[sp], c.sp_bf_topd[sp],
                                      c.sp_bf_dbhmin[sp])
        elseif cm == "FW2" && (cse[1] == 'F' || cse[1] == 'f')
            tcf, mcf, _ = wc_fw2_westside_vol(ceq, d, h, bark)
        elseif cm == "FW2"
            w = cr_fw2_vol(ceq, d, h; bark = bark, topd = 4.5f0, bftopd = 4.5f0, stump = 1f0,
                           iregn = 6, board_cor = 'N', merch_opt = 23, sf_hs = true)
            tcf = w[1]; mcf = w[4] + w[7]
        else
            tcf, mcf, _ = ca_behre_vol(sp, ifor, d, h, bark)
        end
        return (max(tcf, 0f0), d >= (sp == 11 ? 6f0 : 7f0) ? max(mcf, 0f0) : 0f0, bark, true)   # ca/grinit.f DBHMIN
    elseif v isa InlandEmpire                                  # ie: region-6 Behre / FW2 / region-1-2 DVE
        bark = ie_bratio(sp, d); dbhmin = sp == 7 ? 6f0 : 7f0
        if occursin("BEH", eq)
            tcf, mcf, _, _ = ie_behre_vol(sp, Int(s.plot.forest_idx), d, h, bark)
            return (max(tcf, 0f0), d >= dbhmin ? max(mcf, 0f0) : 0f0, bark, true)
        end
        # NATCRS gets the stand's NVEL region (fvsvol.f:90-96 IREGN=KODFOR/100) exactly as compute_volumes! — the Colville
        # (KODFOR 621, REGN 6) merch rules (MEASURED FVSie_g16 374547584489998 2015 Aboveground_Merch_Live 39.78778 live,
        # region-1 rules 39.71698).
        iregn = fvsvol_iregn(s)
        bcor, mopt = iregn == 6 ? ('N', 23) : ('Y', 22)
        w = startswith(eq, "I") ?
            cr_fw2_vol(eq, d, h; bark = bark, topd = 4.5f0, bftopd = 4.5f0, stump = 1f0, iregn = iregn,
                       board_cor = bcor, merch_opt = mopt, sf_hs = true, ht2td = zeros(Float32, 2)) :
            ie_dve_vol(eq, d, h, true)
        return (max(w[1], 0f0), d >= dbhmin ? max(w[4] + w[7], 0f0) : 0f0, bark, true)
    end
    # CI / TT / UT (compute_volumes_{ci,tt,ut}!): MATW r4vol, FW2 (region 4), DVE woodland
    dbhmin = sp == 7 ? 7f0 : 8f0
    bark = v isa CentralIdaho ? ci_bratio(s.coef.species, sp, d) :
           v isa Teton        ? tt_bratio(sp, d) :
                                bark_ratio(s.calib.bark_a, s.calib.bark_b, sp, d)
    if mdl == "MAT"
        mtopp = v isa Teton ? 6f0 : 6f0 * bark                # TT: MTOPS=TOPD (no bark), CI/UT: TOPD·BARK
        tcf, mcf = r4vol_volumes(eq, d, h, mtopp, 0f0)
        return (max(tcf, 0f0), d >= dbhmin ? max(mcf, 0f0) : 0f0, bark, true)
    elseif mdl == "FW2" && !(v isa Teton)
        w = cr_fw2_vol(eq, d, h; bark = bark, topd = 6f0, bftopd = 6f0, stump = 1f0, iregn = 4)
        return (max(w[1], 0f0), d >= dbhmin ? max(w[4] + w[7], 0f0) : 0f0, bark, true)
    end
    # DVE woodland: NATCRS ends CTKFLG = .TRUE. for every equation (fvsvol.f:509-531), so fmsvol.f:141-142's CFTOPK
    # trims a lost-height DVE snag as well (trim = true).
    vol1 = (v isa Utah && se[1] == '3') ? r3d2hv_vol1(eq, d, h) : r4d2h_vol1(eq, d, h)
    return (max(vol1, 0f0), d >= (v isa Teton ? 8f0 : dbhmin) ? max(vol1, 0f0) : 0f0, bark, true)
end

"FMSVOL snag bole for an intact stem: MAX(0.005454154·H, TCF) (fmsvol.f:150). `nothing` off the western layer."
function ffe_west_snag_bole(s::StandState, sp::Int, d::Float32, h::Float32)
    w = ffe_west_nocut(s, sp, d, h); w === nothing && return nothing
    h <= 0f0 && return 0f0
    return max(_ffe_xfloor(s.variant, h), w[1])
end

"""
    ffe_west_snag_vol_at(s, sp, d, htd, htcur) -> Float32 | nothing

FMSVOL(II, XHT=HTIH): the snag's death-form tree (DBHS, HTDEAD) through NATCRS, then CFTOPK at IHT=INT(XHT·100)
when the snag has lost height (the fat lower bole, not a short tree), VOL2HT = MAX(0.005454154·HTDEAD, TCF).
"""
function ffe_west_snag_vol_at(s::StandState, sp::Int, d::Float32, htd::Float32, htcur::Float32; always::Bool = false)
    _ffe_west_vol(s.variant) || return nothing
    htd <= 0f0 && return 0f0
    # BC: every FMSVOL call passes XHT = HTIH/HTIS (> −1) ⇒ the CFVOL-truncated + CFTOPK-trimmed VN, no cone floor.
    s.variant isa BritishColumbia && return (d < 1f0 ? 0f0 : bc_fmsvol(s, sp, d, htd, htcur)[1])
    w = ffe_west_nocut(s, sp, d, htd)
    tcf, mcf, bark, trim = w
    # `always`: fmsvol.f FMSVOL passes XHT=HTIH (> −1) ⇒ LTKIL ⇒ CFTOPK at INT(XHT·100)/100 on EVERY snag — even an
    # un-broken one (HTRUNC then truncates HTDEAD to 0.01 ft, and TCF·VOLTK/VOLT rounds).
    if (htcur < htd || always) && tcf > 0f0 && trim
        # CFTOPK reads the species' STMP/TOPD (WC/PN BLM forests TOPD 5); the other layer variants keep 1 / 4.5.
        stmp, topd = (s.variant isa WestCascades || s.variant isa PacificNorthwest) ?
                     (s.control.sp_stump_ht[sp], s.control.sp_top_diam[sp]) : (1f0, 4.5f0)
        tcf, _ = cr_cftopk(tcf, mcf, d, htd, tcf, bark, unsafe_trunc(Int, htcur * 100f0), stmp, topd)
    end
    return max(0.005454154f0 * htd, tcf)
end

# ---------------------------------------------------------------------------------------------------------------
# Snag height loss (FMSNGHT) parameters. {v}/fmvinit.f HTX(I,1:4) per species (hard first-50% / hard 50-95% / soft
# first-50% / soft 50-95%; TT sets 1:2 and copies to 3:4, fmvinit.f:411-412), extracted from the species SELECT CASE.
# HTR1 = 0.0228, HTR2 = 0.01 in all six; HTXSFT = 2 (IE/KT/CI/TT), 10 (UT, CR — cr/fmvinit.f:132). jl seeded HTX only
# for NE/LS/EM, so these variants' snags never lost height (live IE S248112 LP input snags shrink every year).
const _IE_FM_HTX = NTuple{4,Float32}[(0.9f0, 0.9f0, 0.9f0, 0.9f0), (0.9f0, 0.9f0, 0.9f0, 0.9f0), (0.9f0, 0.9f0, 0.9f0, 0.9f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0)]
const _KT_FM_HTX = NTuple{4,Float32}[(0.9f0, 0.9f0, 0.9f0, 0.9f0), (0.9f0, 0.9f0, 0.9f0, 0.9f0), (0.9f0, 0.9f0, 0.9f0, 0.9f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.1f0, 1.1f0, 1.1f0, 1.1f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0)]
const _CI_FM_HTX = NTuple{4,Float32}[(0.4f0, 0.4f0, 0.4f0, 0.4f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.5f0, 1.5f0, 1.5f0, 1.5f0), (0.9f0, 0.9f0, 0.9f0, 0.9f0), (0.3f0, 0.3f0, 0.3f0, 0.3f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (1.5f0, 1.5f0, 1.5f0, 1.5f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0)]
const _TT_FM_HTX = NTuple{4,Float32}[(0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (1.01f0, 4.61f0, 1.01f0, 4.61f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.001f0, 0.001f0, 0.001f0, 0.001f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0)]
const _UT_FM_HTX = NTuple{4,Float32}[(0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.9f0, 0.9f0, 0.9f0, 0.9f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (1.0f0, 1.0f0, 1.0f0, 1.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0)]

const _CR_FM_HTX = NTuple{4,Float32}[(1.494f0, 1.494f0, 1.494f0, 1.494f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (0.9f0, 0.9f0, 0.9f0, 0.9f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.9f0, 0.9f0, 0.9f0, 0.9f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (1.494f0, 1.494f0, 1.494f0, 1.494f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0978f0, 0.0978f0, 0.0978f0, 0.0978f0), (0.0462f0, 0.0462f0, 0.0462f0, 0.0462f0), (0.0f0, 0.0f0, 0.0f0, 0.0f0)]

# canada/fire/bc/fmvinit.f HTX(I,1:4): species 1-10 the KT values, 11-13/15 (hardwoods) 1.0, 14 (other conifer, as DF) 0.9.
const _BC_FM_HTX = NTuple{4,Float32}[ntuple(_ -> x, 4) for x in
    Float32[0.9, 0.9, 0.9, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.0, 1.0, 1.0, 1.0, 0.9, 1.0]]

"{v}/fmvinit.f HTX table for the western layer (EM keeps its own EM_FM_HTX seeding), else `nothing`."
_ffe_west_htx(v) = v isa SoutheastAlaska ? _AK_FM_HTX :
                   v isa InlandEmpire ? _IE_FM_HTX : v isa Kootenai ? _KT_FM_HTX : v isa CentralIdaho ? _CI_FM_HTX :
                   v isa BritishColumbia ? _BC_FM_HTX :
                   v isa Teton ? _TT_FM_HTX : v isa Utah ? _UT_FM_HTX : v isa CentralRockies ? _CR_FM_HTX : nothing
_snag_htr1(::InlandEmpire) = 0.0228f0
_snag_htr1(::Kootenai) = 0.0228f0
_snag_htr1(::BritishColumbia) = 0.0228f0   # canada/fire/bc/fmvinit.f:110 HTR1 (HTR2 0.01, HTXSFT 2.0)
_snag_htr1(::CentralIdaho) = 0.0228f0
_snag_htr1(::Teton) = 0.0228f0
_snag_htr1(::Utah) = 0.0228f0
_snag_htr1(::CentralRockies) = 0.0228f0
_snag_htr1(::SoutheastAlaska) = 0.02f0      # ak/fmvinit.f:140 HTR1
_snag_htr1(::Klamath) = 0.03406f0           # nc/fmvinit.f:122 HTR1 (50% height loss in 20 yr)
"HTR2 (after-50% snag height-loss rate, {v}/fmvinit.f): 0.01 everywhere but AK (ak/fmvinit.f:141 HTR2=0.02)."
_snag_htr2(v) = v isa SoutheastAlaska ? 0.02f0 : 0.01f0
"HTXSFT (soft-snag height-loss multiplier, {v}/fmvinit.f): UT/CR 10, the default 2 elsewhere."
_snag_htxsft(v) = (v isa Utah || v isa CentralRockies || v isa Klamath || v isa WestSierra ||
                   v isa CentralCalifornia || v isa OregonCoast) ? 10f0 :
                  (v isa SoutheastAlaska || r6_ffe_code(v) === :none) ? 2f0 : 1f0   # bm/ec/so/pn/op/wc HTXSFT 1.0; ak 2.0
# HTR1 of the R6 group (fmvinit.f): EC 0.0228; BM/SO/PN/OP/WC 0.03406.
_snag_htr1(::EastCascades) = 0.0228f0
_snag_htr1(::Union{BlueMountains,SouthCentralOregon,PacificNorthwest,Olympic,WestCascades}) = 0.03406f0

"""
Snag fall FMSFALL ({v}/fmsfall.f) for the western layer: the linear small-snag fall runs below 18" (IE/KT/CI/TT/UT;
EM 12"), with no redcedar exception, and some species fall slower once ≥ 18": BASE = MAX(0.01, BASE·0.32) for
CI ksp 2/8 (ci/fmsfall.f:39), TT 3/8, UT 3/5/8, CR (CASE DEFAULT) 3/17-19 (the shared tt/ut/cr fmsfall.f LDFSP).
Returns (linear_max, slow).
"""
function _ffe_west_fall(v, ksp::Integer)
    v isa EasternMontana && return (12f0, false)
    slow = v isa CentralIdaho ? (ksp == 2 || ksp == 8) :
           v isa Teton        ? (ksp == 3 || ksp == 8) :
           v isa Utah         ? (ksp == 3 || ksp == 5 || ksp == 8) :
           v isa CentralRockies ? (ksp == 3 || 17 <= ksp <= 19) : false
    return (18f0, slow)
end
