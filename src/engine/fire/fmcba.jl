# =============================================================================
# fire/fmcba.jl — FFE per-cycle fuel & cover-type update (FFE chunk F3-state)
#
# Ported from: bin/FVSsn_buildDir/fmcba.f (FMCBA, the deterministic body).
#
# Each cycle, FMCBA establishes the stand's fire/fuels context:
#   - the cover type (the species carrying the most basal area),
#   - percent canopy cover (from per-tree crown areas),
#   - the live herb/shrub surface fuel (by FFE forest type), and
#   - in the first FFE year, the dead surface fuel pools, split into decay classes
#     by each species' share of stand basal area.
# Results live in the per-stand `FireState` (no globals). Gated on `fire.active`, so
# it is a no-op for non-FFE stands; FMCBA itself changes nothing the `.sum` reports
# (it feeds the fire-behavior / consumption chunks F5–F8).
# =============================================================================

# FMCBA's per-record species basal area (the cover type COVTYP = most BA; PRCL = TBA/TOTBA splits the initial dead fuel
# into decay classes). The builds use three forms: bc/cs/em/ie/kt/sn `TBA += BA1·FMPROB`, BA1 = 3.14159·(DBH/24)·(DBH/24);
# ak/bm/ca/ci/cr/ec/ls/nc/ne/oc/on/op/pn/tt/ut/wc/ws `FMTBA += FMPROB·DBH·DBH·0.0054542`; so `TBA += FMPROB·5.454153E-03·
# DBH**2` (integer power = DBH·DBH). jl used the first form everywhere but AK (BM 12827438010497: 2005 initial fuel 1 ULP).
_fmcba_pi_form(v) = v isa Southern || v isa CentralStates || v isa InlandEmpire || v isa EasternMontana ||
                    v isa Kootenai || v isa BritishColumbia
@inline function _fmcba_tba(v, p::Float32, d::Float32)::Float32
    _fmcba_pi_form(v) && return 3.14159f0 * (d / 24f0) * (d / 24f0) * p
    v isa SouthCentralOregon && return p * 5.454153f-3 * (d * d)
    return p * d * d * 0.0054542f0
end

"""
    fmcba!(s) -> StandState

Update the stand's `FireState` cover type, percent cover, big DBH, live fuels, and
(first FFE year) dead-fuel pools (FMCBA, fmcba.f). No-op unless FFE is active.
"""
# so/fmcba.f:893-990 (the SO-FFE snag block FMVINIT leaves "unset" = -1, applied every FMCBA with IF(X.LT.0)): on the
# California forests (KODFOR 500-599, 701) ALLDWN 100 (juniper 150), FALLX 1.235/0.882/0.687 by species group, DECAYX 999
# (hard snags never soften), HTX 1/0/1/0 (height loss stops at 50%), PBSOFT 1.0 / PBSMAL 0.9; on the Oregon forests the
# fire_species_props.csv values (ALLDWN 110/90/100, FALLX = DECAYX = 1), HTX 1 for all four (FMR6HTLS), PBSOFT=PBSMAL=0.
# A SNAGFALL/SNAGDCAY/SNAGBRK/SNAGPBN keyword value is kept (the per-species override dicts / a non-negative PB*).
# MEASURED FVSso_g16 15364795010497 (forest 505) FVS_SnagDet 2020: the inventory SH snags 23.5->17.79 ft live, jl held them.
function _so_snag_params!(s::StandState)
    p = s.fire.params
    kodfor = Int(s.plot.user_forest_code)
    ca = (500 <= kodfor < 600) || kodfor == 701
    nsp = nspecies(s.variant)
    @inbounds for sp in 1:nsp
        k = Int32(sp)
        if ca
            fx, ad = sp in (3, 4, 8, 9, 12, 13, 14, 15, 17, 32) ? (0.882f0, 100f0) :
                     sp in (6, 18, 20) ? (0.687f0, 100f0) : sp == 11 ? (0.687f0, 150f0) : (1.235f0, 100f0)
            haskey(p.snag_fallx_ovr, k) || (p.snag_fallx_ovr[k] = fx)
            haskey(p.snag_alldwn_ovr, k) || (p.snag_alldwn_ovr[k] = ad)
            haskey(p.snag_decayx_ovr, k) || (p.snag_decayx_ovr[k] = 999f0)
            haskey(p.snag_htx, k) || (p.snag_htx[k] = (1f0, 0f0, 1f0, 0f0))
        else
            haskey(p.snag_htx, k) || (p.snag_htx[k] = (1f0, 1f0, 1f0, 1f0))
        end
    end
    p.pb_soft < 0f0 && (p.pb_soft = ca ? 1f0 : 0f0)
    p.pb_smal < 0f0 && (p.pb_smal = ca ? 0.9f0 : 0f0)
    return
end

function fmcba!(s::StandState; load_dead::Bool = true, vtrip::Bool = false)
    fs = s.fire
    (fs === nothing || !fs.active) && return s
    t = s.trees; coef = s.coef
    # FFE is an UNPORTED extension for the ORGANON variants (OregonCoast/Olympic): their species carry
    # no FFE crown-biomass / fuel-coefficient tables in `coef.species` (they hold sp_dbh_min in CONTROL,
    # for volume gating, not `:dbh_min` in coef). Fail with an actionable message rather than a cryptic
    # KeyError deep in the fuel loop below. Growth+volume for these variants is unaffected (FFE only
    # activates under FMIn/SIMFIRE/… keywords). See memory fvsjl-oc-organon-blm-volume.
    if !haskey(coef.species, :dbh_min)
        error("FFE (Fire & Fuels Extension) is not ported for the $(nameof(typeof(s.variant))) " *
              "variant: its ORGANON species carry no FFE crown-biomass/fuel coefficients. Remove the " *
              "FFE keywords (FMIn/SIMFIRE/SNAGINIT/…) from this stand, or run growth+volume only.")
    end
    nsp = length(coef_col(coef, :dbh_min))            # MAXSP
    s.variant isa SouthCentralOregon && _so_snag_params!(s)

    # live herb/shrub fuels are re-set every year. NE (ne/fmcba.f:68) uses a flat constant
    # FULIV=(0.31,0.31) for all forest types; SN uses the FFE-forest-type table with a
    # coastal-plain/piedmont/mountain rough-age × site-index override (FULIV2) where it applies.
    if s.variant isa Northeast
        fs.flive = (0.31f0, 0.31f0)
    elseif s.variant isa CentralStates
        # CS uses the flat FULIV table by FFE forest type (cs/fmcba.f:113); no SN FULIV2 override.
        fs.flive = ffe_live_fuel_loading(coef, ffe_forest_type(s))
    elseif s.variant isa LakeStates
        # LS uses the FULIV table indexed by (IFFEFT, ISZCL) — ls/fmcba.f:138; no SN FULIV2 override.
        fs.flive = ls_live_fuel_loading(s)
    elseif s.variant isa CentralRockies || s.variant isa InlandEmpire || s.variant isa Kootenai ||
           s.variant isa EasternMontana || s.variant isa CentralIdaho ||
           s.variant isa Teton || s.variant isa Utah || s.variant isa BlueMountains ||
           s.variant isa Klamath || s.variant isa WestCascades || s.variant isa PacificNorthwest ||
           s.variant isa EastCascades || s.variant isa SouthCentralOregon ||
           s.variant isa OregonCoast || s.variant isa Olympic || s.variant isa SoutheastAlaska
        # Western (CR/IE/KT/EM/…): live fuel = FULIVE/FULIVI[COVTYP] interpolated by PERCOV — DEFERRED to after
        # the cover-type block below (needs COVTYP + PERCOV). NC additionally needs the top-2 COVCA/COVCAWT.
        # Placeholder here.
        fs.flive = (0f0, 0f0)
    else
        ovr = ffe_live_fuel_override(s)
        fs.flive = ovr === nothing ? ffe_live_fuel_loading(coef, ffe_forest_type(s)) : ovr
    end

    # per-species basal area, total crown area (for percent cover), and the big DBH
    tba = zeros(Float32, nsp)
    totcra = 0f0
    # CR forest-grown crown width is cr_cwcalc (cwcalc.f IWHO=0 — the same CRWDTH FMCBA/FMSSTAGE use), NOT the
    # generic crown_width, which returns the 0.5 default for every CR species ⇒ near-zero crown area ⇒ PERCOV≈0.
    _cr_fm = s.variant isa CentralRockies
    _bm_fm = s.variant isa BlueMountains        # BM CRWDTH via bm_cwcalc (cwcalc.f BMMAP + R6 BF, national library)
    _nc_fm = s.variant isa Klamath              # NC CRWDTH via nc_cwcalc (NCMAP western Bechtold/Crookston library)
    _ws_fm = s.variant isa WestSierra           # WS CRWDTH via ws_cwcalc (WSMAP; R5 forest 511 ⇒ BF=1, same forms as NC)
    _ca_fm = s.variant isa CentralCalifornia    # CA CRWDTH via ca_cwcalc (CAMAP Crookston R6; forest 610 IFOR>5)
    _wc_fm = s.variant isa WestCascades         # WC CRWDTH via wc_cwcalc (WCMAP Crookston R6; forest 618 IFOR=6)
    _pn_fm = s.variant isa PacificNorthwest      # PN CRWDTH via pn_cwcalc (same WCMAP; forest 612 SIUSLAW BF)
    _ec_fm = s.variant isa EastCascades          # EC CRWDTH via ec_cwcalc (ECMAP Crookston R6; forest 608 OKANOGAN BF)
    _so_fm = s.variant isa SouthCentralOregon     # SO CRWDTH via so_cwcalc (SOMAP Crookston R6; forest 601 DESCHUTES BF)
    _oc_fm = s.variant isa OregonCoast            # OC CRWDTH via oc_cwcalc (OCMAP Crookston R6/R1; forest 711 BLM Medford→610 Rogue River BF)
    _op_fm = s.variant isa Olympic                # OP CRWDTH via op_cwcalc (OPMAP Crookston R6; forest 708 BLM Salem→606 Mt Hood BF)
    # IE/KT FMCBA reads CRWDTH(I) filled by ie/cwidth.f → cwcalc.f (IEMAP/KTMAP, IWHO=0) — the same forest-grown
    # value FVS_TreeList reports (`tree_crwdth`). ccfcal MODE=2 (the B1·exp(…) form jl used) agrees for most trees
    # but not the small-tree forms: MEASURED FVSie_g16 11855985010690 2026 FMCBA, ES/AF/GF seedlings CRWDTH 0.5/0.55
    # vs MODE=2 1.2/1.09 ⇒ TOTCRA 32619 vs 32643 ⇒ PERCOV ⇒ FLIVE(2) shrub load 0.16% high. (Without any western
    # crown width IE fell to the eastern `crown_width` 0.5 default ⇒ PERCOV≈1 ⇒ a spurious passive crown fire.)
    # EM CRWDTH via em_cwcalc (em/cwidth.f → em/cwcalc.f IWHO=0: the western Crookston/Bechtold library with
    # BAREA=BA, EL=ELEV, HI=Hopkins index, cwcalc.f:563-576). Without it EM fell to the generic `crown_width`
    # (0.5 ft default) ⇒ PERCOV 0.55 vs live 38.71 ⇒ WMULT 0.5 vs 0.197 ⇒ FWIND 5.0 vs 1.97 + wrong FLIVE.
    # BAREA = the stand BA in every cycle (MEASURED vs FVSem_g16 DEBUG FMCBA: PERCOV cyc1 36.7793 = live 36.78
    # with BA; the NC-style load-time BAREA=1 gives 38.07), unlike NC/WS/CA/OC/OP's cycle-1 clamp.
    _em_fm = s.variant isa EasternMontana
    # CI/TT/UT FMCBA read the common CRWDTH(I) (ci/fmcba.f:340, tt:290, ut:312) — the same forest-grown cwcalc value
    # FVS_TreeList reports (`tree_crwdth`); they fell to the generic crown_width like EM did.
    _citu_fm = s.variant isa CentralIdaho || s.variant isa Teton || s.variant isa Utah ||
               s.variant isa InlandEmpire || s.variant isa Kootenai ||
               s.variant isa SoutheastAlaska   # ak/fmcba.f:195 CWIDTH=CRWDTH(I) (ak_cwcalc, the TreeList CrWidth)
    _ak_fm = s.variant isa SoutheastAlaska
    _west_cw = _cr_fm || _bm_fm || _nc_fm || _ws_fm || _ca_fm || _wc_fm || _pn_fm || _ec_fm || _so_fm || _oc_fm || _op_fm
    _cr_ba = _west_cw ? s.plot.basal_area : 0f0
    # NC CRWDTH (base cwidth.f→cwcalc.f) is computed by CWIDTH at LOAD time, BEFORE the stand BA is
    # accumulated ⇒ the R6-Crookston BAREA term hits cwcalc.f:859 `IF(BAREA.LE.1.) BAREA=1.` (BA=0→1).
    # fmcba in the FIRST FFE cycle reads those load-time CRWDTH; later cycles read the end-of-cycle UPDATE
    # (actual BA). MEASURED vs FVSnc_clean DEBUG FMCBA cyc1: PERCOV=39.01 matches BAREA=1 (jl 39.72), NOT the
    # dense stand BA (jl 44.29). Mirror the load-time clamp for cycle 1.
    # WC differs from NC/WS/CA: its CRWDTH (wc/cwcalc.f) is computed AFTER the stand BA is accumulated, so
    # the Crookston BAREA term uses the actual FFE stand BA (~85 on wct01), NOT the BAREA=1 load-time clamp.
    # MEASURED vs FVSwc_clean cyc0: WF CW 10.31 needs (BA+1)^e with BA≈85 (BA=1 gives 8.74). So WC uses _cr_ba.
    _nc_ba = ((_nc_fm || _ws_fm || _ca_fm || _oc_fm || _op_fm) && s.control.cycle <= Int32(1)) ? 1f0 : _cr_ba

    _cr_el = _west_cw ? s.plot.elevation : 0f0
    _cr_hi = _west_cw ? _cr_hopkins(s.plot.latitude, s.plot.longitude, s.plot.elevation) : 0f0
    _bm_kf = _bm_fm ? bm_kodfor_remap(Int(s.plot.user_forest_code)) : 0   # BM CRWDTH forest BF key (post-FORKOD)
    # CWIDTH=CRWDTH(I) (ie/fmcba.f:244, em/fmcba.f:253): the value CWIDTH last stored (`stored_crwdth`, cwidth.f at load and
    # gradd.f:254), carried by TRIPLE — not one recomputed from a SIMFIRE seam's grown small trees or a thin's residual BA
    # (MEASURED FVSie_g16 11855985010690 2016 burn: records 12/116 CRWDTH 1.0721917 both, jl from the seam DBH 1.0777394 ⇒
    # TOTCRA 30163.06 vs 30151.213 ⇒ PERCOV/WMULT ⇒ Midflame_Wind 1.66760 vs 1.66796, fire kill).
    _old_cw = _stored_crwdth(s.variant)
    # nc/so/ca fmcba.f read CWIDTH = CRWDTH(I), and cwcalc.f routes their Region-5 forests to R5CRWD (a function of
    # sp/D/H only): NC IFOR ≤ 3 or 5 (cwcalc.f:382), SO IFOR 4-9 (:376), CA IFOR ≤ 5 (:385) — the same CRWDTH FVS_TreeList
    # reports (_forest_crwdth). jl ran the R6 Crookston kernels on every forest (MEASURED FVSnc_g16 23660512010900,
    # forest 508→510: inventory PERCOV live 15.09 = the TreeList CrWidth; jl 10.0).
    _ifor = Int(s.plot.forest_idx)
    _r5cw = (_nc_fm && (_ifor <= 3 || _ifor == 5)) || (_so_fm && 4 <= _ifor <= 9) || (_ca_fm && _ifor <= 5)
    cwrec = zeros(Float32, t.n)
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        sp = Int(t.species[i]); d = t.dbh[i]
        d > fs.bigdbh && (fs.bigdbh = d)
        cw = _old_cw ? stored_crwdth(s, i) :
             _r5cw ? _forest_crwdth(s, sp, d, t.height[i], t.crown_pct[i]) :
             _cr_fm ? cr_cwcalc(sp, d, t.height[i], Float32(t.crown_pct[i]), _cr_ba, _cr_el, _cr_hi) :
             _bm_fm ? (t.ffe_oldht[i] > 0f0 ?
                       # CRWDTH(I) as CWIDTH last set it (gradd.f:254 end of cycle / fvs.f:207 load), carried by TRIPLE:
                       # the dims of the FMOLDC-time snapshot, not REGENT's grown small trees (bm/fmcba.f:196 CRWDTH(I))
                       bm_cwcalc(sp, t.ffe_olddbh[i], t.ffe_oldht[i], t.ffe_oldcr[i], _cr_ba, _cr_el, _cr_hi; kodfor = _bm_kf) :
                       bm_cwcalc(sp, d, t.height[i], Float32(t.crown_pct[i]), _cr_ba, _cr_el, _cr_hi; kodfor = _bm_kf)) :
             _nc_fm ? nc_cwcalc(sp, d, t.height[i], Float32(t.crown_pct[i]), _nc_ba, _cr_el, _cr_hi) :
             _ws_fm ? ws_r5crwd(sp, d, t.height[i]) :   # WS: R5CRWD (ws/r5crwd.f), function of sp/D/H only
             _ca_fm ? ca_cwcalc(sp, d, t.height[i], Float32(t.crown_pct[i]), _nc_ba, _cr_el, _cr_hi) :  # CA R6 Crookston (ca/cwcalc.f CAMAP)
             _wc_fm ? wc_cwcalc(sp, d, t.height[i], Float32(t.crown_pct[i]), _nc_ba, _cr_el, _cr_hi; kodfor = Int(s.plot.user_forest_code)) :  # WC R6 Crookston (wc/cwcalc.f WCMAP)
             _pn_fm ? pn_cwcalc(sp, d, t.height[i], Float32(t.crown_pct[i]), _nc_ba, _cr_el, _cr_hi; kodfor = Int(s.plot.user_forest_code)) :  # PN = wc/cwcalc.f (byte-identical) with PN's KODFOR
             _ec_fm ? ec_cwcalc(sp, d, t.height[i], Float32(t.crown_pct[i]), _cr_ba, _cr_el, _cr_hi; kodfor = Int(s.plot.user_forest_code)) :  # EC R6 Crookston (ec/cwcalc.f ECMAP; forest-608 BF)
             _so_fm ? so_cwcalc(sp, d, t.height[i], Float32(t.crown_pct[i]), _cr_ba, _cr_el, _cr_hi) :  # SO R6 Crookston (so/cwcalc.f SOMAP; forest-601 BF)
             _oc_fm ? oc_cwcalc(sp, d, t.height[i], Float32(t.crown_pct[i]), _nc_ba, _cr_el, _cr_hi) :  # OC R6 Crookston (oc/cwcalc.f OCMAP; forest-711→610 BF)
             _op_fm ? op_cwcalc(sp, d, t.height[i], Float32(t.crown_pct[i]), _nc_ba, _cr_el, _cr_hi) :  # OP R6 Crookston (op/cwcalc.f OPMAP; forest-708→606 BF)
             _em_fm ? em_cwcalc(sp, d, t.height[i], Float32(t.crown_pct[i]), s.plot.basal_area, s.plot.elevation,
                                _cr_hopkins(s.plot.latitude, s.plot.longitude, s.plot.elevation)) :   # EM (em/cwcalc.f)
             _citu_fm ? tree_crwdth(s, sp, d, t.height[i], t.crown_pct[i]) :   # CI/TT/UT: CWIDTH=CRWDTH(I) (ci,tt,ut/fmcba.f)
             crown_width(coef, s.species.code2[sp], d, t.height[i], Float32(t.crown_pct[i]), 0,
                         s.plot.latitude, s.plot.longitude, s.plot.elevation)   # forest-grown (CWCALC iwho=0)
        cwrec[i] = cw
    end
    # fmcba.f:189-203 DO I=1,ITRN: TBA(KSP) += BA1·FMPROB(I), TOTCRA += CAREA·FMPROB(I), in FVS's record order. FMMAIN runs
    # after TRIPLE, so in a tripling cycle the walk is the tripled list (originals ×.60, then each record's ×.25/×.15 copies,
    # same DBH/CRWDTH) — jl's list is still untripled here (MEASURED FVSsn private FMCBA trace, 200267456010854 2002: sp74
    # TBA 28.677080 live vs 28.677082 summing the untripled records ⇒ the initial 6-12" fuel 1.0099999 vs 1.01).
    _fm_record_walk(t, vtrip) do i, pr
        pr > 0f0 || return
        sp = Int(t.species[i]); d = t.dbh[i]
        tba[sp] += _fmcba_tba(s.variant, pr, d)
        totcra += 3.1415927f0 * cwrec[i] * cwrec[i] / 4f0 * pr
    end

    # cover type = the species with the most basal area; total BA for the decay split
    bamost = 0f0; covtyp = Int32(0); totba = 0f0
    @inbounds for ksp in 1:nsp
        tba[ksp] > bamost && (bamost = tba[ksp]; covtyp = Int32(ksp))
        totba += tba[ksp]
    end
    # NC (California/westside) uses the TOP TWO cover-type species for the live + initial-dead fuel pools
    # (nc/fmcba.f:298-310): RDPSRT the per-species BA descending → ICT; COVCA(1..2)=ICT(1..2); COVCAWT(j)=
    # FMTBA(ICT(j))/Σ_{i=1,2}FMTBA(ICT(i)). COVTYP is ICT(1) when its BA>0.001 (faithful RDPSRT tie-break).
    covca = (0, 0); covcawt = (0f0, 0f0)
    if _ak_fm
        # ak/fmcba.f:216-217: RDPSRT(MAXSP,FMTBA,ICT) descending; COVTYP=ICT(1) when its BA > 0.001.
        ict = collect(1:nsp)
        rdpsrt!(nsp, tba, ict, true)
        covtyp = tba[ict[1]] > 0.001f0 ? Int32(ict[1]) : Int32(0)
    end
    if s.variant isa Klamath || s.variant isa WestSierra || s.variant isa CentralCalifornia || s.variant isa OregonCoast || s.variant isa Olympic
        ict = collect(1:nsp)
        rdpsrt!(nsp, tba, ict, true)                     # descending indirect sort on tba → ICT
        covtyp = tba[ict[1]] > 0.001f0 ? Int32(ict[1]) : Int32(0)
        xx = 0f0
        @inbounds for i in 1:2
            tba[ict[i]] > 0f0 && (xx += tba[ict[i]])
        end
        covca = (Int(ict[1]), Int(ict[2]))
        covcawt = xx > 0.001f0 ? (tba[ict[1]] / xx, tba[ict[2]] / xx) : (0f0, 0f0)
    end
    # No basal area: FVS fmcba.f sets COVTYP to a VARIANT-SPECIFIC default cover-type species the first
    # year (SN 75, NE 1, CS 48, LS 3 red pine — fmcba.f "COVTYP.EQ.0" block), else keeps the previous
    # cover type (OLDCOVTYP). The default is a SPECIES index feeding DKRCLS/FULIVI; the hardcoded 75 was
    # SN-only and indexed a 0 decay-class on LS (dkr_cls[75]=0) ⇒ an @inbounds OOB write in the cwd fill.
    if covtyp == 0
        covtyp = fs.covtyp != Int32(0) ? fs.covtyp :
                 # IE/KT: bare stand ⇒ COVINI(ITYPE) seral cover species (ie/fmcba.f:279), NOT a fixed species.
                 (s.variant isa InlandEmpire || s.variant isa Kootenai) ? Int32(ie_covini_default(Int(s.plot.habitat_input))) :
                 s.variant isa EasternMontana ? Int32(3)  :   # EM bare-stand default: DF (em/fmcba.f covtyp=3 path)
                 # CI/TT/UT/WC/PN/EC bare stand: COVINI(ITYPE) (the seral cover of the habitat), else the variant's
                 # "NO VALID HABITAT" default (ci:377-387 ICINDX→DF 3; tt:327-337 / ut:351-361 ITYPE→LP 7;
                 # wc:451-462 / pn:425-436 COVINI6(ITYPE)→DF 16; ec:422-431 ITYPE→DF 3). Tables: covini_tables.jl.
                 s.variant isa CentralIdaho ? _covini(CI_COVINI, Int(s.plot.habitat_input), 3) :   # ICINDX (ci/fmcba.f:377-387)
                 s.variant isa Teton ? _covini(TT_COVINI, Int(s.plot.habitat_input), 7) :          # ITYPE (tt/habtyp.f)
                 s.variant isa Utah ? _covini(UT_COVINI, Int(s.plot.habitat_input), 7) :           # ITYPE (ut/habtyp.f)
                 # BM bare stand: COVINI(ITYPE) (bm/fmcba.f:296-300), ITYPE = the PCOML index habtyp.f resolved (79 = its
                 # CWG113 default ⇒ COVINI 4 = grand fir, live herb/shrub 0.30/2.0 — the fixed DF (0.4/2.0) over-loaded FLIVE).
                 s.variant isa BlueMountains ? Int32(_BM_COVINI[(1 <= s.plot.habitat_code <= 92) ? Int(s.plot.habitat_code) : 79]) :
                 s.variant isa Northeast     ? Int32(1)  :
                 s.variant isa CentralStates ? Int32(48) :
                 s.variant isa LakeStates    ? Int32(3)  :
                 # NC bare stand: COVINI5(ITYPE)/COVINI6(ITYPE) by habitat (nc/fmcba.f:337-355); the full
                 # R5/R6 habitat→cover maps are a bare-stand-only path not exercised by nct01 (trees present) —
                 # the fmcba.f "no valid habitat" fallback is Douglas-fir (3), used here until those maps port.
                 # NC (nc/fmcba.f:337-355): R5 forests (KODFOR 5xx or ≥ 705) COVINI5(ITYPE), else COVINI6(ITYPE); DF 3 default.
                 s.variant isa Klamath ?
                     (let kf = Int(s.plot.user_forest_code)
                          ((500 <= kf < 600) || kf >= 705) ? _covini(NC_COVINI5, Int(s.plot.habitat_input), 3) :
                                                             _covini(NC_COVINI6, Int(s.plot.habitat_input), 3)
                      end) :
                 s.variant isa WestCascades ? _covini(WC_COVINI6, Int(s.plot.habitat_input), 16) :
                 s.variant isa PacificNorthwest ? _covini(PN_COVINI6, Int(s.plot.habitat_input), 16) :
                 s.variant isa EastCascades ? _covini(EC_COVINI, Int(s.plot.habitat_input), 3) :
                 # SO (so/fmcba.f:608-616): COVINI(ITYPE) when ITYPE > 0 on an R6 forest (600-699, 799), else PP 10.
                 s.variant isa SouthCentralOregon ?
                     (let kf = Int(s.plot.user_forest_code)
                          ((600 <= kf < 700) || kf == 799) ? _covini(SO_COVINI, Int(s.plot.habitat_input), 10) : Int32(10)
                      end) :
                 s.variant isa SoutheastAlaska ? Int32(11) :   # AK: western hemlock at IY(1) (ak/fmcba.f:236-242), else OLDCOVTYP
                 # WS/CA/OC/OP "NO VALID HABITAT" defaults (ws/fmcba.f:509 PP 10; ca/oc fmcba.f:543 7; op/fmcba.f:436 DF
                 # 16 after COVINI6(ITYPE)). OC's COVINI waits on its habtyp ITYPE port.
                 s.variant isa WestSierra ? _covini(WS_COVINI, Int(s.plot.habitat_input), 10) :     # ws/fmcba.f:500-509
                 # CA (ca/fmcba.f:527-543): IFOR ≥ 6 (R6/BLM) COVINI6(ITYPE), else COVINI5(ITYPE); 7 default.
                 s.variant isa CentralCalifornia ?
                     (Int(s.plot.forest_idx) >= 6 ? _covini(CA_COVINI6, Int(s.plot.habitat_input), 7) :
                                                    _covini(CA_COVINI5, Int(s.plot.habitat_input), 7)) :
                 s.variant isa OregonCoast ? Int32(7) :
                 s.variant isa Olympic ? _covini(OP_COVINI6, Int(s.plot.habitat_input), 16) :
                 # CR (cr/fmcba.f:417-433): COVINI2(ITYPE) on Region-2 forests, COVINI3(ITYPE) on Region 3, else LP 11.
                 s.variant isa CentralRockies ?
                     (Int(s.plot.user_forest_code) ÷ 100 == 2 ? _covini(CR_COVINI2, Int(s.plot.habitat_input), 11) :
                      Int(s.plot.user_forest_code) ÷ 100 == 3 ? _covini(CR_COVINI3, Int(s.plot.habitat_input), 11) :
                      Int32(11)) : Int32(75)
        # CA-FFE top-2 variants (nc:359-360, ws:514-515, ca/oc:548-549): a bare stand carries the ONE cover type,
        # COVCA(1)=COVTYP with weight COVCAWT(1)=1 (COVCA(2)/COVCAWT(2) stay 0 from the FMCBA entry reset). jl kept
        # COVCAWT=(0,0) ⇒ zero live AND initial dead fuel on every bare NC/WS/CA/OC stand. (OP's fmcba.f has no
        # COVCA; its single-COVTYP load equals this one-weight pair.)
        if s.variant isa Klamath || s.variant isa WestSierra || s.variant isa CentralCalifornia ||
           s.variant isa OregonCoast || s.variant isa Olympic
            covca = (Int(covtyp), 0); covcawt = (1f0, 0f0)
        end
    end
    fs.covtyp = covtyp
    fs.percov = (1f0 - fexp(-totcra / 43560f0)) * 100f0   # fmcba.f:228-229 PERCOV=(1.0-EXP(-TOTCRA/43560.))·100 — EXP is expf
    # SO (FCCS/Ottmar) resolves BOTH the live (herb,shrub) and the 11-class dead pool from one COVRINI→FUELINI
    # lookup keyed by COVTYP, the FMSSTAGE structural stage ISSX (with the PERCOV≥60 SE-open→SE-closed bump),
    # and the logging-history model index (so/fmcba.f:639-726). Computed here (once COVTYP+PERCOV are known);
    # the dead pool is reused by the first-year dead-fuel block below.
    so_ini = nothing
    if s.variant isa SouthCentralOregon
        _cls = Int(structure_class(s).class)               # FMSSTAGE IFMST
        _issx = _cls; _logmod = 1
        if _issx == 0
            _issx = 1; _logmod = 2                          # bare/unclassified ⇒ SI + regenerated logging model
        elseif _issx == 2 && fs.percov >= 60f0
            _issx = 3                                       # SE open → SE closed (Ottmar) when cover ≥ 60%
        end
        so_ini = so_fuel_ini(Int(covtyp), _issx, _logmod)
        fs.flive = (so_ini[1], so_ini[2])
    end
    # Western live fuel now that COVTYP + PERCOV are known (fmcba.f:443-449 / ie:283-289)
    s.variant isa CentralRockies && (fs.flive = cr_live_fuel_loading(Int(covtyp), fs.percov))
    (s.variant isa InlandEmpire || s.variant isa Kootenai) && (fs.flive = ie_live_fuel_loading(Int(covtyp), fs.percov))
    s.variant isa EasternMontana && (fs.flive = em_live_fuel_loading(Int(covtyp), fs.percov))
    s.variant isa CentralIdaho && (fs.flive = ci_live_fuel_loading(Int(covtyp), fs.percov))
    s.variant isa Teton && (fs.flive = tt_live_fuel_loading(Int(covtyp), fs.percov))
    s.variant isa Utah && (fs.flive = ut_live_fuel_loading(Int(covtyp), fs.percov))
    s.variant isa BlueMountains && (fs.flive = bm_live_fuel_loading(Int(covtyp), fs.percov))
    s.variant isa Klamath && (fs.flive = nc_live_fuel_loading(covca, covcawt, fs.percov))   # top-2 (nc/fmcba.f:369-378)
    s.variant isa WestSierra && (fs.flive = ws_live_fuel_loading(covca, covcawt, fs.percov))  # top-2 (ws/fmcba.f:519-533)
    s.variant isa CentralCalifornia && (fs.flive = ca_live_fuel_loading(covca, covcawt, fs.percov))  # top-2 (ca/fmcba.f)
    s.variant isa OregonCoast && (fs.flive = oc_live_fuel_loading(covca, covcawt, fs.percov))  # top-2 (oc/fmcba.f ORGANON 50-sp)
    s.variant isa Olympic && (fs.flive = op_live_fuel_loading(covca, covcawt, fs.percov))  # top-2 (op/fmcba.f NWO 39-sp)
    s.variant isa WestCascades && (fs.flive = wc_live_fuel_loading(Int(covtyp), fs.percov))  # single COVTYP (wc/fmcba.f:476-480)
    s.variant isa PacificNorthwest && (fs.flive = pn_live_fuel_loading(Int(covtyp), fs.percov))  # single COVTYP (pn/fmcba.f)
    s.variant isa EastCascades && (fs.flive = ec_live_fuel_loading(Int(covtyp), fs.percov))  # single COVTYP (ec/fmcba.f)
    _ak_fm && (fs.flive = ak_live_fuel_loading(Int(covtyp), fs.percov))   # ak/fmcba.f:253-257 FULIVI/FULIVE by COVTYP

    # dead fuels: loaded once (first FFE year), distributed into decay classes by the species BA share
    # (fmcba.f:375-393). The "hard" (J=2) column comes from ffe_dead_fuel_loading; the "soft" (J=1) column
    # is 0 by default. FUELINIT (hard) / FUELSOFT (soft) override per-size-class values (STFUEL, fmcba.f:320-371).
    # IDC = each species' decay-rate class (DKRCLS).
    if load_dead && !fs.fuels_init
        deffuel = s.variant isa Northeast ? ne_dead_fuel_loading(s) :
                  s.variant isa CentralStates ? cs_dead_fuel_loading(coef, Int(s.plot.forest_type)) :
                  s.variant isa LakeStates ? ls_dead_fuel_loading(s) :
                  s.variant isa CentralRockies ? cr_dead_fuel_loading(Int(covtyp), fs.percov) :  # FUINIE/FUINII × PERCOV
                  (s.variant isa InlandEmpire || s.variant isa Kootenai) ? ie_dead_fuel_loading(Int(covtyp), fs.percov) :
                  s.variant isa EasternMontana ? em_dead_fuel_loading(Int(covtyp), fs.percov) :
                  s.variant isa CentralIdaho ? ci_dead_fuel_loading(Int(covtyp), fs.percov) :
                  s.variant isa Teton ? tt_dead_fuel_loading(Int(covtyp), fs.percov) :
                  s.variant isa Utah ? ut_dead_fuel_loading(Int(covtyp), fs.percov) :
                  s.variant isa BlueMountains ? bm_dead_fuel_loading(Int(covtyp), fs.percov) :
                  s.variant isa Klamath ? nc_dead_fuel_loading(covca, covcawt, fs.percov) :  # top-2 (nc/fmcba.f:421-431)
                  s.variant isa WestSierra ? ws_dead_fuel_loading(covca, covcawt, fs.percov) :  # top-2 (ws/fmcba.f:587-597)
                  s.variant isa CentralCalifornia ? ca_dead_fuel_loading(covca, covcawt, fs.percov) :
                  s.variant isa OregonCoast ? oc_dead_fuel_loading(covca, covcawt, fs.percov) :  # top-2 (oc/fmcba.f ORGANON 50-sp)
                  s.variant isa Olympic ? op_dead_fuel_loading(covca, covcawt, fs.percov) :  # top-2 (op/fmcba.f NWO 39-sp)
                  s.variant isa WestCascades ? wc_dead_fuel_loading(Int(covtyp), fs.percov) :  # single COVTYP (wc/fmcba.f:528-533)
                  s.variant isa PacificNorthwest ? pn_dead_fuel_loading(Int(covtyp), fs.percov) :  # single COVTYP (pn/fmcba.f)
                  s.variant isa EastCascades ? ec_dead_fuel_loading(Int(covtyp), fs.percov) :  # single COVTYP (ec/fmcba.f)
                  _ak_fm ? ak_dead_fuel_loading(Int(s.plot.forest_type)) :   # ak/fmcba.f:270-324 FUINI(FTDEADFU(IFORTP))
                  s.variant isa SouthCentralOregon ? Float32[so_ini[3]...] :  # FCCS/Ottmar STFUEL (so/fmcba.f)
                  ffe_dead_fuel_loading(coef, Int(s.plot.forest_type))
        # Seed the STFUEL override from FIA-DB measured fuel loadings (FVS_STANDINIT FUEL_* → dbsstandin.f
        # FUELINIT, read into plot.ffe_fuel_*) when present AND no explicit FUELINIT/FUELSOFT keyword already set
        # them (the keyword takes precedence, matching FVS where a later-scheduled FUELINIT overrides the DB one).
        isempty(fs.params.stfuel_hard) && !isempty(s.plot.ffe_fuel_hard) && (fs.params.stfuel_hard = copy(s.plot.ffe_fuel_hard))
        isempty(fs.params.stfuel_soft) && !isempty(s.plot.ffe_fuel_soft) && (fs.params.stfuel_soft = copy(s.plot.ffe_fuel_soft))
        ovh = fs.params.stfuel_hard; ovs = fs.params.stfuel_soft
        # CWD(1,ISZ,J,IDC) = CWD + ADD (fmcba.f:597-613): the initial load ADDS to whatever the cut phase already
        # booked this year (FMSCUT crown slash, YARDLOSS downed boles — cuts.f precedes FMMAIN), it does not reset.
        @inbounds for isz in 1:11
            sh = (length(ovh) >= isz && ovh[isz] >= 0f0) ? ovh[isz] : deffuel[isz]   # FUELINIT hard
            ss = (length(ovs) >= isz && ovs[isz] >= 0f0) ? ovs[isz] : 0f0            # FUELSOFT soft
            if totba > 0f0
                for ksp in 1:nsp
                    tba[ksp] > 0f0 || continue
                    idc = ffe_dkr_cls(s, ksp)               # FUELPOOL-overridable decay-rate class
                    prcl = tba[ksp] / totba
                    fs.cwd[isz, 2, idc] += prcl * sh
                    ss > 0f0 && (fs.cwd[isz, 1, idc] += prcl * ss)
                end
            else
                idc = ffe_dkr_cls(s, covtyp)               # FUELPOOL-overridable decay-rate class
                fs.cwd[isz, 2, idc] += sh
                ss > 0f0 && (fs.cwd[isz, 1, idc] += ss)
            end
        end
        # BM decay-rate habitat adjustment (bm/fmcba.f:333-368): applied ONCE at the first FFE year, and only
        # when the user has not set the decay rates with FuelDcay/FuelMult (fs.params.dkr still empty). Store
        # the adjusted matrix into params.dkr so every subsequent fmcwd! reads it (mirrors the persisted DKR).
        if s.variant isa BlueMountains && size(fs.params.dkr, 1) != 11
            fs.params.dkr = bm_adjusted_dkr(bm_itype(Int(s.plot.habitat_code)))   # habtyp.f ITYPE (79 if unmatched)
        end
        # NC decay-rate DCYMLT (nc/fmcba.f:395-414): scale the NC base DKR by the Dunning-code/site-index
        # multiplier at the first FFE year (when the user hasn't set FuelDcay ⇒ params.dkr still empty).
        if s.variant isa Klamath && size(fs.params.dkr, 1) != 11
            _ss = Int(s.plot.site_species)
            _si = (1 <= _ss <= length(s.plot.sp_site_index)) ? s.plot.sp_site_index[_ss] : 0f0
            fs.params.dkr = nc_adjusted_dkr(_si)
        end
        # EC decay-rate habitat adjustment (ec/fmcba.f:457-491): scale the EC base DKR by DKRADJ(TEMP,MOIST,K)
        # at the first FFE year (when the user hasn't set FuelDcay ⇒ params.dkr still empty).
        if s.variant isa EastCascades && size(fs.params.dkr, 1) != 11
            fs.params.dkr = ec_adjusted_dkr(Int(s.plot.habitat_input))
        end
        # WC / PN / OP decay-rate habitat adjustment (wc/fmcba.f:488-521; pn,op/fmcba.f:462-495): WCWMC/WCWMD or
        # PNWMC/PNWMD (the FMR6SDCY tables) by ITYPE.
        if (s.variant isa WestCascades || s.variant isa PacificNorthwest || s.variant isa Olympic) &&
           size(fs.params.dkr, 1) != 11
            _it = Int(s.plot.habitat_input)
            if s.variant isa WestCascades
                _it = clamp(_it, 1, 139); _t = _R6SD_WCWMC[_it]; _m = _R6SD_WCWMD[_it]
            else
                _it = clamp(_it, 1, 75); _t = _R6SD_PNWMC[_it]; _m = _R6SD_PNWMD[_it]
            end
            fs.params.dkr = r6_adjusted_dkr(_FM_DKR_WC, _t, _m;
                                            adj = s.variant isa WestCascades ? _FM_DKRADJ_WC : _FM_DKRADJ_PN)
        end
        # SO decay-rate habitat adjustment (so/fmcba.f:764-836): scale the SO base DKR by DKRADJ(TEMP,MOIST,K)
        # from SOHMC/SOWMD at the first FFE year. The reference stand rides so/habtyp.f's DEFAULT plant
        # association CPS111 = ITYPE 49 (SI 70; same default the SO growth port's SITEAR/SDIDEF ride), giving
        # TEMP=hot/MOIST=dry — explicit-habitat real-FIA stands need the full so/ecocls.f PA decode (documented
        # follow-on, same deferral as growth). habitat_input>0 (a decoded PA) is honored when present.
        if s.variant isa SouthCentralOregon && size(fs.params.dkr, 1) != 11
            _ity = Int(s.plot.habitat_input); _ity <= 0 && (_ity = 49)
            _kf = Int(s.plot.user_forest_code)
            fs.params.dkr = ((500 <= _kf < 600) || _kf == 701) ? so_california_dkr() : so_adjusted_dkr(_ity)
        end
        fs.fuels_init = true
    end
    return s
end

# fmcba.f bare-stand seral cover: COVINI(ITYPE) when 1 ≤ ITYPE ≤ MXVCODE and nonzero, else the variant default.
@inline _covini(tab::AbstractVector, itype::Int, dflt::Integer)::Int32 =
    ((1 <= itype <= length(tab)) && tab[itype] != 0) ? Int32(tab[itype]) : Int32(dflt)

# bm/fmcba.f DATA COVINI(1:92): the seral cover type (species index) per BM plant association (PCOML order).
const _BM_COVINI = Int32[
    9, 9, 3, 3, 3, 3, 3, 3, 3, 3,  3, 3, 3, 9, 9, 9, 9, 8, 8, 8,  8, 8, 9, 9, 9, 9, 9, 9, 9, 9,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7,  7, 7, 7, 7, 7, 7, 5, 5,10,10, 10,10,10,10,10,10,10,10,10,10,
   10,10,10,10,10,10,10, 4, 4, 4,  4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  4, 4, 4, 4, 4, 4, 4, 4, 4, 7,
    7, 7]
