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

"""
    fmcba!(s) -> StandState

Update the stand's `FireState` cover type, percent cover, big DBH, live fuels, and
(first FFE year) dead-fuel pools (FMCBA, fmcba.f). No-op unless FFE is active.
"""
function fmcba!(s::StandState; load_dead::Bool = true)
    fs = s.fire
    (fs === nothing || !fs.active) && return s
    t = s.trees; coef = s.coef
    nsp = length(coef_col(coef, :dbh_min))            # MAXSP

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
           s.variant isa Klamath
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
    _bm_fm = s.variant isa BlueMountains        # BM CRWDTH via bm_cwcalc (BMMAP->cr_cwcalc western library)
    _nc_fm = s.variant isa Klamath              # NC CRWDTH via nc_cwcalc (NCMAP western Bechtold/Crookston library)
    _ws_fm = s.variant isa WestSierra           # WS CRWDTH via ws_cwcalc (WSMAP; R5 forest 511 ⇒ BF=1, same forms as NC)
    _west_cw = _cr_fm || _bm_fm || _nc_fm || _ws_fm
    _cr_ba = _west_cw ? s.plot.basal_area : 0f0
    # NC CRWDTH (base cwidth.f→cwcalc.f) is computed by CWIDTH at LOAD time, BEFORE the stand BA is
    # accumulated ⇒ the R6-Crookston BAREA term hits cwcalc.f:859 `IF(BAREA.LE.1.) BAREA=1.` (BA=0→1).
    # fmcba in the FIRST FFE cycle reads those load-time CRWDTH; later cycles read the end-of-cycle UPDATE
    # (actual BA). MEASURED vs FVSnc_clean DEBUG FMCBA cyc1: PERCOV=39.01 matches BAREA=1 (jl 39.72), NOT the
    # dense stand BA (jl 44.29). Mirror the load-time clamp for cycle 1.
    _nc_ba = ((_nc_fm || _ws_fm) && s.control.cycle <= Int32(1)) ? 1f0 : _cr_ba
    _cr_el = _west_cw ? s.plot.elevation : 0f0
    _cr_hi = _west_cw ? _cr_hopkins(s.plot.latitude, s.plot.longitude, s.plot.elevation) : 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] > 0f0 || continue
        sp = Int(t.species[i]); d = t.dbh[i]
        tba[sp] += 3.14159f0 * (d / 24f0) * (d / 24f0) * t.tpa[i]
        d > fs.bigdbh && (fs.bigdbh = d)
        cw = _cr_fm ? cr_cwcalc(sp, d, t.height[i], Float32(t.crown_pct[i]), _cr_ba, _cr_el, _cr_hi) :
             _bm_fm ? bm_cwcalc(sp, d, t.height[i], Float32(t.crown_pct[i]), _cr_ba, _cr_el, _cr_hi) :
             _nc_fm ? nc_cwcalc(sp, d, t.height[i], Float32(t.crown_pct[i]), _nc_ba, _cr_el, _cr_hi) :
             _ws_fm ? ws_r5crwd(sp, d, t.height[i]) :   # WS: R5CRWD (ws/r5crwd.f), function of sp/D/H only
             crown_width(coef, s.species.code2[sp], d, t.height[i], Float32(t.crown_pct[i]), 0,
                         s.plot.latitude, s.plot.longitude, s.plot.elevation)   # forest-grown (CWCALC iwho=0)
        totcra += 3.1415927f0 * cw * cw / 4f0 * t.tpa[i]
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
    if s.variant isa Klamath || s.variant isa WestSierra
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
                 s.variant isa CentralIdaho ? Int32(3)  :     # CI bare-stand default: DF (cit01 has trees ⇒ unused)
                 s.variant isa Teton ? Int32(3)  :             # TT bare-stand default: DF
                 s.variant isa Utah ? Int32(3)  :              # UT bare-stand default: DF
                 s.variant isa BlueMountains ? Int32(3) :      # BM bare-stand default: DF (bm/fmcba.f COVINI/3)
                 s.variant isa Northeast     ? Int32(1)  :
                 s.variant isa CentralStates ? Int32(48) :
                 s.variant isa LakeStates    ? Int32(3)  :
                 # NC bare stand: COVINI5(ITYPE)/COVINI6(ITYPE) by habitat (nc/fmcba.f:337-355); the full
                 # R5/R6 habitat→cover maps are a bare-stand-only path not exercised by nct01 (trees present) —
                 # the fmcba.f "no valid habitat" fallback is Douglas-fir (3), used here until those maps port.
                 s.variant isa Klamath ? Int32(3) :
                 s.variant isa CentralRockies ? Int32(11) : Int32(75)   # CR: lodgepole pine (fmcba.f:432)
    end
    fs.covtyp = covtyp
    fs.percov = (1f0 - exp(-totcra / 43560f0)) * 100f0
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
                  ffe_dead_fuel_loading(coef, Int(s.plot.forest_type))
        # Seed the STFUEL override from FIA-DB measured fuel loadings (FVS_STANDINIT FUEL_* → dbsstandin.f
        # FUELINIT, read into plot.ffe_fuel_*) when present AND no explicit FUELINIT/FUELSOFT keyword already set
        # them (the keyword takes precedence, matching FVS where a later-scheduled FUELINIT overrides the DB one).
        isempty(fs.params.stfuel_hard) && !isempty(s.plot.ffe_fuel_hard) && (fs.params.stfuel_hard = copy(s.plot.ffe_fuel_hard))
        isempty(fs.params.stfuel_soft) && !isempty(s.plot.ffe_fuel_soft) && (fs.params.stfuel_soft = copy(s.plot.ffe_fuel_soft))
        ovh = fs.params.stfuel_hard; ovs = fs.params.stfuel_soft
        fill!(fs.cwd, 0f0)
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
            fs.params.dkr = bm_adjusted_dkr(Int(s.plot.habitat_code))
        end
        # NC decay-rate DCYMLT (nc/fmcba.f:395-414): scale the NC base DKR by the Dunning-code/site-index
        # multiplier at the first FFE year (when the user hasn't set FuelDcay ⇒ params.dkr still empty).
        if s.variant isa Klamath && size(fs.params.dkr, 1) != 11
            _ss = Int(s.plot.site_species)
            _si = (1 <= _ss <= length(s.plot.sp_site_index)) ? s.plot.sp_site_index[_ss] : 0f0
            fs.params.dkr = nc_adjusted_dkr(_si)
        end
        fs.fuels_init = true
    end
    return s
end
