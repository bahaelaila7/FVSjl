# =============================================================================
# diameter_growth.jl — Southern diameter-growth model (DGF + DGCONS)
#
# Ported from: sn/dgf.f.
#
# `dgcons!` builds the per-species site constant DGCON from site index, slope,
# aspect and the ecological-unit (PCOM) physiographic terms. `dgf!` then evaluates
# ln(DDS) (change in squared inside-bark diameter) for every tree into the scratch
# column `wk2`. The stochastic driver (DGDRIV: calibration + serial correlation)
# turns ln(DDS) into per-tree DG and is ported separately.
#
# Coefficient tables live in dgf_coefficients.jl. The species loop is flattened to
# a single per-tree pass (the Fortran species-chain ordering doesn't affect the
# per-tree value).
# =============================================================================

# Forest-type group of IFORTP → which categorical coefficient applies (dgf.f:453).
function _dgf_forest_group(ifortp::Integer)
    ifortp == 701 || ifortp == 801 || ifortp == 805 ? :nohd :
    ifortp in (104,105,121,124) ? :sfhp :
    ifortp in (141,142,161,162,163,164,166) ? :ylpn :
    ifortp in (165,403,404,405,406,407,409) ? :okpn :
    ifortp in (501,502,503,504,505,510,514,515) ? :upok :
    ifortp in (168,508,601,602,605,606,607,608,702,703,704,705,706,708) ? :lohd :
    ifortp in (103,167,181,401,402,506,507,511,512,513,519,520,802,803,807,809) ? :uphd :
    :none
end

# Parse the SN ecological-unit code (PCOM) → physiographic flag symbol (dgf.f:567).
function _dgf_phys_group(pcom::AbstractString)
    s   = rpad(pcom, 5)
    ch1 = s[1]; c234 = s[1:3]; c45 = s[4:5]; c4 = s[4]; c5 = s[5]
    if ch1 == 'M'
        s[2:4] == "221" && return :pm221
        s[2:4] in ("222","223") && return :pm222
        s[2:4] == "231" && return :pm231
        return :none
    end
    if c234 == "221"
        c45 in ("DD","DE") && return :s231t
        c45 in ("EJ","EG","EN") && return :p222
        return :p221
    elseif c234 in ("222","223")
        return :p222
    elseif c234 == "231"
        c4 in ('A','C','D','I') && return :s231t
        c4 in ('B','E','F','G') && return :s231l
        if c4 == 'H'
            c5 in ('A','B','C') && return :s231l
            c5 in ('D','E','F') && return :p222
            c5 in ('H','I') && return :p232
        end
        return :none
    elseif c234 == "232"
        return c45 == "FF" ? :p234 : :p232
    elseif c234 == "234"
        return :p234
    elseif c234 == "251"
        return :p255
    elseif c234 == "255"
        return c45 == "CH" ? :s231l : :p255
    elseif c234 == "411"
        return :p411
    end
    return :none
end

"""
    dgcons!(state)

Load the site-dependent DG constant `calib.dg_const[sp]` and `calib.atten[sp]`
(DGCONS, dgf.f:564). Call once per stand after site variables are known.
"""
# physiographic group symbol → its CSV coefficient-column name.
const _DG_PHYS_COL = Dict(
    :pm221 => :dg_phys_pm221, :pm222 => :dg_phys_pm222, :pm231 => :dg_phys_pm231,
    :p221  => :dg_phys_p221,  :p222  => :dg_phys_p222,  :s231t => :dg_phys_s231t,
    :s231l => :dg_phys_s231l, :p232  => :dg_phys_p232,  :p234  => :dg_phys_p234,
    :p255  => :dg_phys_p255,  :p411  => :dg_phys_p411)

# forest-type group symbol → its CSV coefficient-column name.
const _DG_FORTYPE_COL = Dict(
    :lohd => :dg_fortype_lowland_hw, :nohd => :dg_fortype_nonoak_hw,
    :okpn => :dg_fortype_oak_pine,   :sfhp => :dg_fortype_softwd_hp,
    :uphd => :dg_fortype_upland_hw,  :upok => :dg_fortype_upland_oak,
    :ylpn => :dg_fortype_yellow_pine)

function dgcons!(s::StandState)
    p, c, sd = s.plot, s.calib, s.coef.species
    prior_obs_count = sd[:dg_prior_obs_count]
    site_coef = sd[:dg_site_index]
    slope_tan = sd[:dg_slope_tan]; slope_cos = sd[:dg_slope_cos_aspect]; slope_sin = sd[:dg_slope_sin_aspect]
    phys = _dgf_phys_group(s.plot.eco_unit)
    phys_coef = haskey(_DG_PHYS_COL, phys) ? sd[_DG_PHYS_COL[phys]] : nothing
    cosa = cos(p.aspect); sina = sin(p.aspect)
    bark_a0 = sd[:bark_intercept]; bark_b0 = sd[:bark_slope]
    @inbounds for sp in 1:MAXSP
        c.atten[sp] = prior_obs_count[sp]
        base = site_coef[sp] * p.sp_site_index[sp] +
               slope_tan[sp] * p.slope +
               slope_cos[sp] * p.slope * cosa +
               slope_sin[sp] * p.slope * sina
        phys_coef === nothing || (base += phys_coef[sp])
        c.dg_const[sp] = base
        c.bark_a[sp] = bark_a0[sp]; c.bark_b[sp] = bark_b0[sp]   # per-stand bark copy
    end
    # Fort Bragg (IFOR==20, dgf.f:636 / bratio.f:106): override the DG attenuation
    # count for longleaf/loblolly and the special bark equations (sp 5,6,8,11,13).
    if p.forest_idx == 20
        c.atten[8] = 2056f0; c.atten[13] = 689f0
        c.bark_a[5]  = 0.1713f0;  c.bark_b[5]  = 0.87459f0
        c.bark_a[6]  = -0.26207f0; c.bark_b[6] = 0.87347f0
        c.bark_a[8]  = -0.43439f0; c.bark_b[8] = 0.91382f0
        c.bark_a[11] = -0.62033f0; c.bark_b[11]= 0.91645f0
        c.bark_a[13] = -0.4671f0;  c.bark_b[13]= 0.90198f0
    end
    # READCORD/READCORH (dgf.f:1168 / htgf.f:332): add ln(COR2)/ln(HCOR2) to DGCON/HTCON before
    # the LSTART calibration. Default off (terms = 1 ⇒ ln 1 = 0, no-op); guarded > 0 like Fortran.
    ctl = s.control
    if ctl.dg_cor2_on
        @inbounds for sp in 1:MAXSP
            ctl.dg_cor2[sp] > 0f0 && (c.dg_const[sp] += log(ctl.dg_cor2[sp]))
        end
    end
    if ctl.htg_cor2_on
        @inbounds for sp in 1:MAXSP
            ctl.htg_cor2[sp] > 0f0 && (c.htg_cor[sp] += log(ctl.htg_cor2[sp]))
        end
    end
    return s
end

"""
    dgf!(state)

Evaluate ln(DDS) for every tree into `scratch.wk[2, i]` (DGF, dgf.f:436), using
the current diameters. `cor` is the per-species calibration correction (0 until
calibration runs). Requires `dgcons!` and AVH/BA already set.
"""
function dgf!(s::StandState, ::Southern)
    p, t, c, sd = s.plot, s.trees, s.calib, s.coef.species
    wk2 = view(s.scratch.wk, 2, :)
    intercept = sd[:dg_intercept]; ln_dbh = sd[:dg_ln_dbh]; dbh_sq = sd[:dg_dbh_squared]
    ln_crown = sd[:dg_ln_crown_pct]; rel_ht = sd[:dg_relative_height]
    stand_ba_c = sd[:dg_stand_basal_area]; point_bal = sd[:dg_point_bal]; planted = sd[:dg_planted]
    ftgrp = _dgf_forest_group(p.forest_type)
    ft_coef = haskey(_DG_FORTYPE_COL, ftgrp) ? sd[_DG_FORTYPE_COL[ftgrp]] : nothing
    kplant = p.managed > 0 ? 1f0 : 0f0
    avh = p.avg_height
    @inbounds for i in 1:t.n
        d = t.dbh[i]
        d <= 0f0 && continue
        sp = t.species[i]
        conspp = c.dg_const[sp] + c.dg_cor[sp]

        relht = 0f0
        if avh > 0f0
            relht = t.height[i] / avh
            relht > 1.5f0 && (relht = 1.5f0)
        end
        ba_v  = p.basal_area <= 0f0 ? 25f0 : p.basal_area
        icr_i = t.crown_pct[i] <= 0 ? 25 : t.crown_pct[i]
        # DGF competition term: pbal = (point BA)·(1 − PCT/100), where PCT is the
        # CROWN-modeled crown ratio (crown_ratio, set by crown_ratio!). pba=PTBAA.
        bal   = (1f0 - t.crown_ratio[i] / 100f0) * ba_v
        pba   = s.density.point_ba[t.plot_id[i]]
        pba <= 0f0 && (pba = ba_v)
        pbal  = pba * (1f0 - t.crown_ratio[i] / 100f0)
        pbal <= 0f0 && (pbal = bal)

        # NB dgf.f uses ALOG here; PROVEN inert (2026-07-05fff: routing log→flog left snt01 bit-exact + suite
        # unchanged ⇒ openlibm log == gfortran ALOG for the real dbh/icr ranges). Kept as Julia `log` — this is
        # the PER-TREE hot loop, so the fpow-companion ccall's ~30% suite cost isn't justified for a zero-diff op
        # (doctrine #8 caveat: only wire the FFI for ops that ACTUALLY differ). The grown-cycle drift is NOT here.
        dds = conspp + intercept[sp] +
              ln_dbh[sp]     * log(d) +
              dbh_sq[sp]     * d * d +
              ln_crown[sp]   * log(Float32(icr_i)) +
              rel_ht[sp]     * relht +
              stand_ba_c[sp] * ba_v +
              point_bal[sp]  * pbal +
              (ft_coef === nothing ? 0f0 : ft_coef[sp]) +
              planted[sp]    * kplant

        # Fort Bragg (IFOR==20, dgf.f:515-537): special longleaf(8)/loblolly(13)
        # diameter-growth equations replace the standard DDS. dg5 is the inside-bark
        # DG; dds = ln(dg5·(2·d_ib + dg5)) re-encodes it for the sqrt growth formula.
        if p.forest_idx == 20 && (sp == 8 || sp == 13)
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
            dib  = d * bark
            site = p.sp_site_index[sp]
            pctf = Float32(t.crown_ratio[i])              # PCT (BA percentile)
            cr   = Float32(icr_i) / 100f0
            dg5 = sp == 8 ?
                dib * (-0.4553f0 * (0.09737f0 - exp(-0.2428f0 * d)) + 0.05574f0 * cr -
                       0.0002965f0 * ba_v - 0.00002481f0 * pba -
                       0.001192f0 * (pctf / 100f0)^(-0.9663f0) +
                       0.0010110f0 * site - 0.007711f0 * relht) :
                dib * (-0.3428f0 * (-0.1741f0 - exp(-0.1328f0 * d)) + 0.1145f0 * cr -
                       0.0001682f0 * ba_v - 0.00003978f0 * pba -
                       0.159400f0 * (pctf / 100f0)^(-0.1299f0) +
                       0.0006204f0 * site + 0.02474f0 * relht)
            dg5 < 0.01f0 && (dg5 = 0.01f0)
            dds = log(dg5 * (2f0 * dib + dg5))
        end

        dds < -9.21f0 && (dds = -9.21f0)
        wk2[i] = dds
    end
    return s
end

"Empirical-Bayes prior variance for DG calibration (dgdriv.f PSIGSQ)."
const DG_PSIGSQ = 0.089827273f0

"""
    calibrate_diameter_growth!(state; scale=1f0)

DGDRIV calibration pass (LSTART, dgdriv.f:150). For each species with enough
measured diameter growth, regress the DGF residuals to get the large-tree
calibration `calib.dg_cor` (COR, empirical-Bayes shrunk) and the attenuation goal
`calib.dg_cor_goal` (WCI). Also seeds the per-tree serial-correlation residual
`trees.old_random` (OLDRN) — measured trees get their residual, calibrated species
fill the rest by regression, and uncalibrated species draw from BACHLO — and the
per-species `calib.vardg` (VARDG). Consumes the RNG, so trees are walked in FVS's
species-sorted order. `scale = YR/FINT` (1 for snt01). Run before `diameter_growth!`.
"""
# GROWTH IDG/IHTG=1/3 (sn/cratet.f:169-180): the input DG/HTG field is a PAST (code 1) or CURRENT
# (code 3) DBH/HT measurement, not the growth increment. Convert it to the increment
# `DG = Q·(DBH − field)` with Q=+1 for code 1 (field = past, current in DBH) and Q=−1 for code 3
# (field = current, past in DBH). A field ≤0 is a MISSING measurement → −1 sentinel (every measured-
# DG consumer filters `≤0`, so the tree falls back to the stand-average backdating); the sentinel is
# essential — without it field=0 would yield the bogus increment DBH−0 = DBH. The diameter
# increment here is OUTSIDE-bark (matching cratet feeding DENSE); the outside→inside BRATIO bark
# correction is applied later, in `calibrate_diameter_growth!` (dgdriv.f:330-333). Height has no
# bark, so the HTG increment is final. Gated: a no-op for the default IDG/IHTG=0.
function apply_growth_input_types!(s::StandState)
    t = s.trees
    idg = s.control.growth_idg
    if idg == 1 || idg == 3
        q = idg == 1 ? 1f0 : -1f0
        @inbounds for i in 1:t.n
            f = t.diam_growth[i]
            t.diam_growth[i] = f <= 0f0 ? -1f0 : q * (t.dbh[i] - f)
        end
    end
    ihtg = s.control.growth_ihtg
    if ihtg == 1 || ihtg == 3
        q = ihtg == 1 ? 1f0 : -1f0
        @inbounds for i in 1:t.n
            f = t.ht_growth[i]
            t.ht_growth[i] = f <= 0f0 ? -1f0 : q * (t.height[i] - f)
        end
    end
    return
end

# Backdate live diameters to the start of the measured-growth period, IN PLACE
# (DENSE/LBKDEN, dense.f:70-128). WK3 = sqrt(d²·r): for a measured-DG tree,
# r = 1−(2·d·gadj − gadj²)/d² (the past inside-bark dbh); unmeasured trees fall back to
# the stand-average BA-growth ratio BAGR. IDG-faithful (dense.f:100-127): exclude only
# MISSING growth (IDG=1/3 missing is already −1 → a genuine 0 is KEPT; IDG=0/2 keep
# input-0 as missing); bark-divide the increment for IDG∈{0,2,3}, SKIP for IDG==1
# (dense.f:101,122). Shared by calibrate_diameter_growth! and init_crown_ratios! so the
# crown's backdated CCF and the DG calibration read ONE backdating — as FVS runs one DENSE.
function _backdate_dbh!(s::StandState)
    t = s.trees; n = t.n
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    idg = s.control.growth_idg
    _cr_bd = s.variant isa CentralRockies; _cr_bd_imod = _cr_bd ? Int(s.plot.model_type) : 0; sd = s.coef.species
    _tt_bd = s.variant isa Teton
    _bm_bd = s.variant isa BlueMountains   # BM bark = bm_bratio (POWER model DIB=BARK1·D^BARK2, not linear a+b·d)
    _bc_bd = s.variant isa BritishColumbia   # BC bark = bc_bratio (constant; shared bark_a/bark_b=0 ⇒ wrong 0.80 floor)
    _ci_bd = s.variant isa CentralIdaho      # CI bark = ci_bratio (POWER DIB=BARK1·D^BARK2) — MISSING branch left the
                                             # DENSE backdating on the linear default (0.9) ⇒ backdated BA 0.13% high
                                             # ⇒ COR fit vs wrong density ⇒ DG low ⇒ ~2% mortality over-kill.
    _ak_bd = s.variant isa SoutheastAlaska   # AK bark = ak_bratio (3-type: power/linear/power); shared bark_a/bark_b unset
    _wc_bd = s.variant isa WestCascades      # WC bark = wc_bratio (POWER a·Dᵇ, all imap=1) — shared linear floors 0.80,
                                             # corrupting the backdated dbh/density ⇒ wrong WF calibration COR (Δ0.098).
    _pn_bd = s.variant isa PacificNorthwest  # PN bark = wc_bratio (POWER, all imap=1) — same #140/CI class as WC
    _ec_bd = s.variant isa EastCascades      # EC bark = wc_bratio (per-species bark_imap POWER/linear) — same watchpoint
    _ca_bd = s.variant isa CentralCalifornia # CA bark = wc_bratio (per-species bark_imap) — same #140/EC class as WC/EC
    _so_bd = s.variant isa SouthCentralOregon # SO bark = so_bratio (so/bratio.f 3-path: CASE1 BARKB / juniper / BRDAT)
    _op_bd = s.variant isa Olympic            # OP bark = op_bratio (op/bratio.f 3-path) — same watchpoint class as WC (WF COR)
    _bk(sp, d) = _cr_bd ? cr_bratio(sd, Int(sp), d, _cr_bd_imod) :
                 _tt_bd ? tt_bratio(Int(sp), Float32(d)) :
                 _bm_bd ? bm_bratio(sd, Int(sp), Float32(d)) :
                 _ci_bd ? ci_bratio(sd, Int(sp), d) :
                 _ak_bd ? ak_bratio(Int(sp), Float32(d)) :
                 _wc_bd ? wc_bratio(sd, Int(sp), Float32(d)) :
                 _pn_bd ? wc_bratio(sd, Int(sp), Float32(d)) :
                 (_ec_bd || _ca_bd) ? wc_bratio(sd, Int(sp), Float32(d)) :
                 _so_bd ? so_bratio(sd, Int(sp), Float32(d)) :
                 _op_bd ? op_bratio(Int(sp), Float32(d)) :        # OP bark = op_bratio (op/bratio.f 3-path) — backdating must match dgf/update
                 _bc_bd ? bc_bratio(Int(sp)) : bark_ratio(bark_a, bark_b, sp, d)
    ismiss = (idg == 1 || idg == 3) ? (g -> g < 0f0) : (g -> g <= 0f0)
    bagr = 0f0; nb = 0f0
    @inbounds for i in 1:n
        g = t.diam_growth[i]; ismiss(g) && continue
        d = t.dbh[i]
        gadj = idg == 1 ? g : g / _bk(t.species[i], d)
        gadj > d && continue
        bagr += 1f0 - (2f0 * d * gadj - gadj * gadj) / (d * d); nb += 1f0
    end
    nb > 0f0 && (bagr /= nb)
    @inbounds for i in 1:n
        d = t.dbh[i]; g = t.diam_growth[i]; r = bagr
        if !ismiss(g)
            gadj = idg == 1 ? g : min(g / _bk(t.species[i], d), d)
            rr = 1f0 - (2f0 * d * gadj - gadj * gadj) / (d * d)
            rr > 0f0 && (r = rr)
        end
        r > 0f0 && (t.dbh[i] = sqrt(d * d * r))
    end
    return s
end

function calibrate_diameter_growth!(s::StandState; scale::Float32 = 1f0, fnmin::Float32 = 5f0)
    t, c = s.trees, s.calib
    sd = s.coef.species
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    sigmar = s.variant isa Olympic ? OP_DG_SIGMAR : sd[:dg_resid_sd]   # OP SIGMAR (op/blkdat.f); OP CSV has no dg_resid_sd column
    _op_cal = s.variant isa Olympic           # OP bark = op_bratio (op/bratio.f) — same watchpoint class as WC (WF COR)
    _cr_cal = s.variant isa CentralRockies; _cr_cal_imod = _cr_cal ? Int(s.plot.model_type) : 0
    _tt_cal = s.variant isa Teton   # TT bark = tt_bratio (PP sp10 IMAP=4 power model, not linear a+b·d)
    _bm_cal = s.variant isa BlueMountains   # BM bark = bm_bratio (POWER model)
    _bc_cal = s.variant isa BritishColumbia   # BC bark = bc_bratio (constant BARK1; shared bark_a/bark_b are 0 → 0.80 floor)
    _ci_cal = s.variant isa CentralIdaho      # CI bark = ci_bratio (POWER) — same missing-branch class as _backdate_dbh!
    _ak_cal = s.variant isa SoutheastAlaska   # AK bark = ak_bratio (3-type); shared bark_a/bark_b unset ⇒ 0.80 floor otherwise
    _wc_cal = s.variant isa WestCascades      # WC bark = wc_bratio (POWER, all imap=1) — shared linear floors 0.80 ⇒ wrong COR
    _pn_cal = s.variant isa PacificNorthwest  # PN bark = wc_bratio (POWER) — apply the WC watchpoint from the start
    _ec_cal = s.variant isa EastCascades      # EC bark = wc_bratio (per-species bark_imap) — WC watchpoint from the start
    _ca_cal = s.variant isa CentralCalifornia # CA bark = wc_bratio (per-species bark_imap) — same #140/EC class (WF COR fix)
    _so_cal = s.variant isa SouthCentralOregon # SO bark = so_bratio (so/bratio.f 3-path)
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    species_sort!(s)

    # Backdate diameters to the start of the measured-growth period (DENSE/LBKDEN,
    # dense.f:70-86): WK3 = sqrt(d²·r). For a tree with measured DG, r=(d−DG/bark)²/d²
    # (so WK3 = past inside-bark-adjusted dbh); unmeasured trees use the stand-average
    # ratio bagr. The calibration DGF predicts from this PAST stand state (past dbh +
    # past BA/PCT/point_ba), EXCEPT the AVHT40 top height (AVH), which stays at the
    # CURRENT stand value — see the dgf! call below.
    saved_dbh = Float32[t.dbh[i] for i in 1:t.n]
    _cr_bd_ccf = 0f0                           # CR: BACKDATED stand CCF (dense.f RELDM1) for the REGENT height calib PCTRED
    _cr_bd_avht = 0f0                          # CR: BACKDATED-window AVHT40 (dense.f AVH) for the same PCTRED (X=AVH·RELDEN/100)
    _bc_bd_ba = 0f0; _bc_bd_relden = 0f0; _bc_bd_pct = Float32[]   # BC: BACKDATED BA/RELDEN/percentile for V2 small-tree HCOR calib (regent.f REGCAL uses the backdated stand)
    _cur_avh = s.plot.avg_height   # current-stand AVHT40 top height (used by the calibration DGF below)
    # NOTRE inflates DEAD-record PROB by FINT/FINTM (cycle-growth period / mortality-observation period) so the
    # recent dead are added back at the right rate to recover the BACKDATED density (notre.f:122-124). FVS keeps
    # that inflation live through calibration and "undoes" it for treelist/FFE uses (FMSSEE/PRTRLS); jl carries
    # the TRUE dead TPA everywhere, so it applies the inflation ONLY here, scoped to the two calibration density
    # passes below, and restores it before returning. Inert when FINTM=FINT (the default).
    _fintr = s.control.growth_fintm > 0f0 ? s.control.growth_fint / s.control.growth_fintm : 1f0
    _livec = t.n
    _saved_dead_tpa = (_fintr != 1f0 && t.ndead > 0) ?
        Float32[t.tpa[j] for j in (_livec + 1):(_livec + t.ndead)] : Float32[]
    if _fintr != 1f0
        @inbounds for j in (_livec + 1):(_livec + t.ndead); t.tpa[j] *= _fintr; end
    end
    # PTBAA (point basal area) in the calibration prediction uses the CURRENT-DBH
    # point BA (dgf.f:493 reads the live PTBAA the last DENSE pass filled at current
    # diameters). The percentile population is the full ITRN, so recently-dead trees
    # (history≠8) count toward PTBAA at their current dbh too (long-dead history-8 are
    # zeroed). Compute it here with the dead partition exposed at current dbh, then
    # restore. (Stand BA/AVH/PCT stay backdated; only this point-BA total is current.)
    # PTBAA uses CURRENT dbh for ALL dead trees (history 6 and 8) — a prism tree
    # contributes a fixed point-BA regardless of dbh, so even a long-dead snag counts.
    # (Only the BACKDATED percentile zeroes history-8; PTBAA does not.)
    let nlive0 = t.n
        t.n = nlive0 + t.ndead
        point_basal_area!(s)
        t.n = nlive0
    end
    cur_point_ba = copy(s.density.point_ba)
    # PTBALT (point BA-in-larger-trees) is CURRENT-dbh in calibration, exactly like PTBAA above.
    # dense.f runs a two-pass backdating: pass-1 (LREDO=T) sets WK5=backdated_dbh²·PROB and feeds
    # PCTILE (the backdated percentile); pass-2 (LREDO=F, dense.f:184) reverts D=DBH(I) so WK5=
    # CURRENT_dbh²·PROB, and PTBAL (dense.f:280, ptbal.f:148 XBALT+=WK5·.005454154·PI/GROSPC) reads
    # that CURRENT WK5. So the calibration DGF's PBAL=PTBALT(I) is current, NOT backdated. jl's
    # compute_density! below overwrites point_bal with the backdated value; capture the current one
    # here (from point_basal_area! above, dead-exposed, same as cur_point_ba) and restore it after.
    # AK-scoped: this is shared-faithful (ptbal.f is common to all western/southern variants) and the
    # same latent gap likely applies cluster-wide, but is validated here only for AK (akt01). On akt01
    # YC (sp3) the omission gave PBAL 52.8 vs live 187.5 (=3×BAF 62.5) ⇒ COR 1.389 vs live 1.434.
    _ak_pbal_fix = s.variant isa SoutheastAlaska
    cur_point_bal = _ak_pbal_fix ? copy(s.density.point_bal) : Float32[]
    # #191: stash the CURRENT-stand RMSQD before backdating so the TT aspen DGFASP calibration prediction uses it
    # (FVS uses current RMSQD in the calibration DGFASP, like the AVH exception below; jl's stand_qmd on the
    # backdated stand would under-predict aspen ⇒ measured>>predicted ⇒ COR falsely BOOSTS aspen DG).
    _TT_CUR_RMSQD[] = stand_qmd(s)    # #195: current RMSQD for the aspen DGFASP calibration (ALL variants: TT/UT/BM/CI/EM/IE aspen dgf! read it; others ignore)
    _backdate_dbh!(s)                         # dense.f:70-128 backdating (IDG-faithful); shared w/ init_crown_ratios!
    # The backdated stand BA/AVH still include the dead trees (kept at current dbh):
    # expose the dead partition for this density pass, then restore. (PTBAA itself is
    # overridden below with the live-only current point_ba; only BA/AVH/PCT use this.)
    nlive = t.n
    saved_dead = Float32[t.dbh[j] for j in (nlive + 1):(nlive + t.ndead)]
    @inbounds for j in (nlive + 1):(nlive + t.ndead)
        t.history[j] == 8 && (t.dbh[j] = 0f0)
    end
    t.n = nlive + t.ndead
    compute_density!(s)                       # past-stand BA/AVH/point_ba/PCT
    # CR REGENT height calib uses the BACKDATED CCF (dense.f RELDM1) for its density modifier PCTRED — NOT the
    # current CCF (regent.f:466 X=AVH·RELDEN/100; RELDEN is the backdated relative density). Capture it HERE with
    # the recently-dead (backdated) trees INCLUDED (t.n = nlive+ndead): they were alive at the period start so
    # dense.f's RELDM1 counts them (history-8 already zeroed above). Live-only omission under-counted the CCF on
    # stands with recent mortality (68 vs live 121 on 1855925743290487); dead-inclusive = 121.13 = live exact.
    # AVH stays CURRENT (not backdated), as in live. (18th-bug fix + factor-2 dead-inclusion correction.)
    _cr_cal && (_cr_bd_ccf = stand_ccf(s); _cr_bd_avht = stand_top_height(s))
    t.n = nlive
    @inbounds for (k, j) in enumerate((nlive + 1):(nlive + t.ndead))
        t.dbh[j] = saved_dead[k]
    end
    s.density.point_ba .= cur_point_ba        # PTBAA from current DBH, live-only (above)
    _ak_pbal_fix && (s.density.point_bal .= cur_point_bal)   # PTBALT current-dbh too (dense.f pass-2; see above)
    # PCT (BA percentile, dense.f pass-1 PCTILE) accumulates the BACKDATED point-BA
    # weights in the CURRENT-dbh rank order (IND is fixed at setup). The percentile
    # POPULATION is the full ITRN including recently-dead trees (history≠8) at their
    # CURRENT dbh — long-dead (history 8) are in the order but zeroed. Only the live
    # trees' crown_ratio is read downstream by dgf!. compute_density! above re-sorted
    # by backdated dbh; recompute crown_ratio here.
    let nlive2 = t.n, ntot = t.n + t.ndead
        rankd = Float32[i <= nlive2 ? saved_dbh[i] : t.dbh[i] for i in 1:ntot]   # current dbh
        wk5   = Float32[i <= nlive2 ? t.dbh[i]^2 * t.tpa[i] :                     # live: backdated
                        (t.history[i] == 8 ? 0f0 : t.dbh[i]^2 * t.tpa[i]) for i in 1:ntot]  # dead: current/0
        # Rank order: FVS's IND from `RDPSRT(ITRN,DBH,IND,.TRUE.)` (gradd.f:186 / dense.f PCTILE) is
        # Scowen's UNSTABLE quicksort — on a current-dbh TIE between a live tree and a same-dbh recently-
        # dead tree, it can order the dead one first, dropping the live tree below 100th percentile. KT
        # needs this exact tie-break (validated: a dead 16.1" ties live tree1, live PCT=89.068 not 100).
        # SN/CR are validated with the stable `sortperm`; keep it for them.
        if s.variant isa Kootenai
            ord = Vector{Int32}(undef, ntot)
            _rdpsrt!(rankd, ord)
        else
            ord = sortperm(rankd; rev = true)
        end
        tot = sum(wk5); cum = 0f0
        if tot > 0f0
            @inbounds for k in ntot:-1:1
                ii = Int(ord[k])
                cum += wk5[ii]
                ii <= nlive2 && (t.crown_ratio[ii] = cum / tot * 100f0)
            end
        end
    end
    # BC V2 small-tree HCOR calibration reads the BACKDATED stand (regent.f REGCAL): capture BA/RELDEN/PCT here.
    _bc_cal && (_bc_bd_ba = s.plot.basal_area; _bc_bd_relden = s.plot.relative_density;
                _bc_bd_pct = Float32[t.crown_ratio[i] for i in 1:t.n])
    if _fintr != 1f0                          # restore TRUE dead TPA (the FINT/FINTM inflation was calibration-
        @inbounds for (k, j) in enumerate((_livec + 1):(_livec + t.ndead)); t.tpa[j] = _saved_dead_tpa[k]; end
    end                                       # only — covers both density passes AND the PCTILE percentile)
    # GROWTH IDG=1/3 bark correction (dgdriv.f:330-333): the measured increment is now the
    # OUTSIDE-bark DBH difference (DBH−past), which DENSE used to backdast above. Convert it to
    # the INSIDE-bark increment with BRATIO at the CURRENT dbh (`saved_dbh`) before the calibration
    # term, so it matches the IDG=0 inside-bark basis. IDG=0/2 supply the inside increment directly.
    if s.control.growth_idg == 1 || s.control.growth_idg == 3
        @inbounds for i in 1:t.n
            if t.diam_growth[i] > 0f0
                # CR: BRATIO = cr_bratio (bark_a/bark_b are 0 ⇒ bark_ratio floors to 0.80, but cr/bratio.f
                # gives ~0.95). Using the 0.80 floor shrinks the measured inside-bark DG ⇒ the DGSCOR
                # cornew drifts ~0.19 more negative, crossing the exp(-2.5)=0.0821 COR out-of-range trap
                # (dgdriv.f:640) ⇒ COR falsely zeroed for measured-DG species (e.g. aspen sp20), which then
                # over-grows where gemdg is explosive on small DBH. 4th CR variant-bark location.
                bk = _cr_cal ? cr_bratio(sd, Int(t.species[i]), saved_dbh[i], _cr_cal_imod) :
                     _tt_cal ? tt_bratio(Int(t.species[i]), saved_dbh[i]) :
                     _bm_cal ? bm_bratio(sd, Int(t.species[i]), saved_dbh[i]) :
                     _ci_cal ? ci_bratio(sd, Int(t.species[i]), saved_dbh[i]) :
                     _ak_cal ? ak_bratio(Int(t.species[i]), saved_dbh[i]) :
                     _wc_cal ? wc_bratio(sd, Int(t.species[i]), saved_dbh[i]) :
                     _pn_cal ? wc_bratio(sd, Int(t.species[i]), saved_dbh[i]) :
                     (_ec_cal || _ca_cal) ? wc_bratio(sd, Int(t.species[i]), saved_dbh[i]) :
                     _so_cal ? so_bratio(sd, Int(t.species[i]), saved_dbh[i]) :
                     _op_cal ? op_bratio(Int(t.species[i]), saved_dbh[i]) :   # OP bark = op_bratio (op/bratio.f) — COR must use it
                     _bc_cal ? bc_bratio(Int(t.species[i])) :
                     bark_ratio(bark_a, bark_b, t.species[i], saved_dbh[i])
                t.diam_growth[i] *= bk
            end
        end
    end

    # FORTYP is computed in the GROW path (per cycle), AFTER the LSTART calibration,
    # so the calibration prediction must NOT include the forest-type term (kuphd etc.
    # are all 0 at LSTART). Zero it for this dgf!, restore after. (dgf.f:453 reads
    # IFORTP, which is 0 until the first STKVAL/FORTYP call in the cycle loop.)
    saved_fortype = s.plot.forest_type
    s.plot.forest_type = 0
    # NE: the BAL competition (ne_badist!) uses the CURRENT stand even though the per-tree prediction is at the
    # backdated dbh — FVS NE calib computes BADIST on the current stand (verified EBAU=52). Stash the current dbh
    # for ne_badist! to read; SN ignores it (point_bal-based, never calls ne_badist!). Cleared right after.
    s.variant isa Northeast && (c.calib_dbh = saved_dbh)
    # AK: the DGF's point-Zeide PRD (SDICAL XMAXPT + SDICLS ZRD) reads the UNCHANGED DBH(I) = CURRENT
    # dbh even during calibration (like PTBALT/PTBAA above) — FVS backdates only DIAM(I), not the DBH
    # array SDICAL/SDICLS sum. Stash the current dbh so ak_point_zeide! uses it (else jl computes PRD on
    # the backdated stand: WS D11.5 gave PRD 0.2257 vs live 0.2645). Cleared right after with calib_dbh.
    s.variant isa SoutheastAlaska && (c.calib_dbh = saved_dbh)
    # AVH (AVHT40 top height) is NOT backdated during calibration: FVS's DENSE backdating pass
    # updates BA/point_ba/PCT at the past dbh, but the calibration DGF's relative-height term
    # reads the CURRENT-stand AVH (like the current point_ba restored at line 347 and the NE
    # current-stand BADIST above). compute_density! above recomputes AVH at the backdated dbh
    # ranking (too low → relht too high → WK2 too high → COR under-shrunk); restore the current
    # value for the prediction. Validated bit-exact vs live FVSsn (COR 1.548, full .sum) on FIA
    # stand 160545945010854 (missing SITE_INDEX, HB/CW/WI increment cores).
    _saved_avh = s.plot.avg_height
    s.plot.avg_height = _cur_avh
    dgf!(s, s.variant)                        # WK2 = DGF prediction at the PAST stand (variant dgf)
    _TT_CUR_RMSQD[] = -1.0f0     # #191/#195: clear the current-RMSQD stash (actual growth uses stand_qmd); unconditional to avoid leaks
    s.plot.avg_height = _saved_avh
    c.calib_dbh = Float32[]
    s.plot.forest_type = saved_fortype
    wk2 = view(s.scratch.wk, 2, :)

    # calibration VMLT (autcor LSTART: new=old=YR, the measurement base period — 5 for SN, 10 for NE,
    # NOT the cycle length FINT, which TIMEINT can change; the projection VMLT below uses FINT). YR is
    # the variant's DG model native period (dgdriv.f VMLTYR; NE's YR=10, blkdat DATA YR/10.0/ — SIGMAR
    # is on a 10-yr basis). Hardcoding 5 over-counts NE's vardg (vrnext(5)≈11 vs vrnext(10)≈29.4 ⇒ ssigma
    # 84% high ⇒ a ~0.5% serial-correlation DG error). SN unchanged (htg_period(Southern)=5).
    # SERLCORR can override the ARMA(1,1) phi/theta; recompute BJRHO only when it does.
    bjrho = _stand_bjrho(s)
    yrcal = Int(htg_period(s.variant))
    _, vmlt = autcor(yrcal, yrcal, bjrho); c.vmlt = vmlt

    # per-species DBH range + endpoint predictions over measured trees. The GST (growth-sample-tree)
    # backdated-DBH floor is VARIANT-specific: SN/NE exclude DBH < 3.0 (sn/ne dgdriv.f:392-396), CS
    # excludes DBH < 5.0 (cs/dgdriv.f:380 `IF(WK3.LT.5.0...)`). A too-low floor over-counts the
    # calibration sample (FN) so a species that FVS leaves uncalibrated (FN < FNMIN=5 ⇒ COR=0) gets a
    # spurious COR — exactly the cst01 WO over-growth (debug-stamped: live FN[WO]=2, COR=0).
    # GST eligibility floor (BKPT): SN/NE hardcode WK3<3.0 (sn/ne dgdriv.f:384), CS/LS WK3<5.0. CR uses the
    # PER-SPECIES BREAK(ISPC) (cr/dgdriv.f:410 BKPT=BREAK(ISPC); e.g. ES=1.0) — the flat 3.0 dropped small GSTs
    # (backdated WK3 2.5" ES) from the calibration, shrinking FN and skewing COR → uniform over-growth.
    gst_min = (s.variant isa CentralStates || s.variant isa LakeStates) ? 5f0 : 3f0   # cs/ls dgdriv.f WK3<5.0
    break_cr = _cr_cal ? sd[:st_break] : nothing
    dn = fill(999f0, MAXSP); dx = zeros(Float32, MAXSP)
    pn = zeros(Float32, MAXSP); px = zeros(Float32, MAXSP)
    @inbounds for i in 1:t.n
        sp = t.species[i]
        bkpt = break_cr === nothing ? gst_min : break_cr[sp]
        (t.dbh[i] < bkpt || t.diam_growth[i] <= 0f0) && continue
        if t.dbh[i] < dn[sp]; dn[sp] = t.dbh[i]; pn[sp] = exp(wk2[i]); end
        if t.dbh[i] > dx[sp]; dx[sp] = t.dbh[i]; px[sp] = exp(wk2[i]); end
    end

    # per-species sums; remember each measured tree's residual
    reslog_t = zeros(Float32, t.n); measured = falses(t.n)
    spopn = zeros(Float32, MAXSP); spopx = zeros(Float32, MAXSP)
    dev = zeros(Float32, MAXSP); devsq = zeros(Float32, MAXSP); fn = zeros(Float32, MAXSP)
    snp = zeros(Float32, MAXSP); snx = zeros(Float32, MAXSP); sny = zeros(Float32, MAXSP)
    snxx = zeros(Float32, MAXSP); snxy = zeros(Float32, MAXSP)
    @inbounds for i in 1:t.n
        sp = t.species[i]; wk3 = t.dbh[i]; dg = t.diam_growth[i]; p = t.tpa[i]
        (wk3 < dn[sp] || wk3 > dx[sp]) && continue
        edds = exp(wk2[i]); spopn[sp] += p; spopx[sp] += edds * p
        dg <= 0f0 && continue
        bark = _cr_cal ? cr_bratio(sd, Int(sp), saved_dbh[i], _cr_cal_imod) :
               _tt_cal ? tt_bratio(Int(sp), saved_dbh[i]) :
               _bm_cal ? bm_bratio(sd, Int(sp), saved_dbh[i]) :
               _ci_cal ? ci_bratio(sd, Int(sp), saved_dbh[i]) :
               _ak_cal ? ak_bratio(Int(sp), saved_dbh[i]) :
               _wc_cal ? wc_bratio(sd, Int(sp), saved_dbh[i]) :
               _pn_cal ? wc_bratio(sd, Int(sp), saved_dbh[i]) :
               (_ec_cal || _ca_cal) ? wc_bratio(sd, Int(sp), saved_dbh[i]) :
               _so_cal ? so_bratio(sd, Int(sp), saved_dbh[i]) :
               _op_cal ? op_bratio(Int(sp), saved_dbh[i]) :   # OP: op_bratio (op/bratio.f) — the TERM bark MUST match
                                                              # (shared bark_a/bark_b=0 floored to 0.80 ⇒ TERM low ⇒ RESLOG −0.098 ⇒ WF COR flips negative)
               _bc_cal ? bc_bratio(Int(sp)) :                 # BC: constant BARK1 (shared bark_a/bark_b=0 ⇒ 0.80 floor, wrong)
               bark_ratio(bark_a, bark_b, sp, saved_dbh[i])   # bark at CURRENT dbh (dgdriv.f:435)
        term = dg * (2f0 * bark * wk3 + dg) * scale
        term <= 0f0 && continue
        reslog = log(term) - wk2[i]
        reslog_t[i] = reslog; measured[i] = true
        fn[sp] += 1f0; dev[sp] += reslog; devsq[sp] += reslog^2
        snp[sp] += p; snx[sp] += p * edds; sny[sp] += p * reslog
        snxx[sp] += p * edds^2; snxy[sp] += p * reslog * edds
    end

    # per-species COR / SIGMA / regression line / VARDG
    slop = zeros(Float32, MAXSP); bnx = zeros(Float32, MAXSP); bny = zeros(Float32, MAXSP)
    calibrated = falses(MAXSP)
    @inbounds for sp in 1:MAXSP
        c.sigma[sp] = sigmar[sp]                      # SIGMA=SIGMAR unless calibrated (dgdriv.f:196)
        # NOCALIB (LDGCAL=.FALSE.) suppresses the self-calibration for a species: skip the COR
        # fit so dg_cor / dg_cor_goal stay 0 (dgdriv.f:567 — correction terms not scaled).
        if s.control.dg_calib_sp[sp] && isct[sp, 1] != 0 && fn[sp] >= fnmin && snp[sp] > 0f0
            bnxv = snx[sp] / snp[sp]; bnyv = sny[sp] / snp[sp]
            csnxy = snxy[sp] - bnxv * bnyv * snp[sp]
            csnxx = snxx[sp] - bnxv * bnxv * snp[sp]
            if csnxx >= 0f0
                slp = csnxy / csnxx
                bpopx = spopx[sp] / spopn[sp]
                sdpred = sqrt(csnxx / (snp[sp] * (1f0 - 1f0 / fn[sp])))
                dist = abs(bpopx - bnxv) / sdpred
                regcor = bnyv + (bpopx - bnxv) * slp
                cornew = dist > 3f0 ? bnyv : dist <= 1f0 ? regcor :
                         bnyv * (dist / 2f0) + regcor * (1f0 - dist / 2f0)
                svar = devsq[sp] - dev[sp]^2 / fn[sp]
                svar_v = (svar / (fn[sp] - 1f0)) / fn[sp]
                psigsq = s.variant isa Northeast ? 0.0898f0 :
                         s.variant isa CentralRockies ? 0.07f0 :
                         s.variant isa Kootenai ? KT_PSIGSQ[sp] :
                         s.variant isa EasternMontana ? EM_PSIGSQ[sp] :
                         s.variant isa Teton ? TT_PSIGSQ[sp] :
                         s.variant isa Utah ? UT_PSIGSQ[sp] :
                         s.variant isa BlueMountains ? BM_PSIGSQ[sp] :
                         s.variant isa BritishColumbia ? BC_PSIGSQ[sp] :
                         s.variant isa CentralIdaho ? CI_PSIGSQ[sp] :
                         s.variant isa InlandEmpire ? IE_PSIGSQ[sp] :
                         s.variant isa Klamath ? NC_PSIGSQ[sp] :
                         s.variant isa Olympic ? 0.0898f0 : DG_PSIGSQ   # OP 0.0898 (op/dgdriv.f DATA PSIGSQ/MAXSP*0.0898/) / NE 0.0898 / SN default
                temp = min(cornew * cornew / psigsq, 72f0)
                wc = 1f0 / (1f0 + exp(-0.5f0 * temp) * sqrt(svar_v / psigsq))
                corv = wc * cornew
                # out-of-range trap (cortem = exp(COR))
                if exp(corv) < 0.0821f0 || exp(corv) > 12.1825f0
                    corv = 0f0
                end
                c.dg_cor[sp] = corv
                # modified residual SD only for calibrated species (dgdriv.f:323-325)
                c.sigma[sp] = sqrt((svar + c.atten[sp] * sigmar[sp]^2) / (fn[sp] + c.atten[sp]))
                # FVS_CalibStats (dgdriv.f:896-962): capture NUMCAL/STDRAT/WC for the CALBSTDB DBS report. STDRAT =
                # SQRT(SVAR/SIGMAR²) with the SAMPLE variance SVAR/(FN−1); WC is the weight computed above.
                c.cal_ntree[sp] = round(Int32, fn[sp])
                c.cal_stdrat[sp] = sigmar[sp] > 0f0 && fn[sp] > 1f0 ?
                                   sqrt((svar / (fn[sp] - 1f0)) / sigmar[sp]^2) : 0f0
                c.cal_wci[sp] = wc
                c.cal_cortem[sp] = exp(corv)      # CORTEM = EXP(COR) at calibration time (before any CORMLT re-scale)
                slop[sp] = slp; bnx[sp] = bnxv; bny[sp] = bnyv; calibrated[sp] = true
            end
        end
        vtemp = exp(c.sigma[sp]^2)
        c.vardg[sp] = (vtemp - 1f0) * vtemp / vmlt
    end

    # seed OLDRN in FVS species-sorted RNG order
    oldrn = t.old_random
    @inbounds for sp in 1:MAXSP
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        if calibrated[sp]
            rx = bny[sp] + (px[sp] - bnx[sp]) * slop[sp]
            rn = bny[sp] + (pn[sp] - bnx[sp]) * slop[sp]
            for k in i1:i2
                i = ind1[k]
                if measured[i]
                    oldrn[i] = reslog_t[i]
                else
                    oldrn[i] = bny[sp] + (exp(wk2[i]) - bnx[sp]) * slop[sp]
                    t.dbh[i] < dn[sp] && (oldrn[i] = rn)
                    t.dbh[i] > dx[sp] && (oldrn[i] = rx)
                end
            end
        else
            bound = s.control.dg_stddev_bound * c.sigma[sp]
            for k in i1:i2
                i = ind1[k]
                z = 0f0
                while true
                    z = bachlo(s.rng, 0f0, c.sigma[sp])
                    z <= bound && break
                end
                oldrn[i] = z
            end
        end
    end

    # COR attenuation goal (WCI) + clamp OLDRN to ±DGSD·SIGMA
    @inbounds for sp in 1:MAXSP
        c.dg_cor_goal[sp] = 0.5f0 * c.dg_cor[sp]
    end
    @inbounds for i in 1:t.n
        lim = s.control.dg_stddev_bound * c.sigma[t.species[i]]
        oldrn[i] > lim && (oldrn[i] = lim)
        oldrn[i] < -lim && (oldrn[i] = -lim)
    end

    # Small-tree height-growth calibration: HCOR_init (regent.f:411-516). For each
    # LHTCAL species (all by default) regress measured-vs-predicted small-tree (dbh<5)
    # height growth: HCOR_init = ln( Σ(HTG·P) / Σ(EDH·P) ) when ≥ NCALHT(5) observations,
    # where EDH = HTCALC-predicted increment (mode 0 age → mode 9 HTGR, ≥0.1) and the
    # measured HTG is scaled by SCALE3 = REGYR/FINTH (regent.f:406/462; REGYR=5, FINTH = the
    # GROWTH keyword's height measurement period, default 5 → SCALE3 = 1; RHCON = 1, no HCOR2). This
    # initial value seeds DIFH = HCOR_init − WCI; the per-cycle attenuation (diameter_growth!)
    # rides it like COR. Uses the CURRENT diameters (saved_dbh) for the dbh<5 filter, since
    # the regent calibration is independent of the large-tree diameter backdating.
    # SN-only: the small-tree height calibration uses the SN ht_curve (HTCALC) coefs + the REGENT
    # model. NE's height-dbh/growth model differs (htdbh + the NE htgf), so its small-tree height
    # calibration is a separate NE piece — skip this block for NE (htg_cor_init stays 0).
    s.variant isa Southern && let bc = (sd[:ht_curve_b1], sd[:ht_curve_b2], sd[:ht_curve_b3], sd[:ht_curve_b4], sd[:ht_curve_b5]),
        montane = !isempty(s.plot.eco_unit) && s.plot.eco_unit[1] == 'M'
        ncalht = 5
        scale3 = s.control.growth_finth > 0f0 ? 5f0 / s.control.growth_finth : 1f0  # REGYR/FINTH
        @inbounds for sp in 1:MAXSP
            i1 = isct[sp, 1]; i1 == 0 && continue
            i2 = isct[sp, 2]
            si = s.plot.sp_site_index[sp]
            snp = 0f0; snx = 0f0; sny = 0f0; nh = 0
            for k in i1:i2
                i = ind1[k]
                saved_dbh[i] >= 5f0 && continue            # large trees excluded
                hh = t.height[i] - t.ht_growth[i]          # start-of-period height (IHTG<2)
                hh < 0.01f0 && continue
                t.ht_growth[i] < 0.001f0 && continue       # no measured height growth
                aget = htcalc_age(bc, sp, si, hh, montane)
                htgr = htcalc_incr(bc, sp, si, aget, montane)
                htgr < 0.1f0 && (htgr = 0.1f0)
                edh = htgr                                  # ·RHCON(=1); EDH≥0.1 already
                p = t.tpa[i]
                term = t.ht_growth[i] * scale3            # TERM = HTG·SCALE3 (regent.f:462)
                snp += p; snx += edh * p; sny += term * p; nh += 1
            end
            nh < ncalht && continue
            cornew = sny / snx                              # (ΣHTG·P)/(ΣEDH·P)
            cornew <= 0f0 && (cornew = 1f-4)
            (cornew < 0.0821f0 || cornew > 12.1825f0) && (cornew = 1f0)
            c.htg_cor_init[sp] = log(cornew)
        end
    end

    # TT UTVAR regent-height calibration (PM/UJ/RM): HCOR = ln(Σ(HTG·SCALE3·P)/Σ(EDH·P)),
    # EDH = POTHTG·PCTRED·VIGOR·RHCON·0.5 (regent.f REGCON 1299-1356). RHCON=1. SCALE3=REGYR/FINTH=10/5=2.
    if s.variant isa Teton
        scale3_tt = s.control.growth_finth > 0f0 ? 10f0 / s.control.growth_finth : 2f0
        avh = stand_top_height(s); relden_c = stand_ccf(s)
        xp = avh * (relden_c / 100f0); xp > 300f0 && (xp = 300f0)
        pctred_c = 1.11436f0 + xp * (-0.011493f0 + xp * (0.43012f-4 + xp * (-0.72221f-7 +
                   xp * (0.5607f-10 - xp * 0.1641f-13))))
        pctred_c > 1f0 && (pctred_c = 1f0); pctred_c < 0.01f0 && (pctred_c = 0.01f0)
        for sp in (4, 11, 12, 13, 16)
            isct[sp, 1] == 0 && continue
            i1 = isct[sp, 1]; i2 = isct[sp, 2]; sitear = s.plot.sp_site_index[sp]
            snx = 0f0; sny = 0f0; nh = 0
            for k in i1:i2
                i = ind1[k]
                hg = t.ht_growth[i]; hg < 0.001f0 && continue      # measured HTG (observed)
                t.height[i] < 0.01f0 && continue
                cr = Float32(t.crown_pct[i]); x = cr / 100f0
                # ★ REGCON POTHTG uses H=0 (evaluated once, NOT the tree's current H): (SJ·1.5−0)/(SJ·1.5)=1
                # ⇒ POTHTG=(SJ/5)·0.83; AND the ·0.5 (UT 10yr→5yr) IS applied. (My earlier current-H/no-0.5
                # form coincidentally matched PM because H≈½·SJ·1.5, but broke UJ where H>½·SJ·1.5.)
                pothtg = (sitear / 5f0) * 0.83f0
                vigor = (150f0 * x * x * x * exp(-6f0 * x)) + 0.3f0; vigor > 1f0 && (vigor = 1f0)
                vigor = 1f0 - ((1f0 - vigor) / 3f0)
                edh = pothtg * pctred_c * vigor * 0.5f0            # ·RHCON=1 · 0.5 (regent.f:1304)
                p = t.tpa[i]; snx += edh * p; sny += hg * scale3_tt * p; nh += 1
            end
            nh < 5 && continue                                     # NCALHT
            cornew = snx > 0f0 ? sny / snx : 1f0; cornew <= 0f0 && (cornew = 1f-4)
            hc = log(cornew)                                        # raw HCOR (regent.f:1364)
            # dgdriv.f:213 attenuation. PM/UJ/RM have NO dgf DG COR (regent-DG) ⇒ dg_cor=0 ⇒ WCI=0 ⇒
            # HCOR_used = CORMLT·HCOR_raw, CORMLT=exp(−0.02773·SFINT), SFINT=YR=10 (validated: 0.758·1.83=1.39).
            cormlt = exp(-0.02773f0 * htg_period(s.variant))
            c.htg_cor_init[sp] = hc; c.htg_cor_small[sp] = cormlt * hc
        end
    end

    # IE regent small-tree HEIGHT calibration (ie/regent.f:1138-1337): compute the RAW regent HCOR into
    # htg_cor_init for NIVAR species via the NIVAR EDH model. Without it, IE NIVAR species had htg_cor_init=0,
    # so the shared attenuation below leaked the diameter COR into the regent height CON ⇒ small-tree
    # over-growth on dense stands where the calibration fires (#171). Inert on iet01 (no measured small-tree HTG).
    s.variant isa InlandEmpire && ie_regent_hcor_init!(s, isct, ind1, saved_dbh)

    # The CS/NE regent HCOR calibration's BALMOD reads the BACKDATED-dbh stand BA (live regent.f BA=177.5,
    # the backdated value, NOT the restored current 242). FVS DENSE (dense.f:79-86) sums the backdated BA over
    # LIVE + the RECENTLY-DEAD records (trees that died within the measurement period, added back at their dbh);
    # `s.plot.basal_area` is live-only (169.5), so sum live+dead here to match CRATET's 177.515 exactly.
    bd_ba_hcor = 0f0; bd_sd2_hcor = 0f0; bd_tpa_hcor = 0f0
    @inbounds for i in 1:(t.n + t.ndead)
        dd = t.dbh[i]; bd_ba_hcor += dd * dd * t.tpa[i] * 0.005454154f0
        bd_sd2_hcor += dd * dd * t.tpa[i]; bd_tpa_hcor += t.tpa[i]
    end
    rmsqd_bd_hcor = bd_tpa_hcor > 0f0 ? sqrt(bd_sd2_hcor / bd_tpa_hcor) : 0f0   # backdated stand QMD (RMSQD basis for LS balmod)
    # The regent HCOR calibration reads the BACKDATED-stand PCT for BALMOD's BAL. It uses the CRATET DENSE
    # percentile, which INCLUDES ALL recently-dead records incl. history-8 (dense.f keeps WK3=dbh; only IMC=9 is
    # zeroed) — UNLIKE the dgf percentile above (which zeroes history-8). So recompute the percentile here over
    # live + ALL dead (backdated dbh), total incl. the history-8 tree, matching CRATET's 177.5 BA basis exactly.
    bd_pct_hcor = zeros(Float32, t.n)
    let nlive2 = t.n, ntot = t.n + t.ndead
        rankd = Float32[i <= nlive2 ? saved_dbh[i] : t.dbh[i] for i in 1:ntot]   # current-dbh rank order
        wk5h  = Float32[t.dbh[i]^2 * t.tpa[i] for i in 1:ntot]                    # backdated BA, ALL dead incl hist-8
        ordh = sortperm(rankd; rev = true); toth = sum(wk5h); cumh = 0f0
        if toth > 0f0
            @inbounds for k in ntot:-1:1
                ii = ordh[k]; cumh += wk5h[ii]
                ii <= nlive2 && (bd_pct_hcor[ii] = cumh / toth * 100f0)
            end
        end
    end
    # restore current diameters + current-stand density (the backdating was local)
    @inbounds for i in 1:t.n; t.dbh[i] = saved_dbh[i]; end
    compute_density!(s)
    # NE small-tree HCOR height calibration (ne/regent.f:411-547). The Southern block above is SN-model-specific
    # (HTCALC ht_curve + SN REGYR=5); NE uses the NC-128 ne_htcalc + BALMOD·RELHTA and REGYR=10. Runs on the
    # CURRENT (restored) dbh/density — regent uses the current dbh, not the DG-backdated one. Each LHTCAL species
    # (all by default) regresses measured-vs-predicted small-tree (dbh<5) height growth: HCOR_init =
    # ln(Σ(HTG·SCALE3·P)/Σ(EDH·P)) with ≥NCALHT(5) obs; SCALE3=REGYR/FINTH=10/FINTH (default 5 ⇒ 2). Stays 0 for
    # a stand with no measured small-tree HTG (net01 carries measured DG, not HTG) ⇒ SN/CS + net01 unaffected.
    if s.variant isa Northeast
        b3ne = sd[:dg_b3]; avh = s.plot.avg_height
        ebau = zeros(Float32, 50); ne_badist!(ebau, s)
        scale3 = s.control.growth_finth > 0f0 ? 10f0 / s.control.growth_finth : 2f0   # REGYR(10)/FINTH(default 5)
        @inbounds for sp in 1:MAXSP
            i1 = isct[sp, 1]; i1 == 0 && continue
            i2 = isct[sp, 2]; si = s.plot.sp_site_index[sp]
            snx = 0f0; sny = 0f0; nh = 0
            for k in i1:i2
                i = ind1[k]
                t.dbh[i] >= 5f0 && continue                       # large trees excluded
                hstart = t.height[i] - t.ht_growth[i]            # start-of-period H for the filter (regent.f:451)
                hstart < 0.01f0 && continue
                t.ht_growth[i] < 0.001f0 && continue             # no measured height growth
                # htcalc.f:389 guard (`IF(HTMAX-H.LE.1.) GO TO 900` ⇒ HTG1=0): a tree at/above the NC-128
                # asymptote (H ≥ HTMAX-1) has no invertible age — the ALOG argument goes negative and yields
                # NaN. FVS guards this INSIDE HTCALC so every call is protected; jl guards per-caller, and this
                # NCALHT calibration path was the one omission. Without it an off-curve small tree (e.g. a
                # 47-ft, 2.9" tree) NaNs `cornew` ⇒ poisons the species-level htg_cor_init ⇒ NaN growth for all
                # small trees of that species (an SN-family crash surfaced only on a rare FIA stand under fire).
                if ne_htcalc_htmax(sp, si) - t.height[i] <= 1f0
                    htgr = 0f0
                else
                    aget = ne_htcalc_age(sp, si, t.height[i])    # HTCALC age/incr on the CURRENT height (regent.f:466)
                    htgr = ne_htcalc_incr(sp, si, aget)
                end
                gmod = ne_balmod(b3ne[sp], ebau, t.dbh[i])       # BALMOD·RELHTA (regent.f:485-491)
                relht = avh > 0f0 ? min(t.height[i] / avh, 1f0) : 0f0
                gmod = 1f0 - (1f0 - gmod) * (1f0 - relht)
                htgr = max(htgr * gmod, 0.1f0)
                edh = max(htgr, 0.1f0)                            # ·RHCON=1
                snx += edh * t.tpa[i]
                sny += t.ht_growth[i] * scale3 * t.tpa[i]        # TERM = HTG·SCALE3
                nh += 1
            end
            nh < 5 && continue                                    # NCALHT
            cornew = sny / snx
            cornew <= 0f0 && (cornew = 1f-4)
            (cornew < 0.0821f0 || cornew > 12.1825f0) && (cornew = 1f0)
            c.htg_cor_init[sp] = log(cornew)
        end
    end
    # CS small-tree HCOR height calibration (cs/regent.f:422-540) — same structure as NE (D22) but cs_htcalc +
    # cs_balmod(b1,b2,b3, BAL, BA, d) with BAL=(1−PCT/100)·BA (PCT = t.crown_ratio, the BA percentile). REGYR=10.
    if s.variant isa CentralStates
        cb1 = sd[:balmod_b1]; cb2 = sd[:balmod_b2]; cb3 = sd[:balmod_b3]
        ba = bd_ba_hcor; avh = s.plot.avg_height    # BACKDATED BA (regent.f reads the backdated stand BA, not current)
        scale3 = s.control.growth_finth > 0f0 ? 10f0 / s.control.growth_finth : 2f0   # REGYR(10)/FINTH(default 5)
        @inbounds for sp in 1:MAXSP
            i1 = isct[sp, 1]; i1 == 0 && continue
            i2 = isct[sp, 2]; si = s.plot.sp_site_index[sp]
            snx = 0f0; sny = 0f0; nh = 0
            for k in i1:i2
                i = ind1[k]
                t.dbh[i] >= 5f0 && continue
                hstart = t.height[i] - t.ht_growth[i]
                hstart < 0.01f0 && continue
                t.ht_growth[i] < 0.001f0 && continue
                aget = cs_htcalc_age(sp, si, t.height[i])            # HTCALC age/incr on CURRENT height
                htgr = cs_htcalc_incr(sp, si, aget)
                bal = (1f0 - bd_pct_hcor[i] / 100f0) * ba            # BAL = (1−PCT/100)·BA (regent.f:450); BACKDATED PCT
                gmod = cs_balmod(cb1[sp], cb2[sp], cb3[sp], bal, ba, t.dbh[i])
                relht = avh > 0f0 ? min(t.height[i] / avh, 1f0) : 0f0
                gmod = 1f0 - (1f0 - gmod) * (1f0 - relht)
                htgr = max(htgr * gmod, 0.1f0)
                edh = max(htgr, 0.1f0)
                snx += edh * t.tpa[i]
                sny += t.ht_growth[i] * scale3 * t.tpa[i]
                nh += 1
            end
            nh < 5 && continue
            cornew = sny / snx
            cornew <= 0f0 && (cornew = 1f-4)
            (cornew < 0.0821f0 || cornew > 12.1825f0) && (cornew = 1f0)
            c.htg_cor_init[sp] = log(cornew)
        end
    end
    # LS small-tree HCOR height calibration (ls/regent.f:419-560) — same structure as CS (D22) but the LS
    # variant hooks: ls_htcalc (MAPLS/LTBHEC via _ls_htcoef) + ls_balmod(sp, d, BA, RMSQD, ...) which reads the
    # stand BA + QMD (not BAL/PCT). regent.f reads the BACKDATED stand BA/QMD (same as cs/regent.f). REGYR=10,
    # SCALE3 = 10/FINTH. Each LHTCAL species (all by default) with ≥NCALHT(5) measured small-tree (dbh<5) HTG
    # observations gets HCOR_init = ln(Σ(HTG·SCALE3·P)/Σ(EDH·P)). Stays 0 with no measured small-tree HTG.
    # This was the port's LS omission: aspen (746) with measured height growth was under-grown ~1.6× (con=1
    # instead of exp(HCOR)) ⇒ dense aspen regen under-thinned/under-grew vs live FVS.
    if s.variant isa LakeStates
        check = sd[:balmod_check]; mb1 = sd[:balmod_b1]; mb2 = sd[:balmod_b2]; mb3 = sd[:balmod_b3]
        mb4 = sd[:balmod_b4]; mc1 = sd[:balmod_c1]; mc2 = sd[:balmod_c2]; bamax1 = sd[:balmod_bamax1]
        ba = bd_ba_hcor; rmsqd = rmsqd_bd_hcor; avh = s.plot.avg_height   # BACKDATED stand BA/QMD (ls/regent.f reads
        scale3 = s.control.growth_finth > 0f0 ? 10f0 / s.control.growth_finth : 2f0   # the backdated DENSE stand; REGYR(10)/FINTH(default 5)
        # FVS regent.f:100 HTGR=0.0 (subroutine entry). The calibration loop's mode-9 HTCALC (regent.f:479) does
        # NOT reset HTGR for an at/above-asymptote tree: htcalc.f:391 `IF(HTMAX-H.LE.1.) GO TO 900` returns
        # WITHOUT setting the HTG1 output arg (label 900 = bare RETURN), so HTGR keeps the PREVIOUS calibration
        # tree's post-line-495 value. This is an FVS uninitialized-carry bug — the *growth* loop guards it with an
        # explicit `HTGR=0.10` (regent.f:208) but the *calibration* loop does not. It only bites species whose
        # small trees (dbh<5) exceed the species asymptote — e.g. tamarack (sp071, HTMAX≈27.9ft@SI25) carries tall
        # skinny dbh<5 records at H=28-34. A faithful drop-in must reproduce the stale carry (setting htgr=0 here
        # inflated cornew 8.888 vs live 4.09 ⇒ 2-3× DG over-growth). The carry propagates across trees AND species
        # (HTGR is never re-zeroed inside DO 100), and is updated for EVERY dbh<5,H>0.01 tree — including
        # measured-HTG=0 trees (regent.f:503 skips them from SNX/SNY only AFTER HTGR is computed at 495/498).
        htgr_carry = 0f0
        @inbounds for sp in 1:MAXSP
            i1 = isct[sp, 1]; i1 == 0 && continue
            i2 = isct[sp, 2]; si = s.plot.sp_site_index[sp]
            snx = 0f0; sny = 0f0; nh = 0
            for k in i1:i2
                i = ind1[k]
                t.dbh[i] >= 5f0 && continue                       # large trees excluded (regent.f:454, pre-HTCALC ⇒ no carry)
                hstart = t.height[i] - t.ht_growth[i]            # start-of-period H for the filter (regent.f:453)
                hstart < 0.01f0 && continue                       # (regent.f:454, pre-HTCALC ⇒ no carry update)
                if ls_htcalc_htmax(sp, si) - t.height[i] <= 1f0  # htcalc.f:391 asymptote guard ⇒ HTCALC leaves HTGR stale
                    htgr = htgr_carry                             # STALE carry (FVS bug), NOT 0
                else
                    aget = ls_htcalc_age(sp, si, t.height[i])    # HTCALC age/incr on the CURRENT height (regent.f:465)
                    htgr = ls_htcalc_incr(sp, si, aget)
                end
                gmod = ls_balmod(sp, t.dbh[i], ba, rmsqd, check, mb1, mb2, mb3, mb4, mc1, mc2, bamax1)
                relht = avh > 0f0 ? min(t.height[i] / avh, 1f0) : 0f0
                gmod = 1f0 - (1f0 - gmod) * (1f0 - relht)
                htgr = max(htgr * gmod, 0.1f0)                    # regent.f:495 (HTGR·GMOD) + 498 (floor 0.1)
                htgr_carry = htgr                                 # FVS keeps HTGR for the next tree (the stale carry)
                t.ht_growth[i] < 0.001f0 && continue              # no measured HTG (regent.f:503, AFTER HTGR ⇒ carry already set)
                edh = htgr                                        # EDH = max(HTGR·RHCON,0.1) = htgr (regent.f:500-501)
                snx += edh * t.tpa[i]
                sny += t.ht_growth[i] * scale3 * t.tpa[i]        # TERM = HTG·SCALE3
                nh += 1
            end
            nh < 5 && continue                                    # NCALHT
            cornew = sny / snx
            cornew <= 0f0 && (cornew = 1f-4)
            (cornew < 0.0821f0 || cornew > 12.1825f0) && (cornew = 1f0)
            c.htg_cor_init[sp] = log(cornew)
        end
    end
    # CR small-tree REGENT height calibration (cr/regent.f REGCAL:445-606). EDH = POTHTG·PCTRED·VIGOR·RHCON
    # (GENGYM potential HTG), aspen/paper-birch (sp20/28) use the Sheppard curve. HCOR_init = ln(Σ(HTG·SCALE3·P)/
    # Σ(EDH·P)) with ≥NCALHT(5) measured dbh<5 HTG. Runs on the CURRENT restored stand (regent uses current dbh).
    # Without this jl held con=1.0 (HCOR=0) ⇒ small-tree height under-grew (sp5 WF con 1.0 vs live 1.047).
    if s.variant isa CentralRockies
        htadj = sd[:st_htadj]; lo = sd[:site_lo]; hi = sd[:site_hi]
        scale3 = s.control.growth_finth > 0f0 ? 10f0 / s.control.growth_finth : 2f0   # REGYR(10)/FINTH(default 5)
        ccf = _cr_bd_ccf; avht = _cr_bd_avht                                            # PCTRED: BACKDATED-window CCF+AVH (dense.f RELDM1+AVH, both dead-inclusive; regent.f:466 X=AVH·RELDEN/100)
        xd = avht * (ccf / 100f0); xd > 300f0 && (xd = 300f0)
        pctred = _CR_AB[1] + xd*(_CR_AB[2] + xd*(_CR_AB[3] + xd*(_CR_AB[4] + xd*(_CR_AB[5] + xd*_CR_AB[6]))))
        pctred > 1f0 && (pctred = 1f0); pctred < 0.01f0 && (pctred = 0.01f0)
        @inbounds for sp in 1:MAXSP
            i1 = isct[sp, 1]; i1 == 0 && continue
            i2 = isct[sp, 2]
            si = s.plot.sp_site_index[sp]
            si > hi[sp] && (si = hi[sp]); si <= lo[sp] && (si = lo[sp] + 0.5f0)
            relsi = (si - lo[sp]) / (hi[sp] - lo[sp]); rsimod = 0.5f0 * (1f0 + relsi)
            pothtg = s.plot.sp_site_index[sp] / (15f0 - 4f0 * relsi) * htadj[sp]        # SITEAR (unclamped), regent.f:534
            ivf = sp in _CR_IVFLAG
            snx = 0f0; sny = 0f0; nh = 0
            for k in i1:i2
                i = ind1[k]
                t.dbh[i] >= 5f0 && continue                       # large trees excluded (regent.f:454)
                hstart = t.height[i] - t.ht_growth[i]             # start-of-period H (IHTG<2, regent.f:534)
                hstart < 0.01f0 && continue
                if sp == 20 || sp == 28                           # aspen/paper birch Sheppard curve (regent.f:542-549)
                    # AG1 = INVERSE Sheppard from the start height H (regent.f:542), NOT birth_age — the
                    # calibration path derives the age from height (the growth path uses ABIRTH; the REGCAL
                    # path inverts the curve). Using birth_age gave a ~30% wrong EDH ⇒ con over-shrunk (0.14
                    # vs live ~0.19) ⇒ aspen regen height under-grew.
                    ag1 = fpow(hstart * 12f0 * 2.54f0 / 26.9825f0, 0.8509f0)
                    ag2 = ag1 + 10f0
                    h2 = 26.9825f0 * fpow(ag2, 1.1752f0) / (2.54f0 * 12f0)
                    edh = (h2 - hstart) * rsimod * 0.75f0
                    edh < 0f0 && (edh = 0f0)
                else
                    xc = Float32(t.crown_pct[i]) / 100f0          # VIGOR from crown ratio (regent.f:531)
                    vigor = 150f0 * fpow(xc, 3f0) * fexp(-6f0 * xc) + 0.3f0
                    vigor > 1f0 && (vigor = 1f0)
                    ivf && (vigor = 1f0 - (1f0 - vigor) / 3f0)
                    edh = pothtg * pctred * vigor                 # ·RHCON(=1), regent.f:536
                end
                t.ht_growth[i] < 0.001f0 && continue              # no measured HTG (regent.f:551)
                snx += edh * t.tpa[i]
                sny += t.ht_growth[i] * scale3 * t.tpa[i]         # TERM = HTG·SCALE3 (regent.f:552)
                nh += 1
            end
            nh < 5 && continue                                    # NCALHT
            cornew = sny / snx
            cornew <= 0f0 && (cornew = 1f-4)
            (cornew < 0.0821f0 || cornew > 12.1825f0) && (cornew = 1f0)
            c.htg_cor_init[sp] = log(cornew)
        end
    end
    # IE PI/JU (sp15/16) small-tree height calib: NO regent HCOR block needed, and it is ALREADY CORRECT.
    # sp15's own regent REGCAL gets N=4 < NCALHT(5) and skips (htg_cor_init=0 — jl matches live). Live's
    # HCOR(15)=0.1436 is a SPILLOVER from the DIAMETER COR (dgdriv.f:183-205: WCI=0.5·COR, HCOR=WCI·(1−CORMLT)),
    # which the htg_cor_small formula below (dg_cor_goal + cormlt_h·(htg_cor_init − dg_cor_goal)) already
    # implements. VERIFIED: jl dg_cor[15]=1.186 / dg_cor_goal=0.593 (=live WCI) ⇒ htg_cor_small[15]=0.14366 ⇒
    # regent CON=1.0·exp=1.15449 == live CON 1.1545 BIT-EXACT. PI/JU is bit-exact-or-cornered: on a PHYSICAL
    # pinyon stand the deterministic HTGR1 (pre-ZZRAN height growth) jl==live to ±0.0001 in the growing regime;
    # per-tree DG + CON also match. The .sum spread is the accepted ZZRAN + tripling residual (ch9) — AMPLIFIED
    # for pinyon because the linear DK maps the ZZRAN-perturbed H→DBH (live draws per-triple-copy, jl once+triples).
    # NOT a growth/calib/crown/mortality bug. See [[fvsjl-ie-variant-port]].
    # BC V2 (LV2ATV) small-tree HEIGHT calibration (regent.f REGCAL:1783-1922). HCOR = ln(Σ(HTG·SCALE3·P)/
    # Σ(EDH·P)) over N≥NCALHT(5) trees with measured HTG>0.001, DBH<5, backdated H≥1.37m. NPER=1 (FINTH=5 ⇒
    # single-subcycle EDH, no density projection): EDH = exp(RHCON + RHLH·ln(H_bd) + RHCCF·RELDEN + RHBAL·BAL),
    # BAL = BA·(100−PCT)·0.0001. HCOR=htg_cor_init; the attenuation below (htg_cor_small) applies it. Without
    # this jl's htg_cor stayed 0 ⇒ V2 small-tree under-grew (essf HCOR(14)=0.5439). VERIFIED via FVSbc RGHCOR14.
    if s.variant isa BritishColumbia
        _bc_zone, _bc_series = bc_stand_zone(s)
        if bc_lv2atv(_bc_zone)
            ba_c = _bc_bd_ba; relden_c = _bc_bd_relden       # BACKDATED stand (REGCAL uses the past-period stand)
            _bc_pct = _bc_bd_pct                              # BACKDATED percentile
            scale3 = s.control.growth_finth > 0f0 ? 5f0 / s.control.growth_finth : 1f0
            regch = bc_v2_regch(s.plot.aspect, s.plot.slope)
            @inbounds for sp in 1:MAXSP
                i1 = isct[sp, 1]; i1 == 0 && continue
                i2 = isct[sp, 2]
                rhcon_sp = bc_v2_rhcon(sp, regch)
                snx = 0f0; sny = 0f0; nh = 0
                for k in i1:i2
                    i = ind1[k]
                    t.dbh[i] >= 5f0 && continue
                    h_bd = t.height[i] - t.ht_growth[i]           # backdated H (IHTG<2)
                    h_bd < 1.37f0 * 3.28084f0 && continue         # H ≥ 1.37 m (regent.f:1809; MtoFT)
                    t.ht_growth[i] < 0.001f0 && continue          # measured HTG present
                    bal = ba_c * (100f0 - _bc_pct[i]) * 0.0001f0
                    edh = exp(rhcon_sp + BC_RG_V2_RHLH[sp]*log(h_bd) +
                              BC_RG_V2_RHCCF[sp]*relden_c + BC_RG_V2_RHBAL[sp]*bal)
                    snx += edh * t.tpa[i]
                    sny += t.ht_growth[i] * scale3 * t.tpa[i]
                    nh += 1
                end
                nh < 5 && continue                                # NCALHT
                cornew = sny / snx
                cornew <= 0f0 && (cornew = 1f-4)
                (cornew < 0.0821f0 || cornew > 12.1825f0) && (cornew = 1f0)
                c.htg_cor_init[sp] = log(cornew)
            end
        end
    end
    return s
end

"""
    diameter_growth!(state, ::Southern; sfint=5f0)

Variant hook: compute each tree's periodic diameter growth into `trees.diam_growth`
(DGDRIV growth mode, sn/dgdriv.f). `DG = sqrt(d_ib² + exp(ln DDS)·frm) − d_ib`,
where `frm = exp(raw)`. While tripling is active (snt01 early cycles, `tripling=true`)
the serial-correlation factor is DETERMINISTIC — `raw = FM·ssigma·rhocp + corr·OLDRN`
(dgdriv.f:90,117, FM=−0.14228); once tripling stops it is the stochastic `dgscor!`
(BACHLO + AR(1)). COR is attenuated each cycle toward its goal (WCI). Trees are
walked in species-sorted order to keep any RNG draws bit-exact. `sfint = IY[icyc+1]−IY[1]`.
"""
const DG_FM = -0.14228f0      # tripling mid-record variance factor (dgdriv.f FM)
const DG_FU =  1.271f0        # tripling upper-record factor (dgdriv.f FU)
const DG_FL = -1.549f0        # tripling lower-record factor (dgdriv.f FL)

function diameter_growth!(s::StandState, ::AbstractVariant; sfint::Float32 = 5f0,
                          tripling::Bool = true)
    t, c = s.trees, s.calib
    sd = s.coef.species
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    # CR uses its own IMAP bark (cr_bratio), not the linear calib.bark_a/b (which are 0 for CR ⇒ the
    # DDS→DG conversion would use the 0.80 floor, understating inside-bark dia ⇒ over-high DG). The
    # DDS→DG bark MUST match the bark cr_gemdg used internally.
    _cr_dg = s.variant isa CentralRockies
    _cr_imodty = _cr_dg ? Int(s.plot.model_type) : 0
    _tt_dg = s.variant isa Teton   # TT bark = tt_bratio (PP sp10 IMAP=4 power model); DDS→DG dib must match
    _bc_dg = s.variant isa BritishColumbia   # BC bark = bc_bratio (constant; calib.bark_a/b=0 ⇒ 0.80 floor otherwise)
    _bm_dg = s.variant isa BlueMountains     # ★#140: BM bark = bm_bratio (POWER); the linear fallback here gave
                                             # ~0.99 (no linear coeffs) vs POWER ~0.86 ⇒ d_ib too large ⇒ DDS→DG
                                             # 7-8% LOW on EVERY tree ⇒ dq10 low ⇒ self-thin under-kill (the #140
                                             # under-thin). Calibration (line 279) + mortality (_mbark) already use
                                             # bm_bratio; this DDS→DG apply site was the missing branch.
    _ak_dg = s.variant isa SoutheastAlaska   # AK bark = ak_bratio (3-type); shared bark_a/bark_b unset ⇒ 0.80 floor otherwise
    _ci_dg = s.variant isa CentralIdaho      # ★ same class as #140: CI bark = ci_bratio (POWER). Linear fallback
                                             # gave a FLAT ~0.90 vs ci_bratio's per-sp/dbh 0.88-0.93 ⇒ DDS→DG off
                                             # ~2% on species where they diverge (net ~0.3%, small since 0.90 ≈ CI
                                             # POWER bark). Calibration (line 280 _ci_bd) + mortality already use
                                             # ci_bratio; this apply site was missed (the "DDS bit-exact" check
                                             # missed the bark-converted DG, exactly as for BM).
    _wc_dg = s.variant isa WestCascades      # ★ same class as #140/CI: WC bark = wc_bratio (POWER, all imap=1). The
    _pn_dg = s.variant isa PacificNorthwest  # PN bark = wc_bratio (POWER) — the DDS→DG conversion watchpoint (from the start)
    _ec_dg = s.variant isa EastCascades      # EC bark = wc_bratio — DDS→DG conversion watchpoint (from the start)
                                             # linear fallback floors 0.80 vs wc_bratio ~0.83-0.90 ⇒ d_ib understated
    _ca_dg = s.variant isa CentralCalifornia # CA bark = wc_bratio (per-species bark_imap) — SAME class: linear-0.80 floor
    _on_dg = s.variant isa Ontario           # ON bark = on_bratio (canada/on/bratio.f, metric H/D). dgdriv.f:201
                                             # D=DBH·BRATIO(ISPC,DBH,HT) uses the ORIGINAL DBH (NOT the grown D that
                                             # dgf.f:362 used to form DDS) — the DDS→DG round-trip watchpoint: the
                                             # linear fallback (bark_a/bark_b=0) would floor to 0.80 vs on_bratio's
                                             # ~0.91-0.96, understating d_ib ⇒ over-high DG. dgf! formed WK2 with
                                             # on_bratio, so this conversion MUST use on_bratio too (same family).
    _so_dg = s.variant isa SouthCentralOregon # SO bark = so_bratio (so/bratio.f 3-path) — DDS→DG conversion
                                             # understated d_ib ⇒ same DDS gave a LARGER DG ⇒ cat01 ~+6% BA/cyc over-growth
                                             # ⇒ DDS→DG (sqrt(d_ib²+DDS)−d_ib) OVER-predicts ~2%/tree ⇒ the multi-cycle
                                             # BA/QMD over-growth (2090 BA +19%). Calibration/mortality/update already
                                             # use wc_bratio; this DDS→DG conversion was the missing branch.
    yr = htg_period(s.variant)   # DG model native period (gradd.f FINT/YR scale): 5 SN, 10 NE
    # DGBND DBH-range bounds are SN-only (NE's DGBND is just the SIZCAP cap, ne/dgbnd.f); `nothing`
    # ⇒ the per-tree bound skips the dlo/dhi adjustment and applies only the size cap.
    dlo_v = haskey(sd, :dg_bound_dbh_lo) ? sd[:dg_bound_dbh_lo] : nothing
    dhi_v = haskey(sd, :dg_bound_dbh_hi) ? sd[:dg_bound_dbh_hi] : nothing
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    oldrn = t.old_random
    nlive = t.n
    # tripling scratch: per-tree upper/lower DG + their serial-correlation residual
    do_trip = tripling && 3 * nlive + t.ndead <= length(t.dbh)
    dgU = do_trip ? Vector{Float32}(undef, nlive) : Float32[]
    dgL = do_trip ? Vector{Float32}(undef, nlive) : Float32[]
    rnU = do_trip ? Vector{Float32}(undef, nlive) : Float32[]
    rnL = do_trip ? Vector{Float32}(undef, nlive) : Float32[]

    # Attenuate COR toward the calibration goal before predicting (dgdriv.f:76-79).
    # The attenuation clock is the cumulative elapsed time SINCE the inventory (FVS
    # SFINT = IY(icyc)−IY(1)), so the first projection cycle uses the FULL COR
    # (cormlt=1 at 0 elapsed years) and it decays thereafter. Use the actual elapsed
    # years from the IY schedule (= current_cycle_year − inventory), NOT `sfint·cycle`
    # — the two are equal only for UNIFORM cycles; a TIMEINT/CYCLEAT non-uniform
    # schedule (e.g. a 10-yr cycle 2) needs the true cumulative time (5, not 10).
    # NB (verified the hard way, -1823 tests): a debug-FVS dgdriv COR dump shows COR
    # ONE CYCLE AHEAD (1.0221 at cyc1) because that WRITE fires AFTER dgdriv updates COR
    # for the NEXT cycle, while dgf already baked the CURRENT (pre-update) COR into WK2.
    # jl's dg_cor at cycle N = the value FVS USES for cycle N's DG (the pre-update one) —
    # the START clock here is correct; do NOT "fix" it to elapsed+sfint.
    elapsed = Float32(current_cycle_year(s) - Int(s.control.cycle_year[1]))
    cormlt = exp(-0.02773f0 * elapsed)
    # The REGENT small-tree height calibration HCOR rides the SAME WCI attenuation as the
    # diameter COR (dgdriv.f:188-194) but on the elapsed-at-END-of-period clock (cycle+1):
    # HCOR = WCI + cormlt_h·DIFH, DIFH = HCOR_init − WCI (set at ICYC=1). This runs for LDGCAL
    # species, and LDGCAL defaults TRUE for ALL species (grinit.f:102; only the unported
    # NOCALIB keyword turns it off), so it always runs here. For a species with NO diameter
    # calibration WCI=0, so the height calibration DECAYS as cormlt_h·HCOR_init (it is NOT held
    # constant — that was the bug); for one with no height calibration HCOR_init=0, reducing to
    # the WCI·(1−cormlt_h) progression. HCOR is SEPARATE from the large-tree HTGF term HTCON
    # (`htg_cor`, from the HCOR2 keyword, 0 for snt01). HCOR_init is computed by the regent
    # regression in `calibrate_diameter_growth!`.
    cormlt_h = exp(-0.02773f0 * (elapsed + sfint))   # elapsed at END of this period (cumulative)
    @inbounds for sp in 1:MAXSP
        c.dg_cor[sp] = c.dg_cor_goal[sp] + cormlt * c.dg_cor_goal[sp]
        c.htg_cor_small[sp] = c.dg_cor_goal[sp] + cormlt_h * (c.htg_cor_init[sp] - c.dg_cor_goal[sp])
    end

    species_sort!(s)
    dgf!(s, s.variant)
    wk2 = view(s.scratch.wk, 2, :)

    # per-cycle ARMA multipliers: AUTCOR(new, old) where `new` = THIS cycle's period and
    # `old` = the PREVIOUS cycle's period (dgdriv.f). For uniform 5-yr cycles both are 5
    # (unchanged); a non-uniform TIMEINT/CYCLEAT schedule (e.g. a 10-yr cycle following a
    # 5-yr one) needs new=10, old=5 — using the base YR for both under-grows the long cycle.
    # The FIRST cycle's `old` is the MEASUREMENT base period = the variant's YR (5 SN, 10 NE),
    # NOT a hardcoded 5 — covmlt=AUTCOR(YR,YR).covar drives CORR; using 5 for NE under-counts it
    # (corr 0.148 vs FVS 0.181 ⇒ a residual serial-correlation DG error). SN unchanged (YR=5).
    cyc = Int(s.control.cycle)
    newp = max(1, cycle_period_at(s.control, cyc))
    # Climate-FVS: scale the large-tree DDS by clgmult's TREEMULT (dgdriv.f:153 CALL CLGMULT(WK4);
    # :217 DDS=EXP(WK2)·WK4). THISYR = IY(ICYC)+FINT/2. Inert unless a CLIMATE keyword activated s.climate.
    (s.climate !== nothing && s.climate.active) &&
        apply_climate_dds!(s, wk2, Float32(current_cycle_year(s)) + Float32(newp) / 2f0)
    # The FIRST projection cycle's `old` period is the DG MEASUREMENT period (dgdriv PVMLT basis) — the GROWTH
    # keyword FINT when overridden from its universal 5-yr default, else the variant native YR (htg_period:
    # 5 SN / 10 NE). Live-stamped: growth_fint10 (GROWTH 10) ⇒ AUTCOR(new=5, old=10) CORR=0.3906, not
    # AUTCOR(5,5)=0.3196 — the ~0.4% cuft growth_fint residual (D2). Default (growth_fint=5) is unchanged for
    # both variants (SN old=5, NE old=10), so every bit-exact scenario stays bit-exact.
    # Use the EXPLICIT DG measurement period (FVS OLDFNT) when one was given — via a GROWTH card OR the
    # FIA-DB DG_MEASURE column (both set `growth_dg_set`). The old `growth_fint != 5f0` test discarded an
    # explicit 5 (indistinguishable from the universal default 5), wrongly falling back to the variant YR:
    # for LS/CS/NE (htg_period 10) an FIA stand's DG_MEASURE=5 became oldp=10, so the first-cycle serial-
    # correlation CORR used AUTCOR(NOLD=10)=0.181 instead of FVS's AUTCOR(NOLD=5)=0.148 — a ~3% over-high
    # tripled-record DG that over-projected the self-thinning QMD (live dgdriv.f DG(I) stamp; see
    # [[fvsjl-ls-morts-growth-projection-bug]]). `growth_dg_set` distinguishes an explicit 5 from the default.
    meas_fint = (s.control.growth_dg_set && s.control.growth_fint > 0f0) ?
                Int(round(s.control.growth_fint)) : Int(htg_period(s.variant))
    oldp = cyc == 0 ? meas_fint : max(1, cycle_period_at(s.control, cyc - 1))
    covmlt, vmlt = autcor(newp, oldp, _stand_bjrho(s))
    pvmlt = c.vmlt > 0f0 ? c.vmlt : vmlt
    corr = covmlt / sqrt(vmlt * pvmlt)
    c.vmlt = vmlt   # FVS dgdriv.f:116 PVMLT=VMLT carry (uniform 5-yr unaffected)
    # BAIMULT (MULTS kind 1): per-species diameter-growth multiplier scaling DDS
    # (dgdriv.f XDGROW=ln(XDMULT) added to ln(DDS) ⇒ DDS·XDMULT).
    cur_year = current_cycle_year(s)   # IY schedule (TIMEINT/CYCLEAT-aware)

    @inbounds for sp in 1:MAXSP
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        xbai = active_multiplier(s.control, :bai, sp, cur_year)
        xdgrow = flog(xbai)                # dgdriv.f:161 XDGROW=ALOG(XDMULT), computed ONCE per species
        vardg = c.vardg[sp]
        # ssigma/rho `log` routed through the gfortran companion (flog) so SSIG matches FVS's ALOG bit-exactly
        # (dgf.f). dgscor's rejection `|FRM|>DGSD·SSIG` is DISCONTINUOUS in SSIG — a libm-log ULP could flip a
        # borderline accept/reject and desync the draw stream — so matching FVS's exact log is the FAITHFUL choice.
        # Proven SAFE + inert (2026-07-05): snt01 bit-exact + full suite unchanged after routing ⇒ the old memory
        # warning "routing ssigma/rho desyncs the RNG" was WRONG (openlibm==gfortran for these vardg ranges), and
        # it DECONFOUNDS the sp33/65 DGSCOR tail: NOT the ssigma log ⇒ it's the upstream WK2 past-dbh/WK3 calib. sqrt stays IEEE.
        evarp1 = (sqrt(1f0 + 4f0 * vardg * pvmlt) + 1f0) / 2f0
        sig1   = sqrt(flog(max(evarp1, 1f0 + eps(Float32))))
        evarp2 = (sqrt(1f0 + 4f0 * vardg * vmlt) + 1f0) / 2f0
        ssigma = sqrt(flog(max(evarp2, 1f0 + eps(Float32))))
        rho = (sig1 > 0f0 && ssigma > 0f0) ?
              flog(1f0 + corr * sqrt((evarp1 - 1f0) * (evarp2 - 1f0))) / (sig1 * ssigma) : 0f0
        rhocp = sqrt(max(1f0 - rho * rho, 0f0))
        frmbase = DG_FM * ssigma * rhocp
        fru = DG_FU * ssigma * rhocp           # upper-triple FRM factor (dgdriv.f:91)
        frl = DG_FL * ssigma * rhocp           # lower-triple FRM factor (dgdriv.f:89)
        for k in i1:i2
            i = ind1[k]
            bark = _cr_dg ? cr_bratio(sd, sp, t.dbh[i], _cr_imodty) :
                   _tt_dg ? tt_bratio(Int(sp), t.dbh[i]) :
                   _bc_dg ? bc_bratio(Int(sp)) :
                   _bm_dg ? bm_bratio(sd, Int(sp), t.dbh[i]) :
                   _ci_dg ? ci_bratio(sd, Int(sp), t.dbh[i]) :
                   _wc_dg ? wc_bratio(sd, Int(sp), t.dbh[i]) :
                   _pn_dg ? wc_bratio(sd, Int(sp), t.dbh[i]) :
                   (_ec_dg || _ca_dg) ? wc_bratio(sd, Int(sp), t.dbh[i]) :
                   _so_dg ? so_bratio(sd, Int(sp), t.dbh[i]) :
                   _on_dg ? on_bratio(Int(sp), t.dbh[i], t.height[i]) :   # dgdriv.f:201 BRATIO(ISPC,DBH,HT), original DBH
                   _ak_dg ? ak_bratio(Int(sp), t.dbh[i]) : bark_ratio(bark_a, bark_b, sp, t.dbh[i])
            d_ib = t.dbh[i] * bark
            # FVS bounds the 5-yr DG (DGBND, dgdriv.f:255-269) THEN scales to the cycle length
            # (gradd.f:79-90, DDS·(FINT/YR)) WITHOUT re-bounding. So DDS here is the 5-yr basis (BAIMULT
            # only); `bsc` bounds the 5-yr DG and then scales by sfint/5. FINT=5 ⇒ identity (no scale).
            # dgdriv.f:161+206: DDS=EXP(WK2 + ALOG(XDMULT)) — the BAIMULT enters in LOG-space BEFORE the exp,
            # NOT as a post-exp `EXP(WK2)*XDMULT` (mathematically equal but a different Float32 op — differs ~36%
            # of the time by 1 ULP for xbai=1.5). xbai=1 ⇒ xdgrow=flog(1)=0 ⇒ fexp(wk2+0)=fexp(wk2), bit-identical
            # to the old ·1.0 ⇒ every non-BAIMULT scenario is untouched (verified 0-diff).
            dds5 = fexp(wk2[i] + xdgrow)                    # YR-yr DDS (BAIMULT: EXP(WK2+ln XDMULT)); YR=5 SN / 10 NE
            # DG bound+scale applied inline (pillar-2: was a per-tree `bsc` closure). `_bsc(dg5)` local macro-
            # style: identical call `_bound_scale(dlo_v, dhi_v, sp, t.dbh[i], d_ib, dg5, sfint, size_cap, yr)`.
            # exp routed through the gfortran companion (fexp, doctrine #8): the tripled records carry
            # OLDRN forward, so a 1-ULP openlibm-vs-libm exp diff COMPOUNDS across cycles (the timeint
            # non-native-cycle volume tail). ssigma/rho log IS now routed (flog, see above) — the old
            # "don't route, it desyncs the RNG" note was DISPROVEN (snt01 bit-exact, suite unchanged); it's
            # inert (openlibm==gfortran for these vardg ranges) and faithful to FVS's ALOG.
            size_cap = s.control.sp_size_cap
            # CR (cr/dgdriv.f:223-225,256-258,270-272) caps EACH DG's spread: GDIF=DG−WKI, GLIM=WKI·0.33,
            # IF GDIF>GLIM DG=WKI+GLIM — WKI is the un-FRM'd central DG (dgdriv.f:213 SQRT(DSQ+DDS)−D). Eastern
            # dgdriv.f has NO GLIM (only DGBND); jl's _bound_scale is DGBND-style, so CR needs this extra cap
            # BEFORE the bound. Without it CR's upper tripled records over-grow (dense-stand BA ~10% high).
            crv = s.variant isa CentralRockies
            wkicr = crv ? (sqrt(d_ib * d_ib + dds5) - d_ib) : 0f0
            if do_trip
                rnpar = oldrn[i]                            # original residual (dgdriv.f:116)
                frmt = frmbase + corr * rnpar; oldrn[i] = frmt
                dgc = sqrt(d_ib * d_ib + dds5 * fexp(frmt)) - d_ib
                crv && (dgc - wkicr > wkicr * 0.33f0) && (dgc = wkicr * 1.33f0)
                t.diam_growth[i] = _bound_scale(dlo_v, dhi_v, sp, t.dbh[i], d_ib, dgc, sfint, size_cap, yr)
                ru = fru + corr * rnpar; rnU[i] = ru
                dgu = sqrt(d_ib * d_ib + dds5 * fexp(ru)) - d_ib
                crv && (dgu - wkicr > wkicr * 0.33f0) && (dgu = wkicr * 1.33f0)
                dgU[i] = _bound_scale(dlo_v, dhi_v, sp, t.dbh[i], d_ib, dgu, sfint, size_cap, yr)
                rl = frl + corr * rnpar; rnL[i] = rl
                dgl = sqrt(d_ib * d_ib + dds5 * fexp(rl)) - d_ib
                crv && (dgl - wkicr > wkicr * 0.33f0) && (dgl = wkicr * 1.33f0)
                dgL[i] = _bound_scale(dlo_v, dhi_v, sp, t.dbh[i], d_ib, dgl, sfint, size_cap, yr)
            else
                if tripling
                    frmt = frmbase + corr * oldrn[i]       # deterministic (dgdriv.f:117)
                    oldrn[i] = frmt
                    frm = fexp(frmt)
                else
                    frm = dgscor!(s.rng, oldrn, i, ssigma, rho, rhocp, wk2[i];
                                  dgsd = s.control.dg_stddev_bound)
                end
                dgc = sqrt(d_ib * d_ib + dds5 * frm) - d_ib
                crv && (dgc - wkicr > wkicr * 0.33f0) && (dgc = wkicr * 1.33f0)
                t.diam_growth[i] = _bound_scale(dlo_v, dhi_v, sp, t.dbh[i], d_ib, dgc, sfint, size_cap, yr)
            end
        end
    end

    # Record tripling is performed LATER (after height growth + mortality), because
    # FVS runs MORTS before TRIPLE (grincr.f) — VARMRT must distribute mortality over
    # the original ITRN records, not the tripled set. Return the per-tree upper/lower
    # DG + serial-correlation residual so `triple_records!` can build the records.
    # htgU/htgL/is_small carry the small-tree (REGENT) tripled-record overrides that
    # `small_tree_growth!` fills after height growth (large-tree records keep is_small
    # = false and inherit the central HTG via copy_tree!).
    htgU = do_trip ? zeros(Float32, nlive) : Float32[]
    htgL = do_trip ? zeros(Float32, nlive) : Float32[]
    is_small = do_trip ? falses(nlive) : BitVector()
    return do_trip ? (nlive = nlive, dgU = dgU, dgL = dgL, rnU = rnU, rnL = rnL,
                      htgU = htgU, htgL = htgL, is_small = is_small) : nothing
end

"""
    triple_records!(state, stash)

TRIPLE (triple.f): split each of the `stash.nlive` original live records into 3
weighted records — central 0.60, upper 0.25 (FU growth), lower 0.15 (FL growth) —
using the DGs stashed by `diameter_growth!`. Run AFTER mortality so the split is of
the post-mortality TPA (matching FVS's MORTS-before-TRIPLE order). Dead records are
pushed to the end. No-op when `stash === nothing`.
"""
function triple_records!(s::StandState, stash)
    stash === nothing && return s
    t = s.trees; nlive = stash.nlive
    dgU = stash.dgU; dgL = stash.dgL; rnU = stash.rnU; rnL = stash.rnL
    htgU = stash.htgU; htgL = stash.htgL; is_small = stash.is_small
    @inbounds for k in t.ndead:-1:1
        copy_tree!(t, 3 * nlive + k, nlive + k)
    end
    @inbounds for i in 1:nlive
        # FVS TRIPLE (triple.f:18) appends BOTH new records per parent contiguously:
        # itfn = ITRN+2i-1 (weight .25), itfn+1 = ITRN+2i (weight .15). This interleaved
        # physical layout (not all-uppers-then-all-lowers) is what TREDEL's swap-from-end
        # walks after a thin, so it must match the oracle's append order exactly.
        u = nlive + 2 * i - 1; l = nlive + 2 * i
        copy_tree!(t, u, i); copy_tree!(t, l, i)
        t.tpa[u] = t.tpa[i] * 0.25f0; t.diam_growth[u] = dgU[i]; t.old_random[u] = rnU[i]
        t.tpa[l] = t.tpa[i] * 0.15f0; t.diam_growth[l] = dgL[i]; t.old_random[l] = rnL[i]
        # the record's period mortality (MortPA) splits with the surviving TPA (0.60/0.25/0.15)
        t.mort_pa[u] = t.mort_pa[i] * 0.25f0; t.mort_pa[l] = t.mort_pa[i] * 0.15f0
        # small-tree records carry per-record height increments (REGENT random effect)
        if is_small[i]
            t.ht_growth[u] = htgU[i]; t.ht_growth[l] = htgL[i]
        end
        t.tpa[i] *= 0.60f0; t.mort_pa[i] *= 0.60f0
        # lineage keys: upper=3K, central=3K+1, lower=3K+2 → species-sort then visits
        # records in the oracle's LNKCHN order (upper,central,lower depth-first), so the
        # untripled DGSCOR consumes the BACHLO stream in the same sequence (RNG-exact).
        kk = t.sort_key[i]
        t.sort_key[u] = 3 * kk; t.sort_key[i] = 3 * kk + 1; t.sort_key[l] = 3 * kk + 2
    end
    t.n = 3 * nlive
    return s
end
