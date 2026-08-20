# =============================================================================
# organon_hook.jl — OC (Oregon Coast) ORGANON live growth hook + StandState copy-back (chunk C7).
#
# Wires the validated ORGANON EXECUTE/GROW orchestration (`organon_execute_swo`, C7 sub-step 1)
# into the shared per-cycle engine, and copies its outputs back into the `StandState` tree records
# at their FVS sites:
#   • DG  → t.diam_growth  = oc_organon_dg  (oc/dgdriv.f:536-557; the shared apply-loop then grows
#                                            the outside-bark DBH by DG/oc_bratio — DBH += ~DGRO)
#   • HTG → t.ht_growth    = HGRO           (oc/htgf.f:96, SCALE·XHT·HGRO·EXP(HTCON) = HGRO defaults)
#   • CR  → t.crown_pct    = ANINT(CR2·100) (oc/crown.f:282)
#   • MORT→ t.tpa         −= MORTEXP·(FINT/5) (oc/morts.f:499)
#
# STATUS (C7 sub-step 2 — StandState apply): `organon_apply_growth!` runs the whole ORGANON growth on
# a `StandState` and applies all four copy-backs for the valid ORGANON trees (IORG=1), bit-exact vs
# FVSoc_clean. `diameter_growth!(::OregonCoast)` now drives it instead of erroring. NON-ORGANON /
# no-big-6 trees (IORG=0) grow FVS-native `oc/dgf.f` (Wykoff, CA-family) — a SEPARATE port still
# pending; until it lands `organon_apply_growth!` leaves those records' DG/HTG at 0 (flagged), so a
# FULL-stand `.sum` is not yet bit-exact — only the ORGANON-tree state is.
#
# MSDI sourcing: the ORGANON max-SDI RVARS(3..5) (=815 on ocmin) is a per-stand ORGANON input; until
# the SITSET path that fills it is ported, `organon_apply_growth!` takes `msdi` (default 0 ⇒ SUBMAX's
# built-in SWO fallback A1). For ocmin pass msdi=815 (the measured RVARS).
# =============================================================================

"organon/execute2.f SI_1/SI_2 = SITE_1−4.5 / SITE_2−4.5 for OC (DF site sp 7, PP site sp 18)."
@inline function _oc_organon_si(s::StandState)
    site1 = s.plot.sp_site_index[7]
    site2 = s.plot.sp_site_index[18]
    # organon/execute2.f:262-267 VERSION=1 SI conversion (mirrors site_setup!/C2)
    if site1 <= 0f0 && site2 > 0f0
        site1 = 1.062934f0*site2
    elseif site2 <= 0f0
        site2 = 0.940792f0*site1
    end
    return site1 - 4.5f0, site2 - 4.5f0
end

"""
    organon_apply_growth!(s; msdi=0f0, cyclg=0) -> OrganonGrowth

Run the ORGANON SWO per-cycle growth on the `StandState` and copy the results into the tree records
for the valid ORGANON trees (IORG=1): grow DBH by the bark-consistent `DG/oc_bratio`, HT by `HGRO`,
set crown `ANINT(CR2·100)`, and reduce TPA by `MORTEXP·(FINT/5)`. Returns the `OrganonGrowth` for
inspection. Deterministic (DGSD=0). Non-ORGANON records are left unchanged (see file header).
"""
function organon_apply_growth!(s::StandState; msdi::Float32 = 0f0, cyclg::Int = 0,
                               fint::Float32 = 5f0)
    t = s.trees
    # oc/dgdriv.f: DGF runs for EVERY tree (WK2 = FVS-native ln DDS), on the ORIGINAL DBH, BEFORE any
    # growth is applied. Then the ORGANON EXECUTE overwrites WK2 for the valid ORGANON trees (IORG=1).
    oc_dgcons!(s)                                    # per-species DGCON (site constants)
    dgf!(s, s.variant)                               # WK2 = FVS-native ln(DDS) for all trees
    wk2 = view(s.scratch.wk, 2, :)
    buf = build_organon_buffer!(s)                   # IORG gate + ORGANON input buffer (original DBH)
    si_1, si_2 = _oc_organon_si(s)
    isp_fvs = Int[Int(t.species[i]) for i in 1:t.n]
    # ACALIB(1,1..18) from setup PREPARE (oc_organon_prepare!) — HTGRO2's minor-species height
    # calibration (organon/htgrowth.f:235-236). Inert (1.0) on FVS/FIA inventory unless a minor
    # ORGANON species had ≥2 measured-height trees; HTGRO1 (big-6) ignores it (htgrowth.f:56).
    calib1 = Float32[s.calib.organon_acalib[1, g] for g in 1:18]
    # ORGANON growth (only when a big-6 tree exists); else the whole stand is FVS-native.
    g = buf.runs ? organon_execute_swo(buf, isp_fvs; si_1=si_1, si_2=si_2,
                       msdi_1=msdi, msdi_2=msdi, msdi_3=msdi, cyclg=cyclg, calib1=calib1) : nothing
    fscale = fint/5f0
    @inbounds for i in 1:t.n
        d0 = t.dbh[i]; d0 <= 0f0 && continue
        sp = isp_fvs[i]
        iorg = buf.iorg[i] == 1
        # DIAMETER: DDS from ORGANON (IORG=1) or the FVS-native DGF (IORG=0); both → DG via the shared
        # sqrt path, DBH grows outside-bark by DG/BARK (oc/dgdriv.f:536-557, update.f).
        dds = (iorg && g !== nothing) ? g.dds[i] : wk2[i]
        dg = oc_organon_dg(sp, d0, dds)
        bark = oc_bratio(sp, d0)
        t.vol_bark[i] = bark             # BRATIO(D_start) for CFTOPK/BFTOPK (vols.f:150); the shared
                                         # apply-loop skips OC so this pre-growth value survives.
        t.diam_growth[i] = dg
        t.dbh[i] = d0 + dg/bark
        # HEIGHT + CROWN: ORGANON HGRO/CR2 for IORG=1; FVS-native HTGF for IORG=0 (oc/htgf.f).
        if iorg && g !== nothing
            t.ht_growth[i] = g.hgro[i]
            t.height[i] += g.hgro[i]
            t.crown_pct[i] = Int32(round(g.cr2[i]*100f0, RoundNearestTiesAway))   # ANINT (oc/crown.f:282)
        else
            pt = Int(t.plot_id[i])
            pccf = (1 <= pt <= length(s.density.point_ccf)) ? s.density.point_ccf[pt] : 0f0
            si = s.plot.sp_site_index[sp]
            h_old = t.height[i]
            htg_large = oc_htgf_native(sp, h_old, t.crown_pct[i], s.plot.avg_height, pccf, si)
            # REGENT: trees with DBH < XMAX(sp) use the small-tree height-age model (SMHTGF), not the
            # large-tree HTGF (which drops RELHT suppression when PCCF<100). DGSD=0 ⇒ deterministic.
            # HCOR (small-tree height calibration) is 0 here; the ratio-estimator calibration is a
            # follow-on (only affects species with ≥5 measured-height small trees). See small_tree_growth.jl.
            is_regen = 1 <= sp <= 50 && d0 < OC_REG_XMAX[sp]
            htg = is_regen ?
                  oc_regent_htg(sp, d0, h_old, t.crown_pct[i], t.crown_ratio[i],
                                s.plot.basal_area, s.plot.avg_height, si, 0f0, htg_large) :
                  htg_large
            t.ht_growth[i] = htg
            t.height[i] += htg
            # REGENT DBH-from-height: trees with DBH < DGMIN(sp) derive DBH growth from the H-D function
            # (regent.f:305; DGMIN≤D<XMAX keeps the large-tree DBH). DBH += DGSM directly (measured live
            # 9987 DBH-INC = DGSM); DK/DKK from HTDBH (Curtis/Arney, LHTDRG=.FALSE. default). Large-tree
            # part of the XWT blend = dg/bark. SCALE2 = YR(5)/FINT (=1 for a 5-yr cycle).
            if is_regen && d0 < OC_REG_DGMIN[sp]
                newd = oc_regent_dbh(sp, d0, h_old, htg, bark, dg/bark, 5f0/fint)
                t.diam_growth[i] = newd - d0
                t.dbh[i] = newd
            end
        end
        # Broken/dead-top trees: grow the NORMAL (NORMHT) height by the same increment as the standing
        # height (update.f:65-67 `NORMHT=INT(REAL(NORMHT)+(HTG*100.+.5))`; op order matched exactly).
        t.norm_ht[i] > 0 &&
            (t.norm_ht[i] = trunc(Int32, Float32(t.norm_ht[i]) + (t.ht_growth[i]*100f0 + 0.5f0)))
        # MORTALITY: ORGANON MORTEXP for every record it grew (oc/morts.f:498-504).
        if g !== nothing
            dead = g.deadexp[i]*fscale
            dead > t.tpa[i] && (dead = t.tpa[i])
            t.mort_pa[i] = dead
            t.tpa[i] -= dead
        end
    end
    return g
end

# --- Live growth hook (oc/dgdriv.f entry) — grow_cycle! integration ----------------------------
# C7: OregonCoast drives the ORGANON engine. Because ORGANON computes DIAMETER, HEIGHT, CROWN and
# MORTALITY together in one EXECUTE, `diameter_growth!(::OregonCoast)` is the SINGLE authority: it
# runs the whole growth and applies all four copy-backs to the StandState HERE, then ZEROS
# `diam_growth`/`ht_growth` so the shared engine's later apply-loop (`DBH += DG/bark`,
# `HT += HTG`) is a no-op — and the OC `height_growth!`/`mortality!`/`crown_ratio_update!`/
# `small_tree_growth!` hooks are no-ops (the work is already done). This avoids the double-apply the
# cooperating-hook split would risk, and is faithful: FVS's `oc/dgdriv.f`→EXECUTE likewise produces
# DG/HTG/CR/MORTEXP in one call, which FVS then copies at dgf/htgf/crown/morts. Returns `nothing`
# (no tripling: INDS(5)=0, DGSD=0 on OC).
function diameter_growth!(s::StandState, ::OregonCoast; tripling::Bool = false,
                          sfint::Float32 = 5f0, kwargs...)
    # MSDI_1/2/3 = SDIDEF(7/18/4) (oc/sitset.f:340-342); inert on ocmin (RD≤RDCC base mortality)
    # but faithful for dense stands. Sourced from the ecoclass-filled sp_sdi_def (site_setup!).
    msdi = (length(s.plot.sp_sdi_def) >= 7 && s.plot.sp_sdi_def[7] > 0f0) ? s.plot.sp_sdi_def[7] : 0f0
    organon_apply_growth!(s; fint=sfint, msdi=msdi)
    t = s.trees
    @inbounds for i in 1:t.n
        t.diam_growth[i] = 0f0    # DBH/HT already grown in organon_apply_growth!; zero so the shared
        t.ht_growth[i]   = 0f0    # grow_cycle! apply-loop (DBH+=DG/bark, HT+=HTG) is inert for OC.
    end
    return nothing
end

# ORGANON did height/mortality/crown inside `diameter_growth!` → these shared hooks are no-ops for OC.
height_growth!(s::StandState, ::OregonCoast; kwargs...) = s
small_tree_growth!(s::StandState, stash, ::OregonCoast; kwargs...) = s
mortality!(s::StandState, ::OregonCoast; kwargs...) = s
crown_ratio_update!(s::StandState, ::OregonCoast; kwargs...) = s
