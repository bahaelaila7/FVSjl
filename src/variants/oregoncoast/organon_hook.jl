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

"""
    organon_thin_update!(th, ba, cyclg) -> th

oc|op/dgdriv.f:304-327 "MANAGE THE THINNING VARIABLES": with OLDBA (grincr.f:283, pre-CUTS BA, th[12]) above the current
(post-thin) BA, shift YST/BART down one slot and set YST(1)=(CYCLG-1)*5, BART(1)=OLDBA-ATBA (ATBA=BA, grincr.f:319),
BABT=OLDBA. Once per cycle (th[13]).
"""
function organon_thin_update!(th::Vector{Float32}, ba::Float32, cyclg::Int)
    if th[13] != Float32(cyclg) + 1f0 && th[12] > ba
        @inbounds for i in 5:-1:2
            th[i] = th[i-1]; th[5+i] = th[5+i-1]
        end
        th[1] = Float32((cyclg - 1) * 5)
        th[6] = th[12] - ba
        th[11] = th[12]
    end
    th[13] = Float32(cyclg) + 1f0
    return th
end

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
                               fint::Float32 = 5f0, msdi2::Float32 = msdi, msdi3::Float32 = msdi)
    t = s.trees
    # oc/dgdriv.f:494-513 COR/HCOR attenuation toward the LSTART calibration goal WCI (shared with OP). DGF below reads
    # COR (CONSPP=DGCON+COR); REGENT reads HCOR. Inert (0) when no species calibrated.
    _organon_cor_attenuate!(s, Float32(current_cycle_year(s) - Int(s.control.cycle_year[1])), fint)
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
    # execute2.f:231-252 CALIB = the cycle-decayed ACALIB (CALH/CALD on): row 1 → HTGRO2, row 3 → DIAMGRO_RUN.
    th = organon_thin_update!(s.calib.organon_thin, s.plot.basal_area, cyclg)
    cal = organon_calib_schedule(s.calib.organon_acalib, cyclg)
    calib1 = Float32[cal[1, g] for g in 1:18]
    calib3 = Float32[cal[3, g] for g in 1:18]
    # ORGANON growth (only when a big-6 tree exists); else the whole stand is FVS-native.
    g = buf.runs ? organon_execute_swo(buf, isp_fvs; si_1=si_1, si_2=si_2,
                       msdi_1=msdi, msdi_2=msdi2, msdi_3=msdi3, cyclg=cyclg, calib1=calib1, calib3=calib3,
                       stor=s.calib.organon_stor, babt=th[11], bart=th[6:10], yt=th[1:5]) : nothing
    fscale = fint/5f0
    yrscale = fint/5f0                               # SCALE=FINT/YR (htgf.f:63, regent.f:197), YR=5 for the ORGANON variants
    # oc/dgdriv.f:463 CALL CLGMULT(WK4) ⇒ :536 DDS=EXP(WK2+XDGROW)·WK4 for EVERY tree (ORGANON-folded or native).
    cw = climate_growth_wk4!(s, Float32(current_cycle_year(s)) + fint / 2f0)
    _misdg = _dm_effects_on(s)                       # oc/dgdriv.f:562 MISDGF (misintoc.f DGPDMR)
    _, _mdgp, _, _mmaxsp = _mis_tables(s.variant)
    sizcap = s.control.sp_size_cap
    dbh0 = s.calib.oc_dbh0; resize!(dbh0, t.n)
    ht0 = s.calib.oc_ht0; cr0 = s.calib.oc_cr0; resize!(ht0, t.n); resize!(cr0, t.n)
    @inbounds for i in 1:t.n; ht0[i] = reinterpret(Int32, t.height[i]); cr0[i] = t.crown_pct[i]; end   # FFE reads PRE-growth HT/ICR
    @inbounds for i in 1:t.n; dbh0[i] = t.dbh[i]; end
    @inbounds for i in 1:t.n
        d0 = t.dbh[i]; d0 <= 0f0 && continue
        sp = isp_fvs[i]
        iorg = buf.iorg[i] == 1
        # DIAMETER: DDS from ORGANON (IORG=1) or the FVS-native DGF (IORG=0); both → DG via the shared
        # sqrt path, DBH grows outside-bark by DG/BARK (oc/dgdriv.f:536-557, update.f).
        dds = (iorg && g !== nothing) ? g.dds[i] : wk2[i]
        # oc/dgdriv.f:536-570: DDS=EXP(WK2+XDGROW)·WK4 (XDGROW=ALOG(XDMULT), the BAIMULT keyword), DG, ×MISDGF (dwarf
        # mistletoe, START-of-cycle DMR), then DGBND (oc/dgbnd.f: the SIZCAP cap only, no DGMAX envelope).
        xdgrow = flog(active_multiplier(s.control, :bai, sp, current_cycle_year(s)))
        dg = oc_organon_dg(sp, d0, dds + xdgrow, cw === nothing ? 1f0 : cw[i])
        _misdg && (dg *= ie_dm_dg_mult(_mdgp, _mmaxsp, sp, Int(t.dmr[i])))
        if (d0 + dg) > sizcap[sp, 1] && sizcap[sp, 3] < 1.5f0
            dg = sizcap[sp, 1] - d0
            dg < 0.01f0 && (dg = 0.01f0)
        end
        bark = oc_bratio(sp, d0)
        dbh_direct = nothing
        t.vol_bark[i] = bark             # BRATIO(D_start) for CFTOPK/BFTOPK (vols.f:150); the shared
                                         # apply-loop skips OC so this pre-growth value survives.
        # HEIGHT (oc/htgf.f): ORGANON HGRO for IORG=1 (htgf.f:93-100 SCALE·XHT·HGRO·EXP(HTCON)); FVS-native HTGF
        # otherwise; then ×MISHGF for every record (htgf.f:266). CROWN: ANINT(CR2·100) for IORG=1 (oc/crown.f:282).
        mhg = mis_hg_mult(s, i)
        h_old = t.height[i]
        crnew_org = 0f0                  # ANINT(CR2·100) of a valid-ORGANON tree, finished after the height is known (below)
        if iorg && g !== nothing
            htg = (yrscale * g.hgro[i]) * mhg           # htgf.f:93 HTG=SCALE*XHT*HGRO*EXP(HTCON), SCALE=FINT/YR; then ×MISHGF
            crnew_org = round(g.cr2[i]*100f0, RoundNearestTiesAway)
        else
            pt = Int(t.plot_id[i])
            pccf = (1 <= pt <= length(s.density.point_ccf)) ? s.density.point_ccf[pt] : 0f0
            si = s.plot.sp_site_index[sp]
            htg = (yrscale * oc_htgf_native(sp, h_old, t.crown_pct[i], s.plot.avg_height, pccf, si)) * mhg   # htgf.f:244 SCALE·HTG
        end
        # REGENT (oc/regent.f:157-370) for EVERY record with D < XMAX(sp), valid ORGANON trees included: the small-tree
        # height (SMHTGF·CON, blended by XWT; XWT=1 for IORG=1 so the ORGANON HTG stands), the 0.1-ft floor and SIZCAP(4);
        # then, for D < DGMIN(sp), the HTDBH-curve diameter step (oc_regent_dg).
        if 1 <= sp <= 50 && d0 < OC_REG_XMAX[sp]
            if !(iorg && g !== nothing)
                si = s.plot.sp_site_index[sp]
                htg = oc_regent_htg(sp, d0, h_old, t.crown_pct[i], t.crown_ratio[i],
                                    s.plot.basal_area, s.plot.avg_height, si, s.calib.htg_cor_small[sp], htg, yrscale)
            end
            htg < 0.1f0 && (htg = 0.1f0)
            if (h_old + htg) > sizcap[sp, 4]
                htg = sizcap[sp, 4] - h_old
                htg < 0.1f0 && (htg = 0.1f0)
            end
            if d0 < OC_REG_DGMIN[sp]
                dg, dbh_direct = oc_regent_dg(sp, d0, h_old, htg, bark, dg, 5f0/fint, iorg && g !== nothing, sizcap)
            end
        end
        t.ht_growth[i] = htg
        t.height[i] = h_old + htg
        if iorg && g !== nothing
            # oc/crown.f:281-288,368-380,426-433 for an IORG=1 tree: CRNEW=ANINT(CR2·100), CHG=CRNEW-ICR, CRNEW=ICR+CHG, ICRI=INT(CRNEW+.5), then
            # (ICR>0) the CRMAX cap CRMAX=(HT·ICR/100+HTG)/(HT+HTG)·100 on the UPDATEd HT, ICRI>CRMAX ⇒ INT(CRMAX+.5) — ORGANON trees included —
            # then ICRI ≤ 95 and ≥ 1 (LORGANON=.TRUE. ⇒ no <10 bump). MEASURED FVSoc_tiered 645182389126144 econ MA 2033: ICR 1,
            # CR2 0.02507 ⇒ ICRI 3, HT 7.954158, HTG 0.12111 ⇒ CRMAX 2.4846 ⇒ ICRI 2 live (jl wrote 3).
            icr_old = Int(t.crown_pct[i])
            icri = trunc(Int, (Float32(icr_old) + (crnew_org - Float32(icr_old))) + 0.5f0)
            if icr_old != 0
                crln = t.height[i] * Float32(icr_old) / 100f0
                crmax = (crln + htg) / (t.height[i] + htg) * 100f0
                Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
            end
            icri > 95 && (icri = 95); icri < 1 && (icri = 1)
            t.crown_pct[i] = Int32(icri)
        end
        t.diam_growth[i] = dg            # YR-basis DG: MORTS reads it so (morts.f:210 G=(DG/BARK)*(FINT/5)); gradd_dg_scale! rescales it after MORTS
        # gradd.f:67-90 (FINT≠YR): DG → FINT basis, BARK=BRATIO(D,HT) at the pre-UPDATE D — the value UPDATE then adds to DBH.
        dg_upd = dg
        if fint != 5f0
            if dg > 0f0
                ddsf = (dg * (2f0 * bark * d0 + dg)) * yrscale
                dg_upd = sqrt((d0 * bark)^2 + ddsf) - bark * d0
            else
                dg_upd = 0f0
            end
        end
        t.dbh[i] = dbh_direct === nothing ? d0 + dg_upd/bark : dbh_direct     # update.f:115 DBH=DBH+DG/BRATIO
        # REGENT (oc/regent.f:285-345) writes DBH(K) itself for a tree with D < DGMIN BEFORE MORTS/FMMAIN run (the rest of the growth
        # is applied at UPDATE), so the pre-UPDATE DBH those passes read is that value, not the start-of-cycle one: MEASURED FVSoc_tiered
        # 15148784010497 plant_cyc cycle 3 MORTS ES014 D=0.10667 (REGENT) vs start 0.10261 ⇒ TOKILL 7.696978 vs 7.696754.
        dbh_direct === nothing || (dbh0[i] = dbh_direct)
        # Broken/dead-top trees: grow the NORMAL (NORMHT) height by the same increment as the standing
        # height (update.f:65-67 `NORMHT=INT(REAL(NORMHT)+(HTG*100.+.5))`; op order matched exactly).
        t.norm_ht[i] > 0 &&
            (t.norm_ht[i] = trunc(Int32, Float32(t.norm_ht[i]) + (t.ht_growth[i]*100f0 + 0.5f0)))
    end
    # MORTALITY: stash ORGANON's MORTEXP; oc/morts.f:498-504 WK2(I)=MORTEXP(I)*(FINT/5.) for EVERY record when ORGANON
    # ran (SMORMT>0). mortality!(::OregonCoast) applies it in the MORTS slot (PROB is untouched until UPDATE in FVS).
    c0 = s.calib
    length(c0.op_mortexp) != t.n && resize!(c0.op_mortexp, t.n)
    @inbounds for i in 1:t.n; c0.op_mortexp[i] = g !== nothing ? g.deadexp[i]*fscale : 0f0; end
    c0.op_org_ran = g !== nothing
    # Stash IORG for the GRADD-CROWN pass (crown_ratio_update!(::OregonCoast)): the FVS-native rank-
    # Weibull crown-change model must run for the IORG=0 trees; the valid-ORGANON (IORG=1) crowns are
    # already set above from CR2. op_iorg is the shared ORGANON per-tree flag (also used by Olympic).
    c = s.calib
    # IORG is a SLOT-indexed COMMON array (oc/dgdriv.f:231-256 sets IORG(1..ITRN) each call and never clears the
    # rest): slots beyond this call's ITRN keep their last value, which CROWN reads for trees ESTAB adds after the
    # call (oc/crown.f:281 IF(LORGANON .AND. IORG(I).EQ.1) with CR2(I)=0 ⇒ ICR=1). Grow-only stash mirrors that.
    length(c.op_iorg) < t.n && resize!(c.op_iorg, t.n)
    @inbounds for i in 1:t.n; c.op_iorg[i] = buf.iorg[i]; end
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
    # inert on ocmin (RD≤RDCC base mortality) but faithful for dense stands; sp_sdi_def from site_setup! (sitset).
    # organon/execute2.f:138-140 MSDI_1/2/3 = RVARS(3/4/5), which oc/sitset.f:327-329 loads with SDIDEF(7) (DF),
    # SDIDEF(4) (WF/GF), SDIDEF(18) (PP) — each its own value (they differ when the ecoclass gives species SDImax).
    sd = length(s.calib.organon_msdi) == 3 ? s.calib.organon_msdi : Float32[s.plot.sp_sdi_def[7], s.plot.sp_sdi_def[4], s.plot.sp_sdi_def[18]]
    # oc/dgdriv.f:177 CYCLG = ICYC-1 = the 0-based jl cycle (EXECUTE's CALIB decay, DG/HG thin + fert clocks).
    organon_apply_growth!(s; fint=sfint, msdi=sd[1], msdi2=sd[2], msdi3=sd[3], cyclg=Int(s.control.cycle))
    # DBH/HT/NORMHT already grown INLINE in organon_apply_growth!; the shared grow_cycle! apply-loop SKIPS OC
    # (simulate.jl) so it does NOT double-apply. We KEEP diam_growth/ht_growth (the applied increment) so the
    # FVS_TreeList DG/HtG columns are populated for OC (they were 0 on every projected cycle when zeroed here).
    return nothing
end

# ORGANON did height/mortality/crown inside `diameter_growth!` → these shared hooks are no-ops for OC.
height_growth!(s::StandState, ::OregonCoast; kwargs...) = s
small_tree_growth!(s::StandState, stash, ::OregonCoast; kwargs...) = s
# crown_ratio_update!(::OregonCoast) is the FVS-native rank-Weibull crown CHANGE for IORG=0 trees —
# see variants/oregoncoast/crown.jl (was a no-op ⇒ IORG=0 crowns frozen ⇒ LP height under-grew).

# MORTALITY snag-booking: ORGANON applied its MORTEXP inline in `diameter_growth!`
# (`organon_apply_growth!`), recording the per-record killed density in `t.mort_pa`. The DBH/HT/CR/TPA
# work is already done, so this hook does NOT re-apply mortality — but the ongoing (non-fire) ORGANON
# mortality still has to reach the FFE standing-snag pool (fmkill.f ITYP=4 → fmsadd/fmsdit), exactly as
# every Wykoff variant's `mortality!` books its `killed` density. Without this, OC booked snags only from
# input-dead + fire, so `snag_summary` decayed to ~0.1× the oracle in the tail and missed the first
# cycle's mortality pulse (ocsnag 1993 14.76 vs oracle 55.26). `book_mortality_snags!` is a no-op when the
# FFE is inactive or the SN biomass coeffs are absent, and guards on `dbh>0` so stale `mort_pa` on empty
# records is skipped. `book_snags=false` on the fire cycle: the fire path owns snag booking there.
function mortality!(s::StandState, ::OregonCoast; fint::Float32 = 5f0, book_snags::Bool = true, kwargs...)
    t = s.trees; n = t.n; c = s.calib
    n == 0 && return _clim_mort_empty!(s, fint)   # ITRN=0 ⇒ oc/morts.f still reaches CLMORTS
    # oc/morts.f:463-466 SMORMT=ΣMORTEXP: the ORGANON rates are used only when EXECUTE ran AND produced some kill.
    smormt = (c.op_org_ran && length(c.op_mortexp) == n) ? sum(view(c.op_mortexp, 1:n)) : 0f0
    if !(smormt > 0f0)
        # No big-6 tree ⇒ EXECUTE never ran (SMORMT=0): oc/morts.f takes the FVS-native SDI/background path, which is
        # ca/morts.f line for line except FINT/5 (YR=5) for FINT/10 — so it runs the shared driver the CA variant
        # uses (PMSC/PMD via IBGMAP, RI=0.5·RI, VARMRT == ca/varmrt.f, MISMRT/CLMORTS/FIXMORT/snags inside). MORTS
        # reads the pre-UPDATE DBH; the OC DG driver grew DBH inline, so run it on the stashed start-of-cycle DBH.
        d0 = c.oc_dbh0
        dend = length(d0) == n ? Float32[t.dbh[i] for i in 1:n] : Float32[]
        isempty(dend) || @inbounds(for i in 1:n; t.dbh[i] = d0[i]; end)
        invoke(mortality!, Tuple{StandState,AbstractVariant}, s, s.variant; fint = fint, book_snags = book_snags)
        isempty(dend) || @inbounds(for i in 1:n; t.dbh[i] = dend[i]; end)
        @inbounds for i in 1:n; t.mort_pa[i] = s.scratch.mort_killed[i]; end
        return s
    end
    length(s.scratch.mort_killed) < n && resize!(s.scratch.mort_killed, n)
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    # oc/morts.f:499 sets WK2(I)=MORTEXP(I)*(FINT/5.) UNclamped (only the local WKI is capped at PROB, :500) and nothing later
    # clamps it unless SIZCAP applies (:713-714); update.f:72-73 caps the kill, dbstrls.f:158 prints WK2 itself. `wk2_raw`
    # keeps that value for FVS_TreeList MortPA via s.scratch.mort_killed (a record with 5-yr MORTEXP·FINT/5 > PROB reports MortPA > TPA).
    wk2_raw = fill(-1f0, n)
    # oc/morts.f:493-504: ORGANON ran ⇒ WK2(I)=MORTEXP(I)·(FINT/5) for every record (the per-record ORGANON kill
    # MORTEXP stashed by the DG driver, fscale already applied). [No-big-6 stands take the FVS-native SDI path, :506-546
    # + VARMRT — not yet ported for OC.]
    if smormt > 0f0
        @inbounds for i in 1:n
            p = t.tpa[i]; p <= 0f0 && continue
            wki = c.op_mortexp[i]
            wki > p && (wki = p)
            killed[i] = wki
            wk2_raw[i] = c.op_mortexp[i]
        end
    end
    # SIZCAP size-limit mortality (oc/morts.f:704-717, reached by the ORGANON path too — it precedes `IF(SMORMT.GT.0.)GO TO 45`): a tree whose
    # end-of-cycle DBH D+G, G=(DG/BARK)*(FINT/5) on the pre-UPDATE D, reaches SIZCAP(1) loses at least P*SIZCAP(2)*FINT/5. MEASURED FVSoc_tiered
    # 645182389126144 TIMEINT 10: DF D21.7 → 23.22 lost 1.4502 of 6.857 TPA live (MORTEXP·2 = 0.031).
    d0v = c.oc_dbh0
    if length(d0v) == n
        sizcap = s.control.sp_size_cap
        @inbounds for i in 1:n
            sp = Int(t.species[i]); (sp < 1 || sp > size(sizcap, 1)) && continue
            p = t.tpa[i]; d = d0v[i]
            bark = oc_bratio(sp, d)
            g = (t.diam_growth[i] / bark) * (fint / 5f0)
            if (d + g) >= sizcap[sp, 1] && trunc(Int, sizcap[sp, 3]) != 1
                kc = (p * sizcap[sp, 2] * fint) / 5f0
                killed[i] < kc && (killed[i] = kc)
                killed[i] > p && (killed[i] = p)
                wk2_raw[i] = -1f0                                      # morts.f:714 WK2 capped at P
            end
        end
    end
    # Climate-FVS mortality (oc/morts.f:801 CALL CLMORTS, after the base kill WK2, before FIXMORT): CLMORTS
    # MAX-combines each record's climate rate into WK2 against the PRE-mortality PROB.
    (s.climate !== nothing && s.climate.active) &&
        apply_climate_mort!(s, killed, Float32(current_cycle_year(s)) + fint / 2f0, fint)
    # MISMRT (gradd.f:96 MISTOE → mistoe.f:522 → mismrt.f:185-191, misintoc.f APMC) MAX-combines the dwarf-mistletoe
    # kill WK2=MAX(WK2,PROB·rate) on the pre-UPDATE DBH (the OC DG driver grew DBH inline; OC never triples, ICL4=0).
    if _ie_mis_variant(s.variant) && _dm_effects_on(s) && any(i -> t.dmr[i] > 0, 1:n)
        d0 = c.oc_dbh0
        dg_end = (length(d0) == n) ? Float32[t.dbh[i] for i in 1:n] : Float32[]
        isempty(dg_end) || @inbounds(for i in 1:n; t.dbh[i] = d0[i]; end)
        ie_dm_mortality_combine!(killed, s, fint, n)
        isempty(dg_end) || @inbounds(for i in 1:n; t.dbh[i] = dg_end[i]; end)
    end
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n
        t.tpa[i] = max(0f0, t.tpa[i] - killed[i])
        t.mort_pa[i] = killed[i]
        wk2_raw[i] > killed[i] && (killed[i] = wk2_raw[i])   # grow_cycle!'s _morts_wk2 reads this buffer for MortPA/OMORT
    end
    return s
end
