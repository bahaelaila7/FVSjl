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
    buf = build_organon_buffer!(s)
    buf.runs || return nothing                      # no big-6 ⇒ FVS-native (oc/dgf.f), nothing here
    si_1, si_2 = _oc_organon_si(s)
    isp_fvs = Int[Int(t.species[i]) for i in 1:t.n]
    g = organon_execute_swo(buf, isp_fvs; si_1=si_1, si_2=si_2,
                            msdi_1=msdi, msdi_2=msdi, msdi_3=msdi, cyclg=cyclg)
    fscale = fint/5f0
    @inbounds for i in 1:t.n
        # mortality applies to every record ORGANON grew (valid + surrogate), oc/morts.f:498-504
        dead = g.deadexp[i]*fscale
        dead > t.tpa[i] && (dead = t.tpa[i])
        t.mort_pa[i] = dead
        t.tpa[i] -= dead
        buf.iorg[i] == 1 || continue                # only valid ORGANON trees get DG/HTG/CR from ORGANON
        dg = oc_organon_dg(isp_fvs[i], t.dbh[i], g.dds[i])
        bark = oc_bratio(isp_fvs[i], t.dbh[i])
        t.diam_growth[i] = dg
        t.dbh[i]    += dg/bark                       # outside-bark DBH growth (update.f)
        t.ht_growth[i] = g.hgro[i]
        t.height[i] += g.hgro[i]
        crnew = round(g.cr2[i]*100f0, RoundNearestTiesAway)   # ANINT (oc/crown.f:282)
        t.crown_pct[i]   = Int32(crnew)
        t.crown_ratio[i] = g.cr2[i]
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
    organon_apply_growth!(s; fint=sfint)
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
