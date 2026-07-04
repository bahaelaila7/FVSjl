# =============================================================================
# fire/snag.jl — snag falldown + decay dynamics (FFE chunk F7-state core)
#
# Ported from: bin/FVSsn_buildDir/fmsfall.f (FMSFALL snag falldown) + the DECAYX
# hard→soft transition (fmvinit.f / fmsnag.f).
#
# When the fire (or ordinary mortality) kills a tree it becomes a standing snag. Each
# year a fraction of the snags fall (transferring to coarse woody debris) and the
# standing-hard snags decay toward soft. These are the per-snag-record rates that drive
# the stateful snag list (the records + per-cycle loop are F7-state's container, which
# builds on these functions). Rates come from the species' snag class (`snag_fallx`,
# `snag_alldwn`, `snag_decayx` in `fire_species_props.csv`).
# =============================================================================

"""
    snag_fall_density(coef, ksp, d, origden, denttl) -> Float32

Density of snags (stems/acre) that fall in one year for a snag record of species `ksp`
and DBH `d` (FMSFALL, fmsfall.f). `origden` is the record's original density, `denttl`
the density still standing. Small snags (< 12" and not redcedar) fall at a linear
`MODRATE·origden`; large snags use the last-5% logic that ramps the final stems down to
zero by the species' `ALLDWN` year.
"""
function snag_fall_density(coef::SpeciesCoefficients, ksp::Integer, d::Float32,
                           origden::Float32, denttl::Float32;
                           fallx::Float32 = coef_col(coef, :snag_fallx)[ksp],
                           alldwn::Float32 = coef_col(coef, :snag_alldwn)[ksp],
                           variant = nothing)::Float32
    # BASE fall rate (fmsfall.f:128/130) is VARIANT-SPECIFIC: SN/CS use −0.001679·d+0.064311; LS uses the
    # "new equation" −0.006·d+0.18 (a much faster fall); NE uses an ALGSLP table (not yet ported — NE keeps
    # the SN form here). The small-snag LINEAR-fall breakpoint also differs: SN/CS = 12" (redcedar ksp2 keeps
    # the last-5% ramp); LS = 18", except cedar/tamarack (ksp 10,11,14) at 12" (fmsfall.f LS:139-145).
    if variant isa LakeStates
        modrate = clamp((-0.006f0 * d + 0.18f0) * fallx, 0.01f0, 1f0)   # FVS clamps MODRATE (not base)
        linear = d < ((ksp == 10 || ksp == 11 || ksp == 14) ? 12f0 : 18f0)
    else
        base = max(0.01f0, -0.001679f0 * d + 0.064311f0)
        modrate = min(1f0, base * fallx)               # FALLX: SNAGFALL-overridable rate correction
        linear = d < 12f0 && ksp != 2                  # small snag (redcedar=2 keeps last-5% logic)
    end
    linear && return modrate * origden
    x = (0.05f0 - 1f0) / (-modrate)                    # year at which 5% remain
    fallm2 = alldwn <= x ? 2f0 : 0.05f0 / (alldwn - x) # final fall rate (last 5%)
    if denttl <= 0.05f0 * origden
        return fallm2 * origden
    end
    dfalln = modrate * origden
    if denttl < dfalln + 0.05f0 * origden              # don't overshoot below 5% in one step
        dfalln = denttl - origden * (0.05f0 - fallm2)
    end
    return dfalln
end

"""
    snag_decay_fraction(coef, ksp) -> Float32

Annual fraction of standing-hard snags of species `ksp` that transition to soft decay
(DECAYX, fmvinit.f — e.g. a 12" tree goes soft in 2/6/10 years for snag class 1/2/3).
"""
@inline snag_decay_fraction(coef::SpeciesCoefficients, ksp::Integer) =
    coef_col(coef, :snag_decayx)[ksp]

"""
    ffe_dkr_cls(s, sp) -> Int

The fuel decay-rate class (DKRCLS, 1-4) for species `sp` — which decay-class column its dead fuel / snag
bole flows into. Prefers the FUELPOOL keyword's per-species override (FireState) over the
`fire_species_props.csv` default.
"""
@inline function ffe_dkr_cls(s::StandState, sp::Integer)::Int
    fs = s.fire
    if fs !== nothing && !isempty(fs.params.dkrcls_ovr)
        ov = get(fs.params.dkrcls_ovr, Int32(sp), Int32(0))
        ov > 0 && return Int(ov)
    end
    return Int(coef_col(s.coef, :dkr_cls)[sp])
end

"""
    add_snag!(fs, sp, dbh, density, year)

Create a standing-dead snag cohort (FMSADD) for `density` stems/acre of species `sp`,
DBH `dbh`, that died in `year`. New snags start fully hard unless the SNAGPSFT keyword set
a per-species initial-soft fraction (PSOFT), in which case that share starts soft. No-op for
non-positive density.
"""
function add_snag!(fs::FireState, sp::Integer, dbh::Float32, density::Float32, year::Integer;
                   bolevol::Float32 = 0f0, height::Float32 = 0f0, yrdead::Integer = year,
                   htcur::Float32 = -1f0, fallvol::Float32 = -1f0)
    density > 0f0 || return
    sn = fs.snags
    # SNAGPSFT: a PSOFT fraction of the new snags is soft at creation (default 0 ⇒ all hard).
    psoft = isempty(fs.params.psoft_ovr) ? 0f0 : get(fs.params.psoft_ovr, Int32(sp), 0f0)
    soft = density * psoft; hard = density - soft
    # `height` = HTDEAD (taper reference); `htcur` = current top (< HTDEAD only for pre-broken SNAGINIT snags
    # or SNAGBRK height-loss). Default htcur = height ⇒ no truncation (ordinary snags).
    hc = htcur > 0f0 ? min(htcur, height) : height
    # `fallvol` = TOTAL-volume bole for the fall→down-wood (CWD1 TVOLI='D'); default (−1) ⇒ reuse `bolevol`.
    fv = fallvol >= 0f0 ? fallvol : bolevol
    push!(sn.sp, Int32(sp));   push!(sn.dbh, dbh)
    push!(sn.den_hard, hard);  push!(sn.den_soft, soft)
    push!(sn.origden, density);  push!(sn.year, Int32(year)); push!(sn.yrdead, Int32(yrdead))
    push!(sn.bolevol, bolevol);  push!(sn.fallvol, fv);  push!(sn.height, height); push!(sn.htcur, hc)
    return
end

"""
    snag_bole_carbon(s) -> Float32

Snag STEM-VOLUME bole carbon in tons C/acre — the faithful FFE Stand-Dead snag basis
(`TOTSNG = (SNVIS+SNVIH)·V2T`, fmdout.f:153): each cohort's death-time stem-volume biomass
(`bolevol`, cuft·V2T) × its still-standing density, × 0.5. This is the bole half of Stand-Dead;
the crown half is CWD2B. (Static here — the snag height-loss that shrinks the bole over time is the
next refinement.) Falls back to Jenkins aboveground for cohorts with `bolevol` unset (e.g. fire snags).
"""
function snag_bole_carbon(s::StandState)::Float32
    fs = s.fire; fs === nothing && return 0f0
    sn = fs.snags; coef = s.coef; c = 0f0
    @inbounds for i in eachindex(sn.sp)
        den = sn.den_hard[i] + sn.den_soft[i]
        den > 0f0 || continue
        b = sn.bolevol[i]
        b <= 0f0 && (b = let (a, _, _) = jenkins_biomass(coef, sn.sp[i], sn.dbh[i]); a end)
        # SNAGBRK: a snag that lost height (htcur < HTDEAD) has a smaller bole. FVS's FMSVOL computes the
        # ORIGINAL tree (dbh, HTDEAD) TRUNCATED at htcur via CFTOPK (the Behre top-kill reduction, fmsvol.f:
        # 101-142), i.e. the fat LOWER bole — NOT a normal short tree of (dbh, htcur). Scale the stored bole by
        # the CFTOPK merch ratio. No-op at the default HTX=0 (htcur ≡ height ⇒ frozen bole, bit-exact).
        if !isempty(fs.params.snag_htx) && sn.htcur[i] < sn.height[i] && sn.height[i] > 0f0
            sp = Int(sn.sp[i]); d = sn.dbh[i]; htd = sn.height[i]
            # merch cubic (mcf_full) + total cubic (vmax) of the death-form tree (d, HTDEAD): LS/NE use the R9
            # Clark volume (v4+v7 merch, v1 total — the same basis as their live-tree merch_cuft_vol / bolevol);
            # SN uses _fm_cuft. (LS _fm_cuft returns 0 — the empty vol_eq path — which silently skipped this.)
            local mcf_full, vmax
            if s.variant isa LakeStates || s.variant isa Northeast
                ifor = Int(s.plot.forest_idx)
                fias = strip(string(coef.code_fia[sp])); fia = isempty(fias) ? 0 : parse(Int, fias)
                dbhmin, topd, scfmind, scftopd, _, _ = s.variant isa LakeStates ? _ls_merch(sp, ifor) : _ne_merch(sp, ifor)
                prod = d >= scfmind ? "01" : "02"; mtopp = d >= scfmind ? scftopd : topd
                v = r9clark_cubic(fia, d, htd, prod, mtopp, topd, 0f0)
                mcf_full = d >= dbhmin ? v[4] + v[7] : 0f0
                vmax = v[1]
            else
                mcf_full = _fm_cuft(s, sp, d, htd; merch = true)
                vmax = _fm_cuft(s, sp, d, htd; merch = false)          # v[1] total cubic (Behre vmax)
            end
            if mcf_full > 0f0
                cc = s.control
                merch_std = (stmp = cc.sp_stump_ht, topd = cc.sp_top_diam, scfstmp = cc.sp_scf_stump,
                             scftop = cc.sp_scf_topd, bftopd = cc.sp_bf_topd, bfstmp = cc.sp_bf_stump)
                bk = bark_ratio(coef, sp, d)
                _, mcf_t, _ = cftopk(merch_std, sp, d, htd, vmax, mcf_full, 0f0, vmax, bk,
                                     round(Int, sn.htcur[i] * 100f0))   # CFTOPK: merch reduced for the broken top
                b *= clamp(mcf_t / mcf_full, 0f0, 1f0)
            end
        end
        c += b * den
    end
    return c * 0.5f0
end

# CWD down-wood size class (1–9) from a stem diameter, matching the FUINI breakpoints
# (<0.25, .25–1, 1–3, 3–6, 6–12, 12–20, 20–35, 35–50, >50 inches).
@inline _cwd_size_class(d::Float32) =
    d < 0.25f0 ? 1 : d < 1f0 ? 2 : d < 3f0 ? 3 : d < 6f0 ? 4 :
    d < 12f0 ? 5 : d < 20f0 ? 6 : d < 35f0 ? 7 : d < 50f0 ? 8 : 9

# Diameter-class breakpoints BP(0:9), inches (fmcwd.f:56). Index j+1 ↔ BP(j).
const _CWD_BP = (0f0, 0.25f0, 1f0, 3f0, 6f0, 12f0, 20f0, 35f0, 50f0, 9999f0)

# Post-burn accelerated snag-fall (FMSNAG/FMSFALL). After a fire, snags present at the burn fall faster
# for PBTIME years: small (<PBSIZE in) snags lose fraction PBSMAL, soft-at-fire snags lose PBSOFT, over
# PBTIME yrs. The PB* params are SNAGPBN-overridable and live on `FireState.params` (FFEParams; SN
# defaults fmvinit.f:1100-1104). update_snags! reads them via `fs.params`.
const _FM_NZERO = 0.01f0    # NZERO: snag density treated as zero; DZERO = NZERO/50 (fmvinit.f:125)

"""
    _cwd_cone_fractions(d, ht) -> NTuple{9,Float32}

Fraction of a fallen bole's volume in each CWD size class 1..9, summing to 1 — the cone-taper
distribution of FVS's FMCWD/CWD1 (fmcwd.f label 1000). The standing stem is modeled as a cone
(DBH `d` in at 4.5 ft, total height `ht` ft); each diameter-class breakpoint BP falls at a height
BPH, and the normalized conic volume `P(h)=r(h)²·(ht−h)/(R1²·ht)` between successive breakpoints
(clipped to the [0.1 ft, ht] bole) is that class's share. Replaces the prior single-class dump
(whole bole into the DBH class), which overloaded the 6–12" class. Normalized to sum 1 so the bole
TOTAL is unchanged (carbon DDW preserved); only the per-class split — FMCFMD's (SMALL,LARGE) input —
is corrected. A bole too short to taper (`ht ≤ 4.6`) or with no taperable volume falls back to its
DBH class.
"""
function _cwd_cone_fractions(d::Float32, ht::Float32, htcur::Float32 = ht)
    f = zeros(Float32, 9)
    d <= 0.1f0 && (d = 0.1f0)
    if ht <= 4.6f0
        f[_cwd_size_class(d)] = 1f0
        return f
    end
    htd = ht
    rhrat = ((htd * 12f0) - 54f0) / (0.5f0 * d)
    # BPH(j) = height (ft) where stem diameter = BP(j); index j+1
    bph = ntuple(j -> max(0.10f0, htd - (0.5f0 * _CWD_BP[j] * rhrat) / 12f0), 10)
    r1 = d * 0.0416666667f0                       # radius (ft) at DBH (= d/12 * 0.5)
    r1 = r1 + 0.10f0 * ((r1 * htd) / (htd - 4.5f0)) # extend cone to the stem base
    r1sq = r1 * r1
    pat(h) = let r2 = r1 * (1f0 - h / htd); (r2 * r2 * (htd - h)) / (r1sq * htd) end
    # FVS CWD1 (fmcwd.f:22,29): the taper is the ORIGINAL tree (HTDEAD = `ht`), but the integration TOP is
    # HIHT(2) = HTIH = the snag's CURRENT height (`htcur`) — a broken/short snag is the fat LOWER bole of the
    # full cone, so its thin top is excluded ⇒ MORE of the bole lands in the large size classes. Normalize by
    # the FULL-cone total (pat(0.10)) NOT the truncated sum, so a truncated snag deposits only its lower
    # fraction (Σf<1) while the bolevol stays the full-tree volume. htcur == ht (default / ordinary snags,
    # SN HTX=0) ⇒ hiht = htd ⇒ Σf = 1 exactly as before (no behavior change).
    loht = 0.10f0; hiht = min(max(htcur, loht), htd)
    @inbounds for j in 1:9
        bphj = bph[j + 1]; bphjm1 = bph[j]        # BPH(j), BPH(j-1)
        (hiht <= bphj || loht > bphjm1) && continue
        hicut = min(hiht, bphjm1); locut = max(loht, bphj)
        locut == hicut && continue
        f[j] = max(0f0, pat(locut) - pat(hicut))
    end
    total_full = pat(loht)                        # P(0.10) − P(htd)=0: the whole-cone volume (normalizer)
    total_full <= 0f0 && (fill!(f, 0f0); f[_cwd_size_class(d)] = 1f0; return f)
    @inbounds for j in 1:9; f[j] /= total_full; end
    return f
end

"""
    update_snags!(s, nyears) -> Float32

Advance every snag cohort `nyears` years (FMSNAG): each year the hard snags decay toward
soft (`snag_decay_fraction`) and a `snag_fall_density` share falls — split proportionally
between the hard and soft pools (fmsnag.f:197-221). The fallen snags transfer into the
coarse-woody-debris pools (`fire.cwd`, CWD1): the fallen aboveground biomass (Jenkins ×
fallen density) is added to the down-wood class for the stem DBH and the species' decay
class. Returns the total density (stems/ac) that fell.
"""
function update_snags!(s::StandState, nyears::Integer; at_year::Union{Nothing,Integer} = nothing)::Float32
    fs = s.fire; (fs === nothing) && return 0f0
    sn = fs.snags; coef = s.coef
    cur = Int(current_cycle_year(s))            # cycle-start year — for the post-burn PBTIME window
    # `eff` is the ACTUAL annual year being stepped. The FFE annual loop (ffe_fuel_update!) advances it
    # year-by-year (at_year = cur, cur+1, …), so a snag CREATED in this cycle but BEFORE the loop — i.e. a
    # FIRE snag (fmburn! runs before ffe_fuel_update!, FVS FMBURN→annual loop) — ages 0,1,2,… across the
    # loop and falls in the years after its death, matching FVS's FMSNAG(IYR−deathyr). Ordinary-mortality
    # snags are created AFTER the loop, so they are never in it this cycle (don't fall) regardless. Default
    # (no at_year) = cur, so all other callers are unchanged.
    eff = at_year === nothing ? cur : Int(at_year)
    fallen = 0f0
    @inbounds for i in eachindex(sn.sp)
        sp = sn.sp[i]
        # Advance the snag only by the years it has actually STOOD since death (capped at this step's
        # nyears): a snag dead `eff − deathyr` years has stood that long (FMSNAG ages each snag by its own
        # (year − deathyr), not a blanket nyears). EXCEPTION — a FIRE snag is created in fmburn! BEFORE this
        # annual loop (its deathyr == this cycle's start `cur`), and FVS's FMSNAG runs AFTER FMBURN in the
        # SAME FMMAIN year, so the fire snag DOES fall in its creation year. Plain `eff−sn.year` gives yrs=0
        # at eff==sn.year and would SKIP that fall — leaving ~5× too many, since the small-snag fall is a
        # CONSTANT modrate·origden/yr (a missing final year matters hugely; fire_carbon 2005 DENIH 74 vs live
        # 15). Add the creation-year step for snags born this cycle (sn.year==cur). SAFE for carbon_snt:
        # ordinary-mortality snags are created AFTER this loop ⇒ never reach here with sn.year==cur ⇒ their
        # falldown stays bit-exact.
        born_now = Int(sn.year[i]) == cur ? 1 : 0
        yrs = clamp(eff - Int(sn.year[i]) + born_now, 0, Int(nyears))
        yrs > 0 || continue
        # a falling snag transfers its BOLE biomass to down wood; the crown is the separate CWD2B path (so
        # don't double-count it). Use the TOTAL-volume `fallvol` (FVS CWD1 TVOLI='D'=total), NOT the merch
        # `bolevol` the Stand-Dead report uses. Fall back to bolevol, then Jenkins, for cohorts with it unset.
        a = sn.fallvol[i]
        a <= 0f0 && (a = sn.bolevol[i])
        a <= 0f0 && (a = let (j, _, _) = jenkins_biomass(coef, sp, sn.dbh[i]); j end)
        idc = ffe_dkr_cls(s, sp)                            # decay-rate class (FUELPOOL-overridable)
        # Distribute the fallen bole down the cone taper across size classes (FMCWD/CWD1) instead of
        # dumping the whole bole into the DBH class. Fractions depend only on (dbh, height) → compute
        # once per cohort. Height unset (0) ⇒ single-class fallback (no behavior change).
        frac = _cwd_cone_fractions(sn.dbh[i], sn.height[i], sn.htcur[i])
        # NO hard→soft density transition for the FALL. FVS's DENIH/DENIS are the snag's INITIAL hard/soft
        # state at CREATION (all ordinary-mortality snags are created HARD → DENIH); the per-snag HARD flag
        # that flips at DKTIME (fmsnag.f:282-285) is a separate DECAY/REPORTING state and does NOT move the
        # fall density. So CWD1's DFIH/DFIS use the initial pool: a mortality snag always falls HARD (×1.00
        # SCNV), never soft. Verified against an instrumented FMSNAG: DFIS=0 every cycle for carbon_snt.
        # (jl previously moved den_hard→den_soft at DKTIME, wrongly applying the 0.80 soft factor to the
        # fall → ~13% DDW bole-fall under-count.) den_soft stays >0 only for snags SEEDED soft.
        for _ in 1:yrs
            denttl = sn.den_hard[i] + sn.den_soft[i]
            denttl > 0f0 || break
            # SNAGFALL per-species overrides of FALLX / ALLDWN (default = CSV value when not overridden).
            fx = get(fs.params.snag_fallx_ovr, Int32(sp), coef_col(coef, :snag_fallx)[sp])
            ad = get(fs.params.snag_alldwn_ovr, Int32(sp), coef_col(coef, :snag_alldwn)[sp])
            dfall = min(denttl, snag_fall_density(coef, sp, sn.dbh[i], sn.origden[i], denttl;
                                                  fallx = fx, alldwn = ad, variant = s.variant))
            dfis = denttl > 0f0 ? sn.den_soft[i] * dfall / denttl : 0f0
            dfih = denttl > 0f0 ? sn.den_hard[i] * dfall / denttl : 0f0
            # Post-burn accelerated fall (FMSNAG fmsnag.f:200-214; rates FMSFALL fmsfall.f:102-119): snags
            # that existed at a fire (died at/before BURNYR) fall faster for PBTIME years — a FLOOR (MAX)
            # on the normal fall. Small (<PBSIZE) snags fall RSMAL≈1−0.1^(1/PBTIME)≈28%/yr; soft-at-fire
            # snags fall RSOFT (PBSOFT=1 ⇒ ~all over PBTIME, via DZERO=NZERO/50); large hard snags are NOT
            # accelerated. Fire-killed snags are hard & mostly small, so the RSMAL·den_hard term drives the
            # post-fire fine down-wood pulse that was missing (fire_early 2005 lt3 1.4 vs live 13.4).
            # BURNYR is the PERSISTENT last actual burn year (FVS keeps it across cycles for the PBTIME
            # window); jl's `fire_year` is the SCHEDULED year and is cleared after firing, so derive the
            # last burn from the accumulated burn_reports instead.
            # Post-burn fall params (SNAGPBN-overridable; defaults = SN fmvinit.f:1100-1104). BURNYR is set
            # at a fire only when the scorch height exceeds PBSCOR (fmburn.f:414) — derive the last qualifying
            # burn from the accumulated burn_reports' scorch (fire_year is the scheduled year, cleared after firing).
            p = fs.params
            byr = 0
            @inbounds for br in fs.burn_reports
                br.scorch > p.pb_scor && Int(br.year) > byr && (byr = Int(br.year))
            end
            if byr > 0 && Int(sn.yrdead[i]) <= byr && 0 <= (cur - byr) <= Int(p.pb_time)
                dzr = (_FM_NZERO / 50f0) / denttl
                rsoft = p.pb_soft < 1f0 ? 1f0 - exp(log(1f0 - p.pb_soft) / p.pb_time) :
                                          1f0 - exp(log(max(1f-9, dzr)) / p.pb_time)
                pbfris = rsoft; pbfrih = 0f0
                if sn.dbh[i] < p.pb_size                         # small snags accelerate (large hard do not)
                    pbfrih = p.pb_smal < 1f0 ? 1f0 - exp(log(1f0 - p.pb_smal) / p.pb_time) :
                                               1f0 - exp(log(max(1f-9, dzr)) / p.pb_time)
                    pbfrih > pbfris && (pbfris = pbfrih)         # fmsnag.f:186-187: bump soft rate to the max
                end
                xs = pbfris * sn.den_soft[i]; xh = pbfrih * sn.den_hard[i]
                dfis < xs && (dfis = xs); dfih < xh && (dfih = xh)
            end
            sn.den_soft[i] -= dfis; sn.den_hard[i] -= dfih
            # Fallen-bole biomass into the down-wood pools, SPLIT by the snag's hard/soft state into the
            # matching CWD pool (FMCWD CWD1, fmcwd.f:K=1 soft DIS → cwd[:,1,:]; K=2 hard DIH → cwd[:,2,:]),
            # with the SCNV density conversion (fmcwd.f:61 SCNV=(0.80 soft,1.00 hard)): a SOFT (decayed)
            # snag's bole contributes 0.80× its volume. The pools decay at DIFFERENT rates (soft/index-1
            # faster ×1.1, hard/index-2 slower; fmcwd.f), so dumping all fallen bole into the hard pool
            # (as before) decayed the soft-snag boles too slowly → they accumulated as the size-5 DDW
            # overshoot. addS → soft pool, addH → hard pool.
            addS = a * dfis * 0.80f0                        # soft-snag fall → soft down-wood (index 1)
            addH = a * dfih                                 # hard-snag fall → hard down-wood (index 2)
            for j in 1:9
                if frac[j] > 0f0
                    fs.cwd[j, 1, idc] += addS * frac[j]
                    fs.cwd[j, 2, idc] += addH * frac[j]
                end
            end
            fallen += dfall
        end
    end
    return fallen
end

"""
    ffe_snag_height_loss!(s, nyears) -> nothing

SNAGBRK snag bole-breakage (FMSNGHT, SN = CASE DEFAULT fmsnght.f:153-164): shrink each snag's CURRENT
height `htcur` toward 0 over `nyears`. No-op unless SNAGBRK set per-species HTX (`snag_htx`) — at the SN
default (HTX=0) snags keep full height and the frozen `bolevol` is used (bit-exact). Per-year loss fraction
= `HTR·HTX[idx]·SFTMULT` (the HTR/HTXSFT scaling cancels the keyword's calibration), so
`htcur ← htcur·(1−lossfrac)^nyears`; the regime index is 1/3 above 0.5·HTD else 2/4, hard (SFTMULT=1) vs
soft (SFTMULT=HTXSFT, once a snag has passed DKTIME). A snag dropping below 1.5 ft becomes fuel (removed).
"""
function ffe_snag_height_loss!(s::StandState, nyears::Integer;
                               at_year::Union{Nothing,Integer} = nothing)
    fs = s.fire; (fs === nothing || isempty(fs.params.snag_htx)) && return
    sn = fs.snags; htxmap = fs.params.snag_htx
    iyr = at_year === nothing ? Int(current_cycle_year(s)) : Int(at_year)
    HTR1 = 0.1f0; HTR2 = 0.01f0; HTXSFT = 2f0    # fmvinit.f:114-115 (HTR1=first-50% rate, HTR2=after-50%)
    @inbounds for i in eachindex(sn.sp)
        (sn.den_hard[i] + sn.den_soft[i]) > 0f0 || continue
        htx = get(htxmap, Int32(sn.sp[i]), nothing); htx === nothing && continue
        htd = sn.height[i]; htc = sn.htcur[i]
        (htd > 0f0 && htc > 0f0) || continue
        # FMSNGHT picks the hard/soft rate from the snag's INITIAL state (FMSNAG calls it with IHRD=1 for the
        # DENIH/hard portion, IHRD=0 for DENIS/soft) — NOT the DKTIME report transition. jl carries one htcur,
        # so use the dominant initial pool: hard rate while den_hard ≥ den_soft (the common all-hard snag).
        soft = sn.den_soft[i] > sn.den_hard[i]
        above = htc > 0.5f0 * htd                       # height regime (>0.5·HTD uses the 50%-loss rate)
        idx = soft ? (above ? 3 : 4) : (above ? 1 : 2)   # HTINDX1/2 hard, 3/4 soft (fmsnght.f:52-59)
        sftmult = soft ? HTXSFT : 1f0
        htr = above ? HTR1 : HTR2                         # first-50% uses HTR1, after-50% uses HTR2 (fmsnght.f:154-159)
        lossfrac = clamp(htr * htx[idx] * sftmult, 0f0, 1f0)
        htnew = htc * (1f0 - lossfrac)^Float32(nyears)
        htnew < 1.5f0 && (htnew = 0f0)                   # fmsnght.f:164 — <1.5 ft ⇒ 'fuel', snag gone
        sn.htcur[i] = htnew
        htnew <= 0f0 && (sn.den_hard[i] = 0f0; sn.den_soft[i] = 0f0)
    end
    return
end

"Total standing snag density (stems/ac) currently in the snag list."
snag_standing_density(fs::FireState) = sum(fs.snags.den_hard) + sum(fs.snags.den_soft)

# Snag-summary DBH class lower bounds (SNPRCL, fmvinit.f:46-51): cumulative thresholds — class i counts
# every snag with DBH ≥ SNPRCL[i], so class 1 (≥0) equals the total.
const _FM_SNPRCL = (0f0, 12f0, 18f0, 24f0, 30f0, 36f0)

"""
    snag_summary(s) -> (; hard, soft)

The FFE snag-summary densities (stems/ac) for the FVS_SnagSum table (FMSSUM, fmssum.f:28-53): `hard`
and `soft` are each a 7-tuple — the standing HARD (`den_hard`) and SOFT (`den_soft`) snag density in the
six cumulative DBH classes (≥0/12/18/24/30/36 in) plus the total in slot 7. FVSjl's per-record
hard/soft split is the current-state equivalent of the Fortran's DENIH/DENIS + HARD flag.
"""
function snag_summary(s::StandState)
    fs = s.fire
    (fs === nothing || !fs.active) && return (hard = ntuple(_ -> 0f0, 7), soft = ntuple(_ -> 0f0, 7))
    sn = fs.snags; thd = zeros(Float32, 7); tsf = zeros(Float32, 7)
    decayx = coef_col(s.coef, :snag_decayx); iyr = Int(current_cycle_year(s))
    dcovr = fs.params.snag_decayx_ovr                       # SNAGDCAY per-species override (empty ⇒ defaults)
    @inbounds for i in eachindex(sn.sp)
        dh = sn.den_hard[i]; ds = sn.den_soft[i]; d = sn.dbh[i]
        # HARD→SOFT decay transition (fmsnag.f:282-284): once a snag has been dead ≥ DKTIME its HARD flag flips,
        # and fmssum.f:9-22 then counts its initially-hard density (DENIH) into the SOFT column. DKTIME = FMSNGDK
        # = DECAYX·(1.24·D + 13.82) for SN (fmsngdk.f DEFAULT, XMOD=1) — the same formula as the falldown TSOFT.
        # The flip is monotonic in (IYR−YRDEAD), so recomputing it at report time matches the persisted flag.
        # KNOWN UPSTREAM BUG (NOT this transition): jl dates periodic-mortality snags at the cycle-START year
        # (mortality.jl `current_cycle_year`), so at the next report they read age≈one full cycle too OLD and
        # over-trip DKTIME vs live (carbon_snag.key 1995: jl 2.9h/39.8s vs live 35.79h/6.91s; totals bit-exact).
        # Live's small/stable soft pool = only the OLD input snags aging past DKTIME. The fix is the snag YRDEAD
        # timing (FVS's annual-loop accounting), coupled to #28 — see docs/audit/BACKLOG.md item 3.
        # Use age = iyr−1−YRDEAD, NOT iyr−YRDEAD: the carbon report (FMCRBOUT, fmmain.f:206) runs BEFORE the
        # annual FMSNAG hard→soft flip (fmsnag.f:282-284), so it reflects the HARD flag as of the PREVIOUS
        # cycle's last FMSNAG year (IY(ICYC)−1). Without the −1 jl over-soften by one year on near-DKTIME snags
        # (verified vs live HARD-flag dump: dbh8.09 dkt5.01 → live age-5<5.01 HARD, jl age-6≥5.01 soft).
        dktime = get(dcovr, Int32(sn.sp[i]), decayx[sn.sp[i]]) * (1.24f0 * d + 13.82f0)
        if Float32(iyr - 1 - Int(sn.yrdead[i])) >= dktime   # TRUE YRDEAD (cycle-end−1 ord. mort.) + report 1yr behind
            ds += dh; dh = 0f0                              # initially-hard snag now reported SOFT (HARD flag false)
        end
        thd[7] += dh; tsf[7] += ds                          # slot 7 = total (all snags)
        for c in 1:6
            d >= _FM_SNPRCL[c] || continue
            thd[c] += dh; tsf[c] += ds
        end
    end
    return (hard = Tuple(thd), soft = Tuple(tsf))
end


"""
    ffe_seed_input_snags!(s) -> StandState

Seed the FFE snag list from the INPUT dead-tree records (the tree-list dead partition `n+1 : n+ndead`,
history codes 6-9) at stand initialization (FMSDIT→FMSADD ITYP=3, fmsdit.f:135). Each becomes a
standing-dead cohort: its STEM-volume bole (`cuft·V2T`) → the Stand-Dead bole, its coarse roots → the
Below-Dead BIOROOT pool. Input snags carry no crown (history ≥7 ⇒ crown already fallen; the records have
`crown_pct = 0`), so there is no CWD2B contribution. Heights/volumes are computed locally here because
`compute_volumes!` covers only the live partition (the dead records arrive with height 0). The snags are
young enough (TSOFT > the inventory age) to be hard, so the static stem bole is exact (no FMSVOL
height-loss yet). Populates the inventory-cycle Stand-Dead / Below-Dead carbon. No-op without FFE / dead records.
"""
function ffe_seed_input_snags!(s::StandState)
    fs = s.fire
    (fs === nothing || !fs.active || s.trees.ndead <= 0) && return s
    t = s.trees; coef = s.coef; c = s.control; sd = coef.species
    v2t = coef_col(coef, :v2t); ifor = Int(s.plot.forest_idx)
    # Input snags PRE-EXIST the inventory: FVS books the input dead trees at IY(1)−FINTM (a measurement
    # period before the start), so the FFE snag falldown (update_snags!, which ages by current−deathyr)
    # ages them from the inventory rather than holding them frozen-full the first cycle. Use the cycle
    # period as FINTM (the common no-GROWTH-keyword case).
    per = round(Int, c.year); per < 1 && (per = 5)
    yr = Int(current_cycle_year(s)) - per
    s.control.merch_init || init_merch_standards!(s)
    @inbounds for i in (t.n + 1):(t.n + t.ndead)
        den = t.tpa[i]; d = t.dbh[i]
        (den > 0f0 && d >= 1f0) || continue
        sp = Int(t.species[i])
        h = t.height[i] > 0f0 ? t.height[i] : max(4.5f0, _htdbh_height(sd, sp, d, ifor))
        # FVS's snag bole (FMDOUT→FMSVOL→CFVOL, fmdout.f:146) is the MERCHANTABLE cubic to the top diameter,
        # NOT the gross total-stem cubic. SN = R8 Clark v[4]; NE = R9 Clark merch v4+v7 (the live-tree
        # `merch_cuft_vol` basis). The SN R8 path returns 0 for NE (empty vol_eq) ⇒ a Jenkins-whole-tree
        # fallback over-count, so NE must use its own R9 merch model (mirrors `ffe_add_snaginit!`).
        if s.variant isa Northeast || s.variant isa LakeStates    # R9 Clark merch (LS vol_eq empty ⇒ R8 path = 0)
            fias = strip(string(coef.code_fia[sp])); fia = isempty(fias) ? 0 : parse(Int, fias)
            dbhmin, topd, scfmind, scftopd, _, _ = s.variant isa LakeStates ? _ls_merch(sp, ifor) : _ne_merch(sp, ifor)
            prod = d >= scfmind ? "01" : "02"; mtopp = d >= scfmind ? scftopd : topd
            v = r9clark_cubic(fia, d, h, prod, mtopp, topd, 0f0)
            mcuft = d >= dbhmin ? v[4] + v[7] : 0f0
        else
            prod, stump, mtopp = d >= c.sp_scf_dbhmin[sp] ?
                ("01", c.sp_scf_stump[sp], c.sp_scf_topd[sp]) : ("02", c.sp_stump_ht[sp], c.sp_top_diam[sp])
            vv, _, _ = _R8CLARK_VOL(s.species.vol_eq[sp], d, h, mtopp, c.sp_top_diam[sp], stump, prod)
            mcuft = vv[4]
        end
        bolevol = mcuft * v2t[sp] / 2000f0
        add_snag!(fs, sp, d, den, yr; bolevol = bolevol, height = h)
        _, _, rbio = jenkins_biomass(coef, sp, d)
        # FVS assumes input snags have been dead 10 years for dead-root decay (fmsadd.f:313-320):
        # XDCAY = (1−CRDCAY)^10. FVSjl was booking the full root biomass (over-counting Below-Dead).
        fs.bioroot += rbio * den * (1f0 - _FM_CRDCAY)^10
    end
    return s
end

"""
    apply_salvage!(s) -> Bool

SALVAGE (act 2520, fmsalv.f): at a scheduled cycle, remove a `PROP` fraction of standing snags within the
DBH/age bounds (and OKSOFT class), reducing den_hard/den_soft. The `PROPLV` proportion of the cut is LEFT
behind and routed to the coarse-woody-debris pools (CWD1 cone taper); the `(1−PROPLV)` remainder is removed
(harvested — its HWP-carbon FATE is a later refinement). All-species (the no-SALVSP default). Returns whether
any salvage fired this cycle. No-op without FFE / a due SALVAGE.
"""
# A snag species `sp` is in the SALVSP list (`isalvs`: 0=all, >0 single, <0 −SPGROUP).
@inline function _salv_included(s::StandState, sp::Int, isalvs::Int)::Bool
    isalvs == 0 && return true
    isalvs > 0 && return sp == isalvs
    g = -isalvs
    return 1 <= g <= length(s.control.sp_groups) && sp in s.control.sp_groups[g]
end

function apply_salvage!(s::StandState)::Bool
    fs = s.fire
    (fs === nothing || !fs.active || isempty(s.control.schedule)) && return false
    yr = Int(current_cycle_year(s)); fvscyc = Int(s.control.cycle) + 1
    sn = fs.snags; coef = s.coef; fired = false
    # SALVSP (act 2501): update the PERSISTENT species cut/leave filter when one is due this cycle.
    for a in s.control.schedule
        a.icflag == Int32(2501) || continue
        (Int(a.year) == yr || (0 < Int(a.year) < 1000 && Int(a.year) == fvscyc)) || continue
        fs.salv_isalvs = Int32(round(a.params[1])); fs.salv_isalvc = Int32(round(a.params[2]))
    end
    isalvs = Int(fs.salv_isalvs); isalvc = Int(fs.salv_isalvc)
    # TOTVOL: total volume of ALL snags before any salvage (fmsalv.f:104-121). Snapshot up front so
    # cutting doesn't shrink the denominator and multiple SALVAGE acts share it. CWDCUT = CUTVOL/TOTVOL.
    totvol = 0f0
    @inbounds for i in eachindex(sn.sp)
        dens = sn.den_hard[i] + sn.den_soft[i]
        dens > 0f0 || continue
        totvol += dens * _salv_snag_vol(coef, sn, i)
    end
    cutvol = 0f0
    for a in s.control.schedule
        a.icflag == Int32(2520) || continue
        (Int(a.year) == yr || (0 < Int(a.year) < 1000 && Int(a.year) == fvscyc)) || continue
        mindb, maxdb, maxag, oksft, prop, proplv = a.params
        oksoft = Int(oksft)
        @inbounds for i in eachindex(sn.sp)
            (sn.den_hard[i] + sn.den_soft[i]) > 0f0 || continue
            # SALVSP filter: cut-list (isalvc=0) cuts only listed species; leave-list (1) leaves them.
            linc = _salv_included(s, Int(sn.sp[i]), isalvs)
            (isalvc == 0 && !linc) && continue
            (isalvc == 1 && linc)  && continue
            d = sn.dbh[i]
            (d >= mindb && d < maxdb) || continue
            (yr - Int(sn.yrdead[i])) <= maxag || continue   # salvage age uses TRUE YRDEAD
            cuth = oksoft != 2 ? prop * sn.den_hard[i] : 0f0   # hard pool cut unless soft-only
            cuts = oksoft != 1 ? prop * sn.den_soft[i] : 0f0   # soft pool cut unless hard-only
            (cuth > 0f0 || cuts > 0f0) || continue
            cutvol += (cuth + cuts) * _salv_snag_vol(coef, sn, i)   # CUTVOL (fmsalv.f:253)
            sn.den_hard[i] = max(0f0, sn.den_hard[i] - cuth)
            sn.den_soft[i] = max(0f0, sn.den_soft[i] - cuts)
            if proplv > 0f0                                    # the left-behind share → down wood (CWD1)
                bole = sn.fallvol[i]
                bole <= 0f0 && (bole = sn.bolevol[i])
                bole <= 0f0 && (bole = let (j, _, _) = jenkins_biomass(coef, sn.sp[i], d); j end)
                idc = ffe_dkr_cls(s, sn.sp[i]); frac = _cwd_cone_fractions(d, sn.height[i], sn.htcur[i])
                addH = bole * cuth * proplv; addS = bole * cuts * proplv * 0.80f0
                for jz in 1:9
                    frac[jz] > 0f0 || continue
                    fs.cwd[jz, 2, idc] += addH * frac[jz]
                    fs.cwd[jz, 1, idc] += addS * frac[jz]
                end
            end
            fired = true
        end
    end
    # Salvaged snags' crown debris-in-waiting → down wood (fmsalv.f:301-340). Because their boles are
    # removed, a CWDCUT = CUTVOL/TOTVOL proportion of EVERY CWD2B year-pool is released to the down-wood
    # pools (P2T; foliage size-0 → litter cwd[10], woody 1-5 → cwd[1-5]) and removed from CWD2B. No /NYRS
    # here (unlike FMCADD's year-1 falldown): the whole pool releases at once when the snag is salvaged.
    if fired && totvol > 0f0 && cutvol > 0f0
        cwdcut = cutvol / totvol
        c2 = fs.cwd2b
        @inbounds for kyr in axes(c2, 3), dkcl in 1:4, sz in 0:5
            pool = c2[dkcl, sz + 1, kyr]
            pool > 0f0 || continue
            down = cwdcut * pool
            fs.cwd[sz == 0 ? 10 : sz, 2, dkcl] += down * _FM_P2T
            c2[dkcl, sz + 1, kyr] = pool - down
        end
    end
    return fired
end

# Per-record snag volume for the SALVAGE CWDCUT ratio (fmsalv.f FMSVOL total volume). Same fallback chain
# as the proplv left-behind bole: total-stem fall volume, then merch bole, then Jenkins stem biomass.
@inline function _salv_snag_vol(coef, sn, i::Int)::Float32
    v = sn.fallvol[i]
    v > 0f0 && return v
    v = sn.bolevol[i]
    v > 0f0 && return v
    j, _, _ = jenkins_biomass(coef, sn.sp[i], sn.dbh[i])
    return j
end

"""
    ffe_add_snaginit!(s) -> StandState

Add the user's SNAGINIT snags to the snag list at the first FFE year (fmsnag.f:90-105, act 2522). Each
request is `(species, DBH-at-death, ht-at-death, age, density)`: a standing-dead cohort that died `age`
years before the inventory (death year = inventory year − age), with the merchantable-cubic stem bole
(same basis as `ffe_seed_input_snags!`) and its coarse roots decayed for `age` years into the Below-Dead
BIOROOT pool. Height-at-death drives the bole + the cone taper of any later falldown. No-op without
requests / FFE.
"""
function ffe_add_snaginit!(s::StandState)
    fs = s.fire
    (fs === nothing || !fs.active || isempty(fs.snaginit)) && return s
    coef = s.coef; c = s.control; sd = coef.species
    v2t = coef_col(coef, :v2t); ifor = Int(s.plot.forest_idx)
    invyr = Int(current_cycle_year(s))
    s.control.merch_init || init_merch_standards!(s)
    @inbounds for (spf, df, htdf, htcf, agef, denf) in fs.snaginit
        sp = Int(spf); d = df; den = denf
        (sp >= 1 && d >= 1f0 && den > 0f0) || continue
        age = max(0, round(Int, agef))
        yr = invyr - age                                     # death year = inventory − AGE (fmsnag.f:99)
        h = htdf > 0f0 ? htdf : max(4.5f0, _htdbh_height(sd, sp, d, ifor))
        htc = htcf > 0f0 ? htcf : h                          # HTIH (current top) → fall-cone truncation
        # Snag stem (bole) volume at death = the variant's MERCH cubic (FMSVOL), × V2T → biomass. NE uses the
        # R9 Clark merch (v4+v7, the same basis as the live-tree `merch_cuft_vol`); SN uses R8 Clark v[4]. Using
        # the SN R8 path for NE returns 0 (the NE vol_eq is not an R8 Clark string) ⇒ snag_bole_carbon then falls
        # back to the full Jenkins ABOVEGROUND (crown+bole) ⇒ the snag carbon was ~8× too high.
        if s.variant isa Northeast || s.variant isa LakeStates
            # LS + NE use the R9 Clark merch (v4+v7) — the same basis as their live-tree merch_cuft_vol; LS
            # vol_eq is EMPTY (R9, not an R8 Clark string), so the R8 path below returns 0 ⇒ snag_bole_carbon
            # falls back to the full Jenkins aboveground (~2× too high). LS reads its own IFOR merch standards.
            fias = strip(string(coef.code_fia[sp])); fia = isempty(fias) ? 0 : parse(Int, fias)
            dbhmin, topd, scfmind, scftopd, _, _ = s.variant isa LakeStates ? _ls_merch(sp, ifor) : _ne_merch(sp, ifor)
            prod = d >= scfmind ? "01" : "02"; mtopp = d >= scfmind ? scftopd : topd
            v = r9clark_cubic(fia, d, h, prod, mtopp, topd, 0f0)
            mcuft = d >= dbhmin ? v[4] + v[7] : 0f0
            tcuft = v[1]                                         # total cubic (fall→CWD1 basis)
        else
            prod, stump, mtopp = d >= c.sp_scf_dbhmin[sp] ?
                ("01", c.sp_scf_stump[sp], c.sp_scf_topd[sp]) : ("02", c.sp_stump_ht[sp], c.sp_top_diam[sp])
            vv, _, _ = _R8CLARK_VOL(s.species.vol_eq[sp], d, h, mtopp, c.sp_top_diam[sp], stump, prod)
            mcuft = vv[4]; tcuft = vv[1]                         # merch (Stand-Dead) / total (fall→CWD1)
        end
        bolevol = mcuft * v2t[sp] / 2000f0
        fallvol = tcuft * v2t[sp] / 2000f0
        add_snag!(fs, sp, d, den, yr; bolevol = bolevol, fallvol = fallvol, height = h, htcur = htc)
        _, _, rbio = jenkins_biomass(coef, sp, d)
        fs.bioroot += rbio * den * (1f0 - _FM_CRDCAY)^age      # dead-root decay over the snag's actual age
    end
    return s
end
