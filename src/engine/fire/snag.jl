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
                           variant = nothing, itype::Integer = 0, kodfor::Integer = 0)::Float32
    # R6 FMSFALL (bm/ec/wc/pn/op fmsfall.f:44-46; so/fmsfall.f Oregon branch): DFALLN = BASE·FALLX·DENTTL — a fraction
    # of the CURRENT density, BASE from FMR6SDCY+FMR6FALL. No SN linear/last-5% ramp and no ALLDWN in this form.
    r6 = variant === nothing ? :none : r6_ffe_code(variant)
    (r6 === :SO && _so_california_fall(kodfor)) && (r6 = :none)
    r6 === :none || return r6_fall_base(r6, ksp, d, itype, kodfor) * fallx * denttl
    # BASE fall rate (fmsfall.f:128/130) is VARIANT-SPECIFIC: SN/CS use −0.001679·d+0.064311; LS uses the
    # "new equation" −0.006·d+0.18 (a much faster fall); NE uses an ALGSLP table (not yet ported — NE keeps
    # the SN form here). The small-snag LINEAR-fall breakpoint also differs: SN/CS = 12" (redcedar ksp2 keeps
    # the last-5% ramp); LS = 18", except cedar/tamarack (ksp 10,11,14) at 12" (fmsfall.f LS:139-145).
    if variant isa LakeStates
        modrate = clamp((-0.006f0 * d + 0.18f0) * fallx, 0.01f0, 1f0)   # FVS clamps MODRATE (not base)
        linear = d < ((ksp == 10 || ksp == 11 || ksp == 14) ? 12f0 : 18f0)
    else
        base = max(0.01f0, -0.001679f0 * d + 0.064311f0)
        if variant !== nothing && _ffe_west_fallform(variant)
            # {v}/fmsfall.f western form: linear below 18" (EM 12"), no redcedar exception; slow-falling species
            # ≥ 18" take BASE = MAX(0.01, BASE·0.32) (CI ksp 2/8, TT 3/8, UT 3/5/8).
            lmax, slow = _ffe_west_fall(variant, ksp)
            (d >= lmax && slow) && (base = max(0.01f0, base * 0.32f0))
            modrate = min(1f0, base * fallx)
            linear = d < lmax
        else
        modrate = min(1f0, base * fallx)               # FALLX: SNAGFALL-overridable rate correction
        # small snag: SN/CS keep the last-5% logic for redcedar (ksp 2).
        linear = d < 12f0 && ksp != 2
        end
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
                   htcur::Float32 = -1f0, fallvol::Float32 = -1f0, slot::Int = 0)
    density > 0f0 || return
    sn = fs.snags
    if slot > 0                                    # FMSADD reuses an emptied record (fmsadd.f:47-56, see fmsadd_bin!)
        psoft = isempty(fs.params.psoft_ovr) ? 0f0 : get(fs.params.psoft_ovr, Int32(sp), 0f0)
        sn.sp[slot] = Int32(sp); sn.dbh[slot] = dbh
        sn.den_soft[slot] = density * psoft; sn.den_hard[slot] = density - density * psoft
        sn.origden[slot] = density; sn.year[slot] = Int32(year); sn.yrdead[slot] = Int32(yrdead)
        sn.bolevol[slot] = bolevol; sn.fallvol[slot] = fallvol >= 0f0 ? fallvol : bolevol
        sn.height[slot] = height; sn.htcur[slot] = htcur > 0f0 ? min(htcur, height) : height
        sn.pbfris[slot] = 0f0; sn.pbfrih[slot] = 0f0
        return
    end
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
    push!(sn.pbfris, 0f0); push!(sn.pbfrih, 0f0)
    return
end

"An FMSADD record that no tree was binned into (fmsadd.f:47-62, :104-111): SPS/YRDEAD set, DEND=DBHS=HTDEAD=0, no density."
function add_empty_snag!(fs::FireState, sp::Integer, year::Integer; yrdead::Integer = year)
    sn = fs.snags
    push!(sn.sp, Int32(sp)); push!(sn.dbh, 0f0); push!(sn.den_hard, 0f0); push!(sn.den_soft, 0f0)
    push!(sn.origden, 0f0); push!(sn.year, Int32(year)); push!(sn.yrdead, Int32(yrdead))
    push!(sn.bolevol, 0f0); push!(sn.fallvol, 0f0); push!(sn.height, 0f0); push!(sn.htcur, 0f0)
    push!(sn.pbfris, 0f0); push!(sn.pbfrih, 0f0)
    return
end

# ── FMSADD (fmsadd.f, identical in all variants) ──────────────────────────────────────────────────────────────────
# One FMSADD call bins that call's new dead trees into snag records by (species SPCL, DBHCL = INT(D/2+1) capped 19,
# HTCL split at MIDHT when the class height range exceeds 20 ft) and gives the records slots in SPECIES-MAJOR order
# (DO SPCL / DO DBHCL / DO HTCL, fmsadd.f:41-110): when more than 3 records exist (GAPS) the first emptied record at or
# after IGAP is reused, else a new one is appended. Every variant routes every snag source (inventory ITYP=3, cut
# ITYP=2, fire/pile-burn ITYP=1, mortality ITYP=4, SNAGINIT) through here: the R6 FMR6HTLS height-loss draws are handed
# out per record (fmsnag.f:139-251), and every variant's FMDOUT/CWD sums run in record order.

mutable struct _FmsaddCtx
    gaps::Bool
    igap::Int
    taken::Vector{Int}
end
_fmsadd_ctx(fs::FireState) = _FmsaddCtx(length(fs.snags.sp) > 3, 1, Int[])   # GAPS = NSNAG > 3 (fmsadd.f:24-25)

"FMSADD record slot for the next class (fmsadd.f:47-62): an emptied, not-yet-taken record from IGAP on, else 0 = append."
function _fmsadd_slot!(fs::FireState, ctx::_FmsaddCtx)::Int
    sn = fs.snags
    if ctx.gaps
        @inbounds for j in ctx.igap:length(sn.sp)
            if sn.den_hard[j] + sn.den_soft[j] <= 0f0 && !(j in ctx.taken)
                push!(ctx.taken, j); ctx.igap = j
                return j
            end
        end
        ctx.gaps = false
    end
    return 0
end

"""
    fmsadd_bin!(s, items, year; yrdead = year, ityp, bolefn)

One FMSADD call. `items` are the call's dead trees in FVS tree order as tuples (sp, d, hsee, ht, hdead, den, htrunc):
`hsee` is FMSSEE's height (class MAXHT/MINHT), `ht` the HT(I) that picks HTCL, `hdead` the tree's HTDEAD contribution
(ITYP=3: MAX(HT,NORMHT·.01), fmsadd.f:161; else HT) and `htrunc` its broken-top height (ITYP=3 only; ≤0 = none). DBHS
and HTDEAD are density-weighted running means; for ITYP=3 the current height is overwritten by each tree (ITRUNC·.01 or
the running HTDEAD), else it is HTDEAD. `bolefn(sp, dbhs, htdead)` gives the record's (bolevol, fallvol) in tons/stem.
"""
function fmsadd_bin!(s::StandState, items, year::Integer; yrdead::Integer = year, ityp::Int, bolefn)
    fs = s.fire; isempty(items) && return s
    mn = Dict{Tuple{Int,Int},Float32}(); mx = Dict{Tuple{Int,Int},Float32}()
    for it in items
        it[6] > 0f0 || continue
        k = (it[1], _snag_dbhcl(it[2]))
        mn[k] = min(get(mn, k, 1000f0), it[3]); mx[k] = max(get(mx, k, 0f0), it[3])
    end
    keys3 = NTuple{3,Int}[]
    for k in sort!(collect(keys(mx)))                      # DO SPCL=1,MAXSP / DO DBHCL=1,19
        mx[k] > 0f0 || continue
        push!(keys3, (k[1], k[2], 1))
        (mx[k] - mn[k]) > 20f0 && push!(keys3, (k[1], k[2], 2))
    end
    rec = Dict{NTuple{3,Int},Vector{Float32}}(k => Float32[0, 0, 0, -1] for k in keys3)   # dend, dbhs, htdead, htih
    for it in items
        den = it[6]; den > 0f0 || continue
        sp = it[1]; d = it[2]; dbhcl = _snag_dbhcl(d); k2 = (sp, dbhcl)
        mid = (mx[k2] - mn[k2]) > 20f0 ? (mx[k2] + mn[k2]) / 2f0 : 0f0
        htcl = (mid <= 0f0 || it[4] < mid) ? 1 : 2
        r = rec[(sp, dbhcl, htcl)]
        tot = r[1] + den
        r[3] = (r[3] * r[1] + it[5] * den) / tot
        r[2] = (r[2] * r[1] + d * den) / tot
        r[1] = tot
        ityp == 3 && (r[4] = it[7] > 0f0 ? it[7] : r[3])
    end
    ctx = _fmsadd_ctx(fs)
    for k in keys3
        r = rec[k]
        j = _fmsadd_slot!(fs, ctx)
        # fmsadd.f:47-62 hands EVERY (SPCL,DBHCL,HTCL) class a record before the trees are binned, so a height class
        # that no tree lands in (HTCL picked on HT(I), the MIDHT split on FMSSEE's heights) still takes a slot — left
        # at DEND=0, DENIH=DENIS=0 (:196 skips it). Later FMSADD calls reuse it as a gap, and it shifts every later
        # record's FMR6HTLS draw (MEASURED FVSbm_g16 645155287126144 salvage: live NSNAG=4 with record 2 empty from
        # 2018, the 2037 snag reused slot 2; jl appended it at 4 ⇒ the height-loss draws landed on other snags).
        if r[1] <= 0f0
            j == 0 && (add_empty_snag!(fs, k[1], year; yrdead = yrdead); push!(ctx.taken, length(fs.snags.sp)))
            continue
        end
        bv, fv = bolefn(k[1], r[2], r[3])
        add_snag!(fs, k[1], r[2], r[1], year; bolevol = bv, fallvol = fv, height = r[3],
                  htcur = ityp == 3 ? r[4] : r[3], yrdead = yrdead, slot = j)
        j == 0 && push!(ctx.taken, length(fs.snags.sp))
    end
    return s
end

"The R6 FMSADD snag bole on a record's class-mean DBHS/HTDEAD: FMSVOL VOL2HT = MAX(0.005454154·H, TCF), bole == fall."
function _r6_snag_bolefn(s::StandState)
    v2t = coef_col(s.coef, :v2t)
    return (sp, d, h) -> (b = max(0.005454154f0 * h, _snag_merch_cuft_on(s, sp, d, h)) * v2t[sp] / 2000f0; (b, b))
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
    sn = fs.snags; c = 0f0
    @inbounds for i in eachindex(sn.sp)
        den = sn.den_hard[i] + sn.den_soft[i]
        den > 0f0 || continue
        c += _snag_bole_tons(s, i) * den
    end
    return c * 0.5f0
end

# The per-stem bole biomass (tons) of snag record `i` as FMDOUT's FMSVOL(I,HTIx)·V2T sees it: the stored death-time bole,
# Jenkins for a record without one, reduced for a broken top (see snag_bole_carbon).
_snag_east_vol(v) = v isa Southern || v isa CentralStates || v isa LakeStates || v isa Northeast

"""
    ffe_east_snag_vol_at(s, sp, d, htdead, xht) -> Float32 (cuft)

FMSVOL(I, XHT) for the eastern family (CS/LS/NE/SN, fmsvol.f:98-153): NATCRS on the snag record's (DBHS, HTDEAD) —
MCF = the DBH-gated merch cubic (LS/NE R9 Clark v4+v7, SN/CS R8 Clark via `_snag_merch_cuft_on`) — then, since
XHT > −1 sets LTKIL, CFTOPK at IHT = INT(XHT·100); VOL2HT = MAX(0.005454154·HTDEAD, MCF). FMDOUT/FMSOUT/FMSALV
call it fresh at every report, so a snag standing below its normal height (inventory ITRUNC/NORMHT, SNAGBRK) is
measured on the fat lower bole of its death-form tree.
"""
function ffe_east_snag_vol_at(s::StandState, sp::Int, d::Float32, htd::Float32, xht::Float32; topkill::Bool = true)::Float32
    coef = s.coef
    local mcf, vmax
    if s.variant isa LakeStates || s.variant isa Northeast
        ifor = Int(s.plot.forest_idx)
        fias = strip(string(coef.code_fia[sp])); fia = isempty(fias) ? 0 : parse(Int, fias)
        dbhmin, topd, scfmind, scftopd, _, _ = s.variant isa LakeStates ? _ls_merch(sp, ifor) : _ne_merch(sp, ifor)
        prod = d >= scfmind ? "01" : "02"; mtopp = d >= scfmind ? scftopd : topd
        v = r9clark_cubic(fia, d, htd, prod, mtopp, topd, 0f0)
        mcf = d >= dbhmin ? v[4] + v[7] : 0f0
        vmax = v[1]
    else
        mcf = _snag_merch_cuft_on(s, sp, d, htd)
        vmax = _fm_cuft(s, sp, d, htd; merch = false)
    end
    if topkill && mcf > 0f0
        cc = s.control
        merch_std = (stmp = cc.sp_stump_ht, topd = cc.sp_top_diam, scfstmp = cc.sp_scf_stump,
                     scftop = cc.sp_scf_topd, bftopd = cc.sp_bf_topd, bfstmp = cc.sp_bf_stump)
        _, mcf, _ = cftopk(merch_std, sp, d, htd, vmax, mcf, 0f0, vmax, bark_ratio(coef, sp, d),
                           unsafe_trunc(Int, xht * 100f0))
    end
    return max(0.005454154f0 * htd, mcf)
end
function _snag_bole_tons(s::StandState, i::Int)::Float32
    fs = s.fire; sn = fs.snags; coef = s.coef
    @inbounds begin
        b = sn.bolevol[i]
        b <= 0f0 && (b = let (a, _, _) = jenkins_biomass(coef, sn.sp[i], sn.dbh[i]); a end)
        # SNAGBRK: a snag that lost height (htcur < HTDEAD) has a smaller bole. FVS's FMSVOL computes the
        # ORIGINAL tree (dbh, HTDEAD) TRUNCATED at htcur via CFTOPK (the Behre top-kill reduction, fmsvol.f:
        # 101-142), i.e. the fat LOWER bole — NOT a normal short tree of (dbh, htcur). Scale the stored bole by
        # the CFTOPK merch ratio. No-op at the default HTX=0 (htcur ≡ height ⇒ frozen bole, bit-exact).
        if _ffe_west_vol(s.variant) && sn.htcur[i] < sn.height[i] && sn.height[i] > 0f0
            sp = Int(sn.sp[i])       # {v}/fmsvol.f(XHT=HTIH): NATCRS(DBHS,HTDEAD) + CFTOPK at the current height
            b = ffe_west_snag_vol_at(s, sp, sn.dbh[i], sn.height[i], sn.htcur[i]) * coef_col(coef, :v2t)[sp] / 2000f0
        elseif _snag_east_vol(s.variant) && sn.height[i] > 0f0
            # fmdout.f:139-147 → fmsvol.f (CS/LS/NE/SN): FMSVOL(I, HTIH) recomputed at every report on (DBHS, HTDEAD)
            sp = Int(sn.sp[i])
            b = ffe_east_snag_vol_at(s, sp, sn.dbh[i], sn.height[i], sn.htcur[i]) * (coef_col(coef, :v2t)[sp] / 2000f0)
        elseif !isempty(fs.params.snag_htx) && sn.htcur[i] < sn.height[i] && sn.height[i] > 0f0
            sp = Int(sn.sp[i]); d = sn.dbh[i]; htd = sn.height[i]
            mcf_full = _fm_cuft(s, sp, d, htd; merch = true)
            vmax = _fm_cuft(s, sp, d, htd; merch = false)          # v[1] total cubic (Behre vmax)
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
    end
    return b
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

# fmcwd.f label 1000 — the shared cone split behind CWD1 (a snag falls), CWD2 (a snag breaks) and CWD3 (a cut tree's
# downed yarding loss): for K=1 (soft: DIS, LOHT(1)) and K=2 (hard: DIH, LOHT(2)) each size class j gets
# DIF = MAX(0, P(LOCUT)−P(HICUT))·TVOLI, ×DEN, of the cone (R1 widened by LOHT(K) when HTD>4.5), NOT renormalized,
# booked only if DIF > 1E-6, as ADD = DIF·V2T·SCNV(K) (SCNV = .80 soft, 1.00 hard) into CWD(1,j,K,DKRCLS), in that
# REAL*4 association. A stem with HTD ≤ 4.5 (RHRAT ≤ 0) puts every breakpoint above its top and adds nothing.
# TVOLI = FMSVL2(…,'D') at (DIAM, HTD), no top-kill (`_fm_tvoli`). One routine for every variant (fmcwd.f is shared).
function _fm_cwd_split!(s::StandState, sp::Int, dbh::Float32, htd::Float32, dis::Float32, dih::Float32,
                        hiht_s::Float32, hiht_h::Float32, loht_s::Float32, loht_h::Float32;
                        tvoli::Float32 = max(0.005454154f0 * htd, _fm_tvoli(s, sp, dbh, htd)))
    (dis + dih) <= 0f0 && return
    diam = dbh <= 0.1f0 ? 0.1f0 : dbh
    rhrat = ((htd * 12f0) - 54f0) / (0.5f0 * diam)
    bph = ntuple(j -> max(0.10f0, htd - (0.5f0 * _CWD_BP[j] * rhrat) / 12f0), Val(10))   # BPH(0:9) → 1:10
    v2t = coef_col(s.coef, :v2t)[sp] / 2000f0              # fmvinit.f V2T = lb/cuft / 2000
    idc = ffe_dkr_cls(s, sp)
    @inbounds for k in 1:2
        den = k == 1 ? dis : dih
        den <= 0f0 && continue
        loht = max(0.10f0, k == 1 ? loht_s : loht_h); hiht = k == 1 ? hiht_s : hiht_h
        r1 = diam * 0.0416666667f0
        htd > 4.5f0 && (r1 = r1 + (loht * ((r1 * htd) / (htd - 4.5f0))))
        r1sq = r1 * r1
        scnv = k == 1 ? 0.80f0 : 1.00f0
        for j in 1:9
            (hiht <= bph[j + 1] || loht > bph[j]) && continue
            hicut = hiht > bph[j] ? bph[j] : hiht
            locut = loht <= bph[j + 1] ? bph[j + 1] : loht
            locut == hicut && continue
            r2 = r1 * (1f0 - (hicut / htd)); p1 = ((r2 * r2) * (htd - hicut)) / (r1sq * htd)
            r2 = r1 * (1f0 - (locut / htd)); p2 = ((r2 * r2) * (htd - locut)) / (r1sq * htd)
            dif = max(0f0, p2 - p1) * tvoli
            dif = dif * den
            dif > 1f-6 && (s.fire.cwd[j, k, idc] += dif * v2t * scnv)
        end
    end
    return
end

# FMSVL2(SP,D,H,-1,VOL,0,'D',.FALSE.) — the fallen/cut stem's volume with no top-kill: the variant's NATCRS TCF on the
# western layer (AK its NVEL F32/DVE/CUR/DEM, `ak_tree_vol`), else the snag-bole cubic of `_snag_merch_cuft_on`.
_fm_tvoli(s::StandState, sp::Int, d::Float32, h::Float32)::Float32 =
    s.variant isa SoutheastAlaska ? ak_tree_vol(s, sp, d, h)[1] :
    _ffe_west_vol(s.variant) ? ffe_west_nocut(s, sp, d, h)[1] : _snag_merch_cuft_on(s, sp, d, h)

# fmcwd.f ENTRY CWD3 (the downed yarding loss of a cut tree): HIHT(2)=HTH, LOHT(2)=.1, hard only. akffe 1993
# THINDBH WH 0.1"×2' — live CWD(1,1,2,2) 2.465E-4 (crown slash only), the renormalized split 4.747E-3.
_cwd3!(s::StandState, sp::Int, dbh::Float32, dih::Float32, hth::Float32; kw...) =
    _fm_cwd_split!(s, sp, dbh, hth, 0f0, dih, 0f0, hth, 1.0f0, 0.10f0; kw...)

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
    # Last qualifying burn year (scorch > PBSCOR) for the post-burn PBTIME fall window — snag- AND
    # year-INDEPENDENT, so compute it ONCE here. (It used to be recomputed inside the per-snag × per-year
    # loops, and `fs.burn_reports::Vector{Any}` boxes `.scorch`/`.year` on every access ⇒ ~1.4 MB per fire
    # run; hoisting drops that to one pass. The `::` asserts de-box the Any-element field reads.)
    byr = 0
    let pbscor = fs.params.pb_scor
        @inbounds for br in fs.burn_reports
            (br.scorch::Float32) > pbscor && Int(br.year::Int) > byr && (byr = Int(br.year::Int))
        end
    end
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
        a = _snag_fall_bole(s, i)
        idc = ffe_dkr_cls(s, sp)                            # decay-rate class (FUELPOOL-overridable)
        # Distribute the fallen bole down the cone taper across size classes (FMCWD/CWD1) instead of
        # dumping the whole bole into the DBH class. Fractions depend only on (dbh, height) → compute
        # once per cohort. Height unset (0) ⇒ single-class fallback (no behavior change).
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
            # fmsfall.f DFALLN is NOT capped at DENTTL: DFIS/DFIH can exceed the pools and are cut only by the
            # fmsnag.f:216-219 DZERO rule — AFTER the LASCO halving (a min() first halved the remainder instead:
            # MEASURED FVSut_g16 42642675010690 aspen snag 2025 DENIH 0.014426 → 0.003766 live, jl 0.007213).
            dfall = snag_fall_density(coef, sp, sn.dbh[i], sn.origden[i], denttl;
                                      fallx = fx, alldwn = ad, variant = s.variant,
                                      itype = _r6_itype(s), kodfor = Int(s.plot.user_forest_code))
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
            if byr > 0
                # fmsnag.f:182-190: in the burn year and the year after ((IYR−BURNYR) ≤ 1) FMSNAG STORES the record's
                # PBFRIS = RSOFT and PBFRIH (RSMAL for DBHS < PBSIZE, raised to PBFRIS for a snag whose HARD flag has
                # flipped); later years reuse them. FMSFALL (fmsfall.f:22-38) gives the rates only while IYR−BURNYR <
                # PBTIME, from this year's DENTTL when PBSOFT/PBSMAL ≥ 1 — expf/logf.
                if eff - byr <= 1
                    rsoft = 0f0; rsmal = 0f0
                    if eff - byr < p.pb_time
                        dzr = (_FM_NZERO / 50f0) / denttl
                        rsoft = p.pb_soft <= 0f0 ? 0f0 :
                                p.pb_soft < 1f0 ? 1f0 - fexp(flog(1f0 - p.pb_soft) / p.pb_time) :
                                                  1f0 - fexp(flog(dzr) / p.pb_time)
                        rsmal = p.pb_smal <= 0f0 ? 0f0 :
                                p.pb_smal < 1f0 ? 1f0 - fexp(flog(1f0 - p.pb_smal) / p.pb_time) :
                                                  1f0 - fexp(flog(dzr) / p.pb_time)
                    end
                    sn.pbfris[i] = rsoft; sn.pbfrih[i] = 0f0
                    if sn.dbh[i] < p.pb_size
                        sn.pbfrih[i] = rsmal
                        (sn.pbfrih[i] < sn.pbfris[i] && !_snag_hard_flag(s, i, eff)) && (sn.pbfrih[i] = sn.pbfris[i])
                        sn.pbfrih[i] > sn.pbfris[i] && (sn.pbfris[i] = sn.pbfrih[i])
                    end
                end
                # fmsnag.f:200-214: UT/TT/CR/BC aspen-cottonwood-birch (LASCO, fmsnag.f:131-155) snags that predate the
                # burn fall at HALF the normal rate for 10 years (still floored by the post-burn rates); every other snag
                # takes the post-burn floor for PBTIME years. MEASURED FVSut_g16 42642675010690 SIMFIRE 2020: 2030
                # Standing_Dead 1.894 live vs 1.502 jl, DDW 19.15 vs 19.49 (fire-killed aspen snags falling at full rate).
                if Int(sn.yrdead[i]) <= byr && 0 <= (eff - byr)
                    if _snag_lasco(s.variant, Int(sn.sp[i])) && (eff - byr) <= 10
                        dfis = dfis * 0.5f0; dfih = dfih * 0.5f0
                        xs = sn.pbfris[i] * sn.den_soft[i]; xh = sn.pbfrih[i] * sn.den_hard[i]
                        dfis < xs && (dfis = xs); dfih < xh && (dfih = xh)
                    elseif (eff - byr) <= p.pb_time
                        xs = sn.pbfris[i] * sn.den_soft[i]; xh = sn.pbfrih[i] * sn.den_hard[i]
                        dfis < xs && (dfis = xs); dfih < xh && (dfih = xh)
                    end
                end
            end
            # fmsnag.f:216-219 (identical in all 24 variant builds): a pool that would keep less than DZERO = NZERO/50
            # falls entirely. (jl applied it only for the R6 variants — MEASURED FVSie_g16 4769882010690 2044: live
            # emptied two records jl still carried at 1E-4/1E-5 stems/ac.)
            let dz = _FM_NZERO / 50f0
                dfis > sn.den_soft[i] - dz && (dfis = sn.den_soft[i])
                dfih > sn.den_hard[i] - dz && (dfih = sn.den_hard[i])
            end
            sn.den_soft[i] -= dfis; sn.den_hard[i] -= dfih
            # Fallen-bole biomass into the down-wood pools, SPLIT by the snag's hard/soft state into the
            # matching CWD pool (FMCWD CWD1, fmcwd.f:K=1 soft DIS → cwd[:,1,:]; K=2 hard DIH → cwd[:,2,:]),
            # with the SCNV density conversion (fmcwd.f:61 SCNV=(0.80 soft,1.00 hard)): a SOFT (decayed)
            # snag's bole contributes 0.80× its volume. The pools decay at DIFFERENT rates (soft/index-1
            # faster ×1.1, hard/index-2 slower; fmcwd.f), so dumping all fallen bole into the hard pool
            # (as before) decayed the soft-snag boles too slowly → they accumulated as the size-5 DDW
            # overshoot. addS → soft pool, addH → hard pool.
            # CWD1(I, DIH, DIS) (fmcwd.f:152-205): HIHT = HTIS/HTIH (the snag's current top), LOHT = 1.0/0.10, and TVOLI =
            # FMSVL2('D') recomputed on (DBHS, HTDEAD) as FVS does — for every variant (_fm_tvoli), not the stored fall bole.
            if sn.height[i] > 0f0
                _fm_cwd_split!(s, Int(sp), sn.dbh[i], sn.height[i], dfis, dfih, sn.htcur[i], sn.htcur[i], 1.0f0, 0.10f0;
                               tvoli = max(0.005454154f0 * sn.height[i], _fm_tvoli(s, Int(sp), sn.dbh[i], sn.height[i])))
            else    # jl-only: a snag with no recorded height (bare add_snag!) books its bole into the DBH class
                kd = _cwd_size_class(sn.dbh[i])
                fs.cwd[kd, 2, idc] += a * dfih; fs.cwd[kd, 1, idc] += a * dfis * 0.80f0
            end
            fallen += dfis + dfih
            # fmsnag.f:226-230: fewer than DZERO left in the record ⇒ it is emptied (the remnant is not added to CWD)
            if sn.den_soft[i] + sn.den_hard[i] <= _FM_NZERO / 50f0
                sn.den_soft[i] = 0f0; sn.den_hard[i] = 0f0
            end
        end
    end
    return fallen
end

# FMSNAG's per-record HARD flag as the snag's year-`iyr` processing sees it: FMSADD creates every record HARD (fmsadd.f:239)
# and FMSNAG flips it at the END of a year's pass once IYR−YRDEAD ≥ DKTIME (fmsnag.f:282-285). So at `iyr` the flag reflects
# the check made at iyr−1 — and only if FMSNAG ran that year, i.e. never before the inventory year (the inventory report
# and the first FMSNAG year still see every input snag HARD).
function _snag_hard_flag(s::StandState, i::Int, iyr::Integer)::Bool
    sn = s.fire.snags
    (iyr - 1) >= Int(s.control.cycle_year[1]) || return true
    sp = Int(sn.sp[i])
    dcx = get(s.fire.params.snag_decayx_ovr, Int32(sp), coef_col(s.coef, :snag_decayx)[sp])
    return Float32(iyr - 1 - Int(sn.yrdead[i])) < _snag_dktime(s, sp, sn.dbh[i], dcx)
end

# The per-stem fall bole (tons) CWD1/CWD2 book: the TOTAL-volume `fallvol`, else `bolevol`, else Jenkins.
@inline function _snag_fall_bole(s::StandState, i::Int)::Float32
    sn = s.fire.snags
    a = sn.fallvol[i]
    a <= 0f0 && (a = sn.bolevol[i])
    a <= 0f0 && (a = let (j, _, _) = jenkins_biomass(s.coef, sn.sp[i], sn.dbh[i]); j end)
    return a
end
# Variants whose fmvinit HTX default is 0 for every species (no height loss), so only SNAGBRK populates `snag_htx`.
_snag_htx0_default(v) = v isa Southern || v isa CentralStates
# fmsnag.f:131-155 LASCO — aspen/cottonwood/paper-birch post-burn snag-fall rule, UT/TT/CR/BC only.
_snag_lasco(v, jsp::Int) = v isa Utah ? (jsp == 6 || jsp == 18 || jsp == 19) :
                           v isa Teton ? (jsp == 6 || jsp == 15) :
                           v isa CentralRockies ? (jsp == 20 || jsp == 21 || jsp == 22 || jsp == 28) :
                           v isa BritishColumbia ? (jsp == 11 || jsp == 12 || jsp == 13 || jsp == 15) : false

# Snag first-50%-height loss rate HTR1 (fmvinit.f). LS=0.1 (faithful) and the SN SNAGBRK keyword's HTX is
# CALIBRATED against this 0.1 (HTR·HTX cancels), so the shared default stays 0.1; only NE, which seeds a RAW
# HTX=1.0 default (ne/fmvinit.f), needs its own HTR1=0.015. (SN/CS default HTX=0 ⇒ inert regardless.)
_snag_htr1(::AbstractVariant) = 0.1f0
_snag_htr1(::Northeast) = 0.015f0
_snag_htr1(::EasternMontana) = 0.0228f0   # em/fmvinit.f HTR1

# CWD2 (fmcwd.f:207-246 → the label-1000 cone split): the stem piece a snag loses to top breakage (between its new
# height LOHT and old height HIHT) goes to down wood, split across the size classes by the same cone taper as the
# CWD1 fall, DIF = MAX(0,P(LOCUT)−P(HICUT))·TVOLI with R1 widened by LOHT (fmcwd.f:347). Not normalized (FVS adds
# the raw cone slice). Enabled per variant as each is validated against live (EM first); the others still drop it.
_ffe_cwd2(v) = _ffe_west_vol(v)   # base fmcwd.f/fmsnag.f — every western-layer variant whose snags lose height

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
    fs = s.fire; fs === nothing && return
    # SN/CS keep HTX=0 (sn/fmvinit.f:1089) yet FMSNAG still calls FMSNGHT for every standing pool: HTSNEW = HTCURR, then
    # fmsnght.f:164 `HTSNEW < 1.5 ⇒ 0` — a snag shorter than 1.5 ft is broken to fuel (CWD2) and its density zeroed
    # (fmsnag.f:262-270) in its first FMSNAG year. So the loop runs with HTX=0 for species without an entry.
    (isempty(fs.params.snag_htx) && !_snag_htx0_default(s.variant)) && return
    sn = fs.snags; htxmap = fs.params.snag_htx
    iyr = at_year === nothing ? Int(current_cycle_year(s)) : Int(at_year)
    # HTR1 (first-50%-height loss rate) is VARIANT-specific (fmvinit.f): SN/CS 0.01, NE 0.015, LS 0.1. HTR2
    # (after-50%) = 0.01 all four. jl formerly hardcoded 0.1 (the LS value) — inert for NE (snag_htx empty)
    # but a latent cross-variant bug; NE now populates snag_htx (=1.0), so its HTR1 must be its own 0.015.
    HTR1 = _snag_htr1(s.variant); HTR2 = _snag_htr2(s.variant); HTXSFT = _snag_htxsft(s.variant)
    ci75 = s.variant isa CentralIdaho
    # PN/WC/BM/EC/OP (+ SO's Oregon forests) FMSNGHT call FMR6HTLS (fmsnght.f:74-93) on EVERY call — one RANN per
    # FMSNGHT, i.e. per snag record per non-empty hard/soft pool — and use its SNHTLS loss whenever the pool's HTX is
    # within [0.99,1.01] (the fmvinit default 1.0). FMSNAG brackets the whole snag loop with RANNGET/RANNPUT
    # (fmsnag.f:113-116/290-293), so the year's draws are rolled back: every year of a cycle replays the SAME
    # sequence, and the main stream is untouched.
    r6 = r6_ffe_code(s.variant)
    (r6 === :SO && _so_california_ht(Int(s.plot.user_forest_code))) && (r6 = :none)
    r6 === :AK && (r6 = :none)          # fmsnght.f: 'AK' falls to CASE DEFAULT (HTR1/HTR2·HTX), no FMR6HTLS draw
    r6save = r6 === :none ? nothing : rannget(s.rng)
    @inbounds for i in eachindex(sn.sp)
        (sn.den_hard[i] + sn.den_soft[i]) > 0f0 || continue
        x2h = 0f0; x2s = 0f0
        if r6 !== :none
            sn.den_hard[i] > 0f0 && (x2h = r6_htls(r6, Int(sn.sp[i]), rann!(s.rng)))   # FMSNGHT(…,IHRD=1,…)
            sn.den_soft[i] > 0f0 && (x2s = r6_htls(r6, Int(sn.sp[i]), rann!(s.rng)))   # FMSNGHT(…,IHRD=0,…)
        end
        htx = get(htxmap, Int32(sn.sp[i]), nothing)
        if htx === nothing
            _snag_htx0_default(s.variant) || continue
            htx = (0f0, 0f0, 0f0, 0f0)
        end
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
        # FMR6HTLS: HTX(KSP,HTINDX) in [0.99,1.01] ⇒ the pool's random loss X2 replaces HTR·HTX·SFTMULT.
        (r6 !== :none && 0.99f0 <= htx[idx] <= 1.01f0) && (lossfrac = soft ? x2s : x2h)
        htnew = htc * (1f0 - lossfrac)^Float32(nyears)
        # fmsnght.f CASE('CI') KSP 1,6 (WP, RC): lose height at HTR1 only while above 75% of HTD, then stop.
        if ci75 && (sn.sp[i] == 1 || sn.sp[i] == 6)
            if htc > 0.75f0 * htd
                lossfrac = clamp(HTR1 * htx[soft ? 3 : 1] * sftmult, 0f0, 1f0)
                htnew = htc * (1f0 - lossfrac)^Float32(nyears)
            else
                htnew = htc
            end
        end
        htnew < 1.5f0 && (htnew = 0f0)                   # fmsnght.f:164 — <1.5 ft ⇒ 'fuel', snag gone
        if _ffe_cwd2(s.variant) && htnew < htc            # CWD2: the broken-off piece → down wood (fmsnag.f:254)
            # CWD2(I, DIH, DIS, OLDHTH, OLDHTS) (fmcwd.f:207-246): HIHT = the old height, LOHT = the new one, the same
            # label-1000 split as CWD1 on TVOLI = FMSVL2('D') of (DBHS, HTDEAD)
            _fm_cwd_split!(s, Int(sn.sp[i]), sn.dbh[i], htd, sn.den_soft[i], sn.den_hard[i], htc, htc, htnew, htnew)
        elseif _snag_htx0_default(s.variant) && htnew < htc
            # SN/CS: the only CWD2 piece at HTX=0 is a <1.5-ft snag broken to fuel — same TVOLI basis as its CWD1 fall
            _fm_cwd_split!(s, Int(sn.sp[i]), sn.dbh[i], htd, sn.den_soft[i], sn.den_hard[i], htc, htc, htnew, htnew;
                           tvoli = max(0.005454154f0 * htd, _fm_tvoli(s, Int(sn.sp[i]), sn.dbh[i], htd)))
        end
        sn.htcur[i] = htnew
        htnew <= 0f0 && (sn.den_hard[i] = 0f0; sn.den_soft[i] = 0f0)
    end
    r6save === nothing || rannput!(s.rng, r6save)          # RANNPUT(SAVESO) (fmsnag.f:290-293)
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
        # YRDEAD TIMING FIXED (#28, mortality.jl:503): ordinary-mortality snags now carry the TRUE FVS death year
        # YRDEAD = IY(ICYC+1)−1 = cycle-END−1 (fmkill.f:140) in `sn.yrdead`, distinct from the cycle-START fall-
        # clock `sn.year`. This RECLASSIFIED the split from wildly-inverted (pre-fix 1995 jl 2.9h/39.8s) to
        # live-tracking (1995 35.8h/6.9s BIT-EXACT; 2000 44.57h/3.46s vs live 44.8/3.3; 2005 66.7h/4.3s). The
        # GRAND TOTAL density is bit-exact every cycle; the ONLY residual is 1–2 knife-edge cohorts at the DKTIME
        # cut. That residual is CORNERED to the grown-DBH Float32 floor (see the classification-faithfulness
        # corner below) — NOT a remaining logic bug.
        # Use age = iyr−1−YRDEAD, NOT iyr−YRDEAD: the carbon report (FMCRBOUT, fmmain.f:206) runs BEFORE the
        # annual FMSNAG hard→soft flip (fmsnag.f:282-284), so it reflects the HARD flag as of the PREVIOUS
        # cycle's last FMSNAG year (IY(ICYC)−1). Without the −1 jl over-soften by one year on near-DKTIME snags
        # (verified vs live HARD-flag dump: dbh8.09 dkt5.01 → live age-5<5.01 HARD, jl age-6≥5.01 soft).
        # CLASSIFICATION-FAITHFULNESS CORNER: with age now a clean integer (all ordinary snags share YRDEAD=cycle-
        # end−1 ⇒ report age exactly 5) and dcx a species constant, the classification `age ≥ dktime` is fully
        # FVS-faithful (true YRDEAD + iyr−1 lag + distributed dktime order below, all live-validated). The ONLY
        # input that can still differ from live is the snag's frozen death-time DBH `d` (grown-cycle accumulated).
        # If `d` were bit-exact the split would be bit-exact — so the 0.233 residual IS the grown-DBH Float32
        # accumulation: the exact flip boundary is d≈8.056 (dktime=age=5), and jl's d=8.0 sits 0.056″ below it
        # while live's same tree may cross it. Same accumulated-growth floor as MYBA/MYSDI (test_dbs_compute).
        # DKTIME must match FVS's EXACT Float32 evaluation order (fmsngdk.f:80, SN CASE DEFAULT):
        # `(1.24·DECAYX·D) + (13.82·DECAYX)`, DECAYX distributed into BOTH terms and multiplied separately —
        # NOT the factored `DECAYX·(1.24·D+13.82)`. At the age≈DKTIME near-tie boundary this sub-ULP order
        # difference flips boundary cohorts' hard/soft classification. (XMOD=1 for SN.)
        dcx = get(dcovr, Int32(sn.sp[i]), decayx[sn.sp[i]])
        dktime = _snag_dktime(s, Int(sn.sp[i]), d, dcx)
        # …and FMSNAG never ran before the inventory year, so the inventory report still sees every snag HARD.
        if (iyr - 1) >= Int(s.control.cycle_year[1]) && Float32(iyr - 1 - Int(sn.yrdead[i])) >= dktime   # TRUE YRDEAD + report 1yr behind
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

# DETAILED snag report DBH class (fmsout.f:150-153): EXCLUSIVE SNPRCL bins 1-6 (< breakpoint), distinct from
# snag_summary's CUMULATIVE (≥) classes. class i = smallest with DBH < SNPRCL[i+1]; class 6 for DBH ≥ 36.
@inline function _snag_detcl(d::Float32)::Int
    @inbounds for c in 1:5
        d < _FM_SNPRCL[c + 1] && return c
    end
    return 6
end

"""
    snag_detail(s) -> Vector of (sp, jcl, death_dbh, hth, hts, vh, vs, tv, yrdied, dh, ds, dt) rows

The DETAILED standing-snag report (fmsout.f → dbsfmdsnag.f / FVS_SnagDet): one row per non-empty
(species, YRDEAD death-year, SNPRCL DBH-class) cohort, with the density-weighted mean DBH, the hard/soft
mean CURRENT height, the hard/soft CURRENT bole cubic volume, the hard/soft/total densities, and the
death year. Uses the SAME hard→soft DKTIME flip as `snag_summary` (an initially-hard snag past DKTIME is
counted soft: its DENIH density/height/volume move to the soft totals — fmsout.f:162-172). Volume is the
per-variant snag stem cubic on the CURRENT height (`_snag_merch_cuft_on`, = FMSVOL LMERCH per variant),
matching the standing-snag bole used in `book_mortality_snags!`. Rows are sorted (sp, yrdead, jcl) for a
deterministic emission order (the SQLite content is order-independent, but the sort keeps tests stable).
"""
function snag_detail(s::StandState)
    fs = s.fire
    (fs === nothing || !fs.active) && return NamedTuple[]
    sn = fs.snags; decayx = coef_col(s.coef, :snag_decayx); iyr = Int(current_cycle_year(s))
    dcovr = fs.params.snag_decayx_ovr
    acc = Dict{Tuple{Int,Int,Int},NTuple{7,Float32}}()   # key (sp,yrdead,jcl) ⇒ (dh,ds,hth,hts,vh,vs,dbh)
    @inbounds for i in eachindex(sn.sp)
        denih = sn.den_hard[i]; denis = sn.den_soft[i]; d = sn.dbh[i]
        (denih + denis) > 0f0 || continue
        d >= _FM_SNPRCL[1] || continue                   # fmsout.f:117 DBHS < SNPRCL(1) skip
        sp = Int(sn.sp[i]); yd = Int(sn.yrdead[i]); jcl = _snag_detcl(d); h = sn.htcur[i]
        vol = _ffe_west_vol(s.variant) ? ffe_west_snag_vol_at(s, sp, d, sn.height[i], h) :   # FMSVOL(XHT=HTIH)
              _snag_east_vol(s.variant) ? ffe_east_snag_vol_at(s, sp, d, sn.height[i], h) :   # fmsout.f:122-133
              _snag_merch_cuft_on(s, sp, d, h)
        ishard = _snag_hard_flag(s, i, iyr)              # HARD(II) (fmsout.f:135; same flip as snag_summary)
        dh  = ishard ? denih : 0f0
        ds  = denis + (ishard ? 0f0 : denih)
        vh  = ishard ? vol * denih : 0f0
        vs  = vol * denis + (ishard ? 0f0 : vol * denih)
        p = get(acc, (sp, yd, jcl), (0f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0))
        acc[(sp, yd, jcl)] = (p[1] + dh, p[2] + ds, p[3] + h * dh, p[4] + h * ds,
                              p[5] + vh, p[6] + vs, p[7] + d * (denih + denis))
    end
    rows = NamedTuple[]
    for (key, v) in acc
        (sp, yd, jcl) = key
        totn = v[1] + v[2]
        totn > 0f0 || continue
        push!(rows, (sp = sp, jcl = jcl, death_dbh = v[7] / totn,
                     hth = v[1] > 0f0 ? v[3] / v[1] : 0f0, hts = v[2] > 0f0 ? v[4] / v[2] : 0f0,
                     vh = v[5], vs = v[6], tv = v[5] + v[6], yrdied = yd,
                     dh = v[1], ds = v[2], dt = totn))
    end
    sort!(rows, by = r -> (r.sp, r.yrdied, r.jcl))
    return rows
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
    c = s.control
    # Input snags PRE-EXIST the inventory: FVS books the input dead trees at YEAR = IY(1)−IFIX(FINTM)
    # (fmsdit.f:135, all 24 variants), FINTM = the MORTALITY observation period (GROWTH field 5 / FIA
    # MORT_MEASURE, default 5 — grinit.f:193), NOT the cycle length. YRDEAD drives the hard→soft DKTIME flip
    # (fmsnag.f:282): with a 10-yr cycle the old `c.year` back-dated them 10 yr, so by 2022 jl reported 14 TPA
    # of EM 196378260020004's <12" input snags soft vs live 0 (all hard until age ≥ DKTIME≈19 yr).
    yr = Int(current_cycle_year(s)) - unsafe_trunc(Int, c.growth_fintm)
    s.control.merch_init || init_merch_standards!(s)
    return _seed_input_snags_fmsadd!(s, yr)   # FMSADD ITYP=3 binning for every variant (fmsadd.f is shared)
end

# FMSADD ITYP=3 (fmsadd.f:98-364, identical in all 24 variants) — the input dead trees are BINNED into snag records
# by (species, DBHCL=INT(D/2+1) capped 19, HTCL split at MIDHT when the class height range exceeds 20 ft), exactly like
# mortality snags (book_mortality_snags!): DBHS and HTDEAD are density-weighted RUNNING means (HTDEAD over
# MAX(HT, NORMHT·.01)), while the current height HTIH/HTIS is OVERWRITTEN by each added tree — ITRUNC·.01 for a
# broken top, else the running HTDEAD (fmsadd.f:345-353). The dead records sit at the top of the tree arrays
# (MAXTRE downward), so FMSADD's I=1..MAXTRE loop visits them in REVERSE input order. The bole (FMSVOL) is on the
# record's class-mean DBHS/HTDEAD; the dead-root biomass stays per tree (×(1−CRDCAY)^10, fmsadd.f:313-320).
# EM live: LP input record 1 DBHCL 3-4 heights 51.5 (jl per-tree mean 53.125). Every variant bins through fmsadd_bin!
# (species-major record order; bole = FMSVOL on the record's DBHS/HTDEAD via _snag_merch_cuft_on, the per-variant
# volume). Tree order = FMSADD's I=1..MAXTRE over the dead block (last input record first).
function _seed_input_snags_fmsadd!(s::StandState, yr::Integer)
    fs = s.fire; t = s.trees; coef = s.coef; sd = coef.species; ifor = Int(s.plot.forest_idx)
    items = Tuple{Int,Float32,Float32,Float32,Float32,Float32,Float32}[]
    @inbounds for i in (t.n + t.ndead):-1:(t.n + 1)
        den = t.tpa[i]; d = t.dbh[i]
        (den > 0f0 && d >= 1f0) || continue
        sp = Int(t.species[i])
        h = t.height[i] > 0f0 ? t.height[i] : max(4.5f0, _htdbh_height(sd, sp, d, ifor; isne = s.variant isa Northeast))
        hd = t.norm_ht[i] > 0 ? max(h, t.norm_ht[i] * 0.01f0) : h
        hsee = t.trunc[i] > 0 ? t.trunc[i] * 0.01f0 : h       # cratet.f:482-488 FMSSEE HS (top-killed ⇒ ITRUNC·.01)
        push!(items, (sp, d, hsee, h, hd, den, t.trunc[i] > 0 ? t.trunc[i] * 0.01f0 : -1f0))
        _, _, rbio = jenkins_biomass(coef, sp, d)
        fs.bioroot += rbio * den * fpowi(1f0 - _FM_CRDCAY, 10)   # REAL**10 ⇒ libgcc __powisf2
    end
    return fmsadd_bin!(s, items, yr; ityp = 3, bolefn = _r6_snag_bolefn(s))
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
    return _fmsalv!(s)
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
        h = htdf > 0f0 ? htdf : max(4.5f0, _htdbh_height(sd, sp, d, ifor; isne = s.variant isa Northeast))
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
        elseif _ffe_west_vol(s.variant)
            mcuft = ffe_west_snag_bole(s, sp, d, h); tcuft = mcuft   # {v}/fmsvol.f MAX(X,TCF) (bole==fall==TCF)
        elseif s.variant isa Klamath
            mcuft = nc_snag_bole_cuft(s, sp, d, h); tcuft = mcuft   # NC total cubic (bole==fall==TCF)
        elseif s.variant isa OregonCoast
            # OC vol_eq are BLM Behre codes ('B…') ⇒ _R8CLARK_VOL returns 0 ⇒ the snag bole was 0 ⇒ the fall
            # fell back to the full Jenkins ABOVEGROUND (crown+bole) ⇒ ~8× too much large down-wood (the NE-class
            # bug). oc/fmsvol.f LMERCH=.FALSE. ⇒ the snag bole + CWD1 fall are the BLM TOTAL cubic (bole==fall==TCF).
            mcuft = oc_tree_cuft(sp, d, h); tcuft = mcuft
        elseif s.variant isa Olympic
            mcuft = op_tree_cuft(sp, d, h); tcuft = mcuft   # OP BLM total cubic (op/fmsvol.f LMERCH=F; same NE-class Jenkins over-book)
        else
            prod, stump, mtopp = d >= c.sp_scf_dbhmin[sp] ?
                ("01", c.sp_scf_stump[sp], c.sp_scf_topd[sp]) : ("02", c.sp_stump_ht[sp], c.sp_top_diam[sp])
            vv, ht1prd, _ = _R8CLARK_VOL(s.species.vol_eq[sp], d, h, mtopp, c.sp_top_diam[sp], stump, prod)
            # FVS FMSVOL→CFVOL MCF = v[4]+v[7] (adds topwood sawtimber-top→merch-top), NOT v[4] alone — same
            # topwood bug fixed in ffe_seed_input_snags! (all snags go through the one FMSVOL). Mirror volume.jl:531.
            mcuft = d >= c.sp_dbh_min[sp] ? vv[4] + vv[7] : 0f0
            (d >= c.sp_dbh_min[sp] && prod == "01" && ht1prd < 10f0) && (mcuft = vv[7])
            tcuft = vv[1]                                        # total (fall→CWD1)
        end
        bolevol = mcuft * v2t[sp] / 2000f0
        fallvol = tcuft * v2t[sp] / 2000f0
        # each SNAGINIT is its own FMSADD(YEAR,-JDO) call (fmsnag.f:99): it may reuse an emptied record
        slot = _fmsadd_slot!(fs, _fmsadd_ctx(fs))
        add_snag!(fs, sp, d, den, yr; bolevol = bolevol, fallvol = fallvol, height = h, htcur = htc, slot = slot)
        _, _, rbio = jenkins_biomass(coef, sp, d)
        # fmsadd.f:386-390: XDCAY=(1-CRDCAY)**PRMS(5) — a REAL exponent (libm powf), 1 when PRMS(5)<=0
        xd = agef > 0f0 ? fpow(1f0 - _FM_CRDCAY, Float32(agef)) : 1f0
        fs.bioroot += rbio * den * xd
    end
    return s
end


# fmsalv.f (one file in every build): SALVAGE with the FFE snag volumes FMSVOL(I,HTIH/HTIS) — on the western layer the
# TOTAL cubic of (DBHS,HTDEAD) trimmed by CFTOPK at the current height, elsewhere the snag's own bole — in cuft (not the
# tons-per-snag fallvol, which weights species by V2T): TOTVOL over every
# snag before any cut, CUTVOL = Σ CUTDIS·ISOFTV + CUTDIH·IHARDV, the PROPLV share left on site through CWD1 (the raw cone
# split), TONRMS += (CUTDIS·ISOFTV + CUTDIH·IHARDV)·V2T·(1−PROPLV) (→ FVS_Fuels Biomass_Removed), then a CWDCUT =
# CUTVOL/TOTVOL share of every CWD2B year-pool falls (DOWN/2000). akffe 2003: Biomass_Removed 1 (jl 0), and the
# released crowns' litter/lt3 fed the fire (Litter_Consumption 2.18447 live vs 2.18676).
function _fmsalv!(s::StandState)::Bool
    fs = s.fire; sn = fs.snags; coef = s.coef; v2t = coef_col(coef, :v2t)
    yr = Int(current_cycle_year(s)); fvscyc = Int(s.control.cycle) + 1
    fired = false
    for a in s.control.schedule
        a.icflag == Int32(2501) || continue
        (Int(a.year) == yr || (0 < Int(a.year) < 1000 && Int(a.year) == fvscyc)) || continue
        fs.salv_isalvs = Int32(round(a.params[1])); fs.salv_isalvc = Int32(round(a.params[2]))
    end
    isalvs = Int(fs.salv_isalvs); isalvc = Int(fs.salv_isalvc)
    west = _ffe_west_vol(s.variant)
    svol(i) = (west && sn.height[i] > 0f0) ?
              ffe_west_snag_vol_at(s, Int(sn.sp[i]), sn.dbh[i], sn.height[i], sn.htcur[i]; always = true) :
              _snag_bole_tons(s, i) / (v2t[sn.sp[i]] / 2000f0)
    totvol = 0f0
    @inbounds for i in eachindex(sn.sp)
        (sn.den_soft[i] + sn.den_hard[i]) > 0f0 || continue
        sn.den_soft[i] > 0f0 && (totvol = totvol + sn.den_soft[i] * svol(i))
        sn.den_hard[i] > 0f0 && (totvol = totvol + sn.den_hard[i] * svol(i))
    end
    cutvol = 0f0
    for a in s.control.schedule
        a.icflag == Int32(2520) || continue
        (Int(a.year) == yr || (0 < Int(a.year) < 1000 && Int(a.year) == fvscyc)) || continue
        mindb, maxdb, maxag, oksft, prop, proplv = a.params
        mindb = max(0f0, mindb); maxdb = min(999f0, maxdb); maxag = max(0f0, maxag)
        (oksft > 2f0 || oksft < 0f0) && (oksft = 0f0)
        oksoft = Int(trunc(oksft)); prop = min(1f0, max(0f0, prop)); proplv = min(1f0, max(0f0, proplv))
        @inbounds for i in eachindex(sn.sp)
            linc = _salv_included(s, Int(sn.sp[i]), isalvs)
            (sn.den_soft[i] + sn.den_hard[i]) <= 0f0 && continue
            (isalvc == 0 && !linc) && continue
            (isalvc == 1 && linc) && continue
            (sn.den_hard[i] <= 0f0 && oksoft == 1) && continue
            (sn.den_soft[i] <= 0f0 && oksoft == 2) && continue
            ((yr - Int(sn.yrdead[i])) > maxag || sn.dbh[i] >= maxdb || sn.dbh[i] < mindb) && continue
            isoftv = sn.den_soft[i] > 0f0 ? svol(i) : 0f0
            ihardv = sn.den_hard[i] > 0f0 ? svol(i) : 0f0
            cutdih = sn.den_hard[i] > 0f0 && oksoft != 2 ? prop * sn.den_hard[i] : 0f0   # all jl snags are HARD
            cutdis = (sn.den_soft[i] > 0f0 && oksoft != 1) ? prop * sn.den_soft[i] : 0f0
            sn.den_soft[i] = sn.den_soft[i] - cutdis; sn.den_hard[i] = sn.den_hard[i] - cutdih
            sn.den_soft[i] <= 0f0 && (sn.den_soft[i] = 0f0); sn.den_hard[i] <= 0f0 && (sn.den_hard[i] = 0f0)
            cutvol = cutvol + (cutdis * isoftv + cutdih * ihardv)
            if sn.height[i] > 0f0                                                              # CWD1(I, DIH, DIS)
                _fm_cwd_split!(s, Int(sn.sp[i]), sn.dbh[i], sn.height[i], cutdis * proplv, cutdih * proplv,
                               sn.htcur[i], sn.htcur[i], 1.0f0, 0.10f0;
                               tvoli = west ? max(0.005454154f0 * sn.height[i],
                                                  _fm_tvoli(s, Int(sn.sp[i]), sn.dbh[i], sn.height[i])) : ihardv + isoftv - (sn.den_soft[i] + cutdis > 0f0 && sn.den_hard[i] + cutdih > 0f0 ? ihardv : 0f0))
            end
            fs.tonrms = fs.tonrms + (cutdis * isoftv + cutdih * ihardv) * (v2t[sn.sp[i]] / 2000f0) * (1f0 - proplv)
            fired = true
        end
    end
    if totvol > 0f0
        cwdcut = cutvol / totvol
        c2 = fs.cwd2b
        @inbounds for kyr in axes(c2, 3), dkcl in 1:4
            for sz in 0:5
                down = cwdcut * c2[dkcl, sz + 1, kyr]
                fs.cwd[sz == 0 ? 10 : sz, 2, dkcl] += down / 2000f0
                c2[dkcl, sz + 1, kyr] = c2[dkcl, sz + 1, kyr] - down
            end
        end
    end
    return fired
end
