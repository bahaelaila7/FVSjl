# =============================================================================
# trees.jl — the per-tree record table, as a Structure-of-Arrays
#
# Ported from: common/ARRAYS.F77  (COMMON /ARRAYS/)
#
# FVS stores every tree attribute as a parallel array indexed 1..ITRN; this is
# already a Structure-of-Arrays, which is exactly what we want for cache-friendly,
# autovectorizable per-tree loops (requirement #4). We keep that layout, give the
# columns readable names, and preallocate every column to MAXTRE so the simulation
# hotpath never allocates (requirement #3).
#
# `n` is the number of active records (Fortran ITRN). Pure scratch columns
# (WK1..WK15, sort indices) live in `Scratch`, not here — this struct holds only
# genuine tree state.
#
# Each field comments its original ARRAYS name for cross-referencing the Fortran.
# =============================================================================

"""
    TreeList

Column-oriented tree records, all preallocated to `MAXTRE`. Only indices `1:n`
are active. Multi-valued attributes (damage codes, pest vars) are stored as
`(k, MAXTRE)` matrices so each tree's slice is contiguous.
"""
mutable struct TreeList
    n::Int                       # live record count                       (ITRN/IREC1)
    ndead::Int                   # dead records, stored at indices n+1:n+ndead (IREC2..)

    # --- identity / classification ---
    species   ::Vector{Int32}    # species index 1..MAXSP                  (ISP)
    plot_id   ::Vector{Int32}    # plot / point number the tree is on      (ITRE)
    tree_id   ::Vector{Int32}    # user tree id                            (IDTREE)
    history   ::Vector{Int32}    # tree history / status code              (IHISTY)
    mort_code ::Vector{Int32}    # mortality code                          (IMC)
    cut_code  ::Vector{Int32}    # cut/removal code                        (KUTKOD)
    special   ::Vector{Int32}    # special handling code                   (ISPECL)
    decay_code::Vector{Int32}    # decay class                             (DECAYCD)
    defect    ::Vector{Int32}    # defect code                             (DEFECT)
    trunc     ::Vector{Int32}    # truncation code                         (ITRUNC)
    norm_ht   ::Vector{Int32}    # normal-height flag                      (NORMHT)
    woodland_stems::Vector{Int32}# woodland stem count                     (WDLDSTEM)

    # --- core dimensions / growth ---
    dbh        ::Vector{Float32} # diameter at breast height (in)          (DBH)
    height     ::Vector{Float32} # total height (ft)                       (HT)
    tpa        ::Vector{Float32} # trees per acre (expansion factor)       (PROB)
    diam_growth::Vector{Float32} # periodic diameter growth (in)           (DG)
    ht_growth  ::Vector{Float32} # periodic height growth (ft)             (HTG)
    crown_pct  ::Vector{Int32}   # crown ratio as integer percent          (ICR)
    crown_ratio::Vector{Float32} # crown ratio as fraction 0..1            (PCT)
    crown_width::Vector{Float32} # crown width (ft)                        (CRWDTH)
    plot_size  ::Vector{Float32} # plot area the record represents (ac)    (PLTSIZ)

    # --- ages ---
    birth_age  ::Vector{Float32} # age at birth, if known                  (ABIRTH)
    age_known  ::Vector{Bool}    # age known (input or a variant site-age dub; growth uses it)
    lbirth     ::Vector{Bool}    # tree age was INPUT (intree.f:190-194 LBIRTH; set only there, TRIPLE copies it) —
                                 # the FVS_TreeList/ATRTList TreeAge gate (dbstrls.f: TREAGE = LBIRTH ? ABIRTH : 0)
    last_diam_year::Vector{Float32} # year of last observed diameter       (YRDLOS)

    # --- volumes ---
    bdft_vol      ::Vector{Float32} # board-foot volume                    (BFV)
    cuft_vol      ::Vector{Float32} # total cubic-foot volume              (CFV)
    merch_cuft_vol::Vector{Float32} # merchantable cubic volume            (MCFV)
    saw_cuft_vol  ::Vector{Float32} # sawtimber cubic volume               (SCFV)
    merch_top_bf  ::Vector{Float32} # board-foot merch top height          (HT2TD[,1])
    merch_top_cf  ::Vector{Float32} # cubic merch top height               (HT2TD[,2])
    cull          ::Vector{Float32} # cull fraction/volume                 (CULL)

    # --- biomass / carbon (per tree) ---
    abvgrd_bio ::Vector{Float32} # above-ground biomass                    (ABVGRD_BIO)
    merch_bio  ::Vector{Float32} # merchantable biomass                    (MERCH_BIO)
    cubsaw_bio ::Vector{Float32} # cubic sawtimber biomass                 (CUBSAW_BIO)
    foliage_bio::Vector{Float32} # foliage biomass                         (FOLI_BIO)
    abvgrd_carb::Vector{Float32} # above-ground carbon                     (ABVGRD_CARB)
    merch_carb ::Vector{Float32} #                                         (MERCH_CARB)
    cubsaw_carb::Vector{Float32} #                                         (CUBSAW_CARB)
    foliage_carb::Vector{Float32}#                                         (FOLI_CARB)
    carbon_frac::Vector{Float32} # carbon fraction by species              (CARB_FRAC)

    # --- per-tree state used by calibration / RNG ---
    mort_pa      ::Vector{Float32} # this period's mortality (trees/ac·GROSPC) for the record (FVS_TreeList MortPA)
    old_crown_pct::Vector{Float32} # crown % previous cycle                (OLDPCT)
    old_random   ::Vector{Float32} # tree's saved random deviate           (OLDRN)
    tree_random  ::Vector{Float32} # per-tree random draw                  (ZRAND)
    sort_key     ::Vector{Float64} # species-sort lineage key (FVS LNKCHN/TRIPLE order)
    # FFE crown-lift previous-cycle state (FMOLDC): height/DBH/crown-% at the END of the previous
    # cycle, used by `compute_crown_lift!` to find how far the crown base rose (OLDHT/OLDCRL/OLDCRW).
    # Carried through record tripling/compaction by `copy_tree!` (children inherit the parent's), so
    # the crown-lift stays aligned across the DG-tripling list growth. 0 ⇒ not yet set (no lift).
    ffe_oldht ::Vector{Float32}    # tree height previous cycle             (OLDHT)
    ffe_olddbh::Vector{Float32}    # DBH previous cycle
    ffe_oldcr ::Vector{Float32}    # crown ratio % previous cycle           (→ OLDCRL, OLDCRW)
    # Bark ratio evaluated at the START-of-cycle DBH and stashed when DG is applied (simulate.jl:438).
    # CFTOPK/BFTOPK must use THIS bark (FVS vols.f:150 `BARK=BRATIO(D)` before line 151's `D=D+DG/BARK`),
    # NOT bark recomputed from the grown DBH — the difference is the broken-top cuft residual. 0 ⇒ not
    # yet grown (cycle-0 LSTART), in which case volume.jl falls back to BRATIO(current DBH) = FVS LSTART.
    vol_bark  ::Vector{Float32}    # BRATIO(D_start) for broken-top volume top-kill
    # Previous cycle's applied diameter growth (KT mortality WK1 — the growth that produced the current
    # DBH, used as the vigor proxy in the Hamilton RIP; kt/morts.f). Snapshotted at the DBH update
    # (simulate.jl) and carried through tripling/compaction. 0 ⇒ cycle-1 (no prior growth → DG override).
    dg_prev   ::Vector{Float32}    # previous-cycle DG                       (WK1)

    # Dwarf mistletoe rating (Hawksworth 0-6), FVS MISCOM IMIST. Seeded from FIA damage codes
    # at setup (western variants only; 0 for non-DM variants). Carried through tripling/compaction.
    dmr       ::Vector{Int32}      # dwarf mistletoe rating 0..6              (IMIST)

    # Birth-cycle height-growth multiplier (live WK4=HTIMLT, estab.f:1063 HTIMLT=FTEMP/(GENTIM+0.0001)).
    # em_esgent!/*_esgent! scale the new tree's HTG by this the cycle it is established. 1.0 = full growth
    # (PLANT/existing, TRAGE≥GENTIM); AUTOES natural regen gets min(TRAGE,GENTIM)/(GENTIM+ε) < 1 (TRAGE=2 ⇒ 0.40
    # at FINT=10). Transient (used only in the birth cycle) but carried through tripling/compaction for safety.
    htimlt    ::Vector{Float32}    # birth-cycle HTG multiplier               (WK4)
    iestat     ::Vector{Int32}   # establishment 'best tree' mortality-immunity date (IESTAT; estab.f:1269 IDSDAT+20,
                                 # morts.f XCHECK; 0 = none). Carried by TRIPLE/COMPRS/TREMOV (triple.f:80, comprs.f:738).
    # Persistent per-tree small-tree height random deviate ZRAND (ie/regent.f:513-537, :1179-1200 — the TTVAR
    # BETA/ZRAND height model). −999 = "draw a fresh one": set at input (intree.f:369/607), for new establishment/
    # sprout records (estab.f:1245/1334/1424, esuckr.f:347) and COMPRS (comprs.f:973); drawn once (BACHLO(0,1),
    # |Z|≤2) in the growth OR the LSTART calibration pass and REUSED until an increment ≤0.1 resets it. Carried by
    # TRIPLE (triple.f:84) and record moves (tremov.f:46/110/153) via copy_tree!.
    zrand      ::Vector{Float32} # (ZRAND)

    # --- multi-valued attributes (k, MAXTRE) ---
    damage::Matrix{Int32}        # 6 damage-agent/severity pairs           (DAMSEV)
    pest_vars::Matrix{Int32}     # 5 pest extension variables              (IPVARS)
    # FFE per-year crown-lift OLDCRW(size 1-5) carried from the cycle it was computed (compute_crown_lift!)
    # so a tree DYING next cycle can add its YRSCYC·OLDCRW to the snag crown (FMSCRO, fmscro.f:147). idx 1=
    # size-1 … 5=size-5 (foliage excluded). 0 ⇒ no lift. Carried through tripling/compaction by copy_tree!.
    ffe_oldcrw::Matrix{Float32}  # 5 woody crown-lift sizes                (OLDCRW)
    # Per-SLOT (not per-record; never copied by copy_tree!) shadow of FVS's ICR storage ABOVE ITRN: the value
    # a vacated live slot keeps after TREDEL (tredel.f moves the last live record into a hole; the end slot's old
    # contents stay), 0 for a slot never used. bm/regent.f reads ICR(K) of a tripled copy's FUTURE slot K before
    # TRIPLE fills it (see small_tree_growth!(::BlueMountains)).
    stale_icr::Vector{Int32}
    # Same per-SLOT shadow for HT: htgf.f:292-307 caps a tripled copy's HTG against HT(ITFN) of its FUTURE slot
    # ITFN=ITRN+2I-1 before TRIPLE (grincr.f:543) writes it — i.e. the height TREDEL left there (0 if never used).
    stale_ht::Vector{Float32}
    # Per-SLOT scratch: the UNCAPPED HTGF increment (htgf.f TEMHTG) of central record I, which htgf.f hands to both
    # tripled copies before capping each against its own stale slot height (consumed by triple_records!).
    temhtg::Vector{Float32}
    # Per-SLOT LBIRTH (the TreeAge gate of dbstrls/dbscuts/dbsatrtls: TREAGE = LBIRTH(I) ? ABIRTH(I) : 0). FVS sets
    # LBIRTH only at input (intree.f:190-194, the IREC1 live slots) and in TRIPLE (triple.f:82 LBIRTH(ITFN)=LBIRTH(I));
    # TREMOV (tremov.f), COMPRS and ESTAB never touch it, so it stays with the SLOT: a record TREDEL moves into a hole
    # reports the hole's flag, and a cohort record booked in a reused slot reports the flag its old occupant left.
    # Never copied by copy_tree!. (The per-record `lbirth` is kept for the input-dead records, intree.f:564.)
    slot_lbirth::Vector{Bool}
    # FFE post-fire crown freeze (fmeff.f:492-506, fmcrow.f:126-127): a scorched record's CROWNW(I,0:5) is set to
    # TCROWN·(1−PROPCR) with GROW(I)=−1, and FMCROW skips recomputing it until GROW has climbed back to 1 (two FMCROW
    # calls, i.e. the rest of the fire cycle and the whole next one). `ffe_grow` is GROW (1 = crowns recomputed as
    # usual), `ffe_crw` the frozen CROWNW. Carried with the record by copy_tree!/TRIPLE (fmtrip.f:40) and TREDEL.
    ffe_grow::Vector{Int8}
    ffe_crw::Matrix{Float32}     # 6 crown sizes 0..5 (foliage, woody 1-5), lb/tree   (CROWNW while GROW<1)
end

function TreeList(maxtre::Int = MAXTRE)
    iz()  = zeros(Int32,   maxtre)
    fz()  = zeros(Float32, maxtre)
    dz()  = zeros(Float64, maxtre)
    TreeList(
        0, 0,
        iz(), iz(), iz(), iz(), iz(), iz(), iz(), iz(), iz(), iz(), iz(), iz(),
        fz(), fz(), fz(), fz(), fz(), iz(), fz(), fz(), fz(),
        fz(), zeros(Bool, maxtre), zeros(Bool, maxtre), fz(),
        fz(), fz(), fz(), fz(), fz(), fz(), fz(),
        fz(), fz(), fz(), fz(), fz(), fz(), fz(), fz(), fz(),
        fz(),                                  # mort_pa
        fz(), fz(), fz(), dz(), fz(), fz(), fz(),
        fz(),                                  # vol_bark
        fz(),                                  # dg_prev
        iz(),                                  # dmr
        ones(Float32, maxtre),                 # htimlt (WK4, default 1.0 = full birth-cycle growth)
        iz(),                                  # iestat
        fill(-999f0, maxtre),                  # zrand (−999 = draw on next use)
        zeros(Int32, 6, maxtre), zeros(Int32, 5, maxtre),
        zeros(Float32, 5, maxtre),              # ffe_oldcrw
        zeros(Int32, maxtre),                   # stale_icr
        zeros(Float32, maxtre),                 # stale_ht
        fill(-1f0, maxtre),                     # temhtg (−1 = not set by an HTGF cap pass this cycle)
        zeros(Bool, maxtre),                    # slot_lbirth
        ones(Int8, maxtre),                     # ffe_grow (fminit.f:972 GROW=1)
        zeros(Float32, 6, maxtre),              # ffe_crw
    )
end

@inline ntrees(t::TreeList) = t.n

# Per-tree vector fields copied by `copy_tree!` (every Vector field of TreeList).
const _TREE_VEC_FIELDS = (
    :species, :plot_id, :tree_id, :history, :mort_code, :cut_code, :special,
    :decay_code, :defect, :trunc, :norm_ht, :woodland_stems,
    :dbh, :height, :tpa, :diam_growth, :ht_growth, :crown_pct, :crown_ratio,
    :crown_width, :plot_size, :birth_age, :age_known, :lbirth, :last_diam_year,
    :bdft_vol, :cuft_vol, :merch_cuft_vol, :saw_cuft_vol, :merch_top_bf,
    :merch_top_cf, :cull, :abvgrd_bio, :merch_bio, :cubsaw_bio, :foliage_bio,
    :abvgrd_carb, :merch_carb, :cubsaw_carb, :foliage_carb, :carbon_frac,
    :mort_pa, :old_crown_pct, :old_random, :tree_random, :sort_key,
    :ffe_oldht, :ffe_olddbh, :ffe_oldcr, :vol_bark, :dg_prev, :dmr, :htimlt, :iestat, :zrand, :ffe_grow)

# Unrolled, type-stable copy of every per-tree vector field. The old `for f in _TREE_VEC_FIELDS`
# loop passed a RUNTIME Symbol to `getfield(t, f)`, whose result type is `Any` — so each copied
# value was boxed and dispatched dynamically (per field, per call: the tripling hot path's scattered
# allocations). This @generated body expands to `t.species[dst]=t.species[src]; …` with CONSTANT
# symbols → fully type-stable, zero allocation. Auto-derived from _TREE_VEC_FIELDS so it stays in sync.
@generated function _copy_tree_vecs!(t::TreeList, dst::Int, src::Int)
    assigns = [:(getfield(t, $(QuoteNode(f)))[dst] = getfield(t, $(QuoteNode(f)))[src])
               for f in _TREE_VEC_FIELDS]
    quote
        @inbounds begin
            $(assigns...)
        end
        nothing
    end
end

"""
    copy_tree!(t, dst, src)

Copy every per-tree attribute from record `src` to record `dst` (all vector
fields plus the `damage`/`pest_vars` matrix columns). Used by record tripling.
"""
@inline function copy_tree!(t::TreeList, dst::Int, src::Int)
    _copy_tree_vecs!(t, dst, src)
    @inbounds begin
        for k in 1:6; t.damage[k, dst]    = t.damage[k, src];    end
        for k in 1:5; t.pest_vars[k, dst] = t.pest_vars[k, src]; end
        for k in 1:5; t.ffe_oldcrw[k, dst] = t.ffe_oldcrw[k, src]; end
        for k in 1:6; t.ffe_crw[k, dst] = t.ffe_crw[k, src]; end
    end
    return t
end

"""
    compact_live!(t)

Drop live records whose TPA has fallen to ≤0 (e.g. removed by a thin), compacting
the remaining live records (order preserved) and shifting the dead partition to
follow. The FVS `TREDEL` after `CUTS`: removed records must not linger, else the
per-tree serial-correlation RNG sequence in the next growth pass diverges.
"""
function compact_live!(t::TreeList)
    w = 0
    @inbounds for i in 1:t.n
        if t.tpa[i] > 0f0
            w += 1
            w != i && copy_tree!(t, w, i)
        end
    end
    if t.ndead > 0 && w < t.n
        @inbounds for k in 1:t.ndead
            copy_tree!(t, w + k, t.n + k)
        end
    end
    t.n = w
    return t
end

"""
    tredel_compact!(t)

FVS `TREDEL` (tredel.f): delete tpa≤0 records by SWAP-FROM-END — fill the
smallest-index vacancy with the largest-index survivor, repeat until the vacancy
and survivor pointers cross. Unlike `compact_live!` this does NOT preserve order;
it reproduces the oracle's exact post-thin physical record layout, which any
physical-order-dependent pass (mortality kill distribution) then walks identically.
`onmove(ivac, irec)`, when given, is called after each record move — the hook through which
extensions that keep their own per-record arrays follow TREDEL (tredel.f:95 RDTDEL).
"""
function tredel_compact!(t::TreeList; thresh::Float32 = 0f0, onmove = nothing)
    n = t.n; ndel = 0
    @inbounds for i in 1:n; t.tpa[i] <= thresh && (ndel += 1); end
    ndel == 0 && return t
    iv = 1; ir = n
    @inbounds while true
        while iv <= n && t.tpa[iv] > thresh; iv += 1; end
        while ir >= 1 && t.tpa[ir] <= thresh; ir -= 1; end
        iv >= ir && break
        copy_tree!(t, iv, ir)                       # TREMOV(IVAC,IREC)
        onmove === nothing || onmove(iv, ir)        # extension per-record moves (RDTDEL, tredel.f:95)
        t.tpa[ir] = 0f0; iv += 1; ir -= 1
    end
    newn = n - ndel
    # FVS leaves the vacated slots newn+1:n holding their old records (moved-from or deleted); remember their ICR
    # before jl slides the dead partition down over them (FVS keeps the dead at MAXTRE, not here).
    @inbounds for k in (newn + 1):n; t.stale_icr[k] = t.crown_pct[k]; t.stale_ht[k] = t.height[k]; end
    if t.ndead > 0
        @inbounds for k in 1:t.ndead; copy_tree!(t, newn + k, n + k); end
    end
    t.n = newn
    # A TREDEL removal triggers Fortran's species-sort REBUILD (spesrt→lnkchn→setup),
    # which re-lists each species in ASCENDING PHYSICAL record index — and crucially
    # does NOT re-run REASS (the post-TRIPLE U,C,L interleave, grincr.f:553), so the
    # tripled lineage order is discarded after a thin/comcup and replaced by physical
    # order. Reset sort_key to the compacted physical position so species_sort! (and
    # thus the DGSCOR per-tree RNG assignment) tracks the oracle post-removal.
    @inbounds for i in 1:newn; t.sort_key[i] = Float64(i); end
    return t
end

"""
    comcup!(t)

COMCUP (comcup.f, top half): each cycle, delete live records whose expansion factor
`PROB` (tpa) has fallen to ≤ 1e-5 — suppressed trees whittled to ~0 by mortality.
Uses the same swap-from-end TREDEL as a thin. Keeping them is harmless to the `.sum`
(they contribute ~0), but they would consume an extra per-tree DGSCOR random deviate
next cycle, drifting the RNG from the oracle. (The COMPRESS-keyword compression in the
bottom half of comcup.f is a management option — not ported here.)
"""
comcup!(t::TreeList; onmove = nothing) = tredel_compact!(t; thresh = 1f-5, onmove = onmove)

"""
    permute_records!(t, lo, perm)

Reorder the contiguous record block `lo:lo+length(perm)-1` so that position `lo+k-1` receives the record that
was at `perm[k]` (every per-tree vector field plus the damage/pest_vars/ffe_oldcrw columns). `perm` must be a
permutation of that block. Used to lay newly established records out in FVS's per-plot booking order.
"""
function permute_records!(t::TreeList, lo::Int, perm::Vector{Int})
    hi = lo + length(perm) - 1
    for f in _TREE_VEC_FIELDS
        v = getfield(t, f)
        tmp = v[perm]
        @inbounds for k in eachindex(perm); v[lo + k - 1] = tmp[k]; end
    end
    t.damage[:, lo:hi]     = t.damage[:, perm]
    t.pest_vars[:, lo:hi]  = t.pest_vars[:, perm]
    t.ffe_oldcrw[:, lo:hi] = t.ffe_oldcrw[:, perm]
    t.ffe_crw[:, lo:hi]    = t.ffe_crw[:, perm]
    return t
end

"""
    spesrt_reorder!(t)

Reset the species-sort lineage key (`sort_key`) to each record's ascending physical
index — the effect of FVS's `SPESRT` rebuild (spesrt.f→lnkchn.f→setup.f), which
re-lists every species in ASCENDING PHYSICAL record order and DISCARDS the post-TRIPLE
`REASS` (upper,central,lower) lineage interleave (grincr.f:553). FVS calls SPESRT
whenever the tree list is modified: after a TREDEL removal (already handled inside
`tredel_compact!`) and — crucially — unconditionally after the establishment model
generates records (esgent.f:49, esnutr.f:138). A stand with active establishment
therefore enters the NEXT growth cycle's DGSCOR in ascending physical order, NOT the
REASS lineage order. (For a NOAUTOES stand no establishment SPESRT fires, so the REASS
interleave survives to the next cycle — which is why the eastern/NOAUTOES stochastic-DG
path keeps needing `species_sort!`'s lineage sort.)
"""
function spesrt_reorder!(t::TreeList)
    @inbounds for i in 1:t.n; t.sort_key[i] = Float64(i); end
    return t
end
