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
    age_known  ::Vector{Bool}    # whether tree age was input              (LBIRTH)
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
    old_crown_pct::Vector{Float32} # crown % previous cycle                (OLDPCT)
    old_random   ::Vector{Float32} # tree's saved random deviate           (OLDRN)
    tree_random  ::Vector{Float32} # per-tree random draw                  (ZRAND)

    # --- multi-valued attributes (k, MAXTRE) ---
    damage::Matrix{Int32}        # 6 damage-agent/severity pairs           (DAMSEV)
    pest_vars::Matrix{Int32}     # 5 pest extension variables              (IPVARS)
end

function TreeList(maxtre::Int = MAXTRE)
    iz()  = zeros(Int32,   maxtre)
    fz()  = zeros(Float32, maxtre)
    TreeList(
        0, 0,
        iz(), iz(), iz(), iz(), iz(), iz(), iz(), iz(), iz(), iz(), iz(), iz(),
        fz(), fz(), fz(), fz(), fz(), iz(), fz(), fz(), fz(),
        fz(), zeros(Bool, maxtre), fz(),
        fz(), fz(), fz(), fz(), fz(), fz(), fz(),
        fz(), fz(), fz(), fz(), fz(), fz(), fz(), fz(), fz(),
        fz(), fz(), fz(),
        zeros(Int32, 6, maxtre), zeros(Int32, 5, maxtre),
    )
end

@inline ntrees(t::TreeList) = t.n

# Per-tree vector fields copied by `copy_tree!` (every Vector field of TreeList).
const _TREE_VEC_FIELDS = (
    :species, :plot_id, :tree_id, :history, :mort_code, :cut_code, :special,
    :decay_code, :defect, :trunc, :norm_ht, :woodland_stems,
    :dbh, :height, :tpa, :diam_growth, :ht_growth, :crown_pct, :crown_ratio,
    :crown_width, :plot_size, :birth_age, :age_known, :last_diam_year,
    :bdft_vol, :cuft_vol, :merch_cuft_vol, :saw_cuft_vol, :merch_top_bf,
    :merch_top_cf, :cull, :abvgrd_bio, :merch_bio, :cubsaw_bio, :foliage_bio,
    :abvgrd_carb, :merch_carb, :cubsaw_carb, :foliage_carb, :carbon_frac,
    :old_crown_pct, :old_random, :tree_random)

"""
    copy_tree!(t, dst, src)

Copy every per-tree attribute from record `src` to record `dst` (all vector
fields plus the `damage`/`pest_vars` matrix columns). Used by record tripling.
"""
@inline function copy_tree!(t::TreeList, dst::Int, src::Int)
    @inbounds begin
        for f in _TREE_VEC_FIELDS
            getfield(t, f)[dst] = getfield(t, f)[src]
        end
        for k in 1:6; t.damage[k, dst]    = t.damage[k, src];    end
        for k in 1:5; t.pest_vars[k, dst] = t.pest_vars[k, src]; end
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
