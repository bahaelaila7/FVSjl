# =============================================================================
# state.jl — StandState: the explicit replacement for FVS's COMMON blocks
#
# The original FVS keeps ~52 COMMON blocks (~1,300 globals) that every subroutine
# mutates. That is what makes the faithful port unreadable and single-threaded.
# Here each COMMON block becomes a plain mutable struct, grouped onto one
# `StandState` that is passed as the first argument to (almost) every function.
#
# Consequences:
#   * no globals               (requirement #1)
#   * one state per stand/thread → no shared mutable state → trivially threadable
#                               (requirement #5)
#   * construct a state, call a function, assert on it → testable
#                               (requirement #6)
#
# Each sub-struct names the COMMON block it comes from; each field carries a short
# meaning + its Fortran name so the port can be cross-checked (requirement #2).
# Defaults reproduce the Fortran BLOCK DATA / initializer values exactly.
# =============================================================================

# ---------------------------------------------------------------------------
# ScheduledActivity — a parsed management activity (thinning/harvest) with its
# calendar date, the CUTS method code (icflag), and up to 6 numeric parameters.
# Populated by the THIN* keyword handlers; consumed by `cuts!` each cycle.
struct ScheduledActivity
    year::Int32                  # calendar year the activity fires
    icflag::Int32                # CUTS method code (8=THINDBH, …)
    params::NTuple{6,Float32}    # method parameters (post-date keyword fields)
end

# Control — COMMON /CONTRL/ + /CONCHR/ : simulation control & flags
# ---------------------------------------------------------------------------
mutable struct Control
    faithful::Bool               # bit-exact Fortran mode (true) vs apply fixes (false)

    # character (/CONCHR/)
    cf_cruise_type::String       # cubic-ft volume cruise type           (CFCTYPE)
    bf_cruise_type::String       # board-ft volume cruise type           (BFCTYPE)
    variant_code::String         # 2-char variant designator             (VARACD)
    sdi_method::String           # "ZEIDE"/"REINEKE"/blank                (CALCSDI)
    tree_format::String          # tree-record format string             (TREFMT)
    keyword_file::String         # keyword file name                     (KWDFIL)
    title::String                # run title                             (ITITLE)
    species_used::Vector{String} # 4-char per-species "used" flags       (IUSED)
    group_names::Vector{String}  # 10-char species-group names           (NAMGRP)
    ptg_names::Vector{String}    # 10-char point-group names             (PTGNAME)

    # logical flags
    auto_thin::Bool              # auto-thinning on                      (LAUTON)
    background_density::Bool     #                                       (LBKDEN)
    do_bdft_vol::Bool            # board-ft volumes requested            (LBVOLS)
    do_cuft_vol::Bool            # cubic volumes requested               (LCVOLS)
    ldcor2::Bool                 #                                       (LDCOR2)
    dub_dg::Bool                 # diameter-growth dubbing               (LDUBDG)
    is_fia::Bool                 # FIA input data                        (LFIA)
    fire_on::Bool                # fire model active                     (LFIRE)
    lflag::Bool                  # general scratch flag                  (LFLAG)
    mort_on::Bool                # mortality active                      (LMORT)
    permanent_plots::Bool        #                                       (LPERM)
    lrcor2::Bool                 #                                       (LRCOR2)
    site_set::Bool               # site index supplied                   (LSITE)
    started::Bool                # simulation started                    (LSTART)
    summary_on::Bool             #                                       (LSUMRY)
    tripling::Bool               # record tripling on                    (LTRIP)
    zeide_sdi::Bool              # use Zeide SDI                         (LZEIDE)
    mort_data::Bool              #                                       (MORDAT)
    no_tripling::Bool            #                                       (NOTRIP)
    fia_nvb::Bool                # FIA NVB biomass                       (LFIANVB)
    metric::Bool                 # metric output mode                    (LMTRIC)

    dg_calib_sp::Vector{Bool}    # per-species DG calibration            (LDGCAL)
    leave_species::Vector{Bool}  # per-species leave flag                (LEAVESP)
    ht_drag_sp::Vector{Bool}     # per-species height drag               (LHTDRG)

    # integer scalars
    error_code::Int32            # error/warning code                    (ICCODE)
    cut_algorithm::Int32         #                                       (ICFLAG)
    icl1::Int32; icl2::Int32; icl3::Int32; icl4::Int32; icl5::Int32; icl6::Int32  # cut limits
    cycle::Int32                 # current cycle number                  (ICYC)
    dg_calib_opt::Int32          # DG calibration option                 (IDG)
    first_cycle_flag::Int32      #                                       (IFST)
    read_flag::Int32             #                                       (IREAD)
    first_live_rec::Int32        # first active tree record              (IREC1)
    dead_rec_cutoff::Int32       #                                       (IREC2)
    irecnt::Int32                #                                       (IRECNT)
    irecrd::Int32                #                                       (IRECRD)
    start_year::Int32            # start date (year)                     (ISTDAT)
    ithnpa::Int32; ithnpi::Int32; ithnpn::Int32   # thinning point counters
    ntrees_active::Int32         # active tree records (mirror TreeList.n)(ITRN)
    unit_calib::Int32            # calibration output unit               (JOCALB)
    unit_list::Int32             # listing unit                         (JOLIST)
    unit_stand::Int32            # standard output unit (stdout)         (JOSTND)
    unit_summary::Int32          # summary output unit                   (JOSUM)
    unit_tree::Int32             # tree-list output unit                 (JOTREE)
    lstknt::Int32; nstknt::Int32 # stocking counters
    ncycle::Int32                # number of cycles requested            (NCYC)
    n_ptgroups::Int32            # number of point groups                (NPTGRP)
    n_spgroups::Int32            # number of species groups              (NSPGRP)
    nspecies::Int32              # species encountered in stand          (NUMSP)

    sp_begin::Vector{Int32}      # per-species begin index               (IBEGIN)
    sp_ref::Vector{Int32}        # chain-sort reference                   (IREF)
    sp_count_tab::Matrix{Int32}  # species count table (MAXSP,2)          (ISCT)
    sp_count::Vector{Int32}      # count per species                     (KOUNT)
    sp_ptr::Vector{Int32}        # pointer per species                   (KPTR)
    method_b::Vector{Int32}      #                                       (METHB)
    method_c::Vector{Int32}      #                                       (METHC)
    input_fields::Vector{Int32}  # input number fields (6)               (INS)
    ptgroup_assign::Matrix{Int32}# point-group assignments (30,52)        (IPTGRP)
    spgroup_assign::Matrix{Int32}# species-group assignments (30,92)      (ISPGRP)
    cycle_year::Vector{Int32}    # year at start of each cycle (MAXCY1)   (IY)
    table_flags::Vector{Int32}   # output table flags (7)                 (ITABLE)

    # real scalars
    auto_eff::Float32            # auto-thin cutting efficiency          (AUTEFF)
    auto_max::Float32            # auto-thin upper limit                 (AUTMAX)
    auto_min::Float32            # auto-thin lower limit                 (AUTMIN)
    ba_max::Float32              # max attainable BA                     (BAMAX)
    ba_min::Float32              #                                       (BAMIN)
    bf_min::Float32              # min harvest BF/acre                   (BFMIN)
    cc_coef::Float32             # canopy overlap coefficient            (CCCOEF)
    cc_coef2::Float32            #                                       (CCCOEF2)
    cf_min::Float32              # min harvest merch cuft/acre           (CFMIN)
    dbh_sdi::Float32             # DBH breakpoint for SDI mortality      (DBHSDI)
    dbh_stage::Float32           # min DBH for Reineke/Curtis SDI        (DBHSTAGE)
    dbh_zeide::Float32           # min DBH for Zeide SDI                 (DBHZEIDE)
    dg_sd::Float32               # DG variance bound (std dev)           (DGSD)
    zeide_dr016::Float32         #                                       (DR016)
    zeide_dr016_at::Float32      #                                       (ATDR016)
    zeide_dr016_old::Float32     #                                       (ODR016)
    cut_eff::Float32             # cutting effectiveness                 (EFF)
    mort_period::Float32         # mortality observation period (yr)     (FINTM)
    thin_ba_wt::Float32          #                                       (PBAWT)
    thin_ccf_wt::Float32         #                                       (PCCFWT)
    thin_tpa_wt::Float32         #                                       (PTPAWT)
    scf_min::Float32             #                                       (SCFMIN)
    special_wt::Float32          #                                       (SPCLWT)
    tcf_min::Float32             #                                       (TCFMIN)
    total_wt::Float32            #                                       (TCWT)
    total_removal::Float32       # total removal (TPA)                   (TRM)
    year::Float32                # current year (float)                  (YR)

    sp_dbh_min::Vector{Float32}  # min DBH for merch cubic vol           (DBHMIN)
    sp_form_class::Vector{Float32}# form class per species               (FRMCLS)
    sp_rcor2::Vector{Float32}    #                                       (RCOR2)
    sp_size_cap::Matrix{Float32} # size cap (MAXSP,4)                     (SIZCAP)
    sp_stump_ht::Vector{Float32} #                                       (STMP)
    sp_top_diam::Vector{Float32} #                                       (TOPD)
    sp_scf_dbhmin::Vector{Float32}#                                      (SCFMIND)
    sp_scf_topd::Vector{Float32} #                                       (SCFTOPD)
    sp_scf_stump::Vector{Float32}#                                       (SCFSTMP)

    schedule::Vector{ScheduledActivity}  # parsed THIN*/harvest activities (cuts!)
    years_cut::Set{Int32}                # years a thin has already been applied (idempotent cuts!)
end

function Control()
    s4(n)   = fill("    ", n)
    Control(
        true,                                                  # faithful
        " ", " ", "SN", "       ", repeat(' ',160), repeat(' ',250), repeat(' ',72),
        s4(MAXSP), fill(" "^10, 30), fill(" "^10, 30),
        false,false,false,false,false,false,false,false,false,false,false,false,
        false,false,false,false,false,false,false,false,false,
        zeros(Bool,MAXSP), zeros(Bool,MAXSP), zeros(Bool,MAXSP),
        Int32(0),Int32(0),                                      # error_code, cut_algorithm
        Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0), # icl1..6
        Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0), # cycle..start_year (9)
        Int32(0),Int32(0),Int32(0),Int32(0),                    # ithnpa,ithnpi,ithnpn,ntrees_active
        Int32(0),Int32(7),Int32(6),Int32(8),Int32(9),           # unit_calib,list,stand,summary,tree
        Int32(0),Int32(0),                                      # lstknt,nstknt
        Int32(0),Int32(0),Int32(0),Int32(0),                    # ncycle,n_ptgroups,n_spgroups,nspecies
        zeros(Int32,MAXSP), zeros(Int32,MAXSP), zeros(Int32,MAXSP,2),
        zeros(Int32,MAXSP), zeros(Int32,MAXSP), zeros(Int32,MAXSP), zeros(Int32,MAXSP),
        zeros(Int32,6), zeros(Int32,30,52), zeros(Int32,30,92),
        zeros(Int32,MAXCY1), zeros(Int32,7),
        0.0f0,60.0f0,45.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,
        0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.98f0,5.0f0,
        0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,
        zeros(Float32,MAXSP), zeros(Float32,MAXSP), zeros(Float32,MAXSP),
        zeros(Float32,MAXSP,4), zeros(Float32,MAXSP), zeros(Float32,MAXSP),
        zeros(Float32,MAXSP), zeros(Float32,MAXSP), zeros(Float32,MAXSP),
        ScheduledActivity[], Set{Int32}(),                      # schedule, years_cut
    )
end

# ---------------------------------------------------------------------------
# PlotData — COMMON /PLOT/ + /PLTCHR/ : stand/plot-level data
# ---------------------------------------------------------------------------
mutable struct PlotData
    # character (/PLTCHR/)
    mgmt_id::String              # management id                        (MGMID)
    ecoregion::String           # Bailey's ecoregion                   (ECOREG)
    pv_ref::String              # PV reference code                     (CPVREF)
    stand_id::String            # stand identification                  (NPLT)
    db_control::String          # database control number              (DBCN)
    eco_unit::String            # SN ecological unit code (e.g. 231Dd)  (PCOM)

    # integer scalars
    stand_age::Int32            # input stand age                       (IAGE)
    aspect_deg::Int32           # input aspect (deg)                    (IASPEC)
    computed_age::Int32         #                                       (ICAGE)
    county::Int32               # FIA county                           (ICNTY)
    forecast_interval::Int32    # integer years/cycle                  (IFINT)
    forest_idx::Int32           # forest subscript                     (IFOR)
    forest_type::Int32          #                                       (IFORTP)
    geo_location::Int32         # 1=N 2=C 3=S                          (IGL)
    inv_forest_type::Int32      #                                       (IIFORTP)
    model_type::Int32           #                                       (IMODTY)
    physio_region::Int32        #                                       (IPHREG)
    points_inv::Int32           #                                       (IPTINV)
    site_species::Int32         # species of max BA                    (ISISP)
    slope_raw::Int32            #                                       (ISLOP)
    n_small::Int32              # records < 3.0" DBH                   (ISMALL)
    state::Int32                # FIA state                           (ISTATE)
    stocking_class::Int32       #                                       (ISTCL)
    size_class::Int32           #                                       (ISZCL)
    habitat_input::Int32        #                                       (ITYPE)
    sp_format_default::Int32    #                                       (JSPINDEF)
    user_forest_code::Int32     #                                       (KODFOR)
    habitat_code::Int32         #                                       (KODTYP)
    managed::Int32              # 0=unmanaged 1=managed               (MANAGD)
    nonstockable::Int32         #                                       (NONSTK)
    n_site_trees::Int32         #                                       (NSITET)
    stand_origin::Int32         # 0=natural 1=plantation              (ISTDORG)

    point_ids::Vector{Int32}    # subplot id vector (MAXPLT)           (IPVEC)
    sp_format::Vector{Int32}    # per-species code format              (JSPIN)
    valid_habitat::Vector{Int32}# valid habitat codes (122)           (JTYPE)

    # real scalars
    aspect::Float32             # aspect (radians)                    (ASPECT)
    at_qmd::Float32             # QMD after thinning                  (ATAVD)
    at_avg_ht::Float32          #                                     (ATAVH)
    at_ba::Float32              #                                     (ATBA)
    at_ccf::Float32             #                                     (ATCCF)
    at_max_sdi::Float32         #                                     (ATSDIX)
    at_tpa::Float32             #                                     (ATTPA)
    avg_height::Float32         # current avg stand height            (AVH)
    basal_area::Float32         # ft²/acre                            (BA)
    baf::Float32                # basal area factor                   (BAF)
    min_dbh_var_plot::Float32   #                                     (BRK)
    before_max_sdi::Float32     #                                     (BTSDIX)
    cov_mult::Float32           #                                     (COVMLT)
    cov_year::Float32           #                                     (COVYR)
    elevation::Float32          # hundreds of feet                    (ELEV)
    cycle_length::Float32       #                                     (FINT)
    fixed_plot_inv::Float32     #                                     (FPA)
    gross_space::Float32        #                                     (GROSPC)
    old_avg_ht::Float32         #                                     (OLDAVH)
    old_ba::Float32             #                                     (OLDBA)
    old_tpa::Float32            #                                     (OLDTPA)
    old_qmd::Float32            #                                     (ORMSQD)
    pi::Float32                 #                                     (PI)
    pct_sdimax_mort_lo::Float32 #                                     (PMSDIL)
    pct_sdimax_mort_hi::Float32 #                                     (PMSDIU)
    relative_density::Float32   # current CCF                         (RELDEN)
    relative_density_prev::Float32 #                                  (RELDM1)
    mai_adj::Float32            #                                     (RMAI)
    qmd::Float32                # quadratic mean diameter             (RMSQD)
    sample_weight::Float32      #                                     (SAMWT)
    sdi_after_cut::Float32      #                                     (SDIAC)
    sdi_after_cut_z::Float32    #                                     (SDIAC2)
    sdi_before_cut::Float32     #                                     (SDIBC)
    sdi_before_cut_z::Float32   #                                     (SDIBC2)
    sdi_max::Float32            # max SDI for stand                   (SDIMAX)
    slope::Float32              # 0..1                                (SLOPE)
    site_index::Float32         # site index of site species          (STNDSI)
    total_fixed_plot::Float32   #                                     (TFPA)
    latitude::Float32           #                                     (TLAT)
    longitude::Float32          # negative = west                     (TLONG)
    total_tpa::Float32          # = sum(tpa[1:n])                     (TPROB)
    var_mult::Float32           #                                     (VMLT)
    var_mult_year::Float32      #                                     (VMLTYR)

    sp_ba_rank::Vector{Float32} #                                     (BARANK)
    sp_relden::Vector{Float32}  # CCF contribution by species         (RELDSP)
    sp_sdi_def::Vector{Float32} # max SDI by species                  (SDIDEF)
    sp_site_index::Vector{Float32} # site index by species            (SITEAR)
    site_trees::Matrix{Float32} # site tree records (MAXSTR,6)        (SITETR)
end

function PlotData()
    PlotData(
        "    ", "    ", " "^10, repeat(' ',26), repeat(' ',40), " "^10,
        Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),
        Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),
        Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),
        zeros(Int32,MAXPLT), zeros(Int32,MAXSP), zeros(Int32,122),
        0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,40.0f0,5.0f0,0.0f0,
        0.0f0,0.0f0,0.0f0,0.0f0,300.0f0,1.0f0,0.0f0,0.0f0,0.0f0,0.0f0,3.14159265f0,
        0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,
        0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,
        zeros(Float32,MAXSP), zeros(Float32,MAXSP), zeros(Float32,MAXSP),
        zeros(Float32,MAXSP), zeros(Float32,MAXSTR,6),
    )
end

# ---------------------------------------------------------------------------
# SpeciesData — variant-loaded species identity (from /PLTCHR/) + per-stand
# derived coefficients. The immutable base coefficient tables live as `const`
# arrays in the variant module (e.g. src/variants/southern/species.jl); this
# struct holds the per-stand mutable pieces. Expanded in C2/C3.
# ---------------------------------------------------------------------------
mutable struct SpeciesData
    alpha::Vector{String}        # 2/4-char alpha codes                  (JSP)
    fia::Vector{String}          # FIA species codes                     (FIAJSP)
    plants::Vector{String}       # PLANTS symbols                        (PLNJSP)
    class_codes::Matrix{String}  # species-tree class codes (MAXSP,3)     (NSP)
    dg_const::Vector{Float32}    # site-dependent DG constant per sp      (DGCON)
    vol_eq::Vector{String}       # per-species NVEL volume equation id    (VEQNNC/VEQNNB)
end
SpeciesData() = SpeciesData(
    fill("    ",MAXSP), fill("    ",MAXSP), fill("      ",MAXSP),
    fill("    ",MAXSP,3), zeros(Float32,MAXSP), fill("           ",MAXSP),
)

# ---------------------------------------------------------------------------
# Calibration — COMMON /CALCOM/ + /HTCAL/ : growth-model calibration (C3)
# ---------------------------------------------------------------------------
mutable struct Calibration
    dg_mult::Vector{Float32}     # per-species DG growth multiplier        (XDMULT)
    htg_mult::Vector{Float32}    # per-species HTG growth multiplier
    dg_const::Vector{Float32}    # site-dependent DG constant              (DGCON)
    dg_cor::Vector{Float32}      # large-tree DG calibration correction    (COR)
    htg_cor::Vector{Float32}     # height-growth calibration constant      (HTCON)
    atten::Vector{Float32}       # prior observation count (Bayes)         (ATTEN)
    sigma::Vector{Float32}       # DG residual standard deviation          (SIGMA)
    vardg::Vector{Float32}       # DG variance                             (VARDG)
    dg_cor_goal::Vector{Float32} # DG calibration attenuation goal          (WCI)
    htg_cor_small::Vector{Float32} # small-tree (REGENT) height calibration (HCOR)
    vmlt::Float32                # ARMA variance multiplier (calibration)  (VMLT)
end
Calibration() = Calibration(ones(Float32,MAXSP), ones(Float32,MAXSP),
    zeros(Float32,MAXSP), zeros(Float32,MAXSP), zeros(Float32,MAXSP),
    zeros(Float32,MAXSP), zeros(Float32,MAXSP), zeros(Float32,MAXSP),
    zeros(Float32,MAXSP), zeros(Float32,MAXSP), 0f0)

# ---------------------------------------------------------------------------
# Density — COMMON /PDEN/ : stand density / SDI scratch (C4). Minimal for now.
# ---------------------------------------------------------------------------
mutable struct Density
    sdi_sum::Float32
    point_ba::Vector{Float32}    # per-point basal area (PTBAA), indexed by subplot
    point_bal::Vector{Float32}   # per-tree BA in larger trees on its point (PTBALT)
    mort_slope::Float32          # Pretzsch self-thinning line slope             (SLPMRT)
    mort_intercept::Float32      # Pretzsch self-thinning line intercept         (CEPMRT)
    tpa_mort::Float32            # last cycle's surviving over-threshold TPA      (TPAMRT)
end
Density() = Density(0.0f0, zeros(Float32, MAXPLT), zeros(Float32, MAXTRE), 0.0f0, 0.0f0, 0.0f0)

# ---------------------------------------------------------------------------
# OutputState — COMMON /OUTCOM/ + /SUMTAB/ : summary table & output controls
# ---------------------------------------------------------------------------
mutable struct OutputState
    summary::Matrix{Float32}     # per-cycle summary rows (MAXCY1, ncols)
end
OutputState() = OutputState(zeros(Float32, MAXCY1, 22))

# ---------------------------------------------------------------------------
# Scratch — COMMON /WORKCM/ + per-tree work columns (WK1..WK15, sort indices).
# Preallocated so the hotpath never allocates (requirement #3).
# ---------------------------------------------------------------------------
mutable struct Scratch
    wk::Matrix{Float32}          # 15 per-tree work columns (15, MAXTRE)  (WK1..WK15)
    idx::Vector{Int32}           # sort index                            (IND)
    idx1::Vector{Int32}          #                                       (IND1)
    idx2::Vector{Int32}          #                                       (IND2)
end
Scratch() = Scratch(zeros(Float32,15,MAXTRE), zeros(Int32,MAXTRE), zeros(Int32,MAXTRE), zeros(Int32,MAXTRE))

# ---------------------------------------------------------------------------
# Extension states — allocated lazily only when the extension is active.
# Filled out in their own chunks (fire C7, econ C8, establishment C4).
# ---------------------------------------------------------------------------
mutable struct Establishment
    active::Bool
end
Establishment() = Establishment(false)

mutable struct DbsState
    enabled::Bool
    db_path::String
    out_db::Any                  # SQLite.DB handle (or nothing); not serialized
    report_flags::Dict{Symbol,Bool}
end
DbsState() = DbsState(false, "FVSOut.db", nothing, Dict{Symbol,Bool}())

mutable struct FireState
    active::Bool
end
FireState() = FireState(false)

mutable struct EconState
    active::Bool
end
EconState() = EconState(false)

# ---------------------------------------------------------------------------
# StandState{V} — the whole simulation state for ONE stand. Parametric on the
# variant so variant hooks dispatch at zero cost. One per thread → no contention.
# ---------------------------------------------------------------------------
mutable struct StandState{V<:AbstractVariant}
    variant::V
    coef::SpeciesCoefficients         # variant coefficients (loaded once from CSV)
    control::Control
    trees::TreeList
    plot::PlotData
    species::SpeciesData
    calib::Calibration
    density::Density
    out::OutputState
    scratch::Scratch
    rng::FVSRng
    estab::Establishment
    dbs::DbsState
    fire::Union{FireState,Nothing}
    econ::Union{EconState,Nothing}
end

"""
    StandState(variant; faithful=true)

Construct a fresh, zeroed stand state for `variant`. Fire/econ are left `nothing`
until their keywords activate them. `faithful` toggles bit-exact Fortran behaviour
(default) vs. the documented bug-fix path (see docs/DIVERGENCES.md).
"""
function StandState(variant::AbstractVariant; faithful::Bool = true)
    ctrl = Control()
    ctrl.faithful = faithful
    ctrl.variant_code = variant_code(variant)
    StandState(
        variant, coefficients(variant), ctrl, TreeList(), PlotData(), SpeciesData(), Calibration(),
        Density(), OutputState(), Scratch(), FVSRng(), Establishment(),
        DbsState(), nothing, nothing,
    )
end
