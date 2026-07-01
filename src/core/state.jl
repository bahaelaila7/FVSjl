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
    aux::Float32                 # 7th param for keywords that need it (THINQFA QFATAR)
end
# 3-arg form (the common case): aux defaults to 0.
ScheduledActivity(year, icflag, params) =
    ScheduledActivity(Int32(year), Int32(icflag), params, 0f0)

# GrowthMultiplier — a keyword growth/mortality multiplier (MULTS, base/mults.f).
# `kind` ∈ (:bai,:htg,:regh,:mort,:regd); applies to `species` (0 = all) from `year`
# onward (the most recent matching one wins). `value` is the per-species multiplier.
# `d1`/`d2` are the DBH window for MORTMULT (morts.f:518: X applied only if d1≤DBH<d2;
# defaults 0/99999 = all trees); unused by the other kinds.
struct GrowthMultiplier
    kind::Symbol
    year::Int32
    species::Int32        # 0 = all species
    value::Float32
    d1::Float32           # MORTMULT DBH window low  (XMDIA1, default 0)
    d2::Float32           # MORTMULT DBH window high (XMDIA2, default 99999)
end
GrowthMultiplier(kind, year, species, value) =
    GrowthMultiplier(kind, year, species, value, 0f0, 99999f0)

# Event-monitor expression AST node (concrete types + evaluator in event_monitor.jl).
abstract type EvNode end

# ConditionalActivity — an IF/THEN/ENDIF block (event monitor): its activities fire
# only in cycles where the condition expression is true. Built by the IF keyword
# handler, evaluated each cycle in `cuts!` via `eval_event` over the parsed AST.
struct ConditionalActivity
    cond::EvNode
    acts::Vector{ScheduledActivity}
    src::String
end

"""
One `ESTUMP` cut record (estump.f) for stump-sprout regeneration. Built per removed
tree of a *sprouting* species, in removal order, and consumed by `esuckr!`:
`species` (FVS species code), `dstmp` (stump DBH = the cut tree's DBH), `prem`
(TPA removed, less standing snags), `plot` (point id), `ishag` (sprout age = the
cycle length IFINT, carried into SPRTHT/ABIRTH).
"""
const CutRecord = @NamedTuple{species::Int32, dstmp::Float32, prem::Float32,
                             plot::Int32, ishag::Int32}

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
    # MORTMSB alternate-mortality (mature-stand breakup): inert by default (QMDMSB=999 ⇒ D10>QMDMSB never).
    msb_qmd::Float32             # QMD threshold above which MSB fires    (QMDMSB; default 999)
    msb_slope::Float32           # self-thin slope (≤−1.605, ≥−10)        (SLPMSB; default 0 = MSB off)
    msb_eff::Float32             # mortality efficiency per pass          (EFFMSB; default 0.90)
    msb_dlo::Float32             # lower DBH (≥) of the kill range        (DLOMSB; default 0)
    msb_dhi::Float32             # upper DBH (<) of the kill range        (DHIMSB; default 999)
    msb_flag::Int32              # 1=from above, 2=from below, 3=throughout (MFLMSB; default 1)
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
    sp_bf_dbhmin::Vector{Float32}#  board-foot min merch DBH              (BFMIND)
    sp_bf_topd::Vector{Float32}  #  board-foot merch top diameter         (BFTOPD)
    sp_bf_stump::Vector{Float32} #  board-foot stump height               (BFSTMP)
    sp_cf_defect::Matrix{Float32}#  cubic-foot defect curve, 9 DBHCLS pts × MAXSP (CFDEFT; MCDEFECT)
    sp_bf_defect::Matrix{Float32}#  board-foot defect curve, 9 DBHCLS pts × MAXSP (BFDEFT; BFDEFECT)
    sp_cf_form0::Vector{Float32} #  cubic log-linear form model intercept (CFLA0; MCFDLN), default 0
    sp_cf_form1::Vector{Float32} #  cubic log-linear form model slope     (CFLA1; MCFDLN), default 1
    sp_bf_form0::Vector{Float32} #  board log-linear form model intercept (BFLA0; BFFDLN), default 0
    sp_bf_form1::Vector{Float32} #  board log-linear form model slope     (BFLA1; BFFDLN), default 1

    schedule::Vector{ScheduledActivity}  # parsed THIN*/harvest activities (cuts!)
    conditionals::Vector{ConditionalActivity} # IF/THEN/ENDIF event-monitor blocks
    years_cut::Set{Int32}                # years a thin has already been applied (idempotent cuts!)
    yardloss_prlost::Float32             # YARDLOSS PRLOST (cuts.f:1461): proportion of harvested merch/saw/
                                         # board volume lost in yarding (0 = inactive). The reported removed
                                         # merch/saw/bdft are scaled by (1−PRLOST); total cubic + BA are not.
    yardloss_prdsng::Float32             # YARDLOSS PRDSNG (cuts.f:1462): of the yarding LOSS, this proportion
                                         # becomes DOWNED snags and (1−PRDSNG) STANDING snags (FFE, cuts.f:1384-85).
    cut_pref::Vector{Int32}              # per-species cut preference (IORDER, set by SPECPREF)
    multipliers::Vector{GrowthMultiplier} # keyword growth/mortality multipliers (MULTS)
    htgstp_events::Vector{ScheduledActivity} # HTGSTOP/TOPKILL top-damage events (htgstp.f);
                                             # icflag = activity (110 HTGSTOP / 111 TOPKILL),
                                             # params = species,HT1,HT2,PRB,AVEPRB,STDPBR
    fixmort_events::Vector{ScheduledActivity} # FIXMORT forced-mortality events (morts.f:781);
                                             # params = species,rate,d1,d2,IP(1=replace/2=add/
                                             # 3=max/4=mult),pointflag
    sp_groups::Vector{Vector{Int32}}         # SPGROUP species groups (ISPGRP members); group N
                                             # = sp_groups[N], referenced in a species field by -N
    sp_group_names::Vector{String}           # SPGROUP group names (NAMGRP), upper-cased
    compute_defs::Vector{Tuple{Int32,String,EvNode}} # COMPUTE event-monitor variable definitions
                                             # (start year, NAME, parsed RHS expression AST)
    compute_vars::Dict{String,Float32}       # current COMPUTE variable values (EVMON user vars)
    volume_events::Vector{ScheduledActivity} # VOLUME/BFVOLUME merch-standard overrides (volkey.f);
                                             # icflag = activity (217 VOLUME cubic / 218 BFVOLUME bd-ft),
                                             # params = species, then the merch standards for that path
    merch_init::Bool                         # whether the per-stand sp_* merch arrays are populated yet
    voleqnum_overrides::Vector{Tuple{Int32,String}} # VOLEQNUM (initre.f:5061): (species idx, NVEL
                                             # cubic equation id) — overrides species.vol_eq after VOLEQDEF
    sp_bf_vol_eq::Vector{String}             # board-foot NVEL equation per species (VEQNNB) — snapshot
                                             # of vol_eq taken BEFORE VOLEQNUM, so board feet keeps the
                                             # default equation when only the cubic eq is overridden
    fertilize_events::Vector{ScheduledActivity} # FERTILIZE (ffin.f): icflag 260, params[1]=efficacy
    ifert_date::Int32                        # year of the last fertilizer application (FFCOM IFFDAT; −1=none)
    ifert_eff::Float32                       # efficacy of that application (FFCOM FFPRMS(4))
    lsprut::Bool                             # stump-sprouting enabled (LSPRUT; SPROUT/NOSPROUT)
    sprout_smult::Float32                    # sprout NUMBER multiplier (SMULT; SPROUT), default 1
    sprout_hmult::Float32                    # sprout HEIGHT multiplier (HMULT; SPROUT), default 1
    # Per-species SPROUT keyword table (esuckr.f activity 450): each entry = (species_code, smult, hmult,
    # dmin, dmax). species_code follows SPDECD: >0 single, 0 all-sproutable, <0 group −code. A parent stump
    # of species s with stump DBH in [dmin,dmax) takes that entry's smult/hmult (last match wins; default 1/1).
    sprout_overrides::Vector{NTuple{5,Float32}}
    cut_log::Vector{CutRecord}               # ESTUMP cut record per removed sprouting tree, fed to ESUCKR
    dg_stddev_bound::Float32                  # DGSD: std-dev bound on stochastic DG variation (DGSTDEV; <1 ⇒ off)
    dg_bjphi::Float32                         # ARMA(1,1) AR parameter for DGSCOR (BJPHI; SERLCORR), default 0.74
    dg_bjthet::Float32                        # ARMA(1,1) MA parameter for DGSCOR (BJTHET; SERLCORR), default 0.42
    age_reset_year::Int32                     # RESETAGE: calendar year of the age reset (−1 = none)   (resage.f)
    age_reset_age::Int32                      # RESETAGE: stand age to set AT that year                (PRMS(1))
    dbs_out_file::String                      # DATABASE/DSNOUT output SQLite file ("" = none)         (DSNOUT)
    dbs_summary::Bool                         # DATABASE SUMMARY ⇒ emit the FVS_Summary table          (ISUMARY)
    dbs_treelist::Bool                        # DATABASE TREELIDB ⇒ emit the FVS_TreeList table        (ITREELIST)
    dbs_compute::Bool                         # DATABASE COMPUTDB ⇒ emit the FVS_Compute table         (ICOMPUTE)
    dbs_cutlist::Bool                         # DATABASE CUTLIST ⇒ emit the FVS_CutList table          (ICUTLIST)
    cutlist_capture::Union{Nothing,Vector{Any}} # active per-cycle cut-record sink (_log_cut!), else nothing
    strclass_on::Bool                         # STRCLASS keyword ⇒ compute the structural stage each cycle (LCALC)
    strclass_thresh::NTuple{6,Float32}        # STRCLASS thresholds: gappct/ssdbh/sawdbh/ccmin/tpamin/pctsmx
    growth_idg::Int32                         # GROWTH: input DIAMETER-growth data type (0=none/incr, 1/3=past DBH, 2=incr) (IDG)
    growth_ihtg::Int32                        # GROWTH: input HEIGHT-growth data type (IHTG)
    growth_fint::Float32                      # GROWTH: DG measurement period (FINT, default 5)
    growth_finth::Float32                     # GROWTH: HTG measurement period (FINTH, default 5)
    growth_fintm::Float32                     # GROWTH: mortality measurement period (FINTM, default 5)
    cycle_lengths::Vector{Int32}              # per-cycle period override (TIMEINT field-1); 0 = use uniform `year`  (IY pre-cumulation)
    cycleat_years::Vector{Int32}              # CYCLEAT-requested extra cycle-boundary years (calendar)              (IWORK1)
    ncycle_eff::Int32                         # effective cycle count after CYCLEAT insertions (build_cycle_schedule!) (NCYC)
    dg_cor2::Vector{Float32}                  # READCORD/REUSCORD large-tree DG correction terms per sp (default 1)  (COR2)
    htg_cor2::Vector{Float32}                 # READCORH/REUSCORH large-tree HTG correction terms per sp (default 1) (HCOR2)
    regh_cor2::Vector{Float32}                # READCORR/REUSCORR small-tree HTG correction terms per sp (default 1) (RCOR2)
    dg_cor2_on::Bool                          # LDCOR2: apply ln(COR2) to DGCON before calibration                   (LDCOR2)
    htg_cor2_on::Bool                         # LHCOR2: apply ln(HCOR2) to HTCON before calibration                  (LHCOR2)
    regh_cor2_on::Bool                        # LRCOR2: small-tree con multiplied by RHCON=RCOR2                     (LRCOR2)
    carbon_report_on::Bool                    # CARBREPT: emit the Stand Carbon Report                              (LCARBON)
    carbon_method::Int32                      # CARBCALC field 1: 0 = FFE (default), 1 = JENKINS                    (ICMETH)
    carbon_units::Int32                       # CARBCALC field 2: 0 = US tons/acre (default, USA), 1 = metric t/ha,
                                              #                   2 = metric tons/acre (fminit.f:909-914 ICMETRC)
    potfire_report_on::Bool                   # POTFIRE: emit the Potential Fire (FMPOFL) report                    (IPFLMB/E)
    unrecognized_keywords::Set{String}        # keywords seen but neither dispatched nor a KNOWN_NOOP — surfaced
                                              # so a silently-ignored SN semantic (e.g. YARDLOSS) can't hide
end

function Control()
    s4(n)   = fill("    ", n)
    Control(
        true,                                                  # faithful
        " ", " ", "SN", "       ", repeat(' ',160), repeat(' ',250), repeat(' ',72),
        s4(MAXSP), fill(" "^10, 30), fill(" "^10, 30),
        false,false,false,false,false,false,false,false,false,false,false,false,
        false,false,false,false,false,false,false,false,false,
        trues(MAXSP), zeros(Bool,MAXSP), zeros(Bool,MAXSP),     # dg_calib_sp(LDGCAL)=on, leave_species, ht_drag_sp
        Int32(0),Int32(0),                                      # error_code, cut_algorithm
        Int32(0),Int32(0),Int32(0),Int32(2),Int32(0),Int32(0), # icl1..6 (icl4=tripling cycle limit, grinit ICL4=2)
        Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0),Int32(0), # cycle..start_year (9)
        Int32(0),Int32(0),Int32(0),Int32(0),                    # ithnpa,ithnpi,ithnpn,ntrees_active
        Int32(0),Int32(7),Int32(6),Int32(8),Int32(9),           # unit_calib,list,stand,summary,tree
        Int32(0),Int32(0),                                      # lstknt,nstknt
        Int32(0),Int32(0),Int32(0),Int32(0),                    # ncycle,n_ptgroups,n_spgroups,nspecies
        zeros(Int32,MAXSP), zeros(Int32,MAXSP), zeros(Int32,MAXSP,2),
        zeros(Int32,MAXSP), zeros(Int32,MAXSP), zeros(Int32,MAXSP), zeros(Int32,MAXSP),
        zeros(Int32,6), zeros(Int32,30,52), zeros(Int32,30,92),
        zeros(Int32,MAXCY1), zeros(Int32,7),
        0.0f0,60.0f0,45.0f0,0.0f0,0.0f0,0.0f0,1.0f0,1.0f0,0.0f0,  # …,cc_coef(CCCOEF)=1, cc_coef2(CCCOEF2)=1 (grinit.f:269-270)
        0.0f0,0.0f0,0.0f0,                                      # dbh_sdi, dbh_stage, dbh_zeide
        999f0, 0f0, 0.90f0, 0f0, 999f0, Int32(1),               # MORTMSB: msb_qmd/slope/eff/dlo/dhi/flag (off by default)
        0.0f0,0.0f0,0.0f0,0.0f0,1.0f0,5.0f0,                    # dg_sd, zeide_dr016/at/old, cut_eff(EFF)=1.0, mort_period=5
        0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,0.0f0,
        zeros(Float32,MAXSP), zeros(Float32,MAXSP), zeros(Float32,MAXSP),
        zeros(Float32,MAXSP,4), zeros(Float32,MAXSP), zeros(Float32,MAXSP),
        zeros(Float32,MAXSP), zeros(Float32,MAXSP), zeros(Float32,MAXSP),
        zeros(Float32,MAXSP), zeros(Float32,MAXSP), zeros(Float32,MAXSP), # sp_bf_dbhmin, sp_bf_topd, sp_bf_stump
        zeros(Float32,9,MAXSP), zeros(Float32,9,MAXSP),         # sp_cf_defect (MCDEFECT), sp_bf_defect (BFDEFECT)
        zeros(Float32,MAXSP), ones(Float32,MAXSP),              # sp_cf_form0/1 (MCFDLN): CFLA0=0, CFLA1=1
        zeros(Float32,MAXSP), ones(Float32,MAXSP),              # sp_bf_form0/1 (BFFDLN): BFLA0=0, BFLA1=1
        ScheduledActivity[], ConditionalActivity[], Set{Int32}(), # schedule, conditionals, years_cut
        0f0,                                                    # yardloss_prlost (YARDLOSS, inactive)
        0f0,                                                    # yardloss_prdsng (YARDLOSS downed-snag proportion)
        zeros(Int32, MAXSP),                                    # cut_pref (IORDER)
        GrowthMultiplier[],                                     # multipliers (MULTS)
        ScheduledActivity[],                                    # htgstp_events (HTGSTOP/TOPKILL)
        ScheduledActivity[],                                    # fixmort_events (FIXMORT)
        Vector{Int32}[], String[],                              # sp_groups, sp_group_names (SPGROUP)
        Tuple{Int32,String,EvNode}[], Dict{String,Float32}(),   # compute_defs, compute_vars (COMPUTE)
        ScheduledActivity[], false,                             # volume_events (VOLUME/BFVOLUME), merch_init
        Tuple{Int32,String}[],                                  # voleqnum_overrides (VOLEQNUM)
        String[],                                               # sp_bf_vol_eq (board equation snapshot)
        ScheduledActivity[], Int32(-1), 0f0,                    # fertilize_events, ifert_date, ifert_eff
        true, 1f0, 1f0, NTuple{5,Float32}[], CutRecord[],       # lsprut (FVS esinit.f:50 default ON; NOAUTOES/NOSPROUT turn it off), sprout_smult, sprout_hmult, sprout_overrides, cut_log
        2f0, 0.74f0, 0.42f0,                                    # dg_stddev_bound(DGSD=2), dg_bjphi(0.74), dg_bjthet(0.42)
        Int32(-1), Int32(0),                                    # age_reset_year(none), age_reset_age
        "", false, false, false,                                # dbs_out_file, dbs_summary, dbs_treelist, dbs_compute (DATABASE)
        false, nothing,                                         # dbs_cutlist, cutlist_capture (FVS_CutList)
        false, SS_THRESH_DEFAULT,                               # strclass_on, strclass_thresh (SSTAGE)
        Int32(0), Int32(0), 5f0, 5f0, 5f0,                      # GROWTH: idg, ihtg, fint, finth, fintm (defaults)
        zeros(Int32, MAXCY1), Int32[], Int32(0),                 # cycle_lengths(TIMEINT), cycleat_years(CYCLEAT), ncycle_eff
        ones(Float32, MAXSP), ones(Float32, MAXSP), ones(Float32, MAXSP),  # dg_cor2/htg_cor2/regh_cor2 (COR2/HCOR2/RCOR2 = 1)
        false, false, false,                                     # dg_cor2_on/htg_cor2_on/regh_cor2_on (LDCOR2/LHCOR2/LRCOR2)
        false, Int32(0), Int32(0),                               # carbon_report_on, carbon_method=FFE, carbon_units=US-t/ac (FVS fminit.f defaults)
        false,                                                   # potfire_report_on (POTFIRE)
        Set{String}(),                                           # unrecognized_keywords
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
    htg_cor_init::Vector{Float32}  # initial small-tree HCOR from the regent calibration regression (HCOR @ ICYC=1)
    bark_a::Vector{Float32}      # per-stand bark intercept (BRATIO; Fort Bragg override)
    bark_b::Vector{Float32}      # per-stand bark slope     (BRATIO; Fort Bragg override)
    ht_dbh_aa::Vector{Float32}   # calibrated Wykoff HT-DBH intercept per sp (NOHTDREG/LHTDRG, cratet.f:329) (AA)
    ht_dbh_iabflg::Vector{Int32} # 0=use calibrated AA (Wykoff), 1=use inventory Curtis-Arney HTDBH (cratet.f) (IABFLG)
    vmlt::Float32                # ARMA variance multiplier (calibration)  (VMLT)
    calib_dbh::Vector{Float32}   # transient: CURRENT-stand dbh for the NE calibration BADIST (empty except
                                 # during calibrate's dgf! — FVS NE calib computes BAL on the current stand while
                                 # predicting at the backdated per-tree dbh; ne_badist! reads this when non-empty)
end
Calibration() = Calibration(ones(Float32,MAXSP), ones(Float32,MAXSP),
    zeros(Float32,MAXSP), zeros(Float32,MAXSP), zeros(Float32,MAXSP),
    zeros(Float32,MAXSP), zeros(Float32,MAXSP), zeros(Float32,MAXSP),
    zeros(Float32,MAXSP), zeros(Float32,MAXSP), zeros(Float32,MAXSP),
    zeros(Float32,MAXSP), zeros(Float32,MAXSP),
    zeros(Float32,MAXSP), ones(Int32,MAXSP), 0f0, Float32[])  # ht_dbh_aa=0, iabflg=1, calib_dbh empty

# ---------------------------------------------------------------------------
# Density — COMMON /PDEN/ : stand density / SDI scratch (C4). Minimal for now.
# ---------------------------------------------------------------------------
mutable struct Density
    sdi_sum::Float32
    point_ba::Vector{Float32}    # per-point basal area (PTBAA), indexed by subplot
    point_bal::Vector{Float32}   # per-tree BA in larger trees on its point (PTBALT)
    point_ccf::Vector{Float32}   # per-point crown competition factor (PCCF), indexed by subplot (dense.f:210)
    point_tpa::Vector{Float32}   # per-point trees-per-acre (PTPA), indexed by subplot           (dense.f:211)
    mort_slope::Float32          # Pretzsch self-thinning line slope             (SLPMRT)
    mort_intercept::Float32      # Pretzsch self-thinning line intercept         (CEPMRT)
    tpa_mort::Float32            # last cycle's surviving over-threshold TPA      (TPAMRT)
end
Density() = Density(0.0f0, zeros(Float32, MAXPLT), zeros(Float32, MAXTRE),
                    zeros(Float32, MAXPLT), zeros(Float32, MAXPLT), 0.0f0, 0.0f0, 0.0f0)

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
    # Per-cycle mortality work buffers (VARMRT) — preallocated to MAXTRE and reused each cycle (sliced to
    # the live count), so the hot mortality path allocates nothing. `killed` is zeroed per call, `efftr`
    # is fully rewritten, `temwk2` is written-before-read per index, so reuse is value-safe.
    mort_killed::Vector{Float32}
    mort_efftr::Vector{Float32}
    mort_temwk2::Vector{Float32}
end
Scratch() = Scratch(zeros(Float32,15,MAXTRE), zeros(Int32,MAXTRE), zeros(Int32,MAXTRE), zeros(Int32,MAXTRE),
                    zeros(Float32,MAXTRE), zeros(Float32,MAXTRE), zeros(Float32,MAXTRE))

# ---------------------------------------------------------------------------
# Extension states — allocated lazily only when the extension is active.
# Filled out in their own chunks (fire C7, econ C8, establishment C4).
# ---------------------------------------------------------------------------
mutable struct Establishment
    active::Bool
    idsdat::Int32       # date of disturbance (ESTAB keyword); -9999 = unset (ESNUTR defaults it)
    ntally::Int32       # regen-tally counter (NTALLY)
    es_seed::Float32    # persisted ESDRAW establishment reseed value (estab.f ESDRAW, SAVEd across tallies)
    years_done::Set{Int32}  # establishment years already applied (idempotent ESNUTR)
end
Establishment() = Establishment(false, Int32(-9999), Int32(0), 0f0, Set{Int32}())

mutable struct DbsState
    enabled::Bool
    db_path::String
    out_db::Any                  # SQLite.DB handle (or nothing); not serialized
    report_flags::Dict{Symbol,Bool}
end
DbsState() = DbsState(false, "FVSOut.db", nothing, Dict{Symbol,Bool}())

"""
Standing-dead (snag) records (FFE FMCOM), structure-of-arrays. Each entry is a cohort of
snags created when trees died: `sp` species, `dbh`, `den_hard`/`den_soft` the densities
(stems/ac) still standing in the hard / soft decay states, `origden` the density at
creation, `year` the year of death. Snags fall and decay each cycle (`update_snags!`).
"""
mutable struct SnagList
    sp::Vector{Int32}
    dbh::Vector{Float32}
    den_hard::Vector{Float32}     # DENIH — initially-hard snags still standing
    den_soft::Vector{Float32}     # DENIS — soft (decayed) snags still standing
    origden::Vector{Float32}      # DEND  — original density at creation
    year::Vector{Int32}           # the FALL-clock start (jl's age-based update_snags! basis); for ordinary
                                  # mortality this is the cycle-START year (tuned so the bit-exact StandDead
                                  # falldown holds), which is NOT FVS's YRDEAD — see `yrdead`.
    yrdead::Vector{Int32}         # the TRUE death year = FVS YRDEAD (input: dead-10yr; fire: fire year; ordinary
                                  # mortality: cycle-END−1, fmkill.f:140 IY(ICYC+1)−1). Used for the hard→soft
                                  # DKTIME classification (snag_summary) + the post-burn window — NOT the fall.
    bolevol::Vector{Float32}      # per-tree death-time STEM-volume bole biomass, tons (cuft·V2T) — the
                                  # FFE snag basis (SNVIS·V2T), distinct from whole-tree Jenkins; 0 = unset
    height::Vector{Float32}       # HTDEAD — snag height at death (ft); drives the cone-taper split of a
                                  # fallen bole across CWD size classes (FMCWD/CWD1). 0 = unset (single-class)
    htcur::Vector{Float32}        # HTIH/HTIS — CURRENT snag height (ft), = `height` at creation; only shrinks
                                  # when SNAGBRK sets HTX>0 (FMSNGHT). Drives the recomputed bole then. At
                                  # default (HTX=0) it stays = `height`, so the frozen `bolevol` is used (bit-exact).
end
SnagList() = SnagList(Int32[], Float32[], Float32[], Float32[], Float32[], Int32[], Int32[], Float32[], Float32[], Float32[])

"""
PotFire-report weather-scenario conditions, overridable by the POTF* keywords (POTFMOIS/POTFWIND/
POTFTEMP/POTFSEAS/POTFPAB) for the SEVERE and MODERATE scenarios of the FVS_PotFire report (fmpofl.f).
A field of −1 (moisture/wind/temp/pab) or 0 (season) means "use the hardcoded scenario default".
"""
mutable struct PotFireCond
    mois::NTuple{7,Float32}   # 7 fuel-moisture fractions (1hr/10hr/100hr/3+/duff/woody/herb); −1 = default
    wind::Float32             # 20-ft wind (mi/h); −1 = default
    temp::Float32             # air temperature (°F); −1 = default
    season::Int32             # burn season 1-4; 0 = default
    pab::Float32              # percent area burned; −1 = default
end
PotFireCond() = PotFireCond(ntuple(_ -> -1f0, 7), -1f0, -1f0, Int32(0), -1f0)

"""
Overridable FFE (Fire & Fuels Extension) model parameters — the scalar coefficients that
FMIN-block keywords change from their variant defaults. Defaults are the SN values (fmvinit.f).
Each keyword that tunes a model coefficient (SNAGPBN, FUELMULT, MOISTURE, …) sets the matching
field here; the apply sites read these instead of hardcoded constants. Grows as more FFE keywords
are ported. Kept on FireState (per-stand, no globals).
"""
mutable struct FFEParams
    # post-burn accelerated snag fall (SNAGPBN; fmvinit.f:1100-1104, applied in FMSNAG/FMSFALL)
    pb_soft::Float32   # PBSOFT: total fraction of soft-at-fire snags fallen after PBTIME (1.0)
    pb_smal::Float32   # PBSMAL: total fraction of small (<PBSIZE) snags fallen after PBTIME (0.9)
    pb_size::Float32   # PBSIZE: small-snag DBH breakpoint, in (12.0)
    pb_time::Float32   # PBTIME: post-burn accelerated-fall window, yr (7.0)
    pb_scor::Float32   # PBSCOR: scorch-height threshold (ft) to trigger post-burn fall (fmburn.f:414) (0.0)
    # SNAGFALL (opt 9): per-species overrides of the snag fall-rate correction (FALLX) and the snag age by
    # which the last 5% fall (ALLDWN). Sparse — only overridden species; snag_fall_density prefers these
    # over the fire_species_props.csv defaults. Keyed by species index.
    snag_fallx_ovr::Dict{Int32,Float32}
    snag_alldwn_ovr::Dict{Int32,Float32}
    # FUELMULT (opt 29) / FUELDCAY (opt 16): override of the [11 size, 4 decay-class] fuel decay-rate matrix
    # DKR. Empty (0×0) ⇒ use the `_FM_DKR` default; the keywords lazily copy _FM_DKR then modify it, and
    # fmcwd! reads this when populated.
    dkr::Matrix{Float32}
    # FUELINIT (opt 21) / FUELSOFT (opt 53): per-size-class (1:11) overrides of the initial hard / soft
    # surface-fuel loading (STFUEL, tons/ac). Empty ⇒ defaults (hard = ffe_dead_fuel_loading, soft = 0);
    # a value ≥ 0 in a slot overrides that size class. fmcba! applies these at the first-FFE-year fuel load.
    stfuel_hard::Vector{Float32}
    stfuel_soft::Vector{Float32}
    # DUFFPROD (opt 17): override of the [11 size, 4 decay-class] proportion-of-decay-to-duff matrix PRDUFF.
    # Empty (0×0) ⇒ the uniform 0.02 default (`_FM_PRDUFF`); fmcwd! reads this when populated.
    prduff::Matrix{Float32}
    # FUELPOOL (opt 19): per-species override of the fuel decay-rate class DKRCLS (1-4) — which decay-class
    # column a species' dead fuel / snag bole flows into. Sparse; `ffe_dkr_cls` prefers it over the CSV.
    dkrcls_ovr::Dict{Int32,Int32}
    # POTF* (opt 30/35/36/41/42): PotFire-report SEVERE (idx 1) and MODERATE (idx 2) scenario conditions.
    potf::NTuple{2,PotFireCond}
    # SNAGPSFT (opt 37): per-species proportion of snags that are SOFT at creation (PSOFT, default 0 = all
    # hard). Sparse; add_snag! splits new density into hard/soft by this fraction.
    psoft_ovr::Dict{Int32,Float32}
    # SNAGDCAY (opt 11): per-species DECAYX override (snag decay-rate multiplier; DKTIME/TSOFT =
    # DECAYX·(1.24·D+13.82)). Sparse — only overridden species; the snag soft-decay transition + the
    # crown-fall TSOFT prefer these over the fire_species_props.csv defaults (0.07/0.21/0.35).
    snag_decayx_ovr::Dict{Int32,Float32}
    # SNAGBRK (opt 10): per-species snag height-LOSS coefficients HTX(1..4) = (hard>0.5HTD, hard<0.5HTD,
    # soft>0.5HTD, soft<0.5HTD), calibrated from the keyword's YRS50/YRS30 (fmin.f:538/546/557/566). Sparse;
    # EMPTY ⇒ HTX=0 = no height loss (the SN default, fmvinit.f:1089) ⇒ snag bole stays frozen at death height.
    # When set, FMSNGHT shrinks the snag height each year (snag.jl) and the bole volume is recomputed from it.
    snag_htx::Dict{Int32,NTuple{4,Float32}}
end
FFEParams() = FFEParams(1.0f0, 0.9f0, 12.0f0, 7.0f0, 0.0f0, Dict{Int32,Float32}(), Dict{Int32,Float32}(),
                        Matrix{Float32}(undef, 0, 0), Float32[], Float32[], Matrix{Float32}(undef, 0, 0),
                        Dict{Int32,Int32}(), (PotFireCond(), PotFireCond()), Dict{Int32,Float32}(),
                        Dict{Int32,Float32}(), Dict{Int32,NTuple{4,Float32}}())

mutable struct FireState
    active::Bool                       # FFE enabled (FMIN keyword)
    covtyp::Int32                      # cover type = species with the most basal area (COVTYP)
    percov::Float32                    # percent canopy cover (PERCOV)
    bigdbh::Float32                    # largest DBH seen in the stand (BIGDBH)
    flive::NTuple{2,Float32}           # live herb / shrub surface fuel, tons/ac (FLIVE)
    cwd::Array{Float32,3}              # surface fuel [size class 1:11, dead(2)/soft(1), decay class 1:4] (CWD)
    fuels_init::Bool                   # dead-fuel pools loaded yet (first FFE year only)
    # scheduled SIMFIRE event (fire_year = 0 ⇒ none) and its conditions
    fire_year::Int32                   # calendar year of the simulated fire (SIMFIRE date)
    swind::Float32                     # 20-ft wind (mi/h)               (SWIND, PRMS1)
    fmois::Int32                       # fuel-moisture dryness model 1–4 (FMOIS, PRMS2)
    atemp::Float32                     # air temperature (°F)            (ATEMP, PRMS3)
    mortcode::Int32                    # 1 = FFE estimates mortality     (MKODE, PRMS4)
    psburn::Float32                    # percent of the stand burned     (PSBURN, PRMS5)
    burnseas::Int32                    # burn season 1–4                 (BURNSEAS, PRMS6)
    flmult::Float32                    # flame-length multiplier         (FLAMEADJ)
    crburn::Float32                    # crown-fire fraction             (FLAMEADJ)
    snags::SnagList                    # standing-dead snag cohorts
    bioroot::Float32                   # dead coarse-root biomass pool, tons/ac (BIOROOT) — accrues at
                                       # tree death, decays at CRDCAY each cycle; → the Below-Dead column
    cwd2b::Array{Float32,3}            # crown debris-in-waiting [decay 1:4, crown size 0:5→idx 1:6,
                                       # year-to-fall 1:60] (CWD2B); the un-fallen part is the Stand-Dead
                                       # crown, it flows to `cwd` (down wood) as it falls (FMSCRO)
    crown_lift_annual::Matrix{Float32} # per-YEAR crown-lift down-wood input [cwd size 1:9, decay 1:4],
                                       # computed once per cycle (FMSDIT OLDCRW) and added each year in the
                                       # fuel loop (FMCADD crown-lift term); zero unless compute_crown_lift!
                                       # has run this cycle (so the non-carbon-report path is unaffected)
    burn_reports::Vector{Any}          # one record per SIMFIRE event (year + moistures + wind + flame +
                                       # scorch + fuel models/weights + consumption) for the FVS_BurnReport /
                                       # FVS_Consumption / FVS_Mortality DBS tables (pushed by fmburn!)
    hwp_fate::Dict{NTuple{3,Int},Float32}  # harvested-wood-products carbon FATE accumulator (fmscut.f:151):
                                       # key (cut_year, product 1=pulp/2=saw, group 1=sw/2=hw) → harvested
                                       # merch BIOMASS (tons/ac). Drives the FVS_Hrv_Carbon table via the
                                       # FAPROP year-since-harvest decay curves. State per stand (no globals).
    fire_smlg::NTuple{2,Float32}       # the (SMALL,LARGE) down-wood an actual SIMFIRE burns on this cycle:
                                       # start-of-cycle + the fire-year's single annual fuel step (FVS runs the
                                       # surface-fuel loop ANNUALLY interleaved with FMBURN, so the fire sees
                                       # cycle-start + 1 yr, not the period-end left by ffe_fuel_update!). (-1,-1)=unset.
    params::FFEParams                  # overridable FFE model coefficients (SNAGPBN/FUELMULT/… keyword tuning)
    moisture_ovr::Vector{Tuple{Int32,NTuple{7,Float32}}}  # MOISTURE keyword: (date, 7 fuel-moisture % —
                                       # 1hr/10hr/100hr/3+/duff/live-woody/live-herb). A fire in the matching
                                       # cycle uses these instead of the FMMOIS dryness-model table (fmburn.f:367).
    snaginit::Vector{NTuple{5,Float32}}  # SNAGINIT keyword: user-added snags (species, DBH-at-death, ht-at-death,
                                       # age, density stems/ac). Seeded at the first FFE year (fmsnag.f:90-105).
    salv_isalvs::Int32                 # SALVSP: salvage species selector (0=all, >0 species, <0 −SPGROUP); persists
    salv_isalvc::Int32                 # SALVSP: 0=cut-list (cut only ISALVS) / 1=leave-list (leave ISALVS, cut rest)
    fuelmodl::Vector{Tuple{Int32,Vector{Tuple{Int32,Float32}}}}  # FUELMODL: (date, [(standard-model#, weight)])
                                       # forced fuel models — used in place of FMCFMD auto-selection at the matching cycle
    fueltret::Vector{Tuple{Int32,Float32}}  # FUELTRET: (date, DPMOD) — fuel-bed DEPTH multiplier from the harvest/
                                       # treatment type, applied to the selected fuel model for ~5 yr after the date (fmusrfm.f)
    defulmod::Dict{Int32,Tuple{Matrix{Float32},Matrix{Float32},Float32,Float32}}  # DEFULMOD: model# → overridden
                                       # (load[2,4], sav[2,4], depth, mext); standard_fuel_model returns this in place of the table
    fire_schedule::Vector{NTuple{7,Float32}}  # every scheduled SIMFIRE as a separate event (FVS OPNEW: one
                                       # activity per keyword) — (year, swind, fmois, atemp, mortcode, psburn,
                                       # burnseas), conditions resolved with FVS defaults at parse time. A fire is
                                       # popped + its conditions loaded into the scalars above when its year falls in
                                       # the current cycle (so >1 SIMFIRE, e.g. fire_repeat, each fire at its own date).
end
FireState() = FireState(false, Int32(0), 0f0, 0f0, (0f0, 0f0), zeros(Float32, 11, 2, 4), false,
                        Int32(0), 20f0, Int32(1), 70f0, Int32(1), 100f0, Int32(1), 1f0, 0f0, SnagList(), 0f0,
                        zeros(Float32, 4, 6, 60), zeros(Float32, 9, 4), Any[], Dict{NTuple{3,Int},Float32}(),
                        (-1f0, -1f0), FFEParams(), Tuple{Int32,NTuple{7,Float32}}[], NTuple{5,Float32}[],
                        Int32(0), Int32(0), Tuple{Int32,Vector{Tuple{Int32,Float32}}}[],
                        Tuple{Int32,Float32}[],
                        Dict{Int32,Tuple{Matrix{Float32},Matrix{Float32},Float32,Float32}}(),
                        NTuple{7,Float32}[])

"""
One ECON harvest cost or revenue record (HRVVRCST / HRVRVN): `amount` per `unit`,
applied to harvested trees with DBH in `[dbh_lo, dbh_hi)`. `sp` is the species the
record applies to (0 = all). Unit codes (ECNCOM.F77): 1=per tree, 2=per MBF, 3=per CCF.
"""
struct EconCostRev
    amount::Float32
    unit::Int32
    dbh_lo::Float32
    dbh_hi::Float32
    sp::Int32
end
EconCostRev(amount, unit, lo, hi) = EconCostRev(Float32(amount), Int32(unit), Float32(lo), Float32(hi), Int32(0))

"ECON economic-analysis state (no globals): discount rate, cost/revenue keyword tables, accumulated streams."
mutable struct EconState
    active::Bool
    discount_rate::Float32                # discount/interest rate (fraction); STRTECON field 2 / 100 (default 0, ecinit.f:15)
    ann_cost::Float32                     # ANNUCST total annual management cost ($/ac/yr)
    hrv_cost::Vector{EconCostRev}         # HRVVRCST variable harvest costs
    hrv_rev::Vector{EconCostRev}          # HRVRVN harvest revenues (per species)
    base_year::Int32                      # ECON analysis start year (−1 until set)
    harvests::Vector{NTuple{3,Float32}}   # accumulated (year, cost, revenue) per harvest
    cycle_cost::Float32                   # this cycle's harvest cost so far (accrued per cut tree)
    cycle_rev::Float32                    # this cycle's harvest revenue so far
    # log-graded HRVRVN (units 4/5) — FVS_EconHarvestValue report (echarv.f). REPORT-ONLY (does not feed PNV).
    tree_log_bf::Dict{Int,Dict{Int,Float32}}        # per-tree per-log-DIB-class BF (idib=>bf), refreshed each
                                                    # compute_volumes! when a unit-4/5 HRVRVN record exists; else empty
    log_grade_rev::Dict{NTuple{3,Int32},Float32}    # (year, speciesIdx, dibClass_x10) => Σ removed board feet
                                                    # (×price/1000 at emit time = the FVS_EconHarvestValue rows)
    tree_log_ft3::Dict{Int,Dict{Int,Float32}}       # unit-5 (FT3_100_LOG) analog of tree_log_bf: per-tree per-log-DIB
                                                    # gross cubic feet (R9LGCFT Smalian, renormalized to VOL(4)+VOL(7))
    log_grade_ft3::Dict{NTuple{3,Int32},Float32}    # (year, speciesIdx, dibClass_x10) => Σ removed cubic feet
                                                    # (×price/100 at emit time = the cubic FVS_EconHarvestValue rows)
end
EconState() = EconState(false, 0.0f0, 0f0, EconCostRev[], EconCostRev[], Int32(-1),
                        NTuple{3,Float32}[], 0f0, 0f0,
                        Dict{Int,Dict{Int,Float32}}(), Dict{NTuple{3,Int32},Float32}(),
                        Dict{Int,Dict{Int,Float32}}(), Dict{NTuple{3,Int32},Float32}())

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
