# =============================================================================
# Western Root Disease (WRD) extension — port of FVS rd/*.f
# =============================================================================
#
# Scope of THIS file (Chunk −1, infrastructure + keyword reader + inert seam):
#   * RD parameter constants        (rd/RDPARM.F77)
#   * RootDiseaseState               (the reduced-path members of /RDCOM/ /RDADD/
#                                     /RDCRY/ /RDARRY/ — commons ported as Julia state)
#   * rd_init_defaults!              (rd/rdinit.f — the scalar / small-array defaults)
#   * kw_rdin!                       (rd/rdin.f entry RDKEY — the RRTYPE / RRINIT /
#                                     SAREA / BBCLEAR / END subset)
#   * rd_active                      (rd/rdatv.f gate  L = RRTINV .OR. RRMAN)
#   * root_disease_setup! / _mn2! / _treg!  — the fvs.f/grincr.f/gradd.f engine seams,
#                                     wired GATED + INERT (early-return unless active);
#                                     the per-cycle MORTALITY body is Chunk 0 (pending).
#
# NOT in this file (Chunk 0 and later — deliberately deferred, see the port report):
#   * rdblk1.f host-species coefficient tables (HABFAC/PKILLS/PNINF/PCOLO/DECFN/…)
#   * rdsetp/rdinoc/rdarea/rdinf/rdmort/rdend/rdgrow mortality + growth-loss bodies
#   * the RD RNG (rdrani/rdrann/rdranp) + random center placement (rdcloc)
# These require a bit-exact FVS growth baseline for validation, which the KT engine
# does not yet provide for the reference stand (a cornered growth straddle exceeds
# the RD mortality signal). See docs and the port report for the measured boundary.
#
# All references are to /workspace/ForestVegetationSimulator/rd/.

# --- rd/RDPARM.F77 PARAMETERs -------------------------------------------------
const RD_DSO    = 1        # dead standing outside center
const RD_DSII   = 2        # dead standing inside, infected
const RD_DSIU   = 3        # dead standing inside, uninfected
const RD_IRRTRE = 1500     # max FVS tree records RD dimensions for
const RD_ITOTSP = 40       # total species across all WRD variants
const RD_ITOTRR = 4        # number of modeled root diseases

# Disease index → name (rd/rdin.f RRTYPE echo; 1=P-annosus 2=S-annosus 3=Armillaria
# 4=Phellinus). RRTYPE 3 in the turnkey oracle echoes "ARMILLARIA".
const RD_DISEASE_NAMES = ("ANNOSUS-P", "ANNOSUS-S", "ARMILLARIA", "PHELLINUS")
const RD_DISEASE_CODES = ('P', 'S', 'A', 'W')   # rd/rdpr.f CHTYPE

# AbstractRootDiseaseState is forward-declared in core/state.jl (StandState field type).

"""
    RootDiseaseState

Julia port of the FVS Western Root Disease common-block state, reduced to the
members exercised by the center-init (RRINIT) reduced mortality path. Populated by
`rd_init_defaults!` then `kw_rdin!`. `active` mirrors the rd/rdatv.f gate.
"""
mutable struct RootDiseaseState <: AbstractRootDiseaseState
    # ---- activation / gate (rd/RDADD.F77 /RDADD/, rd/rdatv.f) ----
    iroot::Int32              # RD model on flag (rdinit 0 → RDIN sets 1)
    rrman::Bool               # manual (keyword) init used  → RRINIT
    rrtinv::Bool              # treelist init used          → RRTREIN (not in chunk 0)
    lrtype::Bool              # RRTYPE keyword was given
    lbbon::Bool               # bark beetles default-on unless BBCLEAR (rdin END handler)
    bbclear::Bool             # BBCLEAR seen → suppress default bark beetles

    # ---- disease selection (rd/rdin.f RRTYPE) ----
    minrr::Int32              # lowest active disease number
    maxrr::Int32              # highest active disease number

    # ---- stand geometry (rd/RDCOM.F77) ----
    sarea::Float32            # stand area, acres
    dimen::Float32            # stand square dimension, ft = sqrt(SAREA)*208.7
    irhab::Int32              # habitat flag (1 default; 2 after IRGEN(8) year)
    istep::Int32              # time-slot counter (= ICYC+1); advanced in rdmn2

    # ---- per-disease init parameters (dim RD_ITOTRR) ----
    ncents::Vector{Int32}     # number of disease centers
    parea::Vector{Float32}    # diseased patch area, acres
    ooarea::Vector{Float32}   # previous-cycle PAREA
    prkill::Vector{Float32}   # infected trees/acre in diseased area (RRINIT fld3)
    prun::Vector{Float32}     # uninfected trees/acre in diseased area (RRINIT fld4)
    rrincs::Vector{Float32}   # initial proportion of roots infected (RRINIT fld5)
    lparea::Vector{Bool}      # PAREA set by keyword (do not default to 0.25·SAREA)
    lonect::Vector{Int32}     # 0 unassigned / 1 stand-is-one-center / 2 multiplot
    ipcflg::Vector{Int32}     # 0 random centers / 1 centers read explicitly

    # ---- mortality-curve constants (rd/rdinit.f; RRTYPE-independent) ----
    xxinf::NTuple{5,Float32}  # DBH knots of the years-to-kill curve (TEMP16)
    yyinf::NTuple{5,Float32}  # years-to-kill knots (TEMP17)
    nninf::Int32              # number of active knots (=3)
    stcut::NTuple{5,Float32}  # stump size-class breakpoints
    xminkl::Vector{Float32}   # min time-to-kill offset per disease (=0)
    xminlf::Vector{Float32}   # min inoculum lifespan per disease ({1,1,1,20})

    RootDiseaseState() = rd_init_defaults!(new())
end

"""
    rd_init_defaults!(rd) -> rd

Port of rd/rdinit.f: the RD common-block default initialization (the members the
reduced center-init mortality path reads). Values verified against the live FVSkt
oracle's RDIN option echo.
"""
function rd_init_defaults!(rd::RootDiseaseState)
    rd.iroot  = Int32(0)
    rd.rrman  = false
    rd.rrtinv = false
    rd.lrtype = false
    rd.lbbon  = true
    rd.bbclear = false

    rd.minrr = Int32(1)       # rdinit MINRR=1, MAXRR=2 (annosus) — overwritten by RRTYPE
    rd.maxrr = Int32(2)

    rd.sarea = 100.0f0
    rd.dimen = 2087.0f0       # sqrt(100)*208.7
    rd.irhab = Int32(1)
    rd.istep = Int32(1)

    n = RD_ITOTRR
    rd.ncents = fill(Int32(20), n)      # rdinit NCENTS=20
    rd.parea  = fill(25.0f0, n)         # rdinit PAREA=25
    rd.ooarea = zeros(Float32, n)
    rd.prkill = fill(0.5f0, n)          # rdinit PRKILL=0.5
    rd.prun   = fill(0.5f0, n)          # rdinit PRUN=0.5
    rd.rrincs = fill(0.1f0, n)          # rdinit RRINCS=0.1
    rd.lparea = fill(false, n)
    rd.lonect = fill(Int32(0), n)
    rd.ipcflg = fill(Int32(0), n)

    rd.xxinf = (0.0f0, 3.9f0, 35.4f0, 0.0f0, 0.0f0)     # TEMP16
    rd.yyinf = (0.0f0, 5.0f0, 40.0f0, 0.0f0, 0.0f0)     # TEMP17
    rd.nninf = Int32(3)
    rd.stcut = (0.0f0, 12.0f0, 24.0f0, 48.0f0, 100.0f0)
    rd.xminkl = fill(0.0f0, n)
    rd.xminlf = Float32[1.0, 1.0, 1.0, 20.0]            # TEMP15
    return rd
end

"""
    rd_active(rd) -> Bool

Port of rd/rdatv.f: L = RRTINV .OR. RRMAN. True once RD has been activated by an
RRINIT (manual) or RRTREIN (treelist) keyword.
"""
rd_active(rd::RootDiseaseState) = rd.rrtinv || rd.rrman
rd_active(::Nothing) = false

# -----------------------------------------------------------------------------
# kw_rdin! — RDIN keyword block reader (rd/rdin.f entry RDKEY subset)
# -----------------------------------------------------------------------------
"""
    kw_rdin!(s, rec, kr)

Parse the `RDin … End` block (rd/rdin.f). Chunk-0 subset: RRTYPE, RRINIT, SAREA,
BBCLEAR, END (RRDOUT/BBOUT recognized no-ops). Populates `s.root_disease`.
"""
function kw_rdin!(s::StandState, rec, kr::KeywordReader)
    s.root_disease === nothing && (s.root_disease = RootDiseaseState())
    rd = s.root_disease
    rd.iroot = Int32(1)                     # rdin.f: entry sets IROOT=1

    while true
        r = read_keyword!(kr)
        (r.status == KW_EOF || r.status == KW_STOP) && break
        k = strip(r.name)
        isempty(k) && continue
        if k == "END"
            rd_in_end!(rd)                  # rd/rdin.f option 9 (END)
            break
        elseif k == "RRTYPE"                # rd/rdin.f option 42
            rd_in_rrtype!(rd, r)
        elseif k == "RRINIT"                # rd/rdin.f option 5
            rd_in_rrinit!(rd, r)
        elseif k == "SAREA"                 # rd/rdin.f option 10
            rd_in_sarea!(rd, r)
        elseif k == "BBCLEAR"               # rd/rdin.f option 41 — suppress default bark beetles
            rd.bbclear = true; rd.lbbon = false
        elseif k == "RRDOUT" || k == "BBOUT" || k == "RRECHO" || k == "SMCOUT"
            # report-request sub-keywords — recognized, no state effect for chunk 0
        else
            # RD sub-keyword not yet ported (SPREAD/CARRY/PSTUMP/RRTREIN/…); record it
            (!isempty(k) && isletter(first(k))) && push!(s.control.unrecognized_keywords, k)
        end
    end
    return nothing
end

# rd/rdin.f option 42 — RRTYPE: pick the active disease number(s) from ARRAY.
function rd_in_rrtype!(rd::RootDiseaseState, r)
    rd.lrtype = true
    minrr = 0; maxrr = 0
    @inbounds for idi in 1:min(RD_ITOTRR, length(r.values))
        if r.present[idi]
            ipoint = Int(nint(r.values[idi]))
            (ipoint < minrr || minrr == 0) && (minrr = ipoint)
            ipoint > maxrr && (maxrr = ipoint)
        end
    end
    minrr <= 0 && (minrr = 1)
    minrr = clamp(minrr, 1, RD_ITOTRR)
    maxrr = clamp(maxrr, 1, RD_ITOTRR)
    maxrr < minrr && (maxrr = minrr)
    # only one disease type modeled at a time, except the P+S annosus pair (1,2)
    if maxrr != minrr && !(minrr == 1 && maxrr == 2)
        maxrr = minrr
    end
    rd.minrr = Int32(minrr)
    rd.maxrr = Int32(maxrr)
    return nothing
end

# rd/rdin.f option 5 — RRINIT: manual center-init. Fields (ARRAY):
#   1 = center-placement flag (<1 ⇒ random), 2 = #centers, 3 = infected TPA/acre,
#   4 = uninfected TPA/acre, 5 = proportion roots infected, 6 = disease area acres,
#   7 = disease type (else MINRR).
function rd_in_rrinit!(rd::RootDiseaseState, r)
    idi = Int(rd.minrr)
    (length(r.present) >= 7 && r.present[7]) && (idi = clamp(Int(nint(r.values[7])), 1, RD_ITOTRR))
    rd.rrman = true

    rndcenters = !(r.present[1] && r.values[1] >= 1.0f0)   # ARRAY(1) < 1 ⇒ random
    if rndcenters
        (r.present[2]) && (rd.ncents[idi] = Int32(min(100, Int(nint(r.values[2])))))
        if r.present[6]
            rd.parea[idi] = Float32(r.values[6]); rd.lparea[idi] = true
        end
        rd.ipcflg[idi] = Int32(0)
        (rd.ncents[idi] == 1 && rd.parea[idi] == rd.sarea) && (rd.lonect[idi] = Int32(1))
    else
        # explicit-center branch (IPCFLG=1) reads center geometry records — Chunk 1.
        (r.present[2]) && (rd.ncents[idi] = Int32(min(100, Int(nint(r.values[2])))))
        rd.parea[idi] = -1.0f0; rd.lparea[idi] = true; rd.ipcflg[idi] = Int32(1)
    end
    # field 3/4 = infected/uninfected TPA in diseased area; field 5 = initial root-infection proportion
    (r.present[3] && r.values[3] >= 0.0f0) && (rd.prkill[idi] = Float32(r.values[3]))
    (r.present[4] && r.values[4] >= 0.0f0) && (rd.prun[idi]   = Float32(r.values[4]))
    (r.present[5] && 0.0f0 <= r.values[5] <= 1.0f0) && (rd.rrincs[idi] = Float32(r.values[5]))
    return nothing
end

# rd/rdin.f option 10 — SAREA: stand area (acres). DIMEN = sqrt(SAREA)*208.7.
function rd_in_sarea!(rd::RootDiseaseState, r)
    if r.present[1] && r.values[1] > 0.0f0
        rd.sarea = Float32(r.values[1])
        rd.dimen = sqrt(rd.sarea) * 208.7f0
    end
    return nothing
end

# rd/rdin.f option 9 — END: defaults not-set patch areas and finalizes LONECT.
# (The default-bark-beetle auto-enable is deferred to Chunk 7; noted via lbbon.)
function rd_in_end!(rd::RootDiseaseState)
    @inbounds for idi in Int(rd.minrr):Int(rd.maxrr)
        rd.lparea[idi] || (rd.parea[idi] = 0.25f0 * rd.sarea)
        (rd.lonect[idi] == 0 && rd.ncents[idi] == 1 && rd.parea[idi] == rd.sarea) &&
            (rd.lonect[idi] = Int32(1))
    end
    return nothing
end

# -----------------------------------------------------------------------------
# Engine seams (rd/rdmn1.f, rd/rdmn2.f, rd/rdtreg.f) — wired GATED + INERT.
# The per-cycle mortality/spread/growth-loss bodies are Chunk 0 (pending); until
# then these are no-ops, so a stand with no RD keyword is byte-identical and a
# stand WITH an RD keyword parses/initializes but applies no mortality yet.
# -----------------------------------------------------------------------------

"""
    root_disease_setup!(s)

fvs.f RDMN1 init seam (called once at stand setup). Inert unless RD is active.
Chunk 0 will call the ported RDSETP here (center placement + initial infection).
"""
function root_disease_setup!(s::StandState)
    rd = s.root_disease
    rd_active(rd) || return nothing
    # Chunk 0 (pending): RDSETP → RDCLOC/RDAREA, RDIPRP, per-record PROBI/PROPI, RDINOC(true).
    return nothing
end

"""
    root_disease_mn2!(s, fint)

grincr.f RDMN2 seam (each cycle, before tripling). Advances the time-slot counter.
Inert (no mortality) until Chunk 0.
"""
function root_disease_mn2!(s::StandState, fint::Real)
    rd = s.root_disease
    (rd_active(rd) && rd.iroot != 0) || return nothing
    # Chunk 0 (pending): ISTEP+=1; IYEAR+=fint; zero RRKILL; RDSUM; RDSHST.
    return nothing
end

"""
    root_disease_treg!(s, fint)

gradd.f RDTREG seam (each cycle, after growth). The main per-cycle RD driver:
Chunk 0 will run RDCNTL (RDSETP→spread→RDMORT) + RDEND (→WK2 mortality) + RDGROW
(growth loss) + RDINOC. Inert until then.
"""
function root_disease_treg!(s::StandState, fint::Real)
    rd = s.root_disease
    (rd_active(rd) && rd.iroot != 0) || return nothing
    # Chunk 0 (pending): the reduced RDCNTL→RDMORT mortality path.
    return nothing
end
