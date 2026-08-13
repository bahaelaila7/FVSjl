# =============================================================================
# newspred.jl (britishcolumbia) — BC NEWSPRED / NISI spatial dwarf-mistletoe (#196).
# Port of canada/newmist (the New & Improved Spread & Intensification model).
# See docs/BC_NEWSPRED_PORT_PLAN.md for the C0-C7 chunk plan + validation vehicle.
#
# CHUNK C1 (this file, so far): DMCOM state → MistletoeState + the DMINIT default
# tables (DMDMR crown-third distribution, DMOPAQ species opacity) + the keyword
# params (NEWSPRED/DMAUTO/MISTPRT). NOT yet wired into the engine — the spread
# (C4), life-history (C5), and stand-coupling/mortality (C6) land in later chunks.
# Until C6, a NEWSPRED-active stand runs DM-FREE (the current validated baseline);
# C1 must stay .sum-INERT.
# =============================================================================

# --- DMCOM PARAMETERs (canada/newmist/DMCOM.F77) ---
const DM_CRTHRD = 3          # crown thirds
# life-history compartments (DMINF 3rd index): IMMATURE / LATENT / SUPPRESSED / ACTIVE / DEAD
const DM_IMMAT  = 1
const DM_LATENT = 2
const DM_SUPRSD = 3
const DM_ACTIVE = 4
const DM_DEAD   = 5
const DM_NPOOL  = 5
const DM_MAXDMR = 6          # DMR (dwarf-mistletoe rating) classes 0..6

# DMDMR[dmr(0:6), crownthird(1:3)] — default initial DM distribution across crown
# thirds for a tree of rating `dmr` (dminitbc.f TPDMR DATA, column-major). Row index
# here is dmr+1 (Julia 1-based; dmr 0 → row 1). E.g. dmr6→(2,2,2), dmr3→(0,1,2),
# dmr1→(0,0,1). Overridable by the DMCRTHRD keyword.
const DM_DMDMR = Int32[
#   ct1 ct2 ct3     dmr
    0   0   0   ;#  0
    0   0   1   ;#  1
    0   0   2   ;#  2
    0   1   2   ;#  3
    0   2   2   ;#  4
    1   2   2   ;#  5
    2   2   2    #  6
]

# DMOPAQ[sp] — relative crown opacity per BC species (dminitbc.f TPOPAQ DATA; April
# 1993 Model Review Workshop Report Table 4.1). BC MAXSP=15 species order:
# PW LW FD BG HW CW PL SE BL PY EP AT AC OC OH. Overridable by the DMOPQ keyword.
const DM_OPAQ = Float32[1.2, 0.9, 1.5, 1.8, 2.0, 1.8, 1.0, 1.7, 1.8, 1.0, 2.0, 2.0, 2.0, 1.5, 2.0]

"""
NISI spatial dwarf-mistletoe state (canada/newmist DMCOM), one per BC stand. `nothing`
on the stand until a NEWSPRED/MISTOE keyword activates it. Per-tree×crown-third×
compartment infection pools (`dminf`) + per-tree DMR (`dmr`); the spatial-grid /
shade / trajectory arrays (DMRDMX/CShd/Shd1/…) are added with the spread core (C4).
"""
mutable struct MistletoeState <: AbstractMistletoeState
    active::Bool                       # MISTOE keyword seen (DM extension on)
    newmod::Bool                       # NEWSPRED — use the NISI spatial spread model (misin.f opt 12)
    prtmis::Bool                       # MISTPRT — emit the DM reports (misin.f opt 6)
    dmrmin::Float32                    # MISTPRT field-1 min DMR to report (default 1.0)
    dmalpha::Float32                   # DMAUTO like-class autocorrelation decay (misin.f opt 24; −999 = unset)
    dmbeta::Float32                    # DMAUTO unlike-class decay (−999 = unset)
    # per-DMR-class crown-third distribution (DMDMR, keyword-overridable copy of DM_DMDMR)
    dmdmr::Matrix{Int32}               # (7, 3)
    opaq::Vector{Float32}              # per-species opacity (copy of DM_OPAQ, DMOPQ-overridable)
    # per-tree state (sized to the stand's live-record count when DMINIT runs)
    dmr::Vector{Int32}                 # per-tree dwarf-mistletoe rating 0..6
    dminf::Array{Float32,3}            # (tree, crownthird 1:3, compartment 1:5) infection pools
    rnseed::Int64                      # DMRNSD spatial-model RNG seed (own stream; never FFI'd)
end

MistletoeState() = MistletoeState(false, false, false, 1.0f0, -999f0, -999f0,
                                  copy(DM_DMDMR), copy(DM_OPAQ),
                                  Int32[], Array{Float32,3}(undef, 0, DM_CRTHRD, DM_NPOOL),
                                  0)

# --- C1 keyword handlers (misin.f) — recognize the DM keywords the YSM stand uses.
# BC-only: for other variants these keywords stay in `unrecognized_keywords` (unchanged
# behaviour). Until C6 wires the model, these only set flags/state and are .sum-INERT.
_dm_state!(s) = (s.mistletoe === nothing && (s.mistletoe = MistletoeState()); s.mistletoe::MistletoeState)

function kw_mistoe!(s::StandState, rec)
    if s.variant isa BritishColumbia
        _dm_state!(s).active = true                         # MISTOE: turn the DM extension on
    else
        push!(s.control.unrecognized_keywords, "MISTOE")
    end
    return
end

function kw_newspred!(s::StandState, rec)
    if s.variant isa BritishColumbia
        _dm_state!(s).newmod = true                         # misin.f opt 12: NEWMOD=.TRUE. (NISI spatial model)
    else
        push!(s.control.unrecognized_keywords, "NEWSPRED")
    end
    return
end

function kw_dmauto!(s::StandState, rec)
    if s.variant isa BritishColumbia
        ms = _dm_state!(s)
        # misin.f opt 24: field1=date, field2=DMALPHA, field3=DMBETA. A supplied >0 value is an
        # error sentinel (→ −999); a blank stays −999 (detected/defaulted later in DMOPTS).
        if length(rec.present) >= 2 && rec.present[2]
            a = Float32(rec.values[2]); ms.dmalpha = a > 0f0 ? -999f0 : a
        end
        if length(rec.present) >= 3 && rec.present[3]
            b = Float32(rec.values[3]); ms.dmbeta = b > 0f0 ? -999f0 : b
        end
    else
        push!(s.control.unrecognized_keywords, "DMAUTO")
    end
    return
end

function kw_mistprt!(s::StandState, rec)
    if s.variant isa BritishColumbia
        ms = _dm_state!(s); ms.prtmis = true                # misin.f opt 6: request DM reports
        (length(rec.present) >= 1 && rec.present[1]) && (ms.dmrmin = Float32(rec.values[1]))  # min DMR to report (default 1.0)
    else
        push!(s.control.unrecognized_keywords, "MISTPRT")
    end
    return
end
