# =============================================================================
# organon_interface.jl — OC (Oregon Coast) FVS↔ORGANON boundary marshalling (chunk C1).
#
# Ported from:
#   • oc/orgspc.f            — ORGSPC / DATA OSPMAP (50 FVS-index → ORGANON FIA code map)
#   • oc/dgdriv.f:217-285    — per-tree eligibility IORG(I), the big-6 stand gate, and the
#                              /ORGANON/ input-buffer fill (SPECIES,DBH1,HT1OR,CR1,SCR1B,
#                              EXPAN1,MGEXP,USER) that FVS hands to ORGANON's EXECUTE.
#
# SCOPE: C1 is the boundary MARSHALLING ONLY. It does NOT run ORGANON growth (that is the
# unported organon/ engine — chunks C3-C6) and it does NOT do the ORGANON PREPARE setup
# calibration / HT-CR dubbing (chunk C2). `build_organon_buffer!` produces the exact tree list
# that would be passed to ORGANON EXECUTE for the current cycle, given the CURRENT StandState.
#
# MEASURED against the live FVSoc_clean oracle (scoped `DEBUG / CRATET DGDRIV`, stand S248112 /
# ocmin, 27 records): the dgdriv `FOR EXECUTE` per-tree dump (I,TREENO,PTNO,DBH1,HT1OR,CR1,
# SCR1B,EXPAN1,MGEXP,USER) and the CRATET `SPECIES` dump, with NVALID=17, NBIG6=17, gate=RUN.
# See docs/OC_VARIANT_PORT_AUDIT.md (C1 verdict).
#
# The two FVS fill seams (oc/cratet.f:226 setup, oc/dgdriv.f:268 per-cycle) share this exact
# ORGSPC map + IORG/big-6 gate + DBH1/HT1OR/CR1 marshalling; they differ only in EXPAN1 (dgdriv
# = PROB, cratet = PROB·PI) and USER (dgdriv = ISPECL, cratet = KUTKOD). This ports the DGDRIV
# (per-cycle growth-entry) form, which is the C1 target.
# =============================================================================

"""
    OC_OSPMAP

`oc/orgspc.f` `DATA OSPMAP` — the 50-element map from FVS species index (1..50) to the ORGANON
SWO species FIA code. EVERY FVS species maps to some FIA code so ORGANON gets stand density right,
but only the 18 in `OC_ORGANON_VALID` are grown by ORGANON; the other 32 pass as a surrogate FIA
code and grow FVS-native.
"""
const OC_OSPMAP = Int32[
#    PC   IC   RC   GF   RF   SH   DF   WH   MH   WB
    242,  81, 242,  17, 202, 202, 202, 263, 263, 122,
#    KP   LP   CP   LM   JP   SP   WP   PP   MP   GP
    122, 122, 122, 122, 122, 117, 122, 122, 122, 122,
#    WJ   BR   GS   PY   OS   LO   CY   BL   EO   WO
    231,  17, 122, 231, 231, 805, 805, 805, 805, 815,
#    BO   VO   IO   BM   BU   RA   MA   GC   DG   FL
    818, 805, 805, 312, 312, 351, 361, 431, 492, 312,
#    WN   TO   SY   AS   CW   WI   CN   CL   OH   RW
    312, 631, 312, 312, 312, 920, 492, 492, 492, 122,
]

"""
    OC_ORGANON_VALID

The 18 FVS species indices that map to a VALID ORGANON SWO species (grown by ORGANON):
2=IC 3=RC 4=GF/WF 7=DF 8=WH 16=SP 18=PP 24=PY 27=CY 30=WO 31=BO 34=BM 36=RA 37=MA 38=GC 39=DG
42=TO 46=WI. (oc/dgdriv.f:234, oc/cratet.f:193.)
"""
const OC_ORGANON_VALID = (2, 3, 4, 7, 8, 16, 18, 24, 27, 30, 31, 34, 36, 37, 38, 39, 42, 46)

"""
    OC_ORGANON_BIG6

The 5 FVS species indices that are ORGANON "big-6" trees (IC, GF/WF, DF, SP, PP): 2, 4, 7, 16, 18.
The stand runs ORGANON only if it has ≥1 big-6 tree (oc/dgdriv.f:227, oc/cratet.f:186). (The "big
6" name counts WF/GF as one; FVS-OC maps both to GF=4, so 5 FVS indices cover the six.)
"""
const OC_ORGANON_BIG6 = (2, 4, 7, 16, 18)

"ORGSPC (oc/orgspc.f): FVS species index → ORGANON FIA code."
@inline orgspc(inspec::Integer) = @inbounds OC_OSPMAP[inspec]

"""
    OrganonBuffer

The FVS→ORGANON input tree list for one cycle — the Julia analogue of the `/ORGANON/` common-block
input arrays (`common/ORGANON.F77:28-35`) filled at `oc/dgdriv.f:268-285`. Only records `1:ntrees`
are active. `runs` is the big-6 stand-gate decision: `false` ⇒ no big-6 tree, all `iorg=0`,
`mortexp=0`, ORGANON EXECUTE is skipped and the stand reverts to FVS-native growth.
"""
mutable struct OrganonBuffer
    ntrees ::Int                 # IREC1 / ITRN — active record count
    nbig6  ::Int                 # count of big-6 eligible trees (dgdriv NBIG6)
    nvalid ::Int                 # count of valid ORGANON trees (Σ iorg)
    runs   ::Bool                # big-6 gate: true ⇒ ORGANON EXECUTE runs
    treeno ::Vector{Int32}       # TREENO(I) = I
    ptno   ::Vector{Int32}       # PTNO(I)   = ITRE(I) (point number)
    species::Vector{Int32}       # SPECIES(I) = ORGSPC(ISP(I)) (ORGANON FIA code)
    iorg   ::Vector{Int32}       # IORG(I) 0/1 valid-ORGANON-tree flag
    dbh1   ::Vector{Float32}     # DBH1(I)   = max(DBH,0.1)
    ht1or  ::Vector{Float32}     # HT1OR(I)  = (HT<4.6 ? 4.6 : HT)
    cr1    ::Vector{Float32}     # CR1(I)    = ICR/100
    scr1b  ::Vector{Float32}     # SCR1B(I)  = 0 (shadow crown ratio)
    expan1 ::Vector{Float32}     # EXPAN1(I) = PROB(I)
    mgexp  ::Vector{Float32}     # MGEXP(I)  = 0 (management expansion factor)
    user   ::Vector{Int32}       # USER(I)   = ISPECL(I)
    mortexp::Vector{Float32}     # MORTEXP(I) (zeroed here; filled by ORGANON mortality, C6)
end

function OrganonBuffer(n::Int)
    iz() = zeros(Int32, n); fz() = zeros(Float32, n)
    OrganonBuffer(0, 0, 0, false, iz(), iz(), iz(), iz(), fz(), fz(), fz(), fz(), fz(), fz(), iz(), fz())
end

"""
    build_organon_buffer!(s) -> OrganonBuffer

Port of `oc/dgdriv.f:217-285`: set the per-tree ORGANON eligibility flag `IORG`, count the big-6
trees, decide the stand gate, and (when the gate opens) fill the `/ORGANON/` input buffer from the
current `StandState.trees`.

Eligibility (`iorg=1`): `HT > 4.5 AND DBH >= 0.1 AND species ∈ OC_ORGANON_VALID`.
Big-6 gate: `nbig6 = Σ(HT>4.5 AND DBH>=0.1 AND species ∈ OC_ORGANON_BIG6)`; if `nbig6 == 0` all
`iorg` are forced 0, `mortexp=0`, and `runs=false` (revert to FVS-native — the buffer past `species`
is left unfilled, exactly as FVS's `GO TO 261` skips the fill loop).

This does NOT run ORGANON growth or the PREPARE calibration/dubbing — it only marshals the tree
list. (C1.)
"""
function build_organon_buffer!(s::StandState)
    t = s.trees
    n = t.n
    buf = OrganonBuffer(n)
    buf.ntrees = n
    nbig6 = 0
    # --- IORG eligibility + big-6 count + surrogate/actual FIA species (dgdriv.f:217-248) ---
    @inbounds for i in 1:n
        sp = Int(t.species[i])
        if t.height[i] > 4.5f0 && t.dbh[i] >= 0.1f0
            (sp in OC_ORGANON_BIG6) && (nbig6 += 1)
            buf.iorg[i] = (sp in OC_ORGANON_VALID) ? Int32(1) : Int32(0)
        else
            buf.iorg[i] = Int32(0)
        end
        buf.species[i] = orgspc(sp)                     # ORGSPC: actual for valid, surrogate otherwise
    end
    buf.nbig6 = nbig6
    # --- big-6 stand gate (dgdriv.f:254-260) ---
    if nbig6 == 0
        @inbounds for i in 1:n
            buf.iorg[i]    = Int32(0)
            buf.mortexp[i] = 0f0
        end
        buf.nvalid = 0
        buf.runs = false
        return buf
    end
    buf.runs = true
    # --- /ORGANON/ input-buffer fill (dgdriv.f:268-285) ---
    nvalid = 0
    @inbounds for i in 1:n
        buf.treeno[i] = Int32(i)
        buf.ptno[i]   = t.plot_id[i]
        d = t.dbh[i];  d < 0.1f0 && (d = 0.1f0)
        buf.dbh1[i]   = d
        # HT1OR: floor to 4.6 ONLY when HT>0 (oc/cratet.f:234); a MISSING height (HT==0) is passed
        # as 0.0 so ORGANON PREPARE flags it MISSHT and dubs it (chunk C2). Flooring 0→4.6 here was
        # the C1 tree-20 residual; this matches the oracle bit-exact.
        h = t.height[i]; (h > 0f0 && h < 4.6f0) && (h = 4.6f0)
        buf.ht1or[i]  = h
        buf.cr1[i]    = Float32(t.crown_pct[i]) / 100f0
        buf.scr1b[i]  = 0f0
        buf.expan1[i] = t.tpa[i]
        buf.mgexp[i]  = 0f0
        buf.user[i]   = t.special[i]
        nvalid += Int(buf.iorg[i])
    end
    buf.nvalid = nvalid
    return buf
end

# --- Per-cycle growth entry (oc/dgdriv.f) — C1 boundary only -----------------------------------
# On OregonCoast() the growth entry builds the ORGANON input tree list (the C1 marshalling), then
# STOPS: the ORGANON growth engine (organon/ diagro/htgrowth/crngrow/mortality — chunks C3-C6) is
# unported, so we do NOT run growth. Doctrine #5: an unported growth hook errors loudly. The
# `build_organon_buffer!` call is retained so the C1 marshalling still exercises here (and so a
# future C2/C3 can slot the ORGANON EXECUTE call in right after it).
function diameter_growth!(s::StandState, ::OregonCoast; kwargs...)
    build_organon_buffer!(s)   # C1: fill the /ORGANON/ input buffer for this cycle
    error("OregonCoast (OC) growth is the ORGANON SWO engine. Ported so far: C1 (boundary " *
          "marshalling), C2 (PREPARE calibration), C3 (DG_SWO diameter growth — `organon_dg_swo`), " *
          "C4 (HG_SWO height growth — `organon_hg_swo`); all bit-exact. Still UNPORTED: crown (C5), " *
          "mortality (C6), and the GROW/EXECUTE per-cycle orchestration + FVS DDS/HTG copy-back " *
          "(C7). See docs/OC_ORGANON_PORT_PLAN.md.")
end
