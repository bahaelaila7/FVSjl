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
const DM_BPCNT  = DM_CRTHRD + 1   # crown-third breakpoints (4): top, 2 dividers, bottom

# --- DMCOM.F77 spatial-grid PARAMETERs (MESH geometry) ---
# The spread model works in MESH units — a MESH-metre spatial grid (default 2 m).
const DM_MESH   = 2               # MESH grid-cell size (metres)
const DM_FPM    = 3.2808f0        # feet per metre (HT() is feet, model wants MESH)
const DM_MXHT   = 50 ÷ DM_MESH    # max stand height in MESH (25)
const DM_MXTHRX = 14 ÷ DM_MESH    # max lateral seed x-travel in MESH (7) = # sampling rings
const DM_ORIGIN = 20 ÷ DM_MESH    # trajectory origin cell (10)
const DM_TWOPIE = 6.283185f0      # 2π (subtended-angle interception, dmtreg.f)

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

# --- Crown-shape discriminant (dmshap.f) — assigns each tree crown to one of 5 shapes
# (1=sphere/circle, 2=cone/triangle, 3=neiloid, 4=paraboloid, 5=ellipsoid) via Fisher's
# linear discriminant (Moeur crown-shape models). BC's 15 species map to 11 discriminant
# groups (dmshap.f MAPBC): PW LW FD BG HW CW PL SE BL PY EP AT AC OC OH →
# [1,2,3,4,5,6,7,8,9,10,6,6,6,3,6] (EP/AT/AC/OH→RC grp6, OC→DF grp3). Groups: 1=WP 2=WL
# 3=DF 4=GF 5=WH 6=RC 7=LP 8=ES 9=AF 10=PP 11=MH. −99 coeff = shape disallowed for that
# group (its score goes hugely negative → never selected). SCORE(shape) = CONST + BCR·CR +
# BHT·HT + BCL·CL + BRAD·RAD + BDBH·DBH + BTPA·TPA; highest score wins (ties keep lower shape).
const DM_MAPBC = Int32[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 6, 6, 6, 3, 6]

const DM_SHP_CONST = reshape(Float32[
    -13.943f0, -37.395f0, -99f0, -32.104f0, -32.042f0,
    -99f0, -38.278f0, -99f0, -44.875f0, -32.742f0,
    -25.38f0, -31.308f0, -35.906f0, -31.435f0, -28.185f0,
    -19.633f0, -23.785f0, -20.062f0, -22.012f0, -22.962f0,
    -19.633f0, -23.785f0, -20.062f0, -22.012f0, -22.962f0,
    -99f0, -27.573f0, -99f0, -22.434f0, -26.275f0,
    -22.322f0, -37.19f0, -99f0, -38.316f0, -31.988f0,
    -99f0, -46.116f0, -99f0, -49.99f0, -43.456f0,
    -24.002f0, -37.611f0, -99f0, -30.234f0, -36.183f0,
    -20.296f0, -24.195f0, -36.526f0, -26.622f0, -22.397f0,
    -24.002f0, -37.611f0, -99f0, -30.234f0, -36.183f0,
], 5, 11)  # [shape 1:5, group 1:11]
const DM_SHP_BCR = reshape(Float32[
    49.794f0, 100.24f0, -99f0, 86.383f0, 94.554f0,
    -99f0, 123.805f0, -99f0, 134.78f0, 116.554f0,
    52.715f0, 60.761f0, 63.156f0, 61.744f0, 59.922f0,
    45.992f0, 54.464f0, 38.042f0, 51.755f0, 53.598f0,
    45.992f0, 54.464f0, 38.042f0, 51.755f0, 53.598f0,
    -99f0, 59.56f0, -99f0, 53.273f0, 55.799f0,
    73.92f0, 104.413f0, -99f0, 102.38f0, 99.877f0,
    -99f0, 106.784f0, -99f0, 110.115f0, 100.844f0,
    56.463f0, 80.909f0, -99f0, 67.587f0, 78.258f0,
    45.808f0, 57.54f0, 60.76f0, 59.567f0, 57.457f0,
    56.463f0, 80.909f0, -99f0, 67.587f0, 78.258f0,
], 5, 11)  # [shape 1:5, group 1:11]
const DM_SHP_BHT = reshape(Float32[
    0.48348f0, 0.92054f0, -99f0, 0.81786f0, 0.85462f0,
    -99f0, 1.03836f0, -99f0, 1.00712f0, 1.02891f0,
    1.02304f0, 1.13332f0, 1.28061f0, 1.13364f0, 1.16673f0,
    0.72015f0, 0.88508f0, 0.56218f0, 0.84169f0, 0.88962f0,
    0.72015f0, 0.88508f0, 0.56218f0, 0.84169f0, 0.88962f0,
    -99f0, 1.08626f0, -99f0, 0.98588f0, 0.98314f0,
    0.83497f0, 1.16782f0, -99f0, 1.17028f0, 1.15645f0,
    -99f0, 1.60533f0, -99f0, 1.62982f0, 1.59203f0,
    1.25278f0, 1.59728f0, -99f0, 1.38587f0, 1.5867f0,
    0.71361f0, 0.91473f0, 0.98645f0, 0.89231f0, 0.91047f0,
    1.25278f0, 1.59728f0, -99f0, 1.38587f0, 1.5867f0,
], 5, 11)  # [shape 1:5, group 1:11]
const DM_SHP_BRAD = reshape(Float32[
    0.56507f0, 1.1108f0, -99f0, 1.82564f0, 1.32317f0,
    -99f0, 1.08869f0, -99f0, 1.8112f0, 0.79969f0,
    -0.10457f0, 0.28897f0, -0.49614f0, 0.07177f0, -0.04803f0,
    0.13463f0, 0.3542f0, 0.24447f0, 0.32306f0, 0.18187f0,
    0.13463f0, 0.3542f0, 0.24447f0, 0.32306f0, 0.18187f0,
    -99f0, 0.92215f0, -99f0, 0.91886f0, 0.48536f0,
    0.62301f0, 1.6984f0, -99f0, 1.57036f0, 1.08601f0,
    -99f0, 1.12257f0, -99f0, 0.86087f0, 0.8917f0,
    1.64765f0, 1.72301f0, -99f0, 1.83248f0, 1.53816f0,
    0.98573f0, 0.60433f0, 2.22885f0, 0.89341f0, 0.5818f0,
    1.64765f0, 1.72301f0, -99f0, 1.83248f0, 1.53816f0,
], 5, 11)  # [shape 1:5, group 1:11]
const DM_SHP_BCL = reshape(Float32[
    -0.8373f0, -1.5604f0, -99f0, -1.3933f0, -1.4959f0,
    -99f0, -1.9316f0, -99f0, -2.1712f0, -1.8977f0,
    -1.3082f0, -1.4288f0, -1.5055f0, -1.4393f0, -1.4367f0,
    -1.1754f0, -1.3304f0, -0.9574f0, -1.3074f0, -1.325f0,
    -1.1754f0, -1.3304f0, -0.9574f0, -1.3074f0, -1.325f0,
    -99f0, -1.2529f0, -99f0, -1.1295f0, -1.1502f0,
    -1.6455f0, -2.2164f0, -99f0, -2.1892f0, -2.1148f0,
    -99f0, -2.2233f0, -99f0, -2.3029f0, -2.1475f0,
    -1.6686f0, -2.0343f0, -99f0, -1.8193f0, -2.0154f0,
    -1.0722f0, -1.3127f0, -1.3903f0, -1.3079f0, -1.2973f0,
    -1.6686f0, -2.0343f0, -99f0, -1.8193f0, -2.0154f0,
], 5, 11)  # [shape 1:5, group 1:11]
const DM_SHP_BTPA = reshape(Float32[
    0.00437f0, 0.00911f0, -99f0, 0.00603f0, 0.00899f0,
    -99f0, 0.01548f0, -99f0, 0.02342f0, 0.02134f0,
    0.00634f0, 0.00706f0, 0.00673f0, 0.00683f0, 0.00673f0,
    0.00782f0, 0.00841f0, 0.00813f0, 0.00778f0, 0.00816f0,
    0.00782f0, 0.00841f0, 0.00813f0, 0.00778f0, 0.00816f0,
    -99f0, 0.00896f0, -99f0, 0.00736f0, 0.0081f0,
    0.006f0, 0.00406f0, -99f0, 0.00542f0, 0.00613f0,
    -99f0, 0.04348f0, -99f0, 0.04908f0, 0.0541f0,
    0.01644f0, 0.01891f0, -99f0, 0.01809f0, 0.01923f0,
    0.01233f0, 0.01068f0, 0.01136f0, 0.0119f0, 0.0107f0,
    0.01644f0, 0.01891f0, -99f0, 0.01809f0, 0.01923f0,
], 5, 11)  # [shape 1:5, group 1:11]
const DM_SHP_BDBH = reshape(Float32[
    0.05236f0, 0.08377f0, -99f0, -0.17327f0, 0.22765f0,
    -99f0, -0.1477f0, -99f0, 0.25077f0, 0.0121f0,
    -0.40495f0, -0.66563f0, -0.42355f0, -0.50661f0, -0.5478f0,
    0.14819f0, -0.02081f0, 0.14503f0, 0.1476f0, 0.00202f0,
    0.14819f0, -0.02081f0, 0.14503f0, 0.1476f0, 0.00202f0,
    -99f0, -0.81991f0, -99f0, -0.99031f0, -0.33308f0,
    0.21973f0, -0.17842f0, -99f0, 0.00182f0, -0.04777f0,
    -99f0, 0.06301f0, -99f0, 0.41333f0, 0.0243f0,
    -0.61101f0, -0.80165f0, -99f0, -0.74084f0, -0.72773f0,
    -0.39061f0, -0.33475f0, -0.99571f0, -0.4074f0, -0.37839f0,
    -0.61101f0, -0.80165f0, -99f0, -0.74084f0, -0.72773f0,
], 5, 11)  # [shape 1:5, group 1:11]

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
    brkpnt::Matrix{Float32}            # (tree, BPCNT 1:4) crown-third breakpoints in MESH units (DMFBRK)
    idmshp::Vector{Int32}              # per-tree crown shape 1:5 (DMSHAP Fisher discriminant)
    dmrdmx::Array{Float32,3}           # (tree, MESH band 1:MXHT, {RADIUS=1,VOLUME=2}) crown frustum geometry (DMSUM)
    rnseed::Int64                      # DMRNSD spatial-model RNG seed (own stream; never FFI'd)
end

MistletoeState() = MistletoeState(false, false, false, 1.0f0, -999f0, -999f0,
                                  copy(DM_DMDMR), copy(DM_OPAQ),
                                  Int32[], Array{Float32,3}(undef, 0, DM_CRTHRD, DM_NPOOL),
                                  Matrix{Float32}(undef, 0, DM_BPCNT), Int32[],
                                  Array{Float32,3}(undef, 0, DM_MXHT, 2),
                                  0)

# DMRDMX 3rd-index tags (DMCOM RADIUS/VOLUME)
const DM_RADIUS = 1
const DM_VOLUME = 2

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

# --- C1 DMINIT (partial): seed the per-tree initial DMR from the input damage codes.
# misdam.f: damage codes 30-34 are dwarf mistletoe (30 generic / 31 LP / 32 WL / 33 DF /
# 34 PP); IMIST/DMRATE = the SEVERITY that follows the code, capped 0..6. jl already stores
# the codes in `t.damage` (6×MAXTRE: dmg1/sev1/dmg2/sev2/dmg3/sev3). Runs once at BC setup
# when the DM model is active. Sizes dmr/dminf to the live-record count. dminf COMPARTMENT
# seeding (which life-history pool the initial rating enters) is deferred to C5.
# Still .sum-INERT — nothing consumes dmr until C6.
function dm_init!(s::StandState)
    ms = s.mistletoe
    (ms === nothing || !(ms.active || ms.newmod)) && return s
    t = s.trees; n = t.n
    ms.dmr = zeros(Int32, n)
    ms.dminf = zeros(Float32, n, DM_CRTHRD, DM_NPOOL)
    ms.brkpnt = zeros(Float32, n, DM_BPCNT)
    ms.idmshp = zeros(Int32, n)
    ms.dmrdmx = zeros(Float32, n, DM_MXHT, 2)
    @inbounds for i in 1:n, j in (1, 3, 5)
        ag = Int(t.damage[j, i])
        (30 <= ag <= 34) && (ms.dmr[i] = Int32(clamp(Int(t.damage[j+1, i]), 0, 6)))  # last DM code wins (misdam.f)
    end
    return s
end

# --- C4 geometry substrate: DMFBRK (dmfbrk.f) — locate the 4 breakpoints defining each
# tree's crown thirds, in MESH units. BrkPnt[i,1] = TOP of crown, BrkPnt[i,4] = BOTTOM;
# [i,2]/[i,3] divide the thirds. HT is feet, ICR (crown_pct) is percent → MESH. This is the
# vertical spatial grid every spread/intensification loop in dmtreg reads. Deterministic;
# recomputed each cycle before spread. Still engine-INERT (nothing consumes brkpnt until C4/C6).
function dm_fbrk!(s::StandState)
    ms = s.mistletoe
    (ms === nothing || !(ms.active || ms.newmod)) && return s
    t = s.trees; n = t.n
    size(ms.brkpnt, 1) == n || (ms.brkpnt = zeros(Float32, n, DM_BPCNT))
    y = 1f0 / (DM_FPM * DM_MESH)                 # feet → MESH
    x = 0.01f0 * y / Float32(DM_CRTHRD)          # ICR% → MESH-per-crown-third scale
    @inbounds for i in 1:n
        z = t.height[i] * Float32(t.crown_pct[i]) * x   # MESH units per crown third
        ms.brkpnt[i, 1] = t.height[i] * y               # top of crown
        for j in 2:DM_BPCNT
            ms.brkpnt[i, j] = ms.brkpnt[i, j-1] - z
        end
    end
    return s
end

# --- C3 geometry: DMSHAP (dmshap.f) — assign each tree crown one of 5 geometric shapes via
# Fisher's linear discriminant. Feeds DMRDMX (the per-MESH-band frustum volume/radius). The
# stand TPA term is MAX(TPROB, OLDTPA) (PLOT.F77; pre-thin TPA). Deterministic; recomputed each
# cycle. Engine-INERT (nothing consumes idmshp until DMRDMX/spread land).
function dm_shap!(s::StandState)
    ms = s.mistletoe
    (ms === nothing || !(ms.active || ms.newmod)) && return s
    t = s.trees; n = t.n
    length(ms.idmshp) == n || (ms.idmshp = zeros(Int32, n))
    tpa = max(s.plot.total_tpa, s.plot.old_tpa)     # PLOT: TPA = MAX(TPROB, OLDTPA)
    @inbounds for i in 1:n
        sp = Int(t.species[i])
        g = (1 <= sp <= length(DM_MAPBC)) ? Int(DM_MAPBC[sp]) : 0
        g == 0 && continue
        rad = t.crown_width[i] * 0.5f0
        cr  = Float32(t.crown_pct[i]) / 100f0
        cl  = cr * t.height[i]
        best = Int32(1); bestsc = -999f0
        for j in 1:5
            sc = DM_SHP_CONST[j, g] + DM_SHP_BCR[j, g]*cr + DM_SHP_BHT[j, g]*t.height[i] +
                 DM_SHP_BCL[j, g]*cl + DM_SHP_BRAD[j, g]*rad + DM_SHP_BDBH[j, g]*t.dbh[i] +
                 DM_SHP_BTPA[j, g]*tpa
            sc > bestsc && (bestsc = sc; best = Int32(j))    # strict > ⇒ ties keep lower shape (dmshap.f .LE. skip)
        end
        ms.idmshp[i] = best
    end
    return s
end

# --- C3 geometry: DMRDMX (dmsum.f:96-256) — per-tree, per-MESH-band crown RADIUS & VOLUME,
# branching on the DMSHAP shape. Each band J (MESH height class) gets the frustum volume FRUST
# and the projected-area-equivalent radius RAD2, both in MESH units. This is what the spread
# loop reads (Level·HtWt·DMRDMX(VOLUME) for seed production; DMRDMX(RADIUS) for the subtended-
# angle interception). Deterministic; recomputed each cycle after dm_shap!. Engine-INERT.
# Faithful to dmsum.f; the dead J1 debug index and DEBUG writes are omitted. A crown with
# HC<=0 (0% crown) has no infectable volume ⇒ left at 0 (avoids the div-by-HC NaN; live never
# reaches such a tree with DM, crowns are dubbed >0).
const _DM_QRTRPI = 1.04720f0   # dmsum.f frustum constant (π/3)
const _DM_PIE    = 3.14159f0
const _DM_HLFPIE = 1.57080f0
function dm_rdmx!(s::StandState)
    ms = s.mistletoe
    (ms === nothing || !(ms.active || ms.newmod)) && return s
    t = s.trees; n = t.n
    (size(ms.dmrdmx, 1) == n && size(ms.idmshp, 1) == n) || (ms.dmrdmx = zeros(Float32, n, DM_MXHT, 2))
    fill!(ms.dmrdmx, 0f0)
    mscl = DM_FPM * DM_MESH                       # feet per MESH cell
    @inbounds for i in 1:n
        ht = t.height[i]
        hc = Float32(t.crown_pct[i]) * ht / 100f0    # crown length (feet)
        hc <= 0f0 && continue
        shp = Int(ms.idmshp[i])
        shp == 0 && continue
        itop = trunc(Int, ht / (mscl + 0.0001f0)) + 1
        itop > DM_MXHT && (itop = DM_MXHT)
        bot  = ht - hc
        ibot = trunc(Int, bot / mscl) + 1
        base = (shp == 1 || shp == 5) ? bot + hc/2f0 : bot
        rad  = t.crown_width[i] * 0.5f0
        cnop1 = mscl * ibot - 2f0*mscl
        cnop2 = cnop1 + mscl
        for j in ibot:itop
            cnop1 += mscl; cnop2 += mscl
            uplim  = min(ht, cnop2)
            j == itop && (uplim = ht)
            lowlim = max(bot, cnop1)
            h2 = uplim - lowlim
            h2 <= 0f0 && continue                   # degenerate band (no vertical extent) → 0
            h1 = lowlim - base
            # Only FRUST feeds the stored RADIUS/VOLUME (dmsum.f RAD2/FRUST); PAREA→RAD1 is
            # debug-only, so its Y/asin/^1.5 terms are omitted (same output, no NaN risk).
            local frust::Float32
            if shp == 1 || shp == 5                 # sphere / ellipsoid
                base >= lowlim && (h1 = base - uplim)
                cst = _DM_QRTRPI * h2 * rad*rad / (hc*hc)
                frust = cst * (3f0*hc*hc - 12f0*h1*h1 - 12f0*h1*h2 - 4f0*h2*h2)
            elseif shp == 2                         # cone / triangle
                r1 = (1f0 - h1/hc) * rad
                r2 = (1f0 - (h1+h2)/hc) * rad
                frust = (_DM_QRTRPI * h2) * (r1*r1 + r2*r2 + r1*r2)
            elseif shp == 3                         # neiloid
                frust = _DM_PIE*rad*rad*h2 - (_DM_HLFPIE*h2*rad*rad/hc)*(2f0*hc - 2f0*h1 - h2)
            else                                    # shp == 4, paraboloid
                frust = (_DM_HLFPIE*h2*rad*rad/hc)*(2f0*hc - 2f0*h1 - h2)
            end
            rad2 = sqrt(max(frust, 0f0) / (_DM_PIE * h2))
            ms.dmrdmx[i, j, DM_RADIUS] = rad2 / mscl
            ms.dmrdmx[i, j, DM_VOLUME] = frust / mscl^3
        end
    end
    return s
end

# --- C4 neighbour-count PDF: BNDIST + GAMMLN (bndist.f) — the Binomial/Poisson/Negative-Binomial
# family PDF for a population of given mean M and variance V (V≈M→Poisson, V>M→NegBinom, V<M→
# Binomial). Called by DMNB to distribute source trees across sampling rings. Pure deterministic
# (no RNG). Returns (pdf[1:ubound+1], endpos, err); pdf[j+1] = P(count=j). GAMMLN = Numerical-
# Recipes Lanczos log-gamma (float32, faithful to the Fortran cof/SqPI DATA).
const _DM_GAMMLN_COF = Float32[76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.12085003f-2, -0.536382f-5]
const _DM_SQPI = 2.50662827465f0
@inline function dm_gammln(arg::Float32)
    x = arg - 1f0
    tmp = x + 5.5f0
    tmp = tmp - (x + 0.5f0) * log(tmp)
    ser = 1f0
    @inbounds for j in 1:6
        x += 1f0
        ser += _DM_GAMMLN_COF[j] / x
    end
    return -tmp + log(_DM_SQPI * ser)
end

function dm_bndist(m::Float32, v::Float32, ubound::Int)
    tol = 1f-6
    pdf = zeros(Float32, ubound + 1)
    endpos = 1
    (m < tol || v < tol) && return (pdf, endpos, true)      # err: degenerate mean/variance
    local method::Int, t1::Float32, t2::Float32, t3::Float32, k::Float32
    if abs(v - m) < tol                                     # Poisson
        method = 1; t1 = -m; t2 = log(m); t3 = 0f0; k = 0f0
    elseif v > m                                            # Negative Binomial
        method = 2; p = (v / m) - 1f0; k = m / p
        t1 = -k * log(1f0 + m/k); t2 = dm_gammln(k); t3 = log(m / (m + k))
    else                                                    # Binomial
        method = 3; p = 1f0 - (v / m); k = m / p
        t1 = dm_gammln(k + 1f0); t2 = log(p); t3 = log(1f0 - p)
    end
    sum = 0f0; plast = 0f0
    for j in 0:ubound
        jp = j + 1; x = Float32(j)
        z = if method == 1
            t1 + x*t2 - dm_gammln(x + 1f0)
        elseif method == 2
            t1 + dm_gammln(k + x) - dm_gammln(x + 1f0) - t2 + x*t3
        else
            (k - x + 1f0) < 0f0 ? -99f0 :
                t1 - dm_gammln(x + 1f0) - dm_gammln(k - x + 1f0) + x*t2 + (k - x)*t3
        end
        z > -75f0 && (pdf[jp] = exp(z))
        sum += pdf[jp]
        if j == 0
            plast = pdf[jp]
        else
            pnow = pdf[jp]
            if pnow < plast && pnow < tol
                endpos = j; break
            elseif sum >= 1f0
                endpos = j; pdf[jp] -= (sum - 1f0); break
            elseif j == ubound
                endpos = j; break
            end
            plast = pnow
        end
    end
    return (pdf, endpos, false)
end
