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
const DM_MXTHRZ = 26 ÷ DM_MESH    # max vertical seed z-travel in MESH (13)
const DM_ORIGIN = 20 ÷ DM_MESH    # trajectory origin cell (10)
const DM_TWOPIE = 6.283185f0      # 2π (subtended-angle interception, dmtreg.f)
const _DM_PIE   = 3.14159f0       # dmcom PIE
const _DM_SQM2AC = 1f0 / 4046.8564f0   # m² → acres
const DM_DSTLEN = 1000            # max sampling-ring source-count array length (DMCOM DSTLEN)
const DM_TOP1   = 1496            # length of the Shd1 trajectory table (DMCOM TOP1)
const DM_MXTRAJ = 1              # max trajectories per grid cell (DMCOM MXTRAJ)
const DM_XX     = 1              # CShd/DMRDMX 3rd-index: x-position
const DM_ZZ     = 2              # CShd 3rd-index: z-position
# CrArea(i) = cumulative circle area (acres) of radius (MESH·i) m; Dstnce(i) = midpoint dist (m)
# of ring i. Compile-time (dminitbc.f:170-179), i = 1..MXTHRX.
const DM_CRAREA = Float32[_DM_PIE * _DM_SQM2AC * Float32(DM_MESH * i)^2 for i in 1:DM_MXTHRX]
const DM_DSTNCE = Float32[Float32(DM_MESH) * (Float32(i) - 0.5f0) for i in 1:DM_MXTHRX]

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
    dmclmp::Float32                    # clumping (variance/mean ratio) for the neighbour PDF (DMINIT 1.0; DMCLMP kw)
    # per-DMR-class crown-third distribution (DMDMR, keyword-overridable copy of DM_DMDMR)
    dmdmr::Matrix{Int32}               # (7, 3)
    opaq::Vector{Float32}              # per-species opacity (copy of DM_OPAQ, DMOPQ-overridable)
    # per-tree state (sized to the stand's live-record count when DMINIT runs)
    dmr::Vector{Int32}                 # per-tree dwarf-mistletoe rating 0..6
    dminf::Array{Float32,3}            # (tree, crownthird 1:3, compartment 1:5) infection pools
    brkpnt::Matrix{Float32}            # (tree, BPCNT 1:4) crown-third breakpoints in MESH units (DMFBRK)
    idmshp::Vector{Int32}              # per-tree crown shape 1:5 (DMSHAP Fisher discriminant)
    dmrdmx::Array{Float32,3}           # (tree, MESH band 1:MXHT, {RADIUS=1,VOLUME=2}) crown frustum geometry (DMSUM)
    sf::Matrix{Float32}                # (DMR-diff 0:6 → 1:7, ring 1:MXTHRX) autocorrelation scaling (DMINIT/DMAUTO)
    dms0::Float64                      # DMRANN current LCG state (own stream; never FFI'd)
    dmss::Float64                      # DMRANN saved seed (DMRNSD)
    rnseed::Int64                      # (reserved)
end

MistletoeState() = MistletoeState(false, false, false, 1.0f0, -999f0, -999f0, 1.0f0,
                                  copy(DM_DMDMR), copy(DM_OPAQ),
                                  Int32[], Array{Float32,3}(undef, 0, DM_CRTHRD, DM_NPOOL),
                                  Matrix{Float32}(undef, 0, DM_BPCNT), Int32[],
                                  Array{Float32,3}(undef, 0, DM_MXHT, 2),
                                  Matrix{Float32}(undef, 7, DM_MXTHRX),
                                  55329.0, 55329.0, 0)

# DMRDMX 3rd-index tags (DMCOM RADIUS/VOLUME)
const DM_RADIUS = 1
const DM_VOLUME = 2

# --- DMRANN (dmrann.f) — the DM spatial RNG: a MINSTD/Park-Miller LCG (16807, 2^31-1) whose
# output divides by 2^31, giving SEL ∈ (0,1). Float64 state (DMS0), seeded to 55329.0 (dminitbc.f:276).
# Its OWN stream — NEVER FFI'd and NEVER shared with the ZZRAN growth RNG; the spatial draws are an
# accepted realization straddle vs live (the whole-loop call order can't be byte-matched).
@inline function dm_rann!(ms::MistletoeState)
    dms1 = mod(16807.0 * ms.dms0, 2147483647.0)     # DMOD(16807D0*DMS0, 2147483647D0)
    sel = Float32(dms1 / 2147483648.0)              # /2^31
    ms.dms0 = dms1
    return sel
end
# DMRNSD: reset the stream to the saved seed (LSET=.FALSE. path, dminitbc.f:281). The LSET=.TRUE.
# set-path (odd-ify + store) is exposed for completeness though the only corpus seed is the default.
dm_rnsd_reset!(ms::MistletoeState) = (ms.dms0 = ms.dmss; ms)
function dm_rnsd_set!(ms::MistletoeState, seed::Float64)
    (mod(seed, 2.0) == 0.0) && (seed += 1.0)
    ms.dmss = seed; ms.dms0 = seed; ms
end

# --- DMSLOP (dmslop.f) — site-slope offset (in MESH height) between a source tree in sampling
# ring `ring` and its target: Offset = INT((ring-0.5)·cos(2π·Rnd)·SLOPE). + = source above target.
# Draws one DMRANN angle ⇒ stochastic. `slope` = stand slope fraction (PLOT SLOPE).
@inline function dm_slop(ms::MistletoeState, ring::Int, slope::Float32)
    d = Float32(ring) - 0.5f0
    rnd = dm_rann!(ms)
    return trunc(Int, d * cos(DM_TWOPIE * rnd) * slope)
end

# --- DMNDMR (dmndmr.f) — recompute each tree's DMR (dwarf-mistletoe rating 0..6) from the
# infection pools after a spread cycle. Per crown third: x = ACTIVE + SUPRSD + DEAD infection;
# rating k = 2 if x>2, else INT(x) with a stochastic +1 (DMRANN draw ≤ frac). DMR = Σ over the
# 3 thirds (so 0..6, Hawksworth). Biocontrol pools (DMINF_BC, MISBCI) are omitted — unported and
# zero for BC/YSM. Stochastic (fractional rounding). Runs after the spread core each cycle.
function dm_ndmr!(s::StandState)
    ms = s.mistletoe
    (ms === nothing || !(ms.active || ms.newmod)) && return s
    t = s.trees; n = t.n
    length(ms.dmr) == n || (ms.dmr = zeros(Int32, n))
    @inbounds for i in 1:n
        rate = 0
        for j in 1:DM_CRTHRD
            x = ms.dminf[i, j, DM_ACTIVE] + ms.dminf[i, j, DM_SUPRSD] + ms.dminf[i, j, DM_DEAD]
            if x > 2f0
                k = 2
            else
                k = trunc(Int, x)
                frac = x - Float32(k)
                dm_rann!(ms) <= frac && (k += 1)
            end
            rate += k
        end
        ms.dmr[i] = Int32(rate)
    end
    return s
end

# --- DMFINF (dmfinf.f) — build the treelist DM-index: sort records by (species, DMR) via
# OPSORT, then fill Ptr[sp, dmr+1, {FST=1,LST=2}] = the first/last positions of each
# (species, DMR) group in the sorted `index`. The spread loop iterates species×DMR through
# this Ptr. Deterministic (OPSORT unstable tie order). Core split out for testability.
function _dm_build_dm_index(species::AbstractVector{<:Integer}, dmr::AbstractVector{<:Integer}, n::Int)
    maxsp = n == 0 ? 1 : Int(maximum(@view species[1:n]))
    ptr = zeros(Int32, maxsp, 7, 2)      # [species, DMR 0:6→1:7, {1=FST, 2=LST}]
    index = zeros(Int32, n)
    n == 0 && return (ptr, index)
    opsort!(n, species, dmr, index, true)
    k = Int(index[1]); sp = Int(species[k]); dm = Int(dmr[k])
    prsp = sp; prdm = dm
    ptr[sp, dm+1, 1] = 1
    @inbounds for i in 2:n
        k = Int(index[i]); sp = Int(species[k]); dm = Int(dmr[k])
        if sp != prsp
            ptr[prsp, prdm+1, 2] = i - 1
            ptr[sp, dm+1, 1] = i
            prsp = sp; prdm = dm
        elseif dm != prdm
            ptr[sp, prdm+1, 2] = i - 1
            ptr[sp, dm+1, 1] = i
            prdm = dm
        end
    end
    ptr[sp, dm+1, 2] = n
    return (ptr, index)
end
dm_finf(s::StandState) = _dm_build_dm_index(s.trees.species, s.mistletoe.dmr, s.trees.n)

# --- DMFDNS (dmfdns.f) — trees/acre density of each target DMR class (0..6) for species `sp`:
# D[i+1] = Σ PROB over the trees in group (sp, i) via the Ptr range. Deterministic.
function dm_fdns(sp::Int, ptr, index, tpa)
    d = zeros(Float32, 7)
    @inbounds for i in 0:6
        fst = ptr[sp, i+1, 1]
        fst > 0 || continue
        acc = 0f0
        for j in fst:ptr[sp, i+1, 2]
            acc += tpa[index[j]]
        end
        d[i+1] = acc
    end
    return d
end

const DM_TINY = 1f-10   # DMCOM DMTINY

# --- DMSRC (dmsrc.f) — build the packed per-DMR-class cumulative source-selection vectors for
# species `sp`: for each class i, walk its trees accumulating cumulative proportional probability
# y += PROB·(1/(D[i]+tiny)) (rises 0→1 within the class); SrcI[k]=tree, SrcCD[k]=cum-prob,
# SPtr[i+1]=end offset of class i. A later uniform draw picks a source tree by bisecting SrcCD.
# Deterministic. Returns (srci, srccd, sptr) — class i occupies k in (sptr[i]|0)+1 .. sptr[i+1].
function dm_src(sp::Int, d, ptr, index, tpa)
    srci = Int32[]; srccd = Float32[]; sptr = zeros(Int, 7)
    k = 0
    @inbounds for i in 0:6
        fst = ptr[sp, i+1, 1]
        if fst > 0
            x = 1f0 / (d[i+1] + DM_TINY)
            y = 0f0
            for j in fst:ptr[sp, i+1, 2]
                k += 1
                y += tpa[index[j]] * x
                push!(srci, index[j]); push!(srccd, y)
            end
        end
        sptr[i+1] = k
    end
    return (srci, srccd, sptr)
end

# --- DMTLST (dmtlst.f) — the list of TARGET trees of species `sp` with DMR `tdmr` that have
# positive expansion (PROB>0), via the Ptr range. Deterministic. Returns the vector of tree
# record indices (Fortran TLst(0)=count → Julia length(tlst)).
function dm_tlst(sp::Int, tdmr::Int, ptr, index, tpa)
    tlst = Int32[]
    fst = ptr[sp, tdmr+1, 1]
    if fst > 0
        @inbounds for i in fst:ptr[sp, tdmr+1, 2]
            j = index[i]
            tpa[j] > 0f0 && push!(tlst, j)
        end
    end
    return tlst
end

# --- DMSAMP (dmsamp.f) — how many SOURCE trees of a given class to place in a target's sampling
# ring: draw a uniform, bisect the neighbour CDF `cnb` to get BigS = total trees in the ring, then
# binomially thin BigS by x = Prop·(D/TotD) (the source-class fraction) via BigS more draws.
# Stochastic (1 + BigS DMRANN draws). Returns S. `cnb` is the (cnb,End) from dm_nb (0-based cnb[j+1]).
function dm_samp!(ms::MistletoeState, totd::Float32, d::Float32, cnb, prop::Float32)
    rnd = dm_rann!(ms)
    bigs = 0
    @inbounds for j in 0:DM_DSTLEN
        if rnd <= cnb[j+1]
            bigs = j; break
        end
    end
    s = 0
    if bigs > 0
        x = prop * (d / totd)
        for _ in 1:bigs
            dm_rann!(ms) <= x && (s += 1)
        end
    end
    return s
end

# --- DMSLST (dmslst.f) — select `n` SOURCE trees of DMR class `dmrcls` (with replacement, weighted
# by the cumulative SrcCD within the class from dm_src), deduped into (tree, count) pairs. For each
# of n draws, bisect SCD[pFrst..pLast] to pick a tree; accumulate its occurrence count. Stochastic
# (n DMRANN draws). Returns (idxs, knts) — m unique sources = length(idxs). sptr is dm_src's length-7
# vector (sptr[c+1] = end offset of class c).
function dm_slst!(ms::MistletoeState, dmrcls::Int, n::Int, sind, scd, sptr)
    pfrst = dmrcls == 0 ? 1 : sptr[dmrcls] + 1
    plast = sptr[dmrcls+1]
    idxs = Int32[]; knts = Int32[]
    for _ in 1:n
        rnd = dm_rann!(ms)
        @inbounds for j in pfrst:plast
            if rnd <= scd[j]
                tree = sind[j]
                pos = findfirst(==(tree), idxs)
                if pos === nothing
                    push!(idxs, tree); push!(knts, Int32(1))
                else
                    knts[pos] += Int32(1)
                end
                break
            end
        end
    end
    return (idxs, knts)
end

# --- DMBSHD (dmbshd.f) — decode the encoded Shd1 trajectory data for grid cell (iz, ix) into the
# CShd trajectory list: for each of `vcnt` trajectories, a header (weight at cs[i,0,XX], StrtVl,
# EndVal=Shd1/2) then (x,z) pairs cs[i,j,XX]/cs[i,j,ZZ] for j=StrtVl..EndVal. Deterministic (a pure
# table decode). Returns (cs[traj, j+1 (j 0:28), {XX,ZZ}], vlen[traj], vcnt). ShdPtr(iz,ix)=0 → vcnt 0.
function dm_bshd(iz::Int, ix::Int)
    cs = zeros(Int32, DM_MXTRAJ, 29, 2)
    vlen = zeros(Int, DM_MXTRAJ)
    ptr = Int(DM_SHDPTR[iz, ix])
    vcnt = 0
    if ptr != 0 && 1 <= ptr < DM_TOP1
        vcnt = Int(DM_SHD1[ptr]); ptr += 1
        @inbounds for i in 1:vcnt
            cs[i, 0+1, DM_XX] = DM_SHD1[ptr]; ptr += 1
            strtvl = Int(DM_SHD1[ptr]); ptr += 1
            endval = Int(DM_SHD1[ptr]) ÷ 2; ptr += 1
            for j in strtvl:endval
                cs[i, j+1, DM_XX] = DM_SHD1[ptr]; ptr += 1
                cs[i, j+1, DM_ZZ] = DM_SHD1[ptr]; ptr += 1
            end
            vlen[i] = endval
        end
    end
    return (cs, vlen, vcnt)
end

# --- DMCYCL core (dmcycl.f:405-510) — one YEAR of life-history compartment advance for a crown
# third. Base BC/YSM path (all biocontrol pools zero): forward cascade ImmLat=xImm·fprop2
# (immature→latent), LatAct=xLat·fprop, SprAct=xSpr·fprop (latent/suppressed→active), ActSpr=
# xAct·bprop (active→suppressed); subtract from source + add to destination; then survival ·spsurv.
# With the default DMLtRx curves + ALGSLP clamp, fprop=1/bprop=0 for every crown third. The New-
# infection intake (xImm+=New) and the DMCAP saturation happen in the driver loop (need NewSpr/
# NewInt/TVol). Deterministic. Returns the advanced (imm, lat, spr, act, ded).
@inline function dm_cycl_advance(imm::Float32, lat::Float32, spr::Float32, act::Float32, ded::Float32,
                                 fprop::Float32, bprop::Float32, fprop2::Float32, spsurv::Float32)
    immlat = imm * fprop2
    latact = lat * fprop
    spract = spr * fprop
    actspr = act * bprop
    imm -= immlat; lat -= latact; spr -= spract; act -= actspr
    lat += immlat; spr += actspr; act += latact + spract
    imm *= spsurv; lat *= spsurv; spr *= spsurv; act *= spsurv; ded *= spsurv
    return (imm, lat, spr, act, ded)
end
const DM_FLWR = 4      # DMFLWR default: years-to-flower → FProp2 = 1/DMFLWR = 0.25 (dminitbc.f:267)
const DM_CAP  = 3.0f0  # DMCAP default: per-crown-third infection carrying capacity (dminitbc.f:268)

# --- DMFSHD (dmfshd.f) — the per-MESH-band canopy shade field Shade[] that DMADLV reads. STOCHASTIC:
# for each canopy band, simulate each tree's expected count (PROB·SQM2AC·Grid²) at random (x,y) on a
# 121×121 grid (DMRANN Poisson), painting a disk of the crown radius with the species opacity (max,
# not additive); Shade = mean opacity over the inner 100×100 (11..110), then per-MESH 1−(1−sh)^MESH.
# Uses the DM RNG (realization straddle vs live). Returns Shade[1:MXHT].
function dm_fshd!(ms::MistletoeState, species, prob, itrn::Int)
    Grid = 121; LowIn = 11; HighIn = 110; Cells = 100f0
    shade = zeros(Float32, DM_MXHT)
    shdlst = Int[]
    @inbounds for u in 1:DM_MXHT
        for v in 1:itrn                                            # add band u if ANY tree has canopy there
            if prob[v] > 0.01f0 && ms.dmrdmx[v, u, DM_RADIUS] > 0.1f0
                push!(shdlst, u); break
            end
        end
    end
    slstln = min(length(shdlst), DM_MXHT)
    g = zeros(Float32, Grid, Grid)
    @inbounds for uu in 1:slstln
        v = shdlst[uu]
        fill!(g, 0f0)
        for i in 1:itrn
            rad = ms.dmrdmx[i, v, DM_RADIUS] * Float32(DM_MESH)
            rad > 0f0 || continue
            tnumbr = prob[i] * _DM_SQM2AC * Float32(Grid)^2
            n = trunc(Int, tnumbr)
            (dm_rann!(ms) <= (tnumbr - n)) && (n += 1)
            opq = ms.opaq[Int(species[i])]                         # DMOPAQ (raw opacity)
            for _ in 1:n
                x = Float32(trunc(Int, dm_rann!(ms) * Grid) + 1)
                y = Float32(trunc(Int, dm_rann!(ms) * Grid) + 1)
                for s in trunc(Int, x-rad):trunc(Int, x+rad), t in trunc(Int, y-rad):trunc(Int, y+rad)
                    if 0 < s <= Grid && 0 < t <= Grid
                        d = sqrt((Float32(s)-x)^2 + (Float32(t)-y)^2)
                        (d <= rad && g[s, t] < opq) && (g[s, t] = opq)
                    end
                end
            end
        end
        sm = 0f0
        for i in LowIn:HighIn, j in LowIn:HighIn
            sm += g[i, j]
        end
        sh = sm / Cells^2
        sh = 1f0 - (1f0 - sh)^DM_MESH
        shade[v] = clamp(sh, 0f0, 1f0)
    end
    return shade
end

# --- DMADLV (dmadlv.f) — accumulate one "level" of the spread field (SFld, to the target) and
# intensification field (IFld, self) from an infected source `srcind` at MESH height `mshht`,
# infection `level`, `cnt` copies, to a target at MESH `dist`. Walks each DMBSHD-decoded seed
# trajectory: inside the source crown (x≤Rad=DMRDMX radius) → IFld += VecWt·Op (opacity capture);
# at the target distance (CShd x == dist) → SFld += Cnt·VecWt·Op; else en-route shading loss
# VecWt −= VecWt·shade[h]. Op = DMOPQ2 = 1−(1−opaq)^MESH (per-MESH-cell opacity, dmtreg.f:229).
# Deterministic given `shade`. II=0/EB=0 (the DMTREG driver values) ⇒ Shd/Shd0 collapse to shade[h].
function dm_adlv!(ms::MistletoeState, species, srcind::Int, cnt::Int,
                  sfld::AbstractVector{Float32}, ifld::AbstractVector{Float32},
                  mshht::Int, dist::Int, level::Float32, shade::AbstractVector{Float32})
    op = 1f0 - (1f0 - ms.opaq[Int(species[srcind])])^DM_MESH        # DMOPQ2
    lszind = max(1, DM_ORIGIN - mshht + 1)                          # source z-index start
    lfzind = max(1, mshht - DM_ORIGIN + 1)                          # field z-index start
    hfzind = min(DM_MXHT, mshht - DM_ORIGIN + DM_MXTHRZ)            # field z-index end
    lsxind = dist                                                  # HSXInd=Dist ⇒ single x
    u = lszind
    @inbounds for i in lfzind:hfzind
        v = lsxind
        for j in lsxind:lsxind                                      # LFXInd..HFXInd = single (dist)
            if 1 <= u <= DM_MXTHRZ && 1 <= v <= DM_MXTHRX           # guard ShdPtr bounds (inert when in range)
                cs, veclen, n = dm_bshd(u, v)
                for k in 1:n
                    vecwt = Float32(cs[k, 1, DM_XX]) * level        # CShd(k,0,XX)
                    xlast = 0f0; hlast = 0
                    for m in 1:veclen[k]
                        h = mshht + Int(cs[k, m+1, DM_ZZ]) - DM_ORIGIN
                        x = Float32(cs[k, m+1, DM_XX]); xlast = x; hlast = h
                        if 1 <= h <= DM_MXHT
                            rad = ms.dmrdmx[srcind, h, DM_RADIUS]
                            y = x - rad
                            if x <= rad
                                loss = vecwt*op; vecwt -= loss; ifld[h] += loss
                            elseif 0f0 < y < 1f0
                                loss = vecwt*op*y; vecwt -= loss; ifld[h] += loss
                            elseif Int(cs[k, m+1, DM_XX]) == dist
                                loss = vecwt*op; sfld[h] += Float32(cnt)*loss
                            else
                                loss = vecwt*shade[h]; vecwt -= loss
                            end
                        end
                        if m == veclen[k] && i == lfzind && lfzind > 1   # edge downward loop
                            for w in i:-1:1
                                rad = ms.dmrdmx[srcind, w, DM_RADIUS]
                                y = xlast - rad
                                if y <= rad
                                    loss = vecwt*op; vecwt -= loss; ifld[w] += loss
                                elseif 0f0 < y < 1f0
                                    loss = vecwt*op*y; vecwt -= loss; ifld[w] += loss
                                elseif Int(cs[k, m+1, DM_XX]) == dist
                                    loss = vecwt*op; sfld[w] += Float32(cnt)*loss
                                else
                                    loss = vecwt*shade[hlast]; vecwt -= loss
                                end
                            end
                        end
                    end
                end
            end
            v += 1
        end
        u += 1
    end
    return nothing
end

# --- SF autocorrelation scaling matrix (dminitbc.f:190-203) — SF[diff,ring] =
# exp(diff·DMALPH · exp(Dstnce[ring]·DMBETA)); reweights source density by the DMR
# difference between source and target class (spatial autocorrelation). DMALPH default
# −0.5, DMBETA default 0.0 (DMAUTO keyword overrides; −999 = unset → default).
function dm_compute_sf!(ms::MistletoeState)
    alph = ms.dmalpha == -999f0 ? -0.50f0 : ms.dmalpha
    beta = ms.dmbeta  == -999f0 ?  0.0f0  : ms.dmbeta
    size(ms.sf) == (7, DM_MXTHRX) || (ms.sf = Matrix{Float32}(undef, 7, DM_MXTHRX))
    @inbounds for j in 1:DM_MXTHRX
        tmp = exp(DM_DSTNCE[j] * beta)
        for i in 0:6
            ms.sf[i+1, j] = exp(Float32(i) * alph * tmp)
        end
    end
    return ms
end

# --- DMAUTO (dmauto.f) — reweight the per-DMR-class source density `d[0:6]` for a target of
# DMR `trgdmr` in ring `rq` by the autocorrelation matrix SF (via the DMR difference |trgdmr−i|),
# preserving the total density DTot. Deterministic. Returns S[0:6] (as 1-based S[i+1]).
function dm_auto(ms::MistletoeState, trgdmr::Int, rq::Int, d)
    s = zeros(Float32, 7)
    dtot = 0f0
    @inbounds for i in 0:6
        dtot += d[i+1]
    end
    dtot <= 0f0 && return s
    ds = 0f0
    @inbounds for i in 0:6
        ds += d[i+1] * ms.sf[abs(trgdmr - i) + 1, rq]   # DMRDFF(i,trgdmr) = |trgdmr-i|
    end
    dtot /= ds
    @inbounds for i in 0:6
        s[i+1] = dtot * d[i+1] * ms.sf[abs(trgdmr - i) + 1, rq]
    end
    return s
end

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
    dm_compute_sf!(ms)
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

# --- C4 neighbour build: DMNB (dmnb.f) — the cumulative distribution of the number of SOURCE
# trees in sampling-ring `rq`'s annulus (outer disk rq minus inner disk rq-1), for a source
# density `d` (trees/acre). Builds each disk's count PDF via BNDIST (mean = CrArea·d, var =
# DMCLMP·mean), conditions each on ≥1 tree, convolves outer⊛inner into the annulus count, then
# cumulates + normalizes. Deterministic (no RNG). Returns CNB (0-based via CNB[k+1]); only
# 0..End meaningful. Faithful to dmnb.f including the j==0 degenerate inner disk (target only).
function dm_nb(rq::Int, d::Float32, dmclmp::Float32)
    nb = zeros(Float32, DM_DSTLEN + 1, 2)   # nb[count+1, disk]
    emark = Int[1, 1]
    cur = 1
    for j in (rq-1):rq
        if j == 0
            nb[2, cur] = 1f0                  # NB(1,1)=1.0: inner disk = the target tree only, P(1)=1
            emark[cur] = 1
        else
            mu  = DM_CRAREA[j] * d
            var = dmclmp * mu
            pdf, endpos, _ = dm_bndist(mu, var, DM_DSTLEN)
            @inbounds for t in 1:DM_DSTLEN+1
                nb[t, cur] = pdf[t]
            end
            emark[cur] = endpos
        end
        cur = 2                               # (prv stays 1)
    end
    prv = 1; cur = 2
    # condition each disk's distribution on ≥1 tree (drop the P(0) mass, renormalize)
    for i in 1:2
        x = 1f0 / (1f0 - nb[1, i])            # nb[1,i] = P(count=0)
        nb[1, i] = 0f0
        @inbounds for j in 1:emark[i]
            nb[j+1, i] *= x
        end
    end
    cnb = zeros(Float32, DM_DSTLEN + 1)
    topend = emark[cur]
    @inbounds for k in 0:topend               # annulus count = outer(cur) ⊛ inner(prv)
        x = 0f0
        for j in k:min(k + emark[prv], DM_DSTLEN)
            yd = Float64(nb[j+1, cur]) * Float64(nb[(j-k)+1, prv])
            yd > 1.0e-25 && (x += Float32(yd))
        end
        cnb[k+1] = x
    end
    @inbounds for k in 1:topend               # cumulate
        cnb[k+1] += cnb[k]
    end
    x = 1f0 / cnb[topend+1]                    # normalize
    @inbounds for k in 0:topend
        cnb[k+1] *= x
    end
    return (cnb, topend)
end
