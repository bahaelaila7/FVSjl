# =============================================================================
# Douglas-fir Tussock Moth (DFTM) defoliator model — dftm/*.f
# =============================================================================
# Second port of the FVS insect/pathogen event-extension family (after DFB).
# UNLIKE DFB (a stochastic-draw + WK2-mortality model), DFTM is a full
# differential-equation *defoliator* population model: an upper (regional) and
# lower (tree-class) module of coupled G/F state equations integrated over a
# fixed 5-year outbreak base period, feeding percent-defoliation → foliage
# biomass → growth-loss + top-kill + mortality back into the FVS treelist.
#
# FVS structure (dftm/, 31 files):
#   * DFTMIN  (dftmin.f)  — the DFTM keyword-block reader (keywds.f option 7,
#                           analogous to DFB=100). 45 sub-keywords → the TMCOM1
#                           /UPPER/LOWER/TMEGGS commons. Ported as `kw_dftmin!`.
#   * TMINIT  (tminit.f + tminit{ec,em,so,tt}.f) — per-stand init of the DFTM
#                           defaults (the B0/R0/B1 model-coefficient arrays, egg
#                           allocations, biomass params, IDFCOD/IGFCOD species
#                           codes). The variant files differ ONLY in IGFCOD
#                           (grand-/white-fir species number). Ported as
#                           `dftm_defaults!` + `dftm_igfcod`.
#   * TMRANN  (tmrann.f)  — the model's OWN Lehmer/MINSTD LCG: S1 = mod(16807·S0,
#                           2147483647); u = S1/2^31. Seed default 55329 (odd).
#                           IDENTICAL algorithm to DFB's DFBRAN. Entry TMRNSD
#                           reseeds (forces odd). Ported as `dftm_rand!` /
#                           `dftm_seed!`. NEVER FFI'd.
#   * DFTMGO  (dftmgo.f)  — per-cycle outbreak gate (grincr.f:402): TMOTPR prob,
#                           OPGET2(810/811) activity match, ITMETH branch
#                           (1=MANSTART deterministic / 2=CRTSTART prob>0.5 /
#                           3=RANSTART TMRANN<prob), LDF/LGF host presence, then
#                           INSCYC forces a 5-yr (TMBASE) cycle. Sets L=go.
#   * TMOTPR  (tmotpr.f)  — DETERMINISTIC stand-outbreak conditional probability
#                           (3 IPRBMT methods: Heller elev/slope/aspect/topo/
#                           crown; Mika-Moore topo/BA/%GF/ash). Ported as
#                           `dftm_otpr` (pure, unit-tested; NOT engine-wired).
#   * TMSCHD  (tmschd.f)  — RANSCHED/MANSCHED regional-outbreak auto-scheduler
#                           (TMOPS → TMSCHD, before the cycle loop).
#   * TMCOUP  (tmcoup.f, 33 KB) — the coupler (gradd.f:103): classifies host
#                           trees into DFTM size/foliage classes, runs DFTMOD,
#                           maps percent-defoliation → growth-loss + mortality.
#   * DFTMOD  (dftmod.f)  — the population-dynamics integrator: for each phase
#                           K(1..3) and tree class, G0COMP/Z1COMP/GFCOMP/Y1COMP/
#                           Y0COMP advance the X0/X1 state vectors (10 sub-steps),
#                           accumulating DPCENT percent-branch-defoliation.
#   * g0comp/gfcomp/y0comp/y1comp/z1comp/uv1/uv2/dfole8 — the G/F/Y function
#                           kernels of the integrator.
#   * tmbmas/tmbchl/biomas — foliage-biomass assignment (4 IBMTYP methods).
#   * garbel/grclas/grpsum/redist/tmdam/g0comp — classification, redistribution,
#                           damage-code tagging, reporting.
#   * tminitec/em/so/tt.f — variant IGFCOD tables (EC=6, EM=9, SO=12, TT=9,
#                           base/IE/BM/UT=4).
#
# LINKAGE: like DFB, NO shipped binary runs DFTM — every FVS*_buildDir links the
# base/exdftm.f NO-OP stub (its entries TMINIT/DFTMIN/DFTMGO/TMCOUP/TMOPS/…
# all RETURN). A DFTM oracle must be relinked = a variant's .o set with the real
# dftm/*.o swapped in for exdftm.o (see scratchpad/dftm/build_ie_dftm.sh).
#
# THIS PORT (chunk 0 — keyword reader + RNG, INERT seam):
#   * `kw_dftmin!`  — faithful DFTMIN block reader → `s.dftm` (DftmState).
#   * `dftm_rand!` / `dftm_seed!` — the TMRANN LCG + TMRNSD reseed (bit-exact).
#   * `dftm_defaults!` — TMINIT defaults; `dftm_igfcod` variant crosswalk.
#   * `dftm_otpr`  — the deterministic TMOTPR outbreak probability (ported +
#                    unit-tested; NOT engine-wired — awaits the DFTMGO seam).
#   NO per-cycle engine seam is wired: a stand carrying a DftmState projects
#   BYTE-IDENTICALLY to one without DFTM (validated by the full FVSjl suite +
#   the relinked FVSie_dftm oracle running DFTM-off ≡ stock). The population
#   dynamics (TMCOUP/DFTMOD) are the next, much larger, chunk.
#
# Float32 discipline: DFTM's Fortran is REAL (Float32). Coefficient arrays and
# every ported arithmetic result are Float32 to stay bit-identical to the oracle.
# =============================================================================

# -----------------------------------------------------------------------------
# TMINIT (tminit.f) DATA — the model-coefficient defaults (Blue-Mountains site).
# -----------------------------------------------------------------------------
# B0(66) /UPPER/, R0(18) /UPPER/, B1(36) /LOWER/ — the DFTMOD integrator's upper
# and lower module parameters. Transcribed in DATA order.
const DFTM_B0_DEFAULT = Float32[
    0.0,   0.0,        0.0,   0.0,   0.0,   0.0,
    0.0,   0.0,        0.0,   0.0,   0.001, 0.001,
    0.002, 0.003,      0.006, 0.013, 0.035, 0.028,
    0.025, 0.028,      0.031, 0.034, 0.035, 0.035,
    0.0,   0.0,        0.0,   0.0,   0.0,   0.0,
    0.0,   0.0,        0.0,   0.0,   0.001, 0.001,
    0.001, 0.002,      0.003, 0.010, 0.016, 0.042,
    0.005, 0.006,      0.007, 0.021, 0.033, 0.056,
    0.50,  0.62,       0.75,  0.80,  200.0, 200.0,
    150.0, 150.0,      0.0,   0.50,  0.60,  0.85,
    0.90,  0.25,       0.0,   0.0,   0.0,   0.0,
]
const DFTM_R0_DEFAULT = Float32[
    0.92, 0.60, 0.07, 0.0,  0.0,  0.0,
    0.95, 0.7,  0.1,  0.02, 0.0,  0.0,
    0.02, 0.02, 0.02, 0.02, 0.02, 0.02,
]
const DFTM_B1_DEFAULT = Float32[
    0.02,  0.02,        0.02,  0.02,   0.02,   0.02,
    5.4,   6.25,        6.25,  2.71,   2.27,   2.2,
    5.4,   6.25,        6.25,  3.69,   3.29,   3.2,
    1.19,  0.081,       0.1,   0.1147, 0.0886, 0.0625,
    0.197, 52.0626838,  0.02,  0.01,   0.0,    0.50,
    0.0,   0.0,         0.0,   0.0,    0.0,    0.0,
]

const DFTM_DEFAULT_SEED = 55329.0f0   # TMINIT: TMRNSD default seed (odd)

"""
    DftmState

Douglas-fir Tussock Moth model state — the FVS `TMCOM1`/`UPPER`/`LOWER`/`TMEGGS`
commons. Populated by `kw_dftmin!` and seeded with the `TMINIT` defaults by
`dftm_defaults!`. `active` mirrors "a DFTM keyword block was read" (LDFTM); the
future per-cycle engine seam is gated on it, so a stand without a DFTM block
projects byte-identically. `igfcod` is the variant's grand-fir species number.
"""
mutable struct DftmState <: AbstractDftmState
    active::Bool         # a DFTM block was read (engine seam gate)
    # --- host / classification (TMCOM1) ---
    ldf::Bool            # LDF    — simulate on Douglas-fir (default true)
    lgf::Bool            # LGF    — simulate on grand/white fir (default true)
    lredis::Bool         # LREDIS — larval redistribution on (default true)
    nclas::NTuple{2,Int32}  # NCLAS — requested DF/GF classes (default 20,20)
    tmpn1::Float32       # TMPN1  — proportion of classes by tree-difference (0.5)
    weight::NTuple{2,Float32}  # WEIGHT — %new-foliage / biomass weights (1,1)
    ibmtyp::Int32        # IBMTYP — biomass method 1..4 (default 4)
    idfcod::Int32        # IDFCOD — Douglas-fir species number (default 3)
    igfcod::Int32        # IGFCOD — grand/white-fir species number (variant: 4/6/9/12)
    # --- outbreak timing (TMCOM1) ---
    itmsch::Int32        # ITMSCH — 1=MANSCHED (user dates) / 2=RANSCHED (auto)
    itmeth::Int32        # ITMETH — 1=MANSTART / 2=CRTSTART / 3=RANSTART (default 1)
    itmslv::Int32        # ITMSLV — >0: salvage survivors after each outbreak
    iprbmt::Int32        # IPRBMT — outbreak-probability method 1..3 (default 1)
    tmwait::Int32        # TMWAIT — RANSCHED min wait years (default 30)
    tmpast::Int32        # TMPAST — year of last recorded outbreak (default 1492)
    tmevnt::Float32      # TMEVNT — RANSCHED annual event probability (default 0.1)
    prbscl::Float32      # PRBSCL — outbreak-probability scaling factor (default 1)
    topo::Float32        # TOPO   — topographic position 1..3 (default 1)
    tmashd::Float32      # TMASHD — soil ash depth, inches (default 15.93)
    tmdefl::Float32      # TMDEFL — SALVAGE critical defoliation % (default 50)
    # --- report / debug (TMCOM1) ---
    itmrep::Int32        # ITMREP — report level 0/1/2 (default 2)
    tmdebu::Bool         # TMDEBU — DEBUG
    tmdtre::Bool         # TMDTRE — DEBUTREE
    lpunch::Bool         # LPUNCH — PUNCH parameter dump
    jotmdk::Int32        # JOTMDK — PUNCH unit
    lnpv2::Bool          # LNPV2  — NPV2 nucleopolyhedrosis-virus override
    lnpv3::Bool          # LNPV3  — NPV3 override
    lransd::Bool         # LRANSD — RANNSEED reseed applied
    ltmprm::Bool         # LTMPRM — a TMPARMS sub-model parameter was redefined
    lchem::Bool          # LCHEM  — CHEMICAL treatment applied
    # --- egg allocation (TMEGGS) ---
    iegtyp::Int32        # IEGTYP — 1=RANLARVA (REGG) / 2=DETLARVA (DEGG)
    dfegg::NTuple{3,Float32}   # DFEGG  — DETLARVA DF thirds (11,9,7)
    gfegg::NTuple{3,Float32}   # GFEGG  — DETLARVA GF thirds (15,10,7)
    dfregg::NTuple{3,Float32}  # DFREGG — RANLARVA DF mean/within/between (9,2,0)
    gfregg::NTuple{3,Float32}  # GFREGG — RANLARVA GF mean/within/between (11,3,0)
    # --- biomass (TMCOM1) ---
    dffbio::NTuple{2,Float32}  # DFFBIO — DF foliage biomass mean/sd (213.8,64.2)
    gffbio::NTuple{2,Float32}  # GFFBIO — GF foliage biomass mean/sd (227,63.7)
    dfpnew::NTuple{2,Float32}  # DFPNEW — DF %new-foliage mean/sd (26.9,12.6)
    gfpnew::NTuple{2,Float32}  # GFPNEW — GF %new-foliage mean/sd (35.2,7.3)
    # --- integrator coefficients (UPPER/LOWER) ---
    b0::Vector{Float32}  # B0(66)  /UPPER/ — upper-module params
    r0::Vector{Float32}  # R0(18)  /UPPER/ — upper-module rates
    b1::Vector{Float32}  # B1(36)  /LOWER/ — lower-module params
    # --- scheduling (TMCOM1) ---
    mansched_years::Vector{Int32}  # MANSCHED IDT dates (OPNEW 810)
    # --- TMRANN LCG state (COMMON in tmrann.f) ---
    rng_s0::Float64      # S0 — current generator state (double, exact mod)
    rng_ss::Float32      # SS — the reseed default (TMRNSD LSET=false resets to it)
end

"""
    dftm_defaults!(variant) -> DftmState

FVS `TMINIT` (tminit.f): the DFTM per-stand run-time defaults, with the variant's
`IGFCOD` grand-fir species number. Called when the first DFTM keyword is seen.
"""
function dftm_defaults!(variant)
    return DftmState(
        false,                                  # active
        true, true, true,                       # ldf, lgf, lredis
        (Int32(20), Int32(20)), 0.5f0,          # nclas, tmpn1
        (1.0f0, 1.0f0), Int32(4),               # weight, ibmtyp
        Int32(3), Int32(dftm_igfcod(variant)),  # idfcod, igfcod
        Int32(1), Int32(1), Int32(0), Int32(1), # itmsch, itmeth, itmslv, iprbmt
        Int32(30), Int32(1492),                 # tmwait, tmpast
        0.1f0, 1.0f0, 1.0f0, 15.93f0, 50.0f0,   # tmevnt, prbscl, topo, tmashd, tmdefl
        Int32(2), false, false, false, Int32(0),# itmrep, tmdebu, tmdtre, lpunch, jotmdk
        false, false, false, false, false,      # lnpv2, lnpv3, lransd, ltmprm, lchem
        Int32(1),                               # iegtyp
        (11.0f0, 9.0f0, 7.0f0), (15.0f0, 10.0f0, 7.0f0),   # dfegg, gfegg
        (9.0f0, 2.0f0, 0.0f0), (11.0f0, 3.0f0, 0.0f0),     # dfregg, gfregg
        (213.8f0, 64.2f0), (227.0f0, 63.7f0),   # dffbio, gffbio
        (26.9f0, 12.6f0), (35.2f0, 7.3f0),      # dfpnew, gfpnew
        copy(DFTM_B0_DEFAULT), copy(DFTM_R0_DEFAULT), copy(DFTM_B1_DEFAULT),
        Int32[],                                # mansched_years
        Float64(DFTM_DEFAULT_SEED), DFTM_DEFAULT_SEED,     # rng_s0, rng_ss
    )
end

# -----------------------------------------------------------------------------
# TMINIT variant IGFCOD (grand/white-fir species number) — tminit{,ec,em,so,tt}.f
# -----------------------------------------------------------------------------
"""
    dftm_igfcod(variant) -> Int

The FVS species number of grand/white fir for a DFTM-linked variant (the only
value that differs across `tminit{,ec,em,so,tt}.f`): base/IE/BM/UT/SO(base)=4,
EC=6, TT/EM=9, SO(tminitso)=12. IDFCOD (Douglas-fir) is 3 in every DFTM variant.
Returns the base default (4) for variants without a specific table.
"""
function dftm_igfcod(variant)::Int
    variant isa EastCascades  && return 6
    variant isa Teton         && return 9
    return 4
end

# -----------------------------------------------------------------------------
# TMRANN (tmrann.f) — the DFTM model's own double-precision Lehmer/MINSTD LCG.
# -----------------------------------------------------------------------------
"""
    dftm_rand!(d) -> Float32

FVS `TMRANN`: `S1 = DMOD(16807·S0, 2147483647)`; `SEL = REAL(S1/2147483648)`; `S0 = S1`.
`S0`/`S1` are DOUBLE, so the modular step is exact integer arithmetic; only the
returned uniform `SEL` truncates to Float32. IDENTICAL to DFB's DFBRAN. Default
seed 55329. NEVER FFI'd.
"""
@inline function dftm_rand!(d::DftmState)::Float32
    s1 = rem(16807.0 * d.rng_s0, 2147483647.0)   # DMOD (exact, values < 2^31)
    d.rng_s0 = s1
    return Float32(s1 / 2147483648.0)            # SEL = REAL(S1 / 2^31)
end

"""
    dftm_seed!(d, seed, present)

FVS `TMRNSD(LSET, SEED)` (tmrann.f ENTRY). RANNSEED with a value (`present`):
force an ODD seed, store it as both the working state `S0` and the reseed default
`SS`. RANNSEED with no value (`present=false`): reset `S0` to the current `SS`.
"""
function dftm_seed!(d::DftmState, seed::Float32, present::Bool)
    if !present
        d.rng_s0 = Float64(d.rng_ss)             # LSET=false: SEED=SS; S0=SS
        return nothing
    end
    (seed % 2.0f0 == 0.0f0) && (seed += 1.0f0)   # AMOD(SEED,2)==0 → SEED+1 (odd)
    d.rng_ss = seed
    d.rng_s0 = Float64(seed)
    return nothing
end

# -----------------------------------------------------------------------------
# TMOTPR (tmotpr.f) — DETERMINISTIC stand-outbreak conditional probability.
# -----------------------------------------------------------------------------
"""
    dftm_otpr(iprbmt, topo, tmashd; elev, slope, aspect, relden, reldsp3, reldsp4,
              tprob, ba, pgfba, has_df) -> Float32 (PROTBK)

FVS `TMOTPR`: the conditional probability that a stand sustains a DFTM outbreak
given a regional one. Three `IPRBMT` methods:

* method 1 (Heller): a logistic in ELEV·100, SLOPE·100 (+aspect terms), TOPO,
  STDCLO (= min(100, RELDEN/150·100)), HOST (= (RELDSP3+RELDSP4)/RELDEN·100), and
  AVCRDI (= 2·√(435.6·(RELDEN/TPROB)/π)). Returns 0 when the stand has no DF
  (`has_df` false). `elev` is in hundreds of feet, `slope` scaled 0-1, `aspect`
  in radians — exactly the PLOT-common units.
* method 2 (Mika-Moore): logistic in (TOPO-1), TMASHD, PGFBA (proportion of BA in
  grand fir), ln(BA).
* method 3: like method 2 without the ash term.

`pgfba` is Σ(0.005454154·DBH²·PROB) over the grand-fir records / BA. Float32
throughout to mirror the REAL Fortran. DETERMINISTIC — no TMRANN draw. NOT yet
engine-wired (awaits the DFTMGO seam that supplies RELDEN/RELDSP/TPROB/PGFBA);
unit-tested against hand-computed logistics.
"""
function dftm_otpr(iprbmt::Integer, topo::Float32, tmashd::Float32;
                   elev::Float32 = 0.0f0, slope::Float32 = 0.0f0, aspect::Float32 = 0.0f0,
                   relden::Float32 = 0.0f0, reldsp3::Float32 = 0.0f0, reldsp4::Float32 = 0.0f0,
                   tprob::Float32 = 1.0f0, ba::Float32 = 1.0f0, pgfba::Float32 = 0.0f0,
                   has_df::Bool = true)::Float32
    if iprbmt <= 1
        has_df || return 0.0f0
        stdclo = 100.0f0
        relden < 150.0f0 && (stdclo = (relden / 150.0f0) * 100.0f0)
        host = ((reldsp3 + reldsp4) / relden) * 100.0f0
        avcrdi = 2.0f0 * sqrt((435.6f0 * (relden / tprob)) / 3.141593f0)
        arg = -0.852926f0 -
              (0.00023762f0 * elev * 100.0f0) -
              (0.00052207f0 * slope * 100.0f0) +
              (0.00348125f0 * slope * 100.0f0 * cos(aspect)) +
              (0.00733782f0 * slope * 100.0f0 * sin(aspect)) -
              (0.375624f0 * topo) +
              (0.0204993f0 * stdclo) +
              (0.0186693f0 * host) +
              (0.0168552f0 * avcrdi)
        return 1.0f0 / (1.0f0 + exp(-arg))
    elseif iprbmt == 2
        arg = -12.7746f0 - 1.4764f0 * (topo - 1.0f0) - 0.0992f0 * tmashd +
              4.8334f0 * pgfba + 2.5529f0 * log(ba)
        return 1.0f0 / (1.0f0 + exp(-arg))
    else
        arg = -12.3652f0 - 1.4050f0 * (topo - 1.0f0) +
              3.8230f0 * pgfba + 2.2426f0 * log(ba)
        return 1.0f0 / (1.0f0 + exp(-arg))
    end
end

# -----------------------------------------------------------------------------
# DFTMGO (dftmgo.f) — per-cycle outbreak host-threshold gate (deterministic).
# -----------------------------------------------------------------------------
"""
    dftm_go_gate(df_probs, gf_probs; ldf, lgf, nclas) -> NamedTuple

FVS `DFTMGO` host-presence gate (dftmgo.f, statements 60–210), the DETERMINISTIC
tail of the per-cycle outbreak decision that runs AFTER the ITMETH start test has
already committed to an outbreak. Given the DF and GF host `PROB` (trees/acre)
records — in FVS `IND1`/`ISCT` order — it decides whether the tussock-moth model
can actually run on each host and how many DFTM tree-classes to build:

* `IDF` = number of Douglas-fir records (`ISCT(IDFCOD,2)-ISCT(IDFCOD,1)+1`); 0 when
  the species is absent (`ISCT(IDFCOD,1)==0`).
* `CNTDF` = Σ`PROB` over the DF records (serial add in record order — bit-exact).
* `LDF` stays on only if it was requested (not `NODFRUN`), the species is present,
  `NCLAS(1)>0`, and `CNTDF ≥ 0.01`; otherwise it drops to `.FALSE.` (stmt 90).
* Grand fir mirrors this with `IGFCOD`/`NCLAS(2)`/`CNTGF`.
* `L = LDF .OR. LGF` — the overall go/no-go (a false `L` is the "outbreak can not
  be simulated" warning). `NACLAS(k) = MIN(count, NCLAS(k))` when host k is on.

DETERMINISTIC (no TMRANN draw): validated by g16 dump-replay against the relinked
`FVSie_dftm` (`DBGGO*` dumps). NOT engine-wired — awaits the DFTMGO/TMBMAS/TMCOUP
seam that supplies the per-record host `PROB` from the FVSjl treelist.
"""
function dftm_go_gate(df_probs::AbstractVector{Float32}, gf_probs::AbstractVector{Float32};
                      ldf::Bool = true, lgf::Bool = true,
                      nclas::Tuple{<:Integer,<:Integer} = (20, 20))
    idf = 0; cntdf = 0.0f0; ldf_out = ldf
    if ldf
        if isempty(df_probs)                     # ISCT(IDFCOD,1) == 0
            ldf_out = false
        else
            idf = length(df_probs)               # I2 - I1 + 1
            if nclas[1] <= 0 || idf < 1
                ldf_out = false
            else
                @inbounds for p in df_probs; cntdf += p; end
                cntdf < 0.01f0 && (ldf_out = false)
            end
        end
    else
        ldf_out = false
    end
    igf = 0; cntgf = 0.0f0; lgf_out = lgf
    if lgf
        if isempty(gf_probs)
            lgf_out = false
        else
            igf = length(gf_probs)
            if nclas[2] <= 0 || igf < 1
                lgf_out = false
            else
                @inbounds for p in gf_probs; cntgf += p; end
                cntgf < 0.01f0 && (lgf_out = false)
            end
        end
    else
        lgf_out = false
    end
    l = ldf_out || lgf_out
    naclas1 = ldf_out ? min(idf, Int(nclas[1])) : 0
    naclas2 = lgf_out ? min(igf, Int(nclas[2])) : 0
    return (cntdf = cntdf, cntgf = cntgf, idf = idf, igf = igf,
            naclas = (naclas1, naclas2), ldf = ldf_out, lgf = lgf_out, l = l)
end

# -----------------------------------------------------------------------------
# TMBMAS (tmbmas.f) IBMTYP=2 — DETERMINISTIC Hatch–Mika foliage-biomass regressions.
# -----------------------------------------------------------------------------
# Method 2 is the only fully deterministic IBMTYP (1/3/4 draw TMBCHL normal errors
# off the TMRANN stream). It assigns per host tree the nominal-branch foliage
# biomass FBIOMS (grams) and the %-new-foliage PCNEWF, feeding DFTMOD/TMCOUP.
# `dgi` = DG·SCALDG, the 10-yr-scaled diameter increment (SCALDG = 10/period);
# `cr` = ICR/100 (crown ratio 0–1); `pct` = the tree's BA percentile (PCT).
# Equations from Hatch & Mika (Univ. of Idaho, 1978). Float32 throughout.

"""
    dftm_bmas2_df(slope, aspect, ba, tprob, relden, dbh, ht, dgi, cr, pct) -> (fbioms, pcnewf)

Douglas-fir branch of `TMBMAS` method 2 (tmbmas.f, stmt 8). FBIOMS clamped to
[91,400] g; PCNEWF (×100) clamped to [11,42] %.
"""
@inline function dftm_bmas2_df(slope::Float32, aspect::Float32, ba::Float32,
        tprob::Float32, relden::Float32, dbh::Float32, ht::Float32,
        dgi::Float32, cr::Float32, pct::Float32)::NTuple{2,Float32}
    fb = 195.20482f0 -
         249.78151f0 * slope -
         99.368f0    * slope * cos(aspect) -
         129.92132f0 * slope * sin(aspect) +
         22.0528f0   * dbh -
         2.14018f0   * ht +
         79.4349f0   * cr -
         0.62932f0   * ba
    fb > 400.0f0 && (fb = 400.0f0)
    fb < 91.0f0  && (fb = 91.0f0)
    pn = 0.49738f0 +
         0.08057f0  * slope +
         0.27017f0  * slope * cos(aspect) +
         0.32162f0  * slope * sin(aspect) +
         0.03936f0  * log(pct * (0.31830989f0 * atan(((relden / 100.0f0) - 1.5f0) / 1.7f0) + 0.5f0)) -
         0.0043258f0 * ht +
         0.077044f0 * dgi -
         0.078413f0 * tprob / 100.0f0 +
         0.0048268f0 * (tprob * tprob / 10000.0f0)
    pn *= 100.0f0
    pn > 42.0f0 && (pn = 42.0f0)
    pn < 11.0f0 && (pn = 11.0f0)
    return (fb, pn)
end

"""
    dftm_bmas2_gf(slope, aspect, ba, tprob, relden, dbh, ht, dgi, cr, pct) -> (fbioms, pcnewf)

Grand/white-fir branch of `TMBMAS` method 2 (tmbmas.f, stmt 18). FBIOMS is an
exponential in DBH/HT/DGI clamped to [125,400] g; PCNEWF (×100) clamped to
[15,47] %. `ba`/`tprob` are unused here (kept for a uniform DF/GF signature). The
`exp` may carry a documented ≤1-ULP transcendental straddle vs gfortran's `expf`.
"""
@inline function dftm_bmas2_gf(slope::Float32, aspect::Float32, ba::Float32,
        tprob::Float32, relden::Float32, dbh::Float32, ht::Float32,
        dgi::Float32, cr::Float32, pct::Float32)::NTuple{2,Float32}
    fb = exp(4.70244f0 + 0.15833f0 * dbh - 0.01429f0 * ht + 0.17778f0 * dgi)
    fb > 400.0f0 && (fb = 400.0f0)
    fb < 125.0f0 && (fb = 125.0f0)
    pn = 0.1269f0 -
         0.017636f0 * log(pct * (0.31830989f0 * atan(((relden / 100.0f0) - 1.5f0) / 1.7f0) + 0.5f0)) +
         0.037998f0 * dbh -
         0.0062929f0 * ht +
         0.49115f0  * cr -
         0.59622f0  * slope
    pn *= 100.0f0
    pn > 47.0f0 && (pn = 47.0f0)
    pn < 15.0f0 && (pn = 15.0f0)
    return (fb, pn)
end

# -----------------------------------------------------------------------------
# kw_dftmin! (dftmin.f) — DFTM keyword-block reader (keywds.f option 7).
# -----------------------------------------------------------------------------
"""
    kw_dftmin!(s, rec, kr)

Parse the `DFTM … END` block (dftm/dftmin.f). Faithfully sets the TMCOM1/UPPER/
LOWER/TMEGGS-equivalent state on `s.dftm` and consumes sub-keyword records up to
`END`, exactly as FVS `DFTMIN` does. The block is entered by the top-level DFTM
keyword (keywds option 7). No per-cycle engine seam is wired yet, so this is
INERT: a stand carrying a DftmState still projects byte-identically to one
without DFTM.

Sub-keywords ported (state-setting): END, REPORT, DATELIST(report-only), NODFRUN,
NOGFRUN, MANSCHED, DEBUG, RANNSEED, WEIGHT, NUMCLASS, RANLARVA, DETLARVA,
DFTMECHO(file-open, no state; consumes 1 supplemental record), REDIST, NOREDIST,
NPV2, NPV3, TMPARMS, CHEMICAL, DEBUTREE, PUNCH, BIOMASS, DFBIOMAS, GFBIOMAS,
RANSCHED, MANSTART, ASHDEPTH, CRTSTART, RANSTART, PROBMETH, TOPO, SALVAGE.
"""
function kw_dftmin!(s::StandState, rec, kr::KeywordReader)
    s.dftm === nothing && (s.dftm = dftm_defaults!(s.variant))
    d = s.dftm
    d.active = true

    while true
        r = read_keyword!(kr)
        (r.status == KW_EOF || r.status == KW_STOP) && break
        k = strip(r.name)
        isempty(k) && continue
        if k == "END"                       # option 1
            break
        elseif k == "REPORT"                # option 2
            r.present[1] && (d.itmrep = Int32(trunc(Int, r.values[1])))
        elseif k == "DATELIST"              # option 3 — TMDTLS report, no state
            # report-only
        elseif k == "NODFRUN"               # option 4
            d.ldf = false
        elseif k == "NOGFRUN"               # option 5
            d.lgf = false
        elseif k == "MANSCHED"              # option 6 — schedule regional outbreak at IDT (default 1)
            idt = r.present[1] ? Int32(trunc(Int, r.values[1])) : Int32(1)
            d.itmsch = Int32(1)
            push!(d.mansched_years, idt)
        elseif k == "DEBUG"                 # option 7
            d.tmdebu = true
            d.itmrep = Int32(2)
        elseif k == "RANNSEED"              # option 8 — TMRNSD reseed (forces odd)
            d.lransd = true
            dftm_seed!(d, Float32(r.values[1]), r.present[1])
        elseif k == "WEIGHT"                # option 9
            w1 = r.present[1] ? Float32(r.values[1]) : d.weight[1]
            w2 = r.present[2] ? Float32(r.values[2]) : d.weight[2]
            d.weight = (w1, w2)
        elseif k == "NUMCLASS"              # option 10 — DF/GF class counts (Σ ≤ 100)
            i1 = r.present[1] ? Int32(trunc(Int, r.values[1])) : d.nclas[1]
            i2 = r.present[2] ? Int32(trunc(Int, r.values[2])) : d.nclas[2]
            if (Int(i1) + Int(i2)) <= 100
                d.nclas = (i1, i2)
                r.present[3] && (d.tmpn1 = Float32(r.values[3]))
            end
        elseif k == "RANLARVA"              # option 11 — RANLARVA egg-mass (REGG)
            isp = trunc(Int, r.values[1])
            if 1 <= isp <= 2
                d.iegtyp = Int32(1)
                regg = _dftm_regg(d)                        # [DFREGG(1..3); GFREGG(1..3)]
                for i in 2:4
                    r.present[i] && (regg[3*isp + i - 4] = Float32(r.values[i]))
                end
                _dftm_set_regg!(d, regg)
            end
        elseif k == "DETLARVA"              # option 12 — DETLARVA egg-mass (DEGG)
            isp = trunc(Int, r.values[1])
            if 1 <= isp <= 2
                d.iegtyp = Int32(2)
                degg = _dftm_degg(d)                        # [DFEGG(1..3); GFEGG(1..3)]
                for i in 2:4
                    r.present[i] && (degg[3*isp + i - 4] = Float32(r.values[i]))
                end
                _dftm_set_degg!(d, degg)
            end
        elseif k == "DFTMECHO"             # option 13 — post-processor file; reads 1 supplemental record
            read_raw_line!(kr)             # consume the file-name line (READ(IREAD) PNAME)
        elseif k == "REDIST"               # option 14
            r.present[1] && (d.b0[62] = Float32(r.values[1]))
            d.lredis = d.b0[62] != 0.0f0
        elseif k == "NOREDIST"             # option 15
            d.b0[62] = 0.0f0
            d.lredis = false
        elseif k == "NPV2"                 # option 20 — nucleopolyhedrosis virus, B0(9..12)
            d.lnpv2 = true
            d.b0[9]  = r.present[1] ? Float32(r.values[1]) : 0.036f0
            d.b0[10] = r.present[2] ? Float32(r.values[2]) : 0.039f0
            d.b0[11] = r.present[3] ? Float32(r.values[3]) : 0.042f0
            d.b0[12] = r.present[4] ? Float32(r.values[4]) : 0.072f0
        elseif k == "NPV3"                 # option 21 — B0(15..18)
            d.lnpv3 = true
            d.b0[15] = r.present[1] ? Float32(r.values[1]) : 0.036f0
            d.b0[16] = r.present[2] ? Float32(r.values[2]) : 0.039f0
            d.b0[17] = r.present[3] ? Float32(r.values[3]) : 0.042f0
            d.b0[18] = r.present[4] ? Float32(r.values[4]) : 0.072f0
        elseif k == "TMPARMS"              # option 26 — set one B0/R0/B1 element by (I,J)
            i = trunc(Int, r.values[1]); j = trunc(Int, r.values[2])
            if r.present[1] && r.present[2] && r.present[3] && 1 <= j <= 6
                d.ltmprm = true
                if 2 <= i <= 12
                    d.b0[(i - 2) * 6 + j] = Float32(r.values[3])
                elseif 19 <= i <= 21
                    d.r0[(i - 19) * 6 + j] = Float32(r.values[3])
                elseif 22 <= i <= 25
                    d.b1[(i - 21) * 6 + j] = Float32(r.values[3])
                end
            end
        elseif k == "CHEMICAL"            # option 27 — larvicide efficacy → B0(K2)
            d.lchem = true
            iph = 3; ia = trunc(Int, r.values[1])
            (1 <= ia <= 4) && (iph = ia)
            instar = 4; ib = trunc(Int, r.values[2])
            (1 <= ib <= 6) && (instar = ib)
            effic = 0.95f0
            (r.present[3] && 0.0f0 <= r.values[3] <= 1.0f0) && (effic = Float32(r.values[3]))
            k1 = (iph - 1) * 6 + instar
            k2 = k1 + 24
            k3 = instar + 12
            d.b0[k2] = 1.0f0 - (((1.0f0 - effic)^0.1f0) /
                                ((1.0f0 - d.r0[k3]) * (1.0f0 - d.b0[k1])))
        elseif k == "DEBUTREE"            # option 28
            d.tmdtre = true
        elseif k == "PUNCH"              # option 30
            jo = trunc(Int, r.values[1])
            if jo > 0
                d.lpunch = true
                d.jotmdk = Int32(jo)
            end
        elseif k == "BIOMASS"           # option 32 — foliage-biomass method 1..4
            ib = trunc(Int, r.values[1])
            (1 <= ib <= 4) && (d.ibmtyp = Int32(ib))
        elseif k == "DFBIOMAS"          # option 33
            r.present[1] && (d.dffbio = (Float32(r.values[1]), d.dffbio[2]))
            r.present[2] && (d.dffbio = (d.dffbio[1], Float32(r.values[2])))
            r.present[3] && (d.dfpnew = (Float32(r.values[3]), d.dfpnew[2]))
            r.present[4] && (d.dfpnew = (d.dfpnew[1], Float32(r.values[4])))
        elseif k == "GFBIOMAS"          # option 34
            r.present[1] && (d.gffbio = (Float32(r.values[1]), d.gffbio[2]))
            r.present[2] && (d.gffbio = (d.gffbio[1], Float32(r.values[2])))
            r.present[3] && (d.gfpnew = (Float32(r.values[3]), d.gfpnew[2]))
            r.present[4] && (d.gfpnew = (d.gfpnew[1], Float32(r.values[4])))
        elseif k == "RANSCHED"          # option 36 — auto-scheduled regional outbreaks
            d.itmsch = Int32(2)
            r.present[1] && (d.tmwait = Int32(trunc(Int, r.values[1])))
            r.present[2] && (d.tmevnt = Float32(r.values[2]))
            r.present[3] && (d.tmpast = Int32(trunc(Int, r.values[3])))
        elseif k == "MANSTART"          # option 37
            d.itmeth = Int32(1)
        elseif k == "ASHDEPTH"          # option 38
            r.present[1] && (d.tmashd = Float32(r.values[1]))
        elseif k == "CRTSTART"          # option 40
            d.itmeth = Int32(2)
        elseif k == "RANSTART"          # option 41
            d.itmeth = Int32(3)
        elseif k == "PROBMETH"          # option 42
            ip = trunc(Int, r.values[1])
            (0 < ip <= 3) && (d.iprbmt = Int32(ip))
            r.present[2] && (d.prbscl = Float32(r.values[2]))
        elseif k == "TOPO"              # option 43
            (r.present[1] && r.values[1] <= 3.0f0) && (d.topo = Float32(r.values[1]))
        elseif k == "SALVAGE"           # option 45
            d.itmslv = Int32(1)
            r.present[1] && (d.tmdefl = Float32(r.values[1]))
        else
            # DFTM sub-keyword not recognized — record it so it can't hide as a
            # silent gap (mirrors the RD/DFB keyword-dispatch policy).
            (!isempty(k) && isletter(first(k))) && push!(s.control.unrecognized_keywords, k)
        end
    end
    return nothing
end

# Egg arrays are stored as NTuples (TMEGGS); RANLARVA/DETLARVA index them as the
# 6-long EQUIVALENCEd REGG/DEGG (DFREGG(1..3);GFREGG(1..3) and DFEGG;GFEGG).
@inline _dftm_regg(d::DftmState) = Float32[d.dfregg[1], d.dfregg[2], d.dfregg[3],
                                           d.gfregg[1], d.gfregg[2], d.gfregg[3]]
@inline _dftm_degg(d::DftmState) = Float32[d.dfegg[1], d.dfegg[2], d.dfegg[3],
                                           d.gfegg[1], d.gfegg[2], d.gfegg[3]]
@inline function _dftm_set_regg!(d::DftmState, v::Vector{Float32})
    d.dfregg = (v[1], v[2], v[3]); d.gfregg = (v[4], v[5], v[6]); nothing
end
@inline function _dftm_set_degg!(d::DftmState, v::Vector{Float32})
    d.dfegg = (v[1], v[2], v[3]); d.gfegg = (v[4], v[5], v[6]); nothing
end

# =============================================================================
# DFTMOD (dftmod.f) — the Douglas-fir Tussock Moth population-dynamics integrator.
# =============================================================================
# VALIDATION STATUS (dump-replay vs instrumented FVSie_dftm, dense.key):
#   * DFTMOD → DPCENT  — BIT-EXACT (0 ULP), 17 non-empty classes + the empty-GF
#                        class NaN reset.  Test: test_dftm.jl "DFTMOD integrator".
#   * dftm_tree_defol / dftm_mortality / dftm_dgloss / dftm_htgloss_notopkill /
#     dftm_topkill (K1 leader; K2 crown PCKILL/HTGLOS) — BIT-EXACT (0 ULP).
#       Test: test_dftm.jl "TMCOUP damage functions".
#   * dftm_garbel/_grclas/_grpsum/_iqrsrt classification — BIT-EXACT (0 ULP): the
#     RDPSRT-sorted pointer, the ISC sector pointers (incl. the empty class 11,9
#     and the cross-block sector underflow 10,11), and Z4/Z2/Z3, all 18 classes.
#       Test: test_dftm.jl "GARBEL/GRCLAS classification".
#   * dftm_tmbchl / dftm_alloc_eggs + the full TMRANN-stream ordering — BIT-EXACT
#     (0 ULP): TMBCHL vs the pristine driver; the RANLARVA egg X7 in JCLAS2 order
#     (draws 1-78); and the DO-380 per-tree PRTOPK + K≥2 RANDOM continuation
#     (draws 79-98).  Test: test_dftm.jl "TMRANN stream".
#   * STILL engine-inert (no simulate.jl seam wires DFTM into a projection): the
#     INSCYC cycle-forcing hook + the gated TMCOUP coupling seam + the end-to-end
#     .sum-DELTA are the remaining DFTM work (see the DFTM handoff memory note).
# =============================================================================
# The upper (regional, module S(0)) + lower (per tree-class, module S(1)) coupled
# G/F/Y state equations of Overton–Colbert–White, integrated over the FVS 5-year
# outbreak base period in three phases K(3)=1..3, ten occasions KP=0..9 each. The
# per-class % branch defoliation DPCENT(I,1)=Phase-II, DPCENT(I,2)=Phase-III feeds
# TMCOUP. A faithful line-by-line port of dftmod.f + its kernels (g0comp/z1comp/
# gfcomp/y0comp/y1comp/dfole8/redist/uv1/uv2). Float32 throughout.
#
# STATIC-LOCAL FIDELITY: FVS is compiled `-fno-automatic`, so three GFCOMP locals
# are STATIC (persist across the KP-occasion calls, zero-initialised at run start):
#   * G19TOT — reset when KP==0, accumulated over KP=1..6  (→ G1(19)).
#   * G17TOT — reset when KP==0 or KP==9, accumulated       (→ G1(9)).
#   * PHI    — assigned only when (G1(4)+G1(18))>0; otherwise it LEAKS its prior
#              value into G1(13)=X1(4)·(exp(10·G1(12)·PHI)−1). Carried across the
#              whole integration, initialised 0 (matches -fno-automatic zero-init).
#
# NaN-RESET GUARD (dftmod.f DO 6): after X1(IX)+=F1(IX), `IF (X1(IX).GE.0.0) …`
# is FALSE for NaN, so a NaN state is "RESET TO 0.0" exactly like a negative. An
# empty GF class (PROB 0 ⇒ Z3=FBIOMS=0) drives 0/0 NaNs at KP=8 (G1(1) has 2·Z1(2)
# =2·Z3=0 in a denominator) and a NaN DPCENT=100·(1−X6/Z3)=100·(1−0/0). Both are
# reproduced bit-exactly by feeding Z3=0 for that class; DPCENT there = NaN.
# =============================================================================

"""
    DftmInteg

Mutable working state of the DFTMOD integrator — the FVS `ICOND`/`GPASS`/`DFOL`/
`UPPER`/`LOWER`/`LIMITS` commons, plus the three persistent GFCOMP static locals.
`n` = IFIN (number of DFTM tree classes). Built by `dftm_integ` and advanced in
place; `dpcent[:,1]`/`dpcent[:,2]` hold the Phase-II/III % branch defoliation.
"""
mutable struct DftmInteg
    n::Int
    # ICOND (per class 1..n)
    iz6::Vector{Int}
    z4::Vector{Float32}; z5::Vector{Float32}
    z2::Vector{Float32}; z3::Vector{Float32}
    x5::Vector{Float32}; x6::Vector{Float32}; x7::Vector{Float32}
    g19::Vector{Float32}
    # GPASS
    g2::Vector{Float32}; g3::Vector{Float32}
    # DFOL
    dpcent::Matrix{Float32}; g19out::Matrix{Float32}
    # UPPER
    x0::Vector{Float32}; b0::Vector{Float32}; r0::Vector{Float32}
    # LOWER
    x1::Vector{Float32}; z1::Vector{Float32}; b1::Vector{Float32}; r1::Vector{Float32}
    g1::Vector{Float32}; f1::Vector{Float32}; y1::Vector{Float32}
    # LIMITS
    kk::Vector{Int}; kp::Int; ic::Vector{Int}; icount::Int
    # persistent GFCOMP static locals
    phi::Float32; g19tot::Float32; g17tot::Float32
end

"""
    dftm_integ(iz6, z2, z3, x5, x6, x7, z4, z5, b0, r0, b1) -> DftmInteg

Assemble the DFTMOD working state for `n = length(iz6)` classes from the ICOND
class arrays (IZ6 species 1/2, Z2 %new-foliage, Z3 FBIOMS, X5 new biomass, X6 old
biomass, X7 eggs, Z4 PROB, Z5) and the UPPER/LOWER coefficient arrays B0(66)/
R0(18)/B1(36). Inputs are copied (the integrator mutates X5/X6/X7/B1).
"""
function dftm_integ(iz6, z2, z3, x5, x6, x7, z4, z5,
                    b0::AbstractVector{Float32}, r0::AbstractVector{Float32},
                    b1::AbstractVector{Float32})
    n = length(iz6)
    DftmInteg(n,
        Int[Int(v) for v in iz6], Float32.(z4), Float32.(z5),
        Float32.(z2), Float32.(z3), Float32.(x5), Float32.(x6), Float32.(x7),
        zeros(Float32, n),                       # g19
        zeros(Float32, n), zeros(Float32, n),    # g2, g3
        zeros(Float32, n, 2), zeros(Float32, n, 2),
        zeros(Float32, 4), collect(Float32, b0), collect(Float32, r0),
        zeros(Float32, 4), zeros(Float32, 4), collect(Float32, b1), zeros(Float32, 18),
        zeros(Float32, 19), zeros(Float32, 4), zeros(Float32, 3),
        Int[0, 0, 0], 0, Int[0, 0, 0], 0,
        0.0f0, 0.0f0, 0.0f0)
end

# REDIST (redist.f): egg redistribution over the stand → X0(kp), G2(j).
@inline function _dftm_redist!(g::DftmInteg, kp::Int, inum::Int)
    sumegg = 0.0f0; tstem = 0.0f0
    @inbounds for j in 1:inum
        sumegg += g.x7[j] * g.z4[j] * g.z5[j]
        tstem  += g.z4[j] * g.z5[j]
    end
    g.x0[kp] = sumegg / tstem
    b062 = g.b0[62]
    @inbounds for j in 1:inum
        g.g2[j] = g.x7[j] + b062 * (g.x0[kp] - g.x7[j])
    end
    return nothing
end

# G0COMP (g0comp.f): upper-module G functions; no-op in phase 1.
@inline function _dftm_g0comp!(g::DftmInteg)
    g.kk[3] == g.kk[1] && return nothing         # K(3)==K(1)
    inum = g.ic[2] - g.ic[1] + 1
    kp = g.kk[3] - 1
    _dftm_redist!(g, kp, inum)
    @inbounds for j in 1:inum
        g.g3[j] = min((1.0f0 - g.z2[j] / 100.0f0) * g.z3[j], g.x6[j])   # AMIN1
    end
    kpp = g.kk[3] + 56
    @inbounds for i in 1:inum
        g.x6[i] = g.g3[i]
        g.x7[i] = (1.0f0 - g.b0[kpp]) * g.g2[i]
    end
    return nothing
end

# Z1COMP (z1comp.f): load the Z/X vectors + phase-specific B1/R1 for class ICOUNT.
@inline function _dftm_z1comp!(g::DftmInteg)
    c = g.icount
    k3 = g.kk[3]
    g.z1[3] = Float32(k3)
    ip = k3
    g.z1[1] = g.z2[c]; g.z1[2] = g.z3[c]
    g.x1[1] = g.x5[c]; g.x1[2] = g.x6[c]; g.x1[3] = g.x7[c]
    g.z1[4] = Float32(g.iz6[c]); g.x1[4] = 0.0f0
    g.b1[29] = g.z1[3]
    g.b1[34] = g.b0[ip + 48]
    g.b1[35] = g.b0[ip + 52]
    g.b1[36] = g.b0[ip + 62]
    izc = g.iz6[c]
    @inbounds for j in 1:6
        i1 = j + (ip - 1) * 6
        i2 = j + (ip + 3) * 6
        i3 = j + (izc - 1) * 6
        g.b1[j]      = g.r0[j + 12]
        g.r1[j]      = g.b0[i1]
        g.r1[j + 6]  = g.b0[i2]
        g.r1[j + 12] = g.r0[i3]
    end
    return nothing
end

# UV1 (uv1.f): lower-module upper-integral kernel.
@inline function _dftm_uv1(kp::Int, x1::Vector{Float32}, b1::Vector{Float32},
                           alpha::Float32, g1::Vector{Float32})::Float32
    g1[12] <= 0.000001f0 && return 0.0f0
    (kp == 0 || kp > 6) && return 0.0f0
    u = x1[3] * x1[4]
    i = kp + 6
    eg = exp(g1[12])
    t = u
    t *= (b1[19] + b1[20] / g1[12])
    t *= (eg - 1.0f0)
    t *= 0.001f0
    t /= b1[21]
    t *= b1[i]
    t *= (1.0f0 - (g1[6] * eg)^alpha)
    t /= (1.0f0 - g1[6] * eg)
    return t
end

# UV2 (uv2.f): lower-module cross-integral kernel.
@inline function _dftm_uv2(g1::Vector{Float32}, kp::Int, x1::Vector{Float32},
                           b1::Vector{Float32}, alpha::Float32, eta::Float32)::Float32
    g1[12] <= 0.000001f0 && return 0.0f0
    (kp == 0 || kp > 6) && return 0.0f0
    u = x1[3] * x1[4]
    i = kp + 12
    (g1[7] == 0.0f0 && eta == 0.0f0) && return 0.0f0
    eg = exp(g1[12])
    t = u
    t *= (b1[19] + b1[20] / g1[12])
    t *= (eg - 1.0f0)
    t *= 0.001f0
    t /= b1[21]
    t *= b1[i]
    t *= g1[6]^alpha
    t *= exp(alpha * g1[12])
    t *= (1.0f0 - (g1[6] * g1[7] * eg)^eta)
    t /= (1.0f0 - g1[6] * g1[7] * eg)
    return t
end

# GFCOMP (gfcomp.f): the S(1) lower-module G/F kernel, one KP occasion.
function _dftm_gfcomp!(g::DftmInteg)
    kp = g.kp
    x1 = g.x1; z1 = g.z1; b1 = g.b1; r1 = g.r1; g1 = g.g1; f1 = g.f1
    inrange = (kp > 0 && kp <= 6)

    g1[6] = inrange ? (1.0f0 - r1[kp]) * (1.0f0 - r1[kp+6]) * (1.0f0 - b1[kp]) : 0.0f0
    g1[7] = inrange ? (1.0f0 - r1[kp+12]) : 0.0f0
    g1[12] = 0.0f0
    if kp >= 1 && kp < 4
        g1[12] = b1[22]
    elseif kp == 4
        g1[12] = b1[23]
    elseif kp > 4 && kp <= 6
        g1[12] = b1[24]
    end
    g1[14] = max(2.0f0 * x1[4] / b1[26] - 1.0f0, 0.0f0)   # AMAX1

    # G1(15): the alpha search (label 10/15/20)
    if !(kp == 0 || kp > 6) && !(x1[1] <= 0.0f0)
        for ii in 1:11
            alpha = Float32(11 - ii)
            g1[15] = alpha
            _dftm_uv1(kp, x1, b1, alpha, g1) < x1[1] && break   # GOTO 20
        end
    else
        g1[15] = 0.0f0
    end

    # G1(16) (label 20/30)
    g1[16] = 0.0f0
    if !(kp == 0 || kp > 6)
        alpha = min(g1[15] + 1.0f0, 10.0f0)
        x1[1] > 0.0f0 && (g1[16] = _dftm_uv1(kp, x1, b1, alpha, g1))
    end

    # G1(17) (label 30/910/915/920/40)
    if (x1[1] + x1[2]) <= 0.0f0 || (kp == 0 || kp > 6)
        g1[17] = 0.0f0
    else
        g1[17] = 10.0f0 - g1[15]
        alpha = g1[15]
        nn = trunc(Int, 10.0f0 - g1[15]) + 1          # N = IFIX(10-G1(15)); NN=N+1
        for ii in 1:nn
            eta = Float32(ii - 1)
            etapls = eta + 1.0f0
            if (_dftm_uv2(g1, kp, x1, b1, alpha, etapls) -
                (x1[1] - _dftm_uv1(kp, x1, b1, alpha, g1))) > x1[2]
                g1[17] = eta                            # GOTO 915
                break
            end
        end
    end

    # G1(18) (label 40/50)
    g1[18] = 0.0f0
    if !(kp == 0 || kp > 6)
        alpha = g1[15]
        eta = 10.0f0 - alpha
        g1[18] = _dftm_uv2(g1, kp, x1, b1, alpha, eta) -
                 (x1[1] - _dftm_uv1(kp, x1, b1, alpha, g1))
        g1[18] < 0.0f0 && (g1[18] = 0.0f0)
    end

    # G1(19) via the STATIC G19TOT (reset at KP==0)
    kp == 0 && (g.g19tot = 0.0f0)
    inrange && (g.g19tot += (10.0f0 - (g1[15] + g1[17])))
    g1[19] = g.g19tot

    g1[1] = 0.0f0
    if kp == 8
        g1[1] = (z1[2] * z1[1] / 100.0f0) *
                ((x1[1] + x1[2] + z1[2]) / (2.0f0 * z1[2]) - g1[19] * b1[27])
    end
    g1[3] = kp == 8 ? x1[1] : 0.0f0
    g1[4] = inrange ? min(g1[16], x1[1]) : 0.0f0
    g1[5] = (kp >= 1 && kp <= 6) ? min(g1[18], x1[2]) : 0.0f0

    g1[8] = 0.0f0
    if !(kp == 0 || kp > 6)
        g1[8] = g1[6]^10.0f0 - 1.0f0
        if g1[15] != 10.0f0 && g1[7] != 0.0f0
            g1[8] = (g1[6]^10.0f0) * (g1[7]^(10.0f0 - g1[15])) - 1.0f0
        end
        if g1[7] == 0.0f0 && g1[15] != 10.0f0
            g1[8] = -1.0f0
        end
    end

    # G1(9) via the STATIC G17TOT (reset at KP==0 or KP==9)
    (kp == 0 || kp == 9) && (g.g17tot = 0.0f0)
    g.g17tot += g1[17]
    g1[9] = g.g17tot

    g1[10] = 0.0f0
    kp == 0 && (g1[10] = -b1[36])
    kp == 5 && (g1[10] = -(1.0f0 + g1[8]) * b1[30])
    kp == 7 && (g1[10] = -b1[34])
    kp == 8 && (g1[10] = (b1[35] * g1[14] * (1.0f0 - b1[28] * g1[9]) - 1.0f0))

    g1[11] = (g1[8] + g1[10]) * x1[3]

    g1[13] = kp == 0 ? b1[25] : 0.0f0
    if (g1[4] + g1[18]) > 0.0f0
        g.phi = (g1[4] + g1[5]) / (g1[4] + g1[18])           # STATIC PHI (else leaks)
    end
    (x1[3] <= 0.0f0 || (x1[1] + x1[2]) <= 0.0f0) && (g1[12] = 0.0f0)
    inrange && (g1[13] = x1[4] * (exp(10.0f0 * g1[12] * g.phi) - 1.0f0))

    f1[1] = g1[1] - g1[4] - g1[3]
    f1[2] = g1[3] - g1[5]
    f1[3] = g1[11]
    f1[4] = g1[13]
    return nothing
end

# Y1COMP (y1comp.f): write the S(1) states back to the class arrays.
@inline function _dftm_y1comp!(g::DftmInteg)
    c = g.icount
    g.y1[1] = g.x1[1]; g.y1[2] = g.x1[2]; g.y1[3] = g.x1[3]
    g.x5[c] = g.y1[1]; g.x6[c] = g.y1[2]; g.x7[c] = g.y1[3]
    g.g19[c] = g.g1[19]
    return nothing
end

# Y0COMP (y0comp.f): final S(0) redistribution (output X0(4) unused downstream).
@inline function _dftm_y0comp!(g::DftmInteg)
    inum = g.ic[2] - g.ic[1] + 1
    _dftm_redist!(g, 4, inum)
    return nothing
end

# DFOLE8 (dfole8.f): fill DPCENT(:,II) = 100·(1 − X6/Z3) at end of phase II/III.
@inline function _dftm_dfole8!(g::DftmInteg)
    ii = g.kk[3] - 1
    @inbounds for i in 1:g.icount
        g.g19out[i, ii] = g.g19[i]
        g.dpcent[i, ii] = 100.0f0 * (1.0f0 - g.x6[i] / g.z3[i])
    end
    return nothing
end

"""
    dftmod!(g::DftmInteg) -> g

FVS `DFTMOD(1, IFIN)`: run the full three-phase, ten-occasion integration over all
`g.n` tree classes, filling `g.dpcent`. Mirrors dftmod.f exactly (G0COMP per phase,
the Z1COMP→10×GFCOMP→Y1COMP class loop with the DO-6 NaN/negative reset, DFOLE8 at
the end of phases 2 & 3, and Y0COMP). `g.b0[57]` is set to ICOUNT as the Fortran does.
"""
function dftmod!(g::DftmInteg)
    g.kk[1] = 1; g.kk[2] = 3
    g.ic[1] = 1; g.ic[2] = g.n
    g.icount = g.ic[2] - g.ic[1] + 1
    g.b0[57] = Float32(g.icount)          # B0(57)=ICOUNT (redistribution-suppression slot)
    @inbounds for i in 1:g.n
        g.dpcent[i, 1] = 0.0f0; g.dpcent[i, 2] = 0.0f0
    end
    k1 = g.kk[1]; k2 = g.kk[2]; ic1 = g.ic[1]; ic2 = g.ic[2]
    for k3 in k1:k2
        g.kk[3] = k3
        _dftm_g0comp!(g)
        for ic3 in ic1:ic2
            g.ic[3] = ic3
            g.icount = ic3 - ic1 + 1
            _dftm_z1comp!(g)
            for i in 1:10
                g.kp = i - 1
                _dftm_gfcomp!(g)
                @inbounds for ix in 1:4
                    g.x1[ix] = g.x1[ix] + g.f1[ix]
                    g.x1[ix] >= 0.0f0 || (g.x1[ix] = 0.0f0)   # NaN or negative → RESET TO 0.0
                end
            end
            _dftm_y1comp!(g)
        end
        (k3 == 2 || k3 == 3) && _dftm_dfole8!(g)
    end
    _dftm_y0comp!(g)
    return g
end

# =============================================================================
# GARBEL / GRCLAS / GRPSUM (garbel.f, grclas.f, grpsum.f) — DFTM classification.
# =============================================================================
# TMCOUP reduces the host tree records of one species to ≤ NACLAS DFTM tree
# classes. GARBEL classifies on a weighted, standardised score of %-new-foliage
# (PCNEWF→Z2) and foliage biomass (FBIOMS→Z3), building the class sector pointers
# ISC into the (RDPSRT-sorted) pointer array and the per-class Z4=ΣPROB, Z2, Z3.
# Method 1 (NCL1 = round(NCLAS·PN1) classes) picks class bounds at the NCL1−1
# largest score gaps; method 2 (NCL2 = NCLAS−NCL1, used only when there are >5
# more records than classes) splits the largest classes. When NRECS == NCLAS the
# result is per-record identity EXCEPT that the gap-bubble + IQRSRT can leave one
# EMPTY sector (ISC low>high) — the PROB-0, FBIOMS-0 class that seeds DFTMOD's NaN.
# =============================================================================

# IQRSRT (base/iqrsrt.f): ascending integer quicksort, in place over list[1..n].
function _dftm_iqrsrt!(list::AbstractVector{<:Integer}, n::Int)
    n < 2 && return list
    iu = zeros(Int, 33); il = zeros(Int, 33)
    m = 1; i = 1; j = n
    t = 0; tt = 0; k = 0; ij = 0; l = 0
    @inbounds while true
        @label l5
        if i >= j; @goto l70; end
        @label l10
        k = i; ij = (i + j) ÷ 2; t = Int(list[ij])
        if Int(list[i]) <= t; @goto l20; end
        list[ij] = list[i]; list[i] = t; t = Int(list[ij])
        @label l20
        l = j
        if Int(list[j]) >= t; @goto l40; end
        list[ij] = list[j]; list[j] = t; t = Int(list[ij])
        if Int(list[i]) <= t; @goto l40; end
        list[ij] = list[i]; list[i] = t; t = Int(list[ij])
        @goto l40
        @label l30
        list[l] = list[k]; list[k] = tt
        @label l40
        l -= 1
        if Int(list[l]) > t; @goto l40; end
        tt = Int(list[l])
        @label l50
        k += 1
        if Int(list[k]) < t; @goto l50; end
        if k <= l; @goto l30; end
        if l - i <= j - k; @goto l60; end
        il[m] = i; iu[m] = l; i = k; m += 1; @goto l80
        @label l60
        il[m] = k; iu[m] = j; j = l; m += 1; @goto l80
        @label l70
        m -= 1
        m <= 0 && return list
        i = il[m]; j = iu[m]
        @label l80
        if j - i >= 11; @goto l10; end
        if i == 1; @goto l5; end
        i -= 1
        @label l90
        i += 1
        if i == j; @goto l70; end
        t = Int(list[i+1])
        if Int(list[i]) <= t; @goto l90; end
        k = i
        @label l100
        list[k+1] = list[k]; k -= 1
        if t < Int(list[k]); @goto l100; end
        list[k+1] = t
        @goto l90
    end
    return list
end

# GRPSUM (grpsum.f): accumulate the standardised, weighted attribute score into p.
@inline function _dftm_grpsum!(p::Vector{Float32}, lipt::AbstractVector{<:Integer},
                               nrecs::Int, atr::Vector{Float32}, wt::Float32)
    ave = 0.0f0; stdv = 0.0f0; xn = Float32(nrecs)
    @inbounds for ii in 1:nrecs
        i = Int(lipt[ii]); ave += atr[i]; stdv += atr[i] * atr[i]   # ATR(I)**2
    end
    stdv = (stdv - ave * ave / xn) / xn
    stdv = stdv > 0.000000001f0 ? sqrt(stdv) : 1.0f0
    ave = ave / xn
    @inbounds for ii in 1:nrecs
        i = Int(lipt[ii]); p[i] = p[i] + wt * (atr[i] - ave) / stdv
    end
    return nothing
end

"""
    dftm_garbel(lipt, prob, pcnewf, fbioms; w1, w2, nclas, pn1)
        -> (z4, z2, z3, isc1, isc2, lipt_sorted, kode)

FVS `GARBEL` specialised to the DFTM call (KEY=(3,4,0…,IMP=1), attributes PCNEWF
and FBIOMS, both weighted). `lipt` is the species block's record-index pointer
(sorted in place by the weighted score); `prob`/`pcnewf`/`fbioms` are indexed by
record number. Returns the per-class Z4=ΣPROB, Z2 (weighted PCNEWF), Z3 (weighted
FBIOMS), the LOCAL sector pointers `isc1`/`isc2` (1-based into the sorted `lipt`),
the sorted pointer, and `kode` (1 ⇒ ΣPROB < 0.001, classification skipped).
"""
function dftm_garbel(gipt::Vector{Int32}, base::Int, nrecs::Int, prob::Vector{Float32},
                     pcnewf::Vector{Float32}, fbioms::Vector{Float32};
                     w1::Float32 = 1.0f0, w2::Float32 = 1.0f0,
                     nclas::Int = 20, pn1::Float32 = 0.5f0)
    # The species block occupies gipt[base : base+nrecs-1]; `bm1 = base-1` maps a
    # LOCAL sector position jj (1-based) to the global index bm1+jj. GARBEL's sector
    # pointers CAN underflow (jj ≤ 0) when a class is empty — reaching into the prior
    # species block's sorted tail (a real Fortran array-adjacency artifact). Running
    # over the GLOBAL gipt (previous blocks already sorted) reproduces it bit-exactly.
    bm1 = base - 1
    lipt = view(gipt, base:base+nrecs-1)
    maxrec = length(prob)
    p = zeros(Float32, maxrec)
    # STEP1: sum PROB
    cntr = 0.0f0
    @inbounds for ii in 1:nrecs; cntr += prob[Int(lipt[ii])]; end
    if cntr < 0.001f0
        return (zeros(Float32, nclas), zeros(Float32, nclas), zeros(Float32, nclas),
                zeros(Int, nclas), zeros(Int, nclas), 1)
    end
    # STEP2: class counts by method
    ncl1 = trunc(Int, nclas * pn1 + 0.5f0)
    ncl1 < 1 && (ncl1 = 1)
    ncl2 = nclas - ncl1
    if !(nrecs - nclas > 5)
        ncl1 = nclas; ncl2 = 0
    end
    # STEP3: weighted scores (KEY(1)=PCNEWF, KEY(2)=FBIOMS)
    _dftm_grpsum!(p, lipt, nrecs, pcnewf, w1)
    _dftm_grpsum!(p, lipt, nrecs, fbioms, w2)
    # STEP4: sort the block on p (descending, RDPSRT .FALSE. = incoming order). `_rdpsrt!`
    # keys over 1..length(key), so sort a compact permutation of the block positions by
    # the record scores, then reorder the gipt slice in place — identical to RDPSRT on IPT.
    pk = Float32[p[Int(lipt[j])] for j in 1:nrecs]
    perm = Int32.(collect(1:nrecs))
    _rdpsrt!(pk, perm; lseq = false)
    blk = Int32[lipt[Int(perm[j])] for j in 1:nrecs]
    @inbounds for j in 1:nrecs; gipt[bm1+j] = blk[j]; end
    # STEP5: class bounds at the NCL1−1 largest score gaps
    diff = zeros(Float32, nclas)
    isc1 = zeros(Int, nclas); isc2 = zeros(Int, nclas)
    nc1 = ncl1 - 1
    ipt1 = Int(gipt[base])
    @inbounds for i in 2:nrecs
        ipt2 = Int(gipt[bm1+i])
        diffp = p[ipt1] - p[ipt2]
        jfin = ncl1
        for jj in 1:nc1
            if diffp <= diff[jj+1]
                jfin = jj; @goto found
            end
            diff[jj] = diff[jj+1]; isc1[jj] = isc1[jj+1]
        end
        @label found
        diff[jfin] = diffp; isc1[jfin] = i
        ipt1 = ipt2
    end
    isc1[1] = nrecs + 1
    _dftm_iqrsrt!(isc1, ncl1)
    # STEP6a: class lengths; redefine isc1 to point to class END
    i1 = 1
    @inbounds for jj in 1:ncl1
        i2 = isc1[jj] - 1
        isc1[jj] = i2
        isc2[jj] = i2 - i1 + 1
        i1 = i2 + 1
    end
    # STEP6b: method 2 — split the NCL2 largest classes
    if ncl2 != 0
        @inbounds for _k in 1:ncl2
            mx = 0; isel = 0
            for jj in 1:ncl1
                if isc2[jj] > mx; mx = isc2[jj]; isel = jj; end
            end
            ncl1 += 1
            mx = trunc(Int, Float32(mx) / 2.0f0 + 0.5f0)
            isc1[ncl1] = isc1[isel]
            isc1[isel] = isc1[isel] - mx
            isc2[ncl1] = mx
            isc2[isel] = isc2[isel] - mx
        end
        _dftm_iqrsrt!(isc1, ncl1)
    end
    # STEP6c/147: assign sector pointers (start,end)
    i1 = 1
    @inbounds for jj in 1:ncl1
        i2 = isc1[jj]
        isc1[jj] = i1
        isc2[jj] = i2
        i1 = i2 + 1
    end
    # STEP7: per-class Z4 = ΣPROB over the (possibly under/overflowing) sector
    z4 = zeros(Float32, nclas); z2 = zeros(Float32, nclas); z3 = zeros(Float32, nclas)
    for i in 1:nclas
        i1 = isc1[i]; i2 = isc2[i]
        acc = 0.0f0
        for jj in i1:i2
            acc += prob[Int(gipt[bm1+jj])]
        end
        z4[i] = acc
    end
    # GRCLAS: weighted per-class attribute means Z2 (PCNEWF), Z3 (FBIOMS)
    _dftm_grclas!(z2, gipt, bm1, isc1, isc2, nclas, pcnewf, prob, z4)
    _dftm_grclas!(z3, gipt, bm1, isc1, isc2, nclas, fbioms, prob, z4)
    return (z4, z2, z3, isc1, isc2, 0)
end

# GRCLAS (grclas.f): per-class Σ(atr·PROB)/Σ(PROB) into out. `bm1 = base-1` maps a
# LOCAL sector position to the global gipt index (may underflow into the prior block).
@inline function _dftm_grclas!(out::Vector{Float32}, gipt::Vector{Int32}, bm1::Int,
                               isc1::Vector{Int}, isc2::Vector{Int}, nclas::Int,
                               atr::Vector{Float32}, prob::Vector{Float32},
                               z4::Vector{Float32})
    for j in 1:nclas
        i1 = isc1[j]; i2 = isc2[j]
        acc = 0.0f0
        for ii in i1:i2
            i = Int(gipt[bm1+ii]); acc += atr[i] * prob[i]
        end
        xp = z4[j]
        xp > 1.0f-30 && (acc = acc / xp)
        out[j] = acc
    end
    return nothing
end

# =============================================================================
# TMCOUP damage (tmcoup.f) — % branch defoliation → tree defoliation class →
# diameter-growth loss, height-growth loss, top-kill, and mortality.
# =============================================================================
# AMORT(14,9): rows = the 14 tree-defoliation classes (1–7 DF, 8–14 GF), columns
# 1=primary-mortality %, 2=secondary-mortality %, 3=top-kill probability test,
# 5–9=cumulative top-kill-class probabilities. Transcribed COLUMN-MAJOR from the
# tmcoup.f DATA statement. RGLOSS(14) = the static per-class diameter-growth-loss
# fractions. TKBOT/TKTOP = the top-kill-class crown thresholds.
const DFTM_AMORT = permutedims(reshape(Float32[
    # col1 primary mortality
    0.0,0.0,0.0090,0.0280,0.1730,0.4770,0.9230, 0.0,0.0,0.0090,0.0280,0.1730,0.4770,0.9230,
    # col2 secondary mortality
    0.0760,0.0760,0.0760,0.1075,0.1999,0.1144,0.0168, 0.0442,0.0428,0.0452,0.0587,0.0926,0.0778,0.0114,
    # col3 top-kill probability test
    0.1158,0.2081,0.4163,0.4820,0.2888,0.2005,0.0295, 0.0559,0.1271,0.2703,0.3838,0.3441,0.2149,0.0316,
    # col4 (unused)
    0.7570,0.6562,0.4736,0.3825,0.3383,0.2081,0.0306, 0.9000,0.8301,0.6755,0.5295,0.3903,0.2303,0.0339,
    # col5..9 cumulative top-kill class probabilities
    0.0694,0.1360,0.1800,0.2376,0.1376,0.1532,0.0226, 0.0341,0.0877,0.1332,0.1733,0.1933,0.0963,0.0142,
    0.0464,0.0696,0.2125,0.1606,0.0599,0.0384,0.0057, 0.0123,0.0293,0.0891,0.1021,0.0587,0.0401,0.0059,
    0.0,0.0,0.0,0.0,0.0,0.0,0.0, 0.0039,0.0078,0.0317,0.0434,0.0225,0.0274,0.0040,
    0.0,0.0,0.0238,0.0838,0.0775,0.0,0.0, 0.0038,0.0019,0.0122,0.0479,0.0157,0.0,0.0,
    0.0,0.0026,0.0,0.0,0.0138,0.0089,0.0013, 0.0018,0.0005,0.0040,0.0172,0.0539,0.0511,0.0075,
], 14, 9), (1, 2))   # → DFTM_AMORT[itab, col]
const DFTM_RGLOSS = Float32[0.846,0.716,0.648,0.610,0.569,0.388,0.388,
                            0.785,0.693,0.628,0.596,0.559,0.455,0.455]
const DFTM_TKBOT = Float32[0.0, 0.0, 0.10, 0.25, 0.50]
const DFTM_TKTOP = Float32[0.0, 0.10, 0.25, 0.50, 0.99]
# The crown-defoliation equation parameters (tmcoup.f P1/P2/P3).
const DFTM_P1 = -22.579678f0
const DFTM_P2 = -22.996253f0
const DFTM_P3 = -598181961.0f0

"""
    dftm_tree_defol(dpmax, g19max, iz6) -> (t, itab)

FVS `TMCOUP` percent-total-tree-defoliation `T` and AMORT-row index `ITAB` for a
tree class. `T` from the maximum food-shortage days `g19max` (≥5 ⇒ 100%; >0.001 ⇒
linear) else the crown-defoliation logistic in the max branch defoliation `dpmax`.
`ITAB` bins T (1–7), +7 for grand fir (`iz6==2`). A NaN `dpmax` still yields a
finite ITAB via the G19 branch (the empty class ⇒ G19MAX>0 ⇒ T=100 ⇒ ITAB 7/14).
"""
@inline function dftm_tree_defol(dpmax::Float32, g19max::Float32, iz6::Integer)
    local t::Float32
    if g19max > 0.001f0
        t = 97.7911f0 + (g19max / 5.0f0) * 2.2089f0
        g19max >= 5.0f0 && (t = 100.0f0)
    else
        t = (DFTM_P1 / (DFTM_P2 + DFTM_P3 * exp(DFTM_P1 * dpmax / 100.0f0))) * 100.0f0
    end
    itab = 1
    (t > 15.0f0 && t <= 35.0f0) && (itab = 2)
    (t > 35.0f0 && t <= 65.0f0) && (itab = 3)
    (t > 65.0f0 && t <= 85.0f0) && (itab = 4)
    (t > 85.0f0 && t <= 95.0f0) && (itab = 5)
    (t > 95.0f0 && t <= 99.5f0) && (itab = 6)
    t > 99.5f0 && (itab = 7)
    iz6 == 2 && (itab += 7)
    return (t, itab)
end

"""
    dftm_mortality(itab, dpmax, prob, wk2_bg) -> wk2

FVS `TMCOUP` DO-320 mortality: `PRMORT = AMORT(ITAB,1)+AMORT(ITAB,2)`, floored at the
normal rate `PRNORM = WK2/PROB`; for the lightest classes (ITAB 1 or 8) interpolated
back toward PRNORM by `(1 − DPMAX/68.03)`. Returns the DFTM mortality `WK2 = PROB·PRMORT`.
Deterministic (no RNG). `wk2_bg` is the pre-DFTM background periodic mortality.
"""
@inline function dftm_mortality(itab::Int, dpmax::Float32, prob::Float32, wk2_bg::Float32)::Float32
    prmort = DFTM_AMORT[itab, 1] + DFTM_AMORT[itab, 2]
    prnorm = wk2_bg / prob
    prmort < prnorm && (prmort = prnorm)
    if itab == 1 || itab == 8
        delta = (prmort - prnorm) * (1.0f0 - (dpmax / 68.03f0))
        delta >= 0.0f0 && (prmort = prmort - delta)
    end
    return prob * prmort
end

"""
    dftm_dgloss(itab, dg) -> (dg_new, dgloss)

FVS `TMCOUP` diameter-growth loss: `DGLOSS = DG·(1 − RGLOSS(ITAB))` capped at `DG`;
`DG_new = DG − DGLOSS`. Deterministic (no RNG).
"""
@inline function dftm_dgloss(itab::Int, dg::Float32)::NTuple{2,Float32}
    dgloss = dg * (1.0f0 - DFTM_RGLOSS[itab])
    dgloss > dg && (dgloss = dg)
    return (dg - dgloss, dgloss)
end

"""
    dftm_htgloss_notopkill(dpmax, htg, fint) -> (htg_new, htgloss)

FVS `TMCOUP` height-growth loss on the NO-top-kill branch: `HTGLOS = HTG·DPMAX/100·(1/FINT)`;
`HTG_new = HTG − HTGLOS` (then floored at 0 by the caller). Deterministic. The top-kill
branch (`PRTOPK < AMORT(ITAB,3)`) is stochastic — one `TMRANN` draw per tree — and is NOT
covered here; the caller applies it when the RNG stream is threaded through.
"""
@inline function dftm_htgloss_notopkill(dpmax::Float32, htg::Float32, fint::Float32)::NTuple{2,Float32}
    htgloss = htg * dpmax / 100.0f0 * (1.0f0 / fint)
    return (htg - htgloss, htgloss)
end

# -----------------------------------------------------------------------------
# TMBCHL (tmbchl.f) — Batchelor composite-rejection normal draw over TMRANN.
# -----------------------------------------------------------------------------
"""
    dftm_tmbchl!(d, xbar, stdev) -> Float32

FVS `TMBCHL(XBAR,STDEV)`: a normal random variate from N(XBAR,STDEV) by the
Batchelor composite-rejection method, consuming three `TMRANN` uniforms per
attempt (`U`, `R1`, `R2`) and redrawing all three on rejection (`Y = −ln(R1) ≤ Z`).
`U ≤ 2/3` uses the uniform branch (`X = 1.5U`, no transcendental in the returned
value); `U > 2/3` the negative-exponential branch. The sign is `sign(0.5 − R2)`.
Float32 throughout to mirror the REAL Fortran; drives the RANLARVA egg allocation
(and, when reached, the between-outbreak larval-density draws) off the DFTM stream.
"""
function dftm_tmbchl!(d::DftmState, xbar::Float32, stdev::Float32)::Float32
    @inbounds while true
        u  = dftm_rand!(d)
        r1 = dftm_rand!(d)
        r2 = dftm_rand!(d)
        local x::Float32, z::Float32
        if u <= (2.0f0 / 3.0f0)
            x = 1.5f0 * u
            z = 0.5f0 * x * x
        else
            x = 1.0f0 - 0.5f0 * log(3.0f0 * u - 2.0f0)   # ALOG
            z = 0.5f0 * (x - 2.0f0)^2
        end
        y = -log(r1)                                     # −ALOG(R1)
        y <= z && continue                               # reject → GOTO 10 (redraw all three)
        a = 0.5f0 - r2
        signx = a / abs(a)                               # +1 if R2<0.5, −1 if R2>0.5
        x = x * signx
        return x * stdev + xbar
    end
end

"""
    dftm_alloc_eggs!(d, iz6, jclas2, dfmean, dfsd, gfmean, gfsd) -> x7

FVS `TMCOUP` DO-170 RANLARVA egg allocation (`IEGTYP==1`): in the `JCLAS2`
(descending-DBH) class order, draw `TMBCHL(DFMEAN,DFREGG(2))` for Douglas-fir
classes (`IZ6==1`) and `TMBCHL(GFMEAN,GFREGG(2))` for grand-fir classes,
flooring negatives at 0, and store the result at the class index `I=JCLAS2(II)`.
Returns the per-class egg vector `X7`. Threads the DFTM `TMRANN` stream through `d`.
"""
function dftm_alloc_eggs!(d::DftmState, iz6::AbstractVector{<:Integer},
                          jclas2::AbstractVector{<:Integer}, dfmean::Float32, dfsd::Float32,
                          gfmean::Float32, gfsd::Float32)
    n = length(iz6)
    x7 = zeros(Float32, n)
    @inbounds for ii in 1:n
        i = Int(jclas2[ii])
        if iz6[i] == 1
            x7[i] = dftm_tmbchl!(d, dfmean, dfsd)
        else
            x7[i] = dftm_tmbchl!(d, gfmean, gfsd)
        end
        x7[i] < 0.0f0 && (x7[i] = 0.0f0)
    end
    return x7
end

"""
    dftm_topkill!(d, itab, dpmax, ht, htg, icr, normht, itrunc, fint)
        -> (ht_new, htg_new, icr_new, normht_new, itrunc_new, jtrunk, pckill, tkill, htgloss)

FVS `TMCOUP` DO-380 per-tree height/top-kill damage. Draws `PRTOPK` (one TMRANN);
if `PRTOPK < AMORT(ITAB,3)` and `DPMAX > 0.0001`, a top-kill occurs: the cumulative
`AMORT(ITAB,5..9)` picks the top-kill class `K` — `K==1` kills one year of leader
growth (`HTGLOS = HTG/FINT`); `K>1` draws a second uniform `RANDOM`, computes
`PCKILL = TKBOT(K)+RANDOM·(TKTOP(K)−TKBOT(K))`, truncates the crown (setting
`NORMHT`/`ITRUNC` when `JTRUNK==1`), and lowers `HT`/`HTG`. Otherwise the no-top-kill
reduction `HTGLOS = HTG·DPMAX/100/FINT` applies. Threads the TMRANN stream through `d`.
`icr`/`normht`/`itrunc` are the tree's pre-damage crown-ratio %/NORMHT/ITRUNC. Returns
the updated geometry (HTG floored at 0 by the caller as in the Fortran).
"""
function dftm_topkill!(d::DftmState, itab::Int, dpmax::Float32, ht::Float32, htg::Float32,
                       icr::Int, normht::Int, itrunc::Int, fint::Float32)
    icri = icr < 0 ? -icr : icr
    cbase = ht * Float32(100 - icri) / 100.0f0
    crown = ht - cbase
    jtrunk = 0; pckill = 0.0f0; tkill = 0.0f0; htgloss = 0.0f0
    prtopk = dftm_rand!(d)
    if prtopk < DFTM_AMORT[itab, 3] && dpmax > 0.0001f0
        prsum = 0.0f0
        for k in 1:5
            prsum += DFTM_AMORT[itab, 4 + k]
            if prtopk <= prsum
                if k == 1
                    htgloss = htg * (1.0f0 / fint)
                    htg = htg - htgloss
                else
                    random = dftm_rand!(d)
                    pckill = DFTM_TKBOT[k] + random * (DFTM_TKTOP[k] - DFTM_TKBOT[k])
                    htgloss = crown * pckill
                    tkill = htgloss
                    htrunc = ht - htgloss
                    jtrunk = 0
                    (htgloss > 5.0f0 && ht > 20.0f0) && (jtrunk = 1)
                    (Float32(itrunc) / 100.0f0 < htrunc && itrunc > 0) && (jtrunk = 0)
                    if jtrunk == 1
                        normht = trunc(Int, ht * 100.0f0 + 0.5f0)     # IFIX(HT*100+0.5)
                        itrunc = trunc(Int, htrunc * 100.0f0 + 0.5f0)
                    end
                    ht = htrunc
                    htg = htg - htgloss
                end
                break
            end
        end
    else
        htgloss = htg * dpmax / 100.0f0 * (1.0f0 / fint)
        htg = htg - htgloss
    end
    return (ht, htg, icri, normht, itrunc, jtrunk, pckill, tkill, htgloss)
end
