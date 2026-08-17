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
