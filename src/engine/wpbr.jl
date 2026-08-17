# =============================================================================
# White Pine Blister Rust (WPBR) model — wpbr/*.f
# =============================================================================
# Third port of the FVS insect/pathogen event-extension family (after DFB and
# DFTM). WPBR is a *canker* pathology model for 5-needle (white) pines: a stand
# rust index (RI) + per-tree growth index (GI) drive expected-canker generation
# (BRECAN), canker growth (BRCGRO), and canker-status → top-kill + mortality
# (BRTSTA/BRCGRO), feeding growth-loss + mortality back into the FVS treelist.
#
# FVS structure (wpbr/, 43 files):
#   * BRIN    (brin.f)   — the WPBR keyword-block reader (keywds.f option 75,
#                          keyword name "BRUST"; cf. DFB=100, DFTM=7, RD=101).
#                          18 sub-keywords → the BRCOM common. Ported `kw_brin!`.
#   * BRANN   (brann.f)  — the model's OWN double-precision Lehmer/MINSTD LCG:
#                          BRS1 = DMOD(16807·BRS0, 2147483647); SEL = BRS1/2^31.
#                          Default seed BRSS = 55329 (odd). ENTRY BRNSED reseeds
#                          (forces odd). IDENTICAL algorithm to DFB's DFBRAN and
#                          DFTM's TMRANN. Ported `wpbr_rand!`/`wpbr_seed!`. NEVER
#                          FFI'd.
#   * BRBLKD  (brblkd.f + brblkd{cr,ie,so}.f) — block data: RNG seed, host-species
#                          map BRSPM (variant-specific), alpha codes BRSPC, damage
#                          code IBRDAM=36, BRPI, RSF, canker-file format ICFMT.
#   * BRINIT  (brinit.f) — per-stand init of the keyword-settable defaults
#                          (pruning/excising specs, DFACT, RATINV, RIDEF, RIBUS,
#                          canker growth rates BRGRTH/BOGRTH, RESIST/PRPSTK) and a
#                          BRNSED(.FALSE.) RNG reset. Ported `wpbr_defaults!`.
#   * BRSETP  (brsetp.f) — per-tree init (ground diameter BRGD, height-to-crown);
#                          called from MAIN.  [later chunk]
#   * BRTREG  (brtreg.f) — the per-cycle driver (gradd.f:124, "begin simulation"),
#                          analogous to DFTMGO/DFBGO.  [later chunk]
#   * BRECAN/BRCANK/BRCINI/BRCGRO/BRCDEL/BRCREM/BRCSTA/BRTSTA/BRGI/BRTARG/BRIBA/
#     BRIBES/BRICAL — the canker/infection dynamics + growth/mortality kernels.
#     Randomness (CALL BRANN) enters in BRECAN (expected cankers), BRCINI (random
#     canker up/out/girdle), BRCGRO (mortality realization), BRCREM (prune/excise
#     success), BRSTYP (stock-type assignment).  [later chunks]
#
# LINKAGE: like DFB/DFTM, NO shipped binary runs WPBR — every FVS*_buildDir links
# the base/exbrus.f NO-OP stub (its entries BRIN/BRSETP/BRTREG/BRINIT/BRDAM/BRCMPR/
# BRESTB/BRTDEL/BRSOR/BRTRIP/BRPR/BRROUT/BRKEY all RETURN). A WPBR oracle must be
# relinked = a variant's .o set with the real wpbr/*.o swapped in for exbrus.o
# (see scratchpad/wpbr/build_ie_wpbr.sh).
#
# THIS PORT (chunk 0 — keyword reader + RNG + defaults, INERT seam):
#   * `kw_brin!`  — faithful BRIN block reader → `s.wpbr` (WpbrState).
#   * `wpbr_rand!` / `wpbr_seed!` — the BRANN LCG + BRNSED reseed (bit-exact).
#   * `wpbr_defaults!` — BRINIT/BRBLKD defaults; `wpbr_brspm` variant host map.
#   NO per-cycle engine seam is wired: a stand carrying a WpbrState projects
#   BYTE-IDENTICALLY to one without WPBR (validated by the FVSjl suite + the
#   relinked FVSie_wpbr oracle running WPBR-off ≡ stock, and a BRUST keyfile whose
#   block parses but has no engine effect). The canker dynamics (BRTREG/BRECAN/
#   BRCGRO/…) are the next, much larger, chunks.
#
# Float32 discipline: WPBR's Fortran is REAL (Float32) except the RNG state (DOUBLE
# BRS0/BRS1/BRSS). Coefficient arrays and every ported arithmetic result are
# Float32 to stay bit-identical to the oracle; the LCG state is Float64.
# =============================================================================

const WPBR_DEFAULT_SEED = 55329.0f0   # BRBLKD: BRSS default (odd); BRS0 init 55329

# -----------------------------------------------------------------------------
# BRBLKD variant host-species map (BRSPM) + alpha codes (BRSPC).
# -----------------------------------------------------------------------------
# BRSPM(MAXSP): a non-zero value is the WPBR host index (1 or 2) of that FVS
# species; 0 = not a blister-rust host. NBRSP=2 host species per variant.
#   base/NI(IE 11-sp legacy): WP@1                         BRSPC = WP,SP
#   IE (23-sp):               WP@1, LM@13                  BRSPC = WP,LM
#   CR (38-sp):               LP@11, PP@13                 BRSPC = LP,PP
#   SO (33-sp):               WP@1, SP@2                   BRSPC = WP,SP
"""
    wpbr_brspm(variant) -> (brspm::Vector{Int}, brspc::NTuple{2,String})

FVS `BRBLKD` host-species map for a WPBR-linked variant: `brspm[i]` is the WPBR
host index (1/2) of FVS species `i`, 0 if not a host. Faithful transcription of
`brblkd{,ie,cr,so}.f`. Falls back to the base/NI single-host map (WP@1).
"""
function wpbr_brspm(variant)
    if variant isa InlandEmpire
        m = zeros(Int, 23); m[1] = 1; m[13] = 2
        return (m, ("WP", "LM"))
    elseif variant isa SouthCentralOregon
        m = zeros(Int, 33); m[1] = 1; m[2] = 2
        return (m, ("WP", "SP"))
    elseif variant isa CentralRockies
        m = zeros(Int, 38); m[11] = 1; m[13] = 2
        return (m, ("LP", "PP"))
    else
        m = zeros(Int, 11); m[1] = 1
        return (m, ("WP", "SP"))
    end
end

"""
    WpbrState

White Pine Blister Rust model state — the FVS `BRCOM`/`DPCOM`/`RIBES` commons
(keyword-settable subset for chunk 0). Populated by `kw_brin!` and seeded with the
`BRINIT`/`BRBLKD` defaults by `wpbr_defaults!`. `active` mirrors "a BRUST block was
read" (BRYES); the future per-cycle engine seam is gated on it, so a stand without
a BRUST block projects byte-identically. NBRSP=2 host species.
"""
mutable struct WpbrState <: AbstractWpbrState
    active::Bool               # BRYES — a BRUST block was read (engine seam gate)
    # --- BRANN LCG state (DPCOM: DOUBLE BRS0/BRS1/BRSS) ---
    brs0::Float64              # BRS0 — current generator state (exact mod)
    brss::Float64              # BRSS — reseed default (BRNSED LSET=false resets to it)
    # --- host map (BRBLKD) ---
    brspm::Vector{Int}         # BRSPM — FVS-species → WPBR host index (0=non-host)
    brspc::NTuple{2,String}    # BRSPC — host alpha codes
    ibrdam::Int32              # IBRDAM — WPBR damage code (36)
    brpi::Float32              # BRPI  — 3.14159
    # --- pruning / excising specs (BRINIT), stored in FVS internal units (cm) ---
    srate::NTuple{2,Float32}   # SRATE — prune / excise success rate (0.9, 0.5)
    htprpr::Float32            # HTPRPR — max prunable ht as proportion of tree ht (0.50)
    exdmin::Float32            # EXDMIN — min DBH for excising, cm (3.0*2.54)
    htmax::NTuple{2,Float32}   # HTMAX — max prune / excise ht, cm (8*30.48, 6*30.48)
    girmax::Float32            # GIRMAX — max % girdle still excisable (50)
    girmrt::Float32            # GIRMRT — % girdle causing mortality/top-kill (100)
    htmin::Float32             # HTMIN — min bole-canker ht for excising, cm (3*2.54)
    outdst::Float32            # OUTDST — min dist-out for prunable canker, cm (6*2.54)
    outnld::Float32            # OUTNLD — dist-out defining non-lethal canker, cm (24*2.54)
    # --- rust / growth index (BRINIT) ---
    dfact::Matrix{Float32}     # DFACT(NBRSP,4) — stand deviation factor by sp/stock (0.33)
    ldfact::Bool               # LDFACT — DEVFACT/RUSTINDX supplied a static factor
    riaf::Vector{Float32}      # RIAF(NBRSP) — rust index adjustment factor (1.0)
    ratinv::NTuple{2,Float32}  # RATINV — canker inactivation rate branch/bole (0.05, 0.01)
    gidef::Float32             # GIDEF — default growth index (15.0)
    ridef::Float32             # RIDEF — default rust index (0.015)
    rimeth::Int32              # RIMETH — RI assignment method (0)
    ribprp::NTuple{3,Float32}  # RIBPRP — ribes species proportions (0, 0.5, 0.5)
    ribus::Matrix{Float32}     # RIBUS(2,3) — ribes bushes/acre old/new by species
    minri::Float32             # MINRI — RUSTINDX method 3/4 min RI
    maxri::Float32             # MAXRI — max RI
    pkage::Float32             # PKAGE — stand age at RI peak
    pkshp::Float32             # PKSHP — RI curve shape
    rsf::NTuple{3,Float32}     # RSF — ribes species calibration factors (2.3, 1.0, 0.64)
    # --- canker growth rates (BRINIT) ---
    brgrth::Matrix{Float32}    # BRGRTH(NBRSP,4) — branch canker radial growth cm (5.0)
    bogrth::Matrix{Float32}    # BOGRTH(NBRSP,4) — bole canker diameter growth cm (4.5)
    # --- stock (BRINIT) ---
    resist::Matrix{Float32}    # RESIST(NBRSP,4) — stock-type resistance (1, 0.33, 0.17, 0.11)
    prpstk::Matrix{Float32}    # PRPSTK(NBRSP,4) — stock-type proportions (1, 0, 0, 0)
    # --- report / IO control (BRINIT/BRBLKD) ---
    lbrsum::Bool               # LBRSUM — print stand summary stats (true)
    lbrdbh::Bool               # LBRDBH — print DBH-class stats (false)
    lmetric::Bool              # LMETRIC — canker data in cm (false)
    lpatpr::Bool               # LPATPR — pathological pruning (false)
    brtl::Bool                 # BRTL  — BRTLST tree-list output requested
    brcl::Bool                 # BRCL  — BRCLST canker-list output requested
    ckinit::Bool               # CKINIT — CANKDATA random canker generation
    icin::Int32                # ICIN  — canker-data logical unit (55)
    idtout::Int32              # IDTOUT — tree-list unit (56)
    idcout::Int32              # IDCOUT — canker-list unit (57)
    icfmt::String              # ICFMT — canker-data read format
    # --- scheduled activities (OPNEW) recorded for the future engine seam ---
    # each entry = (idt, iact, prms) — NOT acted on in the inert chunk.
    activities::Vector{Tuple{Int,Int,Vector{Float32}}}
end

# NBRSP host-species count (parameter in BRCOM.F77).
const WPBR_NBRSP = 2

"""
    wpbr_defaults!(variant) -> WpbrState

FVS `BRINIT` + `BRBLKD`: the WPBR per-stand run-time defaults, with the variant's
host-species map. Called when the first BRUST keyword is seen. The RNG state is
`BRS0 = BRSS = 55329` (BRINIT calls `BRNSED(.FALSE.)` which resets `BRS0` to the
`BRBLKD` default `BRSS`).
"""
function wpbr_defaults!(variant)
    brspm, brspc = wpbr_brspm(variant)
    n = WPBR_NBRSP
    dfact  = fill(0.33f0, n, 4)
    resist = Matrix{Float32}(undef, n, 4)
    prpstk = Matrix{Float32}(undef, n, 4)
    brgrth = fill(5.0f0, n, 4)
    bogrth = fill(4.5f0, n, 4)
    for i in 1:n
        resist[i, 1] = 1.00f0; resist[i, 2] = 0.33f0; resist[i, 3] = 0.17f0; resist[i, 4] = 0.11f0
        prpstk[i, 1] = 1.0f0;  prpstk[i, 2] = 0.0f0;  prpstk[i, 3] = 0.0f0;  prpstk[i, 4] = 0.0f0
    end
    ribus = Float32[0.0 200.0 200.0; 0.0 75.0 75.0]   # RIBUS(2,3): row1=old, row2=new
    return WpbrState(
        false,                                  # active (BRYES false until BRIN sets true)
        Float64(WPBR_DEFAULT_SEED), Float64(WPBR_DEFAULT_SEED),   # brs0, brss
        brspm, brspc, Int32(36), 3.14159f0,     # brspm, brspc, ibrdam, brpi
        (0.9f0, 0.5f0), 0.50f0,                 # srate, htprpr
        3.0f0 * 2.54f0,                         # exdmin
        (8.0f0 * 30.48f0, 6.0f0 * 30.48f0),     # htmax
        50.0f0, 100.0f0,                        # girmax, girmrt
        3.0f0 * 2.54f0, 6.0f0 * 2.54f0, 24.0f0 * 2.54f0,   # htmin, outdst, outnld
        dfact, false, fill(1.0f0, n),           # dfact, ldfact, riaf
        (0.05f0, 0.01f0), 15.0f0, 0.015f0, Int32(0),   # ratinv, gidef, ridef, rimeth
        (0.0f0, 0.5f0, 0.5f0), ribus,           # ribprp, ribus
        0.0f0, 0.0f0, 0.0f0, 0.0f0,             # minri, maxri, pkage, pkshp
        (2.3f0, 1.0f0, 0.64f0),                 # rsf
        brgrth, bogrth, resist, prpstk,         # brgrth, bogrth, resist, prpstk
        true, false, false, false,              # lbrsum, lbrdbh, lmetric, lpatpr
        false, false, false,                    # brtl, brcl, ckinit
        Int32(55), Int32(56), Int32(57),        # icin, idtout, idcout
        "(I7,1X,I1,1X,F3.0,1X,F5.1,1X,F5.1,1X,F4.0,1X,F4.0)",   # icfmt
        Tuple{Int,Int,Vector{Float32}}[],       # activities
    )
end

# -----------------------------------------------------------------------------
# BRANN (brann.f) — the WPBR model's own double-precision Lehmer/MINSTD LCG.
# -----------------------------------------------------------------------------
"""
    wpbr_rand!(w) -> Float32

FVS `BRANN`: `BRS1 = DMOD(16807·BRS0, 2147483647)`; `SEL = REAL(BRS1/2147483648)`;
`BRS0 = BRS1`. `BRS0`/`BRS1` are DOUBLE, so the modular step is exact integer
arithmetic; only the returned uniform `SEL` truncates to Float32. IDENTICAL to
DFB's DFBRAN and DFTM's TMRANN. Default seed 55329. NEVER FFI'd.
"""
@inline function wpbr_rand!(w::WpbrState)::Float32
    s1 = rem(16807.0 * w.brs0, 2147483647.0)   # DMOD (exact, values < 2^31)
    w.brs0 = s1
    return Float32(s1 / 2147483648.0)          # SEL = REAL(BRS1 / 2^31)
end

"""
    wpbr_seed!(w, seed, present)

FVS `BRNSED(LSET, SEED)` (brann.f ENTRY). BRSEED with a value (`present`): force an
ODD seed, store it as both the working state `BRS0` and the reseed default `BRSS`.
BRSEED with no value (`present=false`): reset `BRS0` to the current `BRSS`.
"""
function wpbr_seed!(w::WpbrState, seed::Float32, present::Bool)
    if !present
        w.brs0 = w.brss                          # LSET=false: SEED=BRSS; BRS0=SEED
        return nothing
    end
    (seed % 2.0f0 == 0.0f0) && (seed += 1.0f0)   # AMOD(SEED,2)==0 → SEED+1 (odd)
    w.brss = Float64(seed)
    w.brs0 = Float64(seed)
    return nothing
end

# -----------------------------------------------------------------------------
# kw_brin! (brin.f) — WPBR keyword-block reader (keywds.f option 75, "BRUST").
# -----------------------------------------------------------------------------
"""
    kw_brin!(s, rec, kr)

Parse the `BRUST … END` block (wpbr/brin.f). Faithfully sets the BRCOM-equivalent
state on `s.wpbr` and consumes sub-keyword records up to `END`, exactly as FVS
`BRIN` does. No per-cycle engine seam is wired yet, so this is INERT: a stand
carrying a WpbrState still projects byte-identically to one without WPBR.

Sub-keywords ported (option order 1..18): PRUNE, PRNSPECS, EXCISE, EXSPECS, RIBES,
INACT, END, COMMENT, GROWRATE, BRTLST, BRCLST, CANFMT, DEVFACT, RUSTINDX, BRSEED,
STOCK, CANKDATA, BROUT. Management activities (PRUNE/EXCISE/…) are recorded in
`w.activities` for the future engine seam but not acted on. Supplemental records
(COMMENT lines, CANFMT format, RUSTINDX method 3/4 params, BRTLST/BRCLST/CANKDATA
file names) are consumed to keep the reader byte-aligned. Canker-data reads
(CANKDATA → BRCANK) are deferred: the keyword is parsed and its supplemental
record consumed, but no canker file is opened in the inert chunk.
"""
function kw_brin!(s::StandState, rec, kr::KeywordReader)
    s.wpbr === nothing && (s.wpbr = wpbr_defaults!(s.variant))
    w = s.wpbr
    w.active = true                              # BRYES = .TRUE.

    while true
        r = read_keyword!(kr)
        (r.status == KW_EOF || r.status == KW_STOP) && break
        k = strip(r.name)
        isempty(k) && continue
        if k == "END"                            # option 7
            break
        elseif k == "PRUNE"                       # option 1 — schedule pruning
            idt = (r.present[1] && r.values[1] > 0.0f0) ? trunc(Int, r.values[1]) : 1
            p1 = r.present[2] ? Float32(r.values[2]) : w.srate[1]
            _wpbr_sched!(w, idt, 1001, Float32[p1, Float32(r.values[3]), Float32(r.values[4]),
                                               (r.present[5] && r.values[5] != 0.0f0) ? 1.0f0 : 0.0f0])
        elseif k == "PRNSPECS"                     # option 2 — prune thresholds
            idt = r.present[1] ? trunc(Int, r.values[1]) : 0
            if idt == 0
                r.present[2] && (w.htprpr = Float32(r.values[2]))
                r.present[3] && (w.htmax = (Float32(r.values[3]) * 30.48f0, w.htmax[2]))
                r.present[4] && (w.outdst = Float32(r.values[4]) * 2.54f0)
                r.present[5] && (w.outnld = Float32(r.values[5]) * 2.54f0)
            else
                _wpbr_sched!(w, idt, 1002, _wpbr_prms(r, 2, 5))
            end
        elseif k == "EXCISE"                       # option 3 — schedule excising
            idt = r.present[1] ? trunc(Int, r.values[1]) : 1
            p1 = r.present[2] ? Float32(r.values[2]) : w.srate[2]
            _wpbr_sched!(w, idt, 1003, Float32[p1])
        elseif k == "EXSPECS"                      # option 4 — excise thresholds
            idt = r.present[1] ? trunc(Int, r.values[1]) : 0
            if idt == 0
                r.present[2] && (w.exdmin = Float32(r.values[2]) * 2.54f0)
                r.present[3] && (w.htmax = (w.htmax[1], Float32(r.values[3]) * 30.48f0))
                r.present[4] && (w.girmax = Float32(r.values[4]))
                r.present[5] && (w.girmrt = Float32(r.values[5]))
                r.present[6] && (w.htmin = Float32(r.values[6]) * 2.54f0)
            else
                _wpbr_sched!(w, idt, 1004, _wpbr_prms(r, 2, 6))
            end
        elseif k == "RIBES"                        # option 5 — ribes populations → RI
            idt = r.present[1] ? trunc(Int, r.values[1]) : 1
            if idt == 0
                r.present[2] && (w.ribus[1, 1] = Float32(r.values[2]))
                r.present[3] && (w.ribus[2, 1] = Float32(r.values[3]))
                r.present[4] && (w.ribus[1, 2] = Float32(r.values[4]))
                r.present[5] && (w.ribus[2, 2] = Float32(r.values[5]))
                r.present[6] && (w.ribus[1, 3] = Float32(r.values[6]))
                r.present[7] && (w.ribus[2, 3] = Float32(r.values[7]))
                # BRIBES(REDFAC) recomputes RIDEF — deferred to the dynamics chunk.
            else
                _wpbr_sched!(w, idt, 1005, _wpbr_prms(r, 2, 7))
            end
        elseif k == "INACT"                        # option 6 — canker inactivation rate
            idt = r.present[1] ? trunc(Int, r.values[1]) : 0
            if idt == 0
                r1 = r.present[2] ? Float32(r.values[2]) : w.ratinv[1]
                r2 = r.present[3] ? Float32(r.values[3]) : w.ratinv[2]
                w.ratinv = (r1, r2)
            else
                _wpbr_sched!(w, idt, 1006, _wpbr_prms(r, 2, 3))
            end
        elseif k == "COMMENT"                      # option 8 — comment lines until END
            while !eof(kr.io)
                rl = read_raw_line!(kr)
                startswith(uppercase(strip(String(rl))), "END") && break
            end
        elseif k == "GROWRATE"                     # option 9 — canker growth rates
            i3 = (r.present[1] && r.values[1] != 0.0f0) ? trunc(Int, r.values[1]) : 0
            i4 = (i3 >= 1 && i3 <= length(w.brspm)) ? w.brspm[i3] : 0
            (i3 < 0 || i3 > length(w.brspm) || (i4 == 0 && i3 != 0)) && continue   # invalid/non-host
            i5 = trunc(Int, r.values[2])
            (i5 < 0 || i5 > 4) && continue
            botmp = (r.present[3] && 0.1f0 <= r.values[3] <= 10.0f0) ? Float32(r.values[3]) : 0.0f0
            brtmp = (r.present[4] && 0.1f0 <= r.values[4] <= 10.0f0) ? Float32(r.values[4]) : 0.0f0
            _wpbr_growrate!(w, i4, i5, botmp, brtmp)
        elseif k == "BRTLST"                       # option 10 — tree-list output
            w.brtl = true
            r.present[2] && (w.idtout = Int32(trunc(Int, r.values[2])))
            r.present[3] && read_raw_line!(kr)     # supplemental file-name record
            _wpbr_sched!(w, r.present[1] ? trunc(Int, r.values[1]) : 0, 1007, Float32[])
        elseif k == "BRCLST"                       # option 11 — canker-list output
            w.brcl = true
            r.present[2] && (w.idcout = Int32(trunc(Int, r.values[2])))
            r.present[3] && read_raw_line!(kr)     # supplemental file-name record
            _wpbr_sched!(w, r.present[1] ? trunc(Int, r.values[1]) : 0, 1008, Float32[])
        elseif k == "CANFMT"                        # option 12 — canker read format
            eof(kr.io) || (w.icfmt = String(read_raw_line!(kr)))
        elseif k == "DEVFACT"                       # option 13 — deviation factors
            idt = r.present[1] ? trunc(Int, r.values[1]) : 0
            ksp = r.present[2] ? trunc(Int, r.values[2]) : 0
            (ksp < 1 || ksp > length(w.brspm) || w.brspm[ksp] == 0) && continue
            # FVS brin.f DEVFACT indexes DFACT(KSP,·) by the raw FVS species code,
            # but DFACT is dimensioned (NBRSP,4) — a latent Fortran array overrun
            # that only coincides with the BR host index for host codes ≤ NBRSP.
            # We index by the in-bounds BR host index (the conceptual dimension);
            # inert-seam ⇒ no .sum effect. See handoff for the dynamics-chunk audit.
            bi = w.brspm[ksp]
            if idt == 0
                w.ldfact = true
                for (fld, col) in ((3, 1), (4, 2), (5, 3), (6, 4))
                    if r.present[fld] && 0.01f0 <= r.values[fld] <= 6.0f0
                        w.dfact[bi, col] = Float32(r.values[fld])
                    end
                end
            else
                _wpbr_sched!(w, idt, 1010, _wpbr_prms(r, 2, 6))
            end
        elseif k == "RUSTINDX"                       # option 14 — rust index settings
            r.present[1] && (w.ridef = Float32(r.values[1]))
            rimeth = w.rimeth
            if r.present[3]
                rm = trunc(Int, r.values[3])
                rimeth = (1 <= rm <= 4) ? Int32(rm) : Int32(0)
            end
            w.rimeth = rimeth
            if rimeth == 1 || rimeth == 2
                r.present[4] && (w.ribprp = (Float32(r.values[4]), w.ribprp[2], w.ribprp[3]))
                r.present[5] && (w.ribprp = (w.ribprp[1], Float32(r.values[5]), w.ribprp[3]))
                r.present[6] && (w.ribprp = (w.ribprp[1], w.ribprp[2], Float32(r.values[6])))
                w.ribprp = _wpbr_normprp(w.ribprp)
            elseif rimeth == 3 || rimeth == 4
                # Gaussian (3) / log (4) exposure-time RI: defaults, then a
                # supplemental record of 4 params (0 keeps the default).
                if rimeth == 3
                    w.minri = 1.238867f-3; w.maxri = 7.222711f-3; w.pkage = 23.34596387f0; w.pkshp = 4.078136888f0
                else
                    w.minri = 9.6f-5; w.maxri = 4.7f-3; w.pkage = 1.66f0; w.pkshp = 1.086f0
                end
                v = _wpbr_read_params(kr, 4)
                v !== nothing && begin
                    v[1] > 0.0f0 && (w.minri = v[1]); v[2] > 0.0f0 && (w.maxri = v[2])
                    v[3] > 0.0f0 && (w.pkage = v[3]); v[4] > 0.0f0 && (w.pkshp = v[4])
                end
            end
        elseif k == "BRSEED"                          # option 15 — reseed the LCG
            wpbr_seed!(w, Float32(r.values[1]), r.present[1])
        elseif k == "STOCK"                           # option 16 — stock-type mix
            idt = r.present[1] ? trunc(Int, r.values[1]) : 0
            itsp = trunc(Int, r.values[2]); icls = trunc(Int, r.values[3])
            (itsp < 1 || itsp > length(w.brspm) || icls < 1 || icls > 4 || w.brspm[itsp] == 0) && continue
            # Same raw-KSP overrun as DEVFACT (brin.f STOCK: PRPSTK/RESIST(ITSP,·));
            # index by the in-bounds BR host index. Inert-seam ⇒ no .sum effect.
            bi = w.brspm[itsp]
            if idt == 0
                r.present[4] && (w.prpstk[bi, icls] = Float32(r.values[4]))
                r.present[5] && (w.resist[bi, icls] = Float32(r.values[5]))
            else
                _wpbr_sched!(w, idt, 1009, _wpbr_prms(r, 2, 5))
            end
        elseif k == "CANKDATA"                        # option 17 — read canker data
            r.present[1] && (w.icin = Int32(trunc(Int, r.values[1])))
            (r.present[2] && trunc(Int, r.values[2]) == 1) && (w.lmetric = true)
            (r.present[3] && trunc(Int, r.values[3]) == 1) && (w.ckinit = true)
            # ICIN=15 → inline data after the keyword; else a file-name record.
            # BRCANK canker ingestion is deferred to the dynamics chunk; here we
            # only consume the single supplemental record when a file is named.
            if !(r.present[1] && trunc(Int, r.values[1]) == 15)
                read_raw_line!(kr)
            end
        elseif k == "BROUT"                            # option 18 — output control
            r.present[1] && (w.lbrsum = trunc(Int, r.values[1]) != 0)
            r.present[2] && (w.lbrdbh = trunc(Int, r.values[2]) != 0)
        else
            # BRUST sub-keyword not recognized — record it so it can't hide as a
            # silent gap (mirrors the DFTM/RD keyword-dispatch policy).
            (!isempty(k) && isletter(first(k))) && push!(s.control.unrecognized_keywords, k)
        end
    end
    return nothing
end

# Record a scheduled management activity (OPNEW) for the future engine seam.
@inline function _wpbr_sched!(w::WpbrState, idt::Integer, iact::Integer, prms::Vector{Float32})
    push!(w.activities, (Int(idt), Int(iact), prms))
    return nothing
end

# Collect present fields lo..hi into a Float32 PRMS vector (0 for absent — the
# activity payload; the inert chunk never reads it, so exact packing is deferred).
@inline function _wpbr_prms(r, lo::Int, hi::Int)
    Float32[r.present[i] ? Float32(r.values[i]) : 0.0f0 for i in lo:hi]
end

# BRIN GROWRATE inner loops (brin.f stmt 900): set BOGRTH/BRGRTH over the selected
# species (i4==0 ⇒ all NBRSP) and stock types (i5==0 ⇒ all 4); only non-zero temps
# overwrite.
function _wpbr_growrate!(w::WpbrState, i4::Int, i5::Int, botmp::Float32, brtmp::Float32)
    sps = i4 == 0 ? (1:WPBR_NBRSP) : (i4:i4)
    sts = i5 == 0 ? (1:4) : (i5:i5)
    for sp in sps, st in sts
        botmp != 0.0f0 && (w.bogrth[sp, st] = botmp)
        brtmp != 0.0f0 && (w.brgrth[sp, st] = brtmp)
    end
    return nothing
end

# RUSTINDX RIBPRP normalization (brin.f label 1405): clamp negatives to 0, then
# spread the 1.0-sum deficit across the three ribes proportions until within 1%.
function _wpbr_normprp(p::NTuple{3,Float32})::NTuple{3,Float32}
    a, b, c = max(p[1], 0.0f0), max(p[2], 0.0f0), max(p[3], 0.0f0)
    while true
        excess = 1.0f0 - a - b - c
        (excess > 0.01f0 || excess < -0.01f0) || break
        a += excess / 3.0f0; b += excess / 3.0f0; c += excess / 3.0f0
    end
    return (a, b, c)
end

# Read a supplemental unformatted record of `n` numeric params (RUSTINDX 3/4).
function _wpbr_read_params(kr::KeywordReader, n::Int)
    eof(kr.io) && return nothing
    rl = read_raw_line!(kr)
    toks = split(String(rl))
    v = zeros(Float32, n)
    for i in 1:min(n, length(toks))
        p = tryparse(Float32, toks[i])
        p !== nothing && (v[i] = p)
    end
    return v
end
