# =============================================================================
# Westwide Pine Beetle (WWPB) — LANDSCAPE model (synthetic PPE harness).
# =============================================================================
# USER-approved 2026-08-21 (see memory fvsjl-wwpb-ppe-oracle-absent): build an
# end-to-end WWPB outbreak in FVSjl. The beetle KERNELS below are ported
# bit-exact from the pristine wwpb/*.f (ground truth, driver-golden-validated);
# only the absent outer orchestration (PPMAIN/ALSTD2/SPLAEX/GPGET) is faithfully
# reconstructed. Full architecture + plan: scratchpad/wwpb/HARNESS_PORT_PLAN.md.
#
# This file builds up the landscape state + kernels chunk-by-chunk. It is
# ADDITIVE and INERT (no simulate.jl seam) until the harness is wired — a stand
# projects byte-identically until then (gate: multicycle 339/11).
#
# CHUNK 1-2 (this commit): BMDBHC (DBH→size-class), the per-stand landscape
# state, and BMSDIT (the FVS→BM tree bridge — bins the FVSjl treelist into the
# NSCL size-class × host/nonhost table). Float32 throughout (WWPB's Fortran REAL).
# =============================================================================

const WWPB_NSCL = 10                       # BMPRM NSCL — number of DBH size classes
const WWPB_NUMRV = 9                        # BMPRM NUMRV — number of driving-variable rating values
const WWPB_PI24 = 3.14159f0 / (24.0f0 * 24.0f0)   # bmsdit.f/bminit.f PI24 = PIE/(24·24), Float32 (PIE=REAL param)

# glibc single-precision transcendentals — match gfortran REAL EXP()/x**y bit-exact
# (WWPB's Fortran is all-REAL, so the operations are expf/powf, not the Float64 forms).
@inline _wwpb_expf(x::Float32)::Float32 = ccall((:expf, "libm.so.6"), Float32, (Float32,), x)
@inline _wwpb_powf(x::Float32, y::Float32)::Float32 = ccall((:powf, "libm.so.6"), Float32, (Float32, Float32), x, y)

# -----------------------------------------------------------------------------
# BMDBHC (bmdbhc.f) — assign a DBH to a size class from the UPSIZ breakpoints.
# INDEX=1; DO I=1,NSCL-1: IF DBH < UPSIZ(I) GOTO 20; INDEX=INDEX+1. (1..NSCL)
# -----------------------------------------------------------------------------
@inline function wwpb_dbh_class(upsiz::Vector{Float32}, dbh::Real)::Int
    idx = 1
    @inbounds for i in 1:(WWPB_NSCL - 1)
        Float32(dbh) < upsiz[i] && return idx
        idx += 1
    end
    return idx
end

# -----------------------------------------------------------------------------
# ISPFLL — standing-dead falldown-rate class (1=fast,2=med,3=slow) per species.
# From bmblkd*.f. IE/N-Rockies 11-sp order (WP L DF GF WH C LP S AF PP OTH):
const WWPB_ISPFLL_NI = Int32[2, 3, 2, 1, 2, 3, 2, 1, 2, 1, 2]

# HSPEC natural host designations (bmblkd*.f commented defaults, set live by the
# HOST keyword): MPB(1)+Ips(3) host = LP (sp7 in the IE list); WPB(2) host = PP
# (sp10). Represented as FVS species *codes* so it is variant-list-independent;
# the bridge maps the stand's species index → code to test host membership.
# Natural pine hosts by 2-letter alpha code: LP lodgepole, PP ponderosa.
const WWPB_HOST_MPB_ALPHA = ("LP",)          # PBSPEC 1 (MPB) + 3 (Ips)
const WWPB_HOST_WPB_ALPHA = ("PP",)          # PBSPEC 2 (WPB)

# -----------------------------------------------------------------------------
# WwpbStand — the per-stand BM landscape state (BMCOM, one stand). Size-class ×
# {host=1, nonhost=2}. NSCL+1 slot is the stand summary. Populated by bmsdit!.
# -----------------------------------------------------------------------------
mutable struct WwpbStand
    # --- from the FVS→BM bridge (bmsdit.f) ---
    tree ::Matrix{Float32}      # TREE(NSCL+1, 2)   — TPA per size class × type
    bah  ::Vector{Float32}      # BAH(NSCL+1)       — host basal area per class + summary
    banh ::Vector{Float32}      # BANH(NSCL+1)      — nonhost basal area
    hts  ::Matrix{Float32}      # HTS(NSCL, 2)      — mean height per class × type
    crs  ::Matrix{Float32}      # CRS(NSCL, 2)      — mean crown % per class × type
    hgs  ::Matrix{Float32}      # HGS(NSCL, 2)      — mean height growth per class × type
    tvol ::Matrix{Float32}      # TVOL(NSCL, 2)     — mean cuft volume per class × type
    otpa ::Matrix{Float32}      # OTPA(NSCL, 2)     — initial TPA snapshot
    isph ::Vector{Int32}        # ISPH(2)           — dominant host/nonhost species index
    iqptyp::Vector{Int32}       # IQPTYP(2)         — dead-wood quality pool type
    # --- beetle ledgers (zeroed at setup; filled by the per-year kernels) ---
    oakill::Matrix{Float32}     # OAKILL(NSCL, 2)   — "other agent" kill density
    pbkill::Vector{Float32}     # PBKILL(NSCL)      — pine-beetle-killed proportion
    allkll::Vector{Float32}     # ALLKLL(NSCL+1)    — Ips-killed stems
    scorch::Vector{Float32}     # SCORCH(NSCL)      — severely-scorched proportion
    tpbk ::Array{Float32,3}     # TPBK(NSCL,2,3)    — mortality ledger [class,type,{fast,slow,beetle}]
    fastk::Vector{Float32}      # FASTK(3)          — fast-kill totals [TPA, host vol+BA, nonhost BA]
    # --- GRF / rating-value state (bmcgrf.f) ---
    grf  ::Vector{Float32}      # GRF(NSCL)         — per-size-class growth-reduction factor
    grfstd::Float32             # GRFSTD            — BA-weighted stand GRF
    rvdnst::Float32             # RVDNST            — stand-density rating value
    bastd::Float32              # BASTD             — total stand BA (host+nonhost)
    dvrv ::Vector{Float32}      # DVRV(NUMRV=9)     — per-driving-variable rating values
    # --- stressor inputs (from the DM/RD/rust/other-beetle/lightning/drought/defoliator models;
    #     0 or 1 when those models are off — a pure beetle outbreak drives via RVDNST) ---
    sdmr ::Vector{Float32}      # SDMR(NSCL)   — dwarf-mistletoe rating
    srr  ::Vector{Float32}      # SRR(NSCL)    — root-disease proportion
    ssr  ::Vector{Float32}      # SSR(NSCL)    — stem-rust proportion
    othatt::Vector{Float32}     # OTHATT(NSCL) — other-beetle attack proportion
    topkll::Vector{Float32}     # TOPKLL(NSCL) — Ips top-kill proportion
    strike::Vector{Float32}     # STRIKE(NSCL) — lightning-strike proportion
    rvdsc::Vector{Float32}      # RVDSC(NSCL)  — drought rating value
    rvdfol::Vector{Float32}     # RVDFOL(NSCL) — defoliator rating value
    # --- BKP (beetle-killing-potential / brood) state (bmcbkp.f) ---
    bkp   ::Float32             # BKP    — stand beetle-killing potential (non-Ips)
    bkpips::Float32             # BKPIPS — Ips beetle-killing potential
    oldbkp::Float32             # OLDBKP — BKP before dispersal (feedback ref)
    final ::Vector{Float32}     # FINAL(3) — last tree killed [bkp used, size class, tpa killed]
    strip ::Vector{Float32}     # STRIP(NSCL) — strip-attack proportion
    pslash::Matrix{Float32}     # PSLASH(MXDWHC=2, MXDWSZ=2) — slash colonized by Ips
    # --- special-tree / attractiveness state (bmcspt.f, bmcnum.f) ---
    spclt ::Matrix{Float32}     # SPCLT(NSCL,2) — proportion of "special" (attractive) trees per class × pass
    pitch ::Vector{Float32}     # PITCH(NSCL)  — pitch-out / strip-kill proportion
    atrphe::Float32             # ATRPHE — proportion of trees with attractant pheromone
    numer ::Vector{Float32}     # NUMER(2) — scoring-equation numerator (main pest, Ips)
    tfood ::Vector{Float32}     # TFOOD(2) — total beetle attractive "food"
    repphe::Float32             # REPPHE — repellent-pheromone factor (default → 1)
    ssbatk::Float32             # SSBATK — special-tree BA already attacked (0 initially)
    dwphos::Matrix{Float32}     # DWPHOS(2, MXDWHZ=3) — downed/standing host volume for Ips [type, size]
    # --- kill / dead-wood state (bmistd.f) ---
    sdwp ::Array{Float32,3}     # SDWP(MXDWPC=3, MXDWHZ+1=4, MXDWAG=5) — standing dead-wood volume pool
    spray::Float32              # SPRAY — proportion of beetles killed by spraying (default 0)
    slp  ::Float32              # SLP — stand slope
    habtyp::Int32               # HABTYP — habitat/ecoclass code (fire model)
end

function WwpbStand()
    z2  = () -> zeros(Float32, WWPB_NSCL, 2)
    return WwpbStand(
        zeros(Float32, WWPB_NSCL + 1, 2),                 # tree
        zeros(Float32, WWPB_NSCL + 1), zeros(Float32, WWPB_NSCL + 1),  # bah, banh
        z2(), z2(), z2(), z2(), z2(),                     # hts, crs, hgs, tvol, otpa
        zeros(Int32, 2), zeros(Int32, 2),                 # isph, iqptyp
        z2(),                                             # oakill
        zeros(Float32, WWPB_NSCL), zeros(Float32, WWPB_NSCL + 1), zeros(Float32, WWPB_NSCL),  # pbkill, allkll, scorch
        zeros(Float32, WWPB_NSCL, 2, 3), zeros(Float32, 3),  # tpbk, fastk
        zeros(Float32, WWPB_NSCL), 1.0f0, 1.0f0, 0.0f0, zeros(Float32, WWPB_NUMRV),  # grf, grfstd, rvdnst, bastd, dvrv
        zeros(Float32, WWPB_NSCL), zeros(Float32, WWPB_NSCL), zeros(Float32, WWPB_NSCL),  # sdmr, srr, ssr
        zeros(Float32, WWPB_NSCL), zeros(Float32, WWPB_NSCL), zeros(Float32, WWPB_NSCL),  # othatt, topkll, strike
        zeros(Float32, WWPB_NSCL), zeros(Float32, WWPB_NSCL),  # rvdsc, rvdfol
        0.0f0, 0.0f0, 0.0f0, zeros(Float32, 3), zeros(Float32, WWPB_NSCL), zeros(Float32, 2, 2),  # bkp, bkpips, oldbkp, final, strip, pslash
        zeros(Float32, WWPB_NSCL, 2), zeros(Float32, WWPB_NSCL), 0.0f0,  # spclt, pitch, atrphe
        zeros(Float32, 2), zeros(Float32, 2), 0.0f0, 0.0f0, zeros(Float32, 2, 3),  # numer, tfood, repphe, ssbatk, dwphos
        zeros(Float32, 3, 4, 5), 0.0f0,                   # sdwp, spray
        0.0f0, Int32(0),                                  # slp, habtyp
    )
end

# -----------------------------------------------------------------------------
# host designation: is FVS species index `sp` a host for beetle `pbspec`?
# Natural defaults keyed on the species' 2-letter alpha code (LP/PP). `alpha`
# is the per-species code vector (coef.code2-style, 2-char upper).
# -----------------------------------------------------------------------------
@inline function wwpb_is_host(pbspec::Integer, sp_alpha::AbstractString)::Bool
    a = uppercase(strip(sp_alpha))
    mpb = a in WWPB_HOST_MPB_ALPHA
    wpb = a in WWPB_HOST_WPB_ALPHA
    (pbspec == 1 || pbspec == 4) && mpb && return true   # MPB host
    (pbspec == 2 || pbspec == 4) && wpb && return true   # WPB host
    (pbspec == 3) && mpb && return true                  # Ips host (= LP)
    return false
end

# -----------------------------------------------------------------------------
# bmsdit! (bmsdit.f) — FVS→BM tree BRIDGE. Bins the FVSjl treelist into the
# size-class × host/nonhost table, sets ISPH/IQPTYP, averages HTS/CRS/HGS/TVOL.
# `sp_alpha(sp)` returns the 2-letter species code for host testing; `ispfll`
# is the falldown-rate map (defaults WWPB_ISPFLL_NI). Faithful to bmsdit.f DO-20.
# -----------------------------------------------------------------------------
function bmsdit!(st::WwpbStand, t, w::WwpbState, sp_alpha::Function;
                 ispfll::Vector{Int32}=WWPB_ISPFLL_NI, habtyp::Integer=0, slope::Real=0.0)
    pbspec = Int(w.pbspec)
    upsiz  = w.upsiz
    st.habtyp = Int32(habtyp)
    st.slp    = Float32(slope)

    # zero the accumulators (bmsdit.f DO-100)
    fill!(st.tree, 0.0f0); fill!(st.bah, 0.0f0); fill!(st.banh, 0.0f0)
    fill!(st.hts, 0.0f0); fill!(st.crs, 0.0f0); fill!(st.hgs, 0.0f0); fill!(st.tvol, 0.0f0)

    # per-species host BA accumulator for ISPH (max-BA dominant)
    totbas = Dict{Int,NTuple{2,Float32}}()   # sp → (host_BA, nonhost_BA)

    @inbounds for i in 1:t.n
        t.tpa[i] > 0.0f0 || continue
        sp = Int(t.species[i])
        d  = t.dbh[i]
        K  = wwpb_dbh_class(upsiz, d)
        batemp = d * d * WWPB_PI24
        x = batemp * t.tpa[i]
        lx = wwpb_is_host(pbspec, sp_alpha(sp))
        ty = lx ? 1 : 2
        # BA (host/nonhost) + summary slot
        if lx
            st.bah[K] += x;  st.bah[WWPB_NSCL+1] += x
        else
            st.banh[K] += x; st.banh[WWPB_NSCL+1] += x
        end
        st.tree[K, ty] += t.tpa[i]
        st.tree[WWPB_NSCL+1, ty] += t.tpa[i]
        st.hts[K, ty]  += t.height[i]        * t.tpa[i]
        st.crs[K, ty]  += Float32(t.crown_pct[i]) * t.tpa[i]
        st.hgs[K, ty]  += t.ht_growth[i]     * t.tpa[i]
        st.tvol[K, ty] += t.cuft_vol[i]      * t.tpa[i]
        cur = get(totbas, sp, (0.0f0, 0.0f0))
        totbas[sp] = lx ? (cur[1] + x, cur[2]) : (cur[1], cur[2] + x)
    end

    # ISPH = species with max host/nonhost BA (bmsdit.f DO-25)
    maxh = 0.0f0; maxn = 0.0f0; st.isph[1] = 0; st.isph[2] = 0
    for (sp, (bh, bn)) in totbas
        bn > maxn && (maxn = bn; st.isph[2] = Int32(sp))
        bh > maxh && (maxh = bh; st.isph[1] = Int32(sp))
    end
    for j in 1:2
        s = Int(st.isph[j])
        st.iqptyp[j] = (s > 0 && s <= length(ispfll)) ? ispfll[s] : Int32(1)
    end

    # OTPA snapshot + size-class averages (bmsdit.f DO-30/DO-40)
    @inbounds for k in 1:WWPB_NSCL
        st.otpa[k, 1] = st.tree[k, 1]; st.otpa[k, 2] = st.tree[k, 2]
        for ty in 1:2
            if st.tree[k, ty] > 1.0f-9
                st.hts[k, ty]  /= st.tree[k, ty]
                st.crs[k, ty]  /= st.tree[k, ty]
                st.hgs[k, ty]  /= st.tree[k, ty]
                st.tvol[k, ty] /= st.tree[k, ty]
            else
                st.hts[k, ty] = 0.0f0; st.crs[k, ty] = 0.0f0
                st.hgs[k, ty] = 0.0f0; st.tvol[k, ty] = 0.0f0
            end
        end
    end
    return st
end

# -----------------------------------------------------------------------------
# bmmort! (bmmort.f) — remove killed trees from the size-class table. Called
# TWICE per year: SLOW=false (fast agents: windthrow+fire OAKILL) then SLOW=true
# (beetle PBKILL + slow OAKILL + Ips ALLKLL). Decrements TREE/BAH/BANH, tallies
# the FASTK + TPBK(,,{fast=1,slow=2,beetle=3}) ledgers. Faithful to bmmort.f DO-800.
# NOTE OAKILL enters as a PROPORTION and is converted to TPA (×TREE) here, used,
# then zeroed; PBKILL/ALLKLL persist (zeroed after BKP emerges). Beetles (slow)
# kill HOST only. Deterministic — no transcendentals.
# -----------------------------------------------------------------------------
function bmmort!(st::WwpbStand, slow::Bool)
    st.bah[WWPB_NSCL+1] = 0.0f0; st.banh[WWPB_NSCL+1] = 0.0f0
    st.tree[WWPB_NSCL+1, 1] = 0.0f0; st.tree[WWPB_NSCL+1, 2] = 0.0f0
    @inbounds for k in 1:WWPB_NSCL
        st.oakill[k, 1] *= st.tree[k, 1]       # proportion → TPA
        st.oakill[k, 2] *= st.tree[k, 2]
        if slow
            prdead = st.tree[k, 1] > 0.0f0 ?
                (st.pbkill[k] + st.oakill[k, 1] + st.allkll[k]) / st.tree[k, 1] : 0.0f0
            st.bah[k]     *= (1.0f0 - prdead)
            st.tree[k, 1] *= (1.0f0 - prdead)
        else
            st.fastk[1] += st.oakill[k, 1] + st.oakill[k, 2]
            st.fastk[2] += st.oakill[k, 1] * st.tvol[k, 1] + st.oakill[k, 2] * st.tvol[k, 2]
            prdead = st.tree[k, 1] > 0.0f0 ? st.oakill[k, 1] / st.tree[k, 1] : 0.0f0
            st.fastk[2] += st.bah[k] * prdead
            st.bah[k]     *= (1.0f0 - prdead)
            st.tree[k, 1] -= st.oakill[k, 1]
            prdead = st.tree[k, 2] > 0.0f0 ? st.oakill[k, 2] / st.tree[k, 2] : 0.0f0
            st.fastk[3] += st.banh[k] * prdead
            st.banh[k]    *= (1.0f0 - prdead)
            st.tree[k, 2] -= st.oakill[k, 2]
        end
        # constrain to positive, zero BA if trees gone
        st.bah[k] = max(st.bah[k], 0.0f0); st.banh[k] = max(st.banh[k], 0.0f0)
        st.tree[k, 1] = max(st.tree[k, 1], 0.0f0); st.tree[k, 2] = max(st.tree[k, 2], 0.0f0)
        st.tree[k, 1] <= 0.0f0 && (st.bah[k] = 0.0f0)
        st.tree[k, 2] <= 0.0f0 && (st.banh[k] = 0.0f0)
        st.bah[WWPB_NSCL+1]  += st.bah[k];  st.banh[WWPB_NSCL+1] += st.banh[k]
        st.tree[WWPB_NSCL+1, 1] += st.tree[k, 1]; st.tree[WWPB_NSCL+1, 2] += st.tree[k, 2]
        if slow
            st.tpbk[k, 1, 3] += st.pbkill[k] + st.allkll[k]
            st.tpbk[k, 1, 2] += st.oakill[k, 1]
            st.tpbk[k, 2, 2] += st.oakill[k, 2]
        else
            st.tpbk[k, 1, 1] += st.oakill[k, 1]
            st.tpbk[k, 2, 1] += st.oakill[k, 2]
        end
        st.oakill[k, 1] = 0.0f0; st.oakill[k, 2] = 0.0f0
    end
    return st
end

# -----------------------------------------------------------------------------
# bmcgrf! (bmcgrf.f) — compute the Growth-Reduction Factor (GRF) per host size
# class by multiplying every stressor rating (DM, root disease, stem rust, other
# beetles, fire, lightning, drought, defoliators), each floored at 0.01, product
# floored at 0.01. Then GRFSTD = BA-weighted mean, and the stand-density rating
# RVDNST = 2 − 1.9/(1 + 9·exp(BASTD·−0.033))^3 (applied to GRFSTD if LCDENS).
# `oldgrf` receives last year's GRF (for reproduction). Faithful to bmcgrf.f.
# The one transcendental (exp) + the real power (^3.0) route through glibc.
# -----------------------------------------------------------------------------
function bmcgrf!(st::WwpbStand, w::WwpbState, oldgrf::Vector{Float32}; lcdens::Bool=true)
    pbspec = Int(w.pbspec)
    fill!(st.dvrv, 0.0f0)
    @inbounds for icls in 1:WWPB_NSCL
        oldgrf[icls] = st.grf[icls] > 0.0f0 ? st.grf[icls] : 1.0f0
        g = 1.0f0
        bah = st.bah[icls]
        # DM (SDMR near 0 ⇒ DMR=6)
        gd = 1.0f0 - (st.sdmr[icls] / 6.5f0); gd <= 0.0f0 && (gd = 0.01f0)
        g *= gd; st.dvrv[6] += gd * bah
        # root disease
        gd = 1.0f0 - st.srr[icls]; gd <= 0.0f0 && (gd = 0.01f0)
        g *= gd; st.dvrv[7] += gd * bah
        # stem rust
        gd = 1.0f0 - st.ssr[icls]; gd <= 0.0f0 && (gd = 0.01f0)
        g *= gd; st.dvrv[8] += gd * bah
        # other beetle attacks
        gd = 1.0f0 - (st.othatt[icls] + st.topkll[icls]); gd <= 0.0f0 && (gd = 0.01f0)
        g *= gd; st.dvrv[4] += gd * bah
        # fire (beetle-species-dependent)
        gd = 1.0f0
        if pbspec == 1 || pbspec == 2
            gd = 1.0f0 - 0.5f0 * st.scorch[icls]
        elseif pbspec == 3
            gd = 1.0f0 - 0.99f0 * st.scorch[icls]
        end
        gd <= 0.0f0 && (gd = 0.01f0)
        g *= gd; st.dvrv[1] += gd * bah
        # lightning
        gd = 1.0f0 - 0.99f0 * st.strike[icls]; gd <= 0.0f0 && (gd = 0.01f0)
        g *= gd; st.dvrv[2] += gd * bah
        # drought
        gd = st.rvdsc[icls]; gd <= 0.0f0 && (gd = 0.01f0)
        g *= gd; st.dvrv[3] += gd * bah
        # defoliators (default 1 the first year so GRF ≠ 0)
        st.rvdfol[icls] <= 0.0f0 && (st.rvdfol[icls] = 1.0f0)
        gd = st.rvdfol[icls]; gd <= 0.0f0 && (gd = 0.01f0)
        g *= gd; st.dvrv[5] += gd * bah
        st.grf[icls] = max(0.01f0, g)
    end

    total = 0.0f0; st.grfstd = 0.0f0
    @inbounds for icls in 1:WWPB_NSCL
        x = st.bah[icls]
        total += x
        st.grfstd += st.grf[icls] * x
    end
    if total > 1.0f-9
        st.grfstd /= total
        @inbounds for idv in 1:WWPB_NUMRV; st.dvrv[idv] /= total; end
    else
        st.grfstd = 1.0f0
        @inbounds for idv in 1:WWPB_NUMRV; st.dvrv[idv] = 1.0f0; end
    end

    # stand-density rating value RVDNST = 2 − 1.9/(1 + 9·exp(BASTD·−0.033))^3
    dncf1 = -0.033f0; dncf2 = 3.0f0
    st.bastd = st.bah[WWPB_NSCL+1] + st.banh[WWPB_NSCL+1]
    base = 1.0f0 + 9.0f0 * _wwpb_expf(st.bastd * dncf1)
    st.rvdnst = 2.0f0 - (1.9f0 / _wwpb_powf(base, dncf2))
    if lcdens
        st.grfstd *= st.rvdnst
        st.dvrv[9] = st.rvdnst
    else
        st.dvrv[9] = 1.0f0
    end
    return st
end

# -----------------------------------------------------------------------------
# wwpb_init_coeffs (bminit.f) — the size-class BA + brood "increase" coefficients
# that BMCBKP/BMCNUM/BMISTD consume. MSBA(i)=MID²·(π/576) (BA of the size-class-
# midpoint DBH); UPBA(i)=UPSIZ(i)²·(π/576); INC(1,i)=RSLOPE·DBHMID+B clamped to
# REPMAX above DBHMAX, INC(2,i)=INC(1,i), INC(3,i)=INC(1,1)·0.1 [Ips]. Defaults
# RSLOPE=1/10, REPMAX=4, DBHMAX=36, REPLAC=6 ⇒ B=1−0.6=0.4 (keyword-overridable).
# Depends only on UPSIZ ⇒ model-level (same for all stands). BIT-EXACT vs
# gfortran-16 (scratchpad/wwpb/driver_bminit.f). Integer UPSIZ arithmetic where
# the Fortran uses INTEGER UPSIZ (FLOAT(UPSIZ+LOW), UPSIZ²).
# -----------------------------------------------------------------------------
function wwpb_init_coeffs(upsiz::Vector{Float32};
                          rslope::Float32=Float32(1.0/10.0), repmax::Float32=4.0f0,
                          dbhmax::Float32=36.0f0, replac::Float32=6.0f0)
    msba = zeros(Float32, WWPB_NSCL); upba = zeros(Float32, WWPB_NSCL)
    inc  = zeros(Float32, 3, WWPB_NSCL)
    # WPBA(MXDWSZ=2) — dead-wood-pool BA from WPSIZ breakpoints (bminit.f:128-136); LOW(1)=3.
    wpsiz = WWPB_WPSIZ_DEFAULT   # 10, 20, 60
    wpba = zeros(Float32, 2)
    @inbounds for i in 1:2
        low = i == 1 ? 3 : Int(wpsiz[i-1])
        mid = Float32(Int(wpsiz[i]) + low) * 0.5f0
        wpba[i] = mid * mid * WWPB_PI24
    end
    b = 1.0f0 - (rslope * replac)
    @inbounds for i in 1:WWPB_NSCL
        iu  = Int(upsiz[i])
        low = i == 1 ? 0 : Int(upsiz[i-1])
        mid = Float32(iu + low) * 0.5f0
        msba[i] = mid * mid * WWPB_PI24
        upba[i] = Float32(iu * iu) * WWPB_PI24
        dbhmid = Float32(iu + low) / 2.0f0
        inc[1, i] = dbhmid < dbhmax ? (rslope * dbhmid + b) : repmax
        inc[2, i] = inc[1, i]
    end
    @inbounds for i in 1:WWPB_NSCL
        inc[3, i] = inc[1, 1] * 0.1f0
    end
    # L2D (bminit.f:141-158) — living size class → dead-wood host zone (0 if ≤3").
    l2d = zeros(Int, WWPB_NSCL); jsiz = 1
    @inbounds for isiz in 1:WWPB_NSCL
        if Int(upsiz[isiz]) <= 3
            l2d[isiz] = 0
        else
            dif = (upsiz[isiz] - wpsiz[jsiz]) + (upsiz[isiz-1] - wpsiz[jsiz])
            if upsiz[isiz] <= wpsiz[jsiz]
                l2d[isiz] = jsiz
            elseif dif <= 0.0f0
                l2d[isiz] = jsiz
            else
                jsiz += 1; l2d[isiz] = jsiz
            end
        end
    end
    return (msba = msba, upba = upba, inc = inc, wpba = wpba, l2d = l2d)
end

# -----------------------------------------------------------------------------
# bmcbkp! (bmcbkp.f) — compute the beetle-killing-potential (BKP, "brood") from
# last year's kills. NORMAL (non-bad-year) path, PBSPEC≠3: for each size class,
# BKP += INC(pbspec,isiz)·MSBA·PBKILL (reproduction from filled trees) + the
# FINAL last-tree term + 0.75·STRIP (strip attacks); zeroes PBKILL/STRIP. Ips
# BKPIPS from TOPKLL/ALLKLL. Then BKP·=REPRD, BKPIPS·=REPRDI (generation mults
# from NBGEN/NIBGEN), the LFDBK negative-feedback cap, and the PSLASH·WPBA slash
# term. OLDBKP=BKP. Faithful to bmcbkp.f. The bad-reproduction-year branch
# (GPGET2(317) scheduler) is NOT yet wired — LBAD defaults false (the outbreak
# path); a future BADREP-keyword chunk adds it. Deterministic (no transcendentals).
# -----------------------------------------------------------------------------
function bmcbkp!(st::WwpbStand, w::WwpbState, coeffs;
                 nbgen::Int=1, nibgen::Int=2, ipson::Bool=false, ipsmin::Int=2,
                 wpba::Vector{Float32}=zeros(Float32, 2), lfdbk::Bool=false, tfdbk::Float32=0.0f0)
    pbspec = Int(w.pbspec)
    msba = coeffs.msba; inc = coeffs.inc
    slinc = 5.0f0
    reprd  = nbgen  == 1 ? 1.0f0 : nbgen  == 2 ? 1.5f0 : nbgen  == 3 ? 2.0f0 : 2.5f0
    reprdi = nibgen == 1 ? 1.0f0 : nibgen == 2 ? 1.5f0 : nibgen == 3 ? 2.0f0 : 2.5f0
    if pbspec == 3
        reprd != 0.0f0 && (reprdi = reprd)
        reprd = 0.0f0
    end
    @inbounds for isiz in 1:WWPB_NSCL
        if pbspec != 3
            if Int(st.final[2]) == isiz
                st.pbkill[isiz] -= st.final[3]
                st.pbkill[isiz] < 0.0f0 && (st.pbkill[isiz] = 0.0f0)
                st.bkp += inc[pbspec, isiz] * st.final[1] * st.final[3]
                st.final[1] = 0.0f0; st.final[2] = 0.0f0; st.final[3] = 0.0f0
            end
            st.bkp += msba[isiz] * st.pbkill[isiz] * inc[pbspec, isiz]
            st.pbkill[isiz] = 0.0f0
            st.bkp += msba[isiz] * inc[pbspec, isiz] * 0.75f0 * st.strip[isiz]
            st.strip[isiz] = 0.0f0
        end
        st.bkpips += inc[3, isiz] * st.topkll[isiz] * st.tree[isiz, 1] * msba[ipsmin]
        st.bkpips += inc[3, isiz] * st.allkll[isiz] * msba[ipsmin]
        st.allkll[isiz] = 0.0f0
    end
    st.bkp    *= reprd
    st.bkpips *= reprdi
    if lfdbk && st.bkp > st.oldbkp && st.oldbkp >= tfdbk
        st.bkp = tfdbk
    end
    reprdi *= slinc
    @inbounds for idtyp in 1:2, idsiz in 1:2
        st.bkpips += st.pslash[idtyp, idsiz] * wpba[idsiz] * reprdi
        st.pslash[idtyp, idsiz] = 0.0f0
    end
    st.oldbkp = pbspec != 3 ? st.bkp : st.bkpips
    return st
end

# -----------------------------------------------------------------------------
# bmcspt! (bmcspt.f) — proportion of "special" (beetle-attractive) trees in a
# size class: pitch-outs, lightning strikes, top-kill/other-attack, scorch,
# attractant pheromone (MPB/WPB pass 1; other-attack+pheromone for Ips pass 2).
# MEASURED (scratchpad/wwpb/driver_bmcspt.f, gfortran-16): the Fortran's nested
# inclusion-exclusion overlap correction (DO 1000..1115) is DEAD CODE — its
# innermost `DO 1115 JK=JK+1,SPTCNT` never executes (JK uninitialized/≥SPTCNT in
# the build), so no overlap is subtracted. Effective, bit-exact behavior:
#   all SP==0 → 0 ;  any SP==1 → 1 ;  else min(Σ SP, 1). Faithful to the RUNNING
# pristine routine (doctrine: measure, don't infer). Deterministic.
# -----------------------------------------------------------------------------
function bmcspt!(st::WwpbStand, w::WwpbState, isiz::Int, ipass::Int; ipson::Bool=false)
    pbspec = Int(w.pbspec)
    xsplt = 0.0f0
    if st.tree[isiz, 1] > 1.0f-9
        sp = (0.0f0, 0.0f0, 0.0f0, 0.0f0, 0.0f0); sptcnt = 0
        if pbspec == 1 && ipass == 1
            sp = (st.pitch[isiz], st.strike[isiz], st.topkll[isiz] + st.othatt[isiz], st.scorch[isiz], st.atrphe); sptcnt = 5
        elseif pbspec == 2 && ipass == 1
            sp = (st.pitch[isiz], st.strike[isiz], st.topkll[isiz] + st.othatt[isiz], st.scorch[isiz], st.atrphe); sptcnt = 5
        elseif pbspec == 3 || ipass == 2
            sp = (st.othatt[isiz], st.atrphe, 0.0f0, 0.0f0, 0.0f0); sptcnt = 2
        end
        topflag = true; botflag = false
        @inbounds for i in 1:sptcnt
            sp[i] == 1.0f0 && (topflag = false)
            sp[i] != 0.0f0 && (botflag = true)
        end
        if topflag && botflag
            @inbounds for ii in 1:sptcnt
                xsplt += sp[ii]          # overlap-correction inner loops are dead (measured)
            end
            xsplt > 1.0f0 && (xsplt = 1.0f0)
        elseif !botflag
            xsplt = 0.0f0                # all proportions 0
        else                             # !topflag ⇒ some proportion == 1
            xsplt = 1.0f0
        end
    end
    st.spclt[isiz, ipass] = xsplt
    return xsplt
end

# -----------------------------------------------------------------------------
# bmcnum! (bmcnum.f) — numerator of the between-stand "scoring" (attractiveness)
# equation. For each pass (main pest; +Ips if IPSON): BAIS = total stand BA; per
# size class compute special trees (bmcspt!) → SPCLT, accumulate SPEC, and above
# the min attack size class LISCMIN accumulate SPAREA (special BA) and BAHG (host
# BA). Ips slash adds DWPHOS·WPBA food. Then
#   NUMER(pass) = (USERA·SPEC + 1)·BAIS·(REPPHE·BAHG + SPAREA + SSBATK)/GRFSTD.
# Faithful to bmcnum.f. USERA defaults 1 (bminit.f:239); REPPHE→1 if 0. `coeffs`
# supplies MSBA/WPBA. Deterministic given SPCLT (bmcspt! is deterministic).
# -----------------------------------------------------------------------------
function bmcnum!(st::WwpbStand, w::WwpbState, coeffs;
                 ipson::Bool=false, usera::NTuple{3,Float32}=(1.0f0, 1.0f0, 1.0f0))
    pbspec = Int(w.pbspec)
    msba = coeffs.msba; wpba = coeffs.wpba
    liscmin = (Int(w.iscmin[pbspec]), 0)
    aspec   = (usera[pbspec], 0.0f0)
    st.numer[1] = 0.0f0; st.numer[2] = 0.0f0
    st.tfood[1] = 0.0f0; st.tfood[2] = 0.0f0
    npass = 1
    if ipson
        npass = 2
        liscmin = (liscmin[1], Int(w.iscmin[3]))
        aspec   = (aspec[1], usera[3])
    end
    @inbounds for ipass in 1:npass
        bait = st.bah[WWPB_NSCL+1] + st.banh[WWPB_NSCL+1]   # BAIS
        bait <= 1.0f-6 && continue
        spec = 0.0f0; bahg = 0.0f0; sparea = 0.0f0
        lmin = liscmin[ipass]
        for isiz in 1:WWPB_NSCL
            bmcspt!(st, w, isiz, ipass; ipson = ipson)
            sptree = st.tree[isiz, 1] * st.spclt[isiz, ipass]
            spec += sptree
            if isiz >= lmin
                sparea += sptree * msba[isiz]
                bahg   += st.bah[isiz]
            end
        end
        st.repphe == 0.0f0 && (st.repphe = 1.0f0)
        if ipass == 2 || pbspec == 3
            for idsiz in 1:2, idtyp in 1:2
                st.tfood[ipass] += st.dwphos[idtyp, idsiz] * wpba[idsiz]
                spec += st.dwphos[idtyp, idsiz]
            end
        end
        st.tfood[ipass] += bahg * (bahg / bait)
        st.numer[ipass] = (aspec[ipass] * spec + 1.0f0) * bait *
                          ((st.repphe * bahg) + sparea + st.ssbatk) / st.grfstd
    end
    return st
end

# -----------------------------------------------------------------------------
# bmatct_single! (bmatct.f, single-stand OUTOFF degenerate case) — the landscape
# BKP redistribution collapses to a self-allocation: with one stockable stand and
# the Outside World off, the between-stand SCORE cancels (PROP=1) and CAREA
# cancels (NEWBKP=CAREA·BKP then /CAREA), leaving the "how much BKP makes it into
# the stand" saturation:
#   BKP = TFOOD·(1 − exp(−ALPHA·BKP/(TFOOD+1e-6))),  ALPHA from drought SDD.
# Faithful to bmatct.f lines 707-746 for BMSTND=1, OUTOFF=T (verified vs
# driver_bmatct.f with SPLAAR/SPLALO/SPLADS/GPGET2 stubbed). The exp routes
# through glibc. The full multi-stand spatial dispersal (SPLADS distances, the
# Outside-World immigration BKPIN) is a later chunk — moot at MXSTND=1/OUTOFF.
# -----------------------------------------------------------------------------
function bmatct_single!(st::WwpbStand, w::WwpbState; sdd::Float32=0.0f0, ipson::Bool=false)
    pbspec = Int(w.pbspec)
    alpha = sdd < 0.0f0 ? (0.365f0 - 0.178f0 * sdd) : (0.365f0 - 0.034f0 * sdd)
    alpha > 0.90f0 && (alpha = 0.90f0)
    alpha < 0.01f0 && (alpha = 0.01f0)
    if pbspec != 3
        nb1 = st.bkp                       # NEWBKP(1) = self BKP (PROP=1, CAREA cancels)
        nb2 = ipson ? st.bkpips : 0.0f0    # NEWBKP(2) accumulated only when Ips is a DV
        if nb1 > 0.0f0 && st.tfood[1] > 0.0f0
            st.bkp = st.tfood[1] * (1.0f0 - _wwpb_expf(-(alpha * nb1 / (st.tfood[1] + 1.0f-6))))
        else
            st.bkp = 0.0f0
        end
        if nb2 > 0.0f0
            st.bkpips = st.tfood[2] * (1.0f0 - _wwpb_expf(-(alpha * nb2 / (st.tfood[2] + 1.0f-6))))
        else
            st.bkpips = 0.0f0
        end
    else
        nb1 = st.bkpips
        if nb1 > 0.0f0
            st.bkpips = st.tfood[1] * (1.0f0 - _wwpb_expf(-(alpha * nb1 / (st.tfood[1] + 1.0f-6))))
        else
            st.bkpips = 0.0f0
        end
    end
    return st
end

# =============================================================================
# MULTI-STAND landscape dispersal (bmatct.f full, MXSTND>1) — the interstand
# spatial BKP redistribution. This is the piece that COLLAPSES to bmatct_single!
# when BMSTND=1 & OUTOFF=T (SCORE self-cancels ⇒ PROP=1). Faithful line-by-line
# translation of wwpb/bmatct.f, validated BIT-EXACT vs the gfortran-16 driver
# golden scratchpad/wwpb/driver_bmatct_multi.f (4-stand landscape, OW-off /
# OW-floating / OW-fixed-constant). Transcendentals route through glibc
# logf/expf/powf; SQRT is IEEE-correctly-rounded (= gfortran sqrtf). The spatial
# data model (SPLAAR area, SPLALO location, SPLADS euclidean distance from the
# PPSPLA XLOC/YLOC/AREA arrays) lives in WwpbLandscape.
# =============================================================================

# geometry-unit constants (bmatct.f PARAMETERs, single-precision).
const _WWPB_AC2SQM = Float32(1)/Float32(640)        # acres → sq.mi
const _WWPB_M2MI   = Float32(1)/Float32(1609.34)    # meters → miles
const _WWPB_MI2M   = 1609.34f0                       # miles → meters
const _WWPB_PIE    = 3.14159f0                        # BMPRM PIE

@inline _wwpb_sqrtf(x::Float32)::Float32 = sqrt(x)   # IEEE correctly-rounded == gfortran REAL SQRT

# -----------------------------------------------------------------------------
# WwpbLandscape — the PPSPLA spatial state (per-stand XLOC/YLOC/AREA) plus the
# ACDONE-cached spatial results (ATTC/ATTBYK/AREAO/AREALP/AREAL/RADLC) and the
# per-stand dispersal-output buffers (BKPOUT/BKPIN/SELFBKP/BKPS). One instance
# per landscape; `bmatct_multi!` mutates it across years. The C-constants
# (CBKPO…CRVOND) are resolved-once landscape averages, held here (BMCOM).
# -----------------------------------------------------------------------------
mutable struct WwpbLandscape
    n     ::Int
    xloc  ::Vector{Float32}   # SPLALO X (meters, PPSPLA)
    yloc  ::Vector{Float32}   # SPLALO Y (meters)
    area  ::Vector{Float32}   # SPLAAR area (acres)
    stock ::Vector{Bool}      # STOCK(i)
    # ACDONE-cached spatial results
    acdone::Bool
    attc  ::Matrix{Float32}   # ATTC(2,n)
    attbyk::Matrix{Float32}   # ATTBYK(2,n)
    areao ::Vector{Float32}   # AREAO(2)
    arealp::Vector{Float32}   # AREALP(2)
    areal ::Float32           # AREAL
    radlc ::Float32           # RADLC
    # resolved Outside-World constants (BMCOM CBKPO…CRVOND, resolved once if <0)
    cbkpo ::Vector{Float32}   # CBKPO(2)
    cbaho ::Vector{Float32}   # CBAHO(2)
    cbaspo::Vector{Float32}   # CBASPO(2)
    cspo  ::Vector{Float32}   # CSPO(2)
    cbao  ::Float32           # CBAO
    crvond::Float32           # CRVOND
    # per-stand dispersal-output buffers (BMPCOM)
    bkpout ::Matrix{Float32}  # BKPOUT(2,n)
    bkpin  ::Matrix{Float32}  # BKPIN(2,n)
    selfbkp::Matrix{Float32}  # SELFBKP(2,n)
    bkps   ::Vector{Float32}  # BKPS(n)
end

function WwpbLandscape(xloc::Vector{Float32}, yloc::Vector{Float32},
                       area::Vector{Float32}, stock::Vector{Bool};
                       cbkpo=Float32[-1,-1], cbaho=Float32[-1,-1],
                       cbaspo=Float32[-1,-1], cspo=Float32[-1,-1],
                       cbao::Float32=-1.0f0, crvond::Float32=-1.0f0)
    n = length(xloc)
    z2 = () -> zeros(Float32, 2, n)
    WwpbLandscape(n, copy(xloc), copy(yloc), copy(area), copy(stock),
                  false, z2(), z2(), zeros(Float32, 2), zeros(Float32, 2),
                  0.0f0, 0.0f0,
                  Float32.(cbkpo), Float32.(cbaho), Float32.(cbaspo), Float32.(cspo),
                  Float32(cbao), Float32(crvond),
                  z2(), z2(), z2(), zeros(Float32, n))
end

# -----------------------------------------------------------------------------
# bmatct_multi! (bmatct.f). `stands[i]` supplies NUMER/TFOOD/BKP/BKPIPS (+ BAH/
# BANH/GRFSTD/SPCLT/TREE for the OW landscape-average branch). ATTRACT params:
# `usera/selfa/userc/urmax` are the 3-element USER*(3) arrays (index PBSPEC ∈ 1:3
# for the main beetle, [3] for Ips). Writes the new BKP/BKPIPS back into each
# stand and the BKPOUT/BKPIN/SELFBKP/BKPS buffers into `ls`.
# -----------------------------------------------------------------------------
function bmatct_multi!(ls::WwpbLandscape, stands::Vector{WwpbStand}, w::WwpbState;
                       usera::NTuple{3,Float32}=(1.0f0,1.0f0,1.0f0),
                       selfa::NTuple{3,Float32}=(1.0f0,1.0f0,1.0f0),
                       userc::NTuple{3,Float32}=(100.0f0,100.0f0,100.0f0),
                       urmax::NTuple{3,Float32}=(15.0f0,15.0f0,15.0f0),
                       outoff::Bool=true, ufloat::Float32=-1.0f0, sdd::Float32=0.0f0,
                       rvod::Float32=1.0f0, stocko::Float32=1.0f0,
                       lbad::Bool=false, ibadbb::Int=1,
                       badrep::NTuple{3,Float32}=(1.0f0,1.0f0,1.0f0),
                       ipson::Bool=false,
                       msba::Vector{Float32}=zeros(Float32, WWPB_NSCL))
    n = ls.n
    pbspec = Int(w.pbspec)
    PIE = _WWPB_PIE
    AC2SQM = _WWPB_AC2SQM; M2MI = _WWPB_M2MI; MI2M = _WWPB_MI2M
    npass = ipson ? 2 : 1

    aspec = zeros(Float32, 2); nzero = zeros(Float32, 2)
    ddwt = zeros(Float32, 2);  rmax = zeros(Float32, 2)
    aspec[1] = usera[pbspec]; nzero[1] = selfa[pbspec]
    ddwt[1]  = userc[pbspec]; rmax[1]  = urmax[pbspec]
    if npass == 2
        aspec[2] = usera[3]; nzero[2] = selfa[3]
        ddwt[2]  = userc[3]; rmax[2]  = urmax[3]
    end
    # (GPGET2 option-processor calls omitted: the ATTRACT-schedule keyword changes
    #  are not modeled here — CBKPO… keep their WwpbLandscape values, as with LOK=F.)

    attc = ls.attc; attbyk = ls.attbyk
    kl = zeros(Float32, 2); kls = zeros(Float32, 2); kav = zeros(Float32, 2)
    ko = zeros(Float32, 2); bkpo = zeros(Float32, 2); attwp = zeros(Float32, 2)
    xbkpo = zeros(Float32, 2)

    # ======================= SPATIAL INIT (once, ACDONE) =======================
    if !ls.acdone
        # DO 8: ATTBYK
        @inbounds for curr in 1:n
            caREA = ls.area[curr] * AC2SQM
            r2 = _wwpb_sqrtf(caREA / PIE)
            for ipass in 1:npass
                r2 > rmax[ipass] && (r2 = rmax[ipass])
                attbyk[ipass, curr] = (PIE / ddwt[ipass]) *
                    (_wwpb_logf(abs(ddwt[ipass]*r2*r2 + nzero[ipass])) -
                     _wwpb_logf(abs(nzero[ipass])))
            end
        end

        if outoff
            ls.acdone = true
        else
            xbao = 0.0f0; xrvo = 0.0f0
            xbaho = zeros(Float32, 2); xspo = zeros(Float32, 2); xbaspo = zeros(Float32, 2)
            fill!(xbkpo, 0.0f0)
            # DO 10: landscape centre
            xmax = -1.0f0; ymax = -1.0f0; xmin = 1.0f35; ymin = 1.0f35
            @inbounds for curr in 1:n
                xp = ls.xloc[curr]; yp = ls.yloc[curr]
                xp > xmax && (xmax = xp); xp < xmin && (xmin = xp)
                yp > ymax && (ymax = yp); yp < ymin && (ymin = yp)
            end
            xmid = (xmin + xmax) / 2.0f0; ymid = (ymin + ymax) / 2.0f0
            # DO 13: AREAL/AREALS/DMAX (+ conditional landscape sums)
            areal = 0.0f0; areals = 0.0f0; dmax = -1.0f0
            @inbounds for curr in 1:n
                caREA = ls.area[curr] * AC2SQM
                areal += caREA
                d = _wwpb_sqrtf(caREA / PIE) * MI2M
                xp = ls.xloc[curr]; yp = ls.yloc[curr]
                d = d + _wwpb_sqrtf((xmid-xp)^2 + (ymid-yp)^2)
                d > dmax && (dmax = d)
                ls.stock[curr] || continue
                areals += caREA
                if (ls.cbao < 0.0f0) && (ufloat != -1.0f0)
                    xbao += caREA * (stands[curr].bah[WWPB_NSCL+1] + stands[curr].banh[WWPB_NSCL+1])
                end
                if (ls.crvond < 0.0f0) && (ufloat != -1.0f0)
                    xrvo += caREA * stands[curr].grfstd
                end
                for ipass in 1:npass
                    if (ls.cbkpo[ipass] < 0.0f0) && (ufloat != -1.0f0)
                        if (pbspec != 3) && (ipass == 1)
                            xbkpo[ipass] += caREA * stands[curr].bkp
                        else
                            xbkpo[ipass] += caREA * stands[curr].bkpips
                        end
                    end
                    # host/special sums (only when C-values missing & not fixed OW)
                    if ((ls.cbaho[ipass] >= 0.0f0) && (ls.cspo[ipass] >= 0.0f0) &&
                        (ls.cbaspo[ipass] >= 0.0f0)) || (ufloat == -1.0f0)
                        # GOTO 12
                    else
                        mn = ipass == 1 ? Int(w.iscmin[pbspec]) : Int(w.iscmin[3])
                        tbaho = 0.0f0; tbaspo = 0.0f0; tspo = 0.0f0
                        for isiz in 1:WWPB_NSCL
                            isiz >= mn && (tbaho += stands[curr].bah[isiz])
                            xx = stands[curr].tree[isiz, 1] * stands[curr].spclt[isiz, ipass]
                            tspo += xx
                            isiz >= mn && (tbaspo += xx * msba[isiz])
                        end
                        xbaho[ipass] += caREA * tbaho
                        xspo[ipass]  += caREA * tspo
                        xbaspo[ipass] += caREA * tbaspo
                    end
                end
            end
            radl = dmax * M2MI
            # DO 20: resolve landscape averages (skip if UFLOAT==-1)
            if ufloat != -1.0f0
                ls.cbao   < 0.0f0 && (ls.cbao = xbao / areals)
                ls.crvond < 0.0f0 && (ls.crvond = (xrvo / areals) / rvod)
                for ipass in 1:npass
                    ls.cbkpo[ipass]  < 0.0f0 && (ls.cbkpo[ipass]  = xbkpo[ipass] / areals)
                    ls.cbaho[ipass]  < 0.0f0 && (ls.cbaho[ipass]  = xbaho[ipass] / areals)
                    ls.cbaspo[ipass] < 0.0f0 && (ls.cbaspo[ipass] = xbaspo[ipass] / areals)
                    ls.cspo[ipass]   < 0.0f0 && (ls.cspo[ipass]   = xspo[ipass] / areals)
                end
            end
            # eqns 5,11,15,16
            radlc = _wwpb_sqrtf(areal / PIE)
            for ipass in 1:npass
                rado = radl + rmax[ipass]
                ls.areao[ipass] = (PIE * rado * rado) - areal
                if radlc > rmax[ipass]
                    ls.arealp[ipass] = areal - (PIE * (radlc - rmax[ipass])^2)
                else
                    ls.arealp[ipass] = areal
                end
            end
            ls.areal = areal; ls.radlc = radlc
            # TX (DO 22)
            tx = zeros(Float32, 2)
            for ipass in 1:npass
                tx[ipass] = (PIE / ddwt[ipass]) *
                    (_wwpb_logf(abs(ddwt[ipass]*rmax[ipass]*rmax[ipass] + nzero[ipass])) -
                     _wwpb_logf(abs(nzero[ipass])))
            end
            # ATTC (DO 30)
            @inbounds for curr in 1:n
                ls.stock[curr] || continue
                caREA = ls.area[curr]
                cdist = 0.5f0 * _wwpb_sqrtf(caREA * AC2SQM)
                sx = zeros(Float32, 2)
                for targ in 1:n
                    if targ == curr
                        for ipass in 1:npass
                            sx[ipass] += attbyk[ipass, curr]
                        end
                    else
                        dx0 = ls.xloc[targ] - ls.xloc[curr]
                        dy0 = ls.yloc[targ] - ls.yloc[curr]
                        dist = _wwpb_sqrtf(dx0*dx0 + dy0*dy0) * M2MI
                        tarea = ls.area[targ] * AC2SQM
                        tdist = 0.5f0 * _wwpb_sqrtf(tarea)
                        dist < (cdist + tdist) && (dist = cdist + tdist)
                        for ipass in 1:npass
                            if dist <= rmax[ipass]
                                dxq = ddwt[ipass] * dist * dist
                                sx[ipass] += (tarea / (dxq + nzero[ipass]))
                            end
                        end
                    end
                end
                for ipass in 1:npass
                    attc[ipass, curr] = tx[ipass] - sx[ipass]
                    attc[ipass, curr] < 0.0f0 && (attc[ipass, curr] = 0.0f0)
                end
            end
            ls.acdone = true
        end
    end

    # ============= OUTSIDE-WORLD YEAR SETUP (KL/KO/BKPO/KAV/ATTWP) =============
    if !outoff
        for ipass in 1:npass
            kl[ipass] = 0.0f0; kls[ipass] = 0.0f0; xbkpo[ipass] = 0.0f0
        end
        @inbounds for curr in 1:n
            caREA = ls.area[curr] * AC2SQM
            for ipass in 1:npass
                xx = caREA * stands[curr].numer[ipass]
                kl[ipass] += xx
                if (ufloat == -1.0f0) && ls.stock[curr]
                    kls[ipass] += xx
                    if (pbspec != 3) && (ipass == 1)
                        xbkpo[ipass] += caREA * stands[curr].bkp
                    else
                        xbkpo[ipass] += caREA * stands[curr].bkpips
                    end
                end
            end
        end
        # note: AREALS is only available inside the spatial-init; recompute here.
        areals_now = 0.0f0
        @inbounds for curr in 1:n
            ls.stock[curr] && (areals_now += ls.area[curr] * AC2SQM)
        end
        for ipass in 1:npass
            kl[ipass] = kl[ipass] / ls.areal
        end
        if ufloat == -1.0f0
            for ipass in 1:npass
                ko[ipass]   = kls[ipass] / areals_now
                bkpo[ipass] = xbkpo[ipass] / areals_now
            end
        else
            bao = ls.cbao; rvond = ls.crvond; rvo = rvond * rvod
            for ipass in 1:npass
                baho  = ls.cbaho[ipass]; baspo = ls.cbaspo[ipass]
                spo   = ls.cspo[ipass];  bkpo[ipass] = ls.cbkpo[ipass]
                ko[ipass] = (aspec[ipass]*spo + 1.0f0) * bao * (baho + baspo) / rvo
            end
            if lbad
                if ibadbb != 3 || (ibadbb == 3 && npass == 1)
                    bkpo[1] = bkpo[1] * badrep[pbspec]
                end
                if ibadbb >= 3 && npass == 2
                    bkpo[2] = bkpo[2] * badrep[3]
                end
            end
        end
        if stocko < 1.0f0
            for ipass in 1:npass
                ko[ipass]   = ko[ipass] * stocko
                bkpo[ipass] = bkpo[ipass] * stocko
            end
        end
        for ipass in 1:npass
            kav[ipass] = ((ls.areao[ipass] * ko[ipass]) + (ls.arealp[ipass] * kl[ipass])) /
                         (ls.areao[ipass] + ls.arealp[ipass])
            attwp[ipass] = (PIE / ddwt[ipass]) * kav[ipass] *
                (_wwpb_logf(abs(ddwt[ipass]*rmax[ipass]*rmax[ipass] + nzero[ipass])) -
                 _wwpb_logf(abs(nzero[ipass])))
        end
    end

    # ========================= DISPERSAL (per source) =========================
    newbkp = zeros(Float32, 2, n)
    score  = zeros(Float32, 2, n)
    # zero the output buffers on stockable stands
    @inbounds for curr in 1:n
        ls.stock[curr] || continue
        for ipass in 1:npass
            newbkp[ipass, curr] = 0.0f0
            ls.bkpout[ipass, curr] = 0.0f0
            ls.bkpin[ipass, curr] = 0.0f0
            ls.selfbkp[ipass, curr] = 0.0f0
        end
    end

    totsc = zeros(Float32, 2); attoj = zeros(Float32, 2); totbkp = zeros(Float32, 2)
    @inbounds for curr in 1:n
        ls.stock[curr] || continue
        caREA = ls.area[curr]
        for ipass in 1:npass
            totsc[ipass] = 0.0f0
            outoff || (attoj[ipass] = ko[ipass] * attc[ipass, curr])
            if (pbspec != 3) && (ipass == 1)
                totbkp[ipass] = caREA * stands[curr].bkp
            else
                totbkp[ipass] = caREA * stands[curr].bkpips
            end
        end
        caREA = caREA * AC2SQM
        cdist = 0.5f0 * _wwpb_sqrtf(caREA)
        # target attractiveness
        for targ in 1:n
            if targ == curr
                for ipass in 1:npass
                    score[ipass, targ] = stands[curr].numer[ipass] * attbyk[ipass, curr]
                    totsc[ipass] += score[ipass, targ]
                end
            else
                dx0 = ls.xloc[targ] - ls.xloc[curr]
                dy0 = ls.yloc[targ] - ls.yloc[curr]
                dist = _wwpb_sqrtf(dx0*dx0 + dy0*dy0) * M2MI
                tarea = ls.area[targ] * AC2SQM
                tdist = 0.5f0 * _wwpb_sqrtf(tarea)
                dist < (cdist + tdist) && (dist = cdist + tdist)
                for ipass in 1:npass
                    if dist <= rmax[ipass]
                        dxq = ddwt[ipass] * dist * dist
                        score[ipass, targ] = tarea * stands[targ].numer[ipass] / (dxq + nzero[ipass])
                        totsc[ipass] += score[ipass, targ]
                    else
                        score[ipass, targ] = 0.0f0
                    end
                end
            end
        end
        # Outside-World flux
        if !outoff
            for ipass in 1:npass
                totsc[ipass] += attoj[ipass]
                if totsc[ipass] > 0.0f0
                    ls.bkpout[ipass, curr] = totbkp[ipass] * (attoj[ipass] / totsc[ipass])
                end
                xx = ko[ipass] * ls.areao[ipass] * attwp[ipass]
                rattjo = xx > 0.0f0 ?
                    (caREA * stands[curr].numer[ipass] * attoj[ipass]) / xx : 0.0f0
                areaox = 2753.04f0 * _wwpb_powf(1.0f0 - _wwpb_expf(-0.06366f0 * ls.radlc), 1.5938f0) *
                         _wwpb_powf(1.0f0 - _wwpb_expf(-0.0687f0 * rmax[ipass]), 0.8827f0)
                ls.bkpin[ipass, curr] = bkpo[ipass] * (areaox / AC2SQM) * rattjo
                newbkp[ipass, curr] += ls.bkpin[ipass, curr]
            end
        end
        # allocate BKP to targets
        for targ in 1:n
            for ipass in 1:npass
                prop = 0.0f0
                totsc[ipass] > 0.0f0 && (prop = score[ipass, targ] / totsc[ipass])
                targ == curr && (ls.selfbkp[ipass, targ] = prop * totbkp[ipass])
                newbkp[ipass, targ] += prop * totbkp[ipass]
            end
        end
    end

    # ========================== SATURATION → BKP =============================
    alpha = sdd < 0.0f0 ? (0.365f0 - 0.178f0*sdd) : (0.365f0 - 0.034f0*sdd)
    alpha > 0.90f0 && (alpha = 0.90f0)
    alpha < 0.01f0 && (alpha = 0.01f0)
    @inbounds for curr in 1:n
        ls.stock[curr] || continue
        caREA = ls.area[curr]
        for ipass in 1:2
            newbkp[ipass, curr] = newbkp[ipass, curr] / caREA
            ls.bkpout[ipass, curr] = ls.bkpout[ipass, curr] / caREA
            ls.bkpin[ipass, curr] = ls.bkpin[ipass, curr] / caREA
            ls.selfbkp[ipass, curr] = ls.selfbkp[ipass, curr] / caREA
        end
        s = stands[curr]
        if pbspec != 3
            if (newbkp[1, curr] > 0.0f0) && (s.tfood[1] > 0.0f0)
                s.bkp = s.tfood[1] * (1.0f0 - _wwpb_expf(-(alpha * newbkp[1, curr] / (s.tfood[1] + 1.0f-6))))
                ls.bkps[curr] = s.bkp / newbkp[1, curr]
            else
                s.bkp = 0.0f0; ls.bkps[curr] = 1.0f0
            end
            if newbkp[2, curr] > 0.0f0
                s.bkpips = s.tfood[2] * (1.0f0 - _wwpb_expf(-(alpha * newbkp[2, curr] / (s.tfood[2] + 1.0f-6))))
            else
                s.bkpips = 0.0f0
            end
        else
            if newbkp[1, curr] > 0.0f0
                s.bkpips = s.tfood[1] * (1.0f0 - _wwpb_expf(-(alpha * newbkp[1, curr] / (s.tfood[1] + 1.0f-6))))
                ls.bkps[curr] = s.bkpips / newbkp[1, curr]
            else
                s.bkpips = 0.0f0; ls.bkps[curr] = 1.0f0
            end
        end
    end
    return ls
end

# glibc single-precision logf for AS245/AS63 (matches gfortran REAL LOG bit-exact).
@inline _wwpb_logf(x::Float32)::Float32 = ccall((:logf, "libm.so.6"), Float32, (Float32,), x)

# ALNGAM (bmcbet.f, ALGORITHM AS245 APPL. STATIST. 1989) — log-gamma via rational
# approximations, single-precision REAL. Faithful port (the AS functions are
# defined IN bmcbet.f, so this IS bit-exact-validatable). Coefficients verbatim.
const _WWPB_ALNGAM_R1 = Float32[-2.6668551f0,-2.4438753f1,-2.1969895f1,1.1166754f1,3.1306054f0,6.0777138f-1,1.1940090f1,3.1469011f1,1.5234687f1]
const _WWPB_ALNGAM_R2 = Float32[-7.8335929f1,-1.4204629f2,1.3751941f2,7.8699492f1,4.1643892f0,4.7066876f1,3.1339921f2,2.6350507f2,4.3340002f1]
const _WWPB_ALNGAM_R3 = Float32[-2.1215957f5,2.3066151f5,2.7464764f4,-4.0262111f4,-2.2966072f3,-1.1632849f5,-1.4602593f5,-2.4235740f4,-5.7069100f2]
const _WWPB_ALNGAM_R4 = Float32[2.7919531791f-1,4.9173176105f-1,6.9291059929f-2,3.3503438150f0,6.0124592597f0]
const _WWPB_ALR2PI = 9.1893853320f-1
function _wwpb_alngam(xvalue::Float32)::Float32
    r1 = _WWPB_ALNGAM_R1; r2 = _WWPB_ALNGAM_R2; r3 = _WWPB_ALNGAM_R3; r4 = _WWPB_ALNGAM_R4
    x = xvalue
    (x >= 1.0f38 || x <= 0.0f0) && return 0.0f0          # ifault 2 / 1
    if x < 1.5f0
        if x < 0.5f0
            alngam = -_wwpb_logf(x); y = x + 1.0f0
            y == 1.0f0 && return alngam
        else
            alngam = 0.0f0; y = x; x = (x - 0.5f0) - 0.5f0
        end
        return alngam + x * ((((r1[5]*y+r1[4])*y+r1[3])*y+r1[2])*y+r1[1]) /
                            ((((y+r1[9])*y+r1[8])*y+r1[7])*y+r1[6])
    end
    if x < 4.0f0
        y = (x - 1.0f0) - 1.0f0
        return y * ((((r2[5]*x+r2[4])*x+r2[3])*x+r2[2])*x+r2[1]) /
                   ((((x+r2[9])*x+r2[8])*x+r2[7])*x+r2[6])
    end
    if x < 12.0f0
        return ((((r3[5]*x+r3[4])*x+r3[3])*x+r3[2])*x+r3[1]) /
               ((((x+r3[9])*x+r3[8])*x+r3[7])*x+r3[6])
    end
    y = _wwpb_logf(x)
    alngam = x * (y - 1.0f0) - 0.5f0 * y + _WWPB_ALR2PI
    x > 5.10f6 && return alngam
    x1 = 1.0f0 / x; x2 = x1 * x1
    return alngam + x1 * ((r4[3]*x2+r4[2])*x2+r4[1]) / ((x2+r4[5])*x2+r4[4])
end

# -----------------------------------------------------------------------------
# AS63 regularized incomplete beta I_x(p,q) (Float32). RECONSTRUCT-category: the
# WWPB ALNGAM/BETAIN are absent from the tree (only declared in bmcbet.f). Ported
# as the canonical Applied-Statistics AS63 algorithm, IDENTICAL arithmetic to the
# Fortran driver stub (scratchpad/wwpb/driver_bmcbet.f) so BMCBET's own logic
# validates bit-exact; the AS function itself is the shared reconstruction. `beta`
# = log B(p,q) = lgammaf(p)+lgammaf(q)−lgammaf(p+q).
# -----------------------------------------------------------------------------
function _wwpb_betain(x::Float32, p::Float32, q::Float32, beta::Float32)::Float32
    acu = 0.1f-14
    (x <= 0.0f0) && return 0.0f0
    (x >= 1.0f0) && return 1.0f0
    psq = p + q
    cx = 1.0f0 - x
    if p < psq * x
        xx = cx; cx = x; pp = q; qq = p; indx = true
    else
        xx = x; pp = p; qq = q; indx = false
    end
    term = 1.0f0; ai = 1.0f0; betain = 1.0f0
    ns = Int(trunc(qq + cx * psq))
    rx = xx / cx
    temp = qq - ai
    ns == 0 && (rx = xx)
    while true
        term = term * temp * rx / (pp + ai)
        betain = betain + term
        temp = abs(term)
        if temp <= acu && temp <= acu * betain
            break
        end
        ai += 1.0f0; ns -= 1
        if ns >= 0
            temp = qq - ai
            ns == 0 && (rx = xx)
        else
            temp = psq; psq += 1.0f0
        end
    end
    betain = betain * _wwpb_expf(pp * _wwpb_logf(xx) + (qq - 1.0f0) * _wwpb_logf(cx) - beta) / pp
    return indx ? (1.0f0 - betain) : betain
end

# -----------------------------------------------------------------------------
# bmcbet! (bmcbet.f) — beta-distribution weights BETA(NSCL) over size classes,
# keyed on ABETA (a=ABETA, b=2). BETA(i) = I_{x_i}(a,b) − I_{x_{i-1}}(a,b) over
# the mapped interval [MINSIZE−0.5, MAXSIZE+0.5]; classes init 1e-5. The BMCBET
# wrapper logic is faithful to the pristine routine (RECONSTRUCT only for the
# absent ALNGAM/BETAIN, shared with the Fortran driver ⇒ bit-exact validation).
# -----------------------------------------------------------------------------
function bmcbet!(abeta::Float32, minsize::Int, maxsize::Int)::Vector{Float32}
    beta = fill(1.0f-5, WWPB_NSCL)
    a = abeta; b = 2.0f0
    minx = Float32(minsize) - 0.5f0
    maxx = Float32(maxsize) + 0.5f0
    logbeta = _wwpb_alngam(a) + _wwpb_alngam(b) - _wwpb_alngam(a + b)
    temp0 = 0.0f0
    @inbounds for i in minsize:maxsize
        x = ((Float32(i) + 0.5f0) - minx) / (maxx - minx)
        x <= 0.0f0 && (x = 0.001f0)
        x == 1.0f0 && (x = 1.0f0 - 0.001f0)
        (x < 0.0f0 || x > 1.0f0) && return beta
        (x == 0.0f0 || x == 1.0f0) && return beta
        temp1 = _wwpb_betain(x, a, b, logbeta)
        beta[i] = temp1 - temp0
        temp0 = temp1
    end
    return beta
end

# -----------------------------------------------------------------------------
# bmistd! (bmistd.f) — the STOCHASTIC within-stand kill allocation: distribute
# the stand BKP into per-size-class beetle kills (PBKILL), strip-kills (STRIP),
# and pitch-outs (PITCH), using the beta-distribution size preference (bmcbet!)
# and the BMRANN random stream (wwpb_rand!) — the exact call order is load-
# bearing. Phases: SPRAY reduce → MXISIZ/ISIZ1/ABETA → bmcbet! → special-tree
# kills (BMRANN) → deterministic group-kill (DO 555) → individual-kill (BMRANN,
# DO 888) → PBKILL proportion→TPA + standing-dead-wood pool. Faithful to
# bmistd.f. NZERO=1e-6. `sarea` = stand area (acres, from SPLAAR/the FVS stand).
# -----------------------------------------------------------------------------
function bmistd!(st::WwpbStand, w::WwpbState, coeffs; sarea::Float32)
    NZERO = 1.0f-6; NUNIT = 1.0f0 - NZERO
    msba = coeffs.msba; l2d = coeffs.l2d
    p = zeros(Float32, WWPB_NSCL)
    st.bkp = st.bkp * (1.0f0 - st.spray)
    if st.bkp >= NZERO
        attp = 1.0f0 / sarea
        fill!(st.pitch, 0.0f0)
        mxisiz = 1
        @inbounds for isiz in WWPB_NSCL:-1:1
            if st.tree[isiz, 1] > 0.0f0; mxisiz = isiz; break; end
        end
        isiz1 = 0
        @inbounds for i in WWPB_NSCL:-1:1
            if msba[i] * st.grf[i] <= st.bkp; isiz1 = i; break; end
        end
        isiz1 <= 0 && (isiz1 = 1)
        abeta = 1.0f0
        if st.bkp > 6.0f0
            abeta = 15.0f0
        elseif st.bkp <= 6.0f0 && st.bkp > 3.6f0
            abeta = 2.5f0 + 4.46f0 * (st.bkp - 3.6f0)
        elseif st.bkp <= 3.6f0 && st.bkp >= 1.6f0
            abeta = 1.2f0 + 0.65f0 * (st.bkp - 1.6f0)
        end
        st.bkp < 1.6f0 && (abeta = 1.0f0)
        beta = bmcbet!(abeta, 1, isiz1)

        # ---- special-tree kills ----
        tbasp = 0.0f0
        @inbounds for isiz in 1:mxisiz
            if st.spclt[isiz, 1] * st.tree[isiz, 1] > NZERO
                tbasp += st.tree[isiz, 1] * st.spclt[isiz, 1] * msba[isiz]
            end
        end
        if tbasp > NZERO
            x = 0.0f0
            @inbounds for isiz in 1:mxisiz
                if st.spclt[isiz, 1] * st.tree[isiz, 1] > NZERO
                    x += msba[isiz] * st.tree[isiz, 1] * st.spclt[isiz, 1] / tbasp
                end
                x > NUNIT && (x = 1.0f0); p[isiz] = x
            end
            spkill = 0
            acres = max(floor(Int, sarea + 0.5f0), 1)
            while spkill < acres && st.bkp > NZERO && tbasp > NZERO
                xr = wwpb_rand!(w); isiz = mxisiz
                @inbounds for i in 1:mxisiz
                    if xr <= p[i] || i == mxisiz; isiz = i; break; end
                end
                host = st.tree[isiz, 1] * (st.spclt[isiz, 1] - st.pbkill[isiz])
                attprp = attp >= host ? host / st.tree[isiz, 1] : attp / st.tree[isiz, 1]
                bkpuse = attprp * st.tree[isiz, 1] * msba[isiz]
                bkpkl = bkpuse * st.grf[isiz]
                if st.bkp >= bkpkl
                    st.final[1] = st.bkp < bkpuse ? st.bkp : bkpuse
                    st.final[2] = Float32(isiz); st.final[3] = attprp * st.tree[isiz, 1]
                    st.pbkill[isiz] += attprp; spkill += 1
                    st.bkp -= bkpuse; st.bkp < NZERO && (st.bkp = 0.0f0)
                    if st.pbkill[isiz] > st.spclt[isiz, 1]
                        st.pbkill[isiz] = st.spclt[isiz, 1]
                        tbasp = 0.0f0
                        @inbounds for i in 1:mxisiz
                            if st.spclt[i, 1] > st.pbkill[i] && st.tree[i, 1] > NZERO
                                tbasp += st.tree[i, 1] * msba[i] * (st.spclt[i, 1] - st.pbkill[i])
                            end
                        end
                        tbasp <= NZERO && break
                        x = 0.0f0
                        @inbounds for i in 1:mxisiz
                            if st.spclt[i, 1] > st.pbkill[i] && st.tree[i, 1] > NZERO
                                x += msba[i] * st.tree[i, 1] * (st.spclt[i, 1] - st.pbkill[i]) / tbasp
                            end
                            x > NUNIT && (x = 1.0f0); p[i] = x
                        end
                    end
                else
                    xf = st.bkp / bkpkl
                    if xf >= 0.75f0
                        st.strip[isiz] += attprp * st.tree[isiz, 1]; st.pitch[isiz] += attprp
                    else
                        st.pitch[isiz] += attprp
                        st.final[1] = st.bkp * 0.15f0; st.final[2] = Float32(isiz)
                        st.final[3] = attprp * st.tree[isiz, 1]
                    end
                    st.bkp = 0.0f0
                end
            end
        end

        # ---- random group/individual kills (only if BKP left) ----
        if st.bkp >= NZERO
            pscale = zeros(Float32, WWPB_NSCL)
            @inbounds for isiz in 1:mxisiz
                pscale[isiz] = msba[isiz] * st.tree[isiz, 1] * (1.0f0 - st.pbkill[isiz]) / st.grf[isiz]
            end
            scale = 0.0f0; avklba = 0.0f0
            @inbounds for isiz in 1:mxisiz
                scale += beta[isiz] * pscale[isiz]; avklba += beta[isiz] * pscale[isiz] * msba[isiz]
            end
            if scale < NZERO
                st.bkp = 0.0f0
            else
                avklba /= scale
                totkl = st.bkp * sarea / avklba
                # group-kill
                while totkl > Float32(2 * mxisiz)
                    @inbounds for isiz in 1:mxisiz
                        if beta[isiz] * pscale[isiz] > NZERO
                            host = st.tree[isiz, 1] * (1.0f0 - st.pbkill[isiz])
                            sckl = Float32(trunc(totkl * beta[isiz] * pscale[isiz] / scale))
                            if sckl >= host * sarea
                                attprp = host / st.tree[isiz, 1]; bkpuse = host * msba[isiz]
                            else
                                attprp = sckl / (st.tree[isiz, 1] * sarea); bkpuse = (sckl / sarea) * msba[isiz]
                            end
                            st.bkp -= bkpuse; st.pbkill[isiz] += attprp
                            st.pbkill[isiz] > NUNIT && (st.pbkill[isiz] = 1.0f0)
                        end
                    end
                    @inbounds for isiz in 1:mxisiz
                        st.pbkill[isiz] >= 1.0f0 && (pscale[isiz] = 0.0f0)
                    end
                    scale = 0.0f0; avklba = 0.0f0
                    @inbounds for isiz in 1:mxisiz
                        scale += beta[isiz] * pscale[isiz]; avklba += beta[isiz] * pscale[isiz] * msba[isiz]
                    end
                    if scale > NZERO
                        avklba /= scale; totkl = st.bkp * sarea / avklba
                    else
                        st.bkp = 0.0f0; totkl = 0.0f0
                    end
                end
                # individual-kill P
                if scale > NZERO && st.bkp > NZERO
                    x = 0.0f0
                    @inbounds for isiz in 1:mxisiz
                        if beta[isiz] * pscale[isiz] > NZERO
                            x += beta[isiz] * pscale[isiz] / scale
                        end
                        x > NUNIT && (x = 1.0f0); p[isiz] = x
                    end
                end
                while st.bkp > NZERO
                    if scale <= NZERO; st.bkp = 0.0f0; break; end
                    xr = wwpb_rand!(w); isiz = mxisiz
                    @inbounds for i in 1:mxisiz
                        if xr <= p[i] || i == mxisiz; isiz = i; break; end
                    end
                    host = st.tree[isiz, 1] * (1.0f0 - st.pbkill[isiz])
                    if attp >= host
                        attprp = host / st.tree[isiz, 1]; bkpuse = host * msba[isiz]
                    else
                        attprp = attp / st.tree[isiz, 1]; bkpuse = attp * msba[isiz]
                    end
                    bkpkl = bkpuse * st.grf[isiz]
                    if st.bkp >= bkpkl
                        st.final[1] = st.bkp < bkpuse ? st.bkp : bkpuse
                        st.final[2] = Float32(isiz); st.final[3] = attprp * st.tree[isiz, 1]
                        st.bkp -= bkpuse; st.bkp < NZERO && (st.bkp = 0.0f0)
                        st.pbkill[isiz] += attprp
                        if st.pbkill[isiz] > NUNIT
                            st.pbkill[isiz] = 1.0f0; pscale[isiz] = 0.0f0; scale = 0.0f0
                            @inbounds for ii in 1:mxisiz; scale += pscale[ii] * beta[ii]; end
                            if scale > NZERO
                                x = 0.0f0
                                @inbounds for ii in isiz:mxisiz
                                    if beta[ii] * pscale[ii] > NZERO
                                        x += beta[ii] * pscale[ii] / scale
                                    end
                                    x > NUNIT && (x = 1.0f0); p[ii] = x
                                end
                            end
                        end
                    else
                        xf = st.bkp / bkpkl
                        if xf >= 0.75f0
                            st.strip[isiz] += attprp * st.tree[isiz, 1]; st.pitch[isiz] += attprp
                        else
                            st.pitch[isiz] += attprp
                            st.final[1] = 0.15f0 * st.bkp; st.final[2] = Float32(isiz)
                            st.final[3] = attprp * st.tree[isiz, 1]
                        end
                        st.bkp = 0.0f0
                    end
                end
            end
        end
    end
    # convert PBKILL proportion → TPA + standing dead wood
    @inbounds for isiz in 1:WWPB_NSCL
        st.pbkill[isiz] <= 0.0f0 && continue
        st.pbkill[isiz] *= st.tree[isiz, 1]
        j = l2d[isiz] + 1; k = max(Int(st.iqptyp[1]), 1)
        st.sdwp[k, j, 1] += st.pbkill[isiz] * st.tvol[isiz, 1]
    end
    return st
end

# -----------------------------------------------------------------------------
# bmkill! (bmkill.f) — the BM→FVS mortality HANDBACK. For each FVS tree record,
# convert the beetle-model per-size-class kill ledger (TPBK) + initial TPA (OTPA)
# into a per-record mortality and reconcile with the existing FVS mortality WK2:
#   MFAST/MSLOW/MBTL = min(TPBK(K,type,{1,2,3})/OTPA(K,type), 1)·PROB;  MSUM=sum;
#   if MSUM>MPRG(=WK2): if MFAST+MBTL<MPRG → MSLOW−=(MSUM−MPRG) else MPRG=MFAST+
#   MBTL; else the FVS-only excess goes to dead wood (not the tree mortality).
#   WK2 = MPRG, bounded PROB−WK2 ≥ 1e-6. Faithful to bmkill.f (the SVS/sanitation/
#   dead-wood-pool bookkeeping — LOKS/BTKL/SDWP DO-200 — is deferred; it does not
#   change the returned WK2 tree mortality). `wk2` is modified in place.
# -----------------------------------------------------------------------------
function bmkill!(st::WwpbStand, w::WwpbState, t, wk2::Vector{Float32}, sp_alpha::Function)
    pbspec = Int(w.pbspec)
    upsiz = w.upsiz
    @inbounds for i in 1:t.n
        t.tpa[i] > 0.0f0 || continue
        sp = Int(t.species[i])
        k = wwpb_dbh_class(upsiz, t.dbh[i])
        lx = wwpb_is_host(pbspec, sp_alpha(sp))
        ty = lx ? 1 : 2
        x = st.otpa[k, ty] > 1.0f-9 ? (1.0f0 / st.otpa[k, ty]) : 1.0f9
        mfast = min(st.tpbk[k, ty, 1] * x, 1.0f0) * t.tpa[i]
        mslow = min(st.tpbk[k, ty, 2] * x, 1.0f0) * t.tpa[i]
        mbtl  = min(st.tpbk[k, ty, 3] * x, 1.0f0) * t.tpa[i]
        msum = mfast + mslow + mbtl
        mprg = wk2[i]
        if msum > mprg
            if (mfast + mbtl) < mprg
                mslow -= (msum - mprg)
            else
                mprg = mfast + mbtl
            end
        end
        wk2[i] = mprg
        (t.tpa[i] - wk2[i]) < 1.0f-6 && (wk2[i] = t.tpa[i] - 1.0f-6)
    end
    return wk2
end

# -----------------------------------------------------------------------------
# wwpb_outbreak_cycle! — the reconstructed single-stand orchestration (the
# PPMAIN/BMDRV outer loop, USER-approved; source ABSENT so this is faithful
# translation, not bit-exact). Runs one FVS cycle's outbreak: load the FVS
# treelist into the BM size-class table (bmsdit!), then the per-year loop over
# [iyr1, iyr2] of the bit-exact beetle kernels, then leave the accumulated TPBK
# ledger for bmkill!. MINIMAL outbreak (no fire/wind/lightning/management/other-
# agent keywords): the stressor models default neutral — RVDSC=1 (bmdrgt.f:132,
# "drought model not run"), SDD=0, all others 0 — so the per-year chain is
# bmcgrf→bmcbkp→bmcnum→bmatct_single→bmistd→bmmort(slow). The outbreak must be
# SEEDED: BKP starts at 0, so the first year's bmcbkp needs a nonzero PBKILL
# (inventory beetle-damage, bmsdit.f MICYC==2 LBMDAM) or Outside-World BKPIN
# (deferred). `seed_pbkill` supplies that initial per-size-class beetle-killed
# proportion. Deterministic given the BMRANN seed in `w`.
# -----------------------------------------------------------------------------
function wwpb_outbreak_cycle!(st::WwpbStand, w::WwpbState, coeffs, t, sp_alpha::Function;
                              sarea::Float32, iyr1::Int, iyr2::Int, ipson::Bool=false,
                              seed_pbkill::Union{Nothing,Vector{Float32}}=nothing,
                              habtyp::Integer=0, slope::Real=0.0)
    bmsdit!(st, t, w, sp_alpha; habtyp=habtyp, slope=slope)
    fill!(st.rvdsc, 1.0f0)                      # drought model not run ⇒ RVDSC=1 (neutral)
    seed_pbkill !== nothing && (st.pbkill .= seed_pbkill)   # inventory-damage seed
    oldgrf = zeros(Float32, WWPB_NSCL)
    @inbounds for _iyr in iyr1:iyr2
        bmcgrf!(st, w, oldgrf)
        bmcbkp!(st, w, coeffs; ipson=ipson)
        bmcnum!(st, w, coeffs; ipson=ipson)
        bmatct_single!(st, w; sdd=0.0f0, ipson=ipson)
        bmistd!(st, w, coeffs; sarea=sarea)
        bmmort!(st, true)                        # slow pass: PBKILL→TREE decrement + TPBK ledger
    end
    return st
end

# -----------------------------------------------------------------------------
# bmdrv_multi! (bmdrv.f, MXSTND>1) — the interstand master driver loop. Faithful
# to bmdrv.f's two-phase per-year structure, but with the bit-exact multi-stand
# BMATCT hoisted to the landscape level (that is the whole point of mode-2): all
# stands run phase-1 {BMCGRF/BMCBKP/BMCNUM} first, THEN one landscape BMATCT
# redistributes BKP across neighbors, THEN all stands run phase-2 {BMISTD/BMMORT}
# on the redistributed BKP. Minimal-outbreak chain (drought/fire/wind/management
# neutral, as in wwpb_outbreak_cycle!): the omitted BMDRGT/BMLITE/BMFIRE/BMOBB/
# BMDFOL/BMQMRT/BMMORT(fast) default to no-op when their keywords are absent.
# Stands must be pre-loaded (bmsdit!) with rvdsc=1 and any inventory seed. At
# MXSTND=1 / OUTOFF=T this collapses to wwpb_outbreak_cycle! (bmatct_multi!→
# bmatct_single! degeneracy) — verified by test. No BMDRV Fortran oracle exists
# (FVSppe links the exbm.f BMDRV stub); this composes driver-golden-validated
# kernels, and the multi-stand BMATCT is itself golden-validated.
# -----------------------------------------------------------------------------
function bmdrv_multi!(ls::WwpbLandscape, stands::Vector{WwpbStand}, w::WwpbState, coeffs;
                      area::Vector{Float32}, iyr1::Int, iyr2::Int, ipson::Bool=false,
                      usera::NTuple{3,Float32}=(1.0f0,1.0f0,1.0f0),
                      selfa::NTuple{3,Float32}=(1.0f0,1.0f0,1.0f0),
                      userc::NTuple{3,Float32}=(100.0f0,100.0f0,100.0f0),
                      urmax::NTuple{3,Float32}=(15.0f0,15.0f0,15.0f0),
                      outoff::Bool=true, ufloat::Float32=-1.0f0, sdd::Float32=0.0f0,
                      rvod::Float32=1.0f0, stocko::Float32=1.0f0,
                      lbad::Bool=false, ibadbb::Int=1,
                      badrep::NTuple{3,Float32}=(1.0f0,1.0f0,1.0f0))
    oldgrf = zeros(Float32, WWPB_NSCL)
    @inbounds for _iyr in iyr1:iyr2
        # phase 1: per-stand GRF / BKP / attractiveness numerator (bmdrv.f DO 40)
        for i in eachindex(stands)
            ls.stock[i] || continue
            bmcgrf!(stands[i], w, oldgrf)
            bmcbkp!(stands[i], w, coeffs; ipson=ipson)
            bmcnum!(stands[i], w, coeffs; ipson=ipson, usera=usera)
        end
        # landscape BKP redistribution (bmdrv.f: CALL BMATCT(IYR))
        bmatct_multi!(ls, stands, w; usera=usera, selfa=selfa, userc=userc, urmax=urmax,
                      outoff=outoff, ufloat=ufloat, sdd=sdd, rvod=rvod, stocko=stocko,
                      lbad=lbad, ibadbb=ibadbb, badrep=badrep, ipson=ipson, msba=coeffs.msba)
        # phase 2: per-stand within-stand dynamics + beetle-kill (bmdrv.f DO 45)
        for i in eachindex(stands)
            ls.stock[i] || continue
            bmistd!(stands[i], w, coeffs; sarea=area[i])
            bmmort!(stands[i], true)
        end
    end
    return ls
end

# -----------------------------------------------------------------------------
# wwpb_main_report (bmout.f MAINOUT DO-20 + IPS-slash loop) — aggregate the
# end-of-cycle WwpbStand into the FVS_BM_Main "MAINOUT" stand-summary variables
# that DBSBMMAIN serializes (dbs/dbsbmmain.f). Returns a NamedTuple keyed by the
# FVS_BM_Main column names. Deterministic given the stand state; the arithmetic
# is BIT-EXACT vs gfortran-16 bmout.f (scratchpad/wwpb/bmdbs/golden_bmmain.f) on
# a controlled state. Sanitation columns (BA_San_Remv … VolRemSalv) are 0 unless
# a WWPB sanitation-harvest keyword ran (MXHRVP path, not yet ported) — matching
# the oracle, which also emits 0 there absent a sanitation cut. The outbreak
# STATE feeding this is jl's reconstructed single-stand harness (cornered, like
# every wwpb_apply! kill), so run_keyfile FVS_BM_Main rows are cornered; the
# aggregation + serialization themselves are transcription-golden.
# -----------------------------------------------------------------------------
function wwpb_main_report(st::WwpbStand, coeffs)
    msba = coeffs.msba
    bak_yr = 0.0f0; tpa_yr = 0.0f0; tpak_yr = 0.0f0
    vol_yr = 0.0f0; volh_yr = 0.0f0; volk_yr = 0.0f0
    ba_sp = 0.0f0; spcl_tpa = 0.0f0
    # per-size-class vectors (bmout.f TREEOUT/VOLOUT), consumed by FVS_BM_Tree/FVS_BM_Vol
    tpa_sc = zeros(Float32, WWPB_NSCL); host_sc = zeros(Float32, WWPB_NSCL)
    tkld_sc = zeros(Float32, WWPB_NSCL); spcl_sc = zeros(Float32, WWPB_NSCL)
    tvol_sc = zeros(Float32, WWPB_NSCL); hvol_sc = zeros(Float32, WWPB_NSCL)
    volk_sc = zeros(Float32, WWPB_NSCL)
    @inbounds for i in 1:WWPB_NSCL
        tpakll = st.pbkill[i] + st.allkll[i]
        prophkld = st.tree[i, 1] > 1.0f-6 ? tpakll / st.tree[i, 1] : 0.0f0
        tpa_sc[i]  = st.tree[i, 1] + st.tree[i, 2]
        host_sc[i] = st.tree[i, 1]
        tkld_sc[i] = tpakll
        spcl_sc[i] = st.spclt[i, 1] * st.tree[i, 1]
        tvol_sc[i] = st.tvol[i, 1] * st.tree[i, 1] + st.tvol[i, 2] * st.tree[i, 2]
        hvol_sc[i] = st.tvol[i, 1] * st.tree[i, 1]
        volk_sc[i] = tpakll * st.tvol[i, 1]
        bak_yr  += st.bah[i] * prophkld
        tpa_yr  += tpa_sc[i]
        tpak_yr += tpakll
        vol_yr  += tvol_sc[i]
        volh_yr += hvol_sc[i]
        volk_yr += volk_sc[i]
        ba_sp   += st.tree[i, 1] * msba[i] * st.spclt[i, 1]
        spcl_tpa += spcl_sc[i]
    end
    ips_slsh = 0.0f0
    @inbounds for j in 1:size(st.dwphos, 1), k in 1:size(st.dwphos, 2)
        ips_slsh += st.dwphos[j, k]
    end
    return (PreDispBKP = st.oldbkp, PostDispBKP = st.bkp, StandRV = st.grfstd,
            StandBA = st.bastd, BAH = st.bah[WWPB_NSCL+1], BA_BtlKld = bak_yr,
            TPA = tpa_yr, TPAH = st.tree[WWPB_NSCL+1, 1], TPA_BtlKld = tpak_yr,
            StandVol = vol_yr, VolHost = volh_yr, VolBtlKld = volk_yr,
            BA_Special = ba_sp, Ips_Slash = ips_slsh, SpclTPA = spcl_tpa,
            BA_San_Remv = 0.0f0, BKP_San_Remv = 0.0f0, TPA_SanRemvLv = 0.0f0,
            TPA_SanRemLvDd = 0.0f0, VolRemSan = 0.0f0, VolRemSalv = 0.0f0,
            tpa_sc = tpa_sc, host_sc = host_sc, tkld_sc = tkld_sc, spcl_sc = spcl_sc,
            tvol_sc = tvol_sc, hvol_sc = hvol_sc, volk_sc = volk_sc,
            # FVS_BM_BKP (bmout.f BKPOUT pass-through): kernel-output BKP state. The
            # landscape-dispersal terms (SELFBKP/TO_LS/FRM_LS/IN_OW/OUT2OW/PER_SURV) come
            # from the source-absent PPMAIN orchestrator ⇒ 0 in the single-stand harness.
            bkp_strp = st.final[1], bkp_strp_sc = st.final[2],
            bkp_dvrv = copy(st.dvrv), bkp_fastk = copy(st.fastk))
end

# -----------------------------------------------------------------------------
# wwpb_apply! — the simulate.jl seam. If a DISPERSE outbreak is active, run one
# cycle's single-stand outbreak on the cycle-start stand and reconcile the beetle
# kills with the FVS mortality already applied (mpb_apply!-style), then write the
# combined mortality back to the treelist. Inert unless w.outbreak (a DISPERSE
# keyword was read). The composed outbreak is faithful-reconstruction (the outer
# harness source is absent); every kernel it calls is bit-exact. Called pre-
# tripling from simulate.jl, like the other insect seams.
# -----------------------------------------------------------------------------
function wwpb_apply!(s, old_tpa::Vector{Float32}, fint::Real)
    w = s.wwpb::WwpbState
    (w.outbreak && Int(w.iyr1) > 0) || return s
    t = s.trees
    t.n > 0 || return s
    coeffs = wwpb_init_coeffs(w.upsiz)
    st = WwpbStand()
    code = s.coef.code_alpha
    spα(sp::Int) = (1 <= sp <= length(code)) ? code[sp] : ""
    sarea = s.plot.gross_space > 0.0f0 ? s.plot.gross_space : 1.0f0
    seed = zeros(Float32, WWPB_NSCL)
    sc = Int(w.seed_class)
    (1 <= sc <= WWPB_NSCL) && (seed[sc] = w.seed_tpa)
    # the beetle model works on the cycle-start stand; save the FVS mortality
    # already applied, restore cycle-start TPA for the outbreak + bmkill PROB.
    fvs_mort = Float32[old_tpa[i] - t.tpa[i] for i in 1:t.n]
    @inbounds for i in 1:t.n; t.tpa[i] = old_tpa[i]; end
    wwpb_outbreak_cycle!(st, w, coeffs, t, spα; sarea = sarea,
                         iyr1 = Int(w.iyr1), iyr2 = Int(w.iyr2), seed_pbkill = seed)
    # PPBMMAIN/PPBMTREE/PPBMVOL: accumulate the per-cycle report (holds the stand summary + per-class
    # vectors) for FVS_BM_Main/Tree/Vol (written at finalize).
    if s.control.dbs_bm_main || s.control.dbs_bm_tree || s.control.dbs_bm_vol || s.control.dbs_bm_bkp
        yr = current_cycle_year(s)
        push!(w.main_rows, (yr, wwpb_main_report(st, coeffs)))
    end
    wk2 = copy(fvs_mort)
    bmkill!(st, w, t, wk2, spα)                     # PROB = t.tpa (= old_tpa now)
    @inbounds for i in 1:t.n
        t.tpa[i] = max(old_tpa[i] - wk2[i], 0.0f0)  # reconciled FVS+beetle mortality
    end
    return s
end
