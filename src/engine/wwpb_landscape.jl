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
const WWPB_PI24 = Float32(3.14159 / (24.0 * 24.0))   # bmsdit.f PI24 = PIE/(24·24), the BA constant

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
