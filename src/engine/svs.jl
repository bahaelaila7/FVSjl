# =============================================================================
# svs.jl — SVS (Stand Visualization System) data path — CHUNK 0
#
# Ported from:
#   base/svkey.f    (SVKEY)    — the SVS keyword flag + _index.svs open
#   vbase/svstart.f (SVSTART)  — the cycle-0 inventory-picture seam (fvs.f:333)
#   base/svgtpl.f   (SVGTPL)   — deterministic subplot layout (IPLGEM=0 branch)
#   base/svestb.f   (SVESTB)   — TPA→objects placement (integer path)
#   base/svgtpt.f   (SVGTPT)   — random point in a rectangle (IPLGEM<2 branch)
#   base/svobol.f   (SVOBOL)   + base/svcrol.f (SVCROL) — circle/circle overlap
#   base/svrann.f   (SVRANN)   — the 2nd Park–Miller stream (divisor 2^31) — in rng.jl
#   base/cwcalc.f   (CWCALC)   — KT western forest-grown crown width (Crookston R1/R6)
#   vbase/svout.f   (SVOUT)    — header (:243-258) + live-tree object record (:429-440)
#
# SCOPE (chunk 0): cycle-0 inventory picture only; IPLGEM=0 (one 208.71-ft square acre,
# subplot ids ignored); INTEGER TPA (FRACAD≡0 ⇒ the fractional lottery is skipped, so the
# ONLY svrann draws are SVGTPT's two-per-object); live green trees (IOBJTP=1). Western
# #TREEFORM = WEST.TRF. Later chunks add: fractional lottery + overlap retries (1),
# multi-subplot layout (2), later cycles (3), cut/mortality removal (4), snags/CWD (5).
#
# Validated bit-exact vs the live relinked FVSkt on stand S248112 (kt0.key/kt0.tre):
# _001.svs is byte-identical (see test/integration/test_svs_chunk0.jl).
# =============================================================================

const _SVS_SIDE_IMPERIAL = 208.7103f0   # sqrt(43560) — square-acre side, feet (svgtpl.f:46)
const _SVS_DIAI_FACTOR   = 0.04166667f0 # DBH(in) → stem radius(ft): 1/12 × 0.5 (svestb.f:308)

# KT crown-width equation code per KT species sequence number (cwcalc.f KTMAP, ISPC 1..11):
#   1 WP 2 WL 3 DF 4 GF 5 WH 6 RC 7 LP 8 ES 9 AF 10 PP 11 OT/MH
const _KT_CWEQN = ("11903","07303","20203","01703","26303","24203",
                   "10803","09303","01903","12203","26405")

"""
    kt_crown_width(ispc, D, H, CR, BA, EL) -> Float32

Forest-grown crown width (ft) for KT species sequence `ispc`, DBH `D`, height `H`, crown
ratio percent `CR`, stand basal area `BA`, elevation `EL` (100s ft). Mirrors CWIDTH→CWCALC
(cwcalc.f) for the Kootenai variant: KT is FS region 1 (KODFOR<601) ⇒ the R6 forest factor
BF=1 and the R6-forest adjustment block is skipped. `CL = CR·H·0.01` is crown length;
`BAREA = max(BA, 1)`. All arithmetic in Float32 to match the Fortran REAL single precision.
"""
function kt_crown_width(ispc::Integer, D::Float32, H::Float32, CR::Float32,
                        BA::Float32, EL::Float32)::Float32
    cl    = CR * H * 0.01f0
    barea = BA <= 1f0 ? 1f0 : BA
    eq    = _KT_CWEQN[ispc]
    cw    = 0f0
    if eq == "11903"           # WP — Crookston (R1)
        b = D >= 1f0 ? D : 1f0
        cw = 1.0405f0 * exp(1.2799f0 + 0.11941f0*log(cl) + 0.42745f0*log(b) - 0.07182f0*log(barea))
        D < 1f0 && (cw *= D / 1f0); cw > 35f0 && (cw = 35f0)
    elseif eq == "07303"       # WL — Crookston (R1)
        b = D >= 1f0 ? D : 1f0
        cw = 1.02478f0 * exp(0.99889f0 + 0.19422f0*log(cl) + 0.59423f0*log(b) -
                             0.09078f0*log(H) - 0.02341f0*log(barea))
        D < 1f0 && (cw *= D / 1f0); cw > 40f0 && (cw = 40f0)
    elseif eq == "20203"       # DF — Crookston (R1)
        b = D >= 1f0 ? D : 1f0
        cw = 1.01685f0 * exp(1.48372f0 + 0.27378f0*log(cl) + 0.49646f0*log(b) -
                             0.18669f0*log(H) - 0.01509f0*log(barea))
        D < 1f0 && (cw *= D / 1f0); cw > 80f0 && (cw = 80f0)
    elseif eq == "01703"       # GF — Crookston (R1)
        b = D >= 1f0 ? D : 1f0
        cw = 1.0303f0 * exp(1.14079f0 + 0.20904f0*log(cl) + 0.38787f0*log(b))
        D < 1f0 && (cw *= D / 1f0); cw > 40f0 && (cw = 40f0)
    elseif eq == "26303"       # WH — Crookston (R1)
        b = D >= 0.1f0 ? D : 0.1f0
        cw = 1.02460f0 * exp(1.3522f0 + 0.24844f0*log(cl) + 0.412117f0*log(b) -
                             0.104357f0*log(H) + 0.03538f0*log(barea))
        D < 0.1f0 && (cw *= D / 0.1f0); cw > 54f0 && (cw = 54f0)
    elseif eq == "24203"       # RC — Crookston (R1)
        b = D >= 1f0 ? D : 1f0
        cw = 1.03597f0 * exp(1.46111f0 + 0.26289f0*log(cl) + 0.18779f0*log(b))
        D < 1f0 && (cw *= D / 1f0); cw > 45f0 && (cw = 45f0)
    elseif eq == "10803"       # LP — Crookston (R1)
        b = D >= 0.7f0 ? D : 0.7f0
        cw = 1.03992f0 * exp(1.58777f0 + 0.30812f0*log(cl) + 0.64934f0*log(b) - 0.38964f0*log(H))
        D < 0.7f0 && (cw *= D / 0.7f0); cw > 40f0 && (cw = 40f0)
    elseif eq == "09303"       # ES — Crookston (R1)
        b = D >= 0.1f0 ? D : 0.1f0
        cw = 1.02687f0 * exp(1.28027f0 + 0.2249f0*log(cl) + 0.47075f0*log(b) - 0.15911f0*log(H))
        D < 0.1f0 && (cw *= D / 0.1f0); cw > 40f0 && (cw = 40f0)
    elseif eq == "01903"       # AF — Crookston (R1)
        b = D >= 0.1f0 ? D : 0.1f0
        cw = 1.02886f0 * exp(1.01255f0 + 0.30374f0*log(cl) + 0.37093f0*log(b) - 0.13731f0*log(H))
        D < 0.1f0 && (cw *= D / 0.1f0); cw > 30f0 && (cw = 30f0)
    elseif eq == "12203"       # PP — Crookston (R1)
        b = D >= 2f0 ? D : 2f0
        cw = 1.02687f0 * exp(1.49085f0 + 0.1862f0*log(cl) + 0.68272f0*log(b) - 0.28242f0*log(H))
        D < 2f0 && (cw *= D / 2f0); cw > 46f0 && (cw = 46f0)
    elseif eq == "26405"       # OT/MH — Crookston (R6) model 2 (BF=1 for KT region 1)
        el = EL < 10f0 ? 10f0 : (EL > 79f0 ? 79f0 : EL)
        b  = D >= 1f0 ? D : 1f0
        cw = 3.7854f0 * (b^0.54684f0) * (H^(-0.12954f0)) * (cl^0.16151f0) *
             ((barea + 1f0)^0.03047f0) * (exp(el)^(-0.00561f0))
        D < 1f0 && (cw *= D / 1f0); cw > 45f0 && (cw = 45f0)
    end
    return cw
end

# --- SVOBOL/SVCROL: two stem circles overlap iff (r1+r2)² > centre-distance² (svcrol.f:15-22) ---
@inline function _svs_circles_overlap(x1::Float32, y1::Float32, r1::Float32,
                                      x2::Float32, y2::Float32, r2::Float32)::Bool
    d2 = (x2 - x1)^2 + (y2 - y1)^2
    rr = r1 + r2
    return rr * rr > d2
end

# --- SVGTPT rectangle branch (IPLGEM<2): two svrann draws map to (x,y) in [x1,x2]×[y1,y2] ---
@inline function _svgtpt_rect(rng::FVSRng, x1::Float32, x2::Float32, y1::Float32, y2::Float32)
    xr = svrann!(rng)
    x  = x1 + (x2 - x1) * xr
    yr = svrann!(rng)
    y  = y1 + (y2 - y1) * yr
    return (x, y)
end

# --- SVS object (SVDATA.F77 /SVOBJ/): one visualization stem. IS2F is the tree-record pointer
#     (remapped by SVTRIP on tripling / SVRMOV on removal); (x,y)=XSLOC/YSLOC persist across cycles.
#     iobjtp mirrors IOBJTP: 1 = live green tree (is2f→tree record), 2 = standing snag
#     (is2f→dead-slot IDEAD, displayed as tree#, with `snag` holding the frozen SVSSnag), 0 = removed.
mutable struct SVSObj
    species::Int          # ISP snapshot (fixed; the remapped record is always the same species)
    is2f::Int             # IS2F: tree-record index (live) or dead-slot IDEAD (snag); 0 ⇒ removed
    x::Float32            # XSLOC(ISVOBJ)
    y::Float32            # YSLOC(ISVOBJ)
    iobjtp::Int           # IOBJTP: 1 live, 2 standing snag, 0 removed
    snag::Any             # the frozen SVSSnag when iobjtp==2 (nothing otherwise)
    SVSObj(species::Integer, is2f::Integer, x::Real, y::Real) =
        new(Int(species), Int(is2f), Float32(x), Float32(y), 1, nothing)
end

"""
    svs_place_objects(s) -> Vector{SVSObj}

SVSTART inventory placement (cycle-0): SVGTPL(IPLGEM=0)+SVESTB(integer)+SVGTPT, advancing the SVS
random stream (seeded here from the main stream at the SVSTART seam). Returns the persistent SVS
object list (IOBJTP=1 live green trees) in placement order — the array whose (x,y) survives every
later cycle while IS2F is remapped by tripling/removal. Chunk assumptions (asserted): IPLGEM=0,
integer TPA. `s` must be at the inventory state (post `setup_growth!`/`compute_volumes!`, pre-growth).
"""
function svs_place_objects(s::StandState)::Vector{SVSObj}
    t = s.trees
    @assert Int(s.control.svs_iplgem) == 0 "svs: only IPLGEM=0 supported"

    # --- SVSTART seed: SVS stream starts at the MAIN stream's current s0 (svstart.f:50) ---
    svs_seed!(s.rng)

    side = _SVS_SIDE_IMPERIAL          # imperial (IMETRIC=0); one square acre per subplot

    objs = SVSObj[]
    # --- SVESTB (integer path): NTOADD(i) = ifix(PROB(i)); FRACAD≡0 ⇒ no lottery START draw ---
    for i in 1:t.n
        prob = t.tpa[i]
        ntoadd = trunc(Int, prob)                         # IFIX(PROB) with NOBTS=0
        # integer-TPA guard. A fractional remainder means the lottery path (chunk 1) is needed.
        @assert abs(prob - ntoadd) < 1f-4 "svs: non-integer TPA on record $i (PROB=$prob) needs the lottery (chunk 1)"
        diai = t.dbh[i] * _SVS_DIAI_FACTOR                # stem radius (ft)
        itre = t.plot_id[i]
        for _ in 1:ntoadd
            x = 0f0; y = 0f0
            ncks = 0
            while true
                (x, y) = _svgtpt_rect(s.rng, 0f0, side, 0f0, side)
                # SVOBOL overlap reject vs already-placed live objects on the SAME subplot (svestb.f:323-348).
                overlap = false
                for o in objs
                    t.plot_id[o.is2f] == itre || continue
                    rj = t.dbh[o.is2f] * _SVS_DIAI_FACTOR
                    if _svs_circles_overlap(x, y, diai, o.x, o.y, rj)
                        overlap = true
                        break
                    end
                end
                if overlap
                    ncks += 1
                    ncks > 40 && break                    # give up after 40 tries; place anyway (svestb.f:343)
                    continue
                end
                break
            end
            push!(objs, SVSObj(Int(t.species[i]), i, x, y))
        end
    end
    return objs
end

"""
    svs_svtrip!(objs, nlive)

SVTRIP (svtrip.f): update the object→record pointers for record tripling. FVS TRIPLE splits every
live base record `i` (1..`nlive`) into (central=i, upper=nlive+2i-1, lower=nlive+2i); this mirrors
svtrip.f, distributing the NR objects that pointed at base `i` as ⌊0.6·NR+0.5⌋ to central, then
⌊0.625·rem+0.5⌋ to upper, the rest to lower (identical to triple_records!'s .60/.25/.15 TPA split).
(x,y) are untouched. No-op for a base with ≤1 object (svtrip.f:44). Faithful to the physical
`nlive+2i-1 / nlive+2i` append order in triple_records! (diameter_growth.jl).
"""
function svs_svtrip!(objs::Vector{SVSObj}, nlive::Integer)
    for i in 1:nlive
        itr = (i, nlive + 2i - 1, nlive + 2i)             # (central, upper, lower) — TRIPLE layout
        # count NR = live objects currently pointing at base record i (svtrip.f:40-42)
        nr = 0
        for o in objs
            (o.iobjtp == 1 && o.is2f == i) && (nr += 1)
        end
        nr <= 1 && continue
        nrt1 = trunc(Int, nr * 0.6f0 + 0.5f0)             # IFIX(FLOAT(NR)*.6+.5)
        nrt2 = max(0, trunc(Int, (nr - nrt1) * 0.625f0 + 0.5f0))
        nrt3 = max(0, nr - nrt1 - nrt2)
        nrt2 <= 0 && continue                             # svtrip.f:52
        nrt = (nrt1, nrt2, nrt3)
        slot = 1; ntoc = nrt1
        for o in objs
            (o.iobjtp == 1 && o.is2f == i) || continue
            o.is2f = itr[slot]
            ntoc -= 1
            if ntoc == 0
                slot += 1
                slot > 3 && break
                ntoc = nrt[slot]
                ntoc == 0 && break
            end
        end
    end
    return objs
end

# Expand a per-parent cycle-start snapshot `snap` (length `nlive_pre`) to the post-growth record set
# (length `N`). If TRIPLE fired (N==3·nlive_pre), each parent i seeds central=i, upper=nlive+2i-1,
# lower=nlive+2i (triple_records! layout) with the SAME cycle-start value (children copy the parent).
# Untripled (N==nlive_pre) ⇒ passthrough.
function _svs_expand(snap::AbstractVector{T}, nlive_pre::Integer, N::Integer)::Vector{T} where {T}
    out = Vector{T}(undef, N)
    if N == nlive_pre
        @inbounds for i in 1:N; out[i] = snap[i]; end
    else                                            # tripled
        @inbounds for i in 1:nlive_pre
            out[i] = snap[i]
            out[nlive_pre + 2i - 1] = snap[i]
            out[nlive_pre + 2i]     = snap[i]
        end
    end
    return out
end

"""
    svs_mortality_snags!(objs, s, snap, iyoftd, ndead) -> Int

The SVS mortality→snag seam (gradd.f:164 SVMORT(0,WK2,·) → svmort.f → svrmov.f → svsnad.f), KT
normal-mortality (ISWTCH=2), standing-dead (FALLDIR=-1), non-fire path. Runs AFTER `grow_cycle!` +
`svs_svtrip!` for a cycle, on the post-tripling record set, reconstructing the oracle SVMORT inputs
from the cycle-start snapshot and the engine's per-record mortality:

  * PROB(r)   = `t.tpa[r] + t.mort_pa[r]`   (cycle-start TPA, pre-mortality — bit-exact vs FVS)
  * REMOVE(r) = `t.mort_pa[r]`              (WK2 periodic mortality)
  * DBH/HT/ICR/CRWDTH(r) = cycle-start dims (`snap`, tripled-expanded; snag freezes these)

Ports SVRMOV faithfully: NOBTS per record (live objects) → integer NTOADD=IFIX(PROB−REMOVE)+FRACAD
→ the per-plot fractional-add lottery (RDPSRT-descending on SCORE=ISP·1000+DBH, cumulative RMVSUM,
systematic sampling with ONE `svrann!` draw per qualifying plot) → NTORMV=NOBTS−NTOADD → the
deficit reassignment (same-species, DBH-within-25%, RDPSRT-neighbour search) → ONE `svrann!` draw
per record with NTORMV>0 (SCORE, consumed but unused on the mortality path) → object removal in
object order (IOBJTP→0) → SVSNAD enough-room snag creation (frozen dims, ISTATUS=1, IYRCOD=iyoftd,
FALLDIR=−1). Converts each removed object in place to a standing snag (iobjtp=2, is2f=IDEAD).

The two svrann draw points keep the SVS stream byte-synced with the oracle (verified: draws #17/#18/#19
= 0.4955/0.1365/0.2479). `ndead` is the running dead-slot counter (Ref); returns snags created.
Gated on SUM(REMOVE)>0 (svmort.f:125) — a mortality-free cycle draws nothing. Single SVS plot per
inventory (ISVINV assert) — the multi-subplot lottery is a later chunk.
"""
function svs_mortality_snags!(objs::Vector{SVSObj}, s::StandState,
                              snap::NamedTuple, iyoftd::Integer,
                              ndead::Base.RefValue{Int})::Int
    t = s.trees
    N = t.n
    nlp = snap.nlive_pre
    # --- reconstruct the oracle SVMORT per-record inputs on the post-tripling set ---
    prob   = Vector{Float32}(undef, N)
    remove = Vector{Float32}(undef, N)
    @inbounds for r in 1:N
        remove[r] = t.mort_pa[r]
        prob[r]   = t.tpa[r] + t.mort_pa[r]
    end
    # SUM>0 gate (svmort.f:125): no mortality ⇒ SVRMOV not called ⇒ no draws.
    sumrmv = 0f0; @inbounds for r in 1:N; sumrmv += remove[r]; end
    sumrmv > 0f0 || return 0
    dbh   = _svs_expand(snap.dbh, nlp, N)
    ht    = _svs_expand(snap.ht,  nlp, N)
    icr   = _svs_expand(snap.icr, nlp, N)
    cw    = _svs_expand(snap.cw,  nlp, N)
    isp   = Int[Int(t.species[r]) for r in 1:N]           # species (unchanged by growth)
    itre  = Int[Int(t.plot_id[r]) for r in 1:N]           # subplot id
    @assert length(unique(itre)) == 1 || maximum(itre) == 1 "svs snag #2: multi-subplot lottery is a later chunk"

    # --- NOBTS: live objects per record (svrmov.f:138-147) ---
    nobts = zeros(Int, N)
    for o in objs
        (o.iobjtp == 1 && o.is2f > 0) || continue
        nobts[o.is2f] += 1
    end
    # --- integer NTOADD + FRACAD (svrmov.f:151-183) ---
    ntoadd = zeros(Int, N)
    fracad = zeros(Float32, N)
    @inbounds for i in 1:N
        xm = prob[i] - remove[i]
        if xm >= 0f0
            (xm + 1f0 > 0f0) && (ntoadd[i] = trunc(Int, xm))     # IFIX(XM)
            fracad[i] = xm - Float32(ntoadd[i])
            fracad[i] < 0.00001f0 && (fracad[i] = 0f0)
        else
            fracad[i] = 0f0
        end
    end
    score = Float32[Float32(isp[i]) * 1000f0 + dbh[i] for i in 1:N]

    # --- per-plot fractional-add lottery (svrmov.f:198-294): one svrann! per qualifying plot ---
    nplot = maximum(itre)
    for iplt in 1:nplot
        ntrplt = Int32[i for i in 1:N if itre[i] == iplt]         # records on this plot (record order)
        n = length(ntrplt)
        n == 0 && continue
        sub = Float32[score[Int(r)] for r in ntrplt]              # RDPSRT(N,SCORE,NTRPLT,.FALSE.) — descending
        subord = collect(Int32, 1:n)
        _rdpsrt!(sub, subord; lseq = false)
        sorted = Int[Int(ntrplt[Int(subord[k])]) for k in 1:n]    # record indices, descending score
        # cumulative RMVSUM in sorted order (svrmov.f:230-240)
        rmvsum = zeros(Float32, N)
        i1 = sorted[1]; rmvsum[i1] = fracad[i1]
        for k in 2:n
            i = sorted[k]; rmvsum[i] = fracad[i] + rmvsum[i1]; i1 = i
        end
        rmvsum[sorted[end]] == 0f0 && continue                   # no fractional removal (svrmov.f:244)
        nta  = 0; for r in ntrplt; nta += ntoadd[Int(r)]; end     # svrmov.f:248-251
        wrmv = 0f0; for r in ntrplt; wrmv += remove[Int(r)]; end  # WPP(IPLT)
        ppp  = 0f0; for r in ntrplt; ppp  += prob[Int(r)];   end  # PPP(IPLT)
        toadd = ppp - wrmv - Float32(nta)                        # svrmov.f:260
        toadd < 0.001f0 && continue                              # svrmov.f:263
        itoadd = trunc(Int, toadd + 0.5f0)
        xint   = rmvsum[sorted[end]] / toadd
        x      = svrann!(s.rng)                                  # <<< RNGA draw (svrmov.f:266)
        start  = xint * x
        xprev  = 0f0; iadd = 0
        for k in 1:n                                             # svrmov.f:276-292
            i = sorted[k]
            if start > xprev && start < rmvsum[i]
                start += xint
                ntoadd[i] += 1
                iadd += 1
                xprev = rmvsum[i]
                iadd >= itoadd && break
            else
                xprev = rmvsum[i]
            end
        end
    end

    # --- NTORMV (svrmov.f:298-300) ---
    ntormv = Int[nobts[i] - ntoadd[i] for i in 1:N]
    # --- deficit reassignment (svrmov.f:337-486, ISWTCH=2: no add-tree branch) ---
    nt2 = collect(Int32, 1:N)
    _rdpsrt!(copy(score), nt2; lseq = true)                     # RDPSRT(ITRN,SCORE,NTRPLT,.TRUE.) — descending
    for kk in 1:N
        i = Int(nt2[kk])
        ntormv[i] < 0 || continue
        idel = -1; ibot = 0; itop = 0
        while true
            (itop + ibot == 2) && break                         # svrmov.f:350
            if idel < 0
                idel = -idel
                (idel + kk > N) && (itop = 1)
                itop == 1 && continue
            else
                idel = -(idel + 1)
                (idel + kk < 1) && (ibot = 1)
                ibot == 1 && continue
            end
            ii = Int(nt2[kk + idel])                            # candidate (svrmov.f:370)
            ntormv[ii] < 1 && continue                          # needs surplus (svrmov.f:388)
            isp[ii] == isp[i] || continue                       # same species (svrmov.f:392)
            (min(dbh[ii], dbh[i]) < max(dbh[ii], dbh[i]) * 0.75f0) && continue  # DBH within 25% (svrmov.f:396)
            ntormv[ii] -= 1; ntormv[i] += 1                     # svrmov.f:406-407
            ntormv[i] < 0 && continue                           # still deficit ⇒ keep searching
            break
        end
    end

    # --- one svrann! per record with NTORMV>0 (svrmov.f:491-498): SCORE, unused on the mort path ---
    @inbounds for ii in 1:N
        ntormv[ii] > 0 && svrann!(s.rng)                        # <<< RNGB draw
    end

    # --- removal loop (svrmov.f:502-544, ISWTCH=2 ELSE branch): objects in object order ---
    torem = Tuple{SVSObj,Int}[]
    for o in objs
        (o.iobjtp == 1 && o.is2f > 0) || continue
        isi = o.is2f
        ntormv[isi] > 0 || continue
        o.iobjtp = 0                                            # IOBJTP=0 (removed → snag)
        ntormv[isi] -= 1
        push!(torem, (o, isi))                                  # ISNADD, in object order
    end
    isempty(torem) && return 0

    # --- SVSNAD enough-room branch (svsnad.f:151-301): create standing-dead snags. Open dead-slots
    #     (increasing IDEAD) are filled from the END of ISNADD backward (svsnad.f:162 NSNGS-ISNAG),
    #     so the last-removed object gets the smallest IDEAD. Standing (FALLDIR=-1), ISTATUS=1. ---
    ns = length(torem)
    for k in 1:ns
        (o, isi) = torem[ns - k + 1]                            # reversed fill order
        ndead[] += 1
        idead = ndead[]
        sn = SVSSnag(isp[isi],           # ISNSP
                     dbh[isi],           # ODIA
                     ht[isi],            # OLEN
                     dbh[isi],           # SNGDIA (= ODIA, standing)
                     ht[isi],            # SNGLEN
                     cw[isi],            # CRNDIA (CRWDTH)
                     icr[isi],           # CRNRTO (ICR)
                     Int(iyoftd),        # IYRCOD = IY(ITCYC+1)-1
                     1,                  # ISTATUS = 1 (green; ages to 2=red next SVSNAGE year)
                     -1f0,               # FALLDIR = -1 (standing)
                     Int(isi))           # OIDTRE (source record; unused for standing display)
        o.iobjtp = 2
        o.snag   = sn
        o.is2f   = idead                                        # displayed tree# = IDEAD (svsnad.f:296)
    end
    return ns
end

"""
    svs_picture(objs, s; year, msg) -> String

Format one SVS picture (SVOUT header svout.f:243-258 + live-tree object records svout.f:429-440)
from the persistent object list `objs` and the CURRENT tree state `s.trees`. Each object is emitted
with its (remapped) record's DBH/HT/ICR/crown-width and its persistent (x,y). Removed objects
(IS2F=0) are skipped (svout.f:391 `I.GT.0`).
"""
function svs_picture(objs::Vector{SVSObj}, s::StandState; year::Integer,
                     msg::AbstractString, ilyear::Integer = year - 1)::String
    t = s.trees; p = s.plot
    ba = p.basal_area; el = p.elevation
    io = IOBuffer()
    stand = strip(String(p.stand_id))
    print(io, "#TITLE Stand=", stand, " Year=", _svs_i4(year), " ", msg, "\n")
    print(io, "#TREEFORM WEST.TRF\n")     # western cluster
    print(io, "#FORMAT 2\n")
    print(io, "#PLOTSIZE 208.71 208.71\n")
    print(io, "#UNITS ENGLISH\n")
    print(io, ";                  trcl  stus             fang\n")
    print(io, ";species        tr#  |crcl|   dbh   ht lang |edia crd  cr    crd  cr    crd  cr    crd  cr ex mk  xloc    yloc  z\n")
    for o in objs
        o.iobjtp == 0 && continue         # removed object (svout.f:391 I=IS2F; IF I.GT.0)
        if o.iobjtp == 2                   # standing snag (SVOUT IOBJTP=2 branch)
            svs_write_snag!(io, o.snag::SVSSnag, s, o.is2f, o.x, o.y;
                            iyear = year, ilyear = ilyear, xmod = 1f0)
            continue
        end
        o.is2f == 0 && continue           # svout.f:391 I=IS2F(ISVOBJ); IF (I.GT.0)
        rec  = o.is2f
        sp2  = rpad(rstrip(String(s.species.code2[o.species])), 2)   # SPCD: 2-char, left-justified
        icr  = abs(t.crown_pct[rec])
        xicr = Float32(icr) * 0.01f0
        cw   = kt_crown_width(o.species, t.dbh[rec], t.height[rec],
                              Float32(t.crown_pct[rec]), ba, el)      # CW=CRWDTH(I)
        crad = cw / 2f0
        _svs_write_tree!(io, sp2, rec, t.dbh[rec], t.height[rec], crad, xicr, o.x, o.y)
    end
    return String(take!(io))
end

"""
    svs_render_cycle0(s; msg="Inventory conditions") -> String

Cycle-0 inventory picture (SVSTART seam): place objects then format. Kept for the chunk-0 test.
"""
function svs_render_cycle0(s::StandState; msg::AbstractString = "Inventory conditions")::String
    objs = svs_place_objects(s)
    return svs_picture(objs, s; year = current_cycle_year(s), msg = msg)
end

# --- SVOUT format 30 (svout.f:429-431/439) — shared by live trees (IOBJTP=1) and snags (IOBJTP=2):
#     (A,T16,I5,I3,2I2,F6.1,F6.0,I2,I4,I2,4(F6.1,1X,F4.2),2I2,2F8.2,I2)
#     A=SPCD  I5=tree#  I3=class(ITC)  I2 I2=(0,IPS)  F6.1=DBH  F6.0=HT  I2=lean(0)  I4=dir(IDIR)
#     I2=edia(0)  4×(F6.1 crownRad, 1X, F4.2 crownRatio)  I2 I2=(ex=1,mk=0)  F8.2 F8.2=xloc,yloc  I2=z(0)
# Live trees pass itc=0/ips=1/idir=0; standing snags pass their status-mapped class (98/94/…), ips=1, idir=0.
function _svs_write_tree!(io::IO, sp2::AbstractString, rec::Integer, dbh::Float32, ht::Float32,
                          crad::Float32, xicr::Float32, x::Float32, y::Float32;
                          itc::Integer = 0, ips::Integer = 1, idir::Integer = 0)
    print(io, rpad(sp2, 15))                       # A + T16 (SPCD in cols 1-2, next field at col 16)
    print(io, _svs_i(rec, 5))                      # I5 tree#
    print(io, _svs_i(itc, 3))                      # I3 tree class (0 live; 98/94/… snag)
    print(io, _svs_i(0, 2), _svs_i(ips, 2))        # crown class, plant status (IPS)
    print(io, _svs_f(dbh, 6, 1))                   # F6.1 dbh
    print(io, _svs_f0(ht, 6))                      # F6.0 ht (trailing '.')
    print(io, _svs_i(0, 2), _svs_i(idir, 4), _svs_i(0, 2))   # lean, felling dir, small-end dia
    for _ in 1:4
        print(io, _svs_f(crad, 6, 1), " ", _svs_f(xicr, 4, 2))   # 4×(crown radius, crown ratio)
    end
    print(io, _svs_i(1, 2), _svs_i(0, 2))          # expansion factor (1), marking status (0)
    print(io, _svs_f(x, 8, 2), _svs_f(y, 8, 2))    # xloc, yloc
    print(io, _svs_i(0, 2))                        # z
    print(io, "\n")
end

# --- minimal Fortran edit-descriptor emitters (all fixed-width, blank-filled) ---
_svs_i(v::Integer, w::Integer) = lpad(string(v), w)
_svs_i4(v::Integer) = lpad(string(v), 4, '0')                       # I4.4
_svs_f(v::Real, w::Integer, d::Integer) =
    lpad(_fmt_fixed(Float64(v), d), w)
"Fortran F<w>.0 — round to integer, print WITH a trailing decimal point (e.g. 73.0 → \"   73.\")."
_svs_f0(v::Real, w::Integer) = lpad(string(round(Int, Float64(v))) * ".", w)

"""
    svs_write_cycle0_files(stem, s; imageno=1, msg="Inventory conditions") -> String

The SVSTART seam (fvs.f:333, cycle-0 inventory picture): write the two SVS files FVS emits with
the default JSVPIC=91 — `<stem>_NNN.svs` (the picture: SVOUT header + per-tree object records) and
`<stem>_index.svs` (`#TREELISTINDEX` + one `"Stand=… Year=…" "<pic>"` index line, svout.f:200/211).
Returns the picture-file path. Gated by the caller on `control.svs_on` (SVS keyword seen).
"""
function svs_write_cycle0_files(stem::AbstractString, s::StandState;
                                imageno::Integer = 1, msg::AbstractString = "Inventory conditions")
    body    = svs_render_cycle0(s; msg = msg)
    picfile = string(stem, "_", lpad(string(imageno), 3, '0'), ".svs")
    open(picfile, "w") do io; write(io, body); end
    stand = strip(String(s.plot.stand_id))
    year  = current_cycle_year(s)
    open(string(stem, "_index.svs"), "w") do io
        print(io, "#TREELISTINDEX\n")
        print(io, "\"Stand=", stand, " Year=", _svs_i4(year), " ", msg, "\" \"", basename(picfile), "\"\n")
    end
    return picfile
end

"""
    svs_project!(stem, s; fint=10f0) -> Vector{String}

Full multi-cycle SVS data path for the standard projection: mirror the engine's SVS seams —
SVSTART (fvs.f:333, inventory picture), GRINCR (grincr.f:277, `IF ICYC>1` "Beginning of cycle"
picture) each later cycle, and MAIN (fvs.f:453, "End of projection") — driving `grow_cycle!`
between them and remapping the persistent object list through SVTRIP on every tripling cycle.

Writes `<stem>_NNN.svs` for each picture and the accumulated `<stem>_index.svs` (`#TREELISTINDEX`
+ one line per picture). Returns the picture-file paths. `s` must be at the inventory state
(post `setup_growth!`/`compute_volumes!`, pre-growth). Cycle count comes from the keyword schedule
(`s.control`); tripling fires only for the first ICL4 cycles, exactly as the engine decides it.
"""
function svs_project!(stem::AbstractString, s::StandState; fint::Float32 = 10f0)
    t = s.trees
    stand = strip(String(s.plot.stand_id))
    ncyc  = Int(s.control.ncycle)

    # --- SVSTART (cyc0 inventory picture) ---
    objs   = svs_place_objects(s)
    index  = Tuple{Int,String,String}[]     # (year, msg, picfile-basename)
    pics   = String[]
    imageno = 1
    invyear = current_cycle_year(s)
    ilyear  = invyear                        # ILYEAR: year of the last-written picture (SVOUT updates it)
    ndead   = Ref(0)                         # NDEAD: running SVS dead-slot (IDEAD) counter
    picfile = string(stem, "_", lpad(string(imageno), 3, '0'), ".svs")
    open(picfile, "w") do io; write(io, svs_picture(objs, s; year = invyear,
                                                     msg = "Inventory conditions", ilyear = ilyear)); end
    push!(pics, picfile)
    push!(index, (invyear, "Inventory conditions", basename(picfile)))

    # --- projection loop (mirrors fvs.f: for ICYC=1..NUMCYCLE) ---
    for icyc in 1:ncyc
        # GRINCR seam (grincr.f:277): IF ICYC>1 emit "Beginning of cycle" BEFORE growing this cycle.
        if icyc > 1
            imageno += 1
            byear = current_cycle_year(s)
            picfile = string(stem, "_", lpad(string(imageno), 3, '0'), ".svs")
            open(picfile, "w") do io; write(io, svs_picture(objs, s; year = byear,
                                                             msg = "Beginning of cycle", ilyear = ilyear)); end
            push!(pics, picfile)
            push!(index, (byear, "Beginning of cycle", basename(picfile)))
            ilyear = byear
        end
        # Snapshot cycle-START dims for the SVMORT seam (the snag freezes these; the lottery/deficit
        # search use DBH). Captured BEFORE grow_cycle! grows the records; expanded to the tripled set.
        nlive_pre = t.n
        ba = s.plot.basal_area; el = s.plot.elevation
        snap = (nlive_pre = nlive_pre,
                dbh = Float32[t.dbh[i]        for i in 1:nlive_pre],
                ht  = Float32[t.height[i]     for i in 1:nlive_pre],
                icr = Float32[Float32(t.crown_pct[i]) for i in 1:nlive_pre],
                cw  = Float32[kt_crown_width(Int(t.species[i]), t.dbh[i], t.height[i],
                                             Float32(t.crown_pct[i]), ba, el) for i in 1:nlive_pre])
        iyoftd = current_cycle_year(s) + round(Int, fint) - 1   # IYOFTD = IY(ITCYC+1)-1 (death year)
        # grow one cycle; TRIPLE (if it fired) split every live base record → remap the objects.
        grow_cycle!(s; fint = fint)
        compute_volumes!(s)
        t.n > nlive_pre && svs_svtrip!(objs, nlive_pre)   # tripling appended records ⇒ SVTRIP
        # SVMORT seam (gradd.f:164): mortality→standing-snag object selection + creation.
        svs_mortality_snags!(objs, s, snap, iyoftd, ndead)
    end

    # --- MAIN seam (fvs.f:453): "End of projection" picture ---
    imageno += 1
    eyear = current_cycle_year(s)
    picfile = string(stem, "_", lpad(string(imageno), 3, '0'), ".svs")
    open(picfile, "w") do io; write(io, svs_picture(objs, s; year = eyear,
                                                     msg = "End of projection", ilyear = ilyear)); end
    push!(pics, picfile)
    push!(index, (eyear, "End of projection", basename(picfile)))

    # --- accumulated _index.svs ---
    open(string(stem, "_index.svs"), "w") do io
        print(io, "#TREELISTINDEX\n")
        for (yr, msg, pf) in index
            print(io, "\"Stand=", stand, " Year=", _svs_i4(yr), " ", msg, "\" \"", pf, "\"\n")
        end
    end
    return pics
end

# =============================================================================
# SVS snag data path (mortality → standing-dead snag), KT non-fire path.
#
# Ported from:
#   common/SVDEAD.F77   — the SVS snag list (per dead-slot IDEAD)
#   base/svsnage.f      (SVSNAGE)  — snag aging: crown-diameter decay, crown-ratio, FMSNGHT
#                                    height loss, FMSNGDK hard→soft decay, status progression
#   fire/vbase/fmsnght.f(FMSNGHT)  — KT falls to CASE DEFAULT (western): per-year top-breakage loss
#   fire/vbase/fmsngdk.f(FMSNGDK)  — KT CASE DEFAULT: years-since-death to become a soft snag
#   {kt}/fmvinit.f      — KT snag params (HTR1/HTR2/HTXSFT, per-species HTX & DECAYX)
#   vbase/svout.f:445-540 (SVOUT)  — snag object emission: status→display class, SNCRDI/2, SNCRTO
#
# SCOPE (this sub-chunk): the standing-dead (FALLDIR=-1), non-fire, non-salvage snag DISPLAY +
# AGING, validated bit-exact vs the LIVE FVSkt snag record (measured via an instrumented svout.f
# dump on kt2c cyc2). The mortality object-SELECTION (svrmov.f) + snag CREATION (svsnad.f) + the
# engine SVMORT seam (gradd.f:164, cycle-start dims, pre-tripling record set) are the NEXT sub-chunk.
# =============================================================================

# --- SVDEAD.F77: one SVS snag record (a dead-slot IDEAD). Dims are frozen at death; SVSNAGE ages them.
mutable struct SVSSnag
    sp::Int          # ISNSP    — species sequence number
    odia::Float32    # ODIA     — DBH at death (frozen)
    olen::Float32    # OLEN     — height at death (frozen reference for crown-ratio + height loss)
    sngdia::Float32  # SNGDIA   — current snag diameter (standing ⇒ = ODIA)
    snglen::Float32  # SNGLEN   — current snag length (height-loss-aged)
    crndia::Float32  # CRNDIA   — crown diameter at death (CRWDTH)
    crnrto::Float32  # CRNRTO   — crown ratio at death (ICR), percent
    iyrcod::Int      # IYRCOD   — year of tree death
    istatus::Int     # ISTATUS  — status code (1 green,2 red,3 hard-grey,4 soft-grey,5/6 burn,90-92 WWPB)
    falldir::Float32 # FALLDIR  — fall direction (-1 ⇒ standing)
    oidtre::Int      # OIDTRE   — source tree id
end

# KT snag height/decay params ({kt}/fmvinit.f). HTX/DECAYX indexed by KT species 1..11
# (1 WP 2 WL 3 DF 4 GF 5 WH 6 RC 7 LP 8 ES 9 AF 10 PP 11 MH). All four HTX(sp,1:4) share one value.
const _KT_SNAG_HTR1   = 0.0228f0
const _KT_SNAG_HTR2   = 0.01f0
const _KT_SNAG_HTXSFT = 2.0f0
const _KT_SNAG_HTX    = (0.9f0,0.9f0,0.9f0,1.1f0,1.1f0,1.1f0,1.1f0,1.1f0,1.1f0,1.0f0,1.0f0)
const _KT_SNAG_DECAYX = (1.1f0,1.1f0,1.1f0,0.9f0,0.9f0,0.9f0,0.9f0,0.9f0,0.9f0,1.0f0,1.0f0)

# FMSNGHT KT CASE DEFAULT (fmsnght.f:153-160), one year: shrink snag height by the top-breakage rate.
# `ihard` picks HTINDX1/SFTMULT=1 (hard) vs HTINDX2/SFTMULT=HTXSFT (soft); regime = htcurr vs 0.5·htd.
function _kt_fmsnght(sp::Integer, htd::Float32, htcurr::Float32, ihard::Bool)::Float32
    htx = _KT_SNAG_HTX[sp]
    sftmult = ihard ? 1f0 : _KT_SNAG_HTXSFT
    htnew = htcurr > 0.5f0*htd ?
        htcurr*(1f0 - _KT_SNAG_HTR1*htx*sftmult) :
        htcurr*(1f0 - _KT_SNAG_HTR2*htx*sftmult)
    htnew < 1.5f0 && (htnew = 0f0)                     # fmsnght.f:164 — <1.5 ft ⇒ 'fuel', snag gone
    return htnew
end

# FMSNGDK KT CASE DEFAULT (fmsngdk.f): years-since-death for a snag of diameter `d` to become soft.
_kt_fmsngdk(sp::Integer, d::Float32, xmod::Float32)::Float32 =
    (1.24f0*_KT_SNAG_DECAYX[sp]*d + 13.82f0*_KT_SNAG_DECAYX[sp])*xmod

"""
    svsnage_standing!(sn, iyear, ilyear, xmod) -> (sndi, snht, sncrdi, sncrto)

SVSNAGE (svsnage.f) for a STANDING (FALLDIR=-1), non-fire, non-salvage snag: age `sn` to `iyear`
and return the display diameter / height / crown-diameter / crown-ratio(%). Mirrors svsnage.f:
crown diameter decays 0.90/yr (`CRNDIA·0.90^ITIDIF`); crown ratio uses the PRE-height-loss length
(`(OLEN/SNGLEN)·(CRNRTO·.01−1)+1`, ×100, ICYC>0 branch); the year loop (svsnage.f:330-404) applies
one FMSNGHT height loss + status progression (green/red→red<2yr else hard-grey) + FMSNGDK hard→soft
per year from `max(ilyear+1, iyrcod+1)` to `iyear`. Standing snags keep their death diameter.
"""
function svsnage_standing!(sn::SVSSnag, iyear::Integer, ilyear::Integer, xmod::Float32)
    itidif = iyear > sn.iyrcod ? iyear - sn.iyrcod : 0
    snht = sn.snglen                                   # svsnage.f:102 (pre-loss length)
    sndi = sn.sngdia
    sncrdi = sn.crndia * 0.9f0^itidif                  # svsnage.f:104
    # crown ratio (ICYC>0 branch, svsnage.f:108-110) uses the PRE-loss SNHL length
    if snht > 0.5f0
        sncrto = ((sn.olen/snht)*(sn.crnrto*0.01f0 - 1f0) + 1f0)*100f0
    else
        sncrto = 0f0
    end
    (sncrto <= 0f0 || sncrdi <= 0f0) && (sncrdi = 0f0; sncrto = 0f0)
    # up-to-date? (svsnage.f:127)
    (iyear <= ilyear || itidif <= 0) && return (sndi, snht, sncrdi, sncrto)
    # standing height-loss + status-progression loop (svsnage.f:330-404), FALLDIR=-1
    for icuryr in max(ilyear+1, sn.iyrcod+1):iyear
        ihard = sn.istatus != 4
        snht = _kt_fmsnght(sn.sp, sn.olen, sn.snglen, ihard)
        sn.snglen = snht
        if sn.istatus == 1 || sn.istatus == 2         # CASE(1,2): green/red hard snags
            sn.istatus = itidif < 2 ? 2 : 3
        end
        if sn.istatus != 4                            # FMSNGDK hard→soft (svsnage.f:398-403)
            dktime = _kt_fmsngdk(sn.sp, sn.sngdia, xmod)
            (icuryr - sn.iyrcod) >= dktime && (sn.istatus = 4)
        end
    end
    sn.sngdia = sn.odia; sndi = sn.odia                # svsnage.f:411 — standing ⇒ no diameter change
    return (sndi, snht, sncrdi, sncrto)
end

# SVOUT standing-snag status → display tree-class ITC (svout.f:505-540, FALLDIR=-1 branch).
function svs_snag_class(istatus::Integer)::Int
    istatus == 2  && return 98    # red tree (recently dead standing)
    istatus == 3  && return 94    # hard grey snag
    istatus == 4  && return 94    # soft grey snag
    istatus == 5  && return 97    # recently burned
    istatus == 6  && return 96    # older burned grey
    istatus == 90 && return 90
    istatus == 91 && return 91
    istatus == 92 && return 92
    return 0                      # status 1 is a never-displayed check state (svout.f:511)
end

"""
    svs_write_snag!(io, sn, s, rec; iyear, ilyear=iyear-1, xmod=1f0)

Emit one standing-snag object record (SVOUT IOBJTP=2 branch): age `sn` via `svsnage_standing!`,
map its status to the display class, and write the fmt-30 record with `tree# = rec` (the IS2F dead
slot), DBH=SNDI, HT=SNHT, crown radius=SNCRDI/2, crown ratio=SNCRTO·.01, at the snag's (x,y).
A snag aged to SNHT ≤ 0 is dropped (svout.f:485). Returns true if emitted.
"""
function svs_write_snag!(io::IO, sn::SVSSnag, s::StandState, rec::Integer, x::Float32, y::Float32;
                         iyear::Integer, ilyear::Integer = iyear - 1, xmod::Float32 = 1f0)::Bool
    (sndi, snht, sncrdi, sncrto) = svsnage_standing!(sn, iyear, ilyear, xmod)
    snht <= 0f0 && return false
    sp2  = rpad(rstrip(String(s.species.code2[sn.sp])), 2)
    itc  = svs_snag_class(sn.istatus)
    crad = sncrdi / 2f0
    xicr = sncrto * 0.01f0
    _svs_write_tree!(io, sp2, rec, sndi, snht, crad, xicr, x, y; itc = itc, ips = 1, idir = 0)
    return true
end

"Format a Float64 with exactly `d` fractional digits (round-half-away, like Fortran F edit)."
function _fmt_fixed(v::Float64, d::Integer)::String
    neg = v < 0
    a = abs(v)
    scale = 10.0^d
    n = floor(Int, a * scale + 0.5)                # round half up (matches Fortran F for our data)
    ip = n ÷ Int(scale)
    fp = n % Int(scale)
    s = d == 0 ? string(ip) : string(ip, ".", lpad(string(fp), d, '0'))
    return neg ? "-" * s : s
end
