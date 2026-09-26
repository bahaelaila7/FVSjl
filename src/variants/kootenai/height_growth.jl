# =============================================================================
# height_growth.jl (kootenai) — KT large-tree height growth (kt/htgf.f, western Wykoff form).
#
# Deterministic per-tree height increment (NO age-curve, unlike SN; NO GENGYM, unlike CR):
#   CON = HTCON(sp) + H2COF·HT² + HGLD(sp)·ln(D) + HGLH·ln(HT)          (htgf.f:107)
#   HTG = exp(CON + HDGCOF·ln(DG)) + BIAS;  max(0.1);  ·SCALE·XHT·MISHGF (htgf.f:114-121)
#   size cap: HT+HTG > SIZCAP(sp,4) ⇒ HTG = SIZCAP−HT (min 0.1)         (htgf.f:134)
# Per-stand habitat coefficients (HTCONS entry, htgf.f:175): IHT = MAPHAB(ITYPE) (ITYPE = the KT
# habitat type in p.habitat_input), then HGHCH=HGHC(IHT), H2COF=HGH2(IHT), HDGCOF=HGLDD(IHT);
# HTCON(sp) = HGHCH + HGSC(sp) (+ ln(HCOR2) when READCORH set LHCOR2). DG(sp) is the computed
# diameter increment (t.diam_growth) — now bit-exact from chunk 3, so ln(DG) lines up.
# =============================================================================

# HGLD(11): ln(D) coefficient per species (htgf.f:69)
const KT_HGLD = Float32[-.04935,-.3899,-.4574,-.09775,-.1555,-.1219,-.2454,-.5720,-.1997,-.5657,-.1219]
# HGSC(11): per-species intercept added to the habitat intercept HGHCH (htgf.f:196)
const KT_HGSC = Float32[-.5342,.1433,.1641,-.6458,-.6959,-.9941,-.6004,.2089,-.5478,.7316,-.9941]
const KT_HGLH = 0.23315f0     # ln(HT) coefficient (htgf.f:71)
const KT_HTBIAS = 0.4809f0    # average residual added after exp (htgf.f:71,114)
# MAPHAB(30): ITYPE (1..30) -> habitat-coefficient index IHT (1..8) (htgf.f:185-186)
const KT_HTMAPHAB = Int[1,1, 2,2,2,2,2,2,2, 3,3,4,5,6, 7,7,7,7, 4,4,1,4,4, 8,8,8, 1,1,1,1]
# 8 habitat-indexed coefficient vectors (htgf.f:187-194)
const KT_HGHC  = Float32[2.03035, 1.72222, 1.19728, 1.81759, 2.14781, 1.76998, 2.21104, 1.74090]  # HGHCH intercept
const KT_HGLDD = Float32[0.62144, 1.02372, 0.85493, 0.75756, 0.46238, 0.49643, 0.37042, 0.34003]  # HDGCOF (ln DG)
const KT_HGH2  = Float32[-13.358f-5, -3.809f-5, -3.715f-5, -2.607f-5, -5.200f-5, -1.605f-5, -3.631f-5, -4.460f-5]  # H2COF (HT²)

"""
    height_growth!(s, ::Kootenai; scale=fint/10)

KT `HTGF` hook — per-tree periodic height increment into `trees.ht_growth`. Large-tree model for all
records; REGENT (chunk 6) later overrides the small trees. Deterministic (no ZZRAN here).
"""
function height_growth!(s::StandState, ::Kootenai; scale::Float32 = 1.0f0)
    p, t, c, ctl = s.plot, s.trees, s.calib, s.control
    itype = Int(p.habitat_input)
    iht = (1 <= itype <= 30) ? KT_HTMAPHAB[itype] : 1
    hghch  = KT_HGHC[iht]
    h2cof  = KT_HGH2[iht]
    hdgcof = KT_HGLDD[iht]
    cur_year = current_cycle_year(s)
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0.0f0
        t.tpa[i] <= 0.0f0 && continue
        sp = Int(t.species[i])
        d = t.dbh[i]; hti = t.height[i]
        (d <= 0.0f0 || hti <= 0.0f0) && continue
        dg = t.diam_growth[i]
        dg <= 0.0f0 && continue                              # ln(DG) undefined; DGBND keeps DG>0 for grown trees
        # HTCON(sp) = HGHCH + HGSC(sp) (+ ln(HCOR2) when READCORH active) (htgf.f:196-197)
        htcon = hghch + KT_HGSC[sp]
        (ctl.htg_cor2_on && ctl.htg_cor2[sp] > 0.0f0) && (htcon += log(ctl.htg_cor2[sp]))
        con = htcon + h2cof * hti * hti + KT_HGLD[sp] * log(d) + KT_HGLH * log(hti)
        htg = exp(con + hdgcof * log(dg)) + KT_HTBIAS
        htg < 0.1f0 && (htg = 0.1f0)
        xht = active_multiplier(ctl, :htg, sp, cur_year)     # XHMULT (MULTS kind 2); MISHGF=1
        htg = htg * scale * xht
        cap = ctl.sp_size_cap[sp, 4]                         # SIZCAP(sp,4)
        if hti + htg > cap
            htg = cap - hti; htg < 0.1f0 && (htg = 0.1f0)
        end
        t.ht_growth[i] = htg
    end
    return s
end

"""
    kt_triple_htg!(s, stash; scale)

kt/htgf.f:139-161 (LTRIP): each tripled copy gets its own large-tree HTG from the COPY's spread DG,
`EXP(CON+HDGCOF·ALOG(DG(copy)))+BIAS` with the central record's CON. The two copies differ in one detail that is
reproduced as written: the upper copy (ITFN) is floored at 0.1 and then scaled by SCALE·XHT, the lower copy
(ITFN+1) is scaled first and floored after. Neither copy gets MISHGF. height_growth! only computes the central
HTG, so without this both copies inherited it (ktt01 2000: tree 1 LP copies HTG 6.180/6.180 vs live 5.184/7.371).
Deterministic (no draw). Small trees are overwritten afterwards by REGENT (is_small), as in FVS.
"""
function kt_triple_htg!(s::StandState, stash; scale::Float32 = 1.0f0)
    stash === nothing && return s
    p, t, ctl = s.plot, s.trees, s.control
    itype = Int(p.habitat_input)
    iht = (1 <= itype <= 30) ? KT_HTMAPHAB[itype] : 1
    hghch = KT_HGHC[iht]; h2cof = KT_HGH2[iht]; hdgcof = KT_HGLDD[iht]
    cur_year = current_cycle_year(s)
    dgU = stash.dgU; dgL = stash.dgL; htgU = stash.htgU; htgL = stash.htgL; htg_copy = stash.htg_copy
    @inbounds for i in 1:stash.nlive
        t.tpa[i] <= 0.0f0 && continue
        sp = Int(t.species[i])
        d = t.dbh[i]; hti = t.height[i]
        (d <= 0.0f0 || hti <= 0.0f0 || t.diam_growth[i] <= 0.0f0) && continue   # central skipped ⇒ copies flat
        htcon = hghch + KT_HGSC[sp]
        (ctl.htg_cor2_on && ctl.htg_cor2[sp] > 0.0f0) && (htcon += log(ctl.htg_cor2[sp]))
        con = htcon + h2cof * hti * hti + KT_HGLD[sp] * log(d) + KT_HGLH * log(hti)
        xht = active_multiplier(ctl, :htg, sp, cur_year)
        cap = ctl.sp_size_cap[sp, 4]
        e(dgc) = dgc > 0.0f0 ? exp(con + hdgcof * log(dgc)) : 0.0f0     # ALOG(0)=-Inf ⇒ EXP term 0
        hu = e(dgU[i]) + KT_HTBIAS; hu < 0.1f0 && (hu = 0.1f0); hu = hu * scale * xht      # htgf.f:141-143
        (hti + hu > cap) && (hu = max(cap - hti, 0.1f0))
        hl = (e(dgL[i]) + KT_HTBIAS) * scale * xht; hl < 0.1f0 && (hl = 0.1f0)            # htgf.f:152-154
        (hti + hl > cap) && (hl = max(cap - hti, 0.1f0))
        htgU[i] = hu; htgL[i] = hl; htg_copy[i] = true
    end
    return s
end

# ---------------------------------------------------------------------------------------------------------------------
# kt/cratet.f missing-height dub (DO 150 per species + DO 145 dead records) and kt/smhtrg.f. KT is the only CRATET that
# regresses the stand's own heights as H = AA + BB·√D (≥3 trees with H>4.5, NORMHT≥0, D≥3, LHTDRG) and switches to
# the small-tree line SMA·D+4.5 below the crossover SOLN; otherwise (or BB<0) it dubs with the blkdat Wykoff defaults
# EXP(HT1+HT2/(D+1))+4.5 above SMDBH and SMA·D+4.5 below. jl ran KT through the generic Wykoff-AA dub, so every
# missing height — live and dead — came out on the wrong curve (ktt01 dead WL D=34.6: 31.76 vs live 127.53).
# ---------------------------------------------------------------------------------------------------------------------
const KT_BLK_HT1 = Float32[5.186772, 5.054504, 4.876799, 5.063864, 4.927262, 4.881348, 4.777827, 5.079578, 4.930150,
                           5.019903, 4.77951]
const KT_BLK_HT2 = Float32[-10.421897, -8.618720, -9.146701, -9.892393, -8.727470, -9.628555, -6.336417, -10.201550,
                           -8.825197, -12.014787, -9.31743]
const KT_SMA   = Float32[4.43, 7.49, 5.71, 5.00, 5.32, 4.69, 7.49, 4.98, 4.74, 3.83, 3.92]
const KT_SMDBH = Float32[4.17, 5.84, 5.50, 4.74, 4.74, 5.21, 2.42, 4.96, 3.43, 5.94, 6.17]

# cratet.f:120-148 / 520-540 — the top-kill (TKILL) tail shared by the live (:120) and dead (:142) paths.
@inline function _kt_tkill_tail!(t, ii::Int, h::Float32, dead::Bool)
    t.norm_ht[ii] = (dead && t.height[ii] > 0f0) ? trunc(Int32, t.height[ii] * 100f0 + 0.5f0) :
                                                 trunc(Int32, h * 100f0 + 0.5f0)
    if t.trunc[ii] == 0
        if t.height[ii] > 0f0
            t.trunc[ii] = trunc(Int32, 80f0 * t.height[ii] + 0.5f0)
        else
            t.trunc[ii] = trunc(Int32, 80f0 * h + 0.5f0); t.height[ii] = h
        end
    else
        if t.height[ii] > 0f0
            t.height[ii] < Float32(t.trunc[ii]) * 0.01f0 && (t.height[ii] = Float32(t.trunc[ii]) * 0.01f0)
        else
            t.height[ii] = Float32(t.trunc[ii]) * 0.01f0
        end
    end
    Float32(t.norm_ht[ii]) * 0.01f0 < t.height[ii] && (t.norm_ht[ii] = trunc(Int32, t.height[ii] * 100f0))
    return
end

function kt_dub_missing_heights!(s::StandState)
    t = s.trees; n = t.n; nd = Int(t.ndead); lhtdrg = s.control.ht_drag_sp
    order = species_major_order(s)                   # IND1 (live)
    htkeep = zeros(Float32, n + nd); nsp = length(KT_SMA); ispfl = zeros(Int, nsp)   # KT MAXSP=11
    # ---- SMHTRG (cratet.f:156): only when EVERY live record is under 3" ----
    if n > 0 && all(i -> t.dbh[i] < 3f0, 1:n)
        fill!(ispfl, 1)
        aa = 0f0; bb = 0f0                           # statics: a K1<3 species leaves them from the previous one
        for sp in 1:nsp
            k1 = 0; syy = 0f0; sxx = 0f0; sxy = 0f0; sx2 = 0f0; ind2 = Int[]
            for i in order
                Int(t.species[i]) == sp || continue
                h = t.height[i]; d = t.dbh[i]
                if h > 0f0
                    k1 += 1; syy += h; sxx += d; sxy += d * h; sx2 += d * d
                end
                (h > 0f0 && t.norm_ht[i] == 0) && continue
                push!(ind2, i)
            end
            any(i -> Int(t.species[i]) == sp, order) || continue       # ISCT(sp,1)≤0 ⇒ GO TO 500
            if k1 < 3 || !lhtdrg[sp]
                ispfl[sp] = 0
            else
                xn = Float32(k1); den = xn * sx2 - sxx * sxx
                if den == 0f0
                    ispfl[sp] = 0
                else
                    bb = (xn * sxy - sxx * syy) / den
                    aa = (sx2 * syy - sxx * sxy) / den
                    bb <= 1f0 && (ispfl[sp] = 0)
                end
            end
            for i in ind2                                               # DO 400
                ispfl[sp] == 0 && continue
                t.height[i] > 0f0 && continue
                d = t.dbh[i]; h = aa + bb * d
                if h < 1f0 || h > 20f0 * d
                    ispfl[sp] = 0
                elseif d > 0.1f0 && h < 4.5f0
                    ispfl[sp] = 0
                else
                    htkeep[i] = h
                end
            end
        end
    end
    # ---- DO 150 ISPC ----
    ierrck = s.control.kt_cratet_ierrck            # -fno-automatic static: never reset (carried across stands)
    soln = 0f0                                     # static SOLN: the IERRCK path re-uses the last one
    k1 = 0; counter = 0f0                          # statics carried into a species with no live records (:274-277)
    # dead records in FVS order: IREC2→MAXTRE = the reverse of jl's dead storage (n+1:n+nd)
    deads = collect((n + nd):-1:(n + 1))
    function dub_h(sp, d, aa, bb)
        if k1 < 3 || !lhtdrg[sp] || counter == 1f0
            return d >= KT_SMDBH[sp] ? exp(aa + bb / (d + 1f0)) + 4.5f0 : KT_SMA[sp] * d + 4.5f0
        end
        h = aa + bb * sqrt(d)
        step1 = bb * bb - 4f0 * KT_SMA[sp] * (4.5f0 - aa)
        if step1 <= 0f0
            ierrck = Int32(1)
        else
            rad = sqrt(step1)
            soln = ((bb - rad) / (2f0 * KT_SMA[sp]))^2
            soln <= 0f0 && (soln = ((bb + rad) / (2f0 * KT_SMA[sp]))^2)
        end
        (d < soln || (ierrck == 1 && d < KT_SMDBH[sp])) && (h = KT_SMA[sp] * d + 4.5f0)
        return h
    end
    for sp in 1:nsp
        recs = [i for i in order if Int(t.species[i]) == sp]
        local aa::Float32, bb::Float32
        if isempty(recs)
            aa = KT_BLK_HT1[sp]; bb = KT_BLK_HT2[sp]
        else
            counter = 0f0; k1 = 0
            syy = 0f0; sxx = 0f0; sxy = 0f0; sx2 = 0f0; ind2 = Int[]
            for i in recs
                h = t.height[i]; nh = t.norm_ht[i]; d = t.dbh[i]
                if h <= 4.5f0 || nh < 0 || d < 3f0
                    (h > 0f0 && nh == 0) && continue
                    push!(ind2, i)
                else
                    k1 += 1; xx = sqrt(d)
                    syy += h; sxx += xx; sxy += xx * h; sx2 += d
                end
            end
            fitted = false
            if !(k1 < 3 || !lhtdrg[sp])
                xn = Float32(k1)
                z1 = xn * sxy - sxx * syy; z2 = xn * sx2 - sxx * sxx
                bb = (z1 == 0f0 || z2 == 0f0) ? 0f0 : (xn * sxy - sxx * syy) / (xn * sx2 - sxx * sxx)
                z1 = sx2 * syy - sxx * sxy
                aa = (z1 == 0f0 || z2 == 0f0) ? 0f0 : (sx2 * syy - sxx * sxy) / (xn * sx2 - sxx * sxx)
                fitted = bb >= 0f0
            end
            if !fitted
                aa = KT_BLK_HT1[sp]; bb = KT_BLK_HT2[sp]; counter = 1f0      # label 100
            end
            for ii in ind2                                                   # DO 130
                d = t.dbh[ii]; tkill = t.norm_ht[ii] < 0
                local h::Float32
                if d <= 0.1f0
                    h = 1.01f0
                else
                    h = dub_h(sp, d, aa, bb)
                    h < 4.5f0 && (h = 4.5f0)
                    htkeep[ii] < 4.5f0 && (htkeep[ii] = 4.5f0)
                end
                if tkill
                    _kt_tkill_tail!(t, ii, h, false)                       # label 120
                else
                    t.height[ii] = ispfl[sp] == 1 ? htkeep[ii] : h
                end
            end
        end
        for ii in deads                                                      # label 141, DO 145
            Int(t.species[ii]) == sp || continue
            d = t.dbh[ii]; tkill = t.norm_ht[ii] < 0
            if t.height[ii] > 0f0
                tkill && _kt_tkill_tail!(t, ii, 0f0, true)                  # :142 with the measured HT
                continue                                                     # :146
            end
            local h::Float32
            if d <= 0.1f0
                h = 1.01f0
            else
                h = dub_h(sp, d, aa, bb)
                h < 4.5f0 && (h = 4.5f0)
                htkeep[ii] < 4.5f0 && (htkeep[ii] = 4.5f0)
            end
            if tkill
                _kt_tkill_tail!(t, ii, h, true)
            else
                t.height[ii] = ispfl[sp] == 1 ? htkeep[ii] : h
            end
        end
    end
    s.control.kt_cratet_ierrck = ierrck
    return s
end
