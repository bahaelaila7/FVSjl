# =============================================================================
# ontario/small_tree_growth.jl — ON small-tree REGENT (canada/on/regent.f, D < XMAX = 12 cm).
#
# ON's REGENT is base/regent.f (same structure the SN/CS/LS variants use), differing only in the
# HEIGHT increment: ON replaces the LS Chapman-Richards HTCALC-diff with the Penner-2010 ONSTHG
# annual small-tree height model (htcalc.f:387-393 MODE0=9 → ONSTHG·YRS). For a record with
# DBH ≥ XMAX (4.72" = 12 cm) REGENT is a no-op (the large-tree htont/dgf increment stands).
#
# The sub-12cm chain (LSTART=false, LESTB=false, non-calibrated — all multipliers = 1 for ON:
# con=RHCON·exp(HCOR)=1, HGADJ=1, XRHGRO=XRDGRO=1, SCALE=FNT/REGYR=1, SCALE2=YR/FNT=1, GMOD=1
# (balmod.f), DGMX=DGMAX·SCALE=5.0):
#   1. HTMAX gate (htcalc.f MODE0=0 = `_on_htcalc_htmax`): if HTMAX-H ≤ 1 → HTGR=0.1.
#   2. else HTGR = ONSTHG(sp,D,H,SI,BAL)·10 (BAL = BA·(100-PCT)·0.01); floor 0.1.
#   3. XWT = (D-XMIN)/(XMAX-XMIN) (0 if D≤XMIN); HTGR = HTGR·(1-XWT) + XWT·HTG_large; floor 0.1.
#   4. random ±0.1·HTGR (BACHLO in [-1,1], active DGSD≥1 — ON DGSD=2.0). SIZCAP(sp,4) height cap.
#   5. DBH incr: HK=H+HTG. HK≤4.5 → DG=0, DBH=D+0.001·HK. Else Wykoff/Curtis-Arney HTDBH inverse
#      (DKK from HK, DK from H) → DGSM=(DKK-DK)·BARK → DDS blend → DGSM=√((D·BARK)²+DDS)-BARK·D →
#      DGGR = DGSM·(1-XWT) + XWT·DG_large; floor 0.1, cap DGMX, DIAM budwidth floor, DGBND.
#
# Coefficients auto-extracted from the Fortran DATA into `_regent_data.jl` (onsthg.f OSPMAP+B00..B95,
# blkdat.f HT1/HT2, htdbh.f IWYKCA/SNDBAL/SNALL, rgntsw.f DIAM). Transcendentals via glibc
# logf/expf/powf (on_logf/on_expf/on_powf) for Float32 bit-exactness.
#
# VALIDATION (vs instrumented FVSon_g16 on the ont_sm small-tree stand, scratchpad/on/): **ONSTHG is
# bit-exact** — per-tree HTG(annual) hex matches for both the exp(-5)-clamped dense case and the
# computed case; and the deterministic HTGR matches bit-exact on records where the [XMIN,XMAX] blend
# does not consume the large-tree HTG (D≤XMIN ⇒ xwt=0). On records with xwt>0 the blended HTGR
# inherits a ~1.6% residual from the large-tree HTG/DG for the SAME sub-12cm records — and that is the
# accepted **#206 OLDRN straddle, NOT a bug**: the large-tree Penner DDS/DIAGR are BIT-EXACT for these
# small trees (on_penner_dds fed the oracle's inputs → DDS+DIAGR hex match), so the large-tree DG
# (jl ~0.0249 vs oracle ~0.0233, consistent-sign) differs only by the DGSD=2.0 OLDRN serial-corr draw
# `DG=√(D²+exp(WK2+OLDRN)·SCALE)−D` — the same cornered class as every western variant's cyc1. ont_sm
# is a degenerate super-dense synthetic stand (BA >1000 m²/ha) built to exercise the branch; ont01
# (all-large) → REGENT no-op, so multicycle 339/11 + ontario suite are byte-identical. Recipe:
# scratchpad/on/REGENT_HANDOFF.md.
# =============================================================================

include("_regent_data.jl")

const ON_REG_XMAX  = 4.72f0   # regent.f:78 XMAX = 12 cm (uniform across all 72 species)
const ON_REG_XMIN  = 3.15f0   # regent.f:78 XMIN = 8 cm  (uniform)
const ON_REG_REGYR = 10.0f0   # regent.f:79
const ON_REG_DGMAX = 5.0f0    # regent.f:76 DGMAX = MAXSP*5.0

"ONSTHG (canada/on/onsthg.f): Penner-2010 annual small-tree height increment (ft)."
@inline function _on_onsthg(sp::Int, d::Float32, h::Float32, si::Float32, bal::Float32)::Float32
    ksp  = ON_ST_OSPMAP[sp]
    hm   = max(0.05f0, h * ON_FTtoM)
    lhm  = on_logf(hm)
    sim  = si * ON_FTtoM
    balm = bal * ON_FT2pACRtoM2pHA
    htg = ON_ST_B00[ksp] + ON_ST_B01[ksp]*lhm + ON_ST_B02[ksp]*hm +
          ON_ST_BSI[ksp]*sim + ON_ST_BBAL[ksp]*balm
    htg = min(max(htg, -5f0), 5f0)
    htg = on_expf(htg)
    htg = min(max(htg, 0.0001f0), ON_ST_B95[ksp])
    return htg * ON_MtoFT
end

"HTDBH inverse (canada/on/htdbh.f MODE=1, H→D): Wykoff or Curtis-Arney, floored at the SNDBAL budwidth."
@inline function _on_htdbh_dbh(sp::Int, h::Float32)::Float32
    local d::Float32
    if ON_REG_IWYKCA[sp] == 0                                  # Wykoff (htdbh.f:344)
        d = (ON_REG_HT2[sp] / (on_logf(h - 4.5f0) - ON_REG_HT1[sp])) - 1f0
    else                                                       # Curtis-Arney (htdbh.f:345-352)
        p2 = ON_REG_SNALL_P2[sp]; p3 = ON_REG_SNALL_P3[sp]; p4 = ON_REG_SNALL_P4[sp]
        hat3 = 4.5f0 + p2 * on_expf(-p3 * on_powf(3f0, p4))
        if h >= hat3
            d = on_expf(on_logf((on_logf(h - 4.5f0) - on_logf(p2)) / (-p3)) * (1f0 / p4))
        else
            db = ON_REG_SNDBAL[sp]
            d = (((h - 4.51f0) * (3f0 - db)) / (4.5f0 + p2*on_expf(-p3*on_powf(3f0, p4)) - 4.51f0)) + db
        end
    end
    d < ON_REG_SNDBAL[sp] && (d = ON_REG_SNDBAL[sp])           # htdbh.f:356 IF(D<DB)D=DB
    return d
end

"HTDBH forward (canada/on/htdbh.f MODE=0, D→H): Wykoff or Curtis-Arney. Used to dub missing/broken-top
heights (cratet.f, LHTDRG=.FALSE. all ON species ⇒ always the default-coefficient curve). No DB floor
(htdbh.f:356 applies it only for MODE≠0)."
@inline function _on_htdbh_height(sp::Int, d::Float32)::Float32
    if ON_REG_IWYKCA[sp] == 0                                  # Wykoff (htdbh.f:334)
        return on_expf(ON_REG_HT1[sp] + ON_REG_HT2[sp] / (d + 1f0)) + 4.5f0
    else                                                       # Curtis-Arney (htdbh.f:336-341)
        p2 = ON_REG_SNALL_P2[sp]; p3 = ON_REG_SNALL_P3[sp]
        p4 = ON_REG_SNALL_P4[sp]; db = ON_REG_SNDBAL[sp]
        if d >= 3f0
            return 4.5f0 + p2 * on_expf(-p3 * on_powf(d, p4))
        else
            return ((4.5f0 + p2 * on_expf(-p3 * on_powf(3f0, p4)) - 4.51f0) * (d - db) / (3f0 - db)) + 4.51f0
        end
    end
end

"DGBND (canada/on/dgbnd.f): cap DG at 6·exp(-0.03·DBH) (DBH≤150), floor 0."
@inline function _on_dgbnd(dbh::Float32, dg::Float32)::Float32
    d = dbh > 150f0 ? 150f0 : dbh
    m = 6f0 * on_expf(-0.03f0 * d)
    dg > m && (dg = m)
    dg < 0f0 && (dg = 0f0)
    return dg
end

"""
    small_tree_growth!(s, stash, ::Ontario; fint) — ON REGENT (canada/on/regent.f).

Overrides DGF/HTGF for records with DBH < XMAX (4.72" = 12 cm) with the Penner small-tree height
increment (ONSTHG) blended with the large-tree HTG over [XMIN,XMAX], then dubbed to a diameter
increment. Records ≥ XMAX are untouched (validated no-op). Runs after `height_growth!`.
"""
function small_tree_growth!(s::StandState, stash, ::Ontario; fint::Float32 = 10.0f0)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    sizcap = s.control.sp_size_cap
    ba = p.basal_area
    species_sort!(s)
    trip = stash !== nothing
    nrec = trip ? 3 : 1
    scale  = fint / ON_REG_REGYR                 # FNT/REGYR = 1 at fint=10
    scale2 = ON_REG_REGYR / fint                 # YR/FNT   = 1 at fint=10
    dgmx = ON_REG_DGMAX * scale                  # 5.0
    random_on = s.control.dg_sd >= 1f0           # DGSD (ON=2.0) ⇒ ±10% random effect
    @inbounds for sp in 1:MAXSP
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        si = p.sp_site_index[sp]
        htmax = _on_htcalc_htmax(sp, si)
        for k3 in i1:i2
            i = ind1[k3]
            d = t.dbh[i]
            (d >= ON_REG_XMAX || t.tpa[i] <= 0f0) && continue
            h = t.height[i]
            pct = t.crown_ratio[i]                       # FVS PCT (BA percentile in larger trees)
            bal = ba * (100f0 - pct) * 0.01f0
            if htmax - h <= 1f0
                htgr_s = 0.1f0
            else
                htgr_s = _on_onsthg(sp, d, h, si, bal) * 10f0 * scale   # ONSTHG·YRS·(con·HGADJ·xrhgro=1)
                htgr_s < 0.1f0 && (htgr_s = 0.1f0)
            end
            xwt = d <= ON_REG_XMIN ? 0f0 : (d - ON_REG_XMIN) / (ON_REG_XMAX - ON_REG_XMIN)
            htgr = htgr_s * (1f0 - xwt) + xwt * t.ht_growth[i]
            htgr < 0.1f0 && (htgr = 0.1f0)
            dg_large = t.diam_growth[i]                  # DG(K) large-tree, for the DGSM blend
            for l in 0:(nrec - 1)
                ran = 0f0
                if random_on
                    while true
                        ran = bachlo(s.rng, 0f0, 1f0)
                        (-1f0 <= ran <= 1f0) && break
                    end
                end
                htg = htgr + ran * 0.1f0 * htgr
                htg < 0.1f0 && (htg = 0.1f0)
                (h + htg) > sizcap[sp, 4] && (htg = max(sizcap[sp, 4] - h, 0.1f0))
                hk = h + htg
                if hk <= 4.5f0
                    dg = 0.001f0 * hk                    # regent.f:296 DBH=D+0.001·HK (DG=0); no DIAM/DGBND
                else
                    dkk = _on_htdbh_dbh(sp, hk)
                    dk  = h <= 4.5f0 ? d : _on_htdbh_dbh(sp, h)
                    bark = on_bratio(sp, d, h)
                    if dk < 0f0 || dkk < 0f0
                        dg = htg * 0.2f0 * bark          # regent.f:362 (·XRDGRO=1)
                    else
                        dgsm = (dkk - dk) * bark
                        dgsm < 0f0 && (dgsm = 0f0)
                        dds = dgsm * (2f0 * bark * d + dgsm) * scale2
                        dgsm = sqrt((d * bark)^2 + dds) - bark * d
                        dgsm < 0f0 && (dgsm = 0f0)
                        dg = dgsm * (1f0 - xwt) + xwt * dg_large   # DGGR blend
                        dg < 0.1f0 && (dg = 0.1f0)
                    end
                    dg < 0f0 && (dg = 0.1f0)             # regent.f:382
                    dg > dgmx && (dg = dgmx)             # regent.f:383
                    (d + dg) < ON_REG_DIAM[sp] && (dg = ON_REG_DIAM[sp] - d)   # DIAM budwidth floor
                    dg = _on_dgbnd(d, dg)               # DGBND
                end
                if l == 0
                    t.diam_growth[i] = dg; t.ht_growth[i] = htg
                elseif l == 1
                    stash.dgU[i] = dg; stash.htgU[i] = htg; stash.is_small[i] = true
                else
                    stash.dgL[i] = dg; stash.htgL[i] = htg
                end
            end
        end
    end
    return s
end
