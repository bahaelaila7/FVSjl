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
