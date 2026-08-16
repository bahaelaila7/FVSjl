# =============================================================================
# small_tree_growth.jl (olympic) — OP small-tree REGENT growth (op/regent.f + op/smhgdg.f).
#
# op/regent.f grows trees with DBH < XMAX(ISPC) (≈4"): SMHGDG gives the 5-yr small-tree height (HG5)
# and diameter (DG5) increments, then REGENT blends the small-tree height with the large-tree HTG via
# XWT=(D−XMN)/(XMX−XMN) (XWT=0 for D≤XMN pure small-tree; XWT=1 for D≥XMX or a valid-ORGANON tree ⇒
# large-tree unchanged), and for D<DGMIN(ISPC) overrides the diameter growth:
#   • HK=H+HTG < 4.5 ⇒ DG=0, DBH := D + 0.001·HK  (set directly; the shared apply-loop adds DG=0).
#   • HK ≥ 4.5      ⇒ DG := sqrt((D·BARK)²+DDS)−D·BARK with DDS=(DG5·BARK)·(2·BARK·D+DG5·BARK)·SCALE2.
# Runs AFTER diameter_growth!/height_growth! (gradd.f DGDRIV→HTGF→REGENT order). op is a 5-yr step
# (REGYR=5, SCALE=SCALE2=1); no small-tree calibration on S248112 (oracle small-tree scale = 1.00 ⇒
# CON=1). Deterministic (DGSD=0 ⇒ no ZZRAN/BACHLO draw). Validated per-tree vs live FVSop_clean.
# =============================================================================

# op/smhgdg.f DATA DMAX/BETA (39 species) + ALPHA(0:9,39).
const OP_SMH_DMAX = Float32[
    1.704,1.496,1.639,1.196,1.515,3.396,2.939,1.540,1.683,1.885,
    1.653,1.798,2.474,2.474,1.798,5.373,2.849,2.790,3.419,1.383,
    3.094,3.094,2.011,2.166,3.094,2.475,3.713,0.986,1.219,0.623,
    0.807,0.586,0.860,1.003,1.890,2.166,2.166,0.000,5.373]
const OP_SMH_BETA = Float32[
    0.247400,0.217500,0.179700,0.205600,0.216800,0.216800,0.282200,0.216800,0.281500,0.170400,
    0.168200,0.216800,0.216800,0.216800,0.236900,0.163500,0.172700,0.182900,0.172700,0.302900,
    0.216800,0.216800,0.216800,0.216800,0.216800,0.216800,0.216800,0.216800,0.216800,0.216800,
    0.216800,0.168200,0.216800,0.216800,0.216800,0.216800,0.216800,0.000000,0.163500]
# ALPHA[k+1, sp] for k=0..9 (10 rows × 39 species).
const OP_SMH_ALPHA = Float32[
  # a0
  2.94450 1.75360 2.35710 2.58390 2.47430 3.82050 0.33760 -2.02160 0.59960 0.04520 1.74000 1.84510 3.80850 3.80850 1.84510 2.44730 2.95270 1.68150 2.95270 2.67620 -1.24210 1.45930 -1.19000 -1.24210 -1.24210 -1.24210 -1.24210 -2.19100 0.37550 1.05270 2.49490 -0.80850 1.51560 -3.83450 3.55210 -1.24210 -1.24210 0.00000 2.44730;
  # a1
  0.00000 0.00000 0.00520 0.00000 0.00000 0.00000 0.00000 0.00630 0.00000 0.00800 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00240 0.01240 0.00000 0.01580 0.01240 0.01240 0.01240 0.01240 0.00000 0.01200 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.01240 0.01240 0.00000 0.00000;
  # a2
  0.00000 0.29280 0.00000 0.04100 0.00000 0.05230 0.00000 0.00000 0.00000 0.00000 0.37180 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.35800 0.00000 0.50010 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000;
  # a3
  0.00680 0.00090 0.00060 0.00200 0.00320 0.00510 0.01010 0.00000 0.00800 0.00710 0.00270 0.01670 0.00230 0.00230 0.01670 0.00980 0.00660 0.00680 0.00660 0.00060 0.00000 0.00850 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00190 0.00490 0.00000 0.00120 0.00000 0.00020 0.00000 0.00000 0.00000 0.00980;
  # a4
  0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.71750 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.41610 0.00000 0.66000 0.78130 0.63820 0.60130 0.60130 0.71910 0.00000 0.00000 0.00000 0.00000 0.00000 1.07010 0.00000 0.73120 0.65980 0.00000 0.00000;
  # a5
  0.00000 -0.04460 -0.42690 -0.01520 -0.89340 -0.41020 0.00000 0.00000 0.00000 0.00000 -0.17120 -1.47370 -0.42650 -0.42650 -1.47370 -0.42900 0.00000 0.00000 0.00000 -0.43090 0.00000 -0.60000 0.00000 0.00000 0.00000 0.00000 0.00000 -3.13210 0.00000 0.00000 -0.20850 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 -0.35750;
  # a6
  -0.18950 -2.03490 -1.22190 -2.20600 -2.27090 -1.69680 0.00000 0.00000 0.00000 0.00000 -2.13590 0.00000 -2.09130 -2.09130 0.00000 -0.17100 -0.47340 0.00000 -0.47340 -1.62050 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 -0.60080 -1.70010 0.00000 -0.54780 0.00000 -0.59320 0.00000 0.00000 0.00000 -0.17100;
  # a7
  0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 -1.04790 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000;
  # a8
  -1.40490 -1.38390 0.00000 -0.59150 -1.06900 -1.40010 0.00000 0.00000 0.00000 0.00000 -0.72660 -0.41030 -1.39320 -1.39320 -0.41030 -0.18790 -0.73940 -0.60490 -0.73940 -0.59300 0.00000 -1.22800 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 -0.74510 -0.79520 0.00000 -0.61230 0.00000 -0.50290 0.00000 0.00000 0.00000 -0.18790;
  # a9
  -0.01680 -0.00330 -0.01700 -0.00090 0.00000 -0.01090 -0.00430 0.00000 0.00000 0.00000 -0.00740 -0.01120 -0.00930 -0.00930 -0.01120 -0.01100 -0.02070 -0.01210 -0.02070 -0.00510 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 -0.01010 -0.01770 -0.00810 0.00000 0.00000 -0.00380 0.00000 0.00000 0.00000 -0.01100]

# op/regent.f DATA (39 species).
const OP_REG_XMAX = Float32[4,4,4,4,4,4,4,4,4,4,3,4,4,4,4,4,10,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4]
const OP_REG_XMIN = Float32[2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
const OP_REG_DGMIN = Float32[3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,7,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3]
const OP_REG_DIAM = Float32[0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.2,0.2,0.3,0.4,0.4,0.4,0.4,0.5,0.3,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.3,0.4,0.4,0.2,0.2,0.2,0.2,0.2,0.2,0.2]
const OP_REG_DGMAX = 5.0f0
const OP_REG_REGYR = 5.0f0

"op/smhgdg.f — 5-yr small-tree (HG5, DG5) for FVS species `ispc` (1..39), NON-redwood (RW handled in regent)."
@inline function op_smhgdg(ispc::Int, h::Float32, d::Float32, cr::Float32, ptbal::Float32,
                           ptba::Float32, avht::Float32, si::Float32)
    relht = avht > 0f0 ? h / avht : 0f0
    relht > 1.5f0 && (relht = 1.5f0)
    ptbal2 = log(ptbal + 2.71f0); ptba2 = log(ptba + 2.71f0); relht2 = sqrt(relht)
    boost = ptba < 100f0 ? 1f0 / (1f0 + exp(-3.1f0 + 0.18f0 * ptba)) : 0f0
    beta = OP_SMH_BETA[ispc]
    z = OP_SMH_ALPHA[1, ispc] + OP_SMH_ALPHA[2, ispc] * ptba + OP_SMH_ALPHA[3, ispc] * ptba2 +
        OP_SMH_ALPHA[4, ispc] * ptbal + OP_SMH_ALPHA[5, ispc] * ptbal2 + OP_SMH_ALPHA[6, ispc] * boost +
        OP_SMH_ALPHA[7, ispc] * cr + OP_SMH_ALPHA[8, ispc] * relht + OP_SMH_ALPHA[9, ispc] * relht2 +
        OP_SMH_ALPHA[10, ispc] * si
    dgs = OP_SMH_DMAX[ispc] / (1f0 + exp(z))
    hg5 = dgs / beta
    local dg5::Float32
    if h < 4.5f0
        dg5 = (h + hg5) >= 4.5f0 ? beta * (h + hg5 - 4.5f0) : 0f0
    else
        dg5 = dgs
    end
    return hg5, dg5
end

function small_tree_growth!(s::StandState, stash, ::Olympic; fint::Float32 = 5.0f0)
    p, t, c = s.plot, s.trees, s.calib
    dens = s.density
    n = t.n; n == 0 && return s
    avh = p.avg_height
    org_ran = length(c.op_iorg) == n
    scale = fint / OP_REG_REGYR          # SCALE = FNT/REGYR (FNT=FINT on the non-estab path)
    scale2 = htg_period(s.variant) / fint   # SCALE2 = YR/FNT (YR=5)
    dgmx = OP_REG_DGMAX * scale
    sizcap = s.control.sp_size_cap
    @inbounds for i in 1:n
        t.tpa[i] <= 0f0 && continue
        ispc = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (ispc < 1 || ispc > 39) && continue
        d >= OP_REG_XMAX[ispc] && continue          # op/regent.f:202 D≥XMX ⇒ pure large-tree (skip)
        ispc == 17 && continue                       # RW small-tree handled separately (no RW in S248112)
        iorg = org_ran && c.op_iorg[i] == 1
        # --- SMHGDG small-tree height/diameter (op/smhgdg.f) ---
        cr = Float32(t.crown_pct[i]) * 0.01f0
        pt = Int(t.plot_id[i])
        ptbal = i <= length(dens.point_bal) ? dens.point_bal[i] : 0f0   # PTBALT(IT) — per-TREE BA-in-larger
        ptba  = (1 <= pt <= length(dens.point_ba)) ? dens.point_ba[pt] : 0f0   # PTBAA(ITRE) — per-point BA
        si = ispc <= length(p.sp_site_index) ? p.sp_site_index[ispc] : 0f0
        hg5, dg5 = op_smhgdg(ispc, h, d, cr, ptbal, ptba, avh, si)
        con = exp(c.htg_cor_small[ispc])            # RHCON=1 (no READCORR); HCOR small-tree calib (0 ⇒ CON=1)
        htgr = hg5 * scale * con                     # XRHGRO=1, WK4=1, ZZRAN=0 (DGSD<1)
        htgr < 0.1f0 && (htgr = 0.1f0)
        # --- HTG blend with the large-tree HTG (op/regent.f:282-304) ---
        xmn = OP_REG_XMIN[ispc]; xmx = OP_REG_XMAX[ispc]
        xwt = (d - xmn) / (xmx - xmn)
        (d <= xmn) && (xwt = 0f0)
        iorg && (xwt = 1f0)                          # valid-ORGANON tree ⇒ keep the large-tree HTG
        htg = htgr * (1f0 - xwt) + xwt * t.ht_growth[i]
        htg < 0.1f0 && (htg = 0.1f0)
        cap4 = sizcap[ispc, 4]
        (cap4 > 0f0 && (h + htg) > cap4) && (htg = max(cap4 - h, 0.1f0))
        t.ht_growth[i] = htg
        # --- diameter (op/regent.f:325-454): only for D<DGMIN and NOT a valid-ORGANON tree ---
        (iorg || d >= OP_REG_DGMIN[ispc]) && continue
        hk = h + htg
        if hk < 4.5f0
            t.diam_growth[i] = 0f0
            t.dbh[i] = d + 0.001f0 * hk              # DBH set directly (shared loop adds DG=0)
        else
            bark = op_bratio(ispc, d)
            dg = dg5 * bark                          # DG5·BARK·XRDGRO (XRDGRO=1)
            dg < 0f0 && (dg = 0.1f0)
            dg > dgmx && (dg = dgmx)
            dds = dg * (2f0 * bark * d + dg) * scale2
            dg = sqrt((d * bark)^2 + dds) - bark * d
            (d + dg) < OP_REG_DIAM[ispc] && (dg = OP_REG_DIAM[ispc] - d)   # regent.f:452 (DBH(K)=D here)
            t.diam_growth[i] = dg
        end
    end
    return s
end
