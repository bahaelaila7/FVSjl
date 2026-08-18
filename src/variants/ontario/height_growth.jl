# =============================================================================
# ontario/height_growth.jl — ON large-tree height growth (canada/on/htgf.f + htont.f).
#
# CHUNK 7. `height_growth!(s, ::Ontario; scale)` fills trees.ht_growth per live tree, mirroring
# htgf.f. ON's height increment is a DIAMETER-HEIGHT DIFFERENCE, not a rate:
#     HTNOW = HTONT(sp, DBH,   RMSQD, BA)          ! Penner height at current DBH / stand density
#     HT10  = HTONT(sp, DBH10, QMD10, BA10)        ! Penner height at grown DBH / 10-yr density
#     HTG   = max(HT10 - HTNOW, 0.1) · SCALE · XHMULT · HCOR2    (BALMOD GM≡1.0 for ON)
# where DBH10 = DBH + DG/BRATIO(sp,DBH,HT) and the 10-yr stand density is rebuilt in RECORD order
# (htgf.f DO 20): BA10 = Σ 0.0054542·PROB·DBH10², QMD10 = √(Σ PROB·DBH10² / Σ PROB), clamped BA10≥BA.
#
# HTONT (canada/on/htont.f) is the Penner (2006) diameter-height model, 27 equations, OSPMAP_P maps
# the 72 FVS species → 27 height species. All metric:
#     HTM(m) = 1.3 + (B0 + BSI·SIM + BDBHQ·DBHQM + BBA·BAM)·(1 − exp(−B1·DBHM^B2)),  HT = HTM·MtoFT
#     DBHM = DBH·INtoCM,  DBHQM = QMD·INtoCM,  BAM = BA·FT2pACRtoM2pHA,  SIM = SITEAR(ISISP)·FTtoM.
# SIM uses the STAND site species (ISISP = p.site_species), constant across trees — NOT per-tree.
# The Sharma–Parton branch (htont.f OSPMAP_SP) is all-zero in the shipped variant ⇒ never used, so it
# is intentionally not ported. All 72 OSPMAP_P entries are >0 ⇒ the "don't grow" (+0.001) branch is
# unreachable on any real ON species; it is handled defensively (returns 0.001 ft, the -finit-local-zero
# value FVS would produce) for a species outside the map.
#
# Transcendentals MUST be glibc single-precision expf/powf (via on_expf/on_powf in diameter_growth.jl):
# native Julia exp/^ drift ~1 ULP; measured glibc = 0-mismatch. VALIDATED bit-exact (Float32-hex) vs
# instrumented FVSon_g16 on ont01: HTNOW 8/8, HT10 8/8, final HTG 8/8, BA10 + QMD10 bit-exact
# (scratchpad/on/replay_htont.jl vs /workspace/.onwork/FVSon_htgdump fort.773).
# =============================================================================

const ON_MtoFT = 3.28084f0   # common/METRIC.F77 (metres → feet)

# --- HTONT Penner diameter-height coefficients (27 eqns), canada/on/htont.f DATA blocks ---
const ON_HT_OSPMAP = Int[
    5, 4, 3, 4, 1, 7, 8,12, 9,14,
   13,11, 9, 9,21,22,26,16,16,25,
   16,16,16,17,23,15,15,20,22,20,
   16,20,15,19,21,21,25,25,16,26,
   26,26,18,15,16,16,24,16,18,16,
   16,16,24,16,16,16,16,16,16,16,
   16,16,16,16,16,16,16,16, 6, 2,
    8,10]
const ON_HT_B0 = Float32[
   18.1847, 5.3531,16.8042,15.8338, 4.8752,
    1.7254,22.042, 14.7619,16.8201, 5.6734,
   17.6292,15.585, 11.9821, 9.1537,17.6005,
   18.0739,14.0772,12.689, 16.6299,15.5719,
   14.6377,18.451, 20.4829,13.8802,17.3966,
   14.7273,15.6942]
const ON_HT_B1 = Float32[
    0.0354, 0.0622, 0.032,  0.00746,0.098,
    0.1515, 0.0175, 0.0264, 0.0277, 0.0572,
    0.0428, 0.0224, 0.038,  0.058,  0.0703,
    0.0724, 0.0689, 0.0589, 0.0658, 0.0313,
    0.0705, 0.0718, 0.0417, 0.0926, 0.0402,
    0.0619, 0.0641]
const ON_HT_B2 = Float32[
    1.044,  1.0429, 1.1956, 1.7776, 0.9494,
    0.7898, 1.2469, 1.0967, 1.2443, 0.9494,
    0.9593, 1.3855, 1.0663, 1.138,  1.0014,
    0.9472, 0.9642, 1.135,  0.8726, 1.2336,
    0.9744, 0.9302, 1.1079, 1.1378, 1.1548,
    1.0194, 0.9536]
const ON_HT_BSI = Float32[
    0.4402, 0.1260, 0.2008, 0.2322, 0.,
    0.,     0.1009, 0.0374, 0.,     0.1163,
    0.2185, 0.04,   0.2047, 0.0602, 0.1437,
    0.1274, 0.0819, 0.14,   0.3082, 0.2515,
    0.3886, 0.3932, 0.1059, 0.0842, 0.1733,
    0.1138, 0.1659]
const ON_HT_BDBHQ = Float32[
    0.,     0.6407, 0.,     0.,     0.5207,
    0.6616, 0.,     0.,     0.1608, 0.5513,
    0.,     0.,     0.,     0.4227, 0.,
    0.,     0.0868, 0.1616, 0.,     0.,
    0.,     0.,     0.1319, 0.,     0.,
    0.3114, 0.1951]
const ON_HT_BBA = Float32[
    0.1191, 0.,     0.1271, 0.,     0.2505,
    0.2669, 0.0997, 0.2879, 0.1117, 0.1014,
    0.1159, 0.1147, 0.0794, 0.1285, 0.0884,
    0.0778, 0.1176, 0.06,   0.1219, 0.1312,
    0.048,  0.,     0.,     0.,     0.0815,
    0.0919, 0.0678]

"""
    on_htont(ispc, d_in, q_in, ba_ft, sim) -> Float32

Ontario Penner diameter-height model (canada/on/htont.f). Returns the predicted total height (ft)
for a tree of species `ispc` (1-based FVS index) at diameter `d_in` (in), given the stand quadratic
mean diameter `q_in` (in), basal area `ba_ft` (ft²/ac) and the pre-computed site term
`sim = SITEAR(ISISP)·FTtoM` (m). Bit-exact vs FVSon_g16.
"""
@inline function on_htont(ispc::Int, d_in::Float32, q_in::Float32, ba_ft::Float32, sim::Float32)::Float32
    ksp = (ispc >= 1 && ispc <= length(ON_HT_OSPMAP)) ? ON_HT_OSPMAP[ispc] : 0
    ksp == 0 && return 0.001f0            # htont.f "don't grow" (+0.001); unreachable on real ON species
    dbhm  = d_in * ON_INtoCM
    bam   = ba_ft * ON_FT2pACRtoM2pHA
    dbhqm = q_in * ON_INtoCM
    htm = 1.3f0 + (ON_HT_B0[ksp] + ON_HT_BSI[ksp]*sim + ON_HT_BDBHQ[ksp]*dbhqm + ON_HT_BBA[ksp]*bam) *
          (1.0f0 - on_expf(-ON_HT_B1[ksp] * on_powf(dbhm, ON_HT_B2[ksp])))
    return htm * ON_MtoFT
end

"""
    height_growth!(s, ::Ontario; scale) — per-tree periodic height increment → trees.ht_growth.

Faithful transcription of canada/on/htgf.f (see file header). `scale = FINT/YR` (the engine passes
`fint/htg_period` = 10/10 = 1 for ON's 10-yr cycle). BALMOD ≡ 1.0 for ON (canada/on/balmod.f).
"""
function height_growth!(s::StandState, ::Ontario; scale::Float32 = 1.0f0)
    p, t, ctl = s.plot, s.trees, s.control
    ba    = p.basal_area
    rmsqd = p.qmd
    isisp = Int(p.site_species)
    sitear = (isisp >= 1 && isisp <= length(p.sp_site_index)) ? p.sp_site_index[isisp] : 0f0
    sim = sitear * ON_FTtoM

    # --- DO 20: rebuild the 10-yr stand density in RECORD order (BA10/QMD10). ---
    ba10 = 0f0; sumdbhsq = 0f0; xtrees = 0f0
    @inbounds for i in 1:t.n
        sp = Int(t.species[i])
        bark  = on_bratio(sp, t.dbh[i], t.height[i])          # BRATIO(ISP,DBH,HT), current dims
        dbh10 = t.dbh[i] + t.diam_growth[i] / bark
        x = t.tpa[i] * dbh10 * dbh10
        ba10     += 0.0054542f0 * x
        xtrees   += t.tpa[i]
        sumdbhsq += x
    end
    qmd10 = xtrees > 0f0 ? sqrt(sumdbhsq / xtrees) : 0f0
    ba10 < ba && (ba10 = ba)

    cur_year = current_cycle_year(s)
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0f0
        sp = Int(t.species[i])
        bark  = on_bratio(sp, t.dbh[i], t.height[i])
        dbh10 = t.dbh[i] + t.diam_growth[i] / bark
        htnow = on_htont(sp, t.dbh[i], rmsqd, ba,   sim)
        ht10  = on_htont(sp, dbh10,    qmd10, ba10, sim)
        htg = ht10 - htnow
        htg < 0.1f0 && (htg = 0.1f0)                          # BALMOD GM=1.0 folded in
        xht  = active_multiplier(ctl, :htg, sp, cur_year)     # XHMULT (HTGMULT keyword; 1 default)
        xht2 = ctl.htg_cor2[sp]                               # HCOR2 (READCORH; 1 default)
        htg = scale * xht * htg * xht2
        cap = ctl.sp_size_cap[sp, 4]                          # SIZCAP(sp,4) height cap (999 default)
        if t.height[i] + htg > cap
            htg = cap - t.height[i]
            htg < 0.1f0 && (htg = 0.1f0)
        end
        t.ht_growth[i] = htg
    end
    return s
end
