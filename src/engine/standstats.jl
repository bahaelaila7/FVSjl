# =============================================================================
# standstats.jl — per-acre expansion (NOTRE) and stand summary statistics
#
# Ported from: base/notre.f (expansion) + base/stats.f / sumout (summary columns).
#
# NOTRE turns each record's tally (PROB) into trees-per-acre using the sampling
# design (variable-radius BAF plots vs fixed plots). The stand statistics
# (TPA, basal area, QMD, ...) are then simple weighted reductions over the trees —
# pure loops, autovectorizable, no allocation.
# =============================================================================

const BA_PER_TREE = 0.005454154f0     # ft² of basal area per inch² of DBH

"""
    finalize_design!(state)

INITRE end-of-keywords finalization (initre.f:360-414). The crucial, easy-to-miss
step: FVS overwrites the variable named `PI` with the plot count `IPTINV` — that is
the divisor NOTRE actually uses (NOT π). Also finalizes the stockable proportion
GROSPC into the reciprocal multiplier NOTRE applies.
"""
function finalize_design!(s::StandState)
    p = s.plot
    p.points_inv <= 0 && (p.points_inv = Int32(1))
    p.pi = Float32(p.points_inv)                      # PI := IPTINV (NOTRE divisor!)
    p.sample_weight < 0f0 && (p.sample_weight = Float32(p.points_inv))
    if p.gross_space < 0f0                             # default (set to -1 at stand start)
        g = (p.pi - Float32(p.nonstockable)) / p.pi
        g > 1f0 && (g = 1f0)
        (p.pi - Float32(p.nonstockable)) <= 0f0 && (g = 1f0)
        p.gross_space = g
    end
    p.gross_space = 1f0 / p.gross_space               # → reciprocal multiplier
    return s
end

"""
    notre!(state)

Compute trees-per-acre (`trees.tpa`) for every live record from the sampling
design. Mirrors notre.f: variable-plot trees (DBH ≥ BRK) scale by BAF/DBH²,
fixed-plot trees (DBH < BRK) by FPA, all times the stockable proportion GROSPC.
"""
function notre!(s::StandState)
    p, t = s.plot, s.trees
    fp  = p.fixed_plot_inv / p.pi
    p.total_fixed_plot > 0f0 && (fp = 1f0 / p.total_fixed_plot)
    vp  = p.baf * 183.3465f0 / p.pi
    fp2 = 0f0
    if p.baf <= 0f0
        vp = 0f0; fp2 = -p.baf / p.pi
    end
    brk = p.min_dbh_var_plot
    # expand live records and the dead partition (n+1:n+ndead) alike — dead trees
    # carry their expanded TPA into the backdated calibration BA.
    @inbounds for i in 1:(t.n + t.ndead)
        pr = t.tpa[i]; d = t.dbh[i]
        pr <= 0f0 && (pr = 1f0)
        pr = d < brk ? pr * fp : pr * vp / (d * d) + pr * fp2
        pr <= 0f0 && (pr = 9.0f-25)
        t.tpa[i] = pr * p.gross_space
    end
    # base species-sort key = original record index (FVS chain order pre-tripling);
    # tripling derives child keys so the per-tree RNG draw order matches the oracle.
    @inbounds for i in 1:(t.n + t.ndead)
        t.sort_key[i] = Float64(i)
    end
    return s
end

"""
    stand_tpa(state)   -> total trees per acre        (TPROB)
    stand_ba(state)    -> basal area ft²/acre         (BA)
    stand_qmd(state)   -> quadratic mean diameter, in (RMSQD)

Weighted reductions over the live trees. Pure (no mutation).
"""
function stand_tpa(s::StandState)
    t = s.trees; tot = 0f0
    @inbounds for i in 1:t.n; tot += t.tpa[i]; end
    return tot
end

# dense.f (identical in every variant build): inside DO 50 ISPC / DO 10 I3 / I=IND1(I3):
#   DP=D*P; WK5(I)=D*DP; TSUMD2=TSUMD2+WK5(I); BATREE=0.005454154*WK5(I); BAT=BAT+BATREE; TPROB=TPROB+P
# BA=BAT, RMSQD=SQRT(TSUMD2/TPROB). Both the D*(D*P) association and the IND1 accumulation order are part of
# the REAL*4 result (BM 30193202010497 cyc1 BA 4261C4DB live vs 4261C4DC for record-order p·K·d²).
# DENSE walk order for BA/RMSQD: FVS IND1. jl's IND1 reconstruction (`_ind1_order`, sort-key lineage) is
# live-measured faithful for BM; on CS kwcov cs_serlcorr it is NOT (IND1 order flips a 2030 BdFt cell vs live while
# the D*(D*P) association alone is inert) ⇒ the CS lineage keys diverge from FVS's LNKCHN chain there (open lead).
# Until that lineage is fixed per variant, only BM walks IND1 here; others keep record order. The D*(D*P)
# association is likewise BM-only: applied to SN it moved test_growth COR 1 ULP off Oracle A and one SN
# test_allspecies cell off live (19765 vs 19766) — SN dense.f is a different source revision (open lead).
# ON (canada/on builds base/dense.f): TSUMD2/BAT in IND1 order with the D*(D*P) association — live-measured on the
# ON fixtures (FVSon_g16 DENSE dump): RMSQD ont_all cyc1 411F20D4 / ont_lite 41266DCA and BA ont_sm cyc0 441BF990
# only with both; record order or p·d² is 1-4 ULP off (ont01 SB HTONT(RMSQD) HTNOW 1 ULP ⇒ HtG).
_dense_order(s::StandState) = (s.variant isa BlueMountains || s.variant isa Ontario) ? _ind1_order(s) : (1:s.trees.n)

function stand_ba(s::StandState)
    t = s.trees; ba = 0f0
    if s.variant isa BlueMountains || s.variant isa InlandEmpire || s.variant isa Ontario
        # dense.f:179-190 — species-major IND1 order, DP=D·P; WK5=D·DP; BATREE=0.005454154·WK5; BAT=BAT+BATREE
        # (live-measured on BM and IE; see _dense_order note for why other variants keep record order).
        @inbounds for i in _ind1_order(s)
            d = t.dbh[i]; ba += BA_PER_TREE * (d * (d * t.tpa[i]))
        end
        return ba
    end
    @inbounds for i in 1:t.n; ba += t.tpa[i] * BA_PER_TREE * t.dbh[i]^2; end
    return ba
end

function stand_qmd(s::StandState)
    t = s.trees; sd2 = 0f0; tpa = 0f0
    _dp = s.variant isa BlueMountains || s.variant isa Ontario   # dense.f DP=D*P; WK5=D*DP (see note above)
    @inbounds for i in _dense_order(s)
        d = t.dbh[i]; p = t.tpa[i]
        sd2 += _dp ? d * (d * p) : p * d^2
        tpa += p
    end
    return tpa > 0f0 ? sqrt(sd2 / tpa) : 0f0
end

"""
    stand_sdimax(s) -> Float32

BA-weighted stand maximum SDI (`SDICAL`, base/sdical.f, pre-CLMAXDEN). General across variants —
the per-species SDImax (`plot.sp_sdi_def`) is variant coefficient data; the averaging is the same
base algorithm. Used by the mortality SDImax cap and the structure-stage PCTSMX demotion (BTSDIX).
"""
# SDICAL(0,XMAX) (sdical.f, byte-identical in every variant build): TREEBA=0.0054542*DBH*DBH*PROB (left-assoc,
# REAL*4) accumulated in IND1 order (SPESRT species groups) into per-species BAXSP and TOTBA; then
# XMAX = Σ_sp SDIDEF(sp)·BAXSP(sp) over species 1..MAXSP, / TOTBA. The per-species-then-species-sum structure
# rounds differently from a per-tree Σ SDIDEF·TREEBA — measured on BM 448369010497 cycle 2: CONST=SDIMAX/K
# live 464FA6A5 vs the per-tree form 464FA6A8.
function stand_sdimax(s::StandState)
    t = s.trees; p = s.plot
    t.n == 0 && return 1f0
    baxsp = s.scratch.sdi_baxsp; fill!(baxsp, 0f0); totba = 0f0
    @inbounds for i in _ind1_order(s)
        tb = 0.0054542f0 * t.dbh[i] * t.dbh[i] * t.tpa[i]
        baxsp[t.species[i]] += tb
        totba += tb
    end
    totba <= 0f0 && return 1f0
    xmax = 0f0
    @inbounds for sp in eachindex(baxsp)
        xmax += p.sp_sdi_def[sp] * baxsp[sp]
    end
    return xmax / totba
end

"""
    stand_top_height(state)

Average height of the largest-diameter 40 trees/acre (AVHT40, the summary "top
height"). Trees are taken in descending-DBH order; the last one is prorated to
hit exactly 40 TPA. (Uses a sort — fine for once-per-cycle stats, not the hotpath.)
"""
# Variants whose AVHT40/DENSE walk FVS's own IND lifecycle: CRATET's IND at cycle 0 ({v}/cratet.f RDPSRT(.FALSE.) on
# IND1 / RDPSRT(.TRUE.) with dead records), then gradd.f:186's fresh RDPSRT(DBH,.TRUE.). The cratet.f sort blocks are
# byte-identical in bm/ci/ut/tt (bm 159-166/270, ci 226-233/337, ut 214-221/325, tt 207-214/318) and ca/so
# (139-142/246), ws (215-218/322); their gradd.f:186 is the same fresh RDPSRT(DBH,IND,.TRUE.). The block is in
# EVERY variant's cratet.f (cr 139-146/250, kt 144-151/259, ie 182-189/293, em 146-153/257, ec 195-202/306,
# nc/wc 136-143/247, ak 150-157/267, bc 124-131, oc 463-470, op 452-459, sn 150-157, ls 135-142), so it is the
# base-code lifecycle, not a per-variant choice.
_fvs_ind_lifecycle(v) = true

function stand_top_height(s::StandState; cratet_ind::Bool = false, legacy_double::Bool = false)
    t = s.trees
    t.n == 0 && return 0f0
    # BM follows FVS's IND lifecycle exactly (dense.f:285-297 / avht40.f walk the CURRENT IND, no own sort):
    # CRATET's IND at cycle 0 (bm_cratet_ind!), a fresh RDPSRT(DBH,.TRUE.) everywhere else (gradd.f:186,
    # cuts.f:302/1840, esnutr.f:129/325). The empirical double sort below stays for the other variants.
    # CI too: ci/cratet.f:230-233/:337 and ci/gradd.f:186 are the same pair of sorts (FIA 753188889290487 cycle-1
    # ATAVH live 71.15 = the IND1-seeded walk; the double sort gave 70.95 ⇒ every CIVAR RELHT/PCTRED off).
    if _fvs_ind_lifecycle(s.variant) && !legacy_double
        idx = view(s.scratch.stat_idx, 1:t.n)
        cratet_ind ? bm_cratet_ind!(s, idx) : _rdpsrt!(view(t.dbh, 1:t.n), idx)
        avh = 0f0; ssumn = 0f0
        for k in 1:t.n
            ii = Int(idx[k]); p = t.tpa[ii]
            ssumn + p > 40f0 && (p = 40f0 - ssumn)
            ssumn += p; avh += t.height[ii] * p
            ssumn >= 40f0 && break
        end
        return ssumn > 0f0 ? avh / ssumn : 0f0
    end
    # avht40.f sorts IND with FVS's RDPSRT (Scowen quickersort, descending DBH) — NOT a stable sort. The
    # tie-break among equal-DBH trees decides WHICH tree lands at the 40-TPA boundary (and so its height
    # enters AVH), so a stable `sortperm!` (ascending-index ties) diverges from live on tie-heavy stands.
    # (The DG `point_basal_area!` also sorts by DBH but its BAL is an order-independent sum ⇒ tie order is
    # inert there; only AVH exposes it.) Use the ported `_rdpsrt!` to match FVS's IND tie-break exactly.
    idx = view(s.scratch.stat_idx, 1:t.n)
    dbhv = view(t.dbh, 1:t.n)
    # FVS cratet.f computes the AVHT40 top-height IND by a DOUBLE sort: IND1 = a fresh RDPSRT(.TRUE.), then
    # RDPSRT(.FALSE.) re-sorts that pre-ordered IND. RDPSRT is unstable, so the `.FALSE.` pass SWAPS equal-DBH
    # ties (the later-read record lands at the 40-tpa boundary), which a single sort does NOT — this is the
    # cycle-0 top-height divergence on tie-heavy stands (equal-DBH, different-height trees at the boundary).
    # NOTE (dig-session #2): the per-cycle FVS path is gradd.f:186 CALL RDPSRT(...,.TRUE.) → dense.f (single
    # sort), while cratet.f cycle-0 empirically double-sorts (dig-session #1). An empirical single-vs-double
    # sweep over 4 tie-heavy stands REFUTED a global per-cycle-single fix: single fixes stand 232271267010854
    # (2003) but REGRESSES the two dig-session #1 stands (1737985937290487 2024/2034, 163925866010854 1976);
    # a 3rd (202594547010854) is sort-INDEPENDENT (genuine small-tree height ULP). The correct tie-break is
    # stand-dependent because RDPSRT is an unstable quicksort on tied DBHs — no global sort choice is bit-exact.
    # Double matches the most stands (2/4 fully), so it stays. Residual TopHt swings on tie-heavy dense stands
    # are the cornered AVHT40 top-height tie-break ULP primitive (density BA/SDI/CCF preserved; converges).
    _rdpsrt!(dbhv, idx)                 # LSEQ=.TRUE. → IND1
    _rdpsrt!(dbhv, idx; lseq = false)   # LSEQ=.FALSE. → re-sort preserving IND1 (swaps ties, matches FVS)
    avh = 0f0; ssumn = 0f0
    for k in 1:t.n
        ii = Int(idx[k])
        p = t.tpa[ii]
        ssumn + p > 40f0 && (p = 40f0 - ssumn)
        ssumn += p
        avh += t.height[ii] * p
        ssumn >= 40f0 && break
    end
    return ssumn > 0f0 ? avh / ssumn : 0f0
end

"""
    point_basal_area!(state)

Fill `density.point_ba[ip]` (PTBAA, per-point basal area) AND `density.point_bal[i]`
(PTBALT, the BA in trees LARGER than tree i on the same point). For each subplot,
trees are taken in descending-DBH order and BA accumulates; PTBALT[i] is the sum
before tree i. The diameter-growth competition term uses PTBALT (= pbal). Per-tree
BA = tpa·0.005454154·DBH²·PI/GROSPC. (Sorts per point — once-per-cycle, not hotpath.)
"""
function point_basal_area!(s::StandState)
    p, t = s.plot, s.trees
    pb = s.density.point_ba; pbal = s.density.point_bal
    fill!(pb, 0f0)
    scale = p.pi / p.gross_space
    npts = 0
    @inbounds for i in 1:t.n
        npts = max(npts, Int(t.plot_id[i]))
        pbal[i] = 0f0
    end
    # ptbal.f SELECT CASE (VARACD): the eastern TWIGS variants CS/LS/NE/ON leave PTBALT = PTBAA = 0 and return
    # (live FVSon_g16 FVS_TreeList_East_Metric PtBAL is 0 on every record).
    (s.variant isa CentralStates || s.variant isa LakeStates || s.variant isa Northeast ||
     s.variant isa Ontario) && return s
    # FVS ptbal.f accumulates PTBALT per point in IND order = RDPSRT(ITRN,DBH,IND,.TRUE.) — Scowen's UNSTABLE
    # Quickersort DBH-descending, NOT a stable sort. Use the ported `_rdpsrt!` so equal-DBH tie-break matches
    # FVS's IND (a stable sortperm! diverges on tie-heavy points; inert for IE which uses PCT not PTBALT, but
    # PTBALT-consuming variants — SN calibration, PN/WC/AK/OP — need the RDPSRT order).
    order = view(s.scratch.stat_idx, 1:t.n)
    _rdpsrt!(view(t.dbh, 1:t.n), order)                                                 # IND: DBH descending, FVS tie-break
    @inbounds for i in order
        ip = Int(t.plot_id[i])
        pbal[i] = pb[ip]                                # BA already accumulated = larger trees
        pb[ip] += t.tpa[i] * BA_PER_TREE * t.dbh[i]^2 * scale
    end
    return s
end

"""
    point_density!(state)

Fill `density.point_ccf[ip]` (PCCF) and `density.point_tpa[ip]` (PTPA) — the per-point crown
competition factor and trees-per-acre (dense.f:210-211). Each tree contributes its open-grown crown
area (the same CCFT as `stand_ccf`) and TPA to its OWN subplot, scaled by PI/GROSPC (`p.pi/p.gross_space`,
the same scale `point_basal_area!` uses) so the point value is the gross per-acre density on that point.
Consumed by the multi-point regen crown ratio (regent.f:178 `CR=0.89722−0.0000461·PCCF`) and the
TCONDMLT point weights (cuts.f:1074 `+PBAWT·PTBAA+PCCFWT·PCCF+PTPAWT·PTPA`). (Once per cycle, not hotpath.)
"""
function point_density!(s::StandState)
    p, t = s.plot, s.trees
    pccf = s.density.point_ccf; ptpa = s.density.point_tpa
    fill!(pccf, 0f0); fill!(ptpa, 0f0)
    pi_f = p.pi; gross = p.gross_space
    kt = s.variant isa Kootenai
    ie = s.variant isa InlandEmpire
    # dense.f accumulates PCCF/PTPA inside DO 50 ISPC / DO 10 I3 / I=IND1(I3) — IND1 (SPESRT) order, not record
    # order; the Float32 sums round differently (BM 30193202010497 cyc1 PCCF 427D72D4 live vs 427D72D3 record-order).
    @inbounds for i in _ind1_order(s)
        ip = Int(t.plot_id[i])
        (1 <= ip <= length(pccf)) || continue
        local ccft
        if kt
            # KT PCCF uses the same per-tree ccfcal polynomial as RELDEN (kt/ccfcal.f), not crown-width area.
            ccft = kt_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]
        elseif ie
            ccft = ie_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]   # ie/ccfcal.f MODE=1
        elseif s.variant isa EasternMontana
            ccft = em_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]   # em/ccfcal.f MODE=1 (for PCCF in dgf! DGPCC term)
        elseif s.variant isa Teton
            ccft = tt_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]   # tt/ccfcal.f MODE=1
        elseif s.variant isa Utah
            ccft = ut_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]   # ut/ccfcal.f MODE=1 (PCCF for dgf! DGPTCC)
        elseif s.variant isa BlueMountains
            ccft = bm_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]   # bm/ccfcal.f MODE=1 (PCCF for dgf! DGPCCF)
        elseif s.variant isa CentralIdaho
            ccft = ci_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]   # ci/ccfcal.f MODE=1 (PCCF for dgf! DGPCCF)
        elseif s.variant isa Klamath
            ccft = nc_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]   # nc/ccfcal.f MODE=1 (PCCF for dgf! DGPCCF)
        elseif s.variant isa WestCascades
            ccft = wc_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]   # wc/ccfcal.f MODE=1 (PCCF, RELDEN)
        elseif s.variant isa PacificNorthwest
            ccft = pn_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]   # pn/ccfcal.f MODE=1
        elseif s.variant isa Olympic
            ccft = pn_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]   # op/ccfcal.f MODE=1 (byte-identical to pn) — PCCF for dgf! DGPCCF term (SP/PP)
        elseif s.variant isa EastCascades
            ccft = ec_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]   # ec/ccfcal.f MODE=1 (per-species)
        elseif s.variant isa SoutheastAlaska
            ccft = ak_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]   # ak/ccfcal.f MODE=1
        elseif s.variant isa OregonCoast
            ccft = oc_tree_ccf(Int(t.species[i]), t.dbh[i], t.height[i]) * t.tpa[i]  # oc/ccfcal.f MODE=1 (R5CRWD)
        elseif s.variant isa CentralCalifornia
            ccft = ca_tree_ccf(Int(t.species[i]), t.dbh[i], t.height[i]) * t.tpa[i]  # ca/ccfcal.f MODE=1 (R5CRWD = OC's)
        elseif s.variant isa SouthCentralOregon
            ccft = so_tree_ccf(Int(t.species[i]), t.dbh[i], t.height[i]) * t.tpa[i]  # so/ccfcal.f MODE=1 — was the generic
                                                                                       # national crown-width path ⇒ PCCF ~100× low (DUBSCR TPCCF 1.4 vs live 153)
        elseif s.variant isa WestSierra
            ccft = ws_ccft(Int(t.species[i]), t.dbh[i], t.height[i], t.tpa[i])  # ws/ccfcal.f MODE=1 (same gap as SO)
        else
            cw  = s.variant isa CentralRockies ?
                  cr_crown_width(Int(t.species[i]), t.dbh[i], Int(p.model_type)) :
                  crown_width(s.coef, s.species.code2[t.species[i]], t.dbh[i], t.height[i], 90, 1,
                              p.latitude, p.longitude, p.elevation)
            ccft = t.dbh[i] > 0.1f0 ? 0.001803f0 * cw * cw * t.tpa[i] : 0.001f0 * t.tpa[i]
        end
        # dense.f:210-211 accumulates each term as `CCFT*PI/GROSPC` — i.e. (ccft·pi)/gross evaluated
        # left-to-right, NOT ccft·(pi/gross) with a precomputed reciprocal-scale. The two differ by ~1
        # Float32 ULP per term; on the dense estab_pccf points that sub-ULP tips a regen-crown INT(CR·100+0.5)
        # boundary. Match FVS's exact op order to deconfound (doctrine #8).
        pccf[ip] += ccft * pi_f / gross
        ptpa[ip] += t.tpa[i] * pi_f / gross
    end
    return s
end

"""
    stand_pct!(state)

Fill `trees.crown_ratio[i]` with PCT, the stand basal-area percentile (PCTILE,
pctile.f via dense.f): trees in descending-DBH order, `PCT[i] = (BA of tree i and
all smaller) / total · 100`. So `1 − PCT/100` is the fraction of stand BA in larger
trees, which the diameter-growth competition term uses. (Despite the field name,
this is FVS's PCT array, not the crown ratio — the crown ratio is `crown_pct`/ICR.)
"""
function stand_pct!(s::StandState; cratet_ind::Bool = false)
    t = s.trees; n = t.n
    n == 0 && return s
    if cratet_ind                                        # BM first grow cycle: CRATET's IND (see bm_cratet_ind!)
        idx = view(s.scratch.stat_idx, 1:n)
        bm_cratet_ind!(s, idx)
        _pctile!(t.crown_ratio, t, idx, n)
        return s
    end
    # PCT is built over FVS's IND = the per-cycle DBH-descending order from gradd.f:186
    # `CALL RDPSRT(ITRN,DBH,IND,.TRUE.)` feeding dense.f/PCTILE. RDPSRT is Scowen's UNSTABLE
    # Quickersort, so equal-DBH ties resolve by the partition order, NOT ascending index. A stable
    # `sortperm!` diverges from FVS on tie-heavy stands: e.g. a dense all-0.1"-seedling regen plot
    # (5 records, one dominant-TPA loblolly) gets its dominant assigned a LOW percentile (own BA
    # fraction) instead of ~100, which inverts VARMRT's self-thinning kill (kills the dominant),
    # collapses the survivor QMD, and makes the morts QMD-convergence loop LIMIT-CYCLE instead of
    # converging — a 2–3× first-cycle mortality error. Use `_rdpsrt!` (single .TRUE. sort) to match.
    idx = view(s.scratch.stat_idx, 1:n)
    _rdpsrt!(view(t.dbh, 1:n), idx)                     # IND: DBH descending, FVS tie-break
    _pctile!(t.crown_ratio, t, idx, n)
    return s
end

# PCTILE (pctile.f) over DENSE's WK5 (dense.f:186-187: DP=D*P; WK5=D*DP), in FVS's exact single-precision
# order — identical in all 24 variant builds. Cumulative from the smallest (IND bottom) up; TOT = the top
# record's cumulative; every other record is divided by PCTIN1 = TOT/100. (NOT ×100/TOT); the top record is set
# to exactly 100. jl's former (D*D)*P and cum/TOT*100 each differed by 1 ULP on a share of records (measured vs
# FVSsn treeszcp_cap cycle 1: 8/27 PCT, 2 EFFTR), which VARMRT's geometric kill amplifies.
function _pctile!(pct::AbstractVector{Float32}, t, idx, n::Int)
    n == 1 && (pct[Int(idx[1])] = 100f0; return pct)   # pctile.f: PERCNT(1)=100, IF(N.LE.1) RETURN
    cum = 0f0
    @inbounds for k in n:-1:1
        ii = Int(idx[k])
        cum += t.dbh[ii] * (t.dbh[ii] * t.tpa[ii])      # WK5 = D*(D*P)
        pct[ii] = cum
    end
    i1 = Int(idx[1])
    tot = pct[i1]
    pct[i1] = tot / 100f0
    tot <= 0f0 && return pct                             # pctile.f: IF(TOT.LE.0.0) RETURN
    pctin1 = pct[i1]
    @inbounds for k in 2:n
        ii = Int(idx[k]); pct[ii] = pct[ii] / pctin1
    end
    pct[i1] = 100f0
    return pct
end

"""
    stand_ccf(state)

Crown competition factor (RELDEN): Σ over trees of the open-grown crown area
(CCFCAL/ccfcal.f): `0.001803·crownwidth²·tpa` (or `0.001·tpa` for DBH ≤ 0.1).
"""
function stand_ccf(s::StandState)
    p, t = s.plot, s.trees
    ccf = 0f0
    if s.variant isa Kootenai
        # KT CCF is a direct per-tree polynomial (kt/ccfcal.f), not the crown-width→area path. Sum over the
        # CURRENT tree list (1:t.n). The dead partition is folded in ONLY when the caller has bumped t.n to
        # include it — i.e. the backdated calibration density pass (compute_density! at t.n=nlive+ndead);
        # the growth-cycle pass runs with t.n=nlive so dead are (correctly) excluded. RELDEN is stored by
        # compute_density! into p.relative_density and read by dgf!, matching FVS's DENSE→DGF flow.
        @inbounds for i in 1:t.n
            ccf += kt_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]
        end
        return ccf
    elseif s.variant isa InlandEmpire
        # IE CCF is the same direct per-species polynomial (ie/ccfcal.f MODE=1); stand CCF = Σ CCFT·P = RELDEN.
        # dense.f:168-229 accumulates it SPECIES-MAJOR in IND1 order into a per-species subtotal RELDSP(ISPC), then
        # RELDT=RELDT+RELDSP(ISPC) — a different Float32 summation order than a flat record-order sum.
        sp_cur = 0; relsp = 0f0
        @inbounds for i in _ind1_order(s)
            sp = Int(t.species[i])
            if sp != sp_cur
                sp_cur == 0 || (ccf += relsp)
                sp_cur = sp; relsp = 0f0
            end
            relsp += ie_tree_ccf(sp, t.dbh[i]) * t.tpa[i]
        end
        sp_cur == 0 || (ccf += relsp)
        return ccf
    elseif s.variant isa Ontario
        # ON CCF = ccfcal.f (LS form) → cwcalc.f open-grown crown WIDTH (IWHO=1, CR=90) via the
        # ISPC→US-code (JSP2) remap, then 0.001803·CW²·P. HI needs stand lat/long/elev.
        lat = p.latitude; long = p.longitude; elev = p.elevation
        @inbounds for i in 1:t.n
            ccf += on_tree_ccf(Int(t.species[i]), t.dbh[i], lat, long, elev) * t.tpa[i]
        end
        return ccf
    elseif s.variant isa BritishColumbia
        # BC CCF is the same direct per-species polynomial (bc/ccfcal.f MODE=1); bc_tree_ccf folds in ×P.
        @inbounds for i in 1:t.n
            ccf += bc_tree_ccf(Int(t.species[i]), t.dbh[i], t.tpa[i])
        end
        return ccf
    elseif s.variant isa EasternMontana
        # EM CCF is the same direct per-species polynomial (em/ccfcal.f MODE=1, Paine-Hann/NI form).
        @inbounds for i in 1:t.n
            ccf += em_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]
        end
        return ccf
    elseif s.variant isa Utah
        # UT CCF is the same direct per-species polynomial (ut/ccfcal.f MODE=1); stand CCF = Σ CCFT·P = RELDEN.
        @inbounds for i in 1:t.n
            ccf += ut_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]
        end
        return ccf
    elseif s.variant isa BlueMountains
        # BM CCF is the same direct per-species polynomial (bm/ccfcal.f MODE=1); stand CCF = Σ CCFT·P = RELDEN.
        @inbounds for i in 1:t.n
            ccf += bm_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]
        end
        return ccf
    elseif s.variant isa Klamath
        # NC CCF = the direct per-species ccfcal polynomial (nc/ccfcal.f MODE=1); stand CCF = Σ CCFT·P = RELDEN.
        @inbounds for i in 1:t.n
            ccf += nc_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]
        end
        return ccf
    elseif s.variant isa WestCascades
        # WC CCF = the direct 16-group ccfcal polynomial (wc/ccfcal.f MODE=1); stand CCF = Σ CCFT·P = RELDEN.
        @inbounds for i in 1:t.n
            ccf += wc_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]
        end
        return ccf
    elseif s.variant isa PacificNorthwest
        @inbounds for i in 1:t.n
            ccf += pn_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]
        end
        return ccf
    elseif s.variant isa Olympic
        # OP CCF = op/ccfcal.f MODE=1, a direct per-species RD1+RD2·D+RD3·D² polynomial. op/ccfcal.f
        # carries the "PN $Id$" header and its INDCCF/RD1/RD2/RD3 DATA blocks are byte-identical to
        # pn/ccfcal.f, so `pn_tree_ccf` reproduces it exactly; stand CCF = Σ CCFT·P.
        @inbounds for i in 1:t.n
            ccf += pn_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]
        end
        return ccf
    elseif s.variant isa EastCascades
        @inbounds for i in 1:t.n
            ccf += ec_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]
        end
        return ccf
    elseif s.variant isa CentralIdaho
        # CI CCF is the same direct per-species polynomial (ci/ccfcal.f MODE=1); stand CCF = Σ CCFT·P = RELDEN,
        # accumulated SPECIES-MAJOR (dense.f RELDSP(ISPC) subtotals, as IE). The flat record-order sum was a few ULP
        # off: REGCAL fixture backdated RELDEN 215.12527 vs live 215.1252 ⇒ CW/OH SNX 610.90 vs 610.91.
        sp_cur = 0; relsp = 0f0
        @inbounds for i in _ind1_order(s)
            sp = Int(t.species[i])
            if sp != sp_cur
                sp_cur == 0 || (ccf += relsp)
                sp_cur = sp; relsp = 0f0
            end
            relsp += ci_tree_ccf(sp, t.dbh[i]) * t.tpa[i]
        end
        sp_cur == 0 || (ccf += relsp)
        return ccf
    elseif s.variant isa Teton
        # TT CCF is the same direct per-species polynomial (tt/ccfcal.f MODE=1); stand CCF = Σ CCFT·P = RELDEN.
        @inbounds for i in 1:t.n
            ccf += tt_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]
        end
        return ccf
    elseif s.variant isa OregonCoast
        # OC CCF is the open-grown crown-width → area form (oc/ccfcal.f MODE=1, R5CRWD); stand CCF = Σ CCFT·P.
        @inbounds for i in 1:t.n
            ccf += oc_tree_ccf(Int(t.species[i]), t.dbh[i], t.height[i]) * t.tpa[i]
        end
        return ccf
    elseif s.variant isa CentralCalifornia
        # CA CCF = the open-grown crown-width → area form (ca/ccfcal.f MODE=1, R5CRWD = byte-identical to OC's).
        @inbounds for i in 1:t.n
            ccf += ca_tree_ccf(Int(t.species[i]), t.dbh[i], t.height[i]) * t.tpa[i]
        end
        return ccf
    elseif s.variant isa SoutheastAlaska
        # AK CCF is the open-grown crown-width → area form (ak/ccfcal.f MODE=1); stand CCF = Σ CCFT·P.
        @inbounds for i in 1:t.n
            ccf += ak_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]
        end
        return ccf
    elseif s.variant isa SouthCentralOregon
        # SO CCF = so/ccfcal.f MODE=1 (RD polynomial + WC-hardwood + SH/WO r6crwd crown-width²);
        # stand CCF = Σ CCFT·P = RELDEN, read by dgf! CONSPP and regent PCTRED.
        @inbounds for i in 1:t.n
            ccf += so_tree_ccf(Int(t.species[i]), t.dbh[i], t.height[i]) * t.tpa[i]
        end
        return ccf
    elseif s.variant isa WestSierra
        # WS CCF = ws/ccfcal.f MODE=1 (crown-width² (RD1+D·RD2)²·0.001803 native; GB/MC/CA specials);
        # stand CCF = Σ CCFT·P = RELDEN, read by the crown-ratio SCALE (ws/crown.f).
        @inbounds for i in 1:t.n
            ccf += ws_ccft(Int(t.species[i]), t.dbh[i], t.height[i], t.tpa[i])
        end
        return ccf
    end
    @inbounds for i in 1:t.n
        sp = t.species[i]
        cw = s.variant isa CentralRockies ?
             cr_crown_width(Int(sp), t.dbh[i], Int(p.model_type)) :
             crown_width(s.coef, s.species.code2[sp], t.dbh[i], t.height[i], 90, 1,
                         p.latitude, p.longitude, p.elevation)
        ccf += t.dbh[i] > 0.1f0 ? 0.001803f0 * cw * cw * t.tpa[i] : 0.001f0 * t.tpa[i]
    end
    return ccf
end

"""
    stand_sdi(s)

Reported `.sum` stand density index (SDICLS, sdical.f), following the SDICALC method flag
`zeide_sdi` (LZEIDE, SN default Zeide) — the SAME flag the SDImax mortality uses, so the two
stay consistent. **Zeide:** Σ TPA·(D/10)^1.605 over `D ≥ DBHZEIDE` (sdical.f:326). **Reineke:**
the `SDI = SPROB·A + B·SDSQ` Taylor form over `D ≥ DBHSTAGE` (sdical.f:281-327). Defaults
(Zeide, threshold 0) reproduce the prior behavior.
"""
# The reported stand SDI (.sum SDI = SDIBC before a thin, SDIAC after): SDICLS(0,0.,999.,1,...) (sdical.f ENTRY
# SDICLS; fvs.f:440, grincr.f:241) — identical in every variant build. Both loops walk IND1 (SPESRT) order.
#   pass 1 (DBH>=DBHSTAGE): SDSQ=SDSQ+(DBH**2.0)*PROB; SPROB=SPROB+PROB  → A,B
#   pass 2: SDIC  = SDIC  + (A+B*(DBH**2.0))*PROB          (DBH>=DBHSTAGE)   — PER TREE, not SPROB*A+B*SDSQ
#           SDIC2 = SDIC2 + PROB*(DBH/10.)**1.605           (DBH>=DBHZEIDE)
# disply.f:332-338 reports SDIC2 when LZEIDE else SDIC. DBH**2.0 / **1.605 are gfortran powf. CROWN's SDIAC is the
# same SDIC pass (`stand_sdi_reineke`), whatever LZEIDE.
function stand_sdi(s::StandState)
    t = s.trees
    t.n == 0 && return 0f0
    ord = _ind1_order(s)
    if s.control.zeide_sdi
        thr = s.control.dbh_zeide; sdi2 = 0f0
        @inbounds for i in ord
            d = t.dbh[i]
            d >= thr && (sdi2 += t.tpa[i] * fpow(d / 10f0, 1.605f0))
        end
        return sdi2
    end
    thr = s.control.dbh_stage; sdsq = 0f0; sprob = 0f0
    @inbounds for i in ord
        d = t.dbh[i]; d < thr && continue
        sdsq += fpow(d, 2f0) * t.tpa[i]; sprob += t.tpa[i]
    end
    sprob == 0f0 && return 0f0
    k10 = fpow(10f0, -1.605f0)                     # == gfortran's folded 10.0**(-1.605) (3CCB6B13, verified)
    a = k10 * (1f0 - 1.605f0 / 2f0) * fpow(sdsq / sprob, 1.605f0 / 2f0)
    b = k10 * (1.605f0 / 2f0) * fpow(sdsq / sprob, 1.605f0 / 2f0 - 1f0)
    sdic = 0f0
    @inbounds for i in ord
        d = t.dbh[i]
        d >= thr && (sdic += (a + b * fpow(d, 2f0)) * t.tpa[i])
    end
    return sdic
end

"""
    stand_sdi_reineke(s)

SDICLS(0,0.,999.,1,SDIC,…) — the STAGE (Reineke-summation) stand SDI that CROWN reads (SDIAC/SDIBC, grincr.f:241,323;
fvs.f:196). sdical.f:260-283 first sums SDSQ=Σ(DBH**2)·PROB and SPROB=ΣPROB over the DBH≥DBHSTAGE records in IND1
order to form the STAGE A/B, then DISCARDS the closed form SPROB·A+B·SDSQ and re-sums SDIC=Σ(A+B·DBH**2)·PROB in the
same IND1 order (:293-332). Float32 addition is order-dependent, so both passes walk `_ind1_order` (as `stand_sdi`,
dense.f's BA/PCCF and SDICAL already do) — the index-order closed form sat ~1e-6 relative off live's SDIAC
(TT S248112 2030 SDIAC 318.5158 vs live 318.5154), enough to flip an INT(CRNEW+0.5) crown.
"""
function stand_sdi_reineke(s::StandState)
    t = s.trees
    ord = _ind1_order(s)
    thr = s.control.dbh_stage; sprob = 0f0; sdsq = 0f0
    @inbounds for i in ord
        d = t.dbh[i]; d < thr && continue
        sdsq += (d * d) * t.tpa[i]; sprob += t.tpa[i]
    end
    sprob == 0f0 && return 0f0
    # sdical.f:281-282 `(10.0**(-1.605))*…*((SDSQ/SPROB)**(1.605/2.))` — all FVS `**` = gfortran powf, route via
    # the companion (doctrine #8) not Julia's openlibm `^`.
    a = fpow(10f0, -1.605f0) * (1f0 - 1.605f0 / 2f0) * fpow(sdsq / sprob, 1.605f0 / 2f0)
    b = fpow(10f0, -1.605f0) * (1.605f0 / 2f0) * fpow(sdsq / sprob, 1.605f0 / 2f0 - 1f0)
    sdic = 0f0
    @inbounds for i in ord
        d = t.dbh[i]
        d >= thr && (sdic += (a + b * (d * d)) * t.tpa[i])
    end
    return sdic
end
