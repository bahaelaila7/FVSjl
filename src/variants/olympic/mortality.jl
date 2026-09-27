# =============================================================================
# mortality.jl (olympic) — OP mortality (op/morts.f).
#
# op/morts.f is NOT a cooperating split: whenever ORGANON ran this cycle (SMORMT = ΣMORTEXP > 0),
# the ORGANON per-record mortality MORTEXP is used for EVERY tree record — IORG=0 and IORG=1 alike
# (op/morts.f:13-19,331-337): `WKI = MORTEXP(I)·(FINT/5)`, capped at PROB. Only when NO valid ORGANON
# tree exists (SMORMT=0) does op fall back to its FVS-native RIP mortality (the per-species BM0..BM5
# logistic + Gould–Harrington small-tree model). "No valid ORGANON tree" ⇔ the stand has NO big-6
# (GF sp 3 / DF sp 16) record (op/dgdriv.f big-6 gate ⇒ op_build_organon_buffer! sets runs=false ⇒
# op_org_ran=false), because ORGANON PREPARE/EXECUTE gate on a big-6.  For S248112 (has DF) ORGANON
# always runs and the MORTEXP-for-all path is validated; a no-big-6 OP stand exercises the RIP path.
#
# THE RIP PATH (op/morts.f:302-506, IMRTMP = MORTMAP(ISPC)):
#   CASE 1: RIP = BM0 + BM1·√D + BM3·CR^.25 + BM4·(XSITE1+4.5) + BM5·BAL       (5-yr logistic)
#   CASE 2: RIP = BM0 + BM1·D + BM4·(XSITE1+4.5) + BM5·(BAL/D)
#   CASE 3: RIP = BM0 + BM1·D + BM2·D² + BM3·CR + BM4·(XSITE2+4.5) + BM5·BAL
#   CASE 4: RIP = BM0 + BM1·D + BM2·D² + BM3·CR + BM4·(XSITE1+4.5) + BM5·BAL
#     (1-4): RIP=1/(1+e^−RIP) [5-yr] → (1−RIP)^0.2 [annual survival] → 1−RIP·CRADJ [annual RIP].
#   CASE 5: RELHT=min(HT/AVH,1.5); RIP=−6.6707+0.5105·ln(5+BA)−1.3183·RELHT; RIP=1/(1+e^RIP); 1−RIP.
#   CASE 6 (RW sp 17): RIP=2.901447+0.578694·D−0.001793·PTBAL; RIP=1/(1+e^RIP).
#   SMALL TREE (D<3", sp≠17, Gould–Harrington; OVERRIDES the case): DBHA=D+BETA·min(HT,4.5);
#     RIP=ALPHA1+ALPHA2·(PTBAL·MVALUES[MCLASS]/√(DBHA+1))+ALPHA3·RELHT; RIP=1/(1+e^RIP); 1−RIP.
#   Then D floored at 0.5 (morts.f:341, BEFORE the cases + small-tree); MORTMULT window X; RIP floor
#   0.001 (morts.f:475); WKI=P·(1−(1−RIP)^FINT)·X; SIZCAP up-cap (morts.f:483, size-cap default 999 ⇒
#   inert); WKI≤P; SDIMAX<5 ⇒ kill all. CR=ICR·0.01, BAL=(1−PCT/100)·BA, PTBAL=PTBALT, XSITE1/2=
#   SITEAR(16/19), CRADJ=1−e^−(25·CR)² for CR≤0.17. All inputs are START-of-cycle (D=DBH pre-growth,
#   BA/PCT/PTBAL/AVH from grow_cycle!'s compute_density!) — MORTS runs before UPDATE applies the DG.
#   DGT/G/RELDBH/AVED (morts.f:280-347) are DEAD in the RIP path (leftover WC/PN scaffolding); omitted.
#
# The IMODTY==3 (SMC) DF override (morts.f:380-387) is UNREACHABLE on the RIP path: it fires only for
# ISPC==16 (DF), which is big-6, so ORGANON always runs for a DF stand ⇒ the MORTEXP path, not here.
# The op/morts.f SDIMAX-cap over-density iterative pass (morts.f:508-547, WKI·=PASS while SDIA≥SDIMAX
# or BAA>550) and CLMORTS (climate mort) are omitted (PASS=1); both are inert for a below-SDImax stand
# with no CLIMATE block — PROVEN inert by the bit-exact A/B (if either fired in the oracle, TPA/SDI would
# diverge; the stand peaks at SDI 519 < SDIMAX and BA 313 < 550, so PASS stays 1). A/B vs FVSop_clean on
# opnb6 (S248112 all-WF/LP/PP/SP, DGSD=0 deterministic, NUMCYCLE=10) is BIT-EXACT on TPA (536→344), BA,
# SDI and the accretion/mortality growth columns every cycle (test_op_native_mortality_rip.jl); the only
# residual is the cubic/board-foot VOLUME columns from cycle 1+ (max BdFt +1.4%, MCuFt +1.0%) plus one QMD
# rounding straddle (2020: 9.9 vs 10.0) — the pre-existing OC/OP ORGANON GMV merch-cubic volume corner
# (identical surviving population, cycle-0 volume bit-exact ⇒ not this mortality port). FIXMORT is wired
# uniformly (apply_fixmort!, inert without a FIXMORT keyword). Deterministic (DGSD=0). FINT=5 ⇒ FINT/5 = 1.
# =============================================================================

# op/morts.f DATA tables (index by FVS species 1..39; op MAXSP=39, order SF,WF,GF,...,__,OT).
const OP_MORTMAP = Int32[
    2,2,2,2,2,2,2,4,4,2,
    4,4,4,4,4,1,6,3,3,3,
    4,4,4,4,4,4,4,5,4,4,
    4,4,4,4,4,4,4,1,1]
const OP_BM0 = Float32[
    -7.60159, -7.60159, -7.60159, -7.60159, -7.60159,
    -7.60159, -7.60159, -1.922689902, -1.922689902, -7.60159,
    -1.050000682, -1.050000682, -1.050000682, -1.050000682, -1.050000682,
    -4.13142, -4.13142, -0.761609, -0.761609, -0.761609,
    -2.976822456, -2.0, -2.0, -2.0, -4.317549852,
    -2.0, -2.0, -0.0, -1.050000682, -1.050000682,
    -1.050000682, -1.050000682, -4.072781265, -3.020345211, -3.020345211,
    -3.020345211, -2.0, -4.13142, -4.13412]
const OP_BM1 = Float32[
    -0.200523, -0.200523, -0.200523, -0.200523, -0.200523,
    -0.200523, -0.200523, -0.136081990, -0.136081990, -0.200523,
    -0.194363402, -0.194363402, -0.194363402, -0.194363402, -0.194363402,
    -1.13736, -1.13736, -0.529366, -0.529366, -0.529366,
    0.0, -0.5, -0.5, -0.5, -0.057696253,
    -0.5, -0.5, -0.0, -0.194363402, -0.194363402,
    -0.194363402, -0.194363402, -0.176433475, 0.0, 0.0,
    0.0, -0.5, -1.13736, -1.13736]
const OP_BM2 = Float32[
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.002479863, 0.002479863, 0.0,
    0.003803100, 0.003803100, 0.003803100, 0.003803100, 0.003803100,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.015, 0.015, 0.015, 0.0,
    0.015, 0.015, 0.0, 0.003803100, 0.003803100,
    0.003803100, 0.003803100, 0.0, 0.0, 0.0,
    0.0, 0.015, 0.0, 0.0]
const OP_BM3 = Float32[
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -3.178123293, -3.178123293, 0.0,
    -3.557300286, -3.557300286, -3.557300286, -3.557300286, -3.557300286,
    -0.823305, -0.823305, -4.74019, -4.74019, -4.74019,
    -6.223250962, -3.0, -3.0, -3.0, 0.0,
    -3.0, -3.0, 0.0, -3.557300286, -3.557300286,
    -3.557300286, -3.557300286, -1.729453975, -8.467882343, -8.467882343,
    -8.467882343, -3.0, -0.823305, -0.823305]
const OP_BM4 = Float32[
    0.0441333, 0.0441333, 0.0441333, 0.0441333, 0.0441333,
    0.0441333, 0.0441333, 0.0, 0.0, 0.0441333,
    0.003971638, 0.003971638, 0.003971638, 0.003971638, 0.003971638,
    0.0307749, 0.0307749, 0.0119587, 0.0119587, 0.0119587,
    0.0, 0.015, 0.015, 0.015, 0.004861355,
    0.015, 0.015, 0.0, 0.003971638, 0.003971638,
    0.003971638, 0.003971638, 0.0, 0.013966388, 0.013966388,
    0.013966388, 0.015, 0.0307749, 0.0307749]
const OP_BM5 = Float32[
    0.00063849, 0.00063849, 0.00063849, 0.00063849, 0.00063849,
    0.00063849, 0.00063849, 0.004684133, 0.004684133, 0.00063849,
    0.005573601, 0.005573601, 0.005573601, 0.005573601, 0.005573601,
    0.00991005, 0.00991005, 0.00756365, 0.00756365, 0.00756365,
    0.0, 0.01, 0.01, 0.01, 0.00998129,
    0.01, 0.01, 0.0, 0.005573601, 0.005573601,
    0.005573601, 0.005573601, 0.012525642, 0.009461545, 0.009461545,
    0.009461545, 0.01, 0.00991005, 0.00991005]
# Gould–Harrington small-tree model (morts.f:105-119).
const OP_BETA = Float32[
    0.247354, 0.217481, 0.179705, 0.205647, 0.216823,
    0.216823, 0.282203, 0.216823, 0.281542, 0.170425,
    0.168227, 0.216823, 0.216823, 0.216823, 0.236925,
    0.163506, 0.216823, 0.182940, 0.172690, 0.302866,
    0.216823, 0.216823, 0.216823, 0.216823, 0.216823,
    0.216823, 0.216823, 0.216823, 0.216823, 0.216823,
    0.216823, 0.216823, 0.216823, 0.216823, 0.216823,
    0.216823, 0.216823, 0.216823, 0.216823]
const OP_MCLASS = Int32[
    1, 2, 2, 2, 2, 2, 3, 2, 3, 2, 4, 4, 3,
    3, 4, 3, 1, 1, 1, 1, 1, 4, 3, 4, 3, 5, 5, 5, 5, 4,
    4, 5, 1, 1, 3, 3, 5, 1, 5]
const OP_MVALUES = Float32[1.00, 1.500, 2.250, 3.375, 5.062]
const OP_ALPHA   = Float32[-4.4384, 0.0053, -0.6001]

# op/morts.f:302-506 — the FVS-native per-tree RIP mortality (no valid ORGANON tree). Fills `killed`.
function _op_rip_mortality!(killed::AbstractVector{Float32}, s::StandState, n::Int, fint::Float32)
    t, p, d = s.trees, s.plot, s.density
    ba = p.basal_area
    xsite1 = length(p.sp_site_index) >= 16 ? p.sp_site_index[16] : 0f0   # SITEAR(16=DF)
    xsite2 = length(p.sp_site_index) >= 19 ? p.sp_site_index[19] : 0f0   # SITEAR(19=WH)
    avh    = p.avg_height
    sizcap = s.control.sp_size_cap
    sdimax = clim_sdical_xmax(s, stand_sdimax(s), fint)   # op/morts.f SDICAL (+ sdical.f:216 CLMAXDEN)
    year   = current_cycle_year(s)
    fscale = fint / 5f0
    is_stale = Int(t.species[n])            # op/morts.f:480 BARK=BRATIO(IS,...): IS is the stale ISP(ITRN)
    @inbounds for i in 1:n
        pr = t.tpa[i]; pr <= 0f0 && continue
        ispc = Int(t.species[i])
        dd = t.dbh[i]
        dd <= 0.5f0 && (dd = 0.5f0)                        # morts.f:341 floor (used by cases + small-tree)
        cr  = Float32(t.crown_pct[i]) * 0.01f0             # CR = ICR·0.01
        bal = (1f0 - t.crown_ratio[i] / 100f0) * ba        # BAL = (1−PCT/100)·BA
        ptbal = i <= length(d.point_bal) ? d.point_bal[i] : 0f0
        cradj = cr <= 0.17f0 ? (1f0 - exp(-(25f0 * cr)^2f0)) : 1f0
        b0 = OP_BM0[ispc]; b1 = OP_BM1[ispc]; b2 = OP_BM2[ispc]
        b3 = OP_BM3[ispc]; b4 = OP_BM4[ispc]; b5 = OP_BM5[ispc]
        rip = 0f0
        imr = OP_MORTMAP[ispc]
        if imr == 1
            rip = b0 + b1 * sqrt(dd) + b3 * cr^0.25f0 + b4 * (xsite1 + 4.5f0) + b5 * bal
            rip = 1f0 / (1f0 + exp(-rip)); rip = (1f0 - rip)^0.2f0; rip = 1f0 - rip * cradj
        elseif imr == 2
            rip = b0 + b1 * dd + b4 * (xsite1 + 4.5f0) + b5 * (bal / dd)
            rip = 1f0 / (1f0 + exp(-rip)); rip = (1f0 - rip)^0.2f0; rip = 1f0 - rip * cradj
        elseif imr == 3
            rip = b0 + b1 * dd + b2 * dd * dd + b3 * cr + b4 * (xsite2 + 4.5f0) + b5 * bal
            rip = 1f0 / (1f0 + exp(-rip)); rip = (1f0 - rip)^0.2f0; rip = 1f0 - rip * cradj
        elseif imr == 4
            rip = b0 + b1 * dd + b2 * dd * dd + b3 * cr + b4 * (xsite1 + 4.5f0) + b5 * bal
            rip = 1f0 / (1f0 + exp(-rip)); rip = (1f0 - rip)^0.2f0; rip = 1f0 - rip * cradj
        elseif imr == 5
            relht = avh > 0f0 ? t.height[i] / avh : 0f0
            relht > 1.5f0 && (relht = 1.5f0)
            rip = -6.6707f0 + 0.5105f0 * log(5f0 + ba) - 1.3183f0 * relht
            rip = 1f0 / (1f0 + exp(rip)); rip = 1f0 - rip
        elseif imr == 6
            rip = 2.901447f0 + 0.578694f0 * dd - 0.001793f0 * ptbal
            rip = 1f0 / (1f0 + exp(rip))
        end
        # Gould–Harrington small-tree model OVERRIDES the case for D<3", sp≠17 (morts.f:431-448).
        if dd < 3f0 && ispc != 17
            relht = avh > 0f0 ? t.height[i] / avh : 0f0
            relht > 1.5f0 && (relht = 1.5f0)
            hbh = t.height[i]; hbh >= 4.5f0 && (hbh = 4.5f0)
            dbha = dd + OP_BETA[ispc] * hbh
            avalue = OP_MVALUES[OP_MCLASS[ispc]]
            rip = ptbal * avalue / sqrt(dbha + 1f0)
            rip = OP_ALPHA[1] + OP_ALPHA[2] * rip
            rip = rip + OP_ALPHA[3] * relht
            rip = 1f0 / (1f0 + exp(rip)); rip = 1f0 - rip
        end
        # MORTMULT window X (morts.f:458; XMMULT/XMDIA1/XMDIA2 ⇒ active_mort_mult, 1 outside window).
        x = active_mort_mult(s.control, ispc, year, dd)
        # (Establishment "best-tree" 20-yr immunity, morts.f:464-469, is not modeled — no IESTAT.)
        rip < 0.001f0 && (rip = 0.001f0)                  # morts.f:475 minimum annual RIP
        wki = pr * (1f0 - (1f0 - rip)^fint) * x
        # SIZCAP up-cap (morts.f:480-488): G on the CURRENT DG and the stale-IS bark; inert at default 999.
        bark = op_bratio(is_stale, dd)
        g = (t.diam_growth[i] / bark) * fscale
        idmflg = round(Int, sizcap[ispc, 3])
        if (dd + g) >= sizcap[ispc, 1] && idmflg != 1
            wki = max(wki, pr * sizcap[ispc, 2] * fscale)
        end
        wki > pr && (wki = pr)
        sdimax < 5f0 && (wki = pr)                        # morts.f:495 climate kill-all
        killed[i] = wki
    end
    return killed
end

function mortality!(s::StandState, ::Olympic; fint::Float32 = 5.0f0, book_snags::Bool = true)
    t, c = s.trees, s.calib
    n = t.n; n == 0 && return _clim_mort_empty!(s, fint)   # ITRN=0 ⇒ morts.f still reaches CLMORTS
    killed = zeros(Float32, n)
    if c.op_org_ran && length(c.op_mortexp) == n
        # ORGANON MORTEXP for ALL records (op/morts.f:331-336): WKI = MORTEXP·(FINT/5), capped at PROB.
        fscale = fint / 5f0
        @inbounds for i in 1:n
            p = t.tpa[i]; p <= 0f0 && continue
            wki = c.op_mortexp[i] * fscale
            wki > p && (wki = p)
            killed[i] = wki
        end
    else
        # No valid ORGANON tree (SMORMT=0): the FVS-native RIP mortality (op/morts.f:302-506).
        _op_rip_mortality!(killed, s, n, fint)
    end
    # FIXMORT (morts.f:582-857) applies after both paths; inert without a FIXMORT keyword.
    @label morts45
    # Climate-FVS mortality (op/morts.f CALL CLMORTS — after the base mortality and TPAMRT, immediately before
    # FIXMORT), THISYR = IY(ICYC)+FINT/2. Inert unless a CLIMATE keyword activated s.climate.
    (s.climate !== nothing && s.climate.active) &&
        apply_climate_mort!(s, killed, Float32(current_cycle_year(s)) + fint / 2f0, fint)
    apply_fixmort!(s, killed, n, fint)
    # MISMRT (gradd.f:96 MISTOE → mismrt.f:185-191, misintop.f APMC): WK2=MAX(WK2,PROB·rate). OP never triples
    # (ICL4=0), so it always lands here. Inert unless a record carries DMR.
    _ie_mis_variant(s.variant) && ie_dm_mortality_combine!(killed, s, fint, n)
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end
