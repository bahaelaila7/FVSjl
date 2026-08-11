# =============================================================================
# mortality.jl (easternmontana) — EM mortality (em/morts.f).
#
# Stand-level (= KT structure): grown-stand DQ10, CONST=SDIMAX/0.02483133, TMD10=CONST·D10^−1.605,
# T85D10=TMD10·PMSDIU, T55D10=TMD10·PMSDIL. SDI self-thinning target TN10 → RN=1−(1−(T−TN10)/T)^(1/FINT).
#   T ≤ T55D10  → TN10=T, RN=0 (below the SDI limit — the emt01 regime; PMSC background dominates).
#   T > T85D10  → TN10=T85D10 (kill to the 85% line).
#   middle (55–85%): iterative linear-fn fit (em/morts.f label 220) — DEFERRED (emt01 is below-SDI).
# Per-tree, ORIGINAL EM species (1-3,7-10,18): RI=0.5/(1+exp(PMSC+PMD·D+PMDSQ·D²)); RIP=RN if SDI limiting
#   (T>TEM and RN>0) else RI; WKI=P·(1−(1−RIP)^FINT). ADDED species (4-6,11-17,19): KT density Hamilton
#   (deferred — not in emt01). Reuses the shared self-thinning RDPSRT + snag booking.
# =============================================================================

const EM_PMSC  = Float32[5.45676, 5.26043, 5.55086, 0.0, 0.2118, 0.0, 3.87794, 6.41265, 5.88697, 5.58766,
                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.47709, 0.0]
const EM_PMD   = Float32[-0.0118233, -0.0097092, -0.0129121, 0.0, 0.0, 0.0, 0.3078, -0.0127328, -0.0333752,
                         -0.0052485, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0395156, 0.0]
const EM_PMDSQ = Float32[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0174, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

@inline _em_orig_species(sp::Int) = sp <= 3 || (7 <= sp <= 10) || sp == 18

# em/morts.f MORCON: IPDG/IPDG2(ITYPE,IFOR) → POT index → GMULT/REIN (added-species Hamilton). IFOR 1-6
# collapses to two mappings: forests 1-4 (Beaverhead/Custer/Deerlodge/Gallatin)=Bitterroot, 5-6
# (Helena/Lewis&Clark)=Lolo. POT[k]=0.05·k+0.20 (0.25..2.90), so no table needed.
const EM_IPDG_BITT  = Int[7,6,6,6,6,6,5,6,5,6,6,6,6,8,7,7,7,7,7,5,3,4,3,4,4,6,5,1,1,6]
const EM_IPDG_LOLO  = Int[7,7,7,7,6,7,5,6,7,6,6,7,7,9,9,8,9,8,8,5,3,5,4,5,5,7,5,2,4,7]
const EM_IPDG2_BITT = Int[30,29,29,29,28,28,27,31,27,27,28,31,32,32,31,31,32,31,31,25,23,24,23,24,24,27,26,18,16,27]
const EM_IPDG2_LOLO = Int[31,29,30,31,29,30,29,33,31,30,29,34,34,35,34,34,35,33,33,27,23,26,26,27,26,31,26,19,23,31]
@inline _em_pot(k::Int) = 0.05f0 * k + 0.20f0    # POT(k), em/morts.f DATA POT

# em/morts.f label-220 iterative linear-fn fit between the 55%/85% SDI lines → TN10 (target tree count at
# D10). IPATH2=true (came from T≤T55D0, T>T55D10) computes the line once (TEM=T) then goes to 230; IPATH2=false
# (55%<T≤85% at DIA0) Newton-iterates TREEIT (≤100) so exp(CEPT+SLP·ln(DIA0)) ≈ T. TN10 = exp(CEPT+SLP·ln(D10)),
# capped at T85D10.
function _em_tn10_iter(tt::Float32, dia0::Float32, d10::Float32, const_::Float32, pmsdil::Float32,
                       pmsdiu::Float32, t85d10::Float32, t55d0::Float32, ipath2::Bool)::Float32
    treeit = tt + 0.1f0 * tt
    slp = 0f0; cept = 0f0; knt = 1
    while true
        tem = ipath2 ? tt : treeit
        d55m = (log(tem) - log(pmsdil * const_)) / (-1.605f0)
        t55m = log(tem)
        d85m = d55m * 1.25f0
        while true                                    # em/morts.f label 221: bump D85M until SLP ≤ −0.5
            d85m > 5f0 && (d85m = 5f0); d85m < 0.125f0 && (d85m = 0.125f0)
            t85m = log(const_ * (exp(d85m)^(-1.605f0)) * pmsdiu)
            slp = (t85m - t55m) / (d85m - d55m)
            (slp > -0.5f0 && d85m < 5f0) ? (d85m += 0.1f0) : break
        end
        cept = t55m - slp * d55m
        (ipath2 || tt <= t55d0) && break              # GOTO 230 (no Newton for the IPATH=2 path)
        tprime = cept + slp * log(dia0)
        diff = tt - exp(tprime)
        (diff <= 5f0 && diff >= -5f0) && break
        treeit += 0.5f0 * diff; knt += 1
        knt > 100 && break
    end
    tn10 = exp(cept + slp * log(d10))
    tn10 >= t85d10 && (tn10 = t85d10)
    return tn10
end

function mortality!(s::StandState, ::EasternMontana; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    ba = p.basal_area
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    # grown-stand sums (morts.f): T (total tpa), DQ10 (QMD of DBH+DG), AVED (BA-weighted mean DBH)
    tt = 0f0; sd2sq = 0f0; sd0sq = 0f0; dsum = 0f0
    @inbounds for i in 1:n
        pr = t.tpa[i]; d = t.dbh[i]; sp = Int(t.species[i])
        bark = bark_ratio(bark_a, bark_b, sp, d)
        g = t.diam_growth[i] / bark
        sd2sq += pr * (d * d + 2f0 * d * g + g * g); sd0sq += pr * d * d; tt += pr; dsum += d * pr
    end
    tt < 1f-6 && return s
    dq10 = sqrt(sd2sq / tt)          # DQ10 = QMD of DBH+DG (post-growth)
    dq0  = sqrt(sd0sq / tt)          # DQ0  = QMD of DBH (pre-growth) = DIA0
    aved = dsum / tt
    # DIA0<0.3 reset (morts.f 374-376)
    if dq0 < 0.3f0; dq10 = 0.3f0 + dq10 - dq0; dq0 = 0.3f0; end
    # SDI self-thinning boundary (morts.f 455-485)
    sdimax = stand_sdimax(s)
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    pmsdil = p.pct_sdimax_mort_lo > 0f0 ? p.pct_sdimax_mort_lo : 0.55f0
    const_ = sdimax / 0.02483133f0
    tmd10 = const_ * dq10^(-1.605f0); tmd10 > 35000f0 && (tmd10 = 35000f0)
    tmd0  = const_ * dq0^(-1.605f0);  tmd0  > 35000f0 && (tmd0  = 35000f0)
    t85d10 = tmd10 * pmsdiu; t55d10 = tmd10 * pmsdil
    t85d0  = tmd0  * pmsdiu; t55d0  = tmd0  * pmsdil
    # TN10 target (morts.f 200-271): the SDI mature-stand-boundary self-thinning.
    local tn10::Float32
    if tt > t85d0
        tn10 = t85d10                                              # kill to the 85% line
    elseif tt > t55d0
        tn10 = abs(t85d0 - tt) <= 5f0 ? t85d10 :
               _em_tn10_iter(tt, dq0, dq10, const_, pmsdil, pmsdiu, t85d10, t55d0, false)
    elseif tt <= t55d10
        tn10 = tt                                                  # below 55% at both — hold (RN=0)
    else
        tn10 = _em_tn10_iter(tt, dq0, dq10, const_, pmsdil, pmsdiu, t85d10, t55d0, true)   # IPATH=2
    end
    tn10 > tt && (tn10 = tt); tn10 < 0.1f0 && (tn10 = 0f0)
    rn = 1f0 - (1f0 - (tt - tn10) / tt)^(1f0 / fint)
    # SDI-in-effect gate (morts.f:650-653): TEM = min(CONST·D10^−1.605, 35000)·PMSDIL. jl previously omitted the
    # 35000 cap here, so on ultra-dense sub-1" cohorts (tiny dq10 → uncapped TMD10 ≫ 35000) TEM ballooned to
    # ~106k > tt ⇒ jl wrongly fell to BACKGROUND mortality instead of the SDI self-thin ⇒ severe under-kill
    # (#137/#140 class: dense-cohort self-thin). tem == t55d10 (the already-capped TMD10·PMSDIL). #140.
    tem = t55d10
    # Added-species (Hamilton) stand values (em/morts.f:366-391 + MORCON): RZ from the BAMAX-limited BA10,
    # GMULT/REIN from IPDG/IPDG2[ITYPE,IFOR]. Only used by the added-species branch below.
    itype = Int(p.habitat_input); (itype < 1 || itype > 30) && (itype = 1)
    ifor = Int(p.forest_idx); (ifor < 1 || ifor > 6) && (ifor = 1)
    bamax = s.control.ba_max > 0f0 ? s.control.ba_max : ((1 <= itype <= 30) ? EM_BAMAXA[itype] : 0f0)
    bamax <= 0f0 && (bamax = 1f0)
    deltba = 0.005454154f0 * dq10 * dq10 * tt - ba
    ba10 = ba + (bamax - ba) / bamax * deltba
    tb = ba10 / (0.005454154f0 * dq10 * dq10)
    ttb = (tt - tb) / tt; ttb > 0.9999f0 && (ttb = 0.9999f0)
    rz = 1f0 - (1f0 - ttb)^0.1f0
    ipdg  = ifor <= 4 ? EM_IPDG_BITT  : EM_IPDG_LOLO
    ipdg2 = ifor <= 4 ? EM_IPDG2_BITT : EM_IPDG2_LOLO
    poten1 = _em_pot(ipdg[itype]);  gmult1 = 0.90f0 / poten1; rein1 = (1f0 - (poten1 / 20f0 + 1f0)^(-1.605f0)) / 0.06821f0
    poten2 = _em_pot(ipdg2[itype]); gmult2 = 2.50f0 / poten2; rein2 = (1f0 - (poten2 + 1f0)^(-1.605f0)) / 0.86610f0
    sqba = sqrt(ba); icyc1 = Int(s.control.cycle) == 0
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    @inbounds for i in 1:n
        sp = Int(t.species[i]); pr = t.tpa[i]; pr <= 0f0 && continue
        d = t.dbh[i]
        if _em_orig_species(sp)
            ri = 0.5f0 * (1f0 / (1f0 + exp(EM_PMSC[sp] + EM_PMD[sp] * d + EM_PMDSQ[sp] * d * d)))
            rip = rn
            (tt <= tem || rn <= 0f0) && (rip = ri)              # background when SDI not yet limiting
            rip > 1f0 && (rip = 1f0)
            wki = pr * (1f0 - (1f0 - rip)^fint)
            wki > pr && (wki = pr)
            sdimax < 5f0 && (wki = pr)
            killed[i] = wki
        else
            # ADDED species (4-6,11-17,19) — KT/IE Hamilton potential-mortality (em/morts.f:661-723).
            bark = bark_ratio(bark_a, bark_b, sp, d)
            reldbh = d / aved
            dd = d <= 0.5f0 ? 0.5f0 : d
            dgi = t.diam_growth[i]
            ip = d <= 5f0 ? 2 : 1
            gmult = ip == 1 ? gmult1 : gmult2
            wk1 = t.dg_prev[i]; oldfnt = 10f0                    # WK1 = previous cycle's applied DG
            dgt = wk1 / oldfnt
            d <= 1f0 && dgt < 0.05f0 && (dgt = 0.05f0)
            (1f0 < d <= 5f0) && dgt < 0.05f0 && (dgt = 0.05f0 * (5f0 - d) / 4f0)
            g = wk1 / (bark * oldfnt)
            wk1 / oldfnt < dgt && (g = dgt / bark)
            (icyc1 || wk1 == 0f0) && dgi > 0.5f0 && (g = dgi / (bark * 10f0))
            g = g * gmult
            rip = 2.76253f0 + 0.222310f0 * sqrt(dd) - 0.0460508f0 * sqba + 11.2007f0 * g -
                  0.554421f0 / dd + EM_PMSC[sp] + 0.246301f0 * reldbh + 6.07129f0 * g / dd
            rip > 70f0 && (rip = 70f0); rip < -70f0 && (rip = -70f0)
            rip = 1f0 / (1f0 + exp(rip))
            rip = rip * (ip == 1 ? rein1 : rein2)                # ·POTENT
            ripp = ba * rz
            ba <= bamax && (ripp += (bamax - ba) * rip)
            ripp /= bamax
            ripp < rip && (ripp = rip); ripp > 1f0 && (ripp = 1f0)
            # em/morts.f:717-723 species-group NI rate: LL(5) full, RM(6) 20%, others 60%.
            smult = sp == 5 ? 1f0 : (sp == 6 ? 0.2f0 : 0.6f0)
            wki = pr * (1f0 - (1f0 - ripp)^fint) * smult
            wki > pr && (wki = pr)
            killed[i] = wki
        end
    end
    # BAMAX residual-BA cap (em/morts.f BA-check, em/morts.f:885-940). When the user did not set BAMAX, residual
    # BA is capped at BAMAX = SDIMAX·0.5454154·PMSDIU (the SDI-DERIVED value: em/sitset.f:159 sets BAMAX=BAMAXA(ITYPE)
    # but does NOT set LBAMAX, so sdical.f:203 overwrites it with SDIMAX·0.5454154·PMSDIU every morts cycle — the
    # habitat BAMAXA only feeds SDIDEF→SDIMAX at setup). NOT the `bamax` (=EM_BAMAXA) used by the BADIST weighting
    # above. Scales every record's kill up by ADJFAC=(BANEW−BAMAX)/BADEAD, iterating ≤100× until residual BA ≤ BAMAX.
    # This own-copy had OMITTED it (shared southern/mortality.jl + NC/UT/TT have it). Inert when residual BA ≤ BAMAX.
    let bamax_cap = sdimax * 0.5454154f0 * pmsdiu
        if sdimax >= 5f0 && bamax_cap > 0f0
            for _ in 1:100
                banew = 0f0; badead = 0f0
                @inbounds for i in 1:n
                    d = t.dbh[i]; sp = Int(t.species[i])
                    bark = bark_ratio(bark_a, bark_b, sp, d)
                    g = t.diam_growth[i] / bark
                    ba_ = 0.0054542f0 * (d + g)^2
                    banew  += ba_ * (t.tpa[i] - killed[i])
                    badead += ba_ * killed[i]
                end
                ((banew - bamax_cap) > 1f0 && badead > 0f0) || break
                adjfac = (banew - bamax_cap) / badead
                @inbounds for i in 1:n
                    wki = killed[i] * (1f0 + adjfac)
                    wki > t.tpa[i] && (wki = t.tpa[i])
                    killed[i] = wki
                end
            end
        end
    end
    # Dwarf-mistletoe mortality (mismrt.f): MAX-combine per-tree DM kill into killed[] before snags/removal,
    # same as the shared N-Rockies path (southern/mortality.jl). This variant has its own mortality! so it
    # must be wired here; inert on stands with no DM ratings (dmr==0 ⇒ per-tree no-op).
    # FIXMORT (morts.f:781): forced-mortality override applied AFTER the BA-check, before the DM combine —
    # same as southern/mortality.jl:499. This variant has its own mortality! so it must be wired here;
    # inert unless a FIXMORT keyword scheduled events (apply_fixmort! returns early on empty).
    apply_fixmort!(s, killed, n, fint)
    _ie_mis_variant(s.variant) && ie_dm_mortality_combine!(killed, s, fint, n)
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end
