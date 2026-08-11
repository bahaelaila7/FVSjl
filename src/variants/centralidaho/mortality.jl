# =============================================================================
# mortality.jl (centralidaho) — CI mortality (ci/morts.f). Chunk 7. Full Hamilton RIP model
# (RIP=2.76253+0.222310·√D−0.0460508·√BA+11.2007·G−0.554421/DD+PMSC+0.246301·RELDBH+6.07129·G/DD)
# + POTENT reinforcement (POT/IPDG) + BA10 self-thin — IDENTICAL structure to IE (ci/morts.f:269/329).
# CI deltas: PMSC/PMD/POT/IPDG/IPDG2 from ci/morts.f; ITYPE=NIHMAP[ICINDX]; IFOR 1-6; bark=ci_bratio;
# BAMAX from CI_BAMAXA[ICINDX] (Zeide SDImax rides in via stand_sdimax for the size cap gate).
# =============================================================================

let
    d = CI_DATADIR
    m = readlines(joinpath(d, "mort_1d.csv"))
    global const CI_MORT_PMSC = Float32[parse(Float32, split(strip(m[1+i]), ',')[1]) for i in 1:19]
    global const CI_MORT_PMD  = Float32[parse(Float32, split(strip(m[1+i]), ',')[2]) for i in 1:19]
    global const CI_MORT_POT = Float32[parse(Float32, strip(l)) for l in readlines(joinpath(d, "mort_pot.csv"))]
    # each CSV line = forest col j (1..11), 30 itype values ⇒ hcat gives [itype(30), ifor(11)]
    readmat(f) = (ls = readlines(joinpath(d, f)); hcat([Int[parse(Int, strip(x)) for x in split(strip(l), ',')] for l in ls]...))
    global const CI_MORT_IPDG  = readmat("mort_ipdg.csv")     # [30, 11]
    global const CI_MORT_IPDG2 = readmat("mort_ipdg2.csv")    # [30, 11]
end
const CI_MORT_MAPFOR = Int[10, 1, 1, 10, 1, 1]               # ci/morts.f:198 DATA MAPFOR — IFOR→MIFOR (IPDG column)

function mortality!(s::StandState, ::CentralIdaho; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    ba = p.basal_area
    icindx = Int(p.habitat_input)
    itype = (1 <= icindx <= 130) ? Int(CI_NIHMAP[icindx]) : 1
    (itype < 1 || itype > 30) && (itype = 1)
    sdimax = stand_sdimax(s)
    # CI (Zeide) BAMAX = SDIMAX·0.5454154·PMSDIU (ci/morts.f — verified live 265.15 = 571.92·0.5454154·0.85),
    # NOT the site BAMAXA. PMSDIU default 0.85 (fraction).
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi / 100f0 : 0.85f0
    bamax = s.control.ba_max > 0f0 ? s.control.ba_max : sdimax * 0.5454154f0 * pmsdiu
    bamax <= 0f0 && (bamax = 1f0)
    tt = 0f0; sd2sq = 0f0; dsum = 0f0; wprob = 0f0
    @inbounds for i in 1:n
        pr = t.tpa[i]; d = t.dbh[i]; sp = Int(t.species[i])
        bark = ci_bratio(s.coef.species, sp, d)
        g = t.diam_growth[i] / bark
        sd2sq += pr * (d * d + 2f0 * d * g + g * g); tt += pr
        wprob += pr; dsum += d * pr
    end
    tt < 1f-6 && return s
    dq10 = sqrt(sd2sq / tt)
    deltba = 0.005454154f0 * dq10 * dq10 * tt - ba
    ba10 = ba + (bamax - ba) / bamax * deltba
    tb = ba10 / (0.005454154f0 * dq10 * dq10)
    ttb = (tt - tb) / tt; ttb > 0.9999f0 && (ttb = 0.9999f0)
    rz = 1f0 - (1f0 - ttb)^0.1f0
    aved = dsum / wprob
    ifor = Int(p.forest_idx); (ifor < 1 || ifor > 6) && (ifor = 1)
    mifor = CI_MORT_MAPFOR[ifor]                              # ci/morts.f:721 MIFOR=MAPFOR(IFOR)
    poten1 = CI_MORT_POT[CI_MORT_IPDG[itype, mifor]]
    poten2 = CI_MORT_POT[CI_MORT_IPDG2[itype, mifor]]
    gmult1 = 0.90f0 / poten1; rein1 = (1f0 - (poten1 / 20f0 + 1f0)^(-1.605f0)) / 0.06821f0
    gmult2 = 2.50f0 / poten2; rein2 = (1f0 - (poten2 + 1f0)^(-1.605f0)) / 0.86610f0
    sqba = sqrt(ba)
    icyc1 = Int(s.control.cycle) == 0
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    sc = s.control.sp_size_cap
    @inbounds for i in 1:n
        sp = Int(t.species[i]); pr = t.tpa[i]; pr <= 0f0 && continue
        d = t.dbh[i]; bark = ci_bratio(s.coef.species, sp, d)
        reldbh = d / aved
        dd = d <= 0.5f0 ? 0.5f0 : d
        dgi = t.diam_growth[i]
        ip = d <= 5f0 ? 2 : 1
        gmult = ip == 1 ? gmult1 : gmult2
        wk1 = t.dg_prev[i]; oldfnt = 10f0
        dgt = wk1 / oldfnt
        d <= 1f0 && dgt < 0.05f0 && (dgt = 0.05f0)
        (1f0 < d <= 5f0) && dgt < 0.05f0 && (dgt = 0.05f0 * (5f0 - d) / 4f0)
        g = wk1 / (bark * oldfnt)
        wk1 / oldfnt < dgt && (g = dgt / bark)
        (icyc1 || wk1 == 0f0) && dgi > 0.5f0 && (g = dgi / (bark * 10f0))
        g = g * gmult
        rip = 2.76253f0 + 0.222310f0 * sqrt(dd) - 0.0460508f0 * sqba + 11.2007f0 * g -
              0.554421f0 / dd + CI_MORT_PMSC[sp] + 0.246301f0 * reldbh + 6.07129f0 * g / dd
        rip > 70f0 && (rip = 70f0); rip < -70f0 && (rip = -70f0)
        rip = 1f0 / (1f0 + exp(rip))
        rip = rip * (ip == 1 ? rein1 : rein2)
        ripp = ba * rz
        ba <= bamax && (ripp += (bamax - ba) * rip)
        ripp /= bamax
        ripp < rip && (ripp = rip); ripp > 1f0 && (ripp = 1f0)
        # ci/morts.f species multiplier: 60% of NI rate for WB/PY/AS/MC/LM/CW/OH, 20% for WJ (juniper)
        smult = (sp == 11 || sp == 12 || sp == 13 || sp == 15 || sp == 16 || sp == 17 || sp == 19) ? 0.6f0 :
                (sp == 14 ? 0.2f0 : 1.0f0)
        wki = pr * (1f0 - (1f0 - ripp)^fint) * smult
        gsc = (dgi / bark) * (fint / 10f0)
        if (d + gsc) >= sc[sp, 1] && trunc(Int, sc[sp, 3]) != 1
            wki = max(wki, pr * sc[sp, 2] * fint / 10f0)
        end
        wki > pr && (wki = pr)
        sdimax < 5f0 && (wki = pr)
        killed[i] = wki
    end
    # Dwarf-mistletoe mortality (mismrt.f): MAX-combine the per-tree DM kill into killed[] before
    # snags/removal, exactly as the shared N-Rockies path does (southern/mortality.jl). CI has its own
    # mortality! so this must be wired here too; inert on stands with no DM ratings (dmr==0 ⇒ no-op).
    # FIXMORT (morts.f:781): forced-mortality override applied AFTER the BA-check, before the DM combine —
    # same as southern/mortality.jl:499. This variant has its own mortality! so it must be wired here;
    # inert unless a FIXMORT keyword scheduled events (apply_fixmort! returns early on empty).
    apply_fixmort!(s, killed, n, fint)
    _ie_mis_variant(s.variant) && ie_dm_mortality_combine!(killed, s, fint, n)
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end
