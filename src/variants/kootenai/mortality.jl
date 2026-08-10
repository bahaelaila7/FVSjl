# KT mortality (kt/morts.f) — Hamilton model coefficients (chunk 7).
const KT_MORT_POT  = Float32[0.25f0, 0.3f0, 0.35f0, 0.4f0, 0.45f0, 0.5f0, 0.55f0, 0.6f0, 0.65f0, 0.7f0, 0.75f0, 0.8f0, 0.85f0, 0.9f0, 0.95f0, 1f0, 1.05f0, 1.1f0, 1.15f0, 1.2f0, 1.25f0, 1.3f0, 1.35f0, 1.4f0, 1.45f0, 1.5f0, 1.55f0, 1.6f0, 1.65f0, 1.7f0, 1.75f0, 1.8f0, 1.85f0, 1.9f0, 1.95f0, 2f0, 2.05f0, 2.1f0, 2.15f0, 2.2f0, 2.25f0, 2.3f0, 2.35f0, 2.4f0, 2.45f0, 2.5f0, 2.55f0, 2.6f0, 2.65f0, 2.7f0, 2.75f0, 2.8f0, 2.85f0, 2.9f0]  # POT(54) 0.25..2.90
const KT_MORT_PMSC = Float32[0f0, -0.17603f0, 0.317888f0, 0.317888f0, 0.607725f0, 1.57976f0, -0.12057f0, 0.94019f0, 0.2118f0, 0.2118f0, 0f0]  # per-species mortality constant B5
const KT_MORT_IPDG  = Int32[
    7 15 14 12 10 9 13 9 7 15 15;
    6 15 14 12 10 9 12 9 7 15 15;
    6 15 14 12 10 9 12 9 7 15 15;
    6 15 14 12 10 9 12 9 7 15 15;
    6 14 13 11 9 9 12 8 6 14 14;
    6 14 13 11 9 9 12 9 7 14 14;
    5 12 11 10 8 7 10 7 5 9 12;
    6 14 14 10 9 8 12 8 6 13 14;
    5 13 12 10 8 8 11 7 7 13 13;
    6 14 14 11 9 8 11 8 6 14 14;
    6 14 14 11 10 9 11 8 6 14 14;
    6 15 14 12 10 9 12 8 7 13 15;
    6 15 14 12 10 9 12 9 7 14 15;
    8 17 17 14 12 11 15 11 9 18 17;
    7 15 14 14 11 10 13 10 9 14 15;
    7 15 14 14 11 10 13 10 8 14 15;
    7 16 15 13 11 10 14 10 9 14 16;
    7 15 15 13 11 10 13 10 8 14 15;
    7 15 15 13 11 11 13 10 8 14 15;
    5 11 11 9 8 7 9 7 5 10 11;
    3 10 9 6 6 5 8 5 3 9 10;
    4 11 10 8 7 6 9 6 5 9 11;
    3 9 8 7 5 4 7 5 4 8 9;
    4 12 11 9 7 7 10 7 5 11 12;
    4 12 11 10 7 7 10 7 5 11 12;
    6 15 15 12 10 9 12 9 7 15 15;
    5 12 11 10 8 8 10 8 5 9 12;
    1 7 5 4 6 3 4 3 2 6 7;
    1 10 10 8 7 6 8 6 4 8 10;
    6 15 15 12 10 9 12 9 7 15 15]  # [itype 1..30, ifor 1..11] -> POT idx (D>5, IP=1)
const KT_MORT_IPDG2 = Int32[
    30 50 45 41 38 38 41 38 31 45 50;
    29 49 45 41 37 37 41 37 29 45 49;
    29 49 45 41 37 37 41 37 30 45 49;
    29 49 45 41 37 37 41 37 31 45 49;
    28 48 44 40 36 36 40 37 29 44 48;
    28 48 44 40 37 36 40 37 30 44 48;
    27 47 43 39 37 36 39 37 29 43 47;
    31 52 48 42 40 38 45 41 33 48 52;
    27 46 43 39 35 36 37 33 31 41 46;
    27 47 43 39 36 35 38 37 30 42 47;
    28 48 44 40 37 37 39 37 29 43 48;
    31 53 49 44 41 41 43 41 34 47 53;
    32 54 50 45 41 41 44 41 34 48 54;
    32 54 50 45 42 41 46 42 35 50 54;
    31 52 49 45 40 40 43 39 34 47 52;
    31 52 49 45 40 40 43 39 34 47 52;
    32 54 49 45 41 41 45 41 35 46 54;
    31 50 47 43 40 40 43 40 33 44 50;
    31 51 47 44 40 40 43 40 33 44 51;
    25 41 39 36 33 33 35 33 27 36 41;
    23 39 36 33 30 32 33 32 23 34 39;
    24 44 40 36 33 32 35 33 26 37 44;
    23 43 39 36 31 32 35 33 26 37 43;
    24 44 40 37 32 34 36 34 27 38 44;
    24 44 39 38 31 33 36 33 26 37 44;
    27 50 45 41 36 38 41 38 31 43 50;
    26 53 41 38 34 35 37 35 26 34 53;
    18 36 33 30 32 26 30 28 19 30 36;
    16 37 34 31 27 28 32 32 23 33 37;
    27 50 45 41 36 38 41 38 31 43 50]  # [itype, ifor] -> POT idx (D<=5, IP=2)
"""
    mortality!(s, ::Kootenai; fint, book_snags)

KT HAMILTON mortality (kt/morts.f). Stand-setup DQ10/RZ/AVED + per-species MORCON (POTEN→GMULT/REIN via
IPDG/IPDG2/POT) + per-tree RIP (logistic) → RIPP (BAMAX/RZ weighting) → WKI = P·(1−(1−RIPP)^FINT). Fills the
shared killed[] buffer, then reuses the shared apply-tail (book_mortality_snags! + TPA removal). NOT VARMRT.
(Establishment/FIXMORT/MORTMULT windows omitted — default X=1; cycle>1 uses the projected DG for the growth
term G since the past-DG WK1 isn't threaded yet — refine when the full-cycle differential needs it.)
"""
function mortality!(s::StandState, ::Kootenai; fint::Float32 = 10.0f0, book_snags::Bool = true)
    p, t = s.plot, s.trees
    n = t.n; n == 0 && return s
    ba = p.basal_area
    itype = Int(p.habitat_input)
    bamax = s.control.ba_max > 0f0 ? s.control.ba_max : ((1 <= itype <= 30) ? KT_BAMAXA[itype] : 0f0)
    bamax <= 0f0 && (bamax = 1f0)
    sdimax = stand_sdimax(s)
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    # stand sums (morts.f:196-210, RAW record order): T, SD2SQ → DQ10; AVED (BA-weighted mean DBH)
    tt = 0f0; sd2sq = 0f0; dsum = 0f0; wprob = 0f0
    @inbounds for i in 1:n
        pr = t.tpa[i]; d = t.dbh[i]; sp = Int(t.species[i])
        bark = bark_ratio(bark_a, bark_b, sp, d)
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
    # MORCON: POTEN → GMULT/REIN per size class (morts.f:646-653)
    ifor = kt_ifor(p)
    (ifor < 1 || ifor > 11) && (ifor = 8)
    it = (1 <= itype <= 30) ? itype : 1
    poten1 = KT_MORT_POT[KT_MORT_IPDG[it, ifor]]
    poten2 = KT_MORT_POT[KT_MORT_IPDG2[it, ifor]]
    gmult1 = 0.90f0 / poten1; rein1 = (1f0 - (poten1 / 20f0 + 1f0)^(-1.605f0)) / 0.06821f0
    gmult2 = 2.50f0 / poten2; rein2 = (1f0 - (poten2 + 1f0)^(-1.605f0)) / 0.86610f0
    sqba = sqrt(ba)
    icyc1 = Int(s.control.cycle) == 0
    killed = @view s.scratch.mort_killed[1:n]; fill!(killed, 0f0)
    sc = s.control.sp_size_cap
    @inbounds for i in 1:n
        sp = Int(t.species[i]); pr = t.tpa[i]; pr <= 0f0 && continue
        d = t.dbh[i]; bark = bark_ratio(bark_a, bark_b, sp, d)
        reldbh = d / aved
        dd = d <= 0.5f0 ? 0.5f0 : d
        dgi = t.diam_growth[i]
        ip = d <= 5f0 ? 2 : 1
        gmult = ip == 1 ? gmult1 : gmult2
        # G growth term (morts.f:259-272). WK1 = previous cycle's applied DG (vigor proxy), OLDFNT = its
        # period (10 for uniform cycles). Cycle 1 (WK1=0) falls to the DG override.
        wk1 = t.dg_prev[i]; oldfnt = 10f0
        dgt = wk1 / oldfnt
        d <= 1f0 && dgt < 0.05f0 && (dgt = 0.05f0)
        (1f0 < d <= 5f0) && dgt < 0.05f0 && (dgt = 0.05f0 * (5f0 - d) / 4f0)
        g = wk1 / (bark * oldfnt)
        wk1 / oldfnt < dgt && (g = dgt / bark)
        (icyc1 || wk1 == 0f0) && dgi > 0.5f0 && (g = dgi / (bark * 10f0))
        g = g * gmult
        rip = 2.76253f0 + 0.222310f0 * sqrt(dd) - 0.0460508f0 * sqba + 11.2007f0 * g -
              0.554421f0 / dd + KT_MORT_PMSC[sp] + 0.246301f0 * reldbh + 6.07129f0 * g / dd
        rip > 70f0 && (rip = 70f0); rip < -70f0 && (rip = -70f0)
        rip = 1f0 / (1f0 + exp(rip))
        rip = rip * (ip == 1 ? rein1 : rein2)                # ·POTENT
        ripp = ba * rz
        ba <= bamax && (ripp += (bamax - ba) * rip)
        ripp /= bamax
        ripp < rip && (ripp = rip); ripp > 1f0 && (ripp = 1f0)
        wki = pr * (1f0 - (1f0 - ripp)^fint)                 # X=1 (no MORTMULT window)
        gsc = (dgi / bark) * (fint / 10f0)
        if (d + gsc) >= sc[sp, 1] && trunc(Int, sc[sp, 3]) != 1
            wki = max(wki, pr * sc[sp, 2] * fint / 10f0)
        end
        wki > pr && (wki = pr)
        sdimax < 5f0 && (wki = pr)
        killed[i] = wki
    end
    # Dwarf-mistletoe mortality (mismrt.f): MAX-combine per-tree DM kill into killed[] before snags/removal,
    # same as the shared N-Rockies path (southern/mortality.jl). This variant has its own mortality! so it
    # must be wired here; inert on stands with no DM ratings (dmr==0 ⇒ per-tree no-op).
    _ie_mis_variant(s.variant) && ie_dm_mortality_combine!(killed, s, fint, n)
    book_snags && book_mortality_snags!(s, killed, n, fint)
    @inbounds for i in 1:n; t.tpa[i] = max(0f0, t.tpa[i] - killed[i]); end
    return s
end
