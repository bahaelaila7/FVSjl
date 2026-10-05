# =============================================================================
# crown_init.jl — cycle-0 (LSTART) crown initialisation shared across variants.
#
# Every variant's crown.f ends with the same `DO 79 I=IREC2,MAXTRE` block: each inventory standing-dead record
# with a missing crown (ICR<=0) is dubbed with the variant's own DUBSCR (a rejection-bounded main-stream BACHLO
# draw), re-expressed on the normal height for top-killed records, and bounded to [10,95]. The draws are consumed
# BEFORE DGDRIV's, so skipping them shifts every later DGSCOR/REGENT draw of the stand (measured on the PN and SO
# WRD fixtures). `dub_dead_crowns!` is that block; the variant supplies DUBSCR as `dub(i) -> CR fraction`.
# =============================================================================

"""
    dub_dead_crowns!(dub, s)

crown.f `DO 79 I=IREC2,MAXTRE` (all variants): dub missing crowns on the cycle-0 dead records. FVS files dead
records from MAXTRE downward in read order, so IREC2→MAXTRE is the REVERSE of jl's dead storage (t.n+1:t.n+ndead).
`dub(i)` returns the variant's ICRI BEFORE the shared top-kill/bounds tail — most variants `INT(CR*100+.5)` after
DUBSCR (use `icri_round`), but the crown-length forms (CR GEMCR, CI/EM/IE CL/H) truncate `INT(CR*100.)`.
"""
@inline icri_round(cr::Real) = trunc(Int, Float32(cr) * 100f0 + 0.5f0)
function dub_dead_crowns!(dub, s::StandState)
    t = s.trees; n = t.n
    t.ndead > 0 || return s
    @inbounds for i in (n + Int(t.ndead)):-1:(n + 1)
        Int(t.crown_pct[i]) > 0 && continue
        icri = Int(dub(i))
        if t.trunc[i] != 0
            hn = Float32(t.norm_ht[i]) / 100f0; hd = hn - Float32(t.trunc[i]) / 100f0
            cl = (Float32(icri) / 100f0) * hn - hd
            icri = trunc(Int, (cl * 100f0 / hn) + 0.5f0)
        end
        icri > 95 && (icri = 95); icri < 10 && (icri = 10)
        t.crown_pct[i] = Int32(icri)
    end
    return s
end

# bm/cratet.f LSTART crown-init: DENSE runs over the FULL inventory (live + HISTORY 6-9 standing-dead records),
# and that dead-inclusive density (BA/AVH/point-CCF) is what CROWN→DUBSCR sees when dubbing the D<1 / missing-CR
# LIVE trees. jl partitions the dead into t.n+1:t.n+ndead, so it computes the dead-inclusive scalars by temporarily
# extending the live range, then restores live-only (the grow cycle recomputes density before use). Measured on
# 504443988126144: dead-inclusive BA 55.56 / AVH 85.07 / TPCCF 84.3 (vs live-only 35/45/63) — matches live DUBSCR.
# SHARED (one implementation per Fortran mechanism): bm/, nc/ and pn/cratet.f run the identical LSTART sequence
# (RDPSRT(.FALSE.) IND → `LBKDEN=IDG.LT.2; CALL DENSE` over live+dead → CROWN), so BM, NC and PN all call this.
function crown_init_lstart_dead_inclusive!(s::StandState)
    t = s.trees
    nlive = t.n
    # bm/cratet.f:189-195 `LBKDEN = IDG.LT.2; CALL DENSE` — CROWN (cratet.f:610) dubs against THAT density, whose
    # live WK3 is BACKDATED to the start of the measured-growth period (dense.f:70-128: measured-DG trees by their
    # own increment, the rest by the stand-average BAGR) whenever ≥1 live record carries a measured increment
    # (SN>0; else WK3=DBH). Not backdating put the DUBSCR BA/TPCCF at the current-DBH values (1127576412290487:
    # 3 LP with past DBH ⇒ live BA 16.486 / TPCCF 80.04 vs jl 48.25 / 199.2 ⇒ every seedling crown mis-dubbed).
    # _backdate_dbh! is the shared IDG-faithful dense.f port (live 1:t.n, in place); restored below.
    lbkden = s.control.growth_idg < 2
    # base/fvs.f:193-196 — `SDICLS(0,0.,999.,1,SDIAC,SDIAC2,…)` runs IMMEDIATELY BEFORE `CALL CRATET`,
    # with the source comment "SDICLS IS CALLED HERE SO CROWNS WILL DUB CORRECTLY IN VARIANTS USING THE
    # WEIBULL DISTRIBUTION". So the LSTART dub's RELSDI = SDIAC/SDIDEF is NOT zero: it is the stand
    # Reineke SDI over the LIVE inventory at its READ diameters (before CRATET's backdating DENSE).
    # jl passed no crown_sdi ⇒ RELSDI=0 ⇒ ACRNEW = C0 (its maximum) ⇒ every Weibull-dubbed inventory
    # crown came out too high. MEASURED on the OC control stand (S248112, a missing-crown LP D=11.5):
    # live SDIAC 196.15 ⇒ ICR 76; jl with SDIAC=0 ⇒ 85.
    sdiac = stand_sdi_reineke(s)
    saved_live = lbkden ? t.dbh[1:nlive] : Float32[]
    # #151: dense.f:83-87 — in the CRATET backdating DENSE, standing-dead records get WK3=DBH EXCEPT IMC(I)==9
    # (HISTORY 8,9, older-dead) which LOAD DBH=0 ⇒ they add 0 to BA/CCF/SDI while their HEIGHT still counts
    # toward AVH (measured on 449747082489998: live DUBSCR AVH 67.34 vs jl 1.01 without the dead heights).
    # notre.f:119-124 expands EVERY inventory-dead record (IREC2..MAXTRE) with VP/FP/FP2×(FINT/FINTM), so the dead
    # PROB that CRATET's DENSE sums is TREE_COUNT×FINT/FINTM (BM 10/5 ⇒ ×2; FVS divides it back only where a true
    # density is needed — snag init cratet.f:565, dbstrls MortPA). jl's shared notre! leaves dead TPA unscaled, so
    # scale them here for this pass (24001521010900: 5 HISTORY-8 dead at PROB 12.036/1.998 ⇒ TPCCF 77.788→77.813).
    saved = Tuple{Int,Float32}[]
    saved_tpa = Float32[]
    if t.ndead > 0
        t.n = nlive + Int(t.ndead)
        fintr = s.control.growth_fintm > 0f0 ? s.control.growth_fint / s.control.growth_fintm : 1f0
        @inbounds for i in (nlive + 1):(nlive + Int(t.ndead)); push!(saved_tpa, t.tpa[i]); t.tpa[i] *= fintr; end
    end
    # AVHT40 top height from REAL DBH/HT over live + all dead (IND = real-DBH sort, not WK3), with the dead at their
    # notre-expanded PROB (24001521010900: the ×2 dead fill the top 40 TPA ⇒ AVH 48.466 = live, 41.596 unscaled).
    # dense.f:285-297 AVH walk of the :195 CRATET DENSE, over the :164-166 IND (bm_cratet166_ind: LNKCHN-seeded,
    # dead included, RDPSRT(.FALSE.) on real dbh) — the AVH the :610 CROWN→DUBSCR dub reads. The generic double
    # sort put different equal-DBH records at the 40-TPA boundary (23899355010900 AVH 72.567 vs live 71.606).
    # That DENSE runs BEFORE the missing-height dub (DO 130 :363, DO 145 :464) ⇒ heights are as READ (missing = 0):
    # use the pre-dub snapshot cratet_ht_in (41135212010497: dubbed 1.01 seedlings ⇒ AVH 33.671 vs live 33.160).
    ind153 = bm_cratet166_ind(s, view(t.dbh, 1:t.n), nlive, t.n)   # cratet.f:150-153 IND (real DBH, dead included)
    dbh_real = t.dbh[1:t.n]                   # real DBH of live + dead, before backdating/zeroing (PTBAL's WK5)
    avht_real = let ntot = t.n, ord = ind153, hin = s.calib.cratet_ht_in
        use_hin = length(hin) == ntot
        avh = 0f0; ssumn = 0f0
        @inbounds for k in 1:ntot
            ii = Int(ord[k]); p = t.tpa[ii]
            ssumn + p > 40f0 && (p = 40f0 - ssumn)
            ssumn += p; avh += (use_hin ? hin[ii] : t.height[ii]) * p
            ssumn >= 40f0 && break
        end
        ssumn > 0f0 ? avh / ssumn : 0f0
    end
    if lbkden                                 # backdate LIVE WK3 only (after the real-DBH AVH ranking)
        t.n = nlive; _backdate_dbh!(s); t.n = nlive + Int(t.ndead)
    end
    if t.ndead > 0
        @inbounds for i in (nlive + 1):(nlive + Int(t.ndead))
            (t.history[i] == 8 || t.history[i] == 9) || continue
            push!(saved, (i, t.dbh[i])); t.dbh[i] = 0f0
        end
    end
    # dense.f walks its BA/PCCF/CCF sums DO 50 ISPC / DO 10 I3 / I=IND1(I3), and SETUP's IND1 at CRATET is species-major
    # over ALL records in READ order — the dead still interleaved at their input positions (cratet.f deletes them
    # later). jl's IND1 (sort_key) files the dead after the live, so on a dead-bearing stand the Float32 sums ran in
    # another order (IE 3356357010690: calibration BA 2 ULP off ⇒ DO-220 DGF WK2 2-8 ULP ⇒ cycle-1 WK1/MORTS kill).
    # Key the pass on the read order, then restore.
    sk_saved = (t.ndead > 0 && length(s.calib.input_seq) == t.n) ? t.sort_key[1:t.n] : nothing
    sk_saved === nothing || @inbounds(for i in 1:t.n; t.sort_key[i] = Float64(s.calib.input_seq[i]); end)
    # CA/SO/WS CCFCAL is R5CRWD — height-dependent (H<4.5 ⇒ CRWDTH=SM·H) — and this DENSE runs before the missing-
    # height dub, so a missing height adds CCF 0 here: swap the as-read heights in for this pass and its RELDEN.
    ht_swap = _cratet_predub_ccf(s.variant) && length(s.calib.cratet_ht_in) == t.n
    ht_saved = ht_swap ? t.height[1:t.n] : Float32[]
    ht_swap && @inbounds(for i in 1:t.n; t.height[i] = s.calib.cratet_ht_in[i]; end)
    compute_density!(s)                    # CRATET DENSE: backdated live (+ dead-inclusive) BA / point-CCF
    # RELDEN of that DENSE: dense.f's RELDSP walk is species-major over the SAME read-order IND1 (dead interleaved), so
    # take it while the read-order keys are still in place (MEASURED FVSso_g16 449489561489998 REGENT-calibration
    # RELDEN 43445D43 live; jl's restored-key (dead-after-live) walk gave 43445D41 ⇒ PCTRED ⇒ HCOR 2 ULP).
    s.calib.cratet_relden = stand_ccf(s)
    sk_saved === nothing || @inbounds(for i in 1:t.n; t.sort_key[i] = sk_saved[i]; end)
    # dense.f:244 `CALL PCTILE(ITRN,IND,WK5,PCT,TOTAL)` — in the BACKDATING pass, PCT is accumulated over
    # **IND**, which cratet.f sorted on the REAL `DBH`, while the per-tree weight `WK5 = D*D*PROB` uses the
    # BACKDATED diameter (dense.f:184 `IF(LBKDEN.AND.LREDO) D = WK3(I)`). compute_density! above derived BOTH
    # from the backdated diameters, so the ORDER was wrong whenever backdating reshuffles near-equal trees —
    # and PCT feeds the PCR crown model directly (KT/IE/EM `b13*P + b14*log(P)`). MEASURED on the KT WRD
    # fixture with blanked crowns: live P rises monotonically with the read diameter (D 7.9→45.694,
    # 8.0→50.489, 8.2→54.573, 8.4→58.704) while jl's backdated order INVERTED it (50.92, 45.98, 32.92, 37.06).
    # Rebuild PCT with the real-diameter IND and the backdated weights.
    if lbkden && t.n > 0                   # _pctile! indexes idx[1] unguarded; an empty stand (no live
                                           # records AND no dead ones) reaches here with t.n == 0
        # That IND is cratet.f:150-153's `IND=IND1; RDPSRT(ITRN,DBH,IND,.FALSE.)` — the one this DENSE runs on —
        # NOT bm_cratet_ind!'s: with standing-dead records that switches to the identity .TRUE. re-sort of
        # cratet.f:257, which only happens AFTER the dead are deleted. The two differ only in how equal-DBH ties
        # fall, which permuted PCT inside every tie group (EM REGCAL fixture, 2 dead + 12-way ties: LL PCT
        # 1.0/7.3/15.5/36.4 vs live 0.1/6.4/19.8/32.5 ⇒ NIVAR BAL ⇒ SNX 568.53 vs 568.28).
        idx = view(s.scratch.stat_idx, 1:t.n); idx .= ind153
        _pctile!(t.crown_ratio, t, idx, t.n)
    end
    # (s.calib.cratet_relden — RELDEN after cratet.f:195 DENSE, backdated + dead-inclusive — set above → REGENT HCOR cal)
    ht_swap && @inbounds(for i in 1:t.n; t.height[i] = ht_saved[i]; end)
    # The rest of that DENSE's state for the EM LSTART REGCAL (em/cratet.f:553 — no DENSE in between): BA (=OLDBA,
    # backdated), AVH (the AVHT40 walk above), the per-point PCCF (dense.f:202 accumulates it only in the backdated
    # pass) and PCT (dense.f:244). RELDM1 is the dense.f:259 interpolation to the FINTH-year start; its "current"
    # term is the second, current-DBH pass, dead-inclusive with every dead record at its DBH (dense.f:184).
    c = s.calib
    c.cratet_ba = s.plot.basal_area; c.cratet_avh = avht_real
    c.cratet_pccf = copy(s.density.point_ccf); c.cratet_pct = t.crown_ratio[1:nlive]
    # The cycle-0 DEAD records keep this DENSE's PCT and PTBALT for good: later DENSEs run over ITRN=IREC1 (live only),
    # so dbstrls.f's dead rows (:308-440) report them. PTBAL (dense.f:280) runs after the SECOND, current-DBH pass:
    # WK5 = DBH·(DBH·PROB) at the READ diameter (dead PROB ×FINT/FINTM), per point in IND (real-DBH RDPSRT) order,
    # XBALT += WK5·.005454154·PI/GROSPC (ptbal.f:144-145).
    let ntot = t.n
        xb = zeros(Float32, MAXPLT); ptb = zeros(Float32, ntot)
        pif = s.plot.pi; gr = s.plot.gross_space
        @inbounds for k in 1:ntot
            ii = Int(ind153[k]); ip = Int(t.plot_id[ii]); (1 <= ip <= MAXPLT) || continue
            ptb[ii] = xb[ip]
            d = dbh_real[ii]
            xb[ip] = xb[ip] + (d * (d * t.tpa[ii])) * 0.005454154f0 * pif / gr
        end
        c.cratet_ptbaa = xb
        c.cratet_live_ptbal = ptb[1:min(nlive, ntot)]
        if t.ndead > 0
            nd = Int(t.ndead)
            c.cratet_dead_pct = t.crown_ratio[(nlive + 1):(nlive + nd)]
            # pctile.f with N=1 sets PERCNT(1) — array ELEMENT 1, an unused slot when the lone record is a dead one filed
            # at MAXTRE — and returns, so that dead record's PCT stays 0 (jl files it at index 1 and _pctile! gave it 100).
            (nlive == 0 && nd == 1) && (c.cratet_dead_pct[1] = 0f0)
            c.cratet_dead_ptbal = ptb[(nlive + 1):(nlive + nd)]
        else
            c.cratet_dead_pct = Float32[]; c.cratet_dead_ptbal = Float32[]
        end
    end
    @inbounds for (i, d) in saved; t.dbh[i] = d; end
    # dense.f:249-252 second pass: RMSQD = SQRT(TSUMD2/TPROB), TSUMD2 += D·(D·P), over IND1 species-major, current DBH
    let bk = lbkden ? t.dbh[1:nlive] : Float32[], tsumd2 = 0f0, tprob = 0f0, bat = 0f0, iseq = c.input_seq
        lbkden && @inbounds(for i in 1:nlive; t.dbh[i] = saved_live[i]; end)
        ord = length(iseq) == t.n ? sortperm(collect(1:t.n); by = j -> (Int(t.species[j]), iseq[j])) :
              sortperm(collect(1:t.n); by = j -> (Int(t.species[j]), j))
        @inbounds for j in ord
            pj = t.tpa[j]; dj = t.dbh[j]
            wk5 = dj * (dj * pj)
            tprob += pj; tsumd2 += wk5
            bat += 0.005454154f0 * wk5                 # dense.f BATREE = 0.005454154*WK5(I); BAT = BAT + BATREE
        end
        c.cratet_rmsqd = tprob > 0f0 ? sqrt(tsumd2 / tprob) : 0f0
        # dense.f:260 OLDBA = TEMP2 = (BA−OLDBA)·RAT + OLDBA: BA = this pass's BAT, OLDBA = the backdated pass's BAT.
        c.cratet_oldba = lbkden ? let fth = s.control.growth_finth > 0f0 ? s.control.growth_finth : 5f0,
                                      fit = s.control.growth_fint > 0f0 ? s.control.growth_fint : 10f0
            (bat - c.cratet_ba) * (fth / fit) + c.cratet_ba
        end : 0f0
        lbkden && @inbounds(for i in 1:nlive; t.dbh[i] = bk[i]; end)
    end
    c.cratet_reldm1 = c.cratet_relden
    if lbkden
        bk = t.dbh[1:nlive]
        @inbounds for i in 1:nlive; t.dbh[i] = saved_live[i]; end
        rcur = stand_ccf(s)
        @inbounds for i in 1:nlive; t.dbh[i] = bk[i]; end
        finth = s.control.growth_finth > 0f0 ? s.control.growth_finth : 5f0
        fint = s.control.growth_fint > 0f0 ? s.control.growth_fint : 10f0
        c.cratet_reldm1 = (rcur - c.cratet_relden) * (finth / fint) + c.cratet_relden
    end
    @inbounds for (k, i) in enumerate((nlive + 1):(nlive + length(saved_tpa))); t.tpa[i] = saved_tpa[k]; end
    t.n = nlive
    lbkden && @inbounds(for i in 1:nlive; t.dbh[i] = saved_live[i]; end)
    s.plot.avg_height = avht_real
    crown_ratio_update!(s, s.variant; lstart = true, crown_sdi = sdiac)   # DUBSCR-dub live D<1 seedlings + Weibull-dub missing-CR overstory
    compute_density!(s; cratet_ind = true)  # restore live-only density; CRATET IND ⇒ cycle-0 PCT/AVH (cratet.f:692 DENSE)
    return s
end


"""
    topkill_icri(t, i, icri) -> Int

crown.f statement 55, identical in every variant:

    55 IF (.NOT.LSTART .OR. ITRUNC(I).EQ.0) GO TO 59
       HN=REAL(NORMHT(I))/100.0
       HD=HN-REAL(ITRUNC(I))/100.0
       CL=(REAL(ICRI)/100.)*HN-HD
       ICRI=INT((CL*100./HN)+.5)

A top-killed INVENTORY record carries its crown ratio against the height it WOULD have had
(`NORMHT`), not the height that is left. The dead top `HD = HN − ITRUNC/100` is subtracted from the
crown length that `ICRI` implies on `HN`, and what remains is re-stated as a fraction of `HN`. Only
at LSTART — a tree that is topped later during the projection is handled by the CRMAX cap instead.

Caller supplies `icri` AFTER the 9052 rounding and BEFORE the statement-59 bounds. FVS divides by
`HN` unguarded because a record with `ITRUNC > 0` always has `NORMHT > 0` by then (cratet.f resolves
it in the missing-height dub); the guard here keeps a not-yet-resolved `NORMHT` from producing Inf.
"""
@inline function topkill_icri(t, i::Integer, icri::Integer)
    t.trunc[i] == 0 && return Int(icri)
    hn = Float32(t.norm_ht[i]) / 100f0
    hn <= 0f0 && return Int(icri)
    hd = hn - Float32(t.trunc[i]) / 100f0
    cl = (Float32(icri) / 100f0) * hn - hd
    return trunc(Int, (cl * 100f0 / hn) + 0.5f0)
end

"""
    crown_isort(s; lstart=false) -> Vector{Int32}

crown.f `DO 11 JJ=1,ITRN; ISORT(IND(JJ)) = ITRN-JJ+1` — the whole-stand descending-DBH RANK the
Weibull crown model draws its X from (largest ⇒ ITRN, smallest ⇒ 1). The KEY is always the CURRENT
`DBH(I)`, never DBH+DG:

* **cycling** — `gradd.f:177-186` calls `UPDATE` (which applies `DBH += DG/BRATIO`) and only THEN
  `RDPSRT(ITRN,DBH,IND,.TRUE.)`, so IND ranks the already-grown diameter. jl's apply-loop in
  `grow!` runs before `crown_ratio_update_fvs!` too and does NOT clear `t.diam_growth` (the
  FVS_TreeList DG column reads it), so `t.dbh` is already the grown value — adding `DG/BRATIO`
  again double-counts this cycle's growth.
* **LSTART** — `cratet.f` sorts on the READ diameter, and WHICH sort depends on the stand:
  `:188 IF(IREC2.EQ.MAXTP1) GO TO 60` skips the second sort when there are no dead records, leaving
  the `:153 RDPSRT(ITRN,DBH,IND,.FALSE.)` whose IND was PRESET from IND1 (species order) — a
  different tie-break among equal diameters than the sequential seed of the `:257
  RDPSRT(ITRN,DBH,IND,.TRUE.)` that a stand WITH dead records ends on. `bm_cratet_ind!` is that
  branch; it is variant-agnostic (every variant's cratet.f has the identical pair of calls).
"""
function crown_isort(s::StandState; lstart::Bool = false)
    t = s.trees; n = t.n
    idx = Vector{Int32}(undef, n)
    n == 0 && return idx
    if lstart
        bm_cratet_ind!(s, idx)
    else
        _rdpsrt!(view(t.dbh, 1:n), idx)
    end
    isort = Vector{Int32}(undef, n)
    @inbounds for jj in 1:n; isort[idx[jj]] = Int32(n - jj + 1); end
    return isort
end

"""
    point_crown_inputs(s) -> (prd, qmdplt, tpccf)

The per-inventory-point inputs crown.f computes for DUBSCR (2021 edits, identical in every variant): PRD = ZRD/XMAXPT
(SDICAL + SDICLS Zeide, live trees, 0 when XMAXPT≤0), QMDPLT = sqrt((PTBAA/PTPA)/0.005454) floored at 1 (1 when the
point has no trees), TPCCF = PCCF(ITRE). Returned as closures over the point index.
"""
function point_crown_inputs(s::StandState)
    xmaxpt, zrd, _ = point_zeide!(s)
    dens = s.density
    prd(pt) = (1 <= pt <= length(xmaxpt) && xmaxpt[pt] > 0f0) ? zrd[pt] / xmaxpt[pt] : 0f0
    function qmdplt(pt)
        baplt = (1 <= pt <= length(dens.point_ba)) ? dens.point_ba[pt] : 0f0
        tpaplt = (1 <= pt <= length(dens.point_tpa)) ? dens.point_tpa[pt] : 0f0
        q = tpaplt > 0f0 ? sqrt((baplt / tpaplt) / 0.005454f0) : 1f0
        q <= 1f0 ? 1f0 : q
    end
    tpccf(pt) = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0f0
    return prd, qmdplt, tpccf
end

# Variants whose CCFCAL crown width (R5CRWD: CRWDTH = SM·H below 4.5 ft / the spline branch) reads HT, so the
# CRATET backdating DENSE ({ca,so}/cratet.f:171, ws/cratet.f:247 — before the missing-height dub) sees a missing
# height as 0 ⇒ that record adds no CCF (CA FIA 302001653489998: four D=0.1 seedlings, live CCFT 0.00 each).
_cratet_predub_ccf(v) = v isa CentralCalifornia || v isa SouthCentralOregon || v isa WestSierra
