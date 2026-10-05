# =============================================================================
# regent.jl (centralidaho) — CI small-tree growth (ci/regent.f). Chunk 6. Overrides DG/HTG for
# small trees (D < XMAX). NIVAR path (CI conifers): per-subcycle height via the HTGRL regression,
# accumulated, then DBH from HTDBH; DG = new−old DBH. Coeffs from data/centralidaho/regent_coeffs_1d.csv.
#   HTGRL = CNST + HTBA·RELHT·PTBAA + PBAL·PTBALI + STBA·BA + RLHT·RELHT + CRSQ·RCR² + CR·RCR
#         + BAL·TBAL + HDM1·RHDM1 + HDM2·RHDM2 + PTBA·PTBAA   (ci/regent.f:641)
#   CON = RHCON·exp(HCOR)  (RHCON≈1, HCOR = c.htg_cor_small); H2 = H1 + HTGRL·SCALE·XRHGRO·CON
#   DK = exp(DHCN + DHHT·ln(HK) + DHCR·ln(RCR))  (HTDBH); DG = DK − D_old
# SMHTGF (stochastic BETA/FINDAG) is the establishment/aspen path — deferred; cit01 small trees are conifers.
# =============================================================================

let
    path = joinpath(CI_DATADIR, "regent_coeffs_1d.csv")
    hdr = split(strip(readlines(path)[1]), ',')
    rows = [split(strip(l), ',') for l in readlines(path)[2:end]]
    ci(n) = findfirst(==(n), hdr)
    col(n) = Float32[parse(Float32, strip(rows[sp][ci(n)])) for sp in 1:19]
    for n in ("CNST","CR","CRSQ","DGMAX","DHCN","DHCR","DHHT","HDM1","HDM2","HTBA","RLHT","STBA",
              "XMAX","XMIN","BAL","PBAL","PTBA")
        @eval global const $(Symbol("CI_RG_" * n)) = $(col(n))
    end
end
const CI_RG_REGYR = 5.0f0
const CI_RG_DIAM = Float32[0.4,0.3,0.3,0.3,0.2,0.2,0.4,0.3,0.3,0.5,0.3,0.3,0.2,0.3,0.2,0.3,0.2,0.2,0.2]  # ci/regent.f DATA DIAM (min DBH)
const CI_RG_AB = Float32[1.11436, -0.011493, 0.43012f-4, -0.72221f-7, 0.5607f-10, -0.1641f-13]  # PCTRED poly (regent.f:448)
# BKPT = DBH BREAKPOINT (ci/regent.f:490-500,816-822 SELECT CASE): the DBH at/above which the small-tree HT-DBH
# diameter increment is NOT applied (ci/regent.f:993 `IF(D .GE. BKPT) GO TO 23`) — for BKPT≤D<XMAX the tree still
# takes the XWT-blended regent HEIGHT growth but KEEPS its large-tree dgf DG. CASE(17,19) CW/OH=1.0; CASE(14,15)
# juniper/MC=99.0 (fully regent, XMAX=99); DEFAULT (CIVAR conifers + aspen sp13)=3.0. jl previously overwrote the
# large-tree DG for ALL D<XMAX (conifers XMAX=5 ⇒ the [3,5) medium-conifer band lost its large-tree DG; QMD/BA
# deficit as trees mature past 3"), mirroring the UT bug (see src/variants/utah/regent.jl UT_RG_BREAK).
const CI_RG_BREAK = Float32[3,3,3,3,3,3,3,3,3,3, 3,3,3, 99,99, 3, 1, 3, 1]
@inline _ci_ut_species(sp::Int) = sp == 13 || sp == 14 || sp == 15 || sp == 17 || sp == 19  # ci/regent.f UTVAR branch

"""ci/regent.f sp17/19 (CW/OH) Curtis-Arney height→DBH inverse (P2=1709.7229,P3=5.8887,P4=−0.2286)."""
@inline function _ci_ut_htdbh(ht::Float32)::Float32
    p2 = 1709.7229f0; p3 = 5.8887f0; p4 = -0.2286f0
    hat3 = 4.5f0 + p2 * fexp(-p3 * fpow(3.0f0, p4))
    ht >= hat3 ? fexp(flog((flog(ht - 4.5f0) - flog(p2)) / (-p3)) * 1.0f0 / p4) :   # regent.f:1059 ALOG(..)*1./P4
                 ((ht - 4.51f0) * 2.7f0 / (hat3 - 4.51f0)) + 0.3f0
end

# ci/blkdat.f DATA HT1/HT2 — the Wykoff intercept/slope REGENT's UTVAR 13/17/19 diameter dub inverts (AX=HT1 when
# IABFLG=1, else the cratet-fitted AA). (The species CSV ht1/ht2 columns are NOT CI's blkdat values for every species.)
# CI_BLK_HT1/HT2 (ci/blkdat.f:256-265 Wykoff HT1/HT2) are defined once in centralidaho/site_index.jl.
const CI_RG_BACON = 0.005454154f0
@inline _ci_ttvar(sp::Int) = sp == 11 || sp == 12 || sp == 16            # ci/regent.f:508 CASE(11,12,16) TTVAR
@inline _ci_bkpt(sp::Int) = (sp == 17 || sp == 19) ? 1.0f0 : (sp == 14 || sp == 15) ? 99.0f0 : 3.0f0   # :490-500

# ci/smhtgf.f for the TTVAR species (11,12,16): the persistent ZRAND (drawn whenever it is -999 — jl's inventory
# default 0 too — no DGSD gate), BETA1/BETA2 from the subcycle point CCF, reset to -999 on a ≤0.1 increment.
function _ci_smhtgf!(s::StandState, i::Int, cri::Float32, tpccf::Float32)::Float32
    t = s.trees
    zr = t.tree_random[i]
    if zr == 0f0 || zr == -999f0
        while true; zr = bachlo(s.rng, 0f0, 1f0); (-2f0 <= zr <= 2f0) && break; end
        t.tree_random[i] = zr
    end
    t.dbh[i] <= 0f0 && return 0f0                                  # :24-27 D≤0 ⇒ 0, no floor/reset
    beta1 = fexp(1.17527f0 - 0.42124f0 * flog(tpccf))
    beta2 = fexp(-2.56002f0 - 0.58642f0 * flog(tpccf))
    htg1 = beta1 + beta2 * cri
    stddev = htg1 * (1.08720f0 - 0.00230f0 * cri)
    htg = htg1 + zr * stddev
    htg > 0.1f0 && return htg
    t.tree_random[i] = -999f0
    return 0.1f0
end

@inline function _ci_zzran(s::StandState, dgsd::Float32)::Float32
    z = 0f0
    if dgsd >= 1.0f0
        while true; z = bachlo(s.rng, 0f0, 1f0); (z > 0.5f0 || z < -2.0f0) || break; end
    end
    return z
end

"""
    small_tree_growth!(s, stash, ::CentralIdaho; fint)

ci/regent.f (growth call, LESTB=.FALSE.) in Fortran shape:

* DO 17 J=1,NPER subcycles (KPER, NPER=ceil(NTYR/5)), walking IND1 (DO 16 ISPC / DO 15 I3). CIVAR (1-10,18) and UTVAR
  (13,14,15,17,19) species run at J=1 only with SCALE=NTYR/REGYR resp. NTYR/YR; the TTVAR species (11,12,16) run
  every subcycle through ci/smhtgf.f (persistent ZRAND) with TPCCF=PCCF·PPCCF ∈[25,300], SCALE=KPER(J)/REGYR and the
  climate WK4. Each record then feeds the next subcycle's BANEXT/RDNEXT (UPDATE DENSITY, :740-780) with its projected
  D2; the large-tree projection (:266-289) carries CCFCAL's ·P.
* DO 30 / DO 25 (:803-1300): per record K (the central, then the two tripled copies — each with a fresh ZZRAN for
  CIVAR/UTVAR; TTVAR copies take the central's values), HTGR = HTGR1 + ZZRAN term, floored 0.1, blended by XWT with
  K's OWN large-tree HTG (ci_triple_htg! for the NI copies), SIZCAP; below BKPT the diameter dub (CIVAR HT-DBH,
  UTVAR juniper/MC/Wykoff-inverse, TTVAR WK5 vs DKK) on the DDS scale, the D+0.001·HK direct set when HK<4.5, the
  DIAM floor and DGBND.
"""
small_tree_growth!(s::StandState, stash, ::CentralIdaho; fint::Float32 = 10.0f0) = _ci_regent!(s, stash; fint)

# ci/regent.f REGENT(LESTB, ITRNIN). The cycling call (LESTB=.FALSE.) and the ESGENT birth-cycle call (LESTB=.TRUE., from
# ci/esgent.f) are ONE routine; the LESTB differences are marked inline (NTYR=FINT−5/LSKIPH, interpolated BANEXT/RDNEXT
# from the post-thin TEMBA/TEMCCF, the DO-13 crown draw of the new records in STORAGE order, PCTRED/DELMAX on TEMCCF/
# TEMAHT, records I<ITRNIN skipped, PTBALI=BA·(1−PCT/100), aspen SITAGE=ABIRTH, XWT=0, no tripling, and the LESTB DBH
# assignment DBH=DK (≥DIAM)+0.001·HK, DG=DBH). `ba_in/relden_in/avh_in/pccf_in/ptba_in` = the gradd.f DENSE the routine
# reads (post-growth, pre-regen for ESGENT); `atba/atccf/atavh` = grincr.f:316-320's post-thin ATBA/ATCCF/ATAVH.
function _ci_regent!(s::StandState, stash; fint::Float32 = 10.0f0, lestb::Bool = false, nstart::Int = 0,
                     ba_in::Float32 = -1f0, relden_in::Float32 = -1f0, avh_in::Float32 = -1f0,
                     pccf_in::Vector{Float32} = Float32[], ptba_in::Vector{Float32} = Float32[],
                     atba::Float32 = -1f0, atccf::Float32 = -1f0, atavh::Float32 = -1f0)
    p, t, c, dens, ctl = s.plot, s.trees, s.calib, s.density, s.control
    sd = s.coef.species
    n = t.n; n == 0 && return s
    lestb && nstart >= n && return s
    cw = lestb ? nothing : clim_wk4(s, Float32(current_cycle_year(s)) + fint / 2f0)   # CLGMULT WK4 (ci/regent.f:709/947); nothing ⇒ 1
    wk4(i) = cw === nothing ? 1.0f0 : cw[i]
    ba = ba_in >= 0f0 ? ba_in : p.basal_area
    relden = relden_in >= 0f0 ? relden_in : p.relative_density
    avh = avh_in >= 0f0 ? avh_in : p.avg_height
    atavh_e = lestb ? (atavh >= 0f0 ? atavh : avh) : avh     # cycling: grincr.f:318 ATAVH=AVH
    kodtyp = Int(p.habitat_code)
    rhdm1 = (500 <= kodtyp < 600) ? 1.0f0 : 0.0f0
    rhdm2 = (600 <= kodtyp < 700) ? 1.0f0 : 0.0f0
    regyr = CI_RG_REGYR; yr = 10.0f0
    ntyr = trunc(Int, fint)                                   # :283 NTYR=INT(FINT)
    lestb && (ntyr -= 5)                                      # :300 IF(LESTB) NTYR=NTYR-5
    lskiph = lestb && ntyr <= 0                               # :301
    iyr = trunc(Int, regyr)
    nper = lskiph ? 0 : ntyr ÷ iyr; (!lskiph && ntyr % iyr != 0) && (nper += 1); (!lestb && nper < 1) && (nper = 1)
    kper = zeros(Int, max(nper, 1)); itot = ntyr; nn = nper
    @inbounds for k in 1:nper
        if nn == 1; kper[k] = itot; break; end
        kper[k] = itot ÷ nn; itot -= kper[k]; nn -= 1
    end
    scale2 = lskiph ? 1.0f0 : yr / Float32(ntyr)              # SCALE2 = YR/NTYR (every sub-model, LSKIPH false)
    dgsd = ctl.dg_sd
    cur_year = current_cycle_year(s)
    slo_a = sd[:site_lo]; shi_a = sd[:site_hi]
    aa = c.ht_dbh_aa; iabflg = c.ht_dbh_iabflg; lhtdrg = ctl.ht_drag_sp
    pccf_of(i) = (pt = Int(t.plot_id[i]); pv = isempty(pccf_in) ? dens.point_ccf : pccf_in;
                  (1 <= pt <= length(pv)) ? pv[pt] : 0f0)
    # ---- BANEXT/RDNEXT large-tree projection (:262-289), storage order DO 6 I=1,ITRN ----
    banext = fill(ba, max(nper, 1)); rdnext = fill(relden, max(nper, 1))
    temba = (lestb && atba > 0f0) ? atba : ba                 # :317-321 TEMBA=ATBA (≤0 ⇒ BA), TEMCCF=ATCCF, TEMAHT=ATAVH
    temccf = (lestb && atccf > 0f0) ? atccf : relden
    temaht = atavh_e
    if lestb                                                  # :392-404 label 8: interpolate from old/new density
        bayr = lskiph ? 0f0 : (ba - temba) / Float32(itot)    # ITOT = the remainder the DO 2 loop left (= KPER(NPER))
        ccfyr = lskiph ? 0f0 : (relden - temccf) / Float32(itot)
        nyr = 5
        @inbounds for j in 1:nper
            rdnext[j] = temccf + Float32(nyr) * ccfyr; banext[j] = temba + Float32(nyr) * bayr
            nyr += kper[j]
        end
    elseif nper > 1
        @inbounds for i in 1:n
            d1 = t.dbh[i]; d1 < 3.0f0 && continue
            is = Int(t.species[i]); pr = t.tpa[i]; pr <= 0f0 && continue
            d2 = d1 + t.diam_growth[i] / ci_bratio(sd, is, d1)
            b1 = CI_RG_BACON * d1 * d1; b2 = CI_RG_BACON * d2 * d2
            c1 = ci_tree_ccf(is, d1) * pr; c2 = ci_tree_ccf(is, d2) * pr     # CCFCAL: CCFT·P
            bi = (b2 - b1) / 10.0f0; ci = (c2 - c1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]; pn = pr * fpowi(0.985f0, k)
                rdnext[j] += Float32(k) * ci / pr * pn
                banext[j] += Float32(k) * bi * pn
            end
        end
    end
    pccfv = isempty(pccf_in) ? dens.point_ccf : pccf_in
    ptbav = isempty(ptba_in) ? dens.point_ba : ptba_in
    if lestb                                                  # :411-427 DO 13 (STORAGE order): crown of every new record
        @inbounds for i in (nstart + 1):n
            pt = Int(t.plot_id[i])
            cri0 = 0.89722f0 - 0.0000461f0 * ((1 <= pt <= length(pccfv)) ? pccfv[pt] : 0f0)
            ran = 0f0
            while true; ran = bachlo(s.rng, 0f0, 1f0); (-1f0 <= ran <= 1f0) && break; end
            cri0 = cri0 + 0.07985f0 * ran
            cri0 > 0.90f0 && (cri0 = 0.90f0); cri0 < 0.20f0 && (cri0 = 0.20f0)
            icr0 = trunc(Int32, cri0 * 100f0 + 0.5f0)
            t.crown_pct[i] = icr0
        end
    end
    wk3 = Float32[t.height[i] for i in 1:n]
    wk5 = Float32[t.dbh[i] for i in 1:n]
    ah = lestb ? temaht : avh; rr = lestb ? temccf : relden   # :434-440 R=TEMCCF, AH=TEMAHT under LESTB
    xd = ah * (rr / 100.0f0); xd > 300.0f0 && (xd = 300.0f0)   # :446 X=AH·(R/100)
    pctred = CI_RG_AB[1] + xd*(CI_RG_AB[2] + xd*(CI_RG_AB[3] + xd*(CI_RG_AB[4] + xd*(CI_RG_AB[5] + xd*CI_RG_AB[6]))))
    pctred > 1.0f0 && (pctred = 1.0f0); pctred < 0.01f0 && (pctred = 0.01f0)
    # IND1: the cycling call walks the lineage IND1; ESGENT calls SPESRT first (esgent.f:44) ⇒ species-major, record order
    ind1 = lestb ? sort(collect(1:n); by = i -> (Int(t.species[i]), i)) : species_major_order(s)
    # ---- DO 17 J=1,NPER ----
    ky = 0
    @inbounds for j in 1:nper
        ky += kper[j]
        baj = banext[j]; rdj = rdnext[j]
        ppccf = relden > 0f0 ? 1.0f0 + (rdj - relden) / relden : 0f0
        for i in ind1
            sp = Int(t.species[i]); (1 <= sp <= 19) || continue
            ttvar = _ci_ttvar(sp); utvar = _ci_ut_species(sp); civar = !(ttvar || utvar)
            ((civar || utvar) && j > 1) && continue                  # :528
            (ttvar && lskiph) && continue                             # :527
            d = t.dbh[i]; h = t.height[i]
            d >= CI_RG_XMAX[sp] && continue
            (lestb && i <= nstart) && continue                        # :591 IF(LESTB.AND.I.LT.ITRNIN)
            pr = t.tpa[i]; pr <= 0f0 && continue
            scale = ttvar ? Float32(kper[j]) / regyr : utvar ? Float32(ntyr) / yr : Float32(ntyr) / regyr
            xrhgro = active_multiplier(ctl, :regh, sp, cur_year)
            xrdgro = active_multiplier(ctl, :regd, sp, cur_year)
            con = fexp(c.htg_cor_small[sp])                            # RHCON(=1)·EXP(HCOR)
            h1 = wk3[i]; d1 = wk5[i]
            cri = Float32(t.crown_pct[i])
            tpccf = pccf_of(i) * ppccf; tpccf > 300f0 && (tpccf = 300f0); tpccf < 25f0 && (tpccf = 25f0)
            local htgrl::Float32
            if civar
                pt = Int(t.plot_id[i])
                ptbaa = (1 <= pt <= length(ptbav)) ? ptbav[pt] : ba
                pct = t.crown_ratio[i]
                ptbali = (lestb ? ba : ptbaa) * (1.0f0 - pct / 100.0f0)          # :619-623
                relht = atavh_e > 0f0 ? h / atavh_e : avh > 0f0 ? h / avh : 1.5f0  # :624-630 ATAVH, else AVH, else 1.5
                relht > 1.5f0 && (relht = 1.5f0)
                tbal = (1.0f0 - pct / 100.0f0) * ba
                iicr = ((Int(t.crown_pct[i]) - 1) ÷ 10) + 1; iicr > 9 && (iicr = 9)
                rcr = Float32(iicr)
                htgrl = CI_RG_CNST[sp] + CI_RG_HTBA[sp]*relht*ptbaa + CI_RG_PBAL[sp]*ptbali +
                        CI_RG_STBA[sp]*ba + CI_RG_RLHT[sp]*relht + CI_RG_CRSQ[sp]*rcr*rcr +
                        CI_RG_CR[sp]*rcr + CI_RG_BAL[sp]*tbal + CI_RG_HDM1[sp]*rhdm1 +
                        CI_RG_HDM2[sp]*rhdm2 + CI_RG_PTBA[sp]*ptbaa
            elseif ttvar
                htgrl = _ci_smhtgf!(s, i, cri, tpccf)
            elseif sp == 13                                           # aspen: FINDAG age, Sheppard curve
                slo = slo_a[sp]; shi = shi_a[sp]
                si = p.sp_site_index[sp]; si > shi && (si = shi); si <= slo && (si = slo + 0.5f0)
                rsimod = 0.5f0 * (1.0f0 + (si - slo) / (shi - slo))
                sitage = lestb ? t.birth_age[i] :                     # :660-672 LESTB ⇒ SITAGE=ABIRTH(I)
                         fpow(h * 2.54f0 * 12.0f0 / 26.9825f0, 1.0f0 / 1.1752f0)
                hite1 = 26.9825f0 * fpow(sitage, 1.1752f0); hite2 = 26.9825f0 * fpow(sitage + 10.0f0, 1.1752f0)
                htgrl = (hite2 - hite1) / (2.54f0 * 12.0f0) * rsimod * con * 0.75f0
            else                                                      # 14,15,17,19: POTHTG·PCTRED·VIGOR·CON
                sj = p.sp_site_index[sp]
                pothtg = ((sj / 5.0f0) * (sj * 1.5f0 - h) / (sj * 1.5f0)) * 0.83f0
                xx = cri / 100.0f0
                vigor = 150.0f0 * fpow(xx, 3.0f0) * fexp(-6.0f0 * xx) + 0.3f0; vigor > 1.0f0 && (vigor = 1.0f0)
                sp == 14 && (vigor = 1.0f0 - (1.0f0 - vigor) / 3.0f0)
                htgrl = pothtg * pctred * vigor * con
            end
            h2 = civar ? h1 + htgrl * scale * xrhgro * con :
                 ttvar ? h1 + htgrl * scale * xrhgro * con * wk4(i) : h1 + htgrl * scale * xrhgro
            wk3[i] = h2
            # UPDATE DENSITY FOR NEXT SUBCYCLE (:717-781)
            local d2::Float32
            upd = true
            if civar
                d2 = d; wk5[i] = d2
            elseif ttvar
                rd = pccf_of(i); hless4 = h2 - 4.5f0
                d2 = 0.000231f0*hless4*cri - 0.00005f0*hless4*rd + 0.001711f0*cri + 0.17023f0*hless4 + 0.3f0
                d2 *= xrdgro; d2 < CI_RG_DIAM[sp] && (d2 = CI_RG_DIAM[sp])
                wk5[i] = d2
                (j >= nper || d >= 3.0f0 || h2 <= 4.5f0) && (upd = false)
            elseif sp == 13
                rd = pccf_of(i)
                d2 = -0.41227f0 + 0.16944f0*h2 + 0.003191f0*cri - 0.00220f0*rd
            else
                d2 = d; wk5[i] = d2
                if (sp == 17 || sp == 19) && d >= _ci_bkpt(sp)
                    upd = false
                elseif h2 <= 4.5f0
                    d2 = d + 0.001f0 * h2; wk5[i] = d2
                end
            end
            if upd && j < nper
                c1 = ci_tree_ccf(sp, d1) * pr; c2 = ci_tree_ccf(sp, d2) * pr
                f = fpowi(0.985f0, ky)
                rdnext[j+1] += Float32(ky) * (c2 - c1) / 10.0f0 * f
                banext[j+1] += (CI_RG_BACON * (d2 * d2 - d1 * d1)) * pr * f
            end
        end
    end
    # ---- DO 30 ISPC / DO 25 I3 (+ L-loop over the tripled records) ----
    @inbounds for i in ind1
        sp = Int(t.species[i]); (1 <= sp <= 19) || continue
        ttvar = _ci_ttvar(sp); utvar = _ci_ut_species(sp); civar = !(ttvar || utvar)
        bkpt = _ci_bkpt(sp)
        scale = ttvar ? yr / fint : utvar ? Float32(ntyr) / yr : Float32(ntyr) / regyr
        xrdgro = active_multiplier(ctl, :regd, sp, cur_year)
        xmx = CI_RG_XMAX[sp]; xmn = CI_RG_XMIN[sp]
        dgmx = ttvar ? fint * CI_RG_DGMAX[sp] : CI_RG_DGMAX[sp] * scale
        d = t.dbh[i]
        cri = Float32(t.crown_pct[i])
        iicr = ((Int(t.crown_pct[i]) - 1) ÷ 10) + 1; iicr > 9 && (iicr = 9)
        rcr = Float32(iicr)
        d >= xmx && continue
        (lestb && i <= nstart) && continue                            # :864 IF(LESTB.AND.I.LT.ITRNIN) GO TO 25
        t.tpa[i] <= 0f0 && continue
        h = t.height[i]; hkw = wk3[i]
        dk5 = wk5[i]
        htgr1 = hkw - h
        (ttvar && d < CI_RG_DIAM[sp]) && (d = CI_RG_DIAM[sp])
        if ttvar && lskiph                                            # :887-905 LSKIPH TTVAR: DBH from SMDGF at H
            if h >= 4.5f0
                tp = pccf_of(i); hless4 = h - 4.5f0
                dkl = 0.000231f0*hless4*cri - 0.00005f0*hless4*tp + 0.001711f0*cri + 0.17023f0*hless4 + 0.3f0
                dkl *= xrdgro; dkl < d && (dkl = d)
                t.dbh[i] = dkl; dgl = dkl; dgl > dgmx && (dgl = dgmx); t.diam_growth[i] = dgl
            end
            continue
        end
        civar && htgr1 < 0.0f0 && (htgr1 = 0.0f0)
        dbh_i = t.dbh[i]                                              # DBH(I) (DGDRIV put it in the copy slots too)
        large_htg = lestb ? 0f0 : t.ht_growth[i]
        nrec = (!lestb && stash !== nothing && !ttvar) ? 3 : 1        # :1276 no tripling under LESTB
        for l in 0:(nrec - 1)
            htgr = htgr1
            if (civar || utvar) && lskiph                             # :883-885 LSKIPH: HTG(K)=0, straight to the DBH
                htgr = 0f0
            elseif !ttvar
                zz = _ci_zzran(s, dgsd)
                if civar
                    htgr = htgr1 + zz * 0.1f0 * scale
                elseif sp == 17 || sp == 19
                    htgr = (htgr1 + zz * 0.2f0 * (Float32(ntyr) / 10.0f0)) * wk4(i)
                else                                                  # 13,14,15
                    htgr = htgr1 + zz * 0.1f0 * (Float32(ntyr) / 10.0f0)
                end
            end
            local htg::Float32
            if (civar || utvar) && lskiph
                htg = 0f0
            else
                htgr < 0.1f0 && (htgr = 0.1f0)                        # Dixon 3/4/09
                xwt = (d <= xmn || lestb) ? 0.0f0 : (d - xmn) / (xmx - xmn)   # :971 XWT=0 under LESTB
                lh = l == 0 ? large_htg : (stash.htg_copy[i] ? (l == 1 ? stash.htgU[i] : stash.htgL[i]) : large_htg)
                htg = htgr * (1.0f0 - xwt) + xwt * lh
                capH = ctl.sp_size_cap[sp, 4]
                if h + htg > capH
                    htg = capH - h; htg < 0.1f0 && (htg = 0.1f0)
                end
            end
            dbhK = dbh_i; dgK = 0f0; below = d < bkpt
            if below
                hk = h + htg
                if !ttvar && hk < 4.5f0
                    dgK = 0.0f0; dbhK = d + 0.001f0 * hk             # :993-996 DBH(K)=D+0.001·HK
                else
                    dk = 0f0; dkk = 0f0
                    if civar
                        dk = fexp(CI_RG_DHCN[sp] + CI_RG_DHHT[sp] * flog(hk) + CI_RG_DHCR[sp] * flog(rcr))
                        dkk = h < 4.5f0 ? d : fexp(CI_RG_DHCN[sp] + CI_RG_DHHT[sp] * flog(h) + CI_RG_DHCR[sp] * flog(rcr))
                    elseif utvar
                        if sp == 14
                            sit = p.sp_site_index[sp]
                            dk = (hk - 4.5f0) * 10.0f0 / (sit - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
                            dkk = (h - 4.5f0) * 10.0f0 / (sit - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0)
                            h < 4.5f0 && (dkk = d)
                        else
                            bx = CI_BLK_HT2[sp]; ax = iabflg[sp] == 1 ? CI_BLK_HT1[sp] : aa[sp]
                            if sp == 15
                                dkk = 3.1020f0 + 0.0210f0 * h; dkk < 0.0f0 && (dkk = d)
                                dk = 3.1020f0 + 0.0210f0 * hk; dk < dkk && (dk = dkk + 0.01f0)
                                if !lhtdrg[sp] || iabflg[sp] == 1
                                    dk = _ci_ut_htdbh(hk)
                                    dkk = h <= 4.5f0 ? d : _ci_ut_htdbh(h)
                                end
                            else                                      # 13,17,19 Wykoff inverse
                                dk = (bx / (flog(hk - 4.5f0) - ax)) - 1.0f0; dk < 0.1f0 && (dk = 0.1f0)
                                dkk = h <= 4.5f0 ? d : (bx / (flog(h - 4.5f0) - ax)) - 1.0f0
                            end
                        end
                    else                                              # TTVAR
                        dk = dk5                                      # :870 DK=WK5(I) (read again by the LESTB DBH)
                        if hk >= 4.5f0
                            rd = pccf_of(i); hless4 = h - 4.5f0
                            dkk = 0.000231f0*hless4*cri - 0.00005f0*hless4*rd + 0.001711f0*cri + 0.17023f0*hless4 + 0.3f0
                            bark = ci_bratio(sd, sp, d)
                            dgK = (dk - dkk) * bark
                            dds = dgK * (2.0f0 * bark * d + dgK) * scale2
                            dgK = sqrt(fpow(d * bark, 2f0) + dds) - bark * d
                        else
                            dgK = 0.0f0
                        end
                    end
                    if lestb                                          # :1152-1166 LESTB DBH assignment
                        if sp == 15
                            dat45 = 3.1020f0 + 0.0210f0 * 4.5f0
                            dbhK = (dat45 > 0f0 && hk >= 4.5f0 && lhtdrg[sp] && iabflg[sp] == 0) ?
                                   dk - dat45 + CI_RG_DIAM[sp] : dk
                        else
                            dbhK = dk
                        end
                        (ttvar && dgK > dgmx) && (dgK = dgmx)
                        if utvar || civar
                            dbhK < CI_RG_DIAM[sp] && (dbhK = CI_RG_DIAM[sp])
                            dbhK = dbhK + 0.001f0 * hk
                            dgK = dbhK
                        end
                    else
                        # non-LESTB (:1178-1271)
                        if civar
                            h < 4.5f0 && (dkk = d)
                            dgK = (dk - dkk) * ci_bratio(sd, sp, d) * xrdgro     # (stale BARK; recomputed below)
                            dgK < 0f0 && (dgK = 0f0)
                        elseif utvar
                            bark = ci_bratio(sd, sp, d)
                            if sp == 15
                                h < 4.5f0 && (dkk = d)
                                if dk < 0f0 || dkk < 0f0
                                    dgK = htg * 0.2f0 * bark * xrdgro; dk = d + dgK
                                else
                                    dgK = (dk - dkk) * bark * xrdgro
                                end
                                (lhtdrg[sp] && iabflg[sp] == 0) && (dgK = 0.1f0 * htg * xrdgro)
                                dgK < 0f0 && (dgK = 0.1f0)
                                dgK > dgmx && (dgK = dgmx)
                            else                                          # 13,14,17,19
                                if dk < 0f0 || dkk < 0f0
                                    dgK = htg * 0.2f0 * bark * xrdgro; dk = d + dgK
                                else
                                    dgK = (dk - dkk) * bark * xrdgro
                                end
                                dgK > dgmx && (dgK = dgmx)
                            end
                        end
                        dgK < 0f0 && (dgK = 0f0)
                        barkK = ci_bratio(sd, sp, dbhK)                  # BRATIO(ISPC,DBH(K),HT(K))
                        if civar
                            dgK = (dk - dkk) * barkK * xrdgro; dgK < 0f0 && (dgK = 0f0)
                        elseif utvar
                            dgK > dgmx && (dgK = dgmx)
                        end
                        dds = dgK * (2.0f0 * barkK * d + dgK) * scale2
                        dgK = sqrt(fpow(d * barkK, 2f0) + dds) - barkK * d
                    end
                    (dbhK + dgK) < CI_RG_DIAM[sp] && (dgK = CI_RG_DIAM[sp] - dbhK)
                end
                dgK = dg_bound(nothing, nothing, sp, dbhK, dgK, ctl.sp_size_cap)   # DGBND
            end
            if l == 0
                t.ht_growth[i] = htg
                if below
                    t.dbh[i] = dbhK; t.diam_growth[i] = dgK
                end
            elseif l == 1
                stash.htgU[i] = htg; stash.is_small[i] = true
                below && (stash.dgU[i] = dgK; stash.dbhU[i] = dbhK)
            else
                stash.htgL[i] = htg
                below && (stash.dgL[i] = dgK; stash.dbhL[i] = dbhK)
            end
        end
        # TTVAR tripling (:1277-1290): both copies take the central's DBH/DG/HT/HTG/ICR
        if ttvar && stash !== nothing
            stash.htgU[i] = t.ht_growth[i]; stash.htgL[i] = t.ht_growth[i]; stash.is_small[i] = true
            stash.dgU[i] = t.diam_growth[i]; stash.dgL[i] = t.diam_growth[i]
            stash.dbhU[i] = t.dbh[i]; stash.dbhL[i] = t.dbh[i]
        end
    end
    return s
end

# ci/esgent.f: SPESRT, CALL REGENT(.TRUE.,ITRNIN) — the cycling REGENT with LESTB (_ci_regent!, which also draws the new
# records' crowns, DO 13) — then HTG·WK4 (WK4 = HTIMLT, 1 here) added to HT, the WK4<1 DBH rescale, and the HHTMAX cap
# (HT=HHTMAX, DBH=2.95). REGENT reads the gradd.f:192 DENSE (post-growth, pre-regen: relden_pre/ba_pre/avh_pre/pccf_pre/
# ptba_pre) and grincr.f's post-thin ATBA/ATCCF/ATAVH. The old port grew the birth cohort with one deterministic
# step: no ZZRAN, no DBH=D+0.001·HK below breast height (MEASURED FVSci_g16 5388215010690 PLANT: DF 2025 DBH 0.10418 /
# HtG 3.172 live vs 0.101 / 3.208 jl).
const CI_HHTMAX = Float32[23.0, 27.0, 21.0, 21.0, 22.0, 20.0, 24.0, 18.0, 18.0, 17.0,
                          27.0, 27.0, 16.0, 6.0, 6.0, 27.0, 16.0, 22.0, 16.0]   # ci/blkdat.f:123-125
function ci_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0,
                    atavh::Float32 = -1f0, atba::Float32 = -1f0, atccf::Float32 = -1f0,
                    relden_pre::Float32 = -1f0, ba_pre::Float32 = -1f0, avh_pre::Float32 = -1f0,
                    pccf_pre::Vector{Float32} = Float32[], ptba_pre::Vector{Float32} = Float32[])
    t = s.trees
    nstart >= t.n && return s
    @inbounds for i in (nstart + 1):t.n                       # estab.f:1249-1253: a new record enters with DG=HTG=0, PCT=0
        t.diam_growth[i] = 0f0; t.ht_growth[i] = 0f0; t.crown_ratio[i] = 0f0
    end
    _ci_regent!(s, nothing; fint, lestb = true, nstart, ba_in = ba_pre, relden_in = relden_pre, avh_in = avh_pre,
                pccf_in = pccf_pre, ptba_in = ptba_pre, atba, atccf, atavh)
    @inbounds for i in (nstart + 1):t.n                       # esgent.f:58-71
        sp = Int(t.species[i])
        t.height[i] = t.height[i] + t.ht_growth[i]
        if 1 <= sp <= length(CI_HHTMAX) && t.height[i] > CI_HHTMAX[sp]
            t.height[i] = CI_HHTMAX[sp]; t.dbh[i] = 2.95f0
        end
    end
    return s
end

# The DG(I) ci/dgdriv.f DO 220 (:797-822) leaves on record i: the measured increment (capped at the inside-bark
# DBH when IDG<2) when DG>0 and HT>4.5; 0 when HT<=4.5; otherwise the dub from the second DGF(WK3) call (:795,
# COR final) — needs the calibration's dub_wk2/dub_wk3 stash. Read by the LSTART REGCAL DO 49.
@inline ci_do220_dg(s::StandState, i::Int, dcur::Float32)::Float32 =
    do220_dg(s, i, dcur, ci_bratio(s.coef.species, Int(s.trees.species[i]), dcur))

"""
    ci_regent_hcor_init!(s, isct, ind1, saved_dbh)

ci/regent.f:1305-1500 — the LSTART small-tree HEIGHT calibration (REGENT(.FALSE.,1) from ci/cratet.f:707). jl had
no port, so every CI species ran HCOR_init = 0. Per species with >= NCALHT(5) sub-5" records carrying a measured
HTG: CORNEW = Σ(HTG·SCALE3·P)/Σ(EDH·P), HCOR = ln(CORNEW), trapped to [0.0821, 12.1825] (CORNEW 1 ⇒ HCOR 0).
EDH by sub-model:
- CIVAR (1-10,18): 2.764559 − 0.009643·BA + 0.025303·RCR² (RCR = the crown class, BA the backdated stand BA), ×RHCON.
  Inside the J loop FVS backdates H by HTG AGAIN (`IF(IHTG.LT.2) H=H-HTG(I)`) and drops the record if H<0.01.
- TTVAR (11,12,16): SMHTGF over the NPER subcycles (BETA on TPCCF around the persistent ZRAND, drawn on the main
  stream), EDH = (HK−H)·RHCON.
- UTVAR (13,14,15,17,19): the 10-yr potential (aspen: Sheppard) ·0.5, the last subcycle's value.
SCALE3 = REGYR/FINTH: the `SELECT CASE (ISPC)` that would give the UT species 10/FINTH reads ISPC after DO 45 has
run it to MAXSP+1, so every species takes CASE DEFAULT. NTYR = IFINTH; PPCCF = (RDJ−TEMCCF)/TEMCCF. Stand values are
the cratet.f:262 backdating DENSE's (crown_init_lstart_dead_inclusive! snapshot). t.dbh is the backdated WK3 here.
"""
function ci_regent_hcor_init!(s::StandState, isct::AbstractMatrix, ind1::AbstractVector,
                              saved_dbh::AbstractVector)
    p, t, c = s.plot, s.trees, s.calib
    t.n == 0 && return s
    s.control.growth_ifinth == 0 && return s           # regent.f:1312 IF(IFINTH.EQ.0) GOTO 100
    ctl = s.control; sd = s.coef.species
    rhcon = ones(Float32, 19)                           # ci REGCON: 1.0, or RCOR2 under READCORR
    if ctl.regh_cor2_on
        @inbounds for sp in 1:min(19, length(ctl.regh_cor2)); ctl.regh_cor2[sp] > 0f0 && (rhcon[sp] = ctl.regh_cor2[sp]); end
    end
    snap = c.cratet_relden > 0f0
    ba = snap ? c.cratet_ba : p.basal_area; relden = snap ? c.cratet_relden : p.relative_density
    avh = snap ? c.cratet_avh : p.avg_height
    pccfv = (snap && !isempty(c.cratet_pccf)) ? c.cratet_pccf : s.density.point_ccf
    temccf = relden                                     # TEMCCF = ATCCF (0 at LSTART) ⇒ RELDEN
    finth = ctl.growth_finth > 0f0 ? ctl.growth_finth : 5f0
    scale3 = CI_RG_REGYR / finth                        # CASE DEFAULT (stale ISPC = MAXSP+1)
    ntyr = Int(ctl.growth_ifinth); iyr = Int(CI_RG_REGYR)
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1); nper < 1 && (nper = 1)
    kper = zeros(Int, nper); itot = ntyr; nn = nper
    @inbounds for k in 1:nper
        if nn == 1; kper[k] = itot; break; end
        kper[k] = itot ÷ nn; itot -= kper[k]; nn -= 1
    end
    banext = fill(ba, nper); rdnext = fill(temccf, nper)
    if nper > 1                                         # DO 49 (regent.f:1342-1360): every record, no D gate
        @inbounds for i in 1:t.n
            d1 = t.dbh[i]; sp = Int(t.species[i]); pr = t.tpa[i]
            d2 = d1 + ci_do220_dg(s, i, saved_dbh[i]) / ci_bratio(sd, sp, d1)
            b1 = 0.005454154f0 * d1 * d1; b2 = 0.005454154f0 * d2 * d2
            c1 = ci_tree_ccf(sp, d1) * pr; c2 = ci_tree_ccf(sp, d2) * pr
            bi = (b2 - b1) / 10.0f0; ci = (c2 - c1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]; pn = pr * fpowi(0.985f0, k)
                rdnext[j] += Float32(k) * ci / pr * pn; banext[j] += Float32(k) * bi * pn
            end
        end
    end
    x = avh * (relden / 100f0); x > 300f0 && (x = 300f0)
    pctred = CI_RG_AB[1] + x*(CI_RG_AB[2] + x*(CI_RG_AB[3] + x*(CI_RG_AB[4] + x*(CI_RG_AB[5] + x*CI_RG_AB[6]))))
    pctred > 1.0f0 && (pctred = 1.0f0); pctred < 0.01f0 && (pctred = 0.01f0)
    ihtg = ctl.growth_ihtg
    slo_a = sd[:site_lo]; shi_a = sd[:site_hi]            # SITERANGE(2,ISPC,SLO,SHI)
    @inbounds for sp in 1:19
        ttvar = sp == 11 || sp == 12 || sp == 16
        utvar = _ci_ut_species(sp)
        civar = !(ttvar || utvar)
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        slo = slo_a[sp]; shi = shi_a[sp]
        snp = 0f0; snx = 0f0; sny = 0f0; nh = 0
        for k in i1:i2
            i = Int(ind1[k])
            hg = t.ht_growth[i]
            h = t.height[i]; ihtg < 2 && (h -= hg)
            (saved_dbh[i] >= 5f0 || h < 0.01f0) && continue
            hg < 0.001f0 && continue
            hk = h; cri = Float32(t.crown_pct[i]); edh = 0f0
            ipccf = Int(t.plot_id[i])
            pccf_i = (1 <= ipccf <= length(pccfv)) ? pccfv[ipccf] : 0f0
            dropped = false
            for j in 1:nper
                rdj = rdnext[j]
                ppccf = temccf <= 0f0 ? 0f0 : (rdj - temccf) / temccf
                tpccf = pccf_i * ppccf
                tpccf > 300f0 && (tpccf = 300f0); tpccf < 25f0 && (tpccf = 25f0)
                if civar
                    iicr = ((Int(t.crown_pct[i]) - 1) ÷ 10) + 1; iicr > 9 && (iicr = 9)
                    rcr = Float32(iicr)
                    ihtg < 2 && (h -= hg)                          # regent.f:1425 — backdated AGAIN, every J
                    (saved_dbh[i] >= 5f0 || h < 0.01f0) && (dropped = true; break)
                    edh = 2.764559f0 - 0.009643f0 * ba + 0.025303f0 * rcr * rcr
                elseif ttvar                                       # ci/smhtgf.f CASE(11,12,16)
                    if t.zrand[i] == -999f0
                        while true
                            z = bachlo(s.rng, 0f0, 1f0)
                            (z < -2f0 || z > 2f0) && continue
                            t.zrand[i] = z; break
                        end
                    end
                    if saved_dbh[i] <= 0f0
                        edh = 0f0                                  # D<=0 ⇒ 0, bypasses the 0.1 floor (GO TO 900)
                    else
                        beta1 = fexp(1.17527f0 - 0.42124f0 * flog(tpccf))
                        beta2 = fexp(-2.56002f0 - 0.58642f0 * flog(tpccf))
                        htg1 = beta1 + beta2 * cri
                        stddev = htg1 * (1.08720f0 - 0.00230f0 * cri)
                        edh = htg1 + t.zrand[i] * stddev
                        edh <= 0.1f0 && (edh = 0.1f0; t.zrand[i] = -999f0)
                    end
                    hk += edh
                else                                               # UTVAR
                    si = p.sp_site_index[sp]
                    si > shi && (si = shi); si <= slo && (si = slo + 0.5f0)
                    relsi = (si - slo) / (shi - slo)
                    rsimod = 0.5f0 * (1f0 + relsi)
                    sj = p.sp_site_index[sp]
                    pothtg = ((sj / 5f0) * (sj * 1.5f0 - h) / (sj * 1.5f0)) * 0.83f0
                    xv = cri / 100f0
                    vigor = 150f0 * fpow(xv, 3f0) * fexp(-6f0 * xv) + 0.3f0; vigor > 1f0 && (vigor = 1f0)
                    sp == 14 && (vigor = 1f0 - (1f0 - vigor) / 3f0)
                    if sp == 13
                        ag1 = fpow(h * 12f0 * 2.54f0 / 26.9825f0, 0.8509f0)
                        h2 = (26.9825f0 * fpow(ag1 + 10f0, 1.1752f0)) / (2.54f0 * 12f0)
                        edh = (h2 - h) * rsimod * rhcon[sp] * 0.75f0
                        edh < 0f0 && (edh = 0f0)
                    else
                        edh = pothtg * pctred * vigor * rhcon[sp]
                    end
                    edh *= 0.5f0
                end
            end
            dropped && continue
            civar && (edh *= rhcon[sp])
            ttvar && (edh = (hk - h) * rhcon[sp])
            pr = t.tpa[i]
            snp += pr; snx += edh * pr; sny += hg * scale3 * pr; nh += 1
        end
        nh < 5 && continue                                         # NCALHT
        snx /= snp; sny /= snp
        cornew = sny / snx
        cornew <= 0f0 && (cornew = 1f-4)
        c.htg_cor_init[sp] = (cornew < 0.0821f0 || cornew > 12.1825f0) ? 0f0 : flog(cornew)
    end
    return s
end
