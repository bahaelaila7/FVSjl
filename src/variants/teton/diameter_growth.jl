# =============================================================================
# diameter_growth.jl (teton) — TT large-tree DDS (tt/dgf.f). Chunk 3.
#
# Hooks (EM/KT pattern):
#   tt_dgcons!(s)          — per-species per-stand DGCON/DGDSQ/DGCCF/ATTEN + bark (tt/dgf.f ENTRY DGCONS).
#   dgf!(s, ::Teton)       — per-tree WK2 = DDS (tt/dgf.f DO 10 body). TT stores outside-bark ln(DDS).
#
# TT DGCONS (differs from EM — a SITE-INDEX-CLASS term, not a per-species habitat table):
#   DGCON(sp) = DGSIC(IDGSIM(ISI,sp),sp)·XSITE + DGFOR(ISPFOR,sp)
#             + (DGSASP·sinθ + DGCASP·cosθ + DGSLOP)·SLOPE + DGSLSQ·SLOPE² + DGEL·ELEV + DGEL2·ELEV²
#   XSITE = SITEAR(ISISP) (ln for BI/MC); ISI=ISMAP(ISISP); ISPFOR=MAPLOC(IFOR,sp); DGDSQ=DGDS(MAPDSQ(IFOR,sp),sp);
#   DGCCF=DGCCFC(IGCCFM(ISC,sp),sp) [ISC=clamp(INT(SITEAR)/10,20-60)/10−1]; ATTEN=IBSERV(ISIC,sp) [sp10/13/16→OBSERV];
#   sp10(PP) adds DGHAB(ICHBCL(ITYPE)+1). θ=ASPECT−0.7854 (sp10/13/16 special).
#   CONSPP = DGCON(sp) + COR(sp) + 0.01·DGCCF(sp)·RELDEN  (RELDEN = stand CCF, p.relative_density).
#
# DDS forms (ttt01 = ES/AS/AF/LP/WB → MAIN + ASPEN only; others error until their test stand exists):
#   MAIN (1:3,5,7:9,17): CONSPP + DGLD·lnD + DGBAL·BAL + CR·DGCR + CR²·DGCRSQ + DGDSQ·D²  (BAL=(1−PCT/100)·BA/100)
#   ASPEN (6,14): DGFASP(D,CR,BARK,SI,RMSQD,BA) → ASPDG;  DDS = ASPDG + ln(COR2) + COR   (CR = raw crown pct)
# =============================================================================

# tt/dgf.f DATA SIGMAR (DG residual SD, per species) — dg_resid_sd in species_coefficients.csv. Prior var for
# the DGSCOR COR calibration (shared calibrate_diameter_growth! reads c.sigma).
# tt/dgdriv.f:112 DATA PSIGSQ — per-species prior variance for the DGSCOR COR calibration.
const TT_PSIGSQ = Float32[0.0408, 0.0586, 0.1556, 0.07, 0.0970, 0.1433, 0.0636, 0.0970, 0.0970, 0.0636,
                          0.07, 0.07, 0.0898, 0.1433, 0.07, 0.0898, 0.0858, 0.07]

# tt/bratio.f BARK1/BARK2/IMAP DATA (18 species). PP (sp10) is IMAP=4 (power: BARK1·D^(BARK2−1)),
# which the linear (a+b·d)/d shared bark_ratio cannot express — so TT bark routes through tt_bratio,
# mirroring CR's cr_bratio dispatch. IMAP 1/2/3 reproduce the linear form for d≥1 (ttt01 stays exact).
const TT_BARK1 = Float32[0.969, 0.969, 0.867, 0.0, 0.956, 0.969, 0.969, 0.956, 0.937, 0.809427,
                         0.0, 0.0, 0.94782, 0.950, 0.892, 0.9, 0.969, 0.892]
const TT_BARK2 = Float32[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.016866,
                         0.0, 0.0, 0.0836, 0.0, -0.086, 0.0, 0.0, -0.086]
const TT_IMAP  = Int[2, 2, 2, 1, 2, 2, 2, 2, 2, 4, 1, 1, 3, 2, 3, 2, 2, 3]

# tt/bratio.f BRATIO(IS,D) — faithful all-IMAP bark ratio. IMAP=4 caps at 0.97 and RETURNS (GO TO 100),
# skipping the 0.80 floor; IMAP 1/2/3 apply the 0.99 cap + 0.80 floor.
@inline function tt_bratio(sp::Int, d::Float32)::Float32
    ieqn = TT_IMAP[sp]; b1 = TT_BARK1[sp]; b2 = TT_BARK2[sp]
    temd = d < 1f0 ? 1f0 : d
    local br::Float32
    if ieqn == 1
        br = (b1 == 0f0 && b2 == 0f0) ? (0.9002f0 - 0.3089f0 * (1f0 / (temd > 19f0 ? 19f0 : temd))) :
             1f0 / (b1 + b2 * d)
    elseif ieqn == 2
        br = b1
    elseif ieqn == 3
        br = b1 + b2 * (1f0 / temd)
    else                                        # IMAP=4 (PP): power model, cap 0.97, no 0.80 floor
        d <= 0f0 && return 0.97f0
        br = (b1 * d^b2) / d
        return br > 0.97f0 ? 0.97f0 : br
    end
    br > 0.99f0 && (br = 0.99f0)
    br < 0.80f0 && (br = 0.80f0)
    return br
end

@inline function _tt_bark_ab(sd, sp::Int)
    # tt/bratio.f IMAP: 1(zero)→0.9002−0.3089/D; 2→const BARK1; 3→BARK1+BARK2/D. (IMAP 1-nonzero/4 = PM/UJ/PP,
    # not in ttt01 — handled in the DIAGR chunk.) Return (bark_a, bark_b) for bark_ratio=(a+b·d)/d.
    b1 = sd[:bark1][sp]; b2 = sd[:bark2][sp]; imap = round(Int, sd[:bark_imap][sp])
    if imap == 2
        return (0f0, b1)                          # BRATIO = BARK1 (constant)
    elseif imap == 3
        return (b2, b1)                           # BRATIO = BARK1 + BARK2/D
    elseif imap == 1 && b1 == 0f0 && b2 == 0f0
        return (-0.3089f0, 0.9002f0)              # BRATIO = 0.9002 − 0.3089/D
    else
        return (0f0, b1 == 0f0 ? 0.95f0 : b1)     # IMAP 1-nonzero/4 (PM/UJ/PP) — placeholder (not in ttt01)
    end
end

# tt/dgf.f ENTRY DGCONS — per-species per-stand DG constants.
function tt_dgcons!(s::StandState)
    c = s.calib; p = s.plot; ctl = s.control; sd = s.coef.species
    isisp = Int(p.site_species); (isisp < 1 || isisp > 18) && (isisp = 3)
    ifor  = Int(p.forest_idx);   (ifor < 1 || ifor > 4) && (ifor = 1)
    itype = Int(p.habitat_input)
    elev = p.elevation; aspect = p.aspect; slope = p.slope
    xss = p.sp_site_index[isisp]                       # SITEAR(ISISP)
    lsi = floor(Int, xss / 10f0)
    isic = lsi < 2 ? 1 : (lsi >= 5 ? 5 : lsi)          # site class 1-5 (tt/dgf.f:663-667)
    isic2 = trunc(Int, xss); isic2 > 59 && (isic2 = 60); isic2 < 20 && (isic2 = 20)
    isc = isic2 ÷ 10 - 1                               # 1-5 (tt/dgf.f:675-679)
    (isc < 1 || isc > 5) && (isc = clamp(isc, 1, 5))
    isi = Int(TT_ISMAP[isisp])
    @inbounds for sp in 1:18
        ksi = Int(TT_IGCCFM[isc, sp]); (ksi < 1 || ksi > 5) && (ksi = 1)
        c.dg_ccf[sp] = TT_DGCCFC[ksi, sp]
        ispfor = Int(TT_MAPLOC[ifor, sp]); (ispfor < 1 || ispfor > 5) && (ispfor = 1)
        ispdsq = Int(TT_MAPDSQ[ifor, sp]); (ispdsq < 1 || ispdsq > 4) && (ispdsq = 1)
        temel = elev; xsite = xss; asptem = aspect
        if sp == 16
            temel > 30f0 && (temel = 30f0); xsite = log(xsite); asptem = aspect
        elseif sp == 13
            xsite = log(xsite); asptem = aspect
        elseif sp == 10
            asptem = aspect
        else
            asptem = aspect - 0.7854f0
        end
        idgs = Int(TT_IDGSIM[isi, sp]); (idgs < 1 || idgs > 5) && (idgs = 1)
        dgcon = TT_DGSIC[idgs, sp] * xsite + TT_DGFOR[ispfor, sp] +
                (TT_DGSASP[sp] * sin(asptem) + TT_DGCASP[sp] * cos(asptem) + TT_DGSLOP[sp]) * slope +
                TT_DGSLSQ[sp] * slope * slope + TT_DGEL[sp] * temel + TT_DGEL2[sp] * temel * temel
        c.dg_dsq[sp] = TT_DGDS[ispdsq, sp]
        if sp == 10
            maphab = itype > 0 ? Int(TT_ICHBCL[itype]) + 1 : 2
            (maphab < 1 || maphab > 7) && (maphab = 2)
            dgcon += TT_DGHAB[maphab]
            c.atten[sp] = TT_OBSERV[sp]
        elseif sp == 13 || sp == 16
            c.atten[sp] = TT_OBSERV[sp]
        else
            c.atten[sp] = Float32(TT_IBSERV[isic, sp])
        end
        (ctl.dg_cor2_on && ctl.dg_cor2[sp] > 0f0) && (dgcon += log(ctl.dg_cor2[sp]))
        c.dg_const[sp] = dgcon
        c.bark_a[sp], c.bark_b[sp] = _tt_bark_ab(sd, sp)
    end
    return s
end

# tt/badist.f BAU = BA strictly ABOVE each DBH class (growth pass: DBH, skip HT<4.5 seedlings). Tree-index
# order for Float32 bit-exactness. Used by dgf CASE(15,18) BAUTBA = BAU(int(D+1))/BA.
function _tt_badist_bau(t)
    bau = zeros(Float32, 41); totba = 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] <= 0f0 && continue
        t.height[i] < 4.5f0 && continue                  # seedlings excluded (SEEDS bucket)
        tdbh = t.dbh[i]; tdbh < 1f0 && (tdbh = 1f0)
        icls = trunc(Int, t.dbh[i] + 1f0); icls > 41 && (icls = 41); icls < 1 && (icls = 1)
        treeba = 0.0054542f0 * tdbh * tdbh * t.tpa[i]
        totba += treeba; bau[icls] += treeba
    end
    bau[1] = totba - bau[1]; bau[1] < 0f0 && (bau[1] = 0f0)
    @inbounds for j in 2:41
        bau[j] = bau[j-1] - bau[j]; bau[j] < 0f0 && (bau[j] = 0f0)
    end
    return bau
end

# tt/dgf.f DO 10 — per-tree WK2 = ln(DDS) (outside-bark). ttt01 exercises MAIN + ASPEN.
function dgf!(s::StandState, ::Teton)
    p, t, c = s.plot, s.trees, s.calib
    dens = s.density
    wk2 = view(s.scratch.wk, 2, :)
    relden = p.relative_density
    ba = p.basal_area
    ba100 = ba / 100f0
    logba = ba > 0f0 ? log(ba) : 0f0   # PP CONSPP term uses raw ALOG(BA) (dgf.f:520; TEMBA clamp is dead code)
    rmsqd = stand_qmd(s)
    bau = any(j -> (sp = Int(t.species[j]); sp == 15 || sp == 18), 1:t.n) ? _tt_badist_bau(t) : nothing
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        sp = Int(t.species[i])
        conspp = c.dg_const[sp] + c.dg_cor[sp] + 0.01f0 * c.dg_ccf[sp] * relden
        if sp <= 3 || sp == 5 || (7 <= sp <= 9) || sp == 17
            cr  = Float32(t.crown_pct[i]) * 0.01f0
            bal = (1f0 - t.crown_ratio[i] / 100f0) * ba100
            dds = conspp + TT_DGLD[sp] * log(d) + TT_DGBAL[sp] * bal +
                  cr * TT_DGCR[sp] + cr * cr * TT_DGCRSQ[sp] + c.dg_dsq[sp] * d * d
            dds < -9.21f0 && (dds = -9.21f0)
            wk2[i] = dds
        elseif sp == 6 || sp == 14
            cr = Float32(t.crown_pct[i])                       # aspen: raw crown pct (÷10 inside DGFASP)
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
            si = p.sp_site_index[sp]
            aspdg = _tt_dgfasp(d, cr, bark, si, rmsqd, ba)
            cor2 = (s.control.dg_cor2_on && s.control.dg_cor2[sp] > 0f0) ? s.control.dg_cor2[sp] : 1f0
            dds = aspdg + log(cor2) + c.dg_cor[sp]
            dds < -9.21f0 && (dds = -9.21f0)
            wk2[i] = dds
        elseif sp == 10
            # tt/dgf.f CASE(10) PP — Wykoff DGHAB form w/ point-BAL (PBAL). CONSPP gets −0.257322·ln(BA).
            cr  = Float32(t.crown_pct[i]) * 0.01f0
            ipccf = Int(t.plot_id[i])
            pbal = (1f0 - t.crown_ratio[i] / 100f0) * dens.point_ba[ipccf]
            csp = conspp - 0.257322f0 * logba
            dds = csp + TT_DGLD[sp] * log(d) + cr * (TT_DGCR[sp] + cr * TT_DGCRSQ[sp]) +
                  c.dg_dsq[sp] * d * d + TT_DGDBAL[sp] * pbal / log(d + 1f0)
            dds < -9.21f0 && (dds = -9.21f0)
            wk2[i] = dds
        elseif sp == 4 || sp == 11 || sp == 12
            # tt/dgf.f CASE(4,11,12) PM/UJ/RM DIAGR form. NOTE: REGENT overrides DG for ALL sizes (XMAX=99),
            # so this dgf DDS is a faithful placeholder (the source-variant eqn) that regent replaces.
            dpp = d < 1f0 ? 1f0 : d
            batem = ba < 1f0 ? 1f0 : ba
            si = p.sp_site_index[sp]
            bark = tt_bratio(sp, d)
            df = 0.25897f0 + 1.03129f0 * dpp - 0.0002025464f0 * batem + 0.00177f0 * si
            (df - dpp) > 1f0 && (df = dpp + 1f0)
            df < dpp && (df = dpp)
            diagr = (df - dpp) * bark
            if diagr <= 0f0
                wk2[i] = -9.21f0
            else
                dds = log(diagr * (2f0 * dpp * bark + diagr)) + conspp
                dds < -9.21f0 && (dds = -9.21f0)
                wk2[i] = dds
            end
        elseif sp == 13 || sp == 16
            # tt/dgf.f CASE(13,16) BI/MC — Wykoff form (surrogate from CR/SO) w/ extra DGPCCF·PCCF + DGBA·BA.
            ald = log(d)
            cr  = Float32(t.crown_pct[i]) * 0.01f0
            bal = (1f0 - t.crown_ratio[i] / 100f0) * ba
            ipccf = Int(t.plot_id[i])
            pccf = (1 <= ipccf <= length(dens.point_ccf)) ? dens.point_ccf[ipccf] : 0f0
            dds = conspp + TT_DGLD[sp] * ald + TT_DGBAL[sp] * bal + cr * (TT_DGCR[sp] + cr * TT_DGCRSQ[sp]) +
                  c.dg_dsq[sp] * d * d + TT_DGDBAL[sp] * bal / log(d + 1f0) +
                  TT_DGPCCF[sp] * pccf + TT_DGBA[sp] * ba
            dds < -9.21f0 && (dds = -9.21f0)
            wk2[i] = dds
        elseif sp == 15 || sp == 18
            # tt/dgf.f CASE(15,18) NC/OH — CR-surrogate DIAGR form. Uses COR+DGCON directly (NOT conspp: no
            # CCF·relden term). BAUTBA = BAU(int(D+1))/BA (BA-above-class). ISTAGF=0 ⇒ no DSTAG.
            icls = trunc(Int, d + 1f0); icls > 41 && (icls = 41); icls < 1 && (icls = 1)
            bautba = (bau !== nothing && ba > 0f0) ? bau[icls] / ba : 0f0
            si = p.sp_site_index[sp]
            dpp = d < 1f0 ? 1f0 : d
            batem = ba < 5f0 ? 5f0 : ba
            bark = tt_bratio(sp, d)
            df = (1.55986f0 + 1.01825f0 * dpp - 0.29342f0 * log(batem) +
                  0.00672f0 * si - 0.00073f0 * bautba) * 1.05f0
            df < dpp && (df = dpp)
            diagr = (df - dpp) * bark
            if diagr <= 0f0
                wk2[i] = -9.21f0
            else
                dds = log(diagr * (2f0 * dpp * bark + diagr)) + c.dg_cor[sp] + c.dg_const[sp]
                dds < -9.21f0 && (dds = -9.21f0)
                wk2[i] = dds
            end
        else
            error("TT dgf! DDS form for sp $sp not yet ported — not in ttt01")
        end
    end
    return s
end

# bin/FVStt_buildDir/dgfasp.f (shared ut/) — large aspen DG → ln(DDS).
@inline function _tt_dgfasp(d::Float32, cr::Float32, bark::Float32, si::Float32, rmsqd::Float32, ba::Float32)::Float32
    rel = rmsqd > 0f0 ? d / rmsqd : 0f0
    aspcr = cr / 10f0
    pot = (0.4755f0 - 3.8336f-6 * d^4.1488f0) + (4.510f-2 * aspcr * d^0.67266f0)
    pot <= 0f0 && (pot = 0.01f0)
    fofr = 1.07528f0 * (1f0 - exp(-1.89022f0 * rel))
    gofad = 2.1963f-1 * (rmsqd + 1f0)^0.73355f0
    baact = ba >= 310f0 ? 305f0 : ba
    valmod = 1f0 - exp(-fofr * gofad * ((310f0 - baact) / 310f0)^0.5f0)
    predgr = pot * valmod * (0.48630f0 + 0.01258f0 * si)
    return log(2f0 * d * bark * predgr + predgr * predgr)
end
