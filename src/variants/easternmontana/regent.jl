# =============================================================================
# regent.jl (easternmontana) — EM small-tree growth, em/regent.f.
#
# small_tree_growth! is the non-ESTAB REGENT in the Fortran's own shape: one subcycle loop (DO 17 J / DO 16 ISPC /
# DO 15 I3 over IND1) that grows every sub-model — EMVAR (SMHTGF/SMDGF, persistent ZRAND), NIVAR (LL, log-space
# HTGRL + AX·(H−4.5)^BX DBH), TTVAR (LM, persistent ZRAND + DLESS3), CRVAR/UTVAR (POTHTG·PCTRED·VIGOR / Sheppard
# aspen, one pass) — with the running RDNEXT/BANEXT(J+1) density feedback; then one assembly loop (DO 30 ISPC /
# DO 25 I3) that draws ZZRAN, blends toward the large-tree HTG over [XMIN,XMAX], sets the D<3 DBH increment and loops
# back at 918 for each tripled copy (NIVAR/CRVAR/UTVAR) or duplicates the central record (EMVAR/TTVAR).
# em_regent_aspen_calib! is the LSTART calibration section (label 40); em_esgent! the ESTAB birth-cycle entry.
# =============================================================================

const EM_RG_XMAX = Float32[3,3,3,3,10,99,3,3,3,3,2,4,2,2,2,2,4,3,2]      # em/regent.f DATA XMAX (OH=2.0)
const EM_RG_XMIN = Float32[1.5,1.5,1.5,1.5,2,90,1.5,1.5,1.5,1.5,0.5,2,0.5,0.5,0.5,0.5,2,1.5,0.5]  # DATA XMIN (OH=0.5)
const EM_RG_DIAM = Float32[0.4,0.3,0.3,0.3,0.3,0.3,0.4,0.3,0.3,0.5,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2]
const EM_RG_RHGL = Float32[-0.2785, -0.0480, 0.0]                 # geo-location (IGL 1..3)
const EM_RG_RHHAB = Float32[-0.2146, -0.0941, -0.4916, -0.3582, 0.0]
const EM_RG_MAPHAB = Int[4,4,4,4,4,4,4,4,4,4,4,4,3,2,5,5,5,5,1,4,4,4,4,4,4,4,1,3,3,4]   # ITYPE→RHHAB idx
const _EM_RG_RSAB0 = -0.10987f0; const _EM_RG_RSAB1 = 0.22157f0; const _EM_RG_RSAB2 = -0.12432f0
const _EM_RG_BH = 0.3740f0; const _EM_RG_BCCF = -0.00391f0; const _EM_RG_BBAL = -0.22957f0
const _EM_RG_AX = 0.0658f0; const _EM_RG_BX = 1.3817f0
const _EM_RG_HSIGMA = 0.59f0; const _EM_RG_REGYR = 5.0f0
# Curtis-Arney/Wykoff height-diameter coefficients (em/blkdat.f HT1/HT2, IABFLG=1 for all sp ⇒ AX=HT1),
# and the CRVAR/UTVAR small-tree DGMAX cap (em/regent.f DATA DGMAX). Used by the height→diameter DK model.
const EM_RG_HT1 = Float32[4.1539,4.1539,4.4161,4.1920,4.76537,3.2,4.5356,4.7537,4.5788,4.414,4.4421,4.4421,4.4421,4.4421,4.4421,4.4421,4.4421,4.1539,4.4421]
const EM_RG_HT2 = Float32[-4.212,-4.212,-6.962,-5.1651,-7.61062,-5.0,-5.692,-8.356,-7.138,-8.907,-6.5405,-6.5405,-6.5405,-6.5405,-6.5405,-6.5405,-6.5405,-4.212,-6.5405]
const EM_RG_DGMAX = Float32[0,0,0,2,99,2,0,0,0,0,2.5,2.5,2.5,2.5,2.5,2.5,2.5,0,2.5]


# EMVAR small-tree models (em/smhtgf.f height + em/smdgf.f diameter; coeffs em/blkdat.f:253-270 + smdgf.f DATA).
# EMVAR = _em_orig_species {1,2,3,7,8,9,10,18}. FVS uses these, NOT the NI exp-form (which is sp5/NIVAR only).
const EM_B0ACCF = Float32[1.17527,1.17527,-4.35709,0,0,0,-0.90086,-0.55052,-4.35709,0.405,0,0,0,0,0,0,0,1.17527,0]
const EM_B1ACCF = Float32[-0.42124,-0.42124,0.67307,0,0,0,0.16996,-0.02858,0.67307,0,0,0,0,0,0,0,0,-0.42124,0]
const EM_B0BCCF = Float32[-2.56002,-2.56002,-2.49682,0,0,0,-1.50963,-2.26007,-2.49682,-1.50963,0,0,0,0,0,0,0,-2.56002,0]
const EM_B1BCCF = Float32[-0.58642,-0.58642,-0.51938,0,0,0,-0.61825,-0.67115,-0.51938,-0.61825,0,0,0,0,0,0,0,-0.58642,0]
const EM_B0ASTD = Float32[1.08720,1.08720,1.13785,0,0,0,1.00749,1.09730,1.13785,0.57707,0,0,0,0,0,0,0,1.08720,0]
const EM_B1BSTD = Float32[-0.00230,-0.00230,-0.00185,0,0,0,-0.00435,-0.00130,-0.00185,0.00055,0,0,0,0,0,0,0,-0.00230,0]
const EM_SDHTCR = Float32[0.000231,0.000231,-0.28654,0,0,0,-0.41227,0.04125,-0.15906,0.000335,0,0,0,0,0,0,0,0.000231,0]
const EM_SDHPCF = Float32[-0.00005,-0.00005,0.13469,0,0,0,0.16944,0.17486,0.15323,-0.00020,0,0,0,0,0,0,0,-0.00005,0]
const EM_SDCR   = Float32[0.001711,0.001711,0.002736,0.001711,0.001711,0.003191,0.003191,-0.002371,0.0,0.002621,0,0,0,0,0,0,0,0.001711,0]
const EM_SDHL4  = Float32[0.17023,0.17023,0.00036,0.17023,0.17023,-0.00220,-0.00220,-0.00070,0.0,0.15622,0,0,0,0,0,0,0,0.17023,0]

# SMHTGF (em/smhtgf.f): small-tree HEIGHT increment (real ft). ZRAND drawn once per tree by the caller.
@inline function _em_smhtgf(sp::Int, cr::Float32, tpccf::Float32, zrand::Float32)::Float32
    beta1 = exp(EM_B0ACCF[sp] + EM_B1ACCF[sp]*log(tpccf))
    beta2 = exp(EM_B0BCCF[sp] + EM_B1BCCF[sp]*log(tpccf))
    htg1 = beta1 + beta2*cr
    stddev = htg1*(EM_B0ASTD[sp] + EM_B1BSTD[sp]*cr)
    htgrth = htg1 + zrand*stddev
    return htgrth > 0.1f0 ? htgrth : 0.1f0
end
# SMDGF (em/smdgf.f): small-tree DBH from height. sp{3,7,8,9} linear; else{1,2,10,18} HLESS4-form. RD=TPCCF.
@inline function _em_smdgf(sp::Int, h::Float32, cr::Float32, rd::Float32)::Float32
    if sp == 3 || sp == 7 || sp == 8 || sp == 9
        return EM_SDHTCR[sp] + EM_SDHPCF[sp]*h + EM_SDCR[sp]*cr + EM_SDHL4[sp]*rd
    else
        hl = h - 4.5f0
        return EM_SDHTCR[sp]*hl*cr + EM_SDHPCF[sp]*hl*rd + EM_SDCR[sp]*cr + EM_SDHL4[sp]*hl + 0.3f0
    end
end

# em_esgent! (the birth-cycle ESTAB entry) still caps its EMVAR pass at XMAX (LL's is handled in REGENT proper).
@inline _em_rg_cap(sp::Int) = sp == 5 ? 3.0f0 : EM_RG_XMAX[sp]
# em/regent.f:322-333 PCTRED from X = AH·(RELDEN/100), used by the calibration pass (em_regent_aspen_calib!).
@inline function _em_rg_pctred(ah::Float32, relden::Float32)::Float32
    x = ah * (relden / 100f0); x > 300f0 && (x = 300f0)
    pr = 1.11436f0 + x*(-0.011493f0 + x*(0.43012f-4 + x*(-0.72221f-7 + x*(0.5607f-10 - x*0.1641f-13))))
    pr > 1.0f0 && (pr = 1.0f0); pr < 0.01f0 && (pr = 0.01f0)
    return pr
end
@inline _em_rg_crvar(sp::Int) = sp == 11 || (13 <= sp <= 16) || sp == 19

# em/regent.f:1428-1488 REGCON entry. Only the NIVAR species (LL, sp5) carries the NI habitat/site constant
# RHCON = REGCH + 1.0667 + RHHAB[MAPHAB[ITYPE]] (+ ln RCOR2 under READCORR); every other sub-model (EMVAR, TTVAR,
# CR/UT) has RHCON = 1.0 (or RCOR2). jl gave every species the NI constant, which the LSTART REGCAL divides into
# CORNEW: on the S248112 CW fixture RHCON 0.4704 for 1.0 ⇒ CORNEW 2.50 vs live 1.004 ⇒ CW height growth ~2×.
function em_regcons!(s::StandState)
    p = s.plot; ctl = s.control
    itype = Int(p.habitat_input); (itype < 1 || itype > 30) && (itype = 1)
    igl = Int(p.geo_location); (igl < 1 || igl > 3) && (igl = 1)
    irhhab = EM_RG_MAPHAB[itype]; (irhhab < 1 || irhhab > 5) && (irhhab = 5)
    regch = EM_RG_RHGL[igl] + (_EM_RG_RSAB0 + _EM_RG_RSAB1 * cos(p.aspect) + _EM_RG_RSAB2 * sin(p.aspect)) * p.slope
    rhcon = ones(Float32, 19)
    @inbounds for sp in 1:19
        rc2 = (ctl.regh_cor2_on && sp <= length(ctl.regh_cor2)) ? ctl.regh_cor2[sp] : 0f0
        if sp == 5
            rhcon[sp] = regch + 1.0667f0 + EM_RG_RHHAB[irhhab]
            rc2 > 0f0 && (rhcon[sp] += log(rc2))
        else
            rc2 > 0f0 && (rhcon[sp] = rc2)
        end
    end
    return rhcon
end

"""
    em_regent_aspen_calib!(s)

EM aspen/PB (sp 12/17, UTVAR Sheppard model) two-part setup, run ONCE per stand from the EM branch of
`setup_growth!` AFTER `calibrate_diameter_growth!` (mirrors the CR/TT/UT `_xx_dub_ages!` hooks). Neither
part was ported, so the tiny-aspen small-tree HEIGHT over-grew ~2× (stand 373796912489998: aspen BA cyc1
oracle 68 → jl 123; the inflated `hk=h+htg` feeds the aspen inverse-Wykoff DK diameter).

(A) ABIRTH dub (em/cratet.f:517 → findag.f:91 → pothtg.f:158-161). For aspen/PB with un-set birth_age
    (ABIRTH<=0) and PROB>0, dub the effective tree age from the CURRENT height by inverting the Sheppard
    height–age curve: EFAGE = (H·2.54·12/26.9825)^(1/1.1752). The growth path (below, lines 275-279) reads
    `t.birth_age` for HITE1/HITE2; without the dub it clamped ABIRTH to 1.0 (age 1) ⇒ huge HITE2−HITE1
    increment. `age_known` is set so `setup_growth!`'s per-cycle GRADD advance (gradd.f:205, ABIRTH+=FINT)
    keeps later cycles aged. A DB AGE column (birth_age>0) skips the dub, matching FVS.

(B) REGENT small-tree HEIGHT self-calibration CORNEW (em/regent.f:1298-1365, the LSTART `CALL REGENT(.T.,..)`
    UTVAR-aspen branch). For each sub-5"-DBH aspen/PB record carrying a measured height increment the predicted
    5-yr increment is EDH = (H2−H)·RSIMOD·RHCON(=1)·0.75·0.5, where H is the backdated start height, AG1=EFAGE(H),
    AG2=AG1+10, H2 = the Sheppard height at AG2 (the aspen closed form is J/subcycle-independent, so the DO-55
    subcycle loop leaves the single value). CORNEW = Σ(HTG·SCALE3·P)/Σ(EDH·P), SCALE3 = REGYR(5)/FINTH (the
    HTG_MEASURE period, default 5; =10 for the target stand ⇒ SCALE3=0.5, cancelling the EDH ·0.5). HCOR = ln(CORNEW),
    trapped to CORNEW∈[0.0821,12.1825]. The RAW HCOR seeds `htg_cor_init`; `diameter_growth!`'s shared per-cycle
    attenuation (`htg_cor_small = dg_cor_goal + cormlt_h·(htg_cor_init−dg_cor_goal)`, dgdriv.f WCI/CORMLT) then
    produces the applied CON = exp(htg_cor_small) the growth path multiplies in (was CON=1 ⇒ ~2× HEIGHT over-grow).
    Verified vs FVSem_g16 (stand 373796912489998 sp12): CORNEW=0.27368, HCOR=−1.29579, cyc1 CON=0.3746.

Runs at the setup call site where `calibrate_diameter_growth!` has already restored the read-in dbh
(southern/diameter_growth.jl:801) and transformed `t.ht_growth` to the measured INCREMENT (lines 247-252),
so `t.dbh`/`t.height` are read-in and `t.ht_growth` is the increment — read directly.
"""
function em_regent_aspen_calib!(s::StandState)
    p, t, c = s.plot, s.trees, s.calib
    n = t.n; n == 0 && return s
    dens = s.density
    slo = s.coef.species[:site_lo]; shi = s.coef.species[:site_hi]
    ihtg = s.control.growth_ihtg
    finth = s.control.growth_finth > 0f0 ? s.control.growth_finth : 5f0   # em/grinit.f:192 FINTH default 5
    scale3 = _EM_RG_REGYR / finth                                          # em/regent.f:1095 SCALE3 = REGYR/FINTH
    # (A) ABIRTH dub (em/cratet.f:517 FINDAG → em/pothtg.f), every record with ABIRTH<=0, on the current HT. POTHTG
    # sets ABIRTH=EFAGE for WB/WL/LP/OS (label 100), DF (200), ES/AF (300), PP (400) and the aspen/PB pair (500,
    # PROB>0 only); the other species keep 0. Climate-FVS reads it (clgmult BIRTHYR = THISYR − ABIRTH), and REGENT's
    # aspen Sheppard model reads it for 12/17. jl dubbed only the aspen pair.
    em_findag_abirth!(s)
    # (B) em/regent.f:1195-1360 REGCAL — the small-tree HEIGHT self-calibration, for EVERY species (jl ran it
    # for the aspen pair only). Per species: EDH is the model's predicted increment over the sub-periods and
    # CORNEW = mean(HTG·SCALE3) / mean(EDH) ⇒ HCOR = ln(CORNEW), gated at NCALHT=5 records and trapped to
    # ±2.5σ in ln(C). The EMVAR branch calls SMHTGF, which DRAWS the persistent per-tree ZRAND — those draws
    # land HERE, before DGDRIV, which is why jl's whole main RNG stream sat ahead of live on EM stands.
    s.control.growth_ifinth == 0 && return s           # regent.f:1093 IF(IFINTH.EQ.0) GOTO 100 (no HCOR)
    rhcon = em_regcons!(s)
    # REGCAL runs inside CRATET (em/cratet.f:553) on the state the cratet.f:182 backdating DENSE left — nothing
    # recomputes density in between — so BA/RELDEN/AVH/RELDM1/PCCF/PCT are that DENSE's (snapshotted by
    # crown_init_lstart_dead_inclusive!), NOT the live-only current values jl held here. S248112 CW fixture:
    # live AVH 61.60 / RELDEN 103.66 ⇒ PCTRED 0.5379; the current 63.44 / 119.87 gave 0.459.
    snap = c.cratet_relden > 0f0
    ba = snap ? c.cratet_ba : p.basal_area; relden = snap ? c.cratet_relden : p.relative_density
    avh = snap ? c.cratet_avh : p.avg_height
    reldm1 = snap ? c.cratet_reldm1 : p.relative_density_prev; reldm1 <= 0f0 && (reldm1 = relden)
    pccfv = (snap && !isempty(c.cratet_pccf)) ? c.cratet_pccf : dens.point_ccf
    pctv = (snap && length(c.cratet_pct) >= n) ? c.cratet_pct : t.crown_ratio
    temccf = relden                                     # regent.f:206-207 TEMCCF = ATCCF, else RELDEN
    ntyr = Int(s.control.growth_ifinth); iyr = Int(_EM_RG_REGYR)   # regent.f:185 LSTART ⇒ NTYR = IFINTH (5 unless DB HTG_MEASURE)
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1); nper < 1 && (nper = 1)
    kper = zeros(Int, nper); itot = ntyr; nn = nper
    @inbounds for k in 1:nper
        if nn == 1; kper[k] = itot; break; end
        kper[k] = itot ÷ nn; itot -= kper[k]; nn -= 1
    end
    banext = fill(ba, nper); rdnext = fill(relden, nper)
    # regent.f:1110-1128 DO 49 — with >1 subcycle (IFINTH>5: FIA HTG_MEASURE 6/7/10/15 are common) the start-of-
    # subcycle BA/CCF are extrapolated from the backdated diameter WK3 and the DGDRIV DO-220 increment, for EVERY
    # record (no D≥3 gate, unlike the growth-side copy) and in storage order.
    if nper > 1 && length(c.dub_wk3) == n
        @inbounds for i in 1:n
            d1 = c.dub_wk3[i]; sp = Int(t.species[i]); pr = t.tpa[i]
            d2 = d1 + em_do220_dg(s, i) / em_bratio(sp, d1)
            b1 = _EM_RG_BACON * d1 * d1; b2 = _EM_RG_BACON * d2 * d2
            c1 = em_tree_ccf(sp, d1) * pr; c2 = em_tree_ccf(sp, d2) * pr
            bi = (b2 - b1) / 10.0f0; ci = (c2 - c1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]; pn = pr * fpowi(0.985f0, k)
                rdnext[j] += Float32(k) * ci / pr * pn
                banext[j] += Float32(k) * bi * pn
            end
        end
    end
    pctred = _em_rg_pctred(avh, relden)
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    species_sort!(s)
    @inbounds for sp in 1:MAXSP
        isct[sp, 1] == 0 && continue
        nivar = sp == 5; ttvar = sp == 4
        crvar = _em_rg_crvar(sp); utvar = (sp == 6 || sp == 12 || sp == 17)
        emvar = !(nivar || ttvar || crvar || utvar)
        si = p.sp_site_index[sp]
        si > shi[sp] && (si = shi[sp]); si <= slo[sp] && (si = slo[sp] + 0.5f0)
        relsi = (shi[sp] > slo[sp]) ? (si - slo[sp]) / (shi[sp] - slo[sp]) : 0f0
        rsimod = 0.5f0 * (1f0 + relsi)
        snp = 0f0; snx = 0f0; sny = 0f0; nrec = 0
        for k in Int(isct[sp, 1]):Int(isct[sp, 2])
            i = Int(ind1[k])
            t.tpa[i] <= 0f0 && continue
            htg = t.ht_growth[i]
            h = t.height[i]; ihtg < 2 && (h -= htg)
            (t.dbh[i] >= 5f0 || h < 0.01f0 || htg < 0.001f0) && continue
            hk = h; d = t.dbh[i]; cr = Float32(t.crown_pct[i])
            pt = Int(t.plot_id[i])
            pccf = (1 <= pt <= length(pccfv)) ? pccfv[pt] : 0f0
            edh = 0f0
            for jj in 1:nper
                baj = banext[jj]; rdj = rdnext[jj]
                bal = baj * (100f0 - pctv[i]) * 0.0001f0
                ppccf = emvar ? (temccf <= 0f0 ? 0f0 : (rdj - temccf) / temccf) : (rdj - reldm1) / reldm1
                tpccf = pccf * ppccf
                tpccf > 300f0 && (tpccf = 300f0); tpccf < 25f0 && (tpccf = 25f0)
                if emvar
                    if t.zrand[i] == -999f0                      # smhtgf.f: persistent per-tree ZRAND
                        while true
                            z = bachlo(s.rng, 0f0, 1f0)
                            (z < -2f0 || z > 2f0) && continue
                            t.zrand[i] = z; break
                        end
                    end
                    edh = _em_smhtgf(sp, cr, tpccf, t.zrand[i]); hk += edh
                    edh <= 0.1f0 && (t.zrand[i] = -999f0)        # smhtgf.f: a floored increment resets ZRAND
                elseif nivar
                    edh = exp(rhcon[sp] + _EM_RG_BH * log(hk) + _EM_RG_BCCF * rdj + _EM_RG_BBAL * bal); hk += edh
                elseif ttvar
                    if t.zrand[i] == -999f0
                        while true
                            z = bachlo(s.rng, 0f0, 1f0)
                            (z < -2f0 || z > 2f0) && continue
                            t.zrand[i] = z; break
                        end
                    end
                    if d <= 0f0
                        edh = 0f0
                    else
                        beta1 = exp(1.17527f0 - 0.42124f0 * log(tpccf))
                        beta2 = exp(-2.56002f0 - 0.58642f0 * log(tpccf))
                        htg1 = beta1 + beta2 * cr
                        stddev = htg1 * (1.08720f0 - 0.00230f0 * cr)
                        edh = htg1 + t.zrand[i] * stddev
                        edh <= 0.1f0 && (edh = 0.1f0; t.zrand[i] = -999f0)
                    end
                    hk += edh
                else                                              # CRVAR / UTVAR — EDH is the LAST sub-period's
                    x = cr / 100f0
                    vigor = 150f0 * x^3 * exp(-6f0 * x) + 0.3f0; vigor > 1f0 && (vigor = 1f0)
                    sp == 6 && (vigor = 1f0 - (1f0 - vigor) / 3f0)
                    if sp == 12 || sp == 17
                        ag1 = (h * 12f0 * 2.54f0 / 26.9825f0)^0.8509f0
                        h2 = (26.9825f0 * (ag1 + 10f0)^1.1752f0) / (2.54f0 * 12f0)
                        edh = (h2 - h) * rsimod * rhcon[sp] * 0.75f0
                    else
                        pothtg = crvar ? si / (15f0 - 4f0 * relsi) :
                                 (sp == 15 || sp == 16) ? (si / 10f0) * (si * 1.5f0 - h) / (si * 1.5f0) :
                                 sp == 6 ? ((si / 5f0) * (si * 1.5f0 - h) / (si * 1.5f0)) * 0.83f0 : 0f0
                        edh = pothtg * pctred * vigor * rhcon[sp]
                    end
                    edh *= 0.5f0
                end
            end
            edh_t = emvar ? (hk - h) * rhcon[sp] : nivar ? (hk - h) : ttvar ? (hk - h) * rhcon[sp] : edh
            pr = t.tpa[i]
            snp += pr; snx += edh_t * pr; sny += htg * scale3 * pr; nrec += 1
        end
        nrec < 5 && continue                                       # NCALHT (regent.f:1352)
        snx /= snp; sny /= snp
        cornew = snx != 0f0 ? sny / snx : 1f0
        cornew <= 0f0 && (cornew = 1f-4)
        c.htg_cor_init[sp] = (cornew < 0.0821f0 || cornew > 12.1825f0) ? 0f0 : log(cornew)
    end
    return s
end

# em/findag.f + em/pothtg.f — the CRATET ABIRTH (effective age) dub. SI50 (DF) / SI100 (WB, LP-PP, WL, OS) come from
# the site species' SITEAR as findag.f builds them; EFAGE inverts each species group's height-age curve.
# em/blkdat.f:209-212 PLNJSP — the USDA PLANTS symbols Climate-FVS reads the per-species viability columns by.
const _EM_PLNJSP = String["PIAL","LAOC","PSME","PIFL2","LALY","JUSC2","PICO","PIEN","ABLA","PIPO",
                          "FRPE","POTR5","POBAT","POBA2","PODEM","POAN3","BEPA","2TN","2TB"]
climate_plant_symbols(::EasternMontana) = _EM_PLNJSP

function em_findag_abirth!(s::StandState)
    t = s.trees; p = s.plot
    isisp = Int(p.site_species); sitisp = (1 <= isisp <= length(p.sp_site_index)) ? p.sp_site_index[isisp] : 0f0
    @inbounds for i in 1:t.n
        t.birth_age[i] > 0f0 && continue
        sp = Int(t.species[i]); h = t.height[i]
        si50 = 0f0; si100 = 0f0
        if sp == 3
            si50 = 1.0096f0 + 0.6279f0 * sitisp
            isisp == 3 && (si50 = sitisp)
        end
        if sp == 1 || (7 <= sp <= 10)
            si100 = sitisp
            isisp == 3 && (si100 = (sitisp - 1.0096f0) / 0.6279f0)
        end
        (sp == 2 || sp == 18) && (si100 = p.sp_site_index[sp])
        efage = NaN32
        if sp == 1 || sp == 7 || sp == 2 || sp == 18                  # pothtg.f label 100
            ccf = 125f0
            a = 9.72443f0 - 0.00091f0 * si100 * ccf - h
            b = -0.23733f0 + 0.0149f0 * si100
            cc = 0.00160f0 - 0.00005f0 * si100
            tem = b * b - 4f0 * a * cc; tem < 0f0 && (tem = 0f0)
            efage = cc != 0f0 ? (-b + sqrt(tem)) / (2f0 * cc) : -a / b
        elseif sp == 3                                                  # label 200
            temht = h - 4.5f0; temht == 1f0 && (temht = 1.1f0)
            temsi = si50 - 4.5f0
            term1 = (42.397f0 * temsi^0.3197f0) / temht - 1f0
            term1 < 0f0 && continue
            efage = exp((log(term1) + 1.0232f0 * log(temsi) - 9.7278f0) / (-1.2934f0))
        elseif sp == 8 || sp == 9                                       # label 300
            h >= si100 && continue
            term1 = log((1f0 - ((h / si100)^(1f0 - 0.302381f0))) / 0.931764f0)
            efage = term1 / (-0.01679f0)
        elseif sp == 10                                                 # label 400
            term1 = ((3.635794f0 * si100^0.916307f0) / h) - 1f0
            term1 <= 0f0 && continue
            term2 = log(term1) - 6.09478f0 + 0.277025f0 * log(si100)
            efage = exp(term2 / (-0.96483f0))
        elseif sp == 12 || sp == 17                                     # label 500 (PROB>0 only)
            (t.tpa[i] <= 0f0 || h <= 0f0) && continue
            efage = (h * 2.54f0 * 12f0 / 26.9825f0)^(1f0 / 1.1752f0)
        end
        isnan(efage) && continue
        t.birth_age[i] = efage
        t.age_known[i] = true
    end
    return s
end

# em/regent.f DATA SLO/SHI — LOCAL to regent (the CR/UT relative-site bounds), NOT the blkdat site ranges.
const EM_RG_SLO = Float32[0,0,0,0,0,5,0,0,0,0,30,30,30,30,30,30,30,0,30]
const EM_RG_SHI = Float32[0,0,0,0,0,15,0,0,0,0,120,70,120,120,120,120,70,0,120]
const _EM_RG_BACON = 0.005454154f0
# em/regent.f:359-379 species → sub-model: 1=EMVAR (SMHTGF/SMDGF), 2=NIVAR (LL), 3=TTVAR (LM), 4=CRVAR, 5=UTVAR.
@inline _em_rg_kind(sp::Int) = sp == 5 ? 2 : sp == 4 ? 3 :
    (sp == 11 || (13 <= sp <= 16) || sp == 19) ? 4 : (sp == 6 || sp == 12 || sp == 17) ? 5 : 1

"""EM `small_tree_growth!` — em/regent.f (non-ESTAB entry, `CALL REGENT(.FALSE.,1)`), ported in the Fortran's own
shape: ONE subcycle loop `DO 17 J / DO 16 ISPC / DO 15 I3 (IND1)` over every sub-model, then ONE assembly loop
`DO 30 ISPC / DO 25 I3` with the tripling loop-back at 918. The sub-models share that single pass, so the order of
the RNG draws (EMVAR/TTVAR ZRAND inside the subcycles, NIVAR/CRVAR/UTVAR ZZRAN in the assembly), the running
RDNEXT/BANEXT(J+1) density feedback, and the Fortran variables that CARRY between trees (BARK, H1, D, HTGR) are
all reproduced — the earlier per-sub-model passes got each of those out of step with live FVS on mixed stands."""
function small_tree_growth!(s::StandState, stash, ::EasternMontana; fint::Float32 = 10.0f0,
                            lestb::Bool = false, itrnin::Int = 1,
                            atba::Float32 = -1f0, atccf::Float32 = -1f0, atavh::Float32 = -1f0)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    n = t.n; n == 0 && return s
    rhcon = em_regcons!(s)                            # regent.f:1480 RHCON (NI constant for LL; 1.0/RCOR2 otherwise)
    ba = p.basal_area; relden = p.relative_density; avh = p.avg_height
    dgsd = s.control.dg_sd; regyr = _EM_RG_REGYR
    yr = Float32(s.control.year)
    scale = yr / fint                                 # regent.f:218 SCALE=YR/FINT
    scale2 = fint / yr                                # regent.f:1024 SCALE2=FINT/YR (DUBSCR crown test)
    cur_year = current_cycle_year(s)
    ntyr = trunc(Int, fint); iyr = Int(regyr)         # regent.f:184-200
    lskiph = false
    if lestb                                          # regent.f:186-187 ESTAB: the rest of the cycle after year 5
        ntyr -= 5; lskiph = ntyr <= 0
    end
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1)
    kper = zeros(Int, 10); itot = ntyr; nn = nper
    @inbounds for k in 1:nper
        nn == 1 && break
        kper[k] = itot ÷ nn; itot -= kper[k]; nn -= 1
    end
    nper > 0 && (kper[nper] = itot)
    # regent.f:234-257 subcycle stand density from the LARGE trees (D≥3). CCFCAL returns CCFT·P (ccfcal.f), and FVS
    # evaluates RDNEXT += ((K·CI)/P)·PN with PN=P·0.985**K (integer power) — the P does NOT cancel in single precision.
    banext = zeros(Float32, 11); rdnext = zeros(Float32, 11)
    temba = atba > 0f0 ? atba : ba                    # regent.f:204-208 TEMBA/TEMCCF/TEMAHT (after-thin values)
    temccf = atccf > 0f0 ? atccf : relden
    temaht = atavh
    if lestb                                          # regent.f:263-276 (label 8): interpolate from the cycle start
        bayr = 0f0; ccfyr = 0f0
        if !lskiph
            bayr = (ba - temba) / Float32(itot); ccfyr = (relden - temccf) / Float32(itot)
        end
        nyr = 5
        @inbounds for j in 1:nper
            rdnext[j] = temccf + Float32(nyr) * ccfyr; banext[j] = temba + Float32(nyr) * bayr
            nyr += kper[j]
        end
    else
        @inbounds for j in 1:nper; banext[j] = ba; rdnext[j] = relden; end
    end
    if !lestb && nper > 1
        @inbounds for i in 1:n
            d1 = t.dbh[i]; d1 < 3.0f0 && continue
            sp = Int(t.species[i]); pr = t.tpa[i]
            d2 = d1 + t.diam_growth[i] / em_bratio(sp, d1)
            b1 = _EM_RG_BACON * d1 * d1; b2 = _EM_RG_BACON * d2 * d2
            c1 = em_tree_ccf(sp, d1) * pr; c2 = em_tree_ccf(sp, d2) * pr
            bi = (b2 - b1) / 10.0f0; ci = (c2 - c1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]; pn = pr * fpowi(0.985f0, k)
                rdnext[j] += Float32(k) * ci / pr * pn
                banext[j] += Float32(k) * bi * pn
            end
        end
    end
    wk3 = Float32[t.height[i] for i in 1:n]           # regent.f:301-302 WK3=HT, WK5=DBH
    wk5 = Float32[t.dbh[i] for i in 1:n]
    if lestb                                          # regent.f:287-299 DO 13: crown for each new record, STORAGE order
        @inbounds for i in itrnin:n
            pcc = (pt = Int(t.plot_id[i]); (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0f0)
            crn = 0.89722f0 - 0.0000461f0 * pcc
            ran = 0f0
            while true
                ran = bachlo(s.rng, 0f0, 1f0); (ran < -1f0 || ran > 1f0) && continue
                break
            end
            crn = crn + 0.07985f0 * ran
            crn > 0.90f0 && (crn = 0.90f0); crn < 0.20f0 && (crn = 0.20f0)
            icn = trunc(Int32, (crn * 100f0) + 0.5f0)
            t.crown_pct[i] = icn                 # ICR only: PCT stays estab.f:1252's 0 (NIVAR BAL = full BAJ)
        end
    end
    ah = lestb ? temaht : avh; r = lestb ? temccf : relden   # regent.f:311-317
    delmax = (ah / 36.0f0) * (0.01232f0 * r - 1.75f0); delmax > 0.0f0 && (delmax = 0.0f0)
    x = ah * (r / 100.0f0); x > 300.0f0 && (x = 300.0f0)             # regent.f:324 X=AH*(R/100.)
    pctred = 1.11436f0 + x*(-0.011493f0 + x*(0.43012f-4 + x*(-0.72221f-7 + x*(0.5607f-10 - x*0.1641f-13))))
    pctred > 1.0f0 && (pctred = 1.0f0); pctred < 0.01f0 && (pctred = 0.01f0)
    # species-major ISCT/IND1 order (DO ISPC / DO I3=I1,I2 / I=IND1(I3)) — every REGENT loop walks it
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    order = Int[]
    @inbounds for sp in 1:MAXSP
        i1 = isct[sp, 1]; i1 == 0 && continue
        for k in i1:isct[sp, 2]
            (1 <= k <= length(ind1)) || continue
            ii = Int(ind1[k]); (1 <= ii <= n) && push!(order, ii)
        end
    end
    pccf_of(i) = (pt = Int(t.plot_id[i]); (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0f0)
    bark_c = NaN32                                    # Fortran BARK — carried between trees (see the CR/UT DGK)
    h1_c = 0f0                                        # Fortran H1 — carried into the assembly loop's NIVAR RELH
    dadj_c = 0f0
    # ---- subcycle loop (regent.f:337-679) ----
    ky = 0
    @inbounds for j in 1:nper
        ky += kper[j]
        baj = banext[j]; rdj = rdnext[j]
        surv = fpowi(0.985f0, ky)
        ppccf = relden > 0.0f0 ? 1.0f0 + (rdj - relden) / relden : 0.0f0
        kfrac = Float32(kper[j]) / regyr
        for i in order
            sp = Int(t.species[i]); kind = _em_rg_kind(sp)
            (kind == 4 || kind == 5) && j > 1 && continue      # regent.f:382 CR/UT: one pass only
            (kind == 1 || kind == 3 || kind == 4) && lskiph && continue   # regent.f:378
            xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
            xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
            hcor = c.htg_cor_small[sp]
            con = kind == 2 ? rhcon[sp] + hcor : rhcon[sp] * fexp(hcor)   # regent.f:386-390
            d = t.dbh[i]; h = t.height[i]
            d >= EM_RG_XMAX[sp] && continue
            lestb && i < itrnin && continue
            pr = t.tpa[i]
            h1 = wk3[i]; d1 = wk5[i]; h1_c = h1
            bal = baj * (100.0f0 - t.crown_ratio[i]) * 0.0001f0
            c1 = em_tree_ccf(sp, d) * pr                        # regent.f:445 CCFCAL(ISPC,D,H1,…)
            b1 = _EM_RG_BACON * d * d
            bark_c = em_bratio(sp, d)       # regent.f:447 BARK=BRATIO(ISPC,D,H)
            rsimod = 0f0; pothtg = 0f0
            if kind == 4 || kind == 5
                si = p.sp_site_index[sp]
                si > EM_RG_SHI[sp] && (si = EM_RG_SHI[sp])
                si <= EM_RG_SLO[sp] && (si = EM_RG_SLO[sp] + 0.5f0)
                relsi = (si - EM_RG_SLO[sp]) / (EM_RG_SHI[sp] - EM_RG_SLO[sp])
                rsimod = 0.5f0 * (1.0f0 + relsi)
                if kind == 4
                    pothtg = p.sp_site_index[sp] / (15.0f0 - 4.0f0 * relsi) * 1.0f0
                elseif sp == 6
                    pothtg = (si / 10.0f0) * (si * 1.5f0 - h1) / (si * 1.5f0)
                end
            end
            pccf = pccf_of(i)
            tpccf = pccf * ppccf; tpccf > 300.0f0 && (tpccf = 300.0f0); tpccf < 25.0f0 && (tpccf = 25.0f0)
            cr = Float32(t.crown_pct[i])
            htgrl = 0f0; h2 = h1
            if kind == 2                                         # NIVAR (regent.f:455-471)
                relh = abs(ah - 4.5f0) < 0.01f0 ? 0.0f0 : (h1 - 4.5f0) / (ah - 4.5f0)
                relh > 1.0f0 && (relh = 1.0f0); relh < 0.0f0 && (relh = 0.0f0)
                dadj_c = delmax * relh * relh - 2.0f0 * delmax * relh + 0.65f0
                htgrl = con + _EM_RG_BH * flog(h1) + _EM_RG_BCCF * rdj + _EM_RG_BBAL * bal
                h2 = h1 + fexp(htgrl) * kfrac * xrhgro
            elseif kind == 3                                     # TTVAR (regent.f:476-508): persistent ZRAND
                if t.zrand[i] == -999f0
                    while true
                        z = bachlo(s.rng, 0.0f0, 1.0f0); (z < -2.0f0 || z > 2.0f0) && continue
                        t.zrand[i] = z; break
                    end
                end
                if d <= 0.0f0
                    htgrl = 0.0f0
                else
                    beta1 = fexp(1.17527f0 - 0.42124f0 * flog(tpccf))
                    beta2 = fexp(-2.56002f0 - 0.58642f0 * flog(tpccf))
                    htg1 = beta1 + beta2 * cr
                    stddev = htg1 * (1.08720f0 - 0.00230f0 * cr)
                    htgrl = htg1 + t.zrand[i] * stddev
                    htgrl <= 0.1f0 && (htgrl = 0.1f0; t.zrand[i] = -999f0)
                end
                h2 = h1 + htgrl * kfrac * xrhgro * con
            elseif kind == 4 || kind == 5                        # CR/UT (regent.f:512-535)
                xv = cr / 100.0f0
                vigor = (150.0f0 * (xv^3) * fexp(-6.0f0 * xv)) + 0.3f0
                vigor > 1.0f0 && (vigor = 1.0f0)
                sp == 6 && (vigor = 1.0f0 - ((1.0f0 - vigor) / 3.0f0))
                if sp == 12 || sp == 17                          # Sheppard aspen/PB
                    ab = Float32(t.birth_age[i])
                    hite1 = 26.9825f0 * fpow(ab, 1.1752f0)
                    hite2 = 26.9825f0 * fpow(ab + 10.0f0, 1.1752f0)
                    htgrl = (hite2 - hite1) / (2.54f0 * 12.0f0) * rsimod * con
                    htgrl = htgrl * 0.75f0
                else
                    htgrl = pothtg * pctred * vigor * con
                end
                h2 = h1 + htgrl * fint / 10.0f0
            else                                                 # EMVAR (regent.f:544-545 → em/smhtgf.f)
                if t.zrand[i] == -999f0
                    while true
                        z = bachlo(s.rng, 0.0f0, 1.0f0); (z < -2.0f0 || z > 2.0f0) && continue
                        t.zrand[i] = z; break
                    end
                end
                local htgrr::Float32
                if d <= 0.0f0
                    htgrr = 0.0f0
                else
                    beta1 = fexp(EM_B0ACCF[sp] + EM_B1ACCF[sp] * flog(tpccf))
                    beta2 = fexp(EM_B0BCCF[sp] + EM_B1BCCF[sp] * flog(tpccf))
                    htg1 = beta1 + beta2 * cr
                    stddev = htg1 * (EM_B0ASTD[sp] + EM_B1BSTD[sp] * cr)
                    htgrr = htg1 + t.zrand[i] * stddev
                    htgrr <= 0.1f0 && (htgrr = 0.1f0; t.zrand[i] = -999f0)
                end
                h2 = h1 + htgrr * kfrac * xrhgro * con
            end
            wk3[i] = h2
            # subcycle diameter + density feedback (regent.f:571-660)
            d2 = d
            if kind == 2
                (j >= nper || d >= 3.0f0 || h2 <= 4.5f0) && continue          # GO TO 14 (no WK5, no density)
                d1v = EM_RG_DIAM[sp] + dadj_c
                h1 > 4.5f0 && (d1v = _EM_RG_AX * fpow(h1 - 4.5f0, _EM_RG_BX) + dadj_c)
                d2v = _EM_RG_AX * fpow(h2 - 4.5f0, _EM_RG_BX) + dadj_c
                dgj = (d2v - d1v) * xrdgro; dgj < 0.0f0 && (dgj = 0.0f0)
                d2 = d + dgj
            elseif kind == 3
                hless4 = h2 - 4.5f0
                dless3 = 0.000231f0 * hless4 * cr - 0.00005f0 * hless4 * pccf + 0.001711f0 * cr + 0.17023f0 * hless4
                d2 = (dless3 + 0.3f0) * xrdgro
                d2 < EM_RG_DIAM[sp] && (d2 = EM_RG_DIAM[sp])
                wk5[i] = d2
                (j >= nper || d >= 3.0f0 || h2 <= 4.5f0) && continue
            elseif kind == 4
                d >= 1.0f0 && continue                                          # GO TO 15
                h2 <= 4.5f0 && (d2 = d + 0.0001f0 * h2)
            elseif kind == 5
                ((sp == 12 || sp == 17) && d >= 3.0f0) && continue
                h2 <= 4.5f0 && (d2 = d + 0.001f0 * h2)
            else
                d2 = _em_smdgf(sp, h2, cr, pccf) * xrdgro                       # SMDGF on RAW PCCF, EVERY height
                d2 < EM_RG_DIAM[sp] && (d2 = EM_RG_DIAM[sp])
            end
            wk5[i] = d2
            (kind == 3 && h2 < 4.5f0) && continue
            if kind == 1
                if j < nper && d < 3.0f0 && h2 > 4.5f0
                    c1e = em_tree_ccf(sp, d1) * pr; c2 = em_tree_ccf(sp, d2) * pr
                    rdnext[j+1] += Float32(ky) * (c2 - c1e) / 10.0f0 * surv
                    banext[j+1] += (_EM_RG_BACON * (d2 * d2 - d1 * d1)) * pr * surv
                end
            else
                c2 = em_tree_ccf(sp, d2) * pr
                rdnext[j+1] += Float32(ky) * (c2 - c1) / 10.0f0 * surv
                banext[j+1] += (_EM_RG_BACON * d2 * d2 - b1) * pr * surv
            end
        end
    end
    # ---- assembly (regent.f DO 30 ISPC / DO 25 I3, :686-1050) ----
    rmai_v = _em_rmai(s)
    ltrip = !lestb && stash !== nothing                # regent.f:1041 no tripling from ESTAB
    @inbounds for i in order
        sp = Int(t.species[i]); kind = _em_rg_kind(sp)
        xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
        xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
        xmx = EM_RG_XMAX[sp]; xmn = EM_RG_XMIN[sp]; diam = EM_RG_DIAM[sp]
        dgmx = kind == 1 ? fint * 0.2f0 : kind == 3 ? fint * EM_RG_DGMAX[sp] :
               kind == 4 ? EM_RG_DGMAX[sp] * fint / 10.0f0 : EM_RG_DGMAX[sp]
        d = t.dbh[i]; cr = Float32(t.crown_pct[i])
        d >= xmx && continue
        lestb && i < itrnin && continue
        h = t.height[i]; hk = wk3[i]
        dk = wk5[i]; htgr = hk - h
        (kind == 1 || kind == 3) && d < diam && (d = diam)             # regent.f:751 D1=DIAM(ISPC)
        if lskiph && kind == 1                                          # regent.f:760-770 (ESTAB, FINT<=5)
            if h >= 4.5f0
                dkl = _em_smdgf(sp, h, cr, pccf_of(i)); dkl < d && (dkl = d)
                t.dbh[i] = dkl; t.diam_growth[i] = dkl > dgmx ? dgmx : dkl
            end
            continue
        elseif lskiph && kind == 3                                      # regent.f:774-785
            if h >= 4.5f0
                hl = h - 4.5f0
                dkl = (0.000231f0*hl*cr - 0.00005f0*hl*pccf_of(i) + 0.001711f0*cr + 0.17023f0*hl + 0.3f0) * xrdgro
                dkl < d && (dkl = d)
                t.dbh[i] = dkl; t.diam_growth[i] = dkl > dgmx ? dgmx : dkl
            end
            continue
        end
        htgr1 = hk - h; htgr1 < 0.0f0 && (htgr1 = 0.0f0)
        icrk = Int(t.crown_pct[i])
        dbh0 = t.dbh[i]                                                 # DBH(K): DGDRIV copied it to the slots
        pccf = pccf_of(i)
        large = (t.ht_growth[i],
                 ltrip ? (stash.htg_copy[i] ? stash.htgU[i] : t.ht_growth[i]) : 0f0,
                 ltrip ? (stash.htg_copy[i] ? stash.htgL[i] : t.ht_growth[i]) : 0f0)
        dg_in = (t.diam_growth[i], ltrip ? stash.dgU[i] : 0f0, ltrip ? stash.dgL[i] : 0f0)
        nrec = (ltrip && kind in (2, 4, 5)) ? 3 : 1
        for l in 0:(nrec - 1)
            dbhk = dbh0; dgk_out = dg_in[l+1]; dbh_set = false
            if lskiph && kind == 2                                   # regent.f:771-773 HTG(K)=0, GO TO 20
                htgr = 0f0
            elseif kind == 1
                hk = wk3[i]; dk = wk5[i]; htgr = hk - h
            elseif kind == 2
                zz = 0.0f0
                while true
                    zz = dgsd >= 1.0f0 ? bachlo(s.rng, 0.0f0, 1.0f0) : 0.0f0
                    (zz > 1.0f0 || zz < -1.5f0) && continue
                    break
                end
                htgr = htgr1 * fexp(zz * _EM_RG_HSIGMA)
            elseif kind == 4 || kind == 5
                zz = 0.0f0
                while true
                    zz = dgsd >= 1.0f0 ? bachlo(s.rng, 0.0f0, 1.0f0) : 0.0f0
                    (zz > 0.5f0 || zz < -2.0f0) && continue
                    break
                end
                kind == 5 && (htgr = (htgr + zz * 0.1f0) * xrhgro)
                kind == 4 && (htgr = (htgr + zz * 0.2f0) * xrhgro)
                htgr < 0.1f0 && (htgr = 0.1f0)
            end
            xwt = (d <= xmn || lestb) ? 0.0f0 : (d - xmn) / (xmx - xmn)   # regent.f:833-835 XWT=0 from ESTAB
            htgk = htgr * (1.0f0 - xwt) + xwt * large[l+1]
            (kind == 4 && htgk < 0.1f0) && (htgk = 0.1f0)
            cap = s.control.sp_size_cap[sp, 4]
            if h + htgk > cap
                htgk = cap - h; htgk < 0.1f0 && (htgk = 0.1f0)
            end
            # diameter (regent.f:852-1012)
            skip23 = (kind == 4 && d >= 1.0f0) || (!(kind == 5 && sp == 6) && d >= 3.0f0)
            if !skip23
                d1v = diam
                if kind == 2
                    relh = abs(ah - 4.5f0) < 0.01f0 ? 0.0f0 : (h1_c - 4.5f0) / (ah - 4.5f0)   # stale H1 (em only)
                    relh > 1.0f0 && (relh = 1.0f0); relh < 0.0f0 && (relh = 0.0f0)
                    dadj_c = delmax * relh * relh - 2.0f0 * delmax * relh + 0.65f0
                    d1v = diam + dadj_c
                    h > 4.5f0 && (d1v = _EM_RG_AX * fpow(h - 4.5f0, _EM_RG_BX) + dadj_c)
                end
                hkk = hk; hk = h + htgk; kind == 1 && (hk = hkk)
                if kind != 3 && hk < 4.5f0
                    kind == 2 && (dbhk = 0.1f0 + diam * 0.01f0 + (h + htgk) * 0.001f0; dbh_set = true)
                    kind == 4 && (dbhk = d + 0.001f0 * (h + htgk); dbh_set = true)
                    dgk_out = 0.0f0
                elseif kind == 1
                    dkk = _em_smdgf(sp, h, cr, pccf)
                    bark_c = em_bratio(sp, d)
                    dg = (dk - dkk) * bark_c
                    dds = dg * (2.0f0 * bark_c * d + dg) * scale
                    dg = sqrt(max((d * bark_c)^2 + dds, 0f0)) - bark_c * d
                    if lestb                                                # regent.f:895-898
                        dbhk = dk; dbh_set = true
                        dg > dgmx && (dg = dgmx)
                    end
                    ((dbhk + dg) < diam && hk >= 4.5f0) && (dg = diam - dbhk)
                    dgk_out = dg
                else
                    dkk = 0f0; dgt = 0f0
                    if kind == 2
                        dk = _EM_RG_AX * fpow(hk - 4.5f0, _EM_RG_BX) + dadj_c
                        dk < diam && (dk = diam)
                        dk = dk + hk * 0.001f0
                    elseif kind == 4 || kind == 5
                        if sp == 6
                            sitear = p.sp_site_index[sp]
                            dk = (hk - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
                            dkk = (h - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0)
                            h < 4.5f0 && (dkk = d)
                        else
                            bx = EM_RG_HT2[sp]                               # em/blkdat.f HT2 (COEFFS)
                            ax = c.ht_dbh_iabflg[sp] == 1 ? EM_RG_HT1[sp] : c.ht_dbh_aa[sp]
                            dk = (bx / (flog(hk - 4.5f0) - ax)) - 1.0f0; dk < 0.1f0 && (dk = 0.1f0)
                            dkk = h <= 4.5f0 ? d : (bx / (flog(h - 4.5f0) - ax)) - 1.0f0
                        end
                    else                                                    # TTVAR
                        if hk >= 4.5f0
                            hless4 = h - 4.5f0
                            dless3 = 0.000231f0 * hless4 * cr - 0.00005f0 * hless4 * pccf + 0.001711f0 * cr + 0.17023f0 * hless4
                            dkk = dless3 + 0.3f0; dkk < diam && (dkk = diam)
                            bark_c = em_bratio(sp, d)
                            dgt = (dk - dkk) * bark_c
                            dds = dgt * (2.0f0 * bark_c * d + dgt)
                            dgt = sqrt(max((d * bark_c)^2 + dds, 0f0)) - bark_c * d
                        else
                            dgt = 0.0f0
                        end
                    end
                    if lestb                                                # regent.f:966-974 ESTAB diameter
                        dbhk = dk; dbh_set = true
                        dg = kind == 3 ? dgt : dgk_out
                        (kind == 3 && dg > dgmx) && (dg = dgmx)
                        (kind == 2 || kind == 4) && (dg = dk)
                        if kind == 4 || kind == 5
                            dbhk < diam && (dbhk = diam)
                            dbhk = dbhk + 0.001f0 * hk
                            dg = dbhk
                        end
                        (dbhk + dg) < diam && (dg = diam - dbhk)
                        dgk_out = dg
                    else
                    dgk = 0.0f0
                    if kind == 2
                        dgk = (dk - d1v) * xrdgro
                    elseif kind == 4 || kind == 5
                        if dk < 0.0f0 || dkk < 0.0f0
                            dgk = htgk * 0.2f0 * bark_c * xrdgro          # stale BARK (regent.f:989)
                            dk = d + dgk
                        else
                            dgk = (dk - dkk) * bark_c * xrdgro            # stale BARK (regent.f:992)
                        end
                        dgk > dgmx && (dgk = dgmx)
                    end
                    dgk < 0.0f0 && (dgk = 0.0f0)
                    bark_c = em_bratio(sp, dbhk)     # regent.f:1002 BRATIO(ISPC,DBH(K),HT(K))
                    dg = kind == 2 ? dgk * bark_c : (kind == 3 ? dgt : dgk)
                    dds = dg * (2.0f0 * bark_c * d + dg) * scale
                    dg = sqrt(max((d * bark_c)^2 + dds, 0f0)) - bark_c * d
                    (dbhk + dg) < diam && (dg = diam - dbhk)
                    dgk_out = dg
                    end
                end
                # label 22: DGBND
                dgk_out = dg_bound(nothing, nothing, sp, dbhk, dgk_out, s.control.sp_size_cap)
                if kind != 1 && !lestb                                       # regent.f:1023-1041 DUBSCR crown
                    d = dbhk
                    bark_c = em_bratio(sp, d)
                    dds2 = (dgk_out * (2.0f0 * bark_c * d + dgk_out)) * scale2
                    dg2 = sqrt(max((d * bark_c)^2 + dds2, 0f0)) - bark_c * d
                    dg2 < 0.0f0 && (dg2 = 0.0f0)
                    dnew = d + dg2
                    if dnew >= 3.0f0
                        if icrk == 0
                            crd = em_dubscr(s.rng, sp, dnew, hk, ba, pccf, avh, rmai_v, dgsd)
                            icrk = trunc(Int, crd * 100.0f0 + 0.5f0)
                        end
                        l == 0 && (t.crown_pct[i] = icrk)
                    end
                end
            end
            # store record K (central l=0 → the tree; copies → the TRIPLE stash)
            if l == 0
                t.ht_growth[i] = htgk
                if !skip23
                    t.diam_growth[i] = dgk_out
                    dbh_set && (t.dbh[i] = dbhk)
                end
            else
                if l == 1
                    stash.htgU[i] = htgk; stash.is_small[i] = true
                    skip23 || (stash.dgU[i] = dgk_out; stash.dbhU[i] = dbhk)
                else
                    stash.htgL[i] = htgk
                    skip23 || (stash.dgL[i] = dgk_out; stash.dbhL[i] = dbhk)
                end
            end
        end
        # EMVAR/TTVAR tripling (regent.f:1043-1056): the two copies are exact duplicates of the central record
        if ltrip && (kind == 1 || kind == 3)
            stash.htgU[i] = t.ht_growth[i]; stash.htgL[i] = t.ht_growth[i]; stash.is_small[i] = true
            stash.dgU[i] = t.diam_growth[i]; stash.dgL[i] = t.diam_growth[i]
        end
    end
    return s
end

# em/blkdat.f HHTMAX — ESGENT's cap on a newly established tree's height (and then DBH=2.95).
const EM_HHTMAX = Float32[23,27,21,27,18,6,24,18,18,17,16,16,16,16,16,16,16,22,16]

"""
    em_esgent!(s, nstart; fint, atba, atccf, atavh)

em/esgent.f: SPESRT, then REGENT(.TRUE.,ITRNIN) for the records established this cycle — the same Fortran-shaped
REGENT as the growth pass, in its ESTAB mode (the crown draw per new record in storage order; FINT−5 years of
subcycling with the density interpolated from the after-thin TEMBA/TEMCCF/TEMAHT; XWT=0; the ESTAB diameter; no
tripling, no DUBSCR). Then per record HTG·WK4 (the establishment height multiplier HTIMLT: 1 for PLANT, 0.4 for
natural regeneration), HT+=HTG, the WK4<1 DBH rescale, and the HHTMAX cap (DBH=2.95). jl formerly grew only the
EMVAR conifers with a local ZRAND, so the next cycle's REGENT redrew it and planted-stand ForTyp/SizeCls diverged.
"""
function em_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0,
                    atba::Float32 = -1f0, atccf::Float32 = -1f0, atavh::Float32 = -1f0)
    t = s.trees
    nstart >= t.n && return s
    species_sort!(s)                                   # esgent.f:47 CALL SPESRT
    small_tree_growth!(s, nothing, s.variant; fint = fint, lestb = true, itrnin = nstart + 1,
                       atba = atba, atccf = atccf, atavh = atavh)
    @inbounds for i in (nstart+1):t.n
        sp = Int(t.species[i])
        htemp = t.height[i] + t.ht_growth[i]
        wk4 = t.htimlt[i]
        t.ht_growth[i] = t.ht_growth[i] * wk4
        t.height[i] = t.height[i] + t.ht_growth[i]
        if wk4 < 1f0
            if t.height[i] < 4.5f0
                t.dbh[i] = 0.1f0 + 0.001f0 * t.height[i]; t.diam_growth[i] = 0f0
            else
                t.dbh[i] = t.dbh[i] * (t.height[i] / htemp)
                t.diam_growth[i] = t.diam_growth[i] * (t.height[i] / htemp)
            end
        end
        if t.height[i] > EM_HHTMAX[sp]
            t.height[i] = EM_HHTMAX[sp]; t.dbh[i] = 2.95f0
        end
    end
    esgent_add_gentim!(s, nstart, fint)                # estab.f:1504 ABIRTH += GENTIM (after ESGENT)
    return s
end
