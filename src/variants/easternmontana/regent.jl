# =============================================================================
# regent.jl (easternmontana) — EM small-tree growth (em/regent.f, NIVAR path).
#
# Small trees (dbh < XMAX[sp]) — for emt01, the 0.1" DF seedling. Mirrors KT's subcycle scaffold
# (REGYR=5 subcycles, per-subcycle density from the large trees) with EM's NIVAR forms:
#   HTGRL = CON + BH·ln(H1) + BCCF·RDJ + BBAL·BAL   (log-space; BH=0.3740, BCCF=−0.00391, BBAL=−0.22957)
#   H2    = H1 + EXP(HTGRL)·(KPER/REGYR)·XRHGRO     (em/regent.f:550 — increment is exp(HTGRL))
#   CON   = RHCON[sp] + HCOR[sp];  RHCON = REGCH + 1.0667 + RHHAB[MAPHAB[ITYPE]]  (em_regcons!)
#   REGCH = RHGL[IGL] + (RSAB0 + RSAB1·cos(ASP) + RSAB2·sin(ASP))·SLOPE
#   DG-dub (per subcycle, D<3 & H2>4.5): D1=DIAM[sp]+DADJ (or AX·(H1−4.5)^BX+DADJ if H1>4.5);
#     D2=AX·(H2−4.5)^BX+DADJ; DGJ=max(D2−D1,0)·XRDGRO; D2=D+DGJ. AX=0.0658, BX=1.3817.
#   DADJ = DELMAX·RELH²−2·DELMAX·RELH+0.65; DELMAX=min((AH/36)·(0.01232·RELDEN−1.75),0); RELH=(H1−4.5)/(AH−4.5).
# Final: HTGR1=WK3−H + ZZRAN·HSIGMA(0.59) (dgsd≥1); XWT blend d∈[XMIN,cap] with the large-tree HTG; size cap.
# NIVAR handles the EMVAR conifers + LL(5) — LL capped at D<3 (regent.f:858) so its D≥3 large trees stay on
# the pure large-tree path (validated: em_LLseed seedlings grow, em_LL large trees bit-exact).
# CRVAR(CO)/TTVAR(LM)/UTVAR(RM,AS,PB) small-tree branches deferred.
# =============================================================================

const EM_RG_XMAX = Float32[3,3,3,3,10,99,3,3,3,3,2,4,2,2,2,2,4,3,3]
const EM_RG_XMIN = Float32[1.5,1.5,1.5,1.5,2,90,1.5,1.5,1.5,1.5,0.5,2,0.5,0.5,0.5,0.5,2,1.5,1.5]
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

# em/regent.f DO-25 CRVAR/UTVAR height→diameter DG (lines 914-1012): the small CR/UT tree whose blended new
# height hk=h+htg has crossed breast height (hk≥4.5) grows DIAMETER from the INVERSE Wykoff height-diameter
# curve DK=(HT2/(ln(hk-4.5)-HT1))-1, NOT the (near-zero) large-tree DG. Juniper(6) uses the SITEAR linear form.
# This is the seam that was MISSING in jl: aspen/PB (UTVAR) and cottonwood/OH (CRVAR) seedlings that grew past
# 4.5' kept the ~0 large-tree DG ⇒ dense multi-species stands under-grew BA ~4× (measured stand 488929762126144:
# jl BA 62 vs oracle 244). Returns the bark-transformed DBH increment (inches), same units as t.diam_growth.
# BARK is the local BRATIO(sp,d,h) — FVS carries the previous DO-25 tree's BARK into the DGK scaling (line 990),
# but line 1001 recomputes it fresh for the DDS transform; for the tiny seedlings that dominate the fix the two
# are ~equal, so a local BRATIO is used (residual cornered to em/regent.f's carried-BARK order-dependence).
@inline function _em_crut_dg(sp::Int, d::Float32, h::Float32, hk::Float32, htg::Float32,
                             sitear::Float32, bark::Float32, dgmx::Float32, scale::Float32,
                             ax::Float32 = EM_RG_HT1[sp])::Float32
    xrdgro = 1.0f0
    local dk::Float32, dkk::Float32
    if sp == 6                                            # RM juniper — SITEAR linear (em/regent.f:915-920)
        dk = (hk - 4.5f0) * 10f0 / (sitear - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
        dkk = (h - 4.5f0) * 10f0 / (sitear - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0); h < 4.5f0 && (dkk = d)
    else                                                  # inverse-Wykoff HD (em/regent.f:922-934). AX (passed in) is
        # the CRATET-fitted HT-DBH intercept AA when IABFLG=0, else the HT1 default. em/regent.f:924-927 selects
        # AX=AA(ISPC) when IABFLG(ISPC)=0 (cratet.f:325-331 fits AA=mean(ln(H−4.5)−HT2/(D+1)) over the species'
        # ≥3 measured D≥3 H>4.5 trees, then sets IABFLG=0 if AA≥0). jl previously hard-coded AX=HT1 here, so every
        # aspen/PB/CO/hardwood small-tree (D<3 UTVAR, D<1 CRVAR) DK used the DEFAULT intercept instead of the
        # stand-fitted one — on FIA aspen stands AA≈4.10 vs HT1=4.4421 ⇒ DK under-predicted ⇒ regent DBH growth
        # ~30% low ⇒ QMD under + the mortality Hamilton vigor G collapsed ⇒ ~2× first-cycle mortality over-kill
        # ⇒ aspen BA under-projected ~30% (stand 31438934010690 y2018 BA 103 vs FVSem_g16 119; jl-fixed = 119).
        # jl already fits this AA in dub_missing_heights! (ht_dbh_aa/ht_dbh_iabflg); this just consumes it, as IE does.
        bx = EM_RG_HT2[sp]
        dk = bx / (log(hk - 4.5f0) - ax) - 1.0f0; dk < 0.1f0 && (dk = 0.1f0)
        dkk = h <= 4.5f0 ? d : bx / (log(h - 4.5f0) - ax) - 1.0f0
    end
    local dgk::Float32
    if dk < 0f0 || dkk < 0f0
        dgk = htg * 0.2f0 * bark * xrdgro
    else
        dgk = (dk - dkk) * bark * xrdgro
    end
    dgk > dgmx && (dgk = dgmx); dgk < 0f0 && (dgk = 0f0)
    dg = dgk                                              # CRVAR/UTVAR: DG(K)=DGK (no *BARK — em/regent.f:1005)
    dds = dg * (2f0 * bark * d + dg) * scale
    dg = sqrt((d * bark)^2 + dds) - bark * d
    (d + dg) < EM_RG_DIAM[sp] && (dg = EM_RG_DIAM[sp] - d)   # DIAM floor (em/regent.f:1010-1012)
    dg < 0f0 && (dg = 0f0)
    return dg
end

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

# Species handled by the current NIVAR-form regent: EMVAR conifers + LL(5) (NIVAR — same HTGRL form; CON
# RHCON+HCOR==RHCON·exp(HCOR) at HCOR=0; BH/BCCF/BBAL are the shared subalpine-fir values). LM/CO/RM/AS/PB deferred.
@inline _em_rg_nivar(sp::Int) = sp == 5   # ONLY LL(5) is NIVAR/exp-form (em/regent.f:364). EMVAR conifers
# {1,2,3,7,8,9,10,18} use SMHTGF+SMDGF (below), NOT the exp-form (was wrongly lumped here ⇒ ~2-3× HTG over-grow).
# Effective regent DIAMETER cap: NIVAR skips the regent for D≥3 (em/regent.f:858 GO TO 23). The conifers
# already have EM_RG_XMAX=3, but LL(5)'s XMAX=10 is only the height-blend range — its regent still caps at 3,
# so D≥3 LL stays on the pure large-tree path (else the XWT DG-override halves the large DG → em_LL breaks).
@inline _em_rg_cap(sp::Int) = sp == 5 ? 3.0f0 : EM_RG_XMAX[sp]
# em/regent.f:322-333 — the CRVAR/UTVAR small-tree HEIGHT-growth DENSITY modifier PCTRED. FVS computes it ONCE
# per stand from X = AH·(RELDEN/100) (top height × relative density), NOT the relative SITE INDEX. jl had wrongly
# fed X = relsi·100 (RELSI·100), which over-stated PCTRED in a dense young stand (measured FVSem_g16 9866226020004
# @cyc2: X=61.49 ⇒ oracle PCTRED 0.554 vs jl relsi·100=75.7 ⇒ 0.757, a 1.37× HTG inflation that COMPOUNDS across
# cycles ⇒ CRVAR (GA/OH) seedlings cross 4.5' too fast ⇒ the inverse-Wykoff DBH jump over-fires ⇒ sub-1" DG
# 0.75 vs oracle 0.11 ⇒ inflated DQ10 ⇒ Hamilton mortality under-kill ⇒ +416 TPA whole-stand). AH=AVH, R=RELDEN.
@inline function _em_rg_pctred(ah::Float32, relden::Float32)::Float32
    x = ah * relden / 100f0; x > 300f0 && (x = 300f0)
    pr = 1.11436f0 + x*(-0.011493f0 + x*(0.43012f-4 + x*(-0.72221f-7 + x*(0.5607f-10 - x*0.1641f-13))))
    pr > 1.0f0 && (pr = 1.0f0); pr < 0.01f0 && (pr = 0.01f0)
    return pr
end
@inline _em_rg_crvar(sp::Int) = sp == 11 || (13 <= sp <= 16) || sp == 19
# UTVAR (UT-variant) small-tree form: RM(6, juniper, XMAX=99 fully-regent), AS(12)/PB(17, aspen Sheppard).
@inline _em_rg_utvar(sp::Int) = sp == 6 || sp == 12 || sp == 17
# Push the regent DG/HTG into the tripling stash so the upper/lower sub-records use the REGENT growth, not
# the stale large-tree dgf DG (TT does this at teton/regent.jl:216). Critical for LM (Wykoff large DG) seedlings.
@inline function _em_rg_stash!(stash, t, i::Int)
    if stash !== nothing && !isempty(stash.dgU) && i <= length(stash.dgU)
        stash.dgU[i] = t.diam_growth[i]; stash.dgL[i] = t.diam_growth[i]
        stash.htgU[i] = t.ht_growth[i]; stash.htgL[i] = t.ht_growth[i]
        !isempty(stash.is_small) && (stash.is_small[i] = true)
    end
end
@inline function _em_dless3(h::Float32, cr::Float32, pccf::Float32)::Float32   # TTVAR DBH-from-height (r:584-587)
    hl = h - 4.5f0
    d = 0.000231f0*hl*cr - 0.00005f0*hl*pccf + 0.001711f0*cr + 0.17023f0*hl + 0.3f0
    return d < EM_RG_DIAM[4] ? EM_RG_DIAM[4] : d
end

# em/regent.f RCON entry: RHCON[sp] = REGCH + 1.0667 + RHHAB[MAPHAB[ITYPE]]. Returns the RHCON vector.
function em_regcons!(s::StandState)
    p = s.plot
    itype = Int(p.habitat_input); (itype < 1 || itype > 30) && (itype = 1)
    igl = Int(p.geo_location); (igl < 1 || igl > 3) && (igl = 1)
    irhhab = EM_RG_MAPHAB[itype]; (irhhab < 1 || irhhab > 5) && (irhhab = 5)
    regch = EM_RG_RHGL[igl] + (_EM_RG_RSAB0 + _EM_RG_RSAB1 * cos(p.aspect) + _EM_RG_RSAB2 * sin(p.aspect)) * p.slope
    base = regch + 1.0667f0 + EM_RG_RHHAB[irhhab]
    rhcon = fill(base, 19)                                        # site constant is species-independent here
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
    # (A) ABIRTH dub from CURRENT height for the aspen/PB pair (pothtg.f:158-161; gated ABIRTH<=0, PROB>0).
    @inbounds for sp in (12, 17), i in 1:n
        Int(t.species[i]) == sp || continue
        (t.tpa[i] <= 0f0 || t.birth_age[i] > 0f0) && continue
        h = t.height[i]; h <= 0f0 && continue
        t.birth_age[i] = (h * 2.54f0 * 12f0 / 26.9825f0)^(1f0 / 1.1752f0)
        t.age_known[i] = true
    end
    # (B) em/regent.f:1195-1360 REGCAL — the small-tree HEIGHT self-calibration, for EVERY species (jl ran it
    # for the aspen pair only). Per species: EDH is the model's predicted increment over the sub-periods and
    # CORNEW = mean(HTG·SCALE3) / mean(EDH) ⇒ HCOR = ln(CORNEW), gated at NCALHT=5 records and trapped to
    # ±2.5σ in ln(C). The EMVAR branch calls SMHTGF, which DRAWS the persistent per-tree ZRAND — those draws
    # land HERE, before DGDRIV, which is why jl's whole main RNG stream sat ahead of live on EM stands.
    rhcon = em_regcons!(s)
    ba = p.basal_area; relden = p.relative_density
    reldm1 = p.relative_density_prev; reldm1 <= 0f0 && (reldm1 = relden)
    temccf = relden                                     # regent.f:206-207 TEMCCF = ATCCF, else RELDEN
    fint = Float32(s.control.year)
    ntyr = Int(round(fint)); iyr = Int(_EM_RG_REGYR)
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1); nper < 1 && (nper = 1)
    kper = zeros(Int, nper); itot = ntyr; nn = nper
    @inbounds for k in 1:nper
        if nn == 1; kper[k] = itot; break; end
        kper[k] = itot ÷ nn; itot -= kper[k]; nn -= 1
    end
    banext = fill(ba, nper); rdnext = fill(relden, nper)
    pctred = _em_rg_pctred(p.avg_height, relden)
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
            pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0f0
            edh = 0f0
            for jj in 1:nper
                baj = banext[jj]; rdj = rdnext[jj]
                bal = baj * (100f0 - Float32(t.crown_ratio[i])) * 0.0001f0
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

function small_tree_growth!(s::StandState, stash, ::EasternMontana; fint::Float32 = 10.0f0)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    n = t.n; n == 0 && return s
    rhcon = em_regcons!(s)
    ba = p.basal_area; relden = p.relative_density; avh = p.avg_height; ah = avh
    dgsd = s.control.dg_sd; regyr = _EM_RG_REGYR
    # subcycle count/lengths (regent.f:186-203, reuse KT logic)
    ntyr = Int(round(fint)); iyr = Int(regyr)
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1); nper < 1 && (nper = 1)
    kper = zeros(Int, nper); itot = ntyr; nn = nper
    @inbounds for i in 1:nper
        if nn == 1; kper[i] = itot; break; end
        kper[i] = itot ÷ nn; itot -= kper[i]; nn -= 1
    end
    # per-subcycle density from the large trees (regent.f:236-256, KT form)
    banext = fill(ba, nper); rdnext = fill(relden, nper)
    if nper > 1
        @inbounds for i in 1:n
            d1 = t.dbh[i]; d1 < 3.0f0 && continue
            sp = Int(t.species[i]); pr = t.tpa[i]
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d1)
            d2 = d1 + t.diam_growth[i] / bark
            b1 = 0.005454154f0 * d1 * d1; b2 = 0.005454154f0 * d2 * d2
            cc1 = em_tree_ccf(sp, d1); cc2 = em_tree_ccf(sp, d2)
            bi = (b2 - b1) / 10.0f0; ci = (cc2 - cc1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]; pn = pr * 0.985f0^k
                rdnext[j] += k * ci / pr * pn; banext[j] += k * bi * pn
            end
        end
    end
    delmax = (avh / 36.0f0) * (0.01232f0 * relden - 1.75f0); delmax > 0.0f0 && (delmax = 0.0f0)
    # per-tree accumulators over subcycles
    wk3 = Float32[t.height[i] for i in 1:n]
    dnow = Float32[t.dbh[i] for i in 1:n]
    @inbounds for j in 1:nper
        baj = banext[j]; rdj = rdnext[j]; kpj = Float32(kper[j])
        for i in 1:n
            sp = Int(t.species[i]); d = t.dbh[i]
            (d >= _em_rg_cap(sp) || t.tpa[i] <= 0.0f0) && continue
            _em_rg_nivar(sp) || continue                          # EMVAR conifers + LL(5)=NIVAR
            con = rhcon[sp] + c.htg_cor_small[sp]                # + HCOR (0 until calibrated)
            h1 = wk3[i]; h1 <= 0f0 && (h1 = 0.1f0)
            pct = Float32(t.crown_ratio[i]); bal = baj * (100.0f0 - pct) * 0.01f0
            htgrl = con + _EM_RG_BH * log(h1) + _EM_RG_BCCF * rdj + _EM_RG_BBAL * bal
            h2 = h1 + exp(htgrl) * (kpj / regyr)                 # XRHGRO=1 (no MULTS)
            wk3[i] = h2
            # per-subcycle DG dub (D<3, H2>4.5)
            if !(j >= nper || d >= 3.0f0) && h2 > 4.5f0
                relh = (h1 - 4.5f0) / (ah - 4.5f0); relh > 1f0 && (relh = 1f0); relh < 0f0 && (relh = 0f0)
                dadj = delmax * relh * relh - 2.0f0 * delmax * relh + 0.65f0
                d1v = h1 > 4.5f0 ? _EM_RG_AX * (h1 - 4.5f0)^_EM_RG_BX + dadj : EM_RG_DIAM[sp] + dadj
                d2v = _EM_RG_AX * (h2 - 4.5f0)^_EM_RG_BX + dadj
                dgj = (d2v - d1v); dgj < 0f0 && (dgj = 0f0)
                dnow[i] += dgj
            end
        end
    end
    # final HTGR1 + ZZRAN + XWT blend + DG override (regent.f:473-560, KT form)
    _sp_order = sortperm(view(t.species, 1:n); alg = Base.Sort.MergeSort)
    @inbounds for oi in 1:n
        i = _sp_order[oi]
        sp = Int(t.species[i]); d = t.dbh[i]
        (d >= _em_rg_cap(sp) || t.tpa[i] <= 0.0f0) && continue
        _em_rg_nivar(sp) || continue
        h = t.height[i]; xmn = EM_RG_XMIN[sp]; xmx = _em_rg_cap(sp)
        htgr1 = wk3[i] - h; htgr1 < 0.0f0 && (htgr1 = 0.0f0)
        zzran = 0.0f0
        if dgsd >= 1.0f0
            while true
                zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                (zzran <= 2.0f0 && zzran >= -2.0f0) && break
            end
        end
        htgr = htgr1 + zzran * _EM_RG_HSIGMA; htgr < 0.15f0 && (htgr = 0.15f0)
        xwt = d <= xmn ? 0.0f0 : (d - xmn) / (xmx - xmn)
        htg = htgr * (1.0f0 - xwt) + xwt * t.ht_growth[i]
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        t.ht_growth[i] = htg
        # DG override from the subcycle-accumulated diameter
        dgnew = dnow[i] - d; dgnew < 0f0 && (dgnew = 0f0)
        dgw = dgnew * (1.0f0 - xwt) + xwt * t.diam_growth[i]
        t.diam_growth[i] = dgw
        _em_rg_stash!(stash, t, i)
    end
    slo = s.coef.species[:site_lo]; shi = s.coef.species[:site_hi]; fint10 = fint/10.0f0
    pctred = _em_rg_pctred(ah, relden)     # em/regent.f:322-333 stand density modifier (AH·RELDEN/100), computed once
    @inbounds for i in 1:n
        sp = Int(t.species[i]); d = t.dbh[i]
        _em_rg_crvar(sp) || continue
        (d >= EM_RG_XMAX[sp] || t.tpa[i] <= 0.0f0) && continue
        h = t.height[i]; sitear = p.sp_site_index[sp]
        si = sitear; si > shi[sp] && (si = shi[sp]); si <= slo[sp] && (si = slo[sp]+0.5f0)
        relsi = (si-slo[sp])/(shi[sp]-slo[sp])
        xcr = Float32(t.crown_pct[i])/100.0f0
        vigor = 150.0f0*xcr^3*exp(-6.0f0*xcr)+0.3f0; vigor>1.0f0 && (vigor=1.0f0)
        pothtg = sitear/(15.0f0-4.0f0*relsi); con = exp(c.htg_cor_small[sp])
        htgr = pothtg*pctred*vigor*con*fint10
        if dgsd >= 1.0f0                                          # CRVAR ZZRAN — ONE draw, add only if in
            zz = bachlo(s.rng, 0.0f0, 1.0f0)                      # [-2,0.5] else skip (em/regent.f:811-813,919)
            (zz <= 0.5f0 && zz >= -2.0f0) && (htgr = htgr + zz * 0.2f0)
        end
        htgr < 0.1f0 && (htgr = 0.1f0)
        xmn = EM_RG_XMIN[sp]; xmx = EM_RG_XMAX[sp]
        xwt = d <= xmn ? 0.0f0 : (d-xmn)/(xmx-xmn)
        largeh = t.ht_growth[i]
        htg = htgr*(1.0f0-xwt)+xwt*largeh; htg < 0.1f0 && (htg=0.1f0)
        t.ht_growth[i] = htg
        # CRVAR small-tree DIAMETER (em/regent.f DO-25, D<1 only — D≥1 keeps large-tree DG, line 855 GOTO 23).
        # hk<4.5: DG=0 (line 884); hk≥4.5: inverse-Wykoff height→diameter DK model. Was previously a near-zero
        # placeholder (d+0.0001·h2) that froze cottonwood/OH seedlings once they grew past breast height.
        if d < 1.0f0
            hk = h + htg
            if hk < 4.5f0
                t.diam_growth[i] = 0.0f0
            else
                bk = bark_ratio(c.bark_a, c.bark_b, sp, d)
                axc = c.ht_dbh_iabflg[sp] == 1 ? EM_RG_HT1[sp] : c.ht_dbh_aa[sp]   # em/regent.f:924-927
                t.diam_growth[i] = _em_crut_dg(sp, d, h, hk, htg, sitear, bk, EM_RG_DGMAX[sp]*fint10, 10f0/fint, axc)
            end
        end
        _em_rg_stash!(stash, t, i)
    end
    # UTVAR (RM6 juniper / AS12,PB17 aspen) — SINGLE-STEP (em/regent.f:507-533,600). CON=1.0 (non-NIVAR).
    @inbounds for i in 1:n
        sp = Int(t.species[i]); d = t.dbh[i]
        _em_rg_utvar(sp) || continue
        (d >= EM_RG_XMAX[sp] || t.tpa[i] <= 0.0f0) && continue
        h = t.height[i]; sitear = p.sp_site_index[sp]
        si = sitear; si > shi[sp] && (si = shi[sp]); si <= slo[sp] && (si = slo[sp]+0.5f0)
        relsi = (si-slo[sp])/(shi[sp]-slo[sp]); rsimod = 0.5f0*(1.0f0+relsi)
        con = exp(c.htg_cor_small[sp])                            # RHCON(=1.0 non-NIVAR)·exp(HCOR)
        if sp == 12 || sp == 17                                   # aspen/PB Sheppard (ABIRTH)
            ab = Float32(t.birth_age[i]); ab < 1.0f0 && (ab = 1.0f0)
            hite1 = 26.9825f0 * ab^1.1752f0
            hite2 = 26.9825f0 * (ab + 10.0f0)^1.1752f0
            htgr = (hite2 - hite1) / (2.54f0 * 12.0f0) * rsimod * con * 0.75f0
        else                                                      # RM juniper — density PCTRED (AH·RELDEN/100)
            xcr = Float32(t.crown_pct[i])/100.0f0
            vigor = 150.0f0*xcr^3*exp(-6.0f0*xcr)+0.3f0; vigor>1.0f0 && (vigor=1.0f0)
            vigor = 1.0f0 - (1.0f0 - vigor)/3.0f0                 # RM ⅔ VIGOR cut (em/regent.f:515)
            pothtg = (si/10.0f0) * (si*1.5f0 - h) / (si*1.5f0)
            htgr = pothtg*pctred*vigor*con
        end
        htgr = htgr * fint10
        if dgsd >= 1.0f0                                          # UTVAR ZZRAN — ·0.1 (not 0.2), bounds [-2,0.5]
            zz = bachlo(s.rng, 0.0f0, 1.0f0)
            (zz <= 0.5f0 && zz >= -2.0f0) && (htgr = htgr + zz * 0.1f0)
        end
        htgr < 0.1f0 && (htgr = 0.1f0)
        xmn = EM_RG_XMIN[sp]; xmx = EM_RG_XMAX[sp]
        xwt = d <= xmn ? 0.0f0 : (d-xmn)/(xmx-xmn)
        htg = htgr*(1.0f0-xwt)+xwt*t.ht_growth[i]; htg < 0.1f0 && (htg=0.1f0)
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        t.ht_growth[i] = htg
        # UTVAR small-tree DIAMETER (em/regent.f DO-25). Juniper(6) always DK (SITEAR linear, XMAX=99); aspen/PB
        # (12,17) D<3 → inverse-Wykoff DK once hk≥4.5 (D≥3 keeps large DG, line 606/857 GOTO 23). hk<4.5 → DG=0
        # (line 881-884). Was previously frozen at the ~0 large-tree DG for hk>4.5 ⇒ aspen seedlings never grew
        # diameter (measured: 0.1" aspen got DG≈0 in jl vs DG≈1.4" in FVSem_g16).
        hk = h + htg
        if sp == 6 || d < 3.0f0
            if hk < 4.5f0
                t.diam_growth[i] = 0.0f0
            else
                bk = bark_ratio(c.bark_a, c.bark_b, sp, d)
                axc = c.ht_dbh_iabflg[sp] == 1 ? EM_RG_HT1[sp] : c.ht_dbh_aa[sp]   # em/regent.f:924-927
                t.diam_growth[i] = _em_crut_dg(sp, d, h, hk, htg, sitear, bk, EM_RG_DGMAX[sp], 10f0/fint, axc)
            end
        end
        _em_rg_stash!(stash, t, i)
    end
    # TTVAR (LM4) — subcycle BETA height + DLESS3 DG (em/regent.f:471-957). CON=1.0. ZRAND redraw/cycle [-2,2].
    lm_present = false
    @inbounds for i in 1:n; (Int(t.species[i]) == 4 && t.tpa[i] > 0f0) && (lm_present = true; break); end
    if lm_present
        wk3lm = Float32[t.height[i] for i in 1:n]; dklm = Float32[t.dbh[i] for i in 1:n]
        zrlm = zeros(Float32, n)
        @inbounds for i in 1:n
            (Int(t.species[i]) == 4 && dgsd >= 1.0f0) || continue
            z = 0f0; while true; z = bachlo(s.rng, 0f0, 1f0); (-2f0 <= z <= 2f0) && break; end
            zrlm[i] = z
        end
        conlm = exp(c.htg_cor_small[4])
        @inbounds for j in 1:nper
            kpj = Float32(kper[j])
            for i in 1:n
                Int(t.species[i]) == 4 || continue
                d = t.dbh[i]; (d >= EM_RG_XMAX[4] || t.tpa[i] <= 0f0) && continue
                h1 = wk3lm[i]; cr = Float32(t.crown_pct[i])
                pt = Int(t.plot_id[i]); pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0f0
                tpccf = pccf; tpccf > 300f0 && (tpccf = 300f0); tpccf < 25f0 && (tpccf = 25f0)
                beta1 = exp(1.17527f0 - 0.42124f0*log(tpccf)); beta2 = exp(-2.56002f0 - 0.58642f0*log(tpccf))
                htg1 = beta1 + beta2*cr; stddev = htg1*(1.08720f0 - 0.00230f0*cr)
                htgrl = htg1 + zrlm[i]*stddev; htgrl < 0.1f0 && (htgrl = 0.1f0)
                h2 = h1 + htgrl*(kpj/regyr)*conlm; wk3lm[i] = h2
                dklm[i] = _em_dless3(h2, cr, pccf)
            end
        end
        @inbounds for i in 1:n
            Int(t.species[i]) == 4 || continue
            d = t.dbh[i]; (d >= EM_RG_XMAX[4] || t.tpa[i] <= 0f0) && continue
            h = t.height[i]; cr = Float32(t.crown_pct[i]); xmn = EM_RG_XMIN[4]; xmx = EM_RG_XMAX[4]
            xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)
            htgr = wk3lm[i] - h; htgr < 0f0 && (htgr = 0f0)
            htg = htgr*(1f0-xwt) + xwt*t.ht_growth[i]
            cap = s.control.sp_size_cap[4,4]; (h+htg > cap) && (htg = max(cap-h, 0.1f0))
            t.ht_growth[i] = htg
            pt = Int(t.plot_id[i]); pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0f0
            dkk = _em_dless3(h, cr, pccf)
            bark = bark_ratio(c.bark_a, c.bark_b, 4, d)
            dgr = (dklm[i] - dkk)*bark
            dds = dgr*(2f0*bark*d + dgr)
            arg = (d*bark)^2 + dds
            dgk = arg > 0f0 ? sqrt(arg) - bark*d : 0f0
            dgk < 0f0 && (dgk = 0f0); dgk > fint*2.0f0 && (dgk = fint*2.0f0)
            t.diam_growth[i] = dgk*(1f0-xwt) + xwt*t.diam_growth[i]
            _em_rg_stash!(stash, t, i)
        end
    end
    # EMVAR conifers {1,2,3,7,8,9,10,18} — SMHTGF (height→wk3e) + SMDGF (diameter→wk5e). em/regent.f EMVAR
    # branch (H2=H1+HTGRTH·(KPER/REGYR)·XRHGRO·CON, CON=RHCON·exp(HCOR)=exp(HCOR) since RHCON=1, regent.f:1484).
    # ZRAND: ONE BACHLO per tree drawn in species order (matches FVS DO 30 ISPC loop). SMHTGF uses CLAMPED TPCCF
    # [25,300]; SMDGF uses RAW point CCF (regent.f:448 vs SMDGF call). CR = crown PERCENT (=FLOAT(ICR), regent.f:451).
    emvar_present = false
    @inbounds for i in 1:n; (_em_orig_species(Int(t.species[i])) && t.tpa[i] > 0f0) && (emvar_present = true; break); end
    if emvar_present
        wk3e = Float32[t.height[i] for i in 1:n]; wk5e = Float32[t.dbh[i] for i in 1:n]
        zre = zeros(Float32, n)
        _emv_order = sortperm(view(t.species, 1:n); alg = Base.Sort.MergeSort)   # species order (stable within sp)
        @inbounds for oi in 1:n
            i = _emv_order[oi]; sp = Int(t.species[i])
            # em/smhtgf.f draws ZRAND UNCONDITIONALLY (no DGSD gate — only dgdriv's OLDRN seeding and DGSCOR
            # are gated). Gating it here made a DGSD<1 run consume fewer draws than live (MEASURED with
            # `DGSTDEV 0.` on the EM RD stand: live's first REGENT ZRAND lands at main-stream rann call 27,
            # jl's at 24 — one BACHLO (3 calls) short, exactly this gate).
            (_em_orig_species(sp) && t.dbh[i] < _em_rg_cap(sp) && t.tpa[i] > 0f0) || continue
            z = 0f0; while true; z = bachlo(s.rng, 0f0, 1f0); (-2f0 <= z <= 2f0) && break; end
            zre[i] = z
        end
        @inbounds for j in 1:nper
            kpj = Float32(kper[j])
            for i in 1:n
                sp = Int(t.species[i]); d = t.dbh[i]
                (_em_orig_species(sp) && d < _em_rg_cap(sp) && t.tpa[i] > 0f0) || continue
                h1 = wk3e[i]; cr = Float32(t.crown_pct[i])
                pt = Int(t.plot_id[i]); pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0f0
                tpccf = pccf; tpccf > 300f0 && (tpccf = 300f0); tpccf < 25f0 && (tpccf = 25f0)
                htgrth = _em_smhtgf(sp, cr, tpccf, zre[i])
                con = exp(c.htg_cor_small[sp])                       # RHCON(=1)·exp(HCOR)
                h2 = h1 + htgrth*(kpj/regyr)*con; wk3e[i] = h2       # XRHGRO=1
                # SMDGF: DBH from grown height (raw PCCF) — ONLY for h2>4.5 (em/regent.f:578 IF(H2.LE.4.5)GO TO 14
                # skips the DBH; the HLESS4-form goes NEGATIVE below 4.5). wk5e stays d until the tree crosses 4.5.
                h2 > 4.5f0 && (wk5e[i] = _em_smdgf(sp, h2, cr, pccf))
            end
        end
        @inbounds for i in 1:n
            sp = Int(t.species[i]); d = t.dbh[i]
            (_em_orig_species(sp) && d < _em_rg_cap(sp) && t.tpa[i] > 0f0) || continue
            h = t.height[i]; cr = Float32(t.crown_pct[i]); xmn = EM_RG_XMIN[sp]; xmx = _em_rg_cap(sp)
            xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)
            htgr = wk3e[i] - h; htgr < 0f0 && (htgr = 0f0)
            htg = htgr*(1f0-xwt) + xwt*t.ht_growth[i]
            cap = s.control.sp_size_cap[sp,4]; (h+htg > cap) && (htg = max(cap-h, 0.1f0))
            t.ht_growth[i] = htg
            pt = Int(t.plot_id[i]); pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0f0
            dkk = h > 4.5f0 ? _em_smdgf(sp, h, cr, pccf) : d          # SMDGF at start height (=d below 4.5, no neg-SMDGF)
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
            dgr = (wk5e[i] - dkk)*bark; dds = dgr*(2f0*bark*d + dgr)
            arg = (d*bark)^2 + dds; dgk = arg > 0f0 ? sqrt(arg) - bark*d : 0f0
            dgk < 0f0 && (dgk = 0f0); dgk > fint*2.0f0 && (dgk = fint*2.0f0)
            t.diam_growth[i] = dgk       # NO XWT blend on DG — FVS blends only HTG (em/regent.f:838); DG is the
                                         # pure SMDGF regent DG (blending the large-tree DG over-grew BA ~20% cyc2+).
            _em_rg_stash!(stash, t, i)
        end
    end
    return s
end

# em/esgent.f (estb/esgent.f:53 CALL REGENT(.TRUE.,ITRNIN)) — grow the JUST-ESTABLISHED regen IN its birth cycle.
# EM was OMITTED from the birth-cycle esgent list (simulate.jl had only CR + TT), so planted/established EM regen
# missed its first-cycle height growth ⇒ the cohort stayed at the ~1' establishment height at the birth-cycle report
# (TopHt 1 vs live 6) ⇒ ~1-cycle height/BA lag ⇒ dense self-thin under-kill (#137). Mirrors the EMVAR branch of
# small_tree_growth! (SMHTGF height + SMDGF diameter) for the new records nstart+1:n over the birth-cycle subperiod
# (GENTIM germination offset, like tt_esgent!). Only EMVAR conifers {1,2,3,7,8,9,10,18}; others left as-established.
function em_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0)
    t = s.trees; c = s.calib; dens = s.density
    nstart >= t.n && return s
    rhcon = em_regcons!(s)                              # RHCON (also refreshes birth-cycle density constants)
    dgsd = s.control.dg_sd
    @inbounds for i in (nstart+1):t.n
        t.tpa[i] <= 0.0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]
        (_em_orig_species(sp) && d < _em_rg_cap(sp)) || continue
        h = t.height[i]; cr = Float32(t.crown_pct[i])
        pt = Int(t.plot_id[i]); pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0.0f0
        tpccf = pccf; tpccf > 300.0f0 && (tpccf = 300.0f0); tpccf < 25.0f0 && (tpccf = 25.0f0)
        zrand = 0.0f0
        if dgsd >= 1.0f0
            while true; zrand = bachlo(s.rng, 0.0f0, 1.0f0); (-2.0f0 <= zrand <= 2.0f0) && break; end
        end
        htgrth = _em_smhtgf(sp, cr, tpccf, zrand)
        con = exp(c.htg_cor_small[sp])                  # RHCON(=1)·exp(HCOR)
        # HTG × per-tree WK4=HTIMLT birth-cycle multiplier (live em/esgent.f:23 HTG=HTG*WK4). PLANT/existing=1.0;
        # AUTOES natural regen=0.40 (#193 — jl formerly used a constant subyr/regyr=1 ⇒ 2.5× seedling over-growth).
        htg = htgrth * t.htimlt[i] * con; htg < 0.0f0 && (htg = 0.0f0)
        cap = s.control.sp_size_cap[sp, 4]; (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        h2 = h + htg
        t.height[i] = h2; t.ht_growth[i] = htg
        if h2 > 4.5f0                                    # DBH only once the tree crosses breast height
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
            d2 = _em_smdgf(sp, h2, cr, pccf)
            dkk = h > 4.5f0 ? _em_smdgf(sp, h, cr, pccf) : d
            dgr = (d2 - dkk) * bark; dds = dgr * (2.0f0 * bark * d + dgr)
            arg = (d * bark)^2 + dds
            dgk = arg > 0.0f0 ? sqrt(arg) - bark * d : 0.0f0
            dgk < 0.0f0 && (dgk = 0.0f0)
            dgk > 0.0f0 && (t.dbh[i] = d + dgk / bark; t.diam_growth[i] = dgk)
        end
    end
    return s
end
