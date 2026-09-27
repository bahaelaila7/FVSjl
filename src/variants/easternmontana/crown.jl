# =============================================================================
# crown.jl (easternmontana) — per-tree CCF (em/ccfcal.f MODE=1).
#
# EM CCF is the direct per-species polynomial (Paine-Hann / NI INT-133 form), summed
# over the tree list = stand CCF (like KT/IE ccfcal, not the CR crown-width→area path).
#   large:  CCFT = RD1 + D·RD2 + D²·RD3
#   small:  CCFT = RDA · D^RDB
# The large/small DBH threshold is species-group dependent (em/ccfcal.f SELECT CASE):
#   sp 5 (LL):                 D≥10 → poly, else RDA·D^RDB           (no 0.001 floor)
#   sp 1-4,6-10,12,17,18:      D≥1  → poly; D>0.1 → RDA·D^RDB; else 0.001
#   sp 11,13-16,19:            D≥10 → poly; D>0.1 → RDA·D^RDB; else 0.001
# Coefficients transcribed verbatim from em/ccfcal.f DATA (species 1=WB..19=OH).
# =============================================================================

const EM_RD1 = Float32[0.0186, 0.0392, 0.0388, 0.01925, 0.03, 0.01925, 0.01925, 0.03, 0.0172,
                       0.0219, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.0204, 0.03]
const EM_RD2 = Float32[0.0146, 0.0180, 0.0269, 0.01676, 0.0216, 0.01676, 0.01676, 0.0173, 0.00876,
                       0.0169, 0.0215, 0.0238, 0.0215, 0.0215, 0.0215, 0.0215, 0.0238, 0.0246, 0.0215]
const EM_RD3 = Float32[0.00288, 0.00207, 0.00466, 0.00365, 0.00405, 0.00365, 0.00365, 0.00259, 0.00112,
                       0.00325, 0.00363, 0.00490, 0.00363, 0.00363, 0.00363, 0.00363, 0.00490, 0.0074, 0.00363]
const EM_RDA = Float32[0.009884, 0.007244, 0.017299, 0.009187, 0.011402, 0.009187, 0.009187, 0.007875,
                       0.011402, 0.007813, 0.011109, 0.008915, 0.011109, 0.011109, 0.011109, 0.011109,
                       0.008915, 0.011109, 0.011109]
const EM_RDB = Float32[1.6667, 1.8182, 1.5571, 1.7600, 1.7560, 1.7600, 1.7600, 1.7360, 1.7560,
                       1.7780, 1.7250, 1.7800, 1.7250, 1.7250, 1.7250, 1.7250, 1.7800, 1.7250, 1.7250]

# em/ccfcal.f MODE=1: per-tree CCF (before the ×P expansion the caller applies).
@inline function em_tree_ccf(sp::Integer, d::Real)::Float32
    # D=0 (the IMC=9 dead the backdating DENSE zeroes, dense.f:75) takes the D≤0.1 branch ⇒ 0.001 (LL: RDA·0**RDB=0)
    d <= 0f0 && return sp == 5 ? 0f0 : 0.001f0
    poly()  = EM_RD1[sp] + d * EM_RD2[sp] + d * d * EM_RD3[sp]
    small() = EM_RDA[sp] * fpow(Float32(d), EM_RDB[sp])       # ccfcal.f RDA*(D**RDB): gfortran powf
    if sp == 5                                   # LL: D≥10 poly, else small (no floor)
        return d >= 10f0 ? poly() : small()
    elseif sp == 11 || (13 <= sp <= 16) || sp == 19   # GA/CW/BA/PW/NC/OH: D≥10 threshold
        return d >= 10f0 ? poly() : (d > 0.1f0 ? small() : 0.001f0)
    else                                          # 1-4,6-10,12,17,18: D≥1 threshold
        return d >= 1f0 ? poly() : (d > 0.1f0 ? small() : 0.001f0)
    end
end

# =============================================================================
# EM crown-ratio update (em/crown.f) — the IE/NI Weibull/DCR change-in-crown model.
# Single PARM(14) vector (IE/NI sp9 form, applied to all species) + habitat intercept
# CRHAB[MAPHAB[ITYPE]] + CRSD=6.35. Mirrors KT's crown_ratio_update! DCR structure
# (crown change = exp(PCR) − exp(DCR-backdated), bounded ±1%/yr; d<3 keeps crown while cycling).
# =============================================================================
# em/crown.f PARM(14): 1-6 density(BA,BA²,lnBA,RELDEN,RELDEN²,lnRELDEN), 7-14(D,D²,lnD,H,H²,lnH,P,lnP)
const EM_CRPARM = Float32[-0.00190, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.23372, 0.0, 0.0, -0.28433, 0.001903, 0.0]
const EM_CRHAB  = Float32[0.09453, -0.07740, 0.07113, 0.2039, 0.06176, 0.1513, 0.09086, 0.1580, 0.09229, 0.01551, 0.0, 0.0, 0.0, 0.0]
const EM_CR_MAPHAB = Int[2,2,2,2,2,2,2,2,2,2,2,2,2,3,4,4,4,4,5,6,6,7,6,1,8,1,9,10,10,6]  # ITYPE→CRHAB idx (13*2,3,4*4,5,6,6,7,6,1,8,1,9,2*10,6)
const _EM_CRSD = 6.35f0

# em/blkdat.f DATA HT1/HT2 — the Wykoff HT-DBH intercept/slope cratet.f dubs missing heights with (and fits AA
# against). EM's species_coefficients.csv ht1/ht2/wykoff_ht2 columns are NOT these (CR placeholders), so the dub
# and the AA fit must read these. Verified vs live FVSem_g16 CRATET debug: sp2 default INTERCEPT 4.1539 / SLOPE −4.212.
const EM_BLK_HT1 = Float32[4.1539, 4.1539, 4.4161, 4.1920, 4.76537, 3.2, 4.5356, 4.7537, 4.5788, 4.414,
                           4.4421, 4.4421, 4.4421, 4.4421, 4.4421, 4.4421, 4.4421, 4.1539, 4.4421]
const EM_BLK_HT2 = Float32[-4.212, -4.212, -6.962, -5.1651, -7.61062, -5.0, -5.692, -8.356, -7.138, -8.907,
                           -6.5405, -6.5405, -6.5405, -6.5405, -6.5405, -6.5405, -6.5405, -4.212, -6.5405]

# em/crown.f DATA — the EMVAR/UTTVAR rank-Weibull crown coefficients (generated from the Fortran).
const EM_CR_WEIBA = Float32[0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
const EM_CR_WEIBB0 = Float32[0.11035, 0.11035, 0.14652, -0.82631, 0.0, 0.0, -0.00359, 0.67059, 0.73693, 0.02663, 0.0, -0.08414, 0.0, 0.0, 0.0, 0.0, -0.08414, 0.11035, 0.0]
const EM_CR_WEIBB1 = Float32[1.10085, 1.10085, 1.09052, 1.06217, 0.0, 0.0, 1.12728, 0.99349, 0.98414, 1.11477, 0.0, 1.14765, 0.0, 0.0, 0.0, 0.0, 1.14765, 1.10085, 0.0]
const EM_CR_WEIBC0 = Float32[0.02774, 0.02774, 1.04746, 3.31429, 0.0, 0.0, 2.60377, -4.25938, -4.16681, 2.95048, 0.0, 2.77500, 0.0, 0.0, 0.0, 0.0, 2.77500, 0.02774, 0.0]
const EM_CR_WEIBC1 = Float32[0.35524, 0.35524, 0.39752, 0.0, 0.0, 0.0, 0.0, 1.35687, 1.33779, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.35524, 0.0]
const EM_CR_C0 = Float32[5.68625, 5.68625, 5.92714, 6.19911, 0.0, 0.0, 5.05870, 7.41093, 7.36476, 5.61047, 0.0, 4.01678, 0.0, 0.0, 0.0, 0.0, 4.01678, 5.68625, 0.0]
const EM_CR_C1 = Float32[-0.04470, -0.04470, -0.03346, -0.02216, 0.0, 0.0, -0.03307, -0.03467, -0.03761, -0.03557, 0.0, -0.01516, 0.0, 0.0, 0.0, 0.0, -0.01516, -0.04470, 0.0]

# em/dubscr.f DATA (19 species) — generated from the Fortran source.
const EM_DUB_BCR0 = Float32[-1.669490, -1.669490, -0.426688, -1.66949, -0.89014, 0.0, -1.669490, -0.426688, -0.426688, -1.669490, 0.0, -0.426688, 0.0, 0.0, 0.0, 0.0, -0.426688, -2.19723, 0.0]
const EM_DUB_BCR1 = Float32[-0.209765, -0.209765, -0.093105, -0.209765, -0.18026, 0.0, -0.209765, -0.093105, -0.093105, -0.209765, 0.0, -0.093105, 0.0, 0.0, 0.0, 0.0, -0.093105, 0.0, 0.0]
const EM_DUB_BCR2 = Float32[0.0, 0.0, 0.022409, 0.0, 0.02233, 0.0, 0.0, 0.022409, 0.022409, 0.0, 0.0, 0.022409, 0.0, 0.0, 0.0, 0.0, 0.022409, 0.0, 0.0]
const EM_DUB_BCR3 = Float32[0.003359, 0.003359, 0.002633, 0.003359, 0.00614, 0.0, 0.003359, 0.002633, 0.002633, 0.003359, 0.0, 0.002633, 0.0, 0.0, 0.0, 0.0, 0.002633, 0.0, 0.0]
const EM_DUB_BCR5 = Float32[0.011032, 0.011032, 0.0, 0.011032, 0.0, 0.0, 0.011032, 0.0, 0.0, 0.011032, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
const EM_DUB_BCR6 = Float32[0.0, 0.0, -0.045532, 0.0, 0.0, 0.0, 0.0, -0.045532, -0.045532, 0.0, 0.0, -0.045532, 0.0, 0.0, 0.0, 0.0, -0.045532, 0.0, 0.0]
const EM_DUB_BCR8 = Float32[0.017727, 0.017727, 0.0, 0.017727, 0.0, 0.0, 0.017727, 0.0, 0.0, 0.017727, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
const EM_DUB_BCR9 = Float32[-0.000053, -0.000053, 0.000022, -0.000053, 0.0, 0.0, -0.000053, 0.000022, 0.000022, -0.000053, 0.0, 0.000022, 0.0, 0.0, 0.0, 0.0, 0.000022, 0.0, 0.0]
const EM_DUB_BCR10 = Float32[0.014098, 0.014098, -0.013115, 0.014098, 0.0, 0.0, 0.014098, -0.013115, -0.013115, 0.014098, 0.0, -0.013115, 0.0, 0.0, 0.0, 0.0, -0.013115, 0.0, 0.0]
const EM_DUB_CRSD = Float32[0.5000, 0.5000, 0.6957, 0.5000, 0.8871, 0.0, 0.6124, 0.6957, 0.6957, 0.4942, 0.0, 0.9310, 0.0, 0.0, 0.0, 0.0, 0.9310, 0.200, 0.0]

"""em/dubscr.f — crown for a record with a missing crown (the EMVAR/UTTVAR d<3 live path and the cycle-0 dead
records). One rejection-bounded BACHLO(0,CRSD) draw when DGSD>=1, then the logistic transform. RMAI is the
MAICAL site MAI (em/maical.f ADJMAI(ISPNUM(ISISP), SITEAR(ISISP), 10), capped 128; grinit default 50)."""
@inline function em_dubscr(rng, sp::Integer, d::Real, h::Real, ba::Real, tpccf::Real, avh::Real,
                           rmai::Real, dgsd::Real)::Float32
    hf = Float32(h)
    cr = EM_DUB_BCR0[sp] + EM_DUB_BCR1[sp]*Float32(d) + EM_DUB_BCR2[sp]*hf + EM_DUB_BCR3[sp]*Float32(ba) +
         EM_DUB_BCR5[sp]*Float32(tpccf) + EM_DUB_BCR6[sp]*(Float32(avh)/hf) + EM_DUB_BCR8[sp]*Float32(avh) +
         EM_DUB_BCR9[sp]*(Float32(ba)*Float32(tpccf)) + EM_DUB_BCR10[sp]*Float32(rmai)
    sd = EM_DUB_CRSD[sp]; fcr = 0f0
    while true
        fcr = dgsd >= 1f0 ? bachlo(rng, 0f0, sd) : 0f0
        abs(fcr) > sd && continue
        break
    end
    abs(cr + fcr) >= 86f0 && (cr = 86f0)
    cr = 1f0 / (1f0 + exp(cr + fcr))
    cr > 0.95f0 && (cr = 0.95f0); cr < 0.05f0 && (cr = 0.05f0)
    return cr
end

# em/maical.f — RMAI = ADJMAI(ISPNUM(ISISP), SITEAR(ISISP)|140 if 0, 10), capped at 128. ISISP defaults to 3.
const EM_MAI_ISPNUM = Int[119, 117, 202, 101, 19, 101, 108, 93, 21, 122, 101, 101, 101, 101, 101, 101, 101, 122, 101]
function _em_rmai(s::StandState)::Float32
    isisp = Int(s.plot.site_species); (isisp < 1 || isisp > 19) && (isisp = 3)
    sssi = s.plot.sp_site_index[isisp]; sssi == 0f0 && (sssi = 140f0)
    r = _adjmai(EM_MAI_ISPNUM[isisp], sssi, 10f0)
    return r > 128f0 ? 128f0 : r
end

function crown_ratio_update!(s::StandState, ::EasternMontana; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t = s.plot, s.trees
    t.n == 0 && return s
    itype = Int(p.habitat_input); it = (1 <= itype <= 30) ? itype : 1
    crcon = EM_CRHAB[clamp(EM_CR_MAPHAB[it], 1, 14)]
    ba = p.basal_area; relden = p.relative_density
    lnba = ba > 0f0 ? log(ba) : 0f0; lnrd = relden > 0f0 ? log(relden) : 0f0
    reldm1 = p.relative_density_prev; oba = p.old_ba
    if reldm1 < 100f0; oba = ba; reldm1 = relden; end
    x1 = (!lstart && oba > 0f0) ? log(oba) : 0f0
    x2 = (!lstart && reldm1 > 0f0) ? log(reldm1) : 0f0
    dgsd = s.control.dg_sd
    ba_a = s.calib.bark_a; ba_b = s.calib.bark_b
    P = EM_CRPARM
    sdiac = crown_sdi                                  # SDIAC (pre-growth Reineke SDI) for RELSDI
    rmai_v = _em_rmai(s)                                # DUBSCR RMAI (em/maical.f)
    # crown.f ISORT: whole-stand descending-DBH rank on the CURRENT DBH (grown at cycling, as-read
    # at LSTART) — shared crown_isort, see crown_init.jl.
    isort = crown_isort(s; lstart = lstart)
    # #158-class species-major RNG order: FVS em/crown.f processes trees SPECIES-MAJOR
    # (`DO 70 ISPC=1,MAXSP; DO 60 I3=I1,I2; I=IND1(I3)`), so the per-tree DUBSCR/NIVAR BACHLO crown draw
    # (line ~90, rejection-sampled with a species-specific SD) is consumed in species order. jl dubbed in raw
    # tree-index order ⇒ on a multi-species seedling cohort the per-tree crown draws were mis-assigned ⇒ wrong
    # crown ⇒ SMHTGF over-/under-grows that seedling's height ⇒ compounding dense-stand divergence. Only the
    # lstart DUB path draws RNG, so iterate species-major there; cycling (no draw) keeps natural order. Within a
    # species, tree-index order = FVS IND1 (stable species bucket). The deterministic per-tree crown update is
    # order-independent, so this is a no-op except on the RNG stream. EM-only (dispatches on ::EasternMontana).
    cur_year = current_cycle_year(s)   # em/crown.f CRNMULT block overwrites CRNMLT/DLOW/DHI per species
    order = species_major_order(s)   # em/crown.f DO 70 ISPC … I=IND1(I3) — LIVE records only; the cycle-0
                                     # dead records get their own DO 79 pass below (they were folded in here).
    @inbounds for i in order
        t.tpa[i] <= 0f0 && continue
        icr = Int(t.crown_pct[i])
        (lstart && icr > 0) && continue
        icr < 0 && (t.crown_pct[i] = Int32(-icr); continue)
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        bark = variant_bratio(s, sp, d, h)
        local icri::Int
        # em/crown.f:322-334 — CRVAR (GA/CW/BA/PW/NC/OH: sp 11,13-16,19; the CR-variant hardwood expansion)
        # use a LINEAR crown-LENGTH model, NOT the NIVAR PARM PCR/DCR logistic below and NOT DUBSCR. For a
        # sub-1" seedling CL clamps to HF ⇒ CR=1.0 ⇒ ICR caps at 95 (measured vs FVSem_g16: seedling ICR=95,
        # VIGOR=150·0.95³·e^-5.7+0.3=0.731; the PARM path gave ICR≈36-62 ⇒ VIGOR clamps to 1.0 ⇒ regent HTG
        # over by 1.0/0.731≈1.37×). No RNG draw on this path (crown.f CRVAR never calls BACHLO). D<3 is NOT
        # frozen for CRVAR — the CL model updates every cycle (the freeze below is the EMVAR/UTTVAR path).
        crvar = sp == 11 || (13 <= sp <= 16) || sp == 19
        lpiju = sp == 6                                    # crown.f:318 LPIJU (RM) — crown-LENGTH form, like CRVAR
        nivar = sp == 5                                    # crown.f:314 NIVAR (LL) — the PARM PCR/DCR model
        cmult, cdlow, cdhi = crn_mult_band(s.control, sp, cur_year; lstart = lstart)
        # em/crown.f's band test is STRICT at the top — `.GE. DLOW .AND. .LT. DHI` (:432/458/467/519) — unlike
        # the westside variants' `.LE. DHI`.
        inband = cdlow <= d < cdhi
        if crvar || lpiju
            htg = t.ht_growth[i]
            hf = h + htg                                  # crown.f:325 HF=H+HTG (H is post-growth at CROWN)
            cl = crvar ? 5.17281f0 + 0.32552f0*hf - 0.01675f0*ba :   # crown.f:328 CRVAR crown-length
                         -0.59373f0 + 0.67703f0*hf                    # crown.f:330 LPIJU crown-length
            cl < 1.0f0 && (cl = 1.0f0)                     # crown.f:331
            cl > hf && (cl = hf)                           # crown.f:332
            crnew = (cl/hf)*100.0f0                        # crown.f:333-334 CR=CL/HF; CRNEW=CR*100
            if !lstart || icr > 0                          # crown.f:450 bounded-change path (±1%/yr)
                chg = crnew - Float32(icr)
                pdifpy = chg/Float32(icr)/fint             # crown.f:452 (no *100; crnew already in %)
                pdifpy > 0.01f0  && (chg = Float32(icr)*0.01f0*fint)
                pdifpy < -0.01f0 && (chg = Float32(icr)*(-0.01f0)*fint)
                icri = trunc(Int, (Float32(icr) + (inband ? chg * cmult : chg)) + 0.5f0)   # crown.f:458-463
                if !lstart && icr != 0                     # crown.f:470-486 crown-length max (cycling, icr>0)
                    crln  = h*Float32(icr)/100.0f0
                    crmax = (crln + htg)/(h + htg)*100.0f0
                    (cmult == 1f0 && icri < 10) && (icri = trunc(Int, crmax + 0.5f0))   # crown.f:485
                    Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
                end
            else                                           # crown.f:465-467 DUB (lstart, icr==0)
                icri = trunc(Int, crnew + 0.5f0)
                inband && (icri = trunc(Int, Float32(icri) * cmult))   # crown.f:467-468
            end
            lstart && (icri = topkill_icri(t, i, icri))   # crown.f stmt 55 (crown_init.jl)
            icri > 95 && (icri = 95)                        # crown.f:590
            (cmult == 1f0 && icri < 10) && (icri = 10)      # stmt 59, crown.f:530
            icri < 1  && (icri = 1)                         # crown.f:595
            t.crown_pct[i] = Int32(icri)
            continue
        end
        # crown.f:336-340 — the size gate differs by class: EMVAR/UTTVAR fall to stmt 58 (DUBSCR) only for a
        # sub-1" record at LSTART; NIVAR falls to 58 whenever D<3.
        small = (crvar || lpiju) ? false : (!nivar ? (d < 1f0 && lstart) : d < 3f0)
        if small
            # crown.f:506-516 stmt 58: NIVAR when cycling ⇒ keep the crown (GO TO 60); EMVAR/UTTVAR with a crown
            # already set ⇒ keep it; otherwise DUBSCR at this record's point, then the shared bounds tail.
            (!lstart && nivar) && continue
            (!nivar && icr != 0) && continue
            pt_i = Int(t.plot_id[i])
            tpccf = (1 <= pt_i <= length(s.density.point_ccf)) ? s.density.point_ccf[pt_i] : 0f0
            icri = icri_round(em_dubscr(s.rng, sp, d, h, ba, tpccf, p.avg_height, rmai_v, dgsd))
            inband && (icri = trunc(Int, Float32(icri) * cmult))   # crown.f:519-520
            icri > 95 && (icri = 95); (cmult == 1f0 && icri < 10) && (icri = 10); icri < 1 && (icri = 1)   # stmt 59
            t.crown_pct[i] = Int32(icri)
            continue
        end
        if !nivar
            # crown.f:341-356 EMVAR/UTTVAR — rank-Weibull crown on the RELDEN SCALE (NOT the NIVAR PARM model,
            # which jl applied to every non-CRVAR species). A/B/C from ACRNEW(RELSDI), X = ISORT/ITRN · SCALE.
            relsdi = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0
            relsdi > 1.5f0 && (relsdi = 1.5f0)
            acrnew = EM_CR_C0[sp] + EM_CR_C1[sp] * relsdi * 100f0
            A = EM_CR_WEIBA[sp]
            B = EM_CR_WEIBB0[sp] + EM_CR_WEIBB1[sp] * acrnew; B < 1f0 && (B = 1f0)
            C = EM_CR_WEIBC0[sp] + EM_CR_WEIBC1[sp] * acrnew; C < 2f0 && (C = 2f0)
            scale = 1f0 - 0.00167f0 * (relden - 100f0)
            scale > 1f0 && (scale = 1f0); scale < 0.30f0 && (scale = 0.30f0)
            x = d > 0f0 ? (Float32(isort[i]) / Float32(t.n)) * scale : rann!(s.rng) * scale
            x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
            crnew = (A + B * (-log(1f0 - x))^(1f0 / C)) * 10f0
            htg = t.ht_growth[i]
            if !lstart || icr > 0                      # crown.f:527-543 stmt 53 tail (shared with CRVAR/LPIJU)
                chg = crnew - Float32(icr)
                pdifpy = chg / Float32(icr) / fint
                pdifpy > 0.01f0  && (chg = Float32(icr) * 0.01f0 * fint)
                pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
                icri = trunc(Int, (Float32(icr) + (inband ? chg * cmult : chg)) + 0.5f0)   # crown.f:458-463
                if !lstart && icr != 0
                    crln  = h * Float32(icr) / 100f0
                    crmax = (crln + htg) / (h + htg) * 100f0
                    (cmult == 1f0 && icri < 10) && (icri = trunc(Int, crmax + 0.5f0))   # crown.f:485
                    Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
                end
            else
                icri = trunc(Int, crnew + 0.5f0)
                inband && (icri = trunc(Int, Float32(icri) * cmult))   # crown.f:467-468
            end
            lstart && (icri = topkill_icri(t, i, icri))   # crown.f stmt 55 (crown_init.jl)
            icri > 95 && (icri = 95); (cmult == 1f0 && icri < 10) && (icri = 10); icri < 1 && (icri = 1)   # stmt 59
            t.crown_pct[i] = Int32(icri)
            continue
        end
        # crown.f:376 NIVAR: a tree whose BACKDATED diameter is <3" had its crown set by REGENT this cycle ⇒ keep it.
        (!lstart && d - t.diam_growth[i]/bark < 3f0) && continue
        xcrcon = crcon + P[1]*ba + P[2]*ba*ba + P[3]*lnba + P[4]*relden + P[5]*relden*relden + P[6]*lnrd
        pp = t.crown_ratio[i]; pp < 0.01f0 && (pp = 0.01f0)
        pcr = xcrcon + P[7]*d + P[8]*d*d + P[9]*log(d) + P[10]*h + P[11]*h*h + P[12]*log(h) + P[13]*pp + P[14]*log(pp)
        exppcr = exp(pcr)
        if lstart
            icri = trunc(Int, icr + (inband ? cmult * exppcr : exppcr)*100f0 + 0.50005f0)   # crown.f:430-435
            dgsd >= 1.0f0 && (icri = trunc(Int, bachlo(s.rng, Float32(icri), _EM_CRSD)))
        else
            dcrcon = crcon + P[1]*oba + P[2]*oba*oba + P[3]*x1 + P[4]*reldm1 + P[5]*reldm1*reldm1 + P[6]*x2
            db = d - t.diam_growth[i]/bark; db <= 0f0 && (db = d)
            hb = h - t.ht_growth[i]; hb <= 0f0 && (hb = h)
            # crown.f:395-398 P=OLDPCT (the previous cycle's PCT, gradd.f:267; the first cycle's is the backdated
            # CRATET percentile, cratet.f:481), falling back to PCT when OLDPCT<=0. The OLDPCT>PCT-after-thinning
            # branch (ONTREM(7)>0) is not carried (as IE/BC). jl used the CURRENT PCT, which held LL crowns at 55
            # where live FVSem_g16 drew them down 55→53→51→49.
            pb = t.old_crown_pct[i]; pb <= 0f0 && (pb = t.crown_ratio[i])
            pb < 0.01f0 && (pb = 0.01f0)
            dcr = dcrcon + P[7]*db + P[8]*db*db + P[9]*log(db) + P[10]*hb + P[11]*hb*hb + P[12]*log(hb) + P[13]*pb + P[14]*log(pb)
            chg = exppcr - exp(dcr)
            if icr > 0
                pdifpy = chg / Float32(icr) / fint * 100f0
                pdifpy > 0.01f0  && (chg = Float32(icr) * 0.01f0 * fint / 100f0)
                pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint / 100f0)
            end
            icri = trunc(Int, Float32(icr) + (inband ? cmult * chg : chg)*100f0 + 0.50005f0)   # crown.f:430-435
        end
        lstart && (icri = topkill_icri(t, i, icri))   # crown.f stmt 55 (crown_init.jl)
        icri > 95 && (icri = 95); (cmult == 1f0 && icri < 5) && (icri = 5)   # stmt 59 NIVAR, crown.f:528
        t.crown_pct[i] = Int32(icri)
    end
    # em/crown.f:545-598 DO 79 — cycle-0 dead records: CRVAR (11,13-16,19) and LPIJU (6) take the crown-LENGTH
    # forms with ICRI=INT(CR*100.) (truncated); everything else DUBSCR at the record's point, ICRI=INT(CR*100+.5).
    lstart && dub_dead_crowns!(s) do i
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        crvar = sp == 11 || (13 <= sp <= 16) || sp == 19
        lpiju = sp == 6
        if crvar || lpiju
            cl = crvar ? 5.17281f0 + 0.32552f0 * h - 0.01675f0 * ba : -0.59373f0 + 0.67703f0 * h
            cl < 1f0 && (cl = 1f0); cl > h && (cl = h)
            return trunc(Int, (cl / h) * 100f0)
        end
        pt = Int(t.plot_id[i])
        tpccf = (1 <= pt <= length(s.density.point_ccf)) ? s.density.point_ccf[pt] : 0f0
        icri_round(em_dubscr(s.rng, sp, d, h, ba, tpccf, p.avg_height, rmai_v, dgsd))
    end
    return s
end

# =============================================================================
# em_cwcalc — EM forest-grown crown width (base/cwidth.f -> cwcalc.f, IWHO=0; the
# western Bechtold/Crookston/Donnelly crown-width library).  EM is Region 1 =>
# KODFOR<601 => BF=1 and the Region-6 forest-BF section is skipped (GO TO 10).
# EMMAP maps FVS EM species 1..19 -> a 5-char CWEQN (FIA code + eqn#).  This is the
# CRWDTH array that COVER's CVCW reads (NOT the em_tree_ccf CCF polynomial above,
# which is the separate CCF path).  Reuses the shared CR kernels (_cr_r6m2 / _cr_bech1
# / _cr_bech2 / _cr_hopkins) plus the R1 (Crookston Region-1) and pure-power helpers
# below.  Math faithful to gfortran: fpow/fexp/flog + left-to-right op order; the
# Hopkins point is the shared WESTERN one.  Dump-replay bit-exact vs FVSem_g16 (all 15
# EM CWEQN forms, incl. small-tree scaling + juniper D>=25 plateau).
# REPORTING-ONLY (COVER); INERT for growth/mortality/CCF/volume.
# =============================================================================
const _EM_CWMAP = ("10105","07303","20203","11301","07204","06602","10803","09303",
                   "01903","12203","74902","74605","74705","74902","74902","74902",
                   "37506","26405","74902")

# Transcendentals routed DIRECTLY through glibc libm (logf/expf/powf), like WSBWE — the
# gfortran-16 oracle's REAL*4 EXP/LOG/** resolve to these, so this is bit-exact vs
# FVSem_g16 (verified 353/353).  The core/fmath shim is built by a *different* gfortran and
# is ~1 ULP off here, so it is intentionally NOT used on this crown-width path.
const _EMCW_LIBM = "libm.so.6"
@inline _emcw_log(x::Float32) = ccall((:logf, _EMCW_LIBM), Float32, (Float32,), x)
@inline _emcw_exp(x::Float32) = ccall((:expf, _EMCW_LIBM), Float32, (Float32,), x)
@inline _emcw_pow(x::Float32, y::Float32) = ccall((:powf, _EMCW_LIBM), Float32, (Float32,Float32), x, y)

# Crookston Region-1 form: mult*EXP(c0 + ccl*ln(CL) + cd*ln(Dm) + ch*ln(H) + cba*ln(BAREA)).
# Only the ln(D) term floors to `dfloor`; CL/H/BAREA use actual values.  Small-tree scale
# ×(D/dfloor).  Terms with a 0 coefficient add +0 (bit-exact no-op).
@inline function _em_r1(mult::Float32, c0::Float32, ccl::Float32, cd::Float32, ch::Float32,
                        cba::Float32, dfloor::Float32, cap::Float32,
                        d::Float32, h::Float32, cl::Float32, barea::Float32)::Float32
    dm = d >= dfloor ? d : dfloor
    v = mult * _emcw_exp(c0 + ccl*_emcw_log(cl) + cd*_emcw_log(dm) + ch*_emcw_log(h) + cba*_emcw_log(barea))
    d < dfloor && (v *= d / dfloor)
    v > cap && (v = cap)
    return v
end

# Pure power form a*D^b (Crookston R6 model-1 07204, Donnelly 37506).  OMIND=1 scaling.
@inline function _em_powf(a::Float32, b::Float32, cap::Float32, d::Float32)::Float32
    dm = d >= 1f0 ? d : 1f0
    v = a * _emcw_pow(dm, b)
    d < 1f0 && (v *= d)
    v > cap && (v = cap)
    return v
end

# Crookston R6 model 2: a·D^b·H^c·CL^dd·(BAREA+1)^e·EXP(EL)^f (BF=1).  OMIND=1 small-tree
# scaling; EL∈[ello,elhi].  e=0/f=0 → the factor is ^0 = ×1.0 (bit-exact no-op).
@inline function _em_r6m2(a::Float32, b::Float32, c::Float32, dd::Float32, e::Float32, f::Float32,
                          d::Float32, h::Float32, cl::Float32, ba1::Float32, el::Float32,
                          ello::Float32, elhi::Float32, cap::Float32)::Float32
    elc = el < ello ? ello : (el > elhi ? elhi : el)
    dm = d >= 1f0 ? d : 1f0
    cw = a * _emcw_pow(dm, b) * _emcw_pow(h, c) * _emcw_pow(cl, dd) *
         _emcw_pow(ba1, e) * _emcw_pow(_emcw_exp(elc), f)
    d < 1f0 && (cw *= d)
    cw > cap && (cw = cap)
    return cw
end

# Bechtold 2004 model 1: a + b·D (MIND=5 small-tree scaling).
@inline function _em_bech1(a::Float32, b::Float32, d::Float32, cap::Float32)::Float32
    dm = d >= 5f0 ? d : 5f0
    cw = a + b * dm
    d < 5f0 && (cw *= d / 5f0)
    cw > cap && (cw = cap)
    return cw
end

# Bechtold 2004 model 2: a + b·D + c·D² + crc·CR + hic·HI.  HI∈[hlo,hhi]; MIND=5 small-tree
# scaling; optional D≥25 plateau (dcap25).
@inline function _em_bech2(a::Float32, b::Float32, c::Float32, crc::Float32, hic::Float32,
                           d::Float32, cr::Float32, hi::Float32, hlo::Float32, hhi::Float32,
                           cap::Float32, dcap25::Bool)::Float32
    hv = hi < hlo ? hlo : (hi > hhi ? hhi : hi)
    dm = d >= 5f0 ? d : 5f0
    cw = a + b * dm + c * dm * dm + crc * cr + hic * hv
    d < 5f0 && (cw *= d / 5f0)
    (dcap25 && d >= 25f0) && (cw = a + b * 25f0 + c * 25f0 * 25f0 + crc * cr + hic * hv)
    cw > cap && (cw = cap)
    return cw
end

function em_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32,
                   el::Float32, hi::Float32)::Float32
    (1 <= sp <= 19) || return 0.5f0
    barea <= 1f0 && (barea = 1f0)          # cwcalc.f: IF(BAREA.LE.1.) BAREA=1.
    cl  = cr * h * 0.01f0                   # CL = CR*H*0.01
    ba1 = barea + 1f0
    eqn = _EM_CWMAP[sp]
    cw = if eqn == "10105";     _em_r6m2(2.2354f0,0.66680f0,-0.11658f0,0.16927f0,0f0,0f0, d,h,cl,ba1,el,1f0,999f0,40f0)
    elseif eqn == "74605";      _em_r6m2(4.7961f0,0.64167f0,-0.18695f0,0.18581f0,0f0,0f0, d,h,cl,ba1,el,1f0,999f0,45f0)
    elseif eqn == "74705";      _em_r6m2(4.4327f0,0.41505f0,-0.23264f0,0.41477f0,0f0,0f0, d,h,cl,ba1,el,1f0,999f0,56f0)
    elseif eqn == "26405";      _em_r6m2(3.7854f0,0.54684f0,-0.12954f0,0.16151f0,0.03047f0,-0.00561f0, d,h,cl,ba1,el,10f0,79f0,45f0)
    elseif eqn == "07303";      _em_r1(1.02478f0,0.99889f0,0.19422f0,0.59423f0,-0.09078f0,-0.02341f0, 1f0,40f0, d,h,cl,barea)
    elseif eqn == "20203";      _em_r1(1.01685f0,1.48372f0,0.27378f0,0.49646f0,-0.18669f0,-0.01509f0, 1f0,80f0, d,h,cl,barea)
    elseif eqn == "10803";      _em_r1(1.03992f0,1.58777f0,0.30812f0,0.64934f0,-0.38964f0,0f0, 0.7f0,40f0, d,h,cl,barea)
    elseif eqn == "09303";      _em_r1(1.02687f0,1.28027f0,0.2249f0,0.47075f0,-0.15911f0,0f0, 0.1f0,40f0, d,h,cl,barea)
    elseif eqn == "01903";      _em_r1(1.02886f0,1.01255f0,0.30374f0,0.37093f0,-0.13731f0,0f0, 0.1f0,30f0, d,h,cl,barea)
    elseif eqn == "12203";      _em_r1(1.02687f0,1.49085f0,0.1862f0,0.68272f0,-0.28242f0,0f0, 2f0,46f0, d,h,cl,barea)
    elseif eqn == "07204";      _em_powf(2.2586f0,0.68532f0,33f0, d)
    elseif eqn == "37506";      _em_powf(5.8980f0,0.4841f0,25f0, d)
    elseif eqn == "11301";      _em_bech1(4.0181f0,0.8528f0, d,25f0)
    elseif eqn == "06602";      _em_bech2(-4.1599f0,1.3528f0,-0.0233f0,0.0633f0,-0.0423f0, d,cr,hi,-37f0,19f0,29f0,true)
    elseif eqn == "74902";      _em_bech2(4.1687f0,1.5355f0,0f0,0f0,0.1275f0, d,cr,hi,-26f0,-2f0,35f0,false)
    else 0f0 end
    # cwcalc.f final CRWDTH clamp (after label 9000).
    cw < 0.5f0 && (cw = 0.5f0)
    cw > 99.9f0 && (cw = 99.9f0)
    return cw
end
