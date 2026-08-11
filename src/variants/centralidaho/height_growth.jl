# =============================================================================
# height_growth.jl (centralidaho) — CI large-tree height growth (ci/htgf.f). Chunk 4.
# NI conifers (DEFAULT: sp 1-10,18) — CON = HTCON(sp) + H2COF·HTI² + HGLD·ln(D) + HGLH·ln(HTI);
#   HTG = exp(CON + HDGCOF·ln(DG)) + BIAS; max(0.1). HTCON(sp)=HGHC(ITYPE)+HGSC(sp).
# Weibull-curve species (11,12,13,16,17,19): COFLM(11,12,16) / COFAS(13,17,19) by ICR crown class.
# WJ(14),MC(15): same exp form but HTCON=0 (HGLD=0). ITYPE=NIHMAP[ICINDX]; HGHC/HGLDD/HGH2 are (30).
# Tail: sp 11-17,19 also ×exp(HTCON) (HCOR2 small-tree calib); ×SCALE·XHT; SIZCAP height cap.
# =============================================================================

const CI_HGLD = Float32[-0.04935,-0.3899,-0.4574,-0.09775,-0.1555,-0.1219,-0.2454,-0.572,-0.1997,-0.5657,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.1219,0.0]
const CI_HGSC = Float32[-0.5342,0.1433,0.1641,-0.6458,-0.6959,-0.9941,-0.6004,0.2089,-0.5478,0.7316,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.9941,0.0]
const CI_HGHC = Float32[2.03035,2.03035,1.72222,1.72222,1.72222,1.72222,1.72222,1.72222,1.72222,1.19728,1.19728,1.81759,2.14781,1.76998,2.21104,2.21104,2.21104,2.21104,1.81759,1.81759,2.03035,1.81759,1.81759,1.7409,1.7409,1.7409,2.03035,2.03035,2.03035,2.03035]
const CI_HGLDD = Float32[0.62144,0.62144,1.02372,1.02372,1.02372,1.02372,1.02372,1.02372,1.02372,0.85493,0.85493,0.75756,0.46238,0.49643,0.37042,0.37042,0.37042,0.37042,0.75756,0.75756,0.62144,0.75756,0.75756,0.34003,0.34003,0.34003,0.62144,0.62144,0.62144,0.62144]
const CI_HGH2 = Float32[-0.00013358,-0.00013358,-3.809e-05,-3.809e-05,-3.809e-05,-3.809e-05,-3.809e-05,-3.809e-05,-3.809e-05,-3.715e-05,-3.715e-05,-2.607e-05,-5.2e-05,-1.605e-05,-3.631e-05,-3.631e-05,-3.631e-05,-3.631e-05,-2.607e-05,-2.607e-05,-0.00013358,-2.607e-05,-2.607e-05,-4.46e-05,-4.46e-05,-4.46e-05,-0.00013358,-0.00013358,-0.00013358,-0.00013358]
const CI_HGLH = 0.23315f0
const CI_HTBIAS = 0.4809f0
const CI_COFLM = reshape(Float32[37.0,85.0,1.77836,-0.51147,1.88795,1.20654,0.57697,3.57635,0.90283,45.0,100.0,1.66674,0.25626,1.45477,1.11251,0.67375,2.17942,0.88103,45.0,90.0,1.6477,0.30546,1.35015,0.94823,0.70453,2.4648,1.00316],9,3)
const CI_COFAS = reshape(Float32[30.0,85.0,2.00995,0.03288,1.81059,1.28612,0.72051,3.00551,1.01433,30.0,85.0,2.00995,0.03288,1.81059,1.28612,0.72051,3.00551,1.01433,35.0,85.0,1.80388,-0.07682,1.70032,1.29148,0.72343,2.91519,0.95244],9,3)

@inline _ci_is_weibull(sp) = sp == 11 || sp == 12 || sp == 13 || sp == 16 || sp == 17 || sp == 19

function height_growth!(s::StandState, ::CentralIdaho; scale::Float32 = 1.0f0)
    p, t, ctl = s.plot, s.trees, s.control
    sd = s.coef.species
    icindx = Int(p.habitat_input)
    itype = (1 <= icindx <= 130) ? Int(CI_NIHMAP[icindx]) : 1
    (itype < 1 || itype > 30) && (itype = 1)
    hghch = CI_HGHC[itype]; h2cof = CI_HGH2[itype]; hdgcof = CI_HGLDD[itype]
    relden = p.relative_density
    iage = Int(p.stand_age)
    cur_year = current_cycle_year(s); icyc = Int(ctl.cycle)
    fint = scale * 10f0
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0.0f0
        t.tpa[i] <= 0.0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; hti = t.height[i]
        (d <= 0.0f0 || hti <= 0.0f0) && continue
        dg = t.diam_growth[i]
        # HTCON(sp): NI (1-10,18) = HGHCH+HGSC(sp); 11-17,19 = 0. +ln(HCOR2) when LHCOR2.
        htcon = (_ci_is_weibull(sp) || sp == 14 || sp == 15) ? 0.0f0 : hghch + CI_HGSC[sp]
        (ctl.htg_cor2_on && ctl.htg_cor2[sp] > 0.0f0) && (htcon += log(ctl.htg_cor2[sp]))
        htg = 0.0f0
        if !_ci_is_weibull(sp)
            dg <= 0.0f0 && continue                            # ln(DG) undefined
            con = htcon + h2cof * hti * hti + CI_HGLD[sp] * log(d) + CI_HGLH * log(hti)
            htg = exp(con + hdgcof * log(dg)) + CI_HTBIAS
            htg < 0.1f0 && (htg = 0.1f0)
        else
            iicr = trunc(Int, Float32(t.crown_pct[i]) / 10.0f0 + 0.5f0)
            iicr > 9 && (iicr = 9)
            k = iicr <= 2 ? 1 : iicr <= 7 ? 2 : 3
            cof = (sp == 11 || sp == 12 || sp == 16) ? CI_COFLM : CI_COFAS   # COFLM for WB/PY/LM, COFAS for AS/CW/OH
            c1=cof[1,k]; c2=cof[2,k]; c3=cof[3,k]; c4=cof[4,k]; c5=cof[5,k]; c6=cof[6,k]; c7=cof[7,k]; c8=cof[8,k]; c9=cof[9,k]
            if hti <= 4.5f0 || (0.1f0 + c1) <= d || (4.5f0 + c2) <= hti
                htg = 0.1f0
            else
                temd = d <= 0.2f0 ? 0.2f0 : d
                y1 = (temd - 0.1f0) / c1; y2 = (hti - 4.5f0) / c2
                fby1 = log(y1 / (1.0f0 - y1)); fby2 = log(y2 / (1.0f0 - y2))
                z = (c4 + c6 * fby2 - c7 * (c3 + c5 * fby1)) * (1.0f0 - c7^2)^(-0.5f0)
                if sp == 13 || sp == 17 || sp == 19                # aspen-group z-adjust (COFAS)
                    zadj = 0.1f0 - 0.10273f0 * z + 0.00273f0 * z * z
                    zadj < 0.0f0 && (zadj = 0.0f0)
                    z = z + zadj
                end
                if iage != 0 && icyc == 0                  # FVS ICYC==1 only (jl cycle 0-based); <=1 double-fires
                    ixage = iage + cur_year - Int(ctl.cycle_year[1])
                    if ixage < 40 && ixage > 10 && d < 9.0f0 && z <= 2.0f0
                        zadj = 0.3564f0 * dg * fint / 10f0
                        closur = relden < 100.0f0 ? 1.0f0 : Float32(t.crown_ratio[i]) / 100.0f0
                        zadj *= closur
                        (iicr == 9 || iicr == 8) && (zadj *= 1.1f0)
                        z += zadj; z > 2.0f0 && (z = 2.0f0)
                    end
                end
                bark = ci_bratio(sd, sp, d)
                dia = d + dg / bark
                if (0.1f0 + c1) > dia
                    psi = c8 * ((dia - 0.1f0) / (0.1f0 + c1 - dia))^c9 * exp(z * ((1.0f0 - c7^2))^0.5f0 / c6)
                    h = ((psi / (1.0f0 + psi)) * c2) + 4.5f0
                    h < hti && (h = hti)
                    htg = h - hti
                else
                    htg = 0.1f0
                end
            end
        end
        xht = active_multiplier(ctl, :htg, sp, cur_year)
        if _ci_is_weibull(sp) || sp == 14 || sp == 15            # 11-17,19 → ×exp(HTCON)
            htg = htg * scale * xht * exp(htcon)
        else
            htg = htg * scale * xht
        end
        cap = ctl.sp_size_cap[sp, 4]
        if hti + htg > cap
            htg = cap - hti; htg < 0.1f0 && (htg = 0.1f0)
        end
        t.ht_growth[i] = htg
    end
    return s
end
