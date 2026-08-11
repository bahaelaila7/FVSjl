# =============================================================================
# height_growth.jl (inlandempire) — IE large-tree height growth (ie/htgf.f).
#
# NI conifers (sp<=12,14,23) — identical form to KT:
#   CON = HTCON(sp) + H2COF*HTI^2 + HGLD(sp)*ln(D) + HGLH*ln(HTI)
#   HTG = exp(CON + HDGCOF*ln(DG)) + BIAS;  max(0.1)
# sp15,16 (PM/RM pinyon-juniper): HTG = 0 (no large-tree height growth).
# sp13,17 (LM/PY) + sp18-22 (aspen/CO/MM/PB/OH): COFLM/COFAS Weibull curve (ie/htgf.f:90-160),
#   K by ICR crown class; degenerate guards -> 0.1; DIA=D+DG/bark; PSI/H -> HTG=H-HTI.
# Tail (all): HTG*SCALE*XHT; sp{13,15:22} also *exp(HTCON) (HCOR2 small-tree calib); *MISHGF(=1);
#   SIZCAP(sp,4) height cap.  Per-stand HTCONS: IHT=MAPHAB(ITYPE) -> HGHCH/H2COF/HDGCOF; HTCON(sp)=
#   HGHCH+HGSC(sp) for NI else 0 (+ln(HCOR2) when LHCOR2).
# =============================================================================

const IE_HGLD = Float32[-.04935,-.3899,-.4574,-.09775,-.1555,-.1219,-.2454,-.5720,-.1997,-.5657,-.1219,
                        -.3899,0.0,-.1997,0,0,0,0,0,0,0,0,-.1219]
const IE_HGSC = Float32[-.5342,.1433,.1641,-.6458,-.6959,-.9941,-.6004,.2089,-.5478,.7316,-.9941,
                        .1433,0.0,-.5478,0,0,0,0,0,0,0,0,-.9941]
const IE_HGLH = 0.23315f0
const IE_HTBIAS = 0.4809f0
# MAPHAB(30): ITYPE -> IHT (1..8). Identical to KT_HTMAPHAB.
const IE_HTMAPHAB = Int[1,1, 2,2,2,2,2,2,2, 3,3,4,5,6, 7,7,7,7, 4,4,1,4,4, 8,8,8, 1,1,1,1]
const IE_HGHC  = Float32[2.03035, 1.72222, 1.19728, 1.81759, 2.14781, 1.76998, 2.21104, 1.74090]
const IE_HGLDD = Float32[0.62144, 1.02372, 0.85493, 0.75756, 0.46238, 0.49643, 0.37042, 0.34003]
const IE_HGH2  = Float32[-13.358f-5, -3.809f-5, -3.715f-5, -2.607f-5, -5.200f-5, -1.605f-5, -3.631f-5, -4.460f-5]
# COFLM/COFAS (9,3): Weibull height-curve coefficients (by ICR crown class K). ie/htgf.f DATA.
const IE_COFLM = reshape(Float32[
    37.0,85.0,1.77836,-0.51147,1.88795,1.20654,0.57697,3.57635,0.90283,
    45.0,100.0,1.66674,0.25626,1.45477,1.11251,0.67375,2.17942,0.88103,
    45.0,90.0,1.64770,0.30546,1.35015,0.94823,0.70453,2.46480,1.00316], 9, 3)
const IE_COFAS = reshape(Float32[
    30.0,85.0,2.00995,0.03288,1.81059,1.28612,0.72051,3.00551,1.01433,
    30.0,85.0,2.00995,0.03288,1.81059,1.28612,0.72051,3.00551,1.01433,
    35.0,85.0,1.80388,-0.07682,1.70032,1.29148,0.72343,2.91519,0.95244], 9, 3)

"""height_growth!(s, ::InlandEmpire; scale) — per-tree periodic height increment into trees.ht_growth."""
function height_growth!(s::StandState, ::InlandEmpire; scale::Float32 = 1.0f0)
    p, t, ctl = s.plot, s.trees, s.control
    itype = Int(p.habitat_input)
    iht = (1 <= itype <= 30) ? IE_HTMAPHAB[itype] : 1
    hghch  = IE_HGHC[iht]
    h2cof  = IE_HGH2[iht]
    hdgcof = IE_HGLDD[iht]
    relden = p.relative_density
    iage = Int(p.stand_age)
    cur_year = current_cycle_year(s)
    icyc = Int(ctl.cycle)
    fint = scale * 10f0                                    # SCALE=FINT/YR, YR=10 => FINT = scale*10
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0.0f0
        t.tpa[i] <= 0.0f0 && continue
        sp = Int(t.species[i])
        d = t.dbh[i]; hti = t.height[i]
        (d <= 0.0f0 || hti <= 0.0f0) && continue
        dg = t.diam_growth[i]
        # HTCON(sp): NI = HGHCH+HGSC(sp), special = 0; +ln(HCOR2) when LHCOR2 (small-tree calib).
        htcon_ni = hghch + IE_HGSC[sp]
        htcon = (sp <= 12 || sp == 14 || sp == 23) ? htcon_ni : 0.0f0
        (ctl.htg_cor2_on && ctl.htg_cor2[sp] > 0.0f0) && (htcon += log(ctl.htg_cor2[sp]))
        con = 0.0f0
        htg = 0.0f0
        if sp <= 12 || sp == 14 || sp == 23
            dg <= 0.0f0 && continue                        # ln(DG) undefined
            con = htcon + h2cof * hti * hti + IE_HGLD[sp] * log(d) + IE_HGLH * log(hti)
            htg = exp(con + hdgcof * log(dg)) + IE_HTBIAS
            htg < 0.1f0 && (htg = 0.1f0)
        elseif sp == 15 || sp == 16
            continue                                       # PI/JU: HTG stays 0
        else
            # LM/PY (13/17) + aspen (18-22) Weibull curve
            iicr = trunc(Int, Float32(t.crown_pct[i]) / 10.0f0 + 0.5f0)
            iicr > 9 && (iicr = 9)
            k = iicr <= 2 ? 1 : iicr <= 7 ? 2 : 3
            cof = (sp == 13 || sp == 17) ? IE_COFLM : IE_COFAS
            c1=cof[1,k]; c2=cof[2,k]; c3=cof[3,k]; c4=cof[4,k]; c5=cof[5,k]; c6=cof[6,k]; c7=cof[7,k]; c8=cof[8,k]; c9=cof[9,k]
            if hti <= 4.5f0 || (0.1f0 + c1) <= d || (4.5f0 + c2) <= hti
                htg = 0.1f0
            else
                temd = d <= 0.2f0 ? 0.2f0 : d
                y1 = (temd - 0.1f0) / c1
                y2 = (hti - 4.5f0) / c2
                fby1 = log(y1 / (1.0f0 - y1))
                fby2 = log(y2 / (1.0f0 - y2))
                z = (c4 + c6 * fby2 - c7 * (c3 + c5 * fby1)) * (1.0f0 - c7^2)^(-0.5f0)
                if sp != 13 && sp != 17
                    zadj = 0.1f0 - 0.10273f0 * z + 0.00273f0 * z * z
                    zadj < 0.0f0 && (zadj = 0.0f0)
                    z = z + zadj
                end
                # young-lodgepole/aspen accelerator (age 10..40, D<9, cycle 1)
                if iage != 0 && icyc == 0                  # FVS ICYC==1 only (jl cycle 0-based); <=1 double-fires
                    ixage = iage + cur_year - Int(ctl.cycle_year[1])
                    if ixage < 40 && ixage > 10 && d < 9.0f0 && z <= 2.0f0
                        zadj = 0.3564f0 * dg * fint / 10f0
                        closur = relden < 100.0f0 ? 1.0f0 : Float32(t.crown_ratio[i]) / 100.0f0
                        zadj = zadj * closur
                        (iicr == 9 || iicr == 8) && (zadj = zadj * 1.1f0)
                        z = z + zadj
                        z > 2.0f0 && (z = 2.0f0)
                    end
                end
                bark = ie_bratio(sp, d)
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
        # tail: SCALE*XHT; sp{13,15:22} also *exp(HTCON); *MISHGF(=1); SIZCAP cap
        xht = active_multiplier(ctl, :htg, sp, cur_year)
        if sp == 13 || (15 <= sp <= 22)
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
