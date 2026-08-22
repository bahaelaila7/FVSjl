# CA (CentralCalifornia) FFE initial surface-fuel loading — ca/fmcba.f DATA FULIVE/FULIVI/FUINIE/FUINII.
# CA is a California/R5+R6 variant with the same TOP-TWO-cover-species structure as WS/NC (ca/fmcba.f COVCA(2)):
# live + initial-dead fuel pools interpolated from the top two cover-type species (COVCA(1..2), weighted by
# COVCAWT), NOT the single dominant COVTYP the interior variants use. Tables indexed by CA species 1..50;
# columns (herb, shrub) for the live tables, and the 11 FFE size classes for the dead tables. Established
# (60% cover) vs initiating (10% cover); PERCOV-interpolated over XCOV=[10,60]. Values transcribed directly
# from ca/fmcba.f (script-extracted, from J. Brown / Ottmar 2000). Reuses `_cr_algslp2`.

const _CA_FULIVE = Float32[    # established (60% cover) (herb, shrub) — ca/fmcba.f:125
    0.20 0.20;
    0.20 0.20;
    0.20 0.20;
    0.15 0.10;
    0.15 0.10;
    0.15 0.10;
    0.20 0.20;
    0.20 0.20;
    0.15 0.20;
    0.20 0.10;
    0.20 0.10;
    0.20 0.10;
    0.20 0.10;
    0.20 0.10;
    0.20 0.25;
    0.20 0.25;
    0.15 0.10;
    0.20 0.25;
    0.20 0.20;
    0.23 0.22;
    0.14 0.35;
    0.15 0.20;
    0.20 0.20;
    0.20 0.20;
    0.20 0.25;
    0.23 0.22;
    0.25 0.25;
    0.23 0.22;
    0.23 0.22;
    0.23 0.22;
    0.23 0.22;
    0.23 0.22;
    0.23 0.22;
    0.20 0.20;
    0.20 0.20;
    0.20 0.20;
    0.20 0.20;
    0.25 0.25;
    0.20 0.20;
    0.20 0.25;
    0.20 0.20;
    0.25 0.25;
    0.20 0.20;
    0.25 0.25;
    0.25 0.25;
    0.25 0.25;
    0.20 0.20;
    0.20 0.20;
    0.23 0.22;
    0.20 0.20;
]
const _CA_FULIVI = Float32[    # initiating (10% cover) (herb, shrub) — ca/fmcba.f:181
    0.40 2.00;
    0.40 2.00;
    0.40 2.00;
    0.30 2.00;
    0.30 2.00;
    0.30 2.00;
    0.40 2.00;
    0.40 2.00;
    0.30 2.00;
    0.40 1.00;
    0.40 1.00;
    0.40 1.00;
    0.40 1.00;
    0.40 1.00;
    0.25 1.00;
    0.25 1.00;
    0.30 2.00;
    0.25 1.00;
    0.40 2.00;
    0.55 0.35;
    0.10 2.06;
    0.30 2.00;
    0.40 2.00;
    0.40 2.00;
    0.25 1.00;
    0.55 0.35;
    0.18 2.00;
    0.55 0.35;
    0.55 0.35;
    0.55 0.35;
    0.55 0.35;
    0.55 0.35;
    0.55 0.35;
    0.40 2.00;
    0.40 2.00;
    0.40 2.00;
    0.40 2.00;
    0.18 2.00;
    0.40 2.00;
    0.25 1.00;
    0.40 2.00;
    0.18 2.00;
    0.40 2.00;
    0.18 1.32;
    0.18 1.32;
    0.18 1.32;
    0.40 2.00;
    0.40 2.00;
    0.55 0.35;
    0.40 2.00;
]
const _CA_FUINIE = Float32[    # established, 11 size classes — ca/fmcba.f:236
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;
    2.2 2.2 5.2 15 20 15 0 0 0 1 35;
    0.7 0.7 3 7 7 0 0 0 0 0.6 25;
    0.7 0.7 3 7 7 0 0 0 0 0.6 25;
    0.7 0.7 3 7 7 0 0 0 0 0.6 25;
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;
    2.2 2.2 5.2 15 20 15 0 0 0 1 35;
    1.1 1.1 2.2 10 10 0 0 0 0 0.6 30;
    0.7 0.7 1.6 2.5 2.5 0 0 0 0 1.4 5;
    0.7 0.7 1.6 2.5 2.5 0 0 0 0 1.4 5;
    0.7 0.7 1.6 2.5 2.5 0 0 0 0 1.4 5;
    0.7 0.7 1.6 2.5 2.5 0 0 0 0 1.4 5;
    0.7 0.7 1.6 2.5 2.5 0 0 0 0 1.4 5;
    0.7 0.7 1.6 2.5 2.5 0 0 0 0 1.4 5;
    0.7 0.7 1.6 2.5 2.5 0 0 0 0 1.4 5;
    1 1 1.6 10 10 10 0 0 0 0.8 30;
    0.9 0.9 1.2 7 8 0 0 0 0 0.6 15;
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;
    0.3 0.7 1.4 0.2 0.1 0 0 0 0 3.9 0;
    0.1 0.2 0.4 0.5 0.8 1 0 0 0 0.1 0;
    1.1 1.1 2.2 10 10 0 0 0 0 0.6 30;
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;
    0.9 0.9 1.2 7 8 0 0 0 0 0.6 15;
    0.3 0.7 1.4 0.2 0.1 0 0 0 0 3.9 0;
    0.2 0.6 2.4 3.6 5.6 0 0 0 0 1.4 16.8;
    0.3 0.7 1.4 0.2 0.1 0 0 0 0 3.9 0;
    0.3 0.7 1.4 0.2 0.1 0 0 0 0 3.9 0;
    0.3 0.7 1.4 0.2 0.1 0 0 0 0 3.9 0;
    0.3 0.7 1.4 0.2 0.1 0 0 0 0 3.9 0;
    0.3 0.7 1.4 0.2 0.1 0 0 0 0 3.9 0;
    0.3 0.7 1.4 0.2 0.1 0 0 0 0 3.9 0;
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;
    0.2 0.6 2.4 3.6 5.6 0 0 0 0 1.4 16.8;
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;
    0.9 0.9 1.2 7 8 0 0 0 0 0.6 15;
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;
    0.2 0.6 2.4 3.6 5.6 0 0 0 0 1.4 16.8;
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;
    0.2 0.6 2.4 3.6 5.6 0 0 0 0 1.4 16.8;
    0.2 0.6 2.4 3.6 5.6 0 0 0 0 1.4 16.8;
    0.2 0.6 2.4 3.6 5.6 0 0 0 0 1.4 16.8;
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;
    0.3 0.7 1.4 0.2 0.1 0 0 0 0 3.9 0;
    0.9 0.9 1.6 3.5 3.5 0 0 0 0 0.6 10;
]
const _CA_FUINII = Float32[    # initiating, 11 size classes — ca/fmcba.f:294
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;
    1.6 1.6 3.6 6 8 6 0 0 0 0.5 12;
    0.5 0.5 2 2.8 2.8 0 0 0 0 0.3 12;
    0.5 0.5 2 2.8 2.8 0 0 0 0 0.3 12;
    0.5 0.5 2 2.8 2.8 0 0 0 0 0.3 12;
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;
    1.6 1.6 3.6 6 8 6 0 0 0 0.5 12;
    0.7 0.7 1.6 4 4 0 0 0 0 0.3 12;
    0.1 0.1 0.2 0.5 0.5 0 0 0 0 0.5 0.8;
    0.1 0.1 0.2 0.5 0.5 0 0 0 0 0.5 0.8;
    0.1 0.1 0.2 0.5 0.5 0 0 0 0 0.5 0.8;
    0.1 0.1 0.2 0.5 0.5 0 0 0 0 0.5 0.8;
    0.1 0.1 0.2 0.5 0.5 0 0 0 0 0.5 0.8;
    0.1 0.1 0.2 0.5 0.5 0 0 0 0 0.5 0.8;
    0.1 0.1 0.2 0.5 0.5 0 0 0 0 0.5 0.8;
    0.6 0.6 0.8 6 6 6 0 0 0 0.4 12;
    0.6 0.7 0.8 2.8 3.2 0 0 0 0 0.3 7;
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;
    0.1 0.1 0 0 0 0 0 0 0 2.9 0;
    0.2 0.4 0.2 0 0 0 0 0 0 0.2 0;
    0.7 0.7 1.6 4 4 0 0 0 0 0.3 12;
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;
    0.6 0.7 0.8 2.8 3.2 0 0 0 0 0.3 7;
    0.1 0.1 0 0 0 0 0 0 0 2.9 0;
    0.1 0.4 5 2.2 2.3 0 0 0 0 0.8 5.6;
    0.1 0.1 0 0 0 0 0 0 0 2.9 0;
    0.1 0.1 0 0 0 0 0 0 0 2.9 0;
    0.1 0.1 0 0 0 0 0 0 0 2.9 0;
    0.1 0.1 0 0 0 0 0 0 0 2.9 0;
    0.1 0.1 0 0 0 0 0 0 0 2.9 0;
    0.1 0.1 0 0 0 0 0 0 0 2.9 0;
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;
    0.1 0.4 5 2.2 2.3 0 0 0 0 0.8 5.6;
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;
    0.6 0.7 0.8 2.8 3.2 0 0 0 0 0.3 7;
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;
    0.1 0.4 5 2.2 2.3 0 0 0 0 0.8 5.6;
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;
    0.1 0.4 5 2.2 2.3 0 0 0 0 0.8 5.6;
    0.1 0.4 5 2.2 2.3 0 0 0 0 0.8 5.6;
    0.1 0.4 5 2.2 2.3 0 0 0 0 0.8 5.6;
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;
    0.1 0.1 0 0 0 0 0 0 0 2.9 0;
    0.5 0.5 1 1.4 1.4 0 0 0 0 0.3 5;
]

@inline function ca_live_fuel_loading(covca::NTuple{2,Int}, covcawt::NTuple{2,Float32}, percov::Float32)::NTuple{2,Float32}
    herb = 0f0; shrub = 0f0
    @inbounds for j in 1:2
        c = covca[j]; (1 <= c <= 50) || continue; wt = covcawt[j]
        herb  += _cr_algslp2(percov, 10f0, 60f0, _CA_FULIVI[c, 1] * wt, _CA_FULIVE[c, 1] * wt)
        shrub += _cr_algslp2(percov, 10f0, 60f0, _CA_FULIVI[c, 2] * wt, _CA_FULIVE[c, 2] * wt)
    end
    return (herb, shrub)
end

function ca_dead_fuel_loading(covca::NTuple{2,Int}, covcawt::NTuple{2,Float32}, percov::Float32)::Vector{Float32}
    out = zeros(Float32, 11)
    @inbounds for isz in 1:11, j in 1:2
        c = covca[j]; (1 <= c <= 50) || continue; wt = covcawt[j]
        out[isz] += _cr_algslp2(percov, 10f0, 60f0, _CA_FUINII[c, isz] * wt, _CA_FUINIE[c, isz] * wt)
    end
    return out
end

# =============================================================================
# ca_cwcalc — CA crown width (ft) for FMCBA's PERCOV. cat01 is forest 610 = R6 (IFOR=6 > 5), so CA uses the
# Crookston R6 CAMAP path (ca/cwcalc.f), NOT R5CRWD (ca/cwcalc.f:385 branches R5CRWD only for IFOR≤5). Without
# it the generic crown_width=0.5 collapses PERCOV≈0 → the wrong initiating-stand fuel loads (measured PERCOV
# 0.27). CA_CWMAP maps CA species 1..50 → a 5-char CWEQN (FIA code + eqn#) from ca/cwcalc.f DATA CAMAP. The five
# cat01 species (DF/WF/SP/LP/PP) use the Crookston-R6 model-2 form (_cr_r6m2); DF/WF/SP coefficients are IDENTICAL
# to NC/WS (20205/01505/11705). LP(10805)/PP(12205) added here. Errors on CA species whose CWEQN is not yet ported.
const _CA_CWMAP = ("04105","08105","24205","01505","02006","02105","20205","26305","26403","10105",
                   "10305","10805","10805","11301","11605","11705","11905","12205","12702","12702",
                   "06405","09204","21104","23104","11605","80102","80502","80702","80702","81505",
                   "81802","82102","83902","31206","31206","35106","36102","63102","35106","31206",
                   "31206","63102","63102","74605","74705","31206","98102","98102","31206","21104")

function ca_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32, el::Float32, hi::Float32)::Float32
    (1 <= sp <= 50) || return 0f0
    eqn = _CA_CWMAP[sp]; cl = cr * h * 0.01f0; ba1 = barea + 1f0
    # Crookston-R6 model 2, BF folded into the leading coef (cwcalc.f CASE(610,710,711) Rogue River — cat01's
    # forest 610; same BF table as OC/oc_cwcalc, whose ref forest is also 610). BF: DF(202)=1.0, WF(015)=1.0
    # (neither in the 610 table), SP(117)=1.048, LP(108)=0.944, PP(122)=0.918. Without BF the per-tree CrWidth was
    # ~5% off (the aggregated FFE crown-biomass masked it); folding it in makes the FVS_TreeList column bit-exact.
    if     eqn == "20205"; return _cr_r6m2(6.0227f0*1.000f0,0.54361f0,-0.20669f0,0.20395f0,-0.00644f0,-0.00378f0, d,h,cl,ba1,el, 1f0,75f0,80f0)  # Douglas-fir (BF 1.0)
    elseif eqn == "01505"; return _cr_r6m2(5.0312f0*1.000f0,0.53680f0,-0.18957f0,0.16199f0, 0.04385f0,-0.00651f0, d,h,cl,ba1,el, 2f0,75f0,35f0)  # white fir (BF 1.0)
    elseif eqn == "11705"; return _cr_r6m2(3.5930f0*1.048f0,0.63503f0,-0.22766f0,0.17827f0, 0.04267f0,-0.00290f0, d,h,cl,ba1,el, 5f0,75f0,56f0)  # sugar pine (BF 1.048)
    elseif eqn == "10805"; return _cr_r6m2(6.6941f0*0.944f0,0.81980f0,-0.36992f0,0.17722f0,-0.01202f0,-0.00882f0, d,h,cl,ba1,el, 1f0,79f0,40f0)  # lodgepole pine (BF 0.944)
    elseif eqn == "12205"; return _cr_r6m2(4.7762f0*0.918f0,0.74126f0,-0.28734f0,0.17137f0,-0.00602f0,-0.00209f0, d,h,cl,ba1,el, 13f0,75f0,50f0) # ponderosa pine (BF 0.918)
    elseif eqn == "09204"; return _nc_donnelly(2.8232f0,0.66326f0, d, 38f0)   # bristlecone/BR (Donnelly R6)
    else
        error("ca_cwcalc: crown-width equation $(eqn) (CA species $(sp)) not yet ported — cat01 exercises only " *
              "DF/WF/SP/LP/PP; the remaining CAMAP equations are a follow-up crown-width chunk.")
    end
end
