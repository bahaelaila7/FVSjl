# EM establishment base height (em/essubh.f — subsequent/planted "tallest" tree height). Ported deterministically:
# HHT = EXP(PN), the mean of FVS's HHT=EXP(PN + EMSQR·DILATE·BNORML·σ) since EMSQR=±DRAW is symmetric (mean≈0);
# the stochastic order-statistic spread is DEFERRED (like CI's intentional ZZRAN deferral). Validated bit-exact
# vs FVSem instrument-replay on em_plant.key (DF: PN=0.156307, HHT=1.18492). See docs/ESTB_SUBSYSTEM_CHUNK_PLAN.md.
#
# Inputs (em/estab.f:470-493): L=ALOG(AGE); B=BAA competition-BA clamp[1,400]; XC=SLO·cos(asp)/XS=SLO·sin(asp);
# E=ELEV(100s ft); IHTSER=height series via the habitat-code bracket search (below); IPHY=physiographic pos
# (default 3, esplt2.f:192/223); IPREP=site-prep (default 1=NONE; per-record variation deferred).

# IHTSER from the habitat code (esplt2.f:231-238 bracket search → MYGRUP → em/estab.f:493 MYHTS).
const _EM_ES_IEND = Int32[269,299,319,335,385,394,399,499,509,515,519,522,523,524,529,564,579,584,589,599,
                          634,637,644,649,659,669,689,699,709,719,739,744,799]
const _EM_ES_MYGRUP = Int32[3,1,4,2,4,3,4,3,8,6,8,7,5,7,8,9,10,6,8,5,13,16,11,14,16,12,15,12,14,15,11,15,14]
const _EM_ES_MYHTS = Int32[1,2,2,2,3,3,3,3,4,4,5,5,5,5,5,5]   # em/estab.f:111 MYHTS/1,3*2,4*3,2*4,6*5/

@inline function em_ihtser(kodtyp::Integer)::Int
    ihab = 16
    @inbounds for j in 1:length(_EM_ES_IEND)
        if kodtyp <= _EM_ES_IEND[j]; ihab = Int(_EM_ES_MYGRUP[j]); break; end
    end
    (ihab < 1 || ihab > 16) && (ihab = 16)
    Int(_EM_ES_MYHTS[ihab])
end

# UHAB(5 series,·)/UPRE(4 prep,·)/UPHY(5 phys,·) — only the species with nonzero rows (em/essubh.f DATA); others 0.
const _EM_UHAB = Dict(2=>(-0.01541f0,-0.03814f0,0.11409f0,0.35334f0,0.0f0),
                      3=>(-0.21858f0,-0.03354f0,0.22756f0,0.51988f0,0.0f0),
                      7=>(-0.29969f0,-0.15449f0,0.04545f0,-0.00601f0,0.0f0),
                      8=>(0.0f0,0.0f0,0.18740f0,0.26511f0,0.0f0),
                      10=>(-0.02287f0,-0.14710f0,0.19278f0,0.13817f0,0.0f0))
const _EM_UPRE = Dict(2=>(0.0f0,-0.11310f0,-0.06246f0,0.009632f0),
                      3=>(0.0f0,0.06961f0,0.19508f0,0.17952f0),
                      7=>(0.0f0,0.11502f0,0.02486f0,0.13080f0),
                      8=>(0.0f0,0.10587f0,0.27072f0,0.16240f0),
                      10=>(0.0f0,0.20729f0,0.18491f0,0.11864f0))
const _EM_UPHY = Dict(1=>(-0.18731f0,-0.48682f0,-0.32160f0,-0.16113f0,0.0f0),
                      3=>(-0.27801f0,-0.20433f0,-0.12317f0,-0.26736f0,0.0f0),
                      7=>(0.32401f0,0.14743f0,0.22165f0,0.24559f0,0.0f0),
                      8=>(0.41120f0,0.01164f0,0.22217f0,0.15834f0,0.0f0))
@inline _em_uhab(sp,ih) = haskey(_EM_UHAB,sp) ? _EM_UHAB[sp][ih] : 0.0f0
@inline _em_upre(sp,ip) = haskey(_EM_UPRE,sp) ? _EM_UPRE[sp][ip] : 0.0f0
@inline _em_uphy(sp,iy) = haskey(_EM_UPHY,sp) ? _EM_UPHY[sp][iy] : 0.0f0

# Deterministic EM subsequent/planted base height (em/essubh.f, HHT=EXP(PN)). agel=ALOG(AGE), baa clamp[1,400],
# xc=SLO·cos(asp), xs=SLO·sin(asp), slo, elev(100s ft), ih=IHTSER(1-5), iy=IPHY(1-5), ip=IPREP(1-4).
function em_essubh_hht(sp::Int, agel::Float32, baa::Float32, xc::Float32, xs::Float32,
                       slo::Float32, elev::Float32, ih::Int, iy::Int, ip::Int)::Float32
    if sp == 1                       # WB (NI WP)
        return exp(-1.51302f0 + 1.24537f0*agel - 0.003052f0*baa + _em_uphy(1,iy))
    elseif sp == 2                   # WL
        return exp(-1.36257f0 + 1.21548f0*agel - 0.003797f0*baa + _em_uhab(2,ih) + _em_upre(2,ip))
    elseif sp == 3                   # DF [validated]
        return exp(-2.16416f0 + 1.28151f0*agel - 0.0031363f0*baa + _em_uhab(3,ih) + _em_upre(3,ip) +
                   _em_uphy(3,iy) - 0.09626f0*xc - 0.23946f0*xs - 0.14589f0*slo)
    elseif sp == 5 || sp == 9        # LL / AF (NI AF)
        return exp(-2.06377f0 + 1.18184f0*agel - 0.0044465f0*baa + 0.06615f0*xc + 0.03085f0*xs - 0.37402f0*slo)
    elseif sp == 7                   # LP (BWAF/BWB4 habitat flags deferred = 0)
        return exp(-0.27105f0 + 1.32027f0*agel - 0.008208f0*baa + _em_upre(7,ip) + _em_uphy(7,iy) + _em_uhab(7,ih) -
                   0.15385f0*xc + 0.04156f0*xs - 0.49186f0*slo - 0.04744f0*elev + 0.0003511f0*elev*elev)
    elseif sp == 8                   # ES
        return exp(-2.93213f0 + 1.43503f0*agel - 0.002504f0*baa + _em_upre(8,ip) + _em_uphy(8,iy) + _em_uhab(8,ih))
    elseif sp == 10                  # PP
        return exp(-1.99480f0 + 1.53946f0*agel - 0.00402f0*baa + _em_uhab(10,ih) + _em_upre(10,ip) - 0.01155f0*elev)
    elseif sp == 18                  # OS (NI OT)
        return exp(-2.42379f0 + 1.52366f0*agel - 0.003256f0*baa)
    elseif sp == 4 || sp == 6        # LM / RM — fixed
        return 0.5f0
    else                             # GA/AS/CW/BA/PW/NC/PB/OH (11-17,19) — fixed
        return 5.0f0
    end
end
