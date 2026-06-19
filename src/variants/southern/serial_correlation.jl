# =============================================================================
# serial_correlation.jl — stochastic diameter-growth serial correlation
#
# Ported from: base/autcor.f (ARMA(1,1) variance/covariance multipliers),
# base/dgscor.f (per-tree auto-correlated DG error), and the OLDRN setup in
# sn/dgdriv.f (calibration pass).
#
# SN runs diameter growth STOCHASTICALLY: grinit sets DGSD=2, so each tree's
# ln(DDS) is perturbed by a bounded-normal serial-correlation factor frm=exp(raw),
# where raw = BACHLO(0,ssigma)*rhocp + rho*OLDRN_prev, drawn from the faithful RNG.
# OLDRN (the prior-period residual) is seeded in the calibration pass — measured
# trees get the regression residual, calibrated species fill non-measured trees by
# regression, and uncalibrated species draw OLDRN from BACHLO. Because this consumes
# the RNG, the per-tree iteration MUST follow FVS's species-sorted order (IND1/ISCT).
# =============================================================================

const DG_BJPHI  = 0.74f0      # Box–Jenkins AR parameter   (BJPHI, sn/grinit.f)
const DG_BJTHET = 0.42f0      # Box–Jenkins MA parameter   (BJTHET)
const DG_DGSD   = 2.0f0       # std-dev bound on DG variance (DGSD, sn/grinit.f)

# ARMA(1,1) autocorrelation series BJRHO[1..40] (autcor.f LSTART init).
function _dg_bjrho()
    r = zeros(Float32, 40)
    r[1] = (1f0 - DG_BJPHI * DG_BJTHET) * (DG_BJPHI - DG_BJTHET) /
           (1f0 + DG_BJTHET * (DG_BJTHET - 2f0 * DG_BJPHI))
    @inbounds for i in 2:40
        r[i] = r[i-1] * DG_BJPHI
    end
    return r
end
const DG_BJRHO = _dg_bjrho()

"""
    autcor(new_period, old_period) -> (cov, vrnext)

ARMA(1,1) variance (`vrnext`=VMLT) and covariance (`cov`=COVMLT) multipliers for
the random DG component across a growth cycle of length `new_period`, preceded by
one of length `old_period` (years). autcor.f.
"""
function autcor(newv::Integer, oldv::Integer)
    nv = Int(newv); ov = Int(oldv)
    var = 0f0
    @inbounds for i in 1:(nv - 1)
        var += DG_BJRHO[i] * Float32(nv - i)
    end
    vrnext = Float32(nv) + 2f0 * var

    l = nv + ov - 1; l > 40 && (l = 40)
    nbig = max(nv, ov); nsml = min(nv, ov)
    t = 0f0; dt = 1f0; covar = 0f0
    @inbounds for i in 1:l
        t += dt
        covar += DG_BJRHO[i] * t
        i == nsml && (dt = 0f0)
        i == nbig && (dt = -1f0)
    end
    return covar, vrnext
end

"""
    dgscor!(rng, oldrn, it, ssig, rho, rhocp, wk2_it) -> frm

DGSCOR: draw the per-tree auto-correlated DG error multiplier `frm=exp(raw)` and
store the raw residual back into `oldrn[it]` (for next cycle). `wk2_it` is the
tree's ln(DDS) prediction; large predictions taper the error to 0. dgscor.f.
"""
@inline function dgscor!(rng::FVSRng, oldrn::AbstractVector{Float32}, it::Integer,
                         ssig::Real, rho::Real, rhocp::Real, wk2_it::Real)::Float32
    frm = 0f0
    if DG_DGSD >= 1f0
        bound = DG_DGSD * Float32(ssig)
        while true
            frm = bachlo(rng, 0f0, Float32(ssig))
            frm = frm * Float32(rhocp) + Float32(rho) * oldrn[it]
            abs(frm) > bound || break
        end
    end
    dds = Float32(wk2_it)
    if dds > 5f0
        frm = 0f0
    elseif dds > 4f0
        frm = (dds - 4f0) * frm
    end
    oldrn[it] = frm
    return exp(frm)
end

"""
    species_sort!(state)

Fill the species-sorted tree index `scratch.idx1` (IND1) and the per-species
range table `scratch.sp_count_tab` (ISCT[sp,1:2]): trees grouped by species in
ascending tree-record order. This is the order growth/calibration consume the RNG,
so it must match FVS. (Pre-regen the within-species order is just ascending index.)
"""
function species_sort!(s::StandState)
    t = s.trees
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    @inbounds for sp in 1:MAXSP
        isct[sp, 1] = 0; isct[sp, 2] = 0
    end
    pos = 0
    @inbounds for sp in 1:MAXSP
        start = pos + 1
        for i in 1:t.n
            t.species[i] == sp || continue
            pos += 1
            ind1[pos] = Int32(i)
        end
        if pos >= start
            isct[sp, 1] = Int32(start); isct[sp, 2] = Int32(pos)
        end
    end
    return s
end
