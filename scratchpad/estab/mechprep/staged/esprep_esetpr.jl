# =============================================================================
# ESTAB site-preparation (MECHPREP / BURNPREP) — faithful transcription
# =============================================================================
# Port of estb/esprep.f (default site-prep probabilities), estb/esetpr.f
# (MECHPREP/BURNPREP keyword → PNONE/PMECH/PBURN/IALN + ZMECH/ZBURN), and the
# estb/estab.f:373-399 "SAMPLE WITHOUT REPLACEMENT" per-plot IPPREP assignment
# from the WK6 site-prep RNG vector.
#
# STATUS: transcription VALIDATED bit-exact at the KERNEL level (ie_esetpr_sample
# reproduces the live FVSie_estabdump IPPREP given the oracle's exact SUMUP+WK6 —
# see test_mechprep.jl / golden_mechprep.jl). NOT wired into the live engine: the
# live end-to-end validation is blocked (see VALIDATION.md — jl's disturbance-tally
# seed stream is desynced from the oracle on the harvest scenario that exercises
# site prep, and PROB1 is prep-invariant for the only available IE fixture).
#
# Fortran field meanings (esprin.f / esetpr.f, measured on FVSie_estabdump):
#   MECHPREP <date> <pct> : ARRAY(1)=date/cycle (ZMECH), ARRAY(2)=% of plots to
#                           mechanically prep (0-100). Sets PMECH=pct/100, IALN(2)=1.
#   BURNPREP <date> <pct> : same, sets PBURN=pct/100, IALN(3)=1, ZBURN=date.
#   pct blank ⇒ PMECH/PBURN=0 (no plots prepped). Activity codes 493(MECH)/491(BURN).
# =============================================================================

"""
    ie_esprep(iser, aspect, slope, ba, elev) -> (pnone, pmech, pburn)

esprep.f — DEFAULT site-prep probabilities when the user supplied NO MECHPREP/
BURNPREP keyword (estab.f:246→ESPREP). `iser` = the ESTOCK habitat series index
(MYTYPE(ITYPE)); `aspect` = ASPECT (radians); `slope` = SLOPE (fraction);
`ba` = stand BA; `elev` = ELEV (100s ft). XPREP(3,ISER) is the per-series
calibration bias (esblkd.f DATA — MUST be ported from the FVSie block data for
this default path to be exact; the keyword-driven path below does NOT use it).
"""
function ie_esprep(iser::Integer, aspect::Real, slope::Real, ba::Real, elev::Real,
                   xprep::AbstractMatrix)
    ca = cos(Float32(aspect)); sa = sin(Float32(aspect)); sl = Float32(slope)
    b  = Float32(ba); el = Float32(elev)
    xp(k) = Float32(xprep[k, iser])
    pn = 1.043151f0 + xp(1) - 0.220954f0*ca*sl + 0.369575f0*sa*sl + 0.769112f0*sl +
         0.260178f0*log(b + 1f0) - 0.029689f0*el
    pnone = 1f0/(1f0 + exp(-pn))
    pn = -1.852031f0 + xp(2) + 0.492668f0*ca*sl + 0.192020f0*sa*sl - 0.966674f0*sl -
         0.085920f0*log(b + 1f0) + 0.024939f0*el
    pmech = 1f0/(1f0 + exp(-pn))
    pn = -15.195303f0 + xp(3) + 0.0519477f0*ca*sl - 0.6135848f0*sa*sl - 0.0890163f0*sl -
         0.377915f0*log(b + 1f0) + 0.5303707f0*el - 0.0049081f0*el*el
    pburn = 1f0/(1f0 + exp(-pn))
    return (pnone, pmech, pburn)
end

"""
    ie_esetpr(pmech_pct, pburn_pct; zmech, zburn) -> NamedTuple

esetpr.f — resolve the MECHPREP/BURNPREP keyword values into the site-prep
mixture. `pmech_pct`/`pburn_pct` = the keyword %-of-plots (0-100, or `nothing`
if that keyword absent). Returns PMECH/PBURN fractions, the IALN flags, and
whether any keyword was supplied. (METH — the "prep all plots" NPRMS==0 branch —
is inert here because esprin.f always stores 1 param, so NPRMS is always ≥1;
verified on FVSie_estabdump: a bare `MECHPREP <date>` gives PMECH=0.)
"""
function ie_esetpr(pmech_pct, pburn_pct; zmech::Real, zburn::Real)
    pmech = 0f0; pburn = 0f0; ialn2 = 0; ialn3 = 0
    if pburn_pct !== nothing
        pburn = Float32(pburn_pct) / 100f0; ialn3 = 1
    end
    if pmech_pct !== nothing
        pmech = Float32(pmech_pct) / 100f0; ialn2 = 1
    end
    return (pmech = pmech, pburn = pburn, ialn2 = ialn2, ialn3 = ialn3,
            zmech = Float32(zmech), zburn = Float32(zburn),
            any_kw = (ialn2 == 1 || ialn3 == 1))
end

"""
    ie_esetpr_normalize(pnone, pmech, pburn, ialn2, ialn3) -> (s1, s2, s3)

estab.f:246-378 — turn (PNONE,PMECH,PBURN) into the normalized cumulative-bucket
weights SUMUP(1..3) used by the sampler. With a user keyword (IALN set) the
mixture is (1-PMECH-PBURN, PMECH, PBURN) [renormalized if PMECH+PBURN>1];
without a keyword the ESPREP triple is used. Then SUMUP = each / Σ.
"""
function ie_esetpr_normalize(pnone::Real, pmech::Real, pburn::Real, ialn2::Integer, ialn3::Integer)
    pm = Float32(pmech); pb = Float32(pburn); pn = Float32(pnone)
    if ialn2 == 1 || ialn3 == 1                     # estab.f:249 user keyword path
        s = pm + pb
        if s > 1f0; pm /= s; pb /= s; end           # estab.f:353-356
        pn = 1f0 - pm - pb                          # estab.f:363
    end
    tot = pn + pm + pb                              # estab.f:375
    return (pn/tot, pm/tot, pb/tot)
end

"""
    ie_esetpr_sample(sumup, wk6, nptids, idup) -> Vector{Int}

estab.f:373-399 "SECTION TO CHOOSE SITE PREPS. SAMPLE WITHOUT REPLACEMENT."
Given the normalized bucket weights `sumup` (=SUMUP(1..3)) and the WK6 site-prep
RNG vector (`wk6`, length nptids·idup, already drawn on the :estab stream), assign
each of the DUPNPT=nptids·idup replicated plots an IPPREP ∈ {1=NONE,2=MECH,3=BURN}.

BIT-EXACT vs FVSie_estabdump given identical (sumup, wk6) — see test_mechprep.jl.
"""
function ie_esetpr_sample(sumup, wk6::AbstractVector, nptids::Integer, idup::Integer)
    dupnpt = Float32(nptids * idup)
    su = Float32[sumup[1], sumup[2], sumup[3]]
    ipprep = Vector{Int}(undef, nptids * idup)
    n = 0
    @inbounds for _ii in 1:idup, _nn in 1:nptids
        n += 1
        draw = Float32(wk6[n]) * (((dupnpt + 1f0) - Float32(n)) / dupnpt)   # estab.f:388
        s = 0f0; sel = 0
        for i in 1:2                                                        # estab.f:389-393
            s += su[i]
            if draw > s; continue; end
            sel = i; break
        end
        if sel == 0                                                        # estab.f:394-395
            sel = 3
            su[3] < 0f0 && (sel = 1)
        end
        ipprep[n] = sel                                                    # estab.f:396
        su[sel] -= 1f0 / dupnpt                                            # estab.f:397 (without replacement)
        su[sel] < 0f0 && (su[sel] = 0f0)
    end
    return ipprep
end
