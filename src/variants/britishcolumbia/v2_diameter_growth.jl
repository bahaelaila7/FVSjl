# =============================================================================
# v2_diameter_growth.jl (britishcolumbia) — BC V2 (LV2ATV) large-tree DDS formula. CHUNK 2b (V2 regime).
#
# The V2 regime (dgf.f LV2ATV=true) applies to NON-ICH/IDF/SBS/SBPS zones (ESSF/MS/PP). Its large-tree DDS
# is the IE-form Wykoff (dgf.f:1940 + 1992-1995): a linear predictor → log. FORMULA + the 5 per-species
# coefficients VALIDATED bit-exact (99.1%, ≤1 ULP) vs instrumented FVSbc on all_BC_essf (ESSFdk/01),
# GIVEN the oracle's CONSPP + DGDSQ. ⚠ The CONSPP habitat-class RESOLUTION (DGHAB/DGFOR/DGDS via MAPHAB/
# MAPLOC/MAPDSQ + SEILTDG site params) is task #133's next step — until then this can't run in-engine.
# =============================================================================

const BC_V2_DGLD   = Float32[0.56445,0.54140,0.56888,0.68810,0.68712,0.58705,0.89503,0.73045,0.86240,0.66101,0.89778,0.89778,0.89778,0.56888,0.89778]
const BC_V2_DGCR   = Float32[1.08338,1.03478,2.06850,1.93969,1.64133,1.29360,1.85558,1.54643,0.52044,1.31618,1.28403,1.28403,1.28403,2.06850,1.28403]
const BC_V2_DGCRSQ = Float32[0.0,0.07509,-0.62361,-0.78258,-0.27244,0.0,-0.36393,-0.26635,0.86236,0.0,0.0,0.0,0.0,-0.62361,0.0]
const BC_V2_DGBAL  = Float32[0.42112,0.43637,0.50202,0.45142,0.0,0.74596,-0.03665,0.25639,0.0,0.0,0.0,0.0,0.0,0.50202,0.0]
const BC_V2_DGDBAL = Float32[-2.08272,-2.03256,-2.11590,-1.76812,-0.80918,-2.28375,-0.43329,-1.18218,-0.51270,-1.25881,-0.66110,-0.66110,-0.66110,-2.11590,-0.66110]

"""
    bc_v2_dds(sp, d_in, bal, cr, conspp, dgdsq) -> ln(DDS) (WK2)

BC V2 large-tree DDS (dgf.f:1940 + 1992-1995), IE-form: `DDS = CONSPP + DGLD·ln(D) + DGBAL·BAL +
CR·(DGCR + CR·DGCRSQ) + DGDSQ·D² + DGDBAL·BAL/ln(D+1)`; then `max(-9.21, log(max(0.001, DDS)))`.
`d_in` inches; `conspp` (DGCON+COR+CCF terms) + `dgdsq` (habitat-class) come from the CONSPP resolution
(bc_v2_dgcons!, task #133 — not yet ported). VALIDATED bit-exact (99.1%, ≤1 ULP) on all_BC_essf.
"""
@inline function bc_v2_dds(sp::Integer, d_in::Real, bal::Real, cr::Real, conspp::Real, dgdsq::Real)
    d = Float32(d_in); ald = log(d)
    dds = Float32(conspp) + BC_V2_DGLD[sp]*ald + BC_V2_DGBAL[sp]*Float32(bal) +
          Float32(cr)*(BC_V2_DGCR[sp] + Float32(cr)*BC_V2_DGCRSQ[sp]) +
          Float32(dgdsq)*d*d + BC_V2_DGDBAL[sp]*Float32(bal)/log(d + 1f0)
    return max(-9.21f0, log(max(0.001f0, dds)))
end

# -----------------------------------------------------------------------------
# V2 CONSPP resolution (dgf.f ENTRY DGCONS, LV2ATV NI-branch 2032-2075). Per-species DGCON via habitat/
# forest maps, then per-tree CONSPP = DGCON + COR + 0.01·DGCCF·RELDEN (SEICCF=0 for the NI branch, which
# is the LSPPOK=false path — i.e. no matched SEI subzone). ITYPE/IFOR = grinit defaults 4/4 (grinit.f:201,207):
# becset's ITYPE override doesn't fire for the DEFAULT-renamed ESSF/MS/PP stands, so ITYPE stays at 4.
# NI-branch DGCON VALIDATED bit-exact (sp14/OC on ESSFdk/01, DGCON=0.7631423). ⚠ SEI-branch (LSPPOK=true,
# JBECY>0 matched subzone) NOT ported — dgf! errors loudly if a present species would need it.
# -----------------------------------------------------------------------------

# CSV loaders (rows = species 1..15; cols = class values). Jagged Vector{Vector} indexed [sp][class].
_bc_v2_mat(name)  = [Float32[parse(Float32, x) for x in r[2:end]] for r in _bc_readcsv("v2dg_$name.csv")[2:end]]
_bc_v2_imat(name) = [Int[parse(Int, x)       for x in r[2:end]] for r in _bc_readcsv("v2dg_$name.csv")[2:end]]
_bc_v2_vec(name)  = Float32[parse(Float32, r[2]) for r in _bc_readcsv("v2dg_$name.csv")[2:end]]

const BC_V2_DGHAB  = _bc_v2_mat("DGHAB")     # [sp][ISPHAB]  habitat-class DG constant
const BC_V2_DGFOR  = _bc_v2_mat("DGFOR")     # [sp][ISPFOR]  forest-location DG constant
const BC_V2_DGDS   = _bc_v2_mat("DGDS")      # [sp][ISPDSQ]  D² coefficient by class
const BC_V2_DGCCFA = _bc_v2_mat("DGCCFA")    # [sp][ISPCCF]  CCF coefficient by class
const BC_V2_OBSERV = _bc_v2_mat("OBSERV")    # [sp][ISPHAB]  DGSCOR shrinkage prior (ATTEN)
const BC_V2_MAPHAB = _bc_v2_imat("MAPHAB")   # [sp][ITYPE] → ISPHAB
const BC_V2_MAPLOC = _bc_v2_imat("MAPLOC")   # [sp][IFOR]  → ISPFOR
const BC_V2_MAPDSQ = _bc_v2_imat("MAPDSQ")   # [sp][IFOR]  → ISPDSQ
const BC_V2_MAPCCF = _bc_v2_imat("MAPCCF")   # [sp][ITYPE] → ISPCCF
const BC_V2_DGEL   = _bc_v2_vec("DGEL")
const BC_V2_DGEL2  = _bc_v2_vec("DGEL2")
const BC_V2_DGSASP = _bc_v2_vec("DGSASP")
const BC_V2_DGCASP = _bc_v2_vec("DGCASP")
const BC_V2_DGSLOP = _bc_v2_vec("DGSLOP")
const BC_V2_DGSLSQ = _bc_v2_vec("DGSLSQ")

const BC_V2_ITYPE = 4    # grinit.f:207 default (becset override inert for DEFAULT-renamed V2 stands)
const BC_V2_IFOR  = 4    # grinit.f:201 default

# bc_lv2atv (V2/V3 regime predicate) is defined in mortality.jl.

"""BC V2 DGCONS — per-species NI-branch DGCON (dgf.f:2061-2074). Stores DGCON→dg_const, OBSERV→atten."""
function bc_v2_dgcons!(s::StandState)
    c = s.calib; ctl = s.control; p = s.plot
    elev = Float32(p.elevation); asp = Float32(p.aspect); slope = Float32(p.slope)
    it = BC_V2_ITYPE; ifr = BC_V2_IFOR
    @inbounds for sp in 1:nspecies(BritishColumbia())
        isphab = BC_V2_MAPHAB[sp][it]
        ispfor = BC_V2_MAPLOC[sp][ifr]
        dgcon = BC_V2_DGHAB[sp][isphab] + BC_V2_DGFOR[sp][ispfor] +
                BC_V2_DGEL[sp]*elev + BC_V2_DGEL2[sp]*elev*elev +
                (BC_V2_DGSASP[sp]*sin(asp) + BC_V2_DGCASP[sp]*cos(asp) + BC_V2_DGSLOP[sp])*slope +
                BC_V2_DGSLSQ[sp]*slope*slope
        (ctl.dg_cor2_on && ctl.dg_cor2[sp] > 0f0) && (dgcon += log(ctl.dg_cor2[sp]))
        c.dg_const[sp] = dgcon
        c.atten[sp]    = BC_V2_OBSERV[sp][isphab]      # ATTEN = OBSERV(ISPHAB, ISPC) (dgf.f:2074)
    end
    return s
end

"""BC V2 `dgf!` body — per-tree DDS via the IE-form (imperial: D inches, BAL=(1−PCT/100)·BA/100)."""
function bc_v2_dgf!(s::StandState)
    p, t, c = s.plot, s.trees, s.calib
    wk2 = view(s.scratch.wk, 2, :)
    relden = p.relative_density; ba100 = p.basal_area / 100f0
    nsp = nspecies(BritishColumbia())
    it = BC_V2_ITYPE; ifr = BC_V2_IFOR
    dgccf = zeros(Float32, nsp); dgdsq = zeros(Float32, nsp)
    @inbounds for sp in 1:nsp
        dgccf[sp] = BC_V2_DGCCFA[sp][BC_V2_MAPCCF[sp][it]]
        dgdsq[sp] = BC_V2_DGDS[sp][BC_V2_MAPDSQ[sp][ifr]]
    end
    @inbounds for i in 1:t.n
        d = t.dbh[i]
        d <= 0f0 && continue
        sp = Int(t.species[i])
        (sp < 1 || sp > nsp) && continue
        conspp = c.dg_const[sp] + c.dg_cor[sp] + 0.01f0 * dgccf[sp] * relden   # SEICCF=0 (NI branch)
        pct = t.crown_ratio[i]                                                 # FVS PCT (trees.jl:51)
        bal = (1f0 - pct / 100f0) * ba100
        cr  = Float32(t.crown_pct[i]) * 0.01f0
        wk2[i] = bc_v2_dds(sp, d, bal, cr, conspp, dgdsq[sp])
    end
    return s
end
