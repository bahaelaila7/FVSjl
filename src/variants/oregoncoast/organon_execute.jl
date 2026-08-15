# =============================================================================
# organon_execute.jl — OC (Oregon Coast) ORGANON per-cycle orchestration + FVS copy-back (chunk C7).
#
# Ties the four validated ORGANON SWO growth components into one per-cycle call in the faithful
# `organon/execute2.f` EXECUTE + `organon/grow.f` GROW order, and ports the FVS-side copy-back seam
# that turns the ORGANON outputs into the FVS tree-record increments:
#   • grow.f growth-1 : DG_SWO           (chunk C3, organon_dg_swo)
#   • grow.f growth-2 : HG_SWO           (chunk C4, organon_hg_swo)
#   • grow.f MORTAL   : PM_SWO/MORTAL_RUN (chunk C6, organon_mortal_swo) — before crown, sets DEADEXP
#   • grow.f CROWGRO  : CROWGRO          (chunk C5, organon_cr_swo)      — consumes DEADEXP survivors
#   • oc/dgdriv.f:446-452 : DGRO → BARK → DIAGR → DDS → WK2  (this file: `oc_bratio` + `oc_organon_dds`)
#   • bin/FVSoc_buildDir/bratio.f : the OC (CA-family) variant bark ratio BRATIO (BARKB(5,29)/JBARK).
#
# GROW sequencing (faithful): DG and HG are computed on the START-of-cycle tree list; MORTALITY runs
# next (needs DGRO for the post-growth BA, produces DEADEXP and reduces EXPAN to survivors); CROWN
# runs last on the survivor stand. So the driver order is dg → hg → mort → cr (mortality BEFORE crown,
# exactly as C5's survivor-expansion finding requires).
#
# The FVS copy-back of DGRO is NOT `DBH += DGRO`: FVS loads `WK2 = ln(DIAGR·(2·DBH·BARK+DIAGR))`
# (DIAGR = DGRO·BARK) then the shared species loop forms the inside-bark increment
# `DG = √((DBH·BARK)² + exp(WK2)) − DBH·BARK` and updates `DBH += DG/BARK`. `oc_organon_dds` produces
# the WK2/DDS FVS loads; `oc_organon_dg` forms the inside-bark DG the StandState apply-loop consumes.
#
# MEASURED bit-exact vs the live FVSoc_clean oracle (scoped `DEBUG 1 / DGDRIV`, stand S248112 / ocmin,
# 17 valid ORGANON trees) — the dgdriv `DGRO,BARK,DIAGR,DDS` dump. See docs/OC_VARIANT_PORT_AUDIT.md.
# =============================================================================

# --- OC bark ratio (bin/FVSoc_buildDir/bratio.f — CA-family BARKB(5,29) + JBARK(50)) -----------
# FVS species index (1..50) → bark-coefficient row (1..29).
const OC_BRATIO_JBARK = Int[
    12, 5,12, 2, 2, 2, 1,13,13,10,
    27,10, 3,10, 4,25,25, 3, 3,26,
    12, 9,29,10, 3,19,18,17,20, 8,
    16,24,21, 6,16,14,23, 7,16,16,
    16,28,16,14,16,16,16,22,16,29]
# BARKB row 1..29: a (col 2), b (col 3), eqn type (col 5). FIA order 202,15,122,116,81,312,431,815,
# 93,108,98,242,263,351,1,818,807,805,801,811,839,981,361,821,117,127,101,631,211.
const OC_BRATIO_A = Float32[0.903563,-0.1593,-0.4448,-0.4448,-0.0549,0.08360,0.15565,-0.30722,0.9,0.9,0.958330,0.949670,0.933710,0.075256,0.933290,-0.26824,-0.17324,-0.19128,-0.75739,-0.78572,0.04817,-0.12791,-0.013484,-0.38289,-0.1429,0.93290,0.93290,-0.26824,0.70120]
const OC_BRATIO_B = Float32[0.989388,0.8911,0.8967,0.8967,0.8374,0.94782,0.90182,0.95956,0.0,0.0,1.0,1.0,1.0,0.94373,1.0,0.95767,0.94403,0.96147,0.93475,0.92472,0.92953,0.96579,0.98155,0.93545,0.8863,0.0,0.0,0.95354,1.04862]
const OC_BRATIO_TYPE = Int[1,2,2,2,2,2,2,2,3,3,1,1,1,2,1,2,2,2,2,2,2,2,2,2,2,3,3,2,1]

"""
    oc_bratio(is, d) -> BRATIO

`bin/FVSoc_buildDir/bratio.f` — the OC (CA-family) bark ratio for FVS species index `is` at DBH `d`.
Eqn 1: DIB=a·Dᵇ; 2: DIB=a+b·D; 3: BRATIO=a. Clamped to [0.80, 0.99]. (`H` is unused in this variant.)
"""
@inline function oc_bratio(is::Integer, d::Float32)
    d <= 0.0f0 && return 0.99f0
    r = OC_BRATIO_JBARK[is]
    a = OC_BRATIO_A[r]; b = OC_BRATIO_B[r]; ty = OC_BRATIO_TYPE[r]
    if ty == 1
        br = (a*fpow(d, b))/d
    elseif ty == 2
        br = (a + b*d)/d
    else
        br = a
    end
    br > 0.99f0 && (br = 0.99f0)
    br < 0.80f0 && (br = 0.80f0)
    return br
end

"""
    oc_organon_dds(is, dbh, dgro) -> (bark, diagr, dds)

`oc/dgdriv.f:446-452` — convert the ORGANON outside-bark diameter growth `dgro` into the FVS WK2
LN(DDS): `BARK=BRATIO`, `DIAGR=DGRO·BARK`, `DDS=ln(DIAGR·(2·DBH·BARK+DIAGR))`, floored at −9.21.
"""
@inline function oc_organon_dds(is::Integer, dbh::Float32, dgro::Float32)
    bark = oc_bratio(is, dbh)
    diagr = dgro*bark
    dds = flog(diagr*(2.0f0*dbh*bark + diagr))
    dds < -9.21f0 && (dds = -9.21f0)
    return bark, diagr, dds
end

"""
    oc_organon_dg(is, dbh, dds) -> DG (inside-bark 5-yr increment)

The shared FVS large-tree DDS→DG conversion for an ORGANON tree (`oc/dgdriv.f:536-557`, with
XDMULT=WK4=1, OLDRN=0, FRM=1): `D=DBH·BARK`, `DG=√(D²+exp(DDS))−D`. The StandState apply-loop then
grows the outside-bark DBH by `DG/BARK`. (DGBND capping is applied by the shared engine.)
"""
@inline function oc_organon_dg(is::Integer, dbh::Float32, dds::Float32)
    bark = oc_bratio(is, dbh)
    d = dbh*bark
    return sqrt(d*d + fexp(dds)) - d
end

"""
    OrganonGrowth

Per-cycle ORGANON SWO growth outputs for the buffer (records 1:n), all deterministic (DGSD=0):
`dgro`/`hgro`/`cr2`/`deadexp` (the C3/C4/C5/C6 per-tree results), `dds` (the FVS WK2 load for the
ORGANON trees), plus the stand quantities `spgrp`,`sba1`,`a1`,`a2` a caller/next cycle needs.
"""
struct OrganonGrowth
    dgro    ::Vector{Float32}
    hgro    ::Vector{Float32}
    cr2     ::Vector{Float32}
    deadexp ::Vector{Float32}
    dds     ::Vector{Float32}
    spgrp   ::Vector{Int32}
    sba1    ::Float32
    a1      ::Float32
    a2      ::Float32
end

"""
    organon_execute_swo(buf, isp_fvs; si_1, si_2, msdi_1, msdi_2, msdi_3, cyclg=0, mort=true)
        -> OrganonGrowth

The ORGANON SWO per-cycle EXECUTE/GROW orchestration (`organon/execute2.f` + `grow.f`, VERSION=1)
off the C1 `/ORGANON/` buffer. Runs, in the faithful GROW order: DG_SWO (C3) → HG_SWO (C4) →
MORTAL_RUN (C6) → CROWGRO (C5), then the `oc/dgdriv.f` DGRO→DDS copy-back. `isp_fvs[i]` is the FVS
species index (1..50) for the bark ratio. `si_1`/`si_2` are SITE_1−4.5 / SITE_2−4.5. Deterministic;
no RNG (ORGANON serial-corr suppressed, `oc/dgdriv.f:549`).
"""
function organon_execute_swo(buf::OrganonBuffer, isp_fvs::AbstractVector{<:Integer};
        si_1::Float32, si_2::Float32, msdi_1::Float32=0.0f0, msdi_2::Float32=0.0f0,
        msdi_3::Float32=0.0f0, cyclg::Int=0, mort::Bool=true,
        calib1::Union{Nothing,Vector{Float32}}=nothing)
    n = buf.ntrees
    dgro, spgrp, sba1, bal1, ball1, a1, a2 =
        organon_dg_swo(buf; si_1=si_1, msdi_1=msdi_1, msdi_2=msdi_2, msdi_3=msdi_3, cyclg=cyclg)
    hgro = organon_hg_swo(buf, dgro, spgrp; si_1=si_1, si_2=si_2, cyclg=cyclg, calib1=calib1)
    deadexp = organon_mortal_swo(buf, dgro, hgro, spgrp, bal1, ball1, a1, a2;
                                 si_1=si_1, cyclg=cyclg, mort=mort)
    cr2 = organon_cr_swo(buf, dgro, hgro, spgrp, deadexp; si_1=si_1, si_2=si_2, cyclg=cyclg)
    dds = zeros(Float32, n)
    @inbounds for i in 1:n
        buf.iorg[i] == 1 || continue                 # only valid ORGANON trees get the DDS copy-back
        _, _, dds[i] = oc_organon_dds(Int(isp_fvs[i]), buf.dbh1[i], dgro[i])
    end
    return OrganonGrowth(dgro, hgro, cr2, deadexp, dds, spgrp, sba1, a1, a2)
end
