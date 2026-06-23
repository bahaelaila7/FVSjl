# =============================================================================
# fire/biomass.jl — Jenkins tree biomass (FFE chunk F1)
#
# Ported from: bin/FVSsn_buildDir/fmcbio.f (FMCBIO) + the BIOGRP table (fmcblk.f)
# + the METRIC.F77 unit constants.
#
# The Jenkins (2003) national-scale biomass equations: per-tree aboveground,
# merchantable, and belowground (root) biomass from DBH and a species biomass
# group. This is the foundation of the FFE carbon pools (FMCADD) and a building
# block for crown/surface fuels; it fills the (otherwise-zero) per-tree biomass
# fields the carbon model consumes.
#
# All biomass is returned in (imperial) tons, matching FMCBIO's KGtoTI scaling.
# =============================================================================

# METRIC.F77 unit conversions used by FMCBIO.
const _IN_TO_CM = 2.54f0
const _KG_TO_TI = 1.102311f0 / 1000f0          # TMtoTI / 1000 (kg → imperial tons)

# Jenkins aboveground coefficients by biomass group 1..10 (fmcbio.f B0A/B1A).
const _JENKINS_B0A = Float32[-2.0336, -2.2304, -2.5384, -2.5356, -2.0773,
                             -2.2094, -1.9123, -2.4800, -2.0127, -0.7152]
const _JENKINS_B1A = Float32[ 2.2592,  2.4435,  2.4814,  2.4349,  2.3323,
                              2.3867,  2.3651,  2.4835,  2.4342,  1.7029]
# Merch + belowground coefficients, indexed [softwood, hardwood] (JGRP) (fmcbio.f).
const _JENKINS_B0M = Float32[-0.3737, -0.3065]
const _JENKINS_B1M = Float32[-1.8055, -5.4240]
const _JENKINS_B0B = Float32[-1.5619, -1.6911]
const _JENKINS_B1B = Float32[ 0.6614,  0.8160]

"""
    jenkins_biomass(coef, sp, dbh) -> (above, merch, root)

Per-tree Jenkins biomass in tons (FMCBIO, fmcbio.f) for species `sp` and DBH `dbh`
(inches). Returns aboveground, merchantable, and belowground (root) biomass.

The species' Jenkins group (`BIOGRP`, `fire_biomass.csv`) selects the aboveground
coefficients; groups 1–5 are softwood and 6–10 hardwood for the merch/root forms.
Equations were fit for trees ≥ 2.5 cm DBH; smaller trees use the 2.5-cm aboveground
biomass scaled linearly by DBH. Merch biomass is non-zero only at/above the species'
merch DBH limit (`DBHMIN`).
"""
@inline function jenkins_biomass(coef::SpeciesCoefficients, sp::Integer, dbh::Float32)
    dcm = dbh * _IN_TO_CM
    dcm > 0f0 || return (0f0, 0f0, 0f0)
    igrp = Int(coef_col(coef, :bio_group)[sp])
    jgrp = igrp > 5 ? 2 : 1                          # softwood (1) vs hardwood (2)
    b0a = _JENKINS_B0A[igrp]; b1a = _JENKINS_B1A[igrp]
    if dcm >= 2.5f0
        above = exp(b0a + b1a * log(dcm))
        root  = above * exp(_JENKINS_B0B[jgrp] + _JENKINS_B1B[jgrp] / dcm)
    else                                             # < 2.5 cm: 2.5-cm value scaled by DBH
        above = exp(b0a + b1a * log(2.5f0)) * (dcm / 2.5f0)
        root  = above * exp(_JENKINS_B0B[jgrp] + _JENKINS_B1B[jgrp] / 2.5f0)
    end
    merch = dbh >= coef_col(coef, :dbh_min)[sp] ?
            above * exp(_JENKINS_B0M[jgrp] + _JENKINS_B1M[jgrp] / dcm) : 0f0
    return (above * _KG_TO_TI, merch * _KG_TO_TI, root * _KG_TO_TI)
end
