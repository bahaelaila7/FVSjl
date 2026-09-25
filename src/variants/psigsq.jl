# dgdriv.f PSIGSQ(ISPC) per variant — see dg_psigsq in southern/diameter_growth.jl for the rationale.
# Included after every variant type exists.
dg_psigsq(::Union{Northeast,CentralStates,LakeStates,Ontario,CentralCalifornia,WestCascades,
                  PacificNorthwest,OregonCoast,Olympic}, sp::Int) = 0.0898f0
dg_psigsq(::CentralRockies, sp::Int) = 0.07f0
dg_psigsq(::Kootenai, sp::Int) = KT_PSIGSQ[sp]
dg_psigsq(::EasternMontana, sp::Int) = EM_PSIGSQ[sp]
dg_psigsq(::Teton, sp::Int) = TT_PSIGSQ[sp]
dg_psigsq(::Utah, sp::Int) = UT_PSIGSQ[sp]
dg_psigsq(::BlueMountains, sp::Int) = BM_PSIGSQ[sp]
dg_psigsq(::BritishColumbia, sp::Int) = BC_PSIGSQ[sp]
dg_psigsq(::CentralIdaho, sp::Int) = CI_PSIGSQ[sp]
dg_psigsq(::InlandEmpire, sp::Int) = IE_PSIGSQ[sp]
dg_psigsq(::Klamath, sp::Int) = NC_PSIGSQ[sp]
dg_psigsq(::WestSierra, sp::Int) = WS_PSIGSQ[sp]
dg_psigsq(::SouthCentralOregon, sp::Int) = SO_PSIGSQ[sp]
dg_psigsq(::EastCascades, sp::Int) = EC_PSIGSQ[sp]
