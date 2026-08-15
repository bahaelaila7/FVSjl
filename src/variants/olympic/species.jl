# =============================================================================
# species.jl — OP (Olympic) species map: FVS species (1..39) → ORGANON FIA code → NWO species group.
#
# Ported from op/orgspc.f (OSPMAP: FVS→ORGANON FIA code) composed with organon/execute2.f
# SPGROUP_RUN VERSION=2/3 (SCODE2: ORGANON FIA code → NWO group 1..11). MEASURED bit-exact against
# the live FVSop_dbg DGDRIV/DG_NWO dump ISPGRP column (stand S248112).
#
# FVS species order (op/orgspc.f, MAXSP=39):
#   1=SF  2=WF  3=GF  4=AF  5=RF  6=SS  7=NF  8=YC  9=IC 10=ES 11=LP 12=JP 13=SP 14=WP 15=PP
#  16=DF 17=RW 18=RC 19=WH 20=MH 21=BM 22=RA 23=MA 24=TO 25=GC 26=AS 27=CW 28=WO 29=J  30=LL
#  31=WB 32=KP 33=PY 34=DG 35=HT 36=CH 37=WI 38=(blank) 39=OT
# =============================================================================

# op/orgspc.f OSPMAP(39): FVS species sequence → valid ORGANON FIA species code.
const OP_OSPMAP = Int32[
    017,017,017,017,202,202,017,242,242,202,   # SF WF GF AF RF SS NF YC IC ES
    202,202,202,202,202,202,263,242,263,263,   # LP JP SP WP PP DF RW RC WH MH
    312,351,361,312,361,312,312,815,231,017,   # BM RA MA TO GC AS CW WO J  LL
    202,202,231,492,492,492,920,492,492]       # WB KP PY DG HT CH WI __ OT

# organon/execute2.f SPGROUP_RUN — SCODE2 (VERSION 2/3, NWO/SMC): ORGANON FIA code → group 1..11.
# 1=DF(202) 2=GF(017) 3=WH(263) 4=RC(242) 5=PY(231) 6=MD(361) 7=BL(312) 8=WO(815) 9=RA(351)
# 10=PD(492) 11=WI(920).  (SWO uses SCODE1/19 groups; RAP uses SCODE3/7 — not this variant.)
const OP_SCODE2 = Int32[202, 17, 263, 242, 231, 361, 312, 815, 351, 492, 920]

"""
    op_organon_fia(sp) -> Int32

op/orgspc.f ORGSPC — FVS species sequence number `sp` (1..39) → valid ORGANON FIA species code.
"""
@inline op_organon_fia(sp::Integer) = OP_OSPMAP[sp]

"""
    op_spgroup_nwo(sp) -> Int

FVS species `sp` (1..39) → NWO ORGANON species GROUP (1..11), i.e. ORGSPC∘SPGROUP_RUN(VERSION=2).
This is the `ISPGRP` that indexes the DG_NWO / HG_NWO / bark coefficient tables. Returns 0 if the
ORGANON FIA code is not an NWO species (should not occur — OSPMAP maps every OP species into SCODE2).
"""
@inline function op_spgroup_nwo(sp::Integer)
    fia = OP_OSPMAP[sp]
    @inbounds for g in 1:length(OP_SCODE2)
        OP_SCODE2[g] == fia && return g
    end
    return 0
end
