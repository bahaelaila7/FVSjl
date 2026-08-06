# UT establishment/planted-tree base height (ut/essubh.f). A FIXED per-species height table (no age/site
# regression, no EMSQR/DILATE) — the same shape as TT/CR. HHT is then clamped to [XMIN, HHTMAX] in the shared
# engine. Species order (ut/blkdat.f): 1=WB 2=LM 3=DF 4=WF 5=BS 6=AS 7=LP 8=ES 9=AF 10=PP 11=PI 12=WJ 13=GO
# 14=PM 15=RM 16=UJ 17=GB 18=NC 19=FC 20=MC 21=BI 22=BE 23=OS 24=OH.
# Values verbatim from ut/essubh.f CASE branches (CASE 11,14,17→0.5; 12,15,16→0.5; 18,19,22→10.0; 20,21→1.0).
const _UT_ESSUBH_HHT = Float32[1.0, 0.5, 2.0, 2.0, 1.0, 5.0, 3.0, 1.5, 0.75, 3.0,
                               0.5, 0.5, 5.0, 0.5, 0.5, 0.5, 0.5, 10.0, 10.0, 1.0,
                               1.0, 10.0, 1.0, 5.0]
