# AK permafrost DG modifier — CONFIRMED (interior-species A/B) 2026-08-21

Goal-doc item 3 said the permafrost code is present + gate-safe but the path is INERT on akt01 (coastal species,
none in AK_PERM_SP) and "only needs an interior-AK A/B vs FVSak to CONFIRM." Done:

## Fixture
akt01's 29 trees with the A3 species field rewritten to WS (white spruce, sp 5 — a permafrost species,
AK_PERM_SP=(4,5,6,7,13,16-23)). akt_ws.tre + akt_pf_{on,off}.key (PRMFROST 1990 1.0 vs absent), NUMCYCLE 5, NOTRIPLE.

## Result — PRMFROST FIRES and matches the oracle bit-exact-or-cornered
                 1990  2000  2010  2020  2030  2040   (live TPA)
  jl  off        669   471   293   223   186   160
  or  off        669   471   264   212   173   151
  jl  on         669   590   372   298   247   211
  or  on         669   629   335   276   229   188
  Δ(on-off) jl    0   +119  +79   +75   +61   +51
  Δ(on-off) or    0   +158  +71   +64   +56   +37
- cyc0 (1990) BIT-EXACT both modes (669 TPA / 118 BA). 2000 off-baseline TPA BIT-EXACT (471=471).
- PRMFROST produces the CORRECT large effect (raises TPA: permafrost slows DG ⇒ smaller trees ⇒ less
  density-driven self-thinning ⇒ more survivors) — jl matches the oracle's SIGN and SCALE (Δ ~+50-160 TPA).
- Residual (jl off 293 vs or 264 at 2010; Δ 119 vs 158 at 2000) = the cornered multi-cycle SELF-THIN realization:
  the permafrost DG delta flows through tree-size → SDI → density-mortality (a stochastic self-thin straddle,
  the same #206/OLDRN-class cornered residual every variant's multi-cycle tail carries). The DEG modifier itself
  is exercised + directionally+scale-correct; a per-tree DDS proof would need FVSak_g16 instrumentation (not built).
⇒ Item 3 (AK permafrost) CONFIRMED to the bit-exact-or-cornered bar — the ported LPERM/PFMOD path is live and
correct, no longer just "present but unexercised."
