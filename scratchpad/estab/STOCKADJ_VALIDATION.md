# STOCKADJ (ESTAB packet keyword) — bit-exact vs FVSie_clean

esnutr.f IACTK 440: STOADJ = PRMS(1). Applied estab.f:578-580 — clamp STOADJ≥0.001, then
PROB1 = logistic(PN+ESB-ESB1)·STOADJ (the AUTOES/ESTOCK re-stocking tally). STOADJ<0.0001
(STOCKADJ 0.0 / NATURAL) takes a SEPARATE estab.f branch (GO TO 137/229/163) reached via
jl's scheduled-NATURAL path — OUT OF SCOPE for this multiply chunk.

Repro: extract_ie.jl builds ie_stockadj.db (stand 12343703010690, 1 stand + 33 trees from
the FIADB). ie_sa_base.key = ESTAB/END (STOADJ=1 default); ie_sa_050.key = ESTAB/STOCKADJ 0.5/END.
Both have THINPRSC 2029 (post-thin LAUTAL AUTOES re-stocking = where STOADJ bites).

Oracle (FVSie_clean, keyfile name on stdin) vs jl (run_keyfile IE) — Δ(sa050−base) live TPA:
  2007 +0/+0 · 2017 +0/+0 · 2027 −88/−88 · 2037 −75/−75 · 2047 −205/−205 · 2057 −172/−172
jl baseline TPA also bit-identical to oracle baseline. BIT-EXACT.

Gate: multicycle 339/11 byte-identical (inert when keyword absent).
