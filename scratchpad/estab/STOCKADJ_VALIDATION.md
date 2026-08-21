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

## TALLY / TALLYONE / TALLYTWO (esin.f 16/11/12) — bit-exact vs FVSie_clean
esnutr.f:163-252: schedule a user establishment tally at a date; NTALLY=IACTK-427 (TALLY/TALLYONE=1, TALLYTWO=2
only with a prior TALLYONE else 1). Wired: kw_estab! pushes ScheduledActivity(date,427/428/429); ie_autoes_establish!
honors a scheduled tally due in [year,next_year) within the KDT+1-IDSDAT≤20 window, overriding the automatic rules.
The ESTAB-END 427 at inv-20 is stale (>20yr) → correctly dropped, no regression.
Repro ie_ty_base.key (ESTAB/END) vs ie_ty_t20.key (ESTAB/TALLY 2020/END) on stand 12343703010690. Δ(TALLY 2020−base)
live TPA: 2027 −88, 2037 −75, 2047 −86, 2057 −198 — jl matches oracle EVERY cycle; jl baseline bit-identical. Gate 339/11.
