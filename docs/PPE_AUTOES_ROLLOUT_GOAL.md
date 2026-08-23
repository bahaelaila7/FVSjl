# ACTIVE GOAL — PPE port + AUTOES fix + ESTAB report-keywords (FVSjl)

USER 2026-08-23: "implement the recovered PPE, do that and all the other actionable items."
Branch: kt-variant-port. Doctrine UNCHANGED: bit-exact-or-cornered per chunk vs the LIVE relinked oracle;
MEASURE (g16 single-.o dump-replay, instrumented .sum byte-identical first); glibc libm ccall for Float32
transcendentals; never FFI an RNG; commit ONLY validated chunks; gate test/integration/test_multicycle.jl
339/11 byte-identical is the hard gate; block-on-signal via foreground/background Bash wait, NOT idle stop-hook cycling.

## STATUS 2026-08-23 (this session)
- ITEM 1 (PPE): ✅ **kernel phase COMPLETE**. 5 golden-anchorable PPBASE kernels bit-exact vs gfortran-16
  goldens — ppe_index_qsort! C11SRT/C26SRT/CH8SRT (4be76cba, 19/19+500 fuzz) + C26BSR/CH8BSR/SPBSRX binary
  search + HXINDX hex-neighbor + ADD1 (0eac7b03, 12,197 baked cases). Pure-kernel golden seam EXHAUSTED: the
  rest (alstd1/2, gp*, splaex, bygrps) are COMMON-mutating orchestration wrappers, not standalone functions —
  they belong to the landscape-wiring phase already covered by ppe_landscape.jl (reconstruction; no e2e oracle).
- ITEM 2 (AUTOES #143): ✅ **RESOLVED bit-exact** (081e0e6c). Two measured bugs — ESAVE body length hard-coded
  135 (ihab-10) instead of ihab-derived 127; kw_estab! forced IDSDAT=-9999 discarding esin.f's -1 ⇒ ingrowth
  fired 1 cyc late. Fixed ⇒ firing ic 1/3/5 + seeds 43303/61997/49053 BIT-FOR-BIT vs FVSie_g16; test_ie_estock
  128/128; gate 339/11. **Follow-on (af3aa7c9): the ~7-15% under-stocked TPA residual — long filed "cornered
  NSTORE regime" — was re-measured post-#143 and found to be THREE real establishment bugs (estab used STAND
  slope/aspect not per-plot PSLO/PASP [variant-general] + ESB calib skipped on cyc-1 ingrowth + stale pre-growth
  BA); fixed ⇒ 2029 TPA BIT-EXACT, rest ~1-2%, residual = 1.2% per-point BA attribution corner.**
- ITEM 3 (ESTAB keywords): ✅ **SURFACE CLOSED**. OUTPUT/HABGROUP inert (823b2a9f); MINPLOTS wired bit-exact
  (660edb9b, seed chains byte-exact at MINPLOTS 50/100/150); MECHPREP/BURNPREP wired LIVE bit-exact (a6fac7e0,
  h_base/h_mech100/h_burn100/h_mech50 FULL rows byte-identical vs FVSie_g16); PASSALL measured-inert; PLOTINFO
  inapplicable (DB pipeline). Only ADDTREES (external program) remains. Gate 339/11.

## ITEM 1 — PPE (Parallel Processing Extension) faithful port  [was wrongly "un-portable"; SOURCE RECOVERED]
The PPE landscape harness source was DELETED from the FVS Fortran repo (commit bc6e2377 "deleted PPE
directories", 2014; commons archived by 48d32e7b, 2018) — NOT absent. RECOVERED 152 files to
scratchpad/ppe/recovered/ppbase/ (ppmain.f 369 / alstd1.f / alstd2.f / splaex.f + gp*/hv*/sp* harness) at rev
bc6e2377^. PPMAIN = standalone PROGRAM: outer master-cycle driver over stands/point-groups that calls the
per-stand engine (GRSTND/GRINCR/HVREPS/GROEND) FVSjl already reimplements.
- Cross-stand-coupled parts (the only PPE-unique behaviors): ALSTD1/2 master cycle · hv* landscape harvest
  budgets · SPLAEX/SPNB spatial neighborhoods (inter-stand disturbance spread) · gp* point-group sampling.
- ORACLE: nothing in the current tree builds PPMAIN. Build FVSppe = compile scratchpad/ppe/recovered/ppbase/*.f
  (gfortran-16 -std=legacy -w -fno-automatic + isoc23 shim) linked with a base-FVS buildDir's *.o → standalone
  landscape executable = live oracle. If the full exe won't link, fall back to per-routine gfortran-16 goldens
  (WWPB-kernel recipe) for ALSTD1/2/SPLAEX/gp*/hv*.
- jl beachhead EXISTS: src/engine/ppe_landscape.jl + wwpb_landscape.jl (behavior-faithful area-weighted multi-
  stand aggregation). Upgrade reconstruction→faithful port: per-chunk bit-exact vs the FVSppe oracle, wired as a
  landscape orchestration layer above the single-stand engine. Additive/INERT (gate 339/11 must hold).

## ITEM 2 — AUTOES establishment-tally desync  [#143/item-1; the keystone lever]
jl does NOT model FVS's AUTOES natural-regen ingrowth tally (estab.f:651+ STOADJ block + per-species ESRANN
draws; documented in jl establishment.jl:227-228). ⇒ jl's establishment ESRANN stream is OFFSET vs the oracle
(measured: jl omits the cyc-1 ingrowth tally ⇒ seed 43303 vs oracle 61997). Consequences: the cornered
"jl est fires ~1 cycle late" seedling-timing residual AND it blocks MECHPREP-live + the model-effect ESTAB
keywords. Investigate estab.f tally scheduling (NTALLY / LAUTAL / NTALLY==99 ingrowth path → GO TO 44;
estab.f:223-259, 646-773) vs jl src/engine/establishment.jl. Port the missing AUTOES ingrowth-tally so the
ESRANN draw sequence matches the oracle bit-for-bit on the under-stocked IE fixture 753189105290487
(FVSie_estabdump oracle, /workspace/.iework/). Validate via per-plot ESDRAW/ITPP/NEWTPP dump-replay bit-exact.
HARD RISK: touching establishment RNG/scheduling can regress the gate — measure-first, back out if it can't hold
339/11. If it resolves: apply the staged MECHPREP kernel (scratchpad/estab/mechprep/staged/) + wire MINPLOTS/
PASSALL/PLOTINFO.

## ITEM 3 — ESTAB report/control keywords  [esin.f opt 6/10/14/18/20]
- OUTPUT (opt 6, IPRINT + JOREGT unit) + HABGROUP (opt 14, prints habitat-group table via ESMSGS) = PURE REPORT
  control, no .sum/RNG effect ⇒ parse into control flags, faithfully INERT in jl (jl emits no establishment .out
  report to gate). Small; wire in kw_estab!, gate 339/11 byte-identical, add a parse test.
- MINPLOTS (opt 20, MINREP) · PASSALL (opt 18, excess-tree passing, IBLK/CONFID) · PLOTINFO (opt 10, per-plot
  input) = MODEL-EFFECT on the establishment plot-replication/RNG ⇒ share the ITEM-2 AUTOES desync blocker; wire
  + validate ONLY after item 2 (else they ride the desynced stream = unvalidatable). MAXPLT is NOT a keyword
  (array-dim parameter).

## DONE / not-remaining (do not re-scope)
WRD/Western Root Disease = COMPLETE (src/engine/root_disease.jl 2683 ln, test_root_disease 1159/1159 confirmed
2026-08-23). NEWSPRED/BC spatial DM = COMPLETE (src/variants/britishcolumbia/newspred.jl + newspred_shd.jl +
mistletoe_report.jl; #196 closed-cornered). WSBWE GENDEFOL (0ffd6466), COVER CVBCAL (e6ea9dc2) = DONE this
session. ADDTREES = external program (un-portable, no in-tree source).

Off-switch (USER's call): touch docs/WESTERN_ROLLOUT_COMPLETE.
