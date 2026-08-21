# DBS extension-output-table emission recipe (oracle A/B unblock) — 2026-08-21

Getting FVSie_clean to EMIT the missing DBS extension tables (item 2: FVS_SnagDet, FVS_StrClass,
FVS_CanProfile, FVS_Climate, FVS_BM_*, ...) for a per-table oracle A/B had stalled THREE times
(FVS_TreeList, FVS_SnagDet x2) on finicky keyfile plumbing. CRACKED for FVS_SnagDet — the recipe:

## The working keyfile structure (validated: emits FVS_SnagDet with 66 real rows)
```
STDIDENT
<CN>
DATABASE
DSNout
<out.db>
DSNin
<in.db>
StandSQL
SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
TreeSQL
SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
END
NUMCYCLE         5
FMIN
SNAGSUM
SNAGOUT
END
DATABASE
SUMMARY
SNAGSUDB
SNAGOUDB
END
PROCESS
STOP
```

## The three keys that were each individually necessary
1. **Two DATABASE blocks**: block 1 carries DSNout + DSNin + the input SQL together; block 2 carries the
   output-table TOGGLES (SUMMARY / SNAGSUDB / SNAGOUDB). Splitting DSNout into its own leading block, or
   putting the toggles with the SQL, gives exit 10/20 or an empty/partial DB. (Same structure that first
   got FVS_Summary to write during the TreeLiDB attempts.)
2. **Toggle values must be BLANK (default 1), not 2**: `SNAGOUDB 2` sets ISDET=2 = REDIRECT mode → renames
   the table (we saw `FVS_Summary2`) / suppresses it. Blank = default = the real table name.
3. **The FFE report keyword that populates the arrays** (the crux):
   - `FMIN ... END` alone activates FFE but writes NEITHER snag table.
   - `SNAGSUM` (FFE keyword, distinct from the DBS `SNAGSUDB`) turns on snag *summary* tracking → FVS_SnagSum.
   - `SNAGOUT` (FFE keyword, fmin.f opt 12: start date/cycle field1=default 1, years field2=default 200)
     sets ISNAGB/ISNAGE — the year window FMSOUT (fmsout.f:90-91) gates on. WITHOUT it FMSOUT RETURNs and
     DBSFMDSNAG is never called, so FVS_SnagDet is never created. WITH it → 66 rows.

## Mapping (dbsin.f DBS option → table → FFE report keyword needed)
- SNAGSUDB (opt 22, ISSUM) → FVS_SnagSum  ← needs FFE `SNAGSUM`
- SNAGOUDB (opt 23, ISDET) → FVS_SnagDet  ← needs FFE `SNAGOUT` (fmsout.f/dbsfmdsnag.f)
- STRCLSDB (opt 24, ISTRCLAS) → FVS_StrClass ← needs the FFE StrClass report keyword (analogous gate — TBD)
- CANProf → FVS_CanProfile (fmpocr.f:245 DBSFMCANPR) ← needs the canopy-profile report
- Each extension table has this pattern: a DBS toggle (enables the SQLite write) AND its extension REPORT
  keyword (populates the arrays + sets the year-window gate). BOTH are required.

## Oracle A/B target now available for FVS_SnagDet
FVS_SnagDet columns: CaseID, StandID, Year, SpeciesFVS/PLANTS/FIA, DBH_Class(1-6), Death_DBH,
Current_Ht_Hard/Soft, Current_Vol_Hard/Soft, Total_Volume, Year_Died, Density_Hard/Soft/Total.
Sample oracle rows (stand 753189105290487, cyc 2029): GF/DBH_Class1 Density_Hard 161.56; RC/1 27.75; ...
jl has the source data (SnagList: sp, dbh, den_hard/den_soft — state.jl:677). Port = aggregate SnagList by
(species, DBH_Class, Year_Died) + serialize via a new write_dbs_snagdet! in io/dbs_output.jl, wire SNAGOUDB.

## CAVEAT for validating the port
FVS_SnagDet bit-exactness requires the underlying FFE snag MODEL (densities/heights/volumes/decay) to be
bit-exact on the test stand. IE FFE is NOT in the validated-FFE set (OC/CA/EC/WS/SO/WC/NC/eastern are) — so
validate FVS_SnagDet on a variant whose FFE snags are already bit-exact (e.g. OC or CA), NOT this IE stand.
The IE stand here only proves the ORACLE-EMISSION recipe; the jl A/B belongs on a validated-FFE variant.

## FVS_SnagDet port — DE-RISKED 2026-08-21 (all pieces in hand, pure serialization)
Confirmed jl's SnagList (state.jl:677) carries EVERY field the table needs — it is a pure serialization port, NOT
a model extension:
  - sp                        → SpeciesFVS/PLANTS/FIA (jl species-code maps)
  - dbh                       → DBH_Class (1-6 binning) + Death_DBH
  - den_hard / den_soft       → Density_Hard / Density_Soft / Density_Total
  - yrdead                    → Year_Died (+ Year = the report cycle)
  - htcur (current snag ht)   → Current_Ht_Hard / Current_Ht_Soft
  - bolevol / fallvol         → Current_Vol_Hard / Current_Vol_Soft / Total_Volume
Write pattern MIRRORS the existing write_dbs_snagsum! (dbs_output.jl:259) — prepare INSERT, loop rows, execute.
EXECUTABLE next-session plan:
  1. Read dbsfmdsnag.f:145-210 (the DO JYR/IDC/JCL body) to match the EXACT aggregation: DBH-class thresholds,
     density-weighting of Current_Ht/Vol, how Death_DBH is a per-class aggregate, the SpeciesFVS/PLANTS/FIA codes.
  2. Add _FVS_SNAGDET_CREATE (17-col schema, dbsfmdsnag.f:104-121) + write_dbs_snagdet!(dbpath,caseid,standid,rows)
     aggregating SnagList by (species, DBH_Class, yrdead).
  3. Wire the SNAGOUDB DBS toggle (keyword_dispatch.jl DATABASE block) → call write_dbs_snagdet! at output.
  4. A/B on a VALIDATED-FFE variant (OC or CA — IE FFE is not in the validated set): run the oracle with the
     recipe above (FMIN/SNAGSUM/SNAGOUT + DBS SNAGOUDB) and jl with the same, diff FVS_SnagDet column-by-column.
  5. Commit ONLY if bit-exact-or-cornered. Same shape for FVS_StrClass (structure_stage.jl has the stage; the
     per-stratum detail — Stratum_N_DBH/Ht/Crown/species — needs jl to expose the strata, a bigger lift).

## FVS_SnagDet — COMPLETE aggregation spec extracted 2026-08-21 (fmsout.f) — port now fully de-risked
FMSOUT builds the per-(species IDC, deathyear-index JYR, dbhclass JCL 1-6) arrays that DBSFMDSNAG binds:
  - JYR = IYR − YRDEAD(II) + 1 (years since death, clamp 1..100).  [jl: IYR − SnagList.yrdead[i] + 1]
  - JCL: DO JCL=1,5 { if DBHS(II) < SNPRCL(JCL+1) → that JCL } else JCL=6.  [SNPRCL = 6 snag-print-class DBH
    breakpoints — LOOK UP in fmcom/blkdat; jl bins SnagList.dbh[i] the same way]
  - Density: TOTDS += DENIS (soft) always; if HARD(II) → TOTDH += DENIH (hard) else TOTDS += DENIH.
    [jl: DENIH=SnagList.den_hard, DENIS=SnagList.den_soft; HARD flag = jl's hard/soft state via yrdead+DKTIME]
  - Height (density-weighted): TOTHTH += HTIH·DENIH (hard), TOTHTS += HTIS·DENIS (soft; +HTIH·DENIH if !HARD).
    [jl: HTIH/HTIS = SnagList.htcur]
  - DBH (density-weighted): TOTDBH += DBHS·(DENIS+DENIH).
  - Volume (summed): TOTVLH += SNVOLH, TOTVLS += SNVOLS.  [SNVOLH/SNVOLS = the snag CUBIC volume computed just
    above line 172 — LOOK UP; NOT necessarily jl's SnagList.bolevol (tons biomass) — this is the OPEN risk: confirm
    the report volume definition matches a jl SnagList field or needs computing]
  Then NORMALIZE: TOTN=TOTDH+TOTDS; TOTDBH/=TOTN; TOTHTH/=TOTDH; TOTHTS/=TOTDS (guard >0). Emit rows where TOTN>0.
  DBSFMDSNAG binds: Year=IYR, Species{FVS/PLANTS/FIA}=JSP/PLNJSP/FIAJSP(IDC), DBH_Class=JCL, Death_DBH=TOTDBH,
  Current_Ht_Hard=TOTHTH, Current_Ht_Soft=TOTHTS, Current_Vol_Hard=INT(TOTVLH), Current_Vol_Soft=INT(TOTVLS),
  Total_Volume=INT(TOTVLH+TOTVLS), Year_Died=IYR−JYR+1, Density_Hard=TOTDH, Density_Soft=TOTDS, Density_Total=TOTN.
REMAINING OPEN (2 lookups + 1 fixture): (a) SNPRCL breakpoints, (b) SNVOLH/SNVOLS volume definition vs jl SnagList
fields, (c) an FFE-bit-exact fixture (AK growth/vol/mort confirmed cyc0 bit-exact 171740bd, but AK FFE snag-model
bit-exactness is UNVERIFIED — verify jl-vs-oracle FVS_SnagSum first, then SnagDet). With those three resolved the
serializer is a direct transcription of the formula above.

## SnagSum/SnagDet REAL blocker isolated 2026-08-21 — it's the FFE snag-REPORT bit-exactness, not the wiring
Wired SNAGSUDB → new s.control.dbs_snagsum → carb_rows collection + write_dbs_snagsum! (state.jl field +
keyword_dispatch DATABASE handler + simulate.jl gate). VERIFIED the wiring emits FVS_SnagSum, gate 339/11 held.
BUT the CONTENTS diverge from the oracle, so per doctrine (commit only validated OUTPUT) the wiring was REVERTED:
  - **AK FFE is UNPORTED in jl** (fmcba.jl:33 errors: ORGANON species carry no FFE crown-biomass coeffs) — AK cannot
    be the SnagDet/SnagSum fixture at all.
  - **EM FFE emits FVS_SnagSum but the snag densities DON'T match** (emt01, 3 cycles):
      oracle 1990:H14.76  2000:H17.38  2010:H15.97
      jl     1990:H0.0    2000:H10.80  2010:H21.46  2020:H26.46/S5.51
    Two divergences: (a) jl reports 0 hard snags at 1990 while the oracle has 14.76 — an INITIAL-SNAG seeding gap (dead
    trees in emt01.tre / first-cycle mortality timing); (b) jl emits an extra 2020 row = a cycle-LABELING offset (jl's
    per-cycle carb_rows collection point vs FMSOUT's report timing differ by ~one cycle).
⇒ The DBS snag tables (FVS_SnagSum AND FVS_SnagDet) are blocked on the per-variant FFE **snag-report** bit-exactness
(initial-snag seeding + report-timing), NOT on the trivial DBS wiring. FFE FIRE is validated cluster-wide, but the FFE
snag-DENSITY report is a separate, un-validated sub-output. To land the DBS snag tables: first make the EM (or CA/CR/
WS/SO/WC) FFE snag report bit-exact vs FVS<v>_clean (align initial-snag seeding + the per-cycle report point), THEN the
SNAGSUDB/SNAGOUDB wiring (verified trivial here) + the SnagDet serializer (FMSOUT formula recorded above) drop in.

## Refinement 2026-08-21 — the snag-report divergence is a REPORT-TIMING offset, not a fundamental blocker
emt01.tre has NO dead-tree records (col-48 I1 = crown-class codes 3-8, not mortality flags). So the oracle's 14.76
hard snags at 1990 = the FIRST cycle's PROJECTED mortality, reported by FMSOUT at the cycle-START label; jl reports
its snags at the cycle-END label (2000) → jl's rows are shifted ~one cycle (hence jl's extra 2020 row and its 1990=0).
The snag DENSITIES themselves are downstream of each side's mortality, which is the accepted cornered DGSCOR/OLDRN
straddle — so FVS_SnagSum/SnagDet contents will be CORNERED-matching (like every .sum multi-cycle column), NOT bit-
exact, which MEETS the doctrine bar. ⇒ the DBS snag tables are NOT fundamentally blocked; the one concrete fix is to
ALIGN jl's snag-report cycle labeling with FMSOUT's (report the cycle's projected snags at the cycle-start year), after
which SnagSum/SnagDet track the oracle cornered. That timing-alignment (in the carb_rows collection point / the
write_dbs_snagsum! year field) is the executable next step — a focused fix, then the reverted SNAGSUDB wiring + the
SnagDet serializer (FMSOUT formula above) drop in and validate cornered on EM/CA/CR/WS/SO/WC (a validated-FFE variant).

## CORRECTION 2026-08-21 — the snag divergence is a real FFE fall-down gap, NOT just timing (my prior refine over-claimed)
Re-reading the actual EM FFE FVS_SnagSum numbers: jl's hard snags ACCUMULATE monotonically (2000 H10.8 → 2010 H21.46
→ 2020 H26.46) while the oracle's stay ~FLAT (1990 H14.76 → 2000 H17.38 → 2010 H15.97). A one-cycle timing shift does
NOT reconcile them (jl-shifted 10.8/21.46/26.46 vs oracle 14.76/17.38/15.97 — jl monotone-rising vs oracle flat). So
jl's EM FFE STANDING-snag report retains snags that the oracle SHEDS (fall-down) — a genuine FFE snag fall-rate/decay
divergence in the snag_summary (current-standing) computation, deeper than the wiring or a labeling offset. (The
down-wood StandDead falldown is separately validated, but the STANDING-snag_summary density evidently isn't.)
⇒ HONEST STATE: the DBS snag tables (FVS_SnagSum/SnagDet) are blocked on validating the per-variant FFE STANDING-snag
report itself (fall-down/decay of the snag_summary), a focused FFE-snag-model dive — NOT the trivial DBS wiring
(verified) nor a mere timing alignment. The emission recipe + FMSOUT aggregation formula + SNAGSUDB wiring are all in
hand and drop in once the snag_summary fall-down is bit-exact-or-cornered on a validated-FFE variant.

## ★ RESUME HANDOFF 2026-08-21 (container restart) — OC oracle A/B is SET UP; run the jl side next
PPE landscape harness DONE + committed (cccf01b4) — unrelated, that item is closed.
DBS FVS_SnagDet: switched the A/B fixture from IE (FFE unvalidated) to **OC (FFE bit-exact-or-cornered)** per the
recipe's own caveat. Oracle emission WORKS on OC:
  - Keyfile: `/workspace/.ocwork/ocsnag.key` (+ `ocsnag.tre`), also in `scratchpad/ocffe/`. Derived from
    ocffe_full.key: added FFE `SNAGSUM`+`SNAGOUT` and a DBS output block (DSNOUT ocsnag_oracle.db / SUMMARY /
    SNAGSUDB / SNAGOUDB / END). Bare TREEDATA (auto-reads <base>.tre).
  - RUN VIA: `./FVSoc_clean --keywordfile=ocsnag.key` (NOT stdin; stdin gave STOP 10 +
    all-zero .sum = trees didn't load).
  - Oracle DB `ocsnag_oracle.db` (also copied to scratchpad/ocffe/): **FVS_SnagSum 12 rows, FVS_SnagDet 379 rows.**
  - sqlite3 CLI is ABSENT — read via Julia `using SQLite, DBInterface` (DataFrames NOT in the project; iterate rows).
  - Oracle FVS_SnagSum hard densities (soft_total=0 all cycles; classes are CUMULATIVE ≥SNPRCL(c), so hTot=h1):
    1993 h1=55.26 · 1998 92.99 · 2003 264.90 · 2008 51.01 · 2013 8.67 · 2018 9.66 · 2023 9.25 · 2028 8.07 ·
    2033 6.87 · 2038 6.01 · 2043 5.30 · 2048 4.71  (REALISTIC fall-down — NOT monotone; a clean bit-exact target.)

jl side matches this structure EXACTLY: `snag_summary(s)` (src/engine/fire/snag.jl:424) returns
`(hard=7-tuple, soft=7-tuple)` where slots 1-6 = CUMULATIVE density `d≥_FM_SNPRCL[c]` and slot 7 = total — the same
cumulative semantics as FVS_SnagSum. `_FM_SNPRCL` (recipe open-item (a)) is ALREADY resolved in jl. So:

NEXT STEP (immediate, ~2 min): run jl OC on `scratchpad/ocffe/ocffe_full.key` capturing `snag_summary(s)` each cycle
(instrument simulate.jl's FFE report point at simulate.jl:916 / summary.jl:268, or a per-cycle callback), and A/B the
hard/soft cumulative densities vs the oracle table above.
  - IF bit-exact-or-cornered on OC ⇒ the snag MODEL is faithful ⇒ port `write_dbs_snagdet!` (pure serialization; full
    aggregation spec already extracted from fmsout.f in this doc: (species IDC, JYR=IYR−yrdead+1 clamp 1..100, JCL by
    SNPRCL), density-weighted Ht/DBH, summed Vol, normalize by TOTN). Resolve the last open lookup: SNVOLH/SNVOLS
    volume definition vs jl SnagList fields (bolevol vs a computed cubic — check against oracle FVS_SnagDet
    Current_Vol_Hard/Soft, which are INT()-truncated). Then wire SNAGOUDB toggle + A/B FVS_SnagDet column-by-column.
  - IF jl snag densities diverge from the oracle on OC ⇒ it's a snag-MODEL divergence (like the prior EM monotone
    finding); characterize + corner, do NOT force the serialization to a divergent target.
Doctrine: commit write_dbs_snagdet! ONLY if the OC A/B is bit-exact-or-cornered. Same two-block recipe reusable for
FVS_StrClass / FVS_CanProfile / FVS_Climate afterward.

## ★★ FVS_SnagDet OC A/B RESULT 2026-08-21 — CORNERED (blocked on OC FFE snag-report bit-exactness)
Ran the A/B: jl OC `snag_summary(s)` per cycle (via ocsnag_jl.key = ocffe_full.key + CARBREPT, captured through
write_sum_file carbon_collect[i][4]) vs oracle FVS_SnagSum/FVS_SnagDet. **The snag DENSITIES DIVERGE** — the DBS
serializer is de-risked (schema + aggregation spec in hand, jl SnagList carries every field) but CANNOT be validated
bit-exact because the underlying per-cycle snag pool differs. Do NOT commit write_dbs_snagdet! (doctrine: commit only
validated output). Precise characterization (measure-don't-infer, the yrDied decomposition is the key tool):

  hTot (total snag density):   Year   Oracle    jl(cyc-start)   jl(cyc-end)
                               1993    55.26     14.76           75.43
                               2003   264.90    194.12           38.28
                               2013     8.67      4.67            3.09
                               2033     6.87      0.64            0.58
  Large-snag classes (h2/h3) match closely; the divergence is in the SMALLEST class + its late-cycle fall-down.

ROOT (oracle FVS_SnagDet yrDied decomposition — the oracle's snag-DATING model jl books differently):
  - Year=1993: yd1988=14.8 (input dead) + **yd1993=40.5** (INVENTORY-YEAR mortality, tiny GF/DF/BR DBH 0.1-0.68).
    jl's 1993 = 14.76 = EXACTLY the yd1988 subtotal ⇒ jl omits the yd1993 current-cycle mortality cohort at the
    1993 report point (jl dates ordinary mortality yrdead=cycle-END−1=1997, and its cycle-start report predates the
    grow_cycle! mortality entirely).
  - Year=1998: yd1988=10.2, **yd1991=35.8 (= SNAGINIT age-2 snags, 1993−2=1991)**, yd1993=29.4, yd1997=17.6.
    SNAGINIT snags first surface at 1998 in the oracle, dated 1991; jl adds/dates them differently.
  - Year=2003: … yd2003=207.5 (SIMFIRE-killed snags).
  Three distinct snag-input timings (inventory mortality / SNAGINIT age-dating / fire) that jl books on a different
  schedule ⇒ a deep FFE snag-report + snag-dating model divergence, downstream of & entangled with the KNOWN cornered
  OC growth/mortality drift (goal-doc: OC .sum mort 213 vs 266 = downstream of cornered growth). NOT a serialization
  bug and NOT the DBS-writer's fault.

VERDICT: DBS FVS_SnagDet port = CORNERED — serializer ready, blocked on the OC FFE snag-report/dating bit-exactness
(an FFE-model investigation, not DBS plumbing). Same conclusion the prior session reached on EM, now confirmed
precisely on the validated-FFE OC variant with the yrDied decomposition. The oracle-EMISSION recipe + serialization
spec remain valid and reusable IF/when the FFE snag-report timing is made bit-exact. Other DBS tables whose source
model IS bit-exact (e.g. FVS_StrClass from structure_stage, FVS_Compute already done) are better serialization
candidates than the snag/fire-derived tables. Artifacts: /workspace/.ocwork/ocsnag*.{key,tre,db},
scratchpad/ocffe/{ocsnag_jl.key, run_snag_ab.jl}.

## ★★ DBS extension-table TRIAGE 2026-08-21 (measured — which of the 7 missing tables are validatable)
After the SnagDet A/B (cornered), triaged the remaining missing DBS tables by whether their SOURCE MODEL can be
validated bit-exact-or-cornered (the serializer is the easy part; the target data is the gate):
  1. **FVS_SnagDet** — CORNERED (this doc above): FFE snag-report/dating diverges even at cyc0 (jl 14.76 vs 55.26).
  2. **FVS_BM_*** (WWPB, 4 tables) — NO ORACLE: WWPB is the source-absent PPE reconstruction (wwpb_outbreak_cycle!),
     so there is nothing to A/B against bit-exact. Un-validatable by construction (same class as PPE item 5).
  3. **FVS_DM_*** (mistletoe/NEWSPRED, 3) — ORACLE DB-CRASH-BLOCKED: FVSbc_clean SIGSEGVs on all DATABASE reads
     (isoc23-shim × SQLite interop), and NEWSPRED's only corpus case is the YSM DB stand ⇒ no A/B vehicle.
  4. **FVS_RD_*** (WRD root disease, 3) — jl WRD is linked-but-DORMANT/partial; source model not in the validated set.
  5. **FVS_StrClass** — 44-col × 3-strata, HEAVILY tree-list-derived (strata by ht/DBH/crown/species-dominance).
     Would corner on growth drift + stratum-boundary straddle, AND needs jl structure_stage to expose per-stratum
     detail (DBH, nom/lg/sm ht, crown base/cover, top-2 species, status per stratum) = a big NEW plumbing lift.
  6. **FVS_CanProfile** (fmpocr DBSFMCANPR) — canopy cover by height; tree-list/crown-derived ⇒ cyc0-exact, multi-
     cycle cornered; reuses jl's existing canopy-profile machinery but still needs a collection hook + oracle A/B.
  7. **FVS_Climate** (dbsclsum ICLIM) — the CLEANEST: per-species Viability/GrowthMult/SiteMult/MxDenMult are
     DETERMINISTIC climate-reader outputs (bit-exact-able); BA/TPA/ViabMort/dClimMort/AutoEstbTPA are tree-list/
     model-derived (cyc0-exact per the validated Climate-FVS, multi-cycle cornered). jl's climate.jl computes all of
     these INLINE but does NOT collect them per-species-per-cycle ⇒ needs a collection hook (like carbon_collect) +
     a CLIMATE-data oracle fixture (CLIMATE kw + climate projection file + DBS CLIMATEDB toggle) + serializer.
VERDICT: model-derived tables (1-4,6) corner or lack an oracle; the classification/reader tables (5,7) are cyc0-
exact-at-best but need substantial NEW jl per-cycle collection plumbing for a multi-cycle-cornered result. The DBS
extension-table item is thus "cornered-or-heavy-plumbing-for-cornered" — NOT a stream of clean bit-exact wins. The
oracle-EMISSION recipe (two-block DATABASE + report keyword + blank toggle) is proven + reusable for any of them.
Best next candidate if pursued: FVS_Climate (add a per-species climate-report collect vector in the climate apply
path, mirror write_dbs_carbon!, A/B on IE with a climate fixture) — cyc0-bit-exact deliverable, multi-cycle cornered.

## ★ FVS_Climate PROBE 2026-08-21 — also needs report-timing/sampling work (not a clean serialization either)
Set up the oracle (fixture EXISTS: /workspace/.iework/climate/clim_iet.key + FVSie_clean): added CLIMREDB (DBS opt 8,
sets ICLIM) to the DATABASE block → climtest.key emits FVS_Climate (429 rows). Emission spec fully extracted
(clauestb.f:178-207, called per cycle during climate auto-establishment): one row per species where SPIMP>0.05 OR
SPVIAB>0.4 (and INDXSPECIES>0), columns SPVIAB(=VSCORE per clgmult.f:117)/SPBA(ΣDBH²·PROB·0.005454154)/SPTPA(ΣPROB)/
SPMORT1/SPMORT2/SPGMULT(=Σtreemult·PROB/ΣPROB, clgmult.f:243-249)/SPSITGM(=xgsite^clgrowmult)/MXDENMLT/POTESTAB.
Feasibility gate PASSED — jl computes every piece per-tree (apply_climate_dds! treemult, species_vscore, clim_mort_rates,
clmaxden, clauestb POTESTAB); the report is a per-species TPA-weighted aggregation of them.
BUT the cyc0 A/B (probe_climrep.jl) does NOT match the oracle:
  - Viability: jl DF/PSME raw spviab 0.95 vs oracle 0.9309; WL 0.879 vs 0.8584; every species off ⇒ a viability
    SAMPLING/transform difference (oracle Viability=VSCORE not raw spviab; jl 0.95>0.5→VSCORE 1.0 ≠ oracle 0.9309 ⇒
    the two are sampling the CSV viability column at a DIFFERENT year, or a species→PLANTS-column mismatch — needs
    a clgmult.f VSCORE-year/column trace).
  - BA/TPA: oracle DF BA=28.82/TPA=157.63 vs jl (pre-growth inventory) 19.45/206.96 ⇒ the oracle's clauestb SPBA/SPTPA
    are sampled POST-growth (report-timing), like the SnagDet FMSOUT-vs-FMCRBOUT timing.
  - Filter: WH (TSHE, viab 0.064) appears in jl but not the oracle ⇒ the SPIMP>0.05 importance filter must be ported.
VERDICT: FVS_Climate = also needs a climate-report-timing + VSCORE-sampling investigation before it can validate — NOT
a clean cyc0 serialization win. Consistent with SnagDet. Oracle + probe saved (/workspace/.iework/climate/{climtest.key,
climtest_oracle.db,probe_climrep.jl}) for a future session that takes on the climate-report-timing chunk. ⇒ CONFIRMED:
the whole DBS extension-table long-tail is "needs per-table model-report-timing work + multi-cycle-cornered", not a
stream of clean bit-exact serializations. The oracle-EMISSION recipe (proven for SnagSum/SnagDet/Climate) is the durable
reusable asset.

## ★★★ BREAKTHROUGH 2026-08-21 — the DBS-table "divergence" is REPORT-COLLECTION TIMING, jl values are CORRECT
Settled the FVS_Climate viability discrepancy DEFINITIVELY. Inline PSME(DF) viability knots for S248112: 1990=0.95,
2030=0.797. The oracle's FVS_Climate row LABELED Year=1990 reports Viability=0.9309 = PSME interpolated at **1995**
(= 1990 + fint). Sampling jl at 1995 reproduces the oracle's "1990" row BIT-EXACT for EVERY species:
  DF 0.9309 · WL 0.8584 · PP 0.4335 · GF 0.8869 · LP 0.8412 · ES 0.8818  (all == oracle).
⇒ **FVS reports the DBS extension tables POST-GROWTH, sampling at the cycle-END year (report_label + fint).** jl's
per-cycle carbon-collect samples at cycle-START (pre-growth), which is why the naive A/B diverged. jl's underlying
climate reads are CORRECT — this is a COLLECTION-POINT bug, not a model divergence. The deterministic climate columns
(Viability/GrowthMult/SiteMult/MxDenMult) are therefore BIT-EXACT-PORTABLE once collected at the post-growth point.
This also explains the SnagDet timing (FMSOUT post-growth) — though SnagDet ALSO has a snag-dating divergence on top,
so it still corners; Climate is clean for the deterministic columns (BA/TPA/mort corner multi-cycle on the tree list).
ACTION: port FVS_Climate collecting per-species at the post-growth cycle-end sampling year (report_year + fint):
climate_report(s, sample_year) → (SPVIAB,SPBA,SPTPA,SPMORT1,SPMORT2,SPGMULT[=Σtreemult·prob/Σprob],SPSITGM[=xgsite^clgrowmult],
MXDENMLT,POTESTAB) with the SPIMP>0.05||SPVIAB>0.4 filter (clauestb.f:196); write_dbs_climate! (dbsclsum.f schema);
wire CLIMREDB toggle. Expect Viability/GM bit-exact, BA/TPA cornered.

## ★★ FVS_Climate climate_report VALIDATED (deterministic cols bit-exact) + 2 gaps found 2026-08-21
Built climate_report(s; report_year, fint) (WIP saved scratchpad/climate/climate_report_wip.jl; validation driver
scratchpad/climate/validate_climrep.jl) mirroring clauestb.f:178-207 + clgmult.f, collecting POST-grow_cycle! at
sample year = report_year + fint/2. A/B vs the oracle "1990" row (post-cyc1 state), all 8 species:
  BIT-EXACT ✓  Viability (WL 0.8584/DF 0.9309/GF 0.8869/LP 0.8412/ES 0.8818/AF 0.7528/PP 0.4335/MM 0.9495 — ALL match),
               SiteMult (1.0), MxDenMult (1.0).
  CORNERED ✓   BA/TPA (WL 13.47 vs 13.49, ES 12.08 vs 12.14 — cyc1 growth straddle, expected).
               GrowthMult: PP 0.7783 EXACT (vscore<0.99 ⇒ deterministic); DF/GF/LP >1.0 (1.0303 vs oracle 1.0) —
               vscore>0.99 ⇒ ps=MAX(xgsite,xrelgr,vscore) (clgmult.f:222-226 clim_treemult) so gm depends on the
               per-tree Leites XDF (birthyr) ⇒ tree-list-dependent ⇒ CORNERED like BA/TPA.
  ⇒ CONFIRMS the timing breakthrough end-to-end: sampling at report_year+fint/2 post-growth reproduces the oracle's
  deterministic climate columns BIT-EXACT.
TWO REMAINING GAPS (genuine unported model pieces, NOT serialization — block a faithful write_dbs_climate!):
  1. **Mort1/Mort2 need the SPCALIB presence-calibration** (clmorts.f:57-75, ICYC=1 branch — SPCALIB(I)=viab·0.9,
     mortality computed RELATIVE to it ⇒ SPMORT1=0 at cycle 1). jl's apply_climate_mort! computes the RAW FYRMORT
     (PP M1=0.2217) with NO SPCALIB ⇒ over-reports at cycle 1 (oracle M1=0). This is the "clmorts.f:92-98 first-cycle
     presence-calibration = chunk C, unported" the climate memory already flagged. Porting SPCALIB fixes the report
     Mort AND the actual multi-cycle climate mortality (a real, .sum-affecting improvement, not just report-cosmetic).
  2. **POTESTAB needs an AutoEstb activity** (clauestb.f:76 `IF(PTREES*AESNTREES>0)`; LAESTB set only by an AutoEstb
     keyword, MYACT=2802). The clim_iet.key CLIMATE block has NO AutoEstb keyword, yet the oracle reports POTESTAB=
     99.44 for DF/GF/ES/MM ⇒ the oracle has an auto-estab source (default AESNTREES? a keyword form jl isn't parsing?)
     UNRESOLVED — jl reports 0. Resolve by DEBUG-dumping clauestb (AESNTREES/LAESTB) in the oracle to find where 99.44
     comes from with no AutoEstb keyword, then match jl's clim_autoestb parsing.
⇒ FVS_Climate = deterministic core VALIDATED + timing solved; a faithful full-table commit awaits SPCALIB (chunk C,
worth porting for the .sum too) + the POTESTAB/AutoEstb-source resolution. climate_report reverted from src (partially
faithful ⇒ not committed per doctrine); WIP preserved. This is the pattern for ALL the model-report DBS tables:
the deterministic columns port cleanly at the post-growth midpoint; the model-derived columns need their model piece
bit-exact first.

## ★ POTESTAB gap RESOLVED + SPCALIB PORTED 2026-08-21 (33290ff0) — FVS_Climate now fully characterized
POTESTAB=99.44 with NO AutoEstb keyword = the clinit.f:40-42 DEFAULTS: AESNTREES=500, NESPECIES=4, AESTOCK=40.
clauestb.f runs every cycle when LCLIMATE and computes POTESTAB with those defaults (line 76 `IF PTREES*AESNTREES>0`
— 500>0 always); LAESTB (set only by an AUTOESTB keyword) gates the ACTUAL establishment, NOT the report. The
NESPECIES=4 default ⇒ the oracle's 4 rows (DF/GF/ES/MM, the top-4 by viability). FIX for the (reverted) WIP
_climate_potestab: use AESNTREES=500/NESPECIES=4/AESTOCK=40 when `isempty(c.autoestb)` instead of returning 0.

**SPCALIB (clmorts chunk C) is now PORTED + committed (33290ff0)** — the hard/valuable part. ⇒ FVS_Climate is
FULLY characterized, every column portable:
  - Viability / SiteMult / MxDenMult — BIT-EXACT (deterministic, sampled at report_year+fint/2 post-growth).
  - ViabMort (SPMORT1) — BIT-EXACT via the now-ported SPCALIB (PP series matches oracle all 11 cycles).
  - dClimMort (SPMORT2) — the transfer-distance DMORT, already ported (apply_climate_mort! ldmort branch).
  - POTESTAB — via clinit defaults (500/4/40), resolved above.
  - GrowthMult (vscore>0.99 species) / BA / TPA — CORNERED on the cyc1 growth straddle (accepted).
REMAINING to ship write_dbs_climate!: re-add climate_report (WIP scratchpad/climate/, + the POTESTAB-default fix),
add the write_dbs_climate! schema/serializer (dbsclsum.f 15-col), a post-growth climate_collect hook in summary.jl
(fires right after grow_cycle!, labels with the pre-advance cycle year), and the CLIMREDB DBS-toggle parse. All
mechanical now that every column is understood + SPCALIB is in src. A future session's clean serialization chunk.

## ★ write_dbs_climate! — remaining pieces precise (2026-08-21)
Correction to "all mechanical": the ONE remaining non-trivial piece is the dClimMort (SPMORT2) column =
Σ(DMORT·CLMRTMLT2·X)/SPWTS per species (clmorts.f:165-268) — a per-tree-weighted aggregation of the transfer-
distance DMORT (tree-list-dependent ⇒ cornered). Cleanest impl: have apply_climate_mort! STORE per-species
c.spmort1 (clim_mort_rates[1], now SPCALIB-calibrated) + c.spmort2 (the ldmort-branch DMORT aggregated) into new
ClimateState fields as it runs, then climate_report reads them (avoids duplicating the transfer-distance loop +
guarantees the report matches the APPLIED mortality). The other columns are done: Viability/SiteMult/MxDenMult
bit-exact, ViabMort bit-exact (SPCALIB in src 33290ff0), POTESTAB via clinit defaults (500/4/40), BA/TPA/GrowthMult
cornered. Then: write_dbs_climate! (dbsclsum.f 15-col INSERT), a climate_collect hook in summary.jl firing right
after grow_cycle! (label=pre-advance cycle year, viab sampled report_year+fint/2), CLIMREDB DBS-toggle parse, and
the run_keyfile write. A/B via /workspace/.iework/climate/climtest.key (CLIMREDB→429 rows) — deterministic cols
bit-exact, tree-list cols within the cyc1 straddle. WIP climate_report: scratchpad/climate/climate_report_wip.jl
(update: mort1→clim_survival_cal(spviab,c.spcalib[sp]); mort2→c.spmort2; POTESTAB defaults when isempty(autoestb)).

## ★★★ FVS_Climate DONE 2026-08-21 (c7e9aeeb) — first missing DBS table ported to a full oracle A/B
Complete end-to-end: climate_report (POST-grow_cycle! collection at report_year+fint/2 midpoint) + write_dbs_climate!
+ ClimateState.spmort1/spmort2 (from apply_climate_mort!) + Control.dbs_climate/CLIMREDB toggle + post-growth
climate_collect hook in write_sum_file + the run_keyfile write. Vs FVSie_clean: Viability+ViabMort BIT-EXACT 147/147,
rest cornered on the DGSD=2.0 OLDRN straddle. test_climate_dbs 13/13, gate 339/11. +2 real bug fixes (SPCALIB 33290ff0,
POTESTAB no-column ranking). ⚠ FALSE-POSITIVE LESSON: an accidental FVSie_clean run into the jl output db gave a bogus
"all-columns bit-exact" (oracle-vs-oracle); ALWAYS rm the jl db + verify jl-only rows before an A/B.

## REUSABLE PATTERN for the remaining deterministic-reader DBS tables (proven by FVS_Climate)
1. Control.dbs_<x> flag + the DBS toggle keyword in kw_database! (the toggle name from dbsin.f's KEYWRD table).
2. A <x>_report(s; report_year, fint) computing the per-<row> report at the POST-grow_cycle! state, sampling any
   time-series at report_year+fint/2 (the FVS midpoint); read model-derived per-species/row values FROM state
   (stored by the model routine that already ran in grow_cycle!) so the report == the applied values.
3. write_dbs_<x>! (schema from dbs<x>.f) + a <x>_collect hook in write_sum_file right after grow_cycle! (label =
   pre-advance cycle year) + the run_keyfile write, gated on the Control flag.
4. A/B on a variant whose underlying MODEL is bit-exact for the table's data; deterministic columns bit-exact,
   tree-list columns cornered on the OLDRN straddle. Test with a fixture keyfile.
REMAINING (status): FVS_SnagDet CORNERED (FFE snag-report/dating diverges even at cyc0). FVS_BM_* NO ORACLE (WWPB
reconstruction). FVS_DM_* DB-CRASH-BLOCKED (NEWSPRED). FVS_RD_* jl-WRD partial. FVS_StrClass 44-col/3-strata heavy +
needs structure_stage to expose per-stratum detail (all tree-list ⇒ fully cornered). FVS_CanProfile canopy-cover-by-
height (fmpocr, tree/crown-derived ⇒ fully cornered, needs a canopy-profile collect). ⇒ FVS_Climate was the cleanest
(2 deterministic columns); the rest are blocked or fully-cornered-tree-list.

## ★ FVS_CanProfile 2026-08-21 — serializer TRIVIAL + ready, but blocked on jl's OC canopy-profile MODEL divergence
Built the full FVS_CanProfile port (write_dbs_canprofile! + canopy_crfill extraction + Control.dbs_canprofile +
CANFPROF FFE keyword parse + post-growth canprof_collect hook + run_keyfile write) — the serializer is trivial: jl's
canopy_bulk_density ALREADY builds the exact CRFILL array (crown fuel by 1-ft layer, lbs/ac-ft) that DBSFMCANPR
serializes; I extracted it as `canopy_crfill(s)`. Emission spec (dbsfmcanpr.f): 1 row per layer I with CRFILL(I)>0 —
Height_ft=I, Height_m=I·0.3048, fuel lbs/ac-ft, kg/m³=·0.45359237/(4046.856422·0.3048). Oracle enable = FFE keyword
`CANFPROF` (fmin.f opt 47, sets ICANPR + year window; NO separate DBS toggle). Oracle vehicle /workspace/.ocwork/
occanpr.key (+.tre) → occanpr_oracle.db (1025 rows). Gate 339/11 held (refactor behavior-preserving).
BLOCKER — the A/B on OC does NOT validate: jl's OC canopy_crfill diverges STRUCTURALLY from the oracle even at cyc0
(NOT cornering). 1993: oracle ht8→11.4, ht12→238.9 (plateau); jl ht5-7 present, ht12→57.7 (~4× low, smooth ramp).
Pre-growth jl profile also diverges ⇒ NOT a timing issue. ROOT (pre-existing OC FFE canopy gap, exposed by this port,
NOT the serializer): (1) OC has NO `fm_canopy_lsw(::OregonCoast)` method ⇒ falls back to `sp<=25` — wrong LSW
softwood set for the 50-species ORGANON OC (includes OC hardwoods / excludes conifers 26-50); (2) the ~4× magnitude
gap also implicates the OC/ORGANON crown_biomass (foliage+½finewoody) feeding adcrwn. ⇒ REVERTED per doctrine
(commit only validated output). To finish: either (a) validate on a variant WITH a proper fm_canopy_lsw + bit-exact
FFE canopy (BM `1:14+17` / CR `1:19+29:37` / SN `≤17+88` / NE `≤25`) — its canopy_crfill is bit-exact so the profile
A/B validates cyc0-exact + OLDRN-cornered; or (b) first fix the OC canopy model (add OC LSW + verify ORGANON crown
biomass). The serializer + wiring are correct + gate-safe; only the OC canopy MODEL blocks the A/B. WIP not kept in
src (reverted); the whole port is ~40 lines re-appliable from this note. NB canopy_crfill(s) extraction is independently
useful (canopy_bulk_density + FVS_CanProfile share it).

## ★★★ FVS_CanProfile DONE 2026-08-21 (4e8ff2c7) — 2nd DBS table validated (on CR, not OC)
The FVS_CanProfile serializer (from the prior note) is now VALIDATED + committed. Two corrections to that note:
  1. TIMING: FVS_CanProfile reports the PRE-growth cycle-START inventory (fmpocr mode 2 alongside FMPOFL), NOT
     post-growth. jl collects canopy_crfill BEFORE grow_cycle! (at the potfire_collect point). The oracle's 229.18
     plateau matches ONLY at pre-growth. ⇒ DBS reports have DIFFERENT timings: Climate=post-growth-midpoint,
     CanProfile=pre-growth-cycle-start, SnagDet/carbon=post-growth. Determine per-table which.
  2. VARIANT: validated on CR (crt01 FFE stand + CANFPROF vs FVScr_clean) — cyc0 68/68 BIT-EXACT, multi-cycle
     cornered on the CR OLDRN straddle. CR has a proper fm_canopy_lsw + validated cr_crown_biomass ⇒ bit-exact
     canopy_crfill. The OC failure was BOTH the wrong timing AND the OC canopy gap below.
★ OC-FFE BUG LEAD surfaced: jl's OC canopy_crfill diverges structurally from the oracle (no fm_canopy_lsw(::Oregon
Coast) ⇒ wrong sp≤25 LSW fallback for the 50-species ORGANON + likely ORGANON crown_biomass ~4× low). Add an OC LSW
method + verify OC crown_biomass to close it (a separate OC-FFE chunk; the CanProfile port doesn't need it).
⇒ TWO DBS tables now DONE (FVS_Climate + FVS_CanProfile). Remaining: SnagDet cornered; BM/DM/RD blocked; StrClass
44-col heavy needs structure_stage per-stratum exposure. The <x>_report + write_dbs_<x>! + collect-hook pattern is
proven twice; the only per-table decision is the collection TIMING (pre vs post growth) + the enabling keyword.

## ★ OC canopy_crfill divergence — REFINED 2026-08-21 (it's crown_biomass, NOT the LSW)
Investigated the OC-FFE canopy bug lead. Read OC fmvinit.f LSW (fire/oc/fmvinit.f:138-320): LSW=TRUE for species
**1-25 + 50** (26-49 FALSE). So the jl AbstractVariant fallback `sp≤25` is NEARLY RIGHT for OC (misses only sp==50,
absent from the test stands) — the LSW is NOT the divergence cause. The real gap is OC/ORGANON **crown_biomass ~4×
LOW**: at the SAME pre-growth timing, oracle 1993 ht12=238.89 vs jl 57.68 (ratio 4.14), and jl spreads crown to
lower heights (ht 5-7) the oracle lacks (crown-base/crown_pct difference). ⇒ closing OC FVS_CanProfile needs an
ORGANON crown-biomass fix (crown_biomass(s,sp,...) for the 18 ORGANON species returns ~1/4 the oracle's foliage+
finewoody) + a crown-base check — a deeper OC-FFE chunk, LATENT (.sum-inert on OC test stands: OC crown fire doesn't
fire, so canopy_bulk_density only affects the report + the un-triggered crown-fire index). A future OC-FFE session:
add fm_canopy_lsw(::OregonCoast)=(1≤sp≤25)||sp==50 (source-faithful, tiny) AND dump FVSoc_g16 fmpocr per-tree
crown-biomass to find the ~4× ORGANON factor. FVS_CanProfile itself is DONE (validated on CR); this OC gap is a
separate OC-FFE model item, not a CanProfile-port gap.

## ★★ OC crown-biomass bug FIXED 2026-08-21 (40298123) — the ~4× gap was eastern-vs-western FMCROWE/FMCROWW
The refined OC-FFE lead is now a validated fix. OC's fmcrow.f (== WC's) routes hardwoods {35,39,40,41,43,44,45,46}
→ FMCROWE (eastern) and all other species → FMCROWW (western), but jl omitted OregonCoast from the cr_crownw
dispatch ⇒ ALL OC species used the eastern FMCROWE ⇒ ~4× under-counted conifer crown fuel. Added OregonCoast to
the cr_crownw path (OC_FFE_ISPMAP already existed; new oc_uses_fmcrowe = the 8 hardwoods). OC cyc0 canopy profile
0/71 → 62/71 BIT-EXACT (ht12 57.7→253.9 vs 238.9; ht20 exact). OC FFE .sum UNCHANGED (surface fire; crown biomass
inert on byram/mort). Gate 339/11. Residual 9 crown-base layers = a smaller secondary OC crown-base/crown_pct detail.
NOTE: **OP (Olympic) is ALSO ORGANON — likely the same bug** (check OP_FFE_ISPMAP + an op_uses_fmcrowe if OP FFE is
exercised). FVS_CanProfile stays validated on CR (68/68, the clean bit-exact case); OC is now 62/71 (was 0).

## OC FVS_CanProfile cyc0 residual (62/68) ROOT-CAUSED 2026-08-21 — small no-height-input DF trees (CORNERED)
Measured the OC canprofile cyc0 (1993) per-layer A/B vs FVSoc_clean (instrumented canopy_crfill per-tree, then
rebuilt the crfill in python from the 23 cyc0 trees). The entire 6-layer residual is EXACTLY two trees: input
records 12 & 13 = small DF (dbh 1.2" and 1.9", crown-code present but **HEIGHT FIELD = 0** in occanpr.tre). FVS
imputes their missing height; jl imputes h=11 / h=13 (both > CANMHT=6 ⇒ included in the canopy profile), but the
oracle EXCLUDES them entirely (oracle crfill is 0.00 at ht6-7 where only these trees reach; ht8-9 = exactly tree-2's
partial/full contribution, no DF-small addition). Confirmed NOT a crown-biomass-value diff: jl's cr_crownw already
does the fmcroww SMWGT small/large blend correctly, and both the suppressed & dominant DF large-tree formulas give
nonzero weight at D=1.2 (hand-computed ~3.3 lb/tree) — so the oracle's exclusion is a TREE-INCLUSION difference, i.e.
the ORGANON missing-height imputation for sub-2" DF gives a height ≤ CANMHT (6 ft) that jl's imputation does not.
CORNERED: (1) `.sum`-INVISIBLE — OC cyc0 .sum is bit-exact (these tiny trees don't move TopHt/QMD/BA/TPA); (2)
report-only (FVS_CanProfile), and the table is already validated BIT-EXACT on CR (68/68); (3) OC fire is surface-
dominated (crown fuel inert on byram/scorch/mort). Fixing needs matching ORGANON's exact org_intree missing-height
model for sub-2" trees = a deep OC-input detail with layered low value. LEFT CORNERED.

## OC redwood (sp50) crown-width GAP surfaced 2026-08-21 (separate, pre-existing)
While validating the fm_canopy_lsw(::OregonCoast) redwood-inclusion fix (8ba9d59c, source-verified), a redwood test
stand hit `oc_cwcalc: crown-width equation 21104 (species 50) not yet ported` (oc/cwcalc.f CASE for sp50 absent).
So OC redwood is under-ported in the crown-width path — a distinct ORGANON-species gap (redwood = a real OC species).
The LSW fix is source-faithful + gate-safe + additive/inert on non-redwood stands; full redwood end-to-end A/B is
blocked until oc_cwcalc gets the sp50 CASE. Noted as the next OC-species lead if redwood coverage is prioritized.

## FVS_SnagDet — FULLY SCOPED 2026-08-21 (aggregation + volume routing pinned; ready to implement)
Read the complete emit path. dbsfmdsnag.f just serializes precomputed arrays SDBH/SHTH/SHTS/SVLH/SVLS/SDH/SDS
(species IDC × years-ago JYR × DBH-class JCL 1:6). The REAL aggregation is **fmsout.f:99-208**:
  Per snag record II (jl SnagList index i), skip if (den_hard+den_soft)≤0 OR dbh<SNPRCL(1)=0:
  - HARD flag = the SAME hard→soft dktime flip as snag_summary (snag.jl:458-462): age=iyr−1−yrdead ≥ dktime ⇒ SOFT.
    dktime = (1.24·dcx·d)+(13.82·dcx), dcx = SNAGDCAY override or coef snag_decayx[sp]. (den_soft always soft;
    den_hard → hard iff NOT flipped, else added to soft.)
  - JYR = iyr − yrdead + 1 (cap 100); YRLAST = max JYR.  Year_Died col = iyr − JYR + 1 = yrdead.
  - JCL (NON-cumulative, unlike SnagSum's cumulative ≥): DO JCL=1,5: if dbh < SNPRCL(JCL+1) → that JCL; else 6.
    SNPRCL = _FM_SNPRCL = (0,12,18,24,30,36).
  - Accumulate by (sp, JYR, JCL):  TOTDS += den_soft(+den_hard if flipped);  TOTDH += den_hard if hard;
    TOTHTS += htcur·den_soft(+htcur·den_hard if flipped);  TOTHTH += htcur·den_hard if hard;
    TOTVLS += vol2ht(htcur)·den_soft(+vol·den_hard if flipped);  TOTVLH += vol2ht(htcur)·den_hard if hard;
    TOTDBH += dbh·(den_soft+den_hard)   [weight = TOTAL orig density, both parts].
  - POST (fmsout.f:180-202): TOTDBH /= (TOTDH+TOTDS);  TOTHTH /= TOTDH (0 if TOTDH=0);  TOTHTS /= TOTDS (0 if 0).
  Emit a row per (sp,JYR,JCL) with TOTDH+TOTDS>0. Columns: DBH_Class=JCL, Death_DBH=TOTDBH, Current_Ht_Hard=TOTHTH,
  Current_Ht_Soft=TOTHTS, Current_Vol_Hard=TOTVLH, Current_Vol_Soft=TOTVLS, Total_Volume=TOTVLH+TOTVLS,
  Year_Died=yrdead, Density_Hard=TOTDH, Density_Soft=TOTDS, Density_Total=TOTDH+TOTDS. SpeciesFVS/PLANTS/FIA =
  coef.code_alpha/code_plants/code_fia[sp] (same as write_dbs_climate!).
**VOLUME ROUTING PINNED**: FMSVOL (fmsvol.f:129-152) dispatches on METHC(ISPC): 6/10→NATCRS, 5/8→OCFVOL, else→CFVOL.
initre.f:6365-6367 `DO ISPC=1,MAXSP; METHC=10` ⇒ **OC snags → NATCRS = R9 Clark cubic**, which jl HAS (r9clark_cubic,
used in snag_bole_carbon for LS/NE). OCFVOL returns 0 for OC (fvsvol.f:569 CASE DEFAULT) — a red herring (never hit
since METHC=10). vol2ht(xht) = MAX(0.005454154·HTDEAD, TCF) where TCF = NATCRS total cubic of (dbh, HTDEAD) TRUNCATED
at xht=htcur via CFTOPK when xht<HTDEAD (LTKIL). ocsnag.key has SNAGBRK ⇒ htcur<height ⇒ the CFTOPK truncation + the
hard/soft-differing-height path ARE exercised (all 4 column groups live). NOTE FVS tracks HTIH (hard) vs HTIS (soft)
separately; jl has a single htcur per cohort — verify jl's htcur applies to both (likely equal per cohort at default).
IMPLEMENTATION (next): (1) snagdet_rows(s;iyr) in snag.jl mirroring the above + a _snag_vol2ht(s,i) via NATCRS/
r9clark_cubic + cftopk (reuse snag_bole_carbon's cftopk setup, grab TCF not MCF); (2) _FVS_SNAGDET_CREATE 17-col +
write_dbs_snagdet! in dbs_output.jl; (3) snagdet_collect in write_sum_file at the _carb_push point (same iyr/snag
state as snag_summary); (4) Control.dbs_snagdet + SNAGOUDB in keyword_dispatch (DATABASE block); (5) create+write in
simulate.jl. VALIDATE column-group by group vs /workspace/.ocwork/ocsnag_oracle.db (FVS_SnagDet, 379 rows, yrs 1993-
2048): density first (snag_summary already ~bit-exact vs SnagSum ⇒ high confidence), then Death_DBH, then Ht (tests
FMSNGHT height-loss), then Vol (tests NATCRS+CFTOPK). Corner any column on the OC grown-DBH Float32 straddle per snag_summary.

## FVS_SnagDet VOLUME-column CORRECTION 2026-08-21 (refines the "NATCRS=R9Clark" note above — it was premature)
Deeper read: NATCRS (fvsvol.f:54+) does NOT compute R9 Clark directly — it calls the National Volume Estimator
Library (VOLINITNVB/VOLINIT) with VOLEQ=VEQNNC(ISPC), the per-species equation number. So the OC snag volume =
NVEL(VEQNNC), which may be a Region-6 Behre/DVEE eqn OR Clark depending on the species' VEQNNC — NOT simply jl's
r9clark_cubic, and DISTINCT from OC's live-tree BLMVOL ('B') path. FURTHER, METHC=10 (→NATCRS) is set in initre.f
inside the CRUISE/CVOLUME keyword handler; whether it runs for a bare ocsnag.key (else METHC stays grinit's 999 →
FMSVOL ELSE → CFVOL Behre) is UNCONFIRMED. ⇒ the SnagDet VOLUME columns (Current_Vol_Hard/Soft, Total_Volume) need
an empirical determination of the OC snag volume routine before they can be made bit-exact — a genuine sub-port
(CFVOL Behre regression or NVEL-VEQNNC), NOT a trivial reuse. TRACTABLE-NOW columns (no volume routine needed):
Density_Hard/Soft/Total (snag_summary basis, ~bit-exact vs SnagSum), Death_DBH, Current_Ht_Hard/Soft (jl htcur +
the FMSNGHT height-loss). RECOMMENDED next-session split: implement + validate SnagDet's density/DBH/height columns
bit-exact-or-cornered first (real, deliverable), then resolve the volume routine (empirical: dump one oracle cohort's
Current_Vol_Hard / Density_Hard = per-tree vol, back-solve vs CFVOL-Behre and NVEL to identify the routine) and add
the two volume columns. The full aggregation spec above is exact; only the vol2ht() primitive is open.

## FVS_SnagDet — IMPLEMENTED + A/B'd 2026-08-21 → REVERTED (blocked on jl snag-LIST model divergence, NOT serialization)
Built the full port (snagdet_rows aggregation faithful to fmsout.f + _snag_vol2ht + write_dbs_snagdet! 17-col +
SNAGOUDB toggle + Control.dbs_snagdet + write_sum_file snagdet_collect at the carbon/snag-summary timing +
simulate.jl write) and A/B'd vs /workspace/.ocwork/ocsnag_oracle.db. **The SERIALIZATION is correct but the port is
NOT bit-exact-or-cornered — the A/B surfaced THREE real divergences in jl's underlying snag LIST, so REVERTED per
"commit only validated":**
  1. **Snag ACCUMULATION**: per-year row counts — oracle GROWS 1993=5·2003=30·2013=14·2023=33·2033=47·2048=51 (ongoing
     mortality snags accumulating as distinct death-year cohorts); jl COLLAPSES to ~4/cycle after 2013 (1993=2·2003=15·
     2013=6·2018=4·…·2048=4). jl's ordinary-mortality snags don't PERSIST across cycles the way FVS's do (falldown too
     fast, or few distinct yrdead cohorts survive) ⇒ jl emits 68 rows vs oracle 379, missing 209/256 (sp,yr,yrdied) keys.
     NOTE the memory's "SnagSum bit-exact" was only the PULSE era (1993-2013, SNAGINIT+fire densities); the post-2013
     ACCUMULATION era (ongoing mortality) diverges — snag_summary likely diverges there too, just wasn't sampled.
  2. **Hard HEIGHT** (Current_Ht_Hard): 2/47 match — jl's htcur (SNAGBRK-reduced current height, FMSNGHT) ≠ FVS HTIH.
  3. **Hard VOLUME** (Current_Vol_Hard): 7/47 — the OC snag volume routine (still unresolved: NATCRS/NVEL vs CFVOL) +
     the height-loss feed both wrong. (Density_Soft/Ht_Soft/Vol_Soft trivially matched — all 0, snags are hard-only here.)
  Density_Hard/Total matched 29/47 on the keys jl HAS (the pulse-era cohorts) — so the aggregation logic is right; the
  DATA (which snags exist, their current heights) is what diverges. ⇒ **SnagDet is blocked on a jl snag-MODEL
  reconciliation** (ordinary-mortality snag persistence/falldown across cycles + the SNAGBRK height-loss model +
  the snag volume routine) — a multi-part model port, NOT the ~15-col serialization the earlier notes assumed. The
  serialization spec above is validated-correct and reusable once the snag list matches. Recommend: FIRST reconcile
  jl's multi-cycle ordinary-mortality snag accumulation vs FVS (fmsfall.f falldown + the per-cycle FMSDIT add) using
  the SnagSum GRAND-TOTAL density across ALL cycles (not just the pulse) as the cheaper 1-number-per-cycle A/B, THEN
  return to SnagDet. Reverted commit-free; tree clean, gate 339/11.

## ★★ SNAG-DENSITY BUG MEASURED 2026-08-21 — jl snag_summary is SYSTEMATICALLY LOW (memory "bit-exact" claim WAS FALSE)
Followed the recommended cheaper diagnostic: FVS_SnagSum grand-total density per cycle, jl vs ocsnag_oracle.db, ALL
12 cycles. jl UNDER-counts at EVERY cycle (Hard_snags_total, oracle→jl): 1993 55.3→14.8 · 1998 93.0→75.4 · 2003
264.9→194.1 · 2008 51.0→38.3 · 2013 8.7→4.7 · 2018 9.7→3.1 · 2023 9.2→1.8 · … · 2048 4.7→0.5 (late cycles ~10× low).
**This CORRECTS the memory/recipe claim that jl snag_summary matched "1993 55.26·2003 264.9·2008 51.0·2013 8.7" — those
were the ORACLE values; jl actually gives 14.8/194/38/4.7. The "bit-exact" was never a real head-to-head (a
measure-your-own-fix violation).** LOCALIZED: the gap is ENTIRELY in DBH-class 1 (small, <12"); classes ≥12" are
BIT-EXACT (1993 oracle 0.61 == jl 0.61 for class2-5). At cyc0 (1993, pre-falldown) the small-snag source is the
SNAGINIT cohort: keyfile `SNAGINIT 10 11 50 40 2 50` = species10, DBH11, HTDEAD50, HTIH40(pre-broken), age2,
density=50/ac. jl's ffe_add_snaginit! calls add_snag!(den=50) but the 1993 report shows only 14.76 (factor 0.295) —
so ~35/ac vanish between creation and the cycle-top report. **NOT tripling** (NOTRIPLE gives identical 14.76). 2 years
of normal falldown (modrate≈0.046/yr ⇒ ×0.91) can't explain 0.295 (that needs ~26 yr). ⇒ a REAL small-snag
creation/falldown bug in the SNAGINIT (or pre-broken-snag/height-loss) path — candidates: (a) jl over-applies
pre-inventory-age falldown while FVS treats the input density as the CURRENT (already-fallen) count; (b) the SNAGBRK
pre-broken (HTIH<HTDEAD) path drops density; (c) species-10 selector/filter drops the cohort so 14.76 is first-cycle
mortality only. NEXT DIAGNOSTIC (cheap, decisive): env-gate a print in ffe_add_snaginit! after add_snag! to confirm
den=50 is created, then a print in update_snags!/snag_summary at the 1993 report to see the density there — isolates
creation-vs-falldown. This bug gates BOTH FVS_SnagSum AND FVS_SnagDet (and feeds FFE Stand-Dead carbon/CWD), so fix it
BEFORE re-attempting the SnagDet serialization. Report-only subsystem; growth/.sum unaffected.

## SNAG-DENSITY bug SHARPENED 2026-08-21 (instrumented ffe_add_snaginit! + snag_summary, env SNAG_DEBUG)
Decisive: at the 1993 (inventory) snag report jl's list has **nrec=2, total 14.76, and the SNAGINIT cohort (den=50,
death-yr 1991) is NOT yet added** — the `SNAGINIT-add den=50` fires AFTER the 1993 report and the cohort first appears
at 1998 (nrec 2→8). ⇒ **jl adds SNAGINIT one cycle LATE (misses the inventory-year FFE snag report)** — a FFE-init
ORDER bug: ffe_add_snaginit! runs after the cyc0 carbon/snag report instead of before it (FVS fmsnag.f adds SNAGINIT
at the first FFE year, BEFORE that year's report). This is the cyc0 gap. The multi-cycle low tail (1998 75 vs 93; late
~10× low) is a SEPARATE falldown/accumulation divergence stacked on top (the 2 mystery cyc0 records + the fast late
decay). NOTE the arithmetic doesn't cleanly close (jl 14.76 non-SNAGINIT vs oracle 55.26≈50+5) — the "2 records" at
1993 need identifying (first-cycle mortality booked early?) before a fix. FIX PATH: move the SNAGINIT add ahead of the
inventory-year snag/carbon report (mind the carbon-push timing), then re-A/B SnagSum grand-total per cycle; expect the
1993 gap to close and the multi-cycle falldown to be the remaining item. Report-only; growth/.sum unaffected; gate 339/11.
