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
