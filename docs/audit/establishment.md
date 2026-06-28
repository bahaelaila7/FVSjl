# Bandaid Audit — Establishment / Regen (`src/engine/establishment.jl`)

Module: ESTAB tree-creation + ESNUTR cycle hook + plot replication (ESSUBH/ESGENT/REGENT-LESTB).
FVS source checked: `bin/FVSsn_buildDir/{estab.f, esnutr.f, esinit.f, esplt2.f, essubh.f, esgent.f, regent.f, blkdat.f}`.

Scope note: the *only* validated scenario (per project memory) is the snt01 **BARE** stand —
`NPTIDS=1`, NATURAL keyword, 5-yr cycle (`FINT=5`), seedling heights floored to `XMIN`. Several
items below are faithful *for that path* but mishandle untested cases; one has no source basis at all.

Faithful items confirmed (counted, not listed): `MINREP=50` (esinit.f:55), the `_ES_HHTMAX` table
(blkdat.f:86-95, matches element-for-element), `dupnpt = nptids·idup` (estab.f:131), the per-replicate
ESRANN reseed `es0 = esdraw|1` (estab.f:528 ESRNSD), the crown coefficients
`0.89722 / 0.0000461 / 0.07985` and the `[0.20,0.90]` clamp + `INT(cr*100+0.5)` using the MAIN RANN
stream (regent.f:178-185), `ibrkup = INT(ptree/10+1)` and `tpa = ptree/brk` (estab.f:589-624),
the species-order `newidx` sort and `sort_key=i` lineage reset (esgent.f:43-44 SPESRT), the height
`else`-branch `BACHLO(0.5,0.25)` loop + `XMIN` floor (estab.f:485-489), and the `idsdat` default
`IY(1)-20` (behavior matches esnutr.f:100 — but the docstring's "(esnutr.f:63)" citation is wrong;
line 63 is the SPECMULT loop, the real default is line 100; behavior is faithful, citation is stale).

---

## FLAG 1 — BANDAID — `gentim` has no FVS basis (establishment.jl:61)

**jl:** `gentim = (Int(yr) + per - idsdat) - per; gentim < 0 && (gentim = 0)`
which algebraically reduces to `gentim = max(yr - idsdat, 0)`, then feeds
`age = per - delay - gentim + trage` (line 81), the age handed to `htcalc_height`.

**Comment claims:** "gentim/delay/trage timing (esnutr/estab/essubh): age = FINT − delay − gentim + trage."

**FVS source — estab.f:317-318:** `GENTIM = FINT-5.0 ; IF (GENTIM.LT.0.0) GENTIM=0.0`.
The per-keyword reset (estab.f:508-512) is `GENTIM = FINT-DELAY-5.0` (else 0). ESSUBH then computes
`AGE = TIME-DELAY-GENTIM+TRAGE` with `TIME=FINT` (estab.f:469, essubh.f:59). **GENTIM depends only on
`FINT` (and DELAY), never on `IDSDAT` or the calendar year.** The jl's `yr - idsdat` is a different
quantity entirely: for the validated BARE stand `idsdat = IY(1)-20 = 1992` (constant) while `yr`
increases each cycle, so jl `gentim` = 20, 25, 30, … whereas the true GENTIM = 0 every cycle (FINT=5).

**Why tests still pass:** EMSQR is unused in ESSUBH (essubh.f:46 `RDANUW=EMSQR`) and for NATURAL
seedlings `htcalc_height(age)` for *any* small age is below `XMIN`, so the wrong `age` is erased by
the `XMIN` floor at line 97. The expression is effectively inert in the tested scenario.

**Severity: BANDAID.** Source explicitly contradicts it. On any path where `htcalc_height(age) > XMIN`
(better sites/species, or `FINT≠5` where the true GENTIM=FINT-5≠0) the regen heights — and thus DBH,
BA, density — will diverge. Faithful replacement: `gentim = max(per - 5, 0)`.

---

## FLAG 2 — GAP — `idup` uses floor where FVS uses ceil (establishment.jl:55)

**jl:** `idup = max(1, fld(_ES_MINREP, nptids))` — floor division.

**FVS source — estab.f:122-130:** the replication loop increments `DUP` until `N = NPTIDS*I ≥ MINREP`
(or `NPTIDS*(I+1) > MAXPLT`); `IDUP = INT(DUP+0.5)`. The exit condition `NPTIDS*I ≥ MINREP` gives the
**smallest** I with `NPTIDS·I ≥ 50`, i.e. `idup = ceil(50/NPTIDS)`, **not** `floor`.

For `NPTIDS=1` (the only tested value) `ceil = floor = 50`, so it matches. For `NPTIDS ∈ {3,4,6,7,…}`
they differ (e.g. NPTIDS=3 → FVS 17, jl 16; NPTIDS=4 → FVS 13, jl 12). `idup` drives `dupnpt`
(the planted-TPA divisor) and the per-replicate record/RNG-draw count, so a wrong `idup` shifts both
the regen record count and the downstream MAIN-RANN crown draws.

**Severity: GAP** — correct only when `NPTIDS` divides 50 (incl. the tested `NPTIDS=1`). Fix: use
`cld(_ES_MINREP, nptids)` (ceil).

Related (same line context, establishment.jl:54): `nptids = IPTINV − NONSTK` is the **no-treedata**
branch only (esplt2.f:73-75, taken when `IPTKNT=0`). For establishment *into a stand that already has
tree records* (`IPTKNT>0`), `NPTIDS` is the count of plot-ID-matched stockable plots (esplt2.f:90-110),
which need not equal `IPTINV−NONSTK`. Faithful for the BARE stand; a GAP for regen-under-canopy.

---

## FLAG 3 — GAP — ESGENT height-growth + REGENT DBH scaling omitted (establishment.jl:102-103, and absence of any first-cycle growth)

**jl:** new trees are added at cycle end and DBH is `dbh = _htdbh_dbh(sd,sp,hht,ifor)` (from the
*establishment* height `hht`), floored to `0.1`, plus `0.001·hht`; the tree is **not** grown this cycle.

**Comment claims:** "REGENT establishment dbh (regent.f:331-334, LESTB branch): DBH = HTDBH⁻¹(**HK**),
floored to the species min **DIAM**, then … 0.001·**HK**."

**FVS source — regent.f:331-335 + esgent.f:48-66:** ESTAB calls `ESGENT → REGENT(.TRUE.)` which **grows
the new trees to end-of-cycle** before assigning DBH. REGENT sets `DBH = DKK = HTDBH⁻¹(HK)` where
`HK = HT+HTG` (the *grown* height, regent.f:284-300), floors to `DIAM(ISPC)` (species budwidth),
adds `0.001·HK`. ESGENT then applies the WK4/HTIMLT scaling (esgent.f:51-62):
`HT = HHT + HTG·WK4`; if `WK4<1` and `HT<4.5` → `DBH = 0.1+0.001·HT`, else `DBH = DBH·(HT/HTEMP)`.

The jl implements **only** the `WK4=0, HT<4.5` corner: for `FINT=5`, `GENTIM=0 ⇒ WK4=HTIMLT=0`
(estab.f:516 `FTEMP/(GENTIM+0.0001)` with `FTEMP=min(TRAGE,0)=0`), so the trees genuinely do not grow
and `DBH = 0.1+0.001·HHT` — which equals the jl's `(_htdbh_dbh floored to 0.1) + 0.001·hht`. That is
why the comment cites `HK`/`DIAM` but the code uses `hht`/`0.1`: they coincide *only* for FINT=5
sub-4.5-ft seedlings.

They diverge when (a) `FINT≠5` → `WK4≠0` → FVS grows the tree by `HTG·WK4` this cycle and scales DBH
(jl does neither), or (b) the established height ≥ 4.5 ft — reachable because several species are
floored to `XMIN ≥ 4.5` (sp 11=5.05, 13=4.70, 38=5.98, 53=4.15, 80=5.98, 81=4.70; blkdat.f:67-76) or a
`PLANT` keyword supplies a tall `TREEHT`. There FVS uses `DKK·(HT/HTEMP)` and the `DIAM(ISPC)` floor,
not `_htdbh_dbh(hht)+0.001·hht` floored to 0.1.

**Severity: GAP.** Faithful for FINT=5 / seedlings <4.5 ft (the tested hickory-dominated BARE regen,
all `XMIN<4.5`); wrong DBH/HT for non-5-yr cycles and for species/PLANTs whose established height ≥4.5 ft.
(Note `ABIRTH = AGEPL+GENTIM`, estab.f:628/707, is also not tracked — a reporting/age GAP.)

---

## FLAG 4 — GAP — PLANT-height branch floors to XMIN, not 0.05 (establishment.jl:97)

**jl:** after both height branches, line 97 applies `hht < es_xmin[sp] && (hht = es_xmin[sp])`
**unconditionally**.

**FVS source — estab.f:475-492:** the two branches floor differently. The **TREEHT≥0.1** branch
(`PLANT` with an explicit height) floors at **0.05** (estab.f:483 `IF(HHT.LT.0.05) HHT=0.05`) — it does
**not** apply the `XMIN` floor. Only the `else`/default branch floors at `XMIN` (estab.f:489). The jl
collapses both into a single `XMIN` floor.

Consequence: a `PLANT` keyword that specifies a small height (say 1 ft) for a high-`XMIN` species
(e.g. sp 38, `XMIN=5.98`) keeps ≈1 ft in FVS but is bumped to 5.98 ft by the jl.

**Severity: GAP** — never exercised by the NATURAL-only BARE test (`treeht<0.1`, else-branch). Also in
this block the jl comment "HTADJ=0" silently drops the `+HTADJ(IPNSPE)` add present in *both* FVS
branches (estab.f:482,488); harmless unless the `HTADJ` (act 442) keyword is used, which
`establishment.jl` does not process at all.

---

## FLAG 5 — GAP — `pccf` hardcoded to 0 (establishment.jl:62, 139)

**jl:** `pccf = 0f0  # point crown competition factor (≈0 for the sparse established plots)` and the
crown formula uses it: `cr = clamp(0.89722 - 0.0000461*pccf + 0.07985*ran_cr, …)`.

**FVS source — regent.f:178:** `CR = 0.89722 - 0.0000461*PCCF(IPCCF)` where `PCCF(IPCCF)` is the
**point** crown-competition factor of the inventory point the tree sits on. For a truly bare plot
`PCCF≈0`, so the coefficient term vanishes and the jl matches. For establishment **into a stocked
stand** (NATURAL/PLANT under an existing canopy) `PCCF>0`, lowering the open-grown crown ratio by
`0.0000461·PCCF` (≈0.01–0.02 at PCCF 200–400) and shifting `ICR` by ~1–2 — a real, source-mandated
adjustment the jl omits.

**Severity: GAP** — justification is an explicit "≈0 … assumption". Faithful for the bare plot, wrong
crown ratio for under-canopy regen.

---

## FLAG 6 — UNVERIFIED — missing WK6 site-prep ESRANN draws (establishment.jl:70-74)

**jl:** per replicate consumes exactly three ESRANN draws (two EMSQR + one `esdraw`/ESAVE) and the
height `BACHLO` draws.

**FVS source — estab.f:202-205:** **before** the plot loop, FVS unconditionally consumes
`IDUP·NPTIDS` ESRANN draws into `WK6` (site-prep sampling), and (NTALLY=1) one extra `ESDRAW` draw +
`ESRNSD` reseed (estab.f:174-179). For `NPTIDS=1, IDUP=50` that is 50 draws the jl never makes. Because
each replicate is reseeded from the *previous* replicate's `ESAVE` (estab.f:528 / jl:124), only the
**first** replicate's draws depend on the post-WK6 generator state — but they do depend on it.

I could not confirm from this file alone whether the jl seeds `s.rng.es0` to FVS's post-WK6/ESDRAW
state elsewhere (esinit/esplt setup). In the tested scenario it is masked anyway: the ESRANN stream
only drives EMSQR (unused, essubh.f:46) and the height `BACHLO` draws, which are erased by the `XMIN`
floor. **Severity: UNVERIFIED** — a latent ESRANN-stream misalignment that would surface only when
established heights are *not* XMIN-floored. To confirm/refute: trace `s.rng.es0` initialization at the
cycle's establishment entry and compare to FVS's `ESRNSD(.TRUE.,ESDRAW)` + 50 `WK6` advances.
