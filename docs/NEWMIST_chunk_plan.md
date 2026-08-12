# NEWMIST (spatial dwarf mistletoe) port — chunk plan

**Status: PLANNED, not started — pending a user go/no-go (see #196 cost/benefit flag).**
This is the fix for the BC YSM multi-cycle under-mortalization (isolated 2026-08-12 to the missing
spatial DM via the DM-free `all_BC` self-thin control). It is a MAJOR port; this plan makes it
executable chunk-by-chunk (doctrine: one chunk/session, commit frequently, validate each chunk).

## Scope
`canada/newmist/` = **~8,761 lines / ~55 routines**, the FVS **spatial** dwarf-mistletoe model activated by
the `NEWSPRED` keyword (distinct from the non-spatial base `mistoe/mistoe.f` already ported for the N-Rockies
cluster via `_ie_mis_variant`). FVSbc links the newmist `mistoe.f` (confirmed: `diff` == `canada/newmist/mistoe.f`).
Currently jl does NOT parse MISTOE/NEWSPRED/DMAUTO (→ `unrecognized_keywords`) and BC is not in any DM dispatch,
so jl runs ZERO DM on these stands.

## Architecture (measured)
- Per-cycle entry: `gradd.f:96 CALL MISTOE` (same site as the base model).
- `MISTOE` (newmist, mistoe.f): `OPFIND` scheduled actions (DMAUTO etc.) → `DMTREG` (spread/intensify) →
  `DMMDMR` (stand-mean DMR) → `RANN` draws (spatial spread — ZZRAN/RNG-straddle class) → `MISINF` (infection)
  → `IF(DMFLAG) MISMRT` (DM-caused mortality). DG-loss via `misdgf`/`dgdriv`.
- Key routines: dmcycl(593 scheduler-driver) / mistoe(558 spread entry) / dmtreg(553) / dmshap(378 crown shape) /
  dmadlv(349 adjust DM levels) / **bndist(342 between-tree distance — the core novel spatial piece)** /
  dminitbc(284 BC init) / dmblkd(256 coefficients block-data) / dmauto(149 DMAUTO keyword) / dmopts(157 keywords) /
  12× dmcw*(148 ea, per-variant crown-width).

## ⚠ Source-read gotchas (found while surveying — save the next session the mis-step)
- **`dmauto.f` is NOT DMR-seeding.** It is the spatial-autocorrelation REWEIGHTING function
  `f = EXP(A·|DMR_s−DMR_t|·EXP(B·D))` over the `SF(DMRDFF(i,TrgDMR), RQ)` matrix (sampling ring `RQ`), called by
  `DMTREG` to redistribute per-DMR-class neighbour densities. The **DMAUTO keyword** just sets its A/B coefficients
  (the `-0.50`). Initial DMR assignment is elsewhere (dminitbc / misintbc / DB damage codes) — find it before chunk 1.
- Core spatial state lives in `DMCOM.F77` (SF matrix, DMRDFF delta-DMR index, sampling-ring geometry MESH). Port the
  common block's meaning, not just arrays.

## Chunks (each: port faithfully → validate → commit)
1. **Keywords + DMR init.** Parse MISTOE/NEWSPRED/DMAUTO/MISTPRT (dmopts) — DMAUTO sets the SF A/B coefs (NOT DMR).
   Find + port the actual initial-DMR assignment (dminitbc.f / misintbc.f / DB damage codes) + coefficients
   (dmblkd.f). `t.dmr` already exists. Wire into keyword_dispatch. NOT `.sum`-validatable alone (DMR seeded but no
   effect yet) — validate the seeded DMR values vs an instrumented newmist dump.
2. **Spatial distance core.** `bndist.f` (between-tree distance) + tree-position/grid handling. The novel piece
   the base model lacks. Unit-test the distance math against instrumented newmist values.
3. **Crown width.** The per-variant `dmcw*` (BC's — verify which; dminitbc may select). Feeds spread.
4. **Spread / intensify.** `DMTREG` + `dmshap` + `dmadlv` + the RANN spread draws + `DMMDMR`/`MISINF`. RNG order
   MUST match live (ZZRAN discipline — never FFI the RNG; derive draw order from mistoe.f).
5. **DG-loss.** `misdgf` DMR→DG multiplier. Check reuse of the ported `ie_dm_dg_mult` vs newmist-specific coef.
6. **DM mortality.** `MISMRT` (mismrt.f) — combine into the BC mortality driver (like `ie_dm_mortality_combine!`).
7. **Output.** dmsum/misprt/dmtlst/dmslst DM reports (optional; `.sum`-inert).
8. **BC wiring + end-to-end.** Dispatch newmist for BC; end-to-end validate.

## Validation
- **Primary target EXISTS** (no DB-capable oracle needed): `tests/FVSbc/YSM-SkyRanch.sum.save` (2022 production
  FVSbc output, WITH MISTOE/NEWSPRED/DMAUTO). Run jl on the YSM stand (units now correct after 42eb555) and diff
  the multi-cycle `.sum` — the target is jl's TPA declining toward the oracle's harder mortality (2077 ~1366).
- **Regression control:** the DM-free `mrun/all_BC.key` (V3, no MISTOE) — jl already MATCHES it (TPA 2089→1253 vs
  2087→1292); the newmist port must leave it BYTE-IDENTICAL (DMR=0 ⇒ inert, like the base model).
- **A/B (optional, currently blocked):** `FVSbc_clean` WITH vs WITHOUT MISTOE would isolate the DM contribution,
  but the relinked oracle SIGSEGVs on ALL DATABASE reads (isoc23-shim × SQLite-C-interop — a relink-infra issue,
  NOT a production-FVS bug; the 2022 oracle read these DBs fine). Either fix the relink (DB-capable FVSbc) or
  convert YSM to inline TREEDATA. Not required — the `.sum.save` is the trustworthy target.

  **DB-capable relink attempts (2026-08-12, both FAILED):** removing the isoc23-shim from relink_bc.sh still
  SIGSEGVs on DB reads (the shim only provides `__isoc23_sscanf`, resolved from libc without it) ⇒ NOT the shim.
  The crash is a deeper gfortran-16 C/Fortran-interop issue (ISO_C_BINDING × SQLite fvsqlite3.c binding, or the
  pre-compiled sqlite3.o ABI). A DB-capable FVSbc relink would need recompiling sqlite3.c/fvsqlite3.c fresh under
  the current toolchain and/or auditing the fsql3 character-passing — deferred (marginal value vs the .sum.save
  target). ⇒ Use the INLINE-TREEDATA route for the A/B if needed, not a DB relink.

## Cost/benefit (why this is a user decision)
~8,800-line port of an OPTIONAL keyword model whose only exercised corpus case is the single YSM DATABASE stand.
The base (non-spatial) mistoe.f is done+validated for the N-Rockies cluster. Weigh before committing a
multi-session effort. If greenlit, execute the chunks above in fresh sessions.
