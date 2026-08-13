# BC NEWSPRED (spatial dwarf-mistletoe / NISI) port plan — #196

USER-greenlit 2026-08-13 ("scope then start porting"). Closes the BC YSM multi-cycle
under-mortalization (jl SDI→1586 vs oracle 926) attributed to the missing spatial DM.

## Model overview
`canada/newmist` = the **NISI** (New & Improved Spread & Intensification) *spatial*
dwarf-mistletoe model. 8,761 lines / 50 `.f` files. Distinct from the base `mistoe.f`
(non-spatial, already ported for N-Rockies via `_ie_mis_variant` + CR's own
`dwarf_mistletoe_model.jl`). BC currently has NO DM model wired.

**What it does:** per crown-third of each tree record it maintains 4 life-history
compartments (IMMATURE / LATENT / NON-FLOWERING / FLOWERING). Infections spread
spatially between neighbouring trees (between-tree distance × crown shape × slope ×
spatial autocorrelation), intensify over cycles, and produce a per-tree DMR (dwarf-
mistletoe rating 0-6). DMR drives HEIGHT + DIAMETER growth multipliers (BHTG/DHTG,
mistoe.f) and extra MORTALITY (MISMRT) — the latter is the YSM gap.

## Call graph (entry → payoff)
```
MISIN0 → DMINIT (dminitbc.f)                      # one-time init of DMCOM commons
FVS cycle → MISTOE (mistoe.f)                     # per-cycle DM driver
   ├─ DMTREG (dmtreg.f)                           # NISI spread & intensification core
   │    ├─ DMOPTS (keywords DMCLMP/DMALPH/DMBETA) → DMNAC
   │    ├─ DMMTRX / DMFBRK / DMFSHD / DMADLV      # spread matrix, bark, shade, adult-lives
   │    ├─ DMNB / DMNTRD / DMNDMR / DMSRC         # neighbour build, target reduction, DMR, source
   │    ├─ DMSLST / DMSLOP / DMSAMP / DMOTHR / DMTLST
   │    └─ BNDIST (bndist.f) + DMSHAP + DMCW*     # spatial substrate (distance/shape/crown-width)
   ├─ DMMDMR (dmmdmr.f)  → per-tree DMR
   ├─ DMCYCL (dmcycl.f)  → life-history compartment advance (immature→latent→…→flowering)
   ├─ DMAUTO (dmauto.f)  → SF spatial-autocorrelation reweighting matrix
   ├─ MISINF (infection) + MISMRT (mortality, USEMRT)   # ← PAYOFF: extra mortality
   └─ HTG/DG multipliers (BHTG/DHTG in mistoe.f)        # ← PAYOFF: growth reduction
Reports: DMSUM / DMTLST / DMSLST / misprt.f
```
Keywords (currently unrecognized in jl): MISTOE, NEWSPRED, DMAUTO, DMCLMP, DMALPH, DMBETA.

## Chunk decomposition (dependency-ordered)
- **C0 — Foundation + validation harness.** Port `DMCOM.F77` common → a Julia
  `MistletoeState` (per-crown-third 4-compartment arrays, per-tree DMR, SF matrix,
  RNG seed `DMRNSD`). Build the VALIDATION VEHICLE (see below) so every later chunk
  validates vs a live oracle. NO model logic yet — just state + harness + a loud
  "NEWSPRED active but unported" error so partial states never silently mis-run
  (doctrine #5).
- **C1 — Keywords + DMINIT.** Parse MISTOE/NEWSPRED/DMAUTO/DMCLMP/DMALPH/DMBETA
  (`dmopts.f`, `dmauto.f` keyword halves) + wire `DMINIT` (dminitbc.f) + `DMRNSD`
  (dmrann.f RNG seed). Validate: keyword echo + init state vs FVSbc DEBUG.
- **C2 — Coefficients + BC crown width.** Port `dmblkd.f` (256-line DATA) + the BC
  crown-width routine (dmcw.f dispatch → the BC-applicable dmcw*). Extract to CSV,
  verify line-for-line.
- **C3 — Spatial substrate.** `bndist.f` (between-tree distance, 342), `dmshap.f`
  (crown shape, 378), `dmslop.f` (slope). Validate per-tree distances/shapes vs g16.
- **C4 — NISI spread & intensification core.** `dmtreg.f` (553) + its helpers
  (dmmtrx/dmfbrk/dmfshd/dmadlv/dmnb/dmntrd/dmndmr/dmsrc/dmslst/dmsamp/dmothr/dmtlst).
  Largest chunk; likely split into C4a (matrix+neighbour build) / C4b (spread+DMR).
- **C5 — Life-history + autocorrelation.** `dmcycl.f` (593, compartment advance) +
  `dmauto.f` (SF reweighting). Validate compartment counts per crown-third vs g16.
- **C6 — PAYOFF: stand coupling.** `mistoe.f` driver + `dmmdmr.f` (DMR) + `misinf`
  (infection) + `mismrt` (MISMRT mortality, USEMRT) + BHTG/DHTG growth multipliers.
  Wire into the BC growth+mortality dispatch. THIS closes the YSM gap. Validate the
  YSM .sum TPA/SDI vs oracle.
- **C7 — Reports.** `dmsum.f` / `dmtlst.f` / `dmslst.f` / `misprt.f` DM summary +
  treelist reports (validation-supporting, low-risk).

## Validation strategy — unblocking the crash-blocked A/B
The clean A/B is blocked because the gfortran-16-relinked `FVSbc_clean` SIGSEGVs on
ALL DATABASE reads (isoc23-shim × SQLite C-interop, dbstreesin.f) — NOT a production
bug. Three vehicles, in priority order:
1. **INLINE TREEDATA reproducer (primary).** Extract YSM's trees from the FVS-ready
   DB into an inline `.tre` + a BC keyfile (STDINFO BEC string) with NEWSPRED/MISTOE
   enabled. The INLINE path in FVSbc_clean WORKS (all_BC/all_BC_essf run inline). This
   yields the oracle `.sum` WITH DM mortality → the C6 payoff target. Build in C0.
2. **FVSbc_g16 instrument-replay (per-chunk).** Full gfortran-16 rebuild of
   `bin/FVSbc_buildDir/*.f` (+ newmist) with the isoc23 shim = the durable
   instrumentable oracle (same recipe as FVS{ut,tt,em}_g16). Add DEBUG/WRITE dumps to
   dmtreg/dmcycl/bndist for per-chunk per-tree validation (DMR, distances, compartments).
   Confirm it reads inline stands (avoids the DB path entirely).
3. **DB-capable relink (fallback).** If an inline reproducer can't capture the YSM
   spatial config, relink FVSbc against a non-shim SQLite to unblock the DB path — only
   if 1+2 prove insufficient.

## Integration notes
- Existing jl base DM: `centralrockies/dwarf_mistletoe_model.jl` (CR) +
  `inlandempire/mistoe_coefficients.jl` (`_ie_mis_variant`, N-Rockies). NEWSPRED is a
  SEPARATE spatial model — do NOT overload the base path; add a `bc_newspred!` dispatch
  gated on the NEWSPRED keyword, wired into the BC cycle only (BC-first; other variants
  have their own dmcw* but are out of scope for #196).
- Doctrine: validate bit-exact-or-cornered per chunk vs the inline oracle / g16; MEASURE
  don't infer; never FFI the RNG (DMRNSD/DMRANN is its own stream — expect a realization
  straddle on the spatial draws, same accepted class as ZZRAN).
- Risk: the spatial neighbour model uses per-tree (x,y) positions FVS synthesizes from
  point/plot layout — confirm jl reproduces the same synthetic positions (bndist input)
  before trusting spread counts; this is the likely first "cornered vs real" fork.

## Status
- **C0 validation foundation DONE 2026-08-13.** Key simplification: jl reads the YSM
  metric DB directly (SQLite.jl + the 42eb555 metric-ingest fix) — NO inline reproducer
  needed for the PAYOFF validation. Set up `/workspace/.bcwork/newspred/ysm271.key`
  (YSM029-271 from FVS-BC.YSM-SkyRanch.db, NOTRIPLE, TIMEINT 5/first-4, NUMCYCLE 24).
  jl(DM-free) vs oracle `YSM-SkyRanch.sum.save`:
    - **cyc0 (2018) BIT-EXACT** (jl 2300/10/311 vs oracle 2300/9/311) → stand loads +
      growth correct; the DM effect is the ENTIRE divergence.
    - **C6 target quantified** — NEWSPRED must kill trees AND suppress growth: by 2077
      jl 1713/79/1586 (TPA/BA/SDI) → oracle 1366/42/926 (−20% TPA, −47% BA, −42% SDI).
  Payoff-validation vehicle = jl DB-run vs `.sum.save` (works today). Inline FVSbc
  reproducer for per-chunk internals DEFERRED to C4 (only chunk likely to need it).
- **C1 (state + keywords) DONE 2026-08-13** (9ae3880 + 3d5644d): MistletoeState + DMINIT
  default tables (DM_DMDMR/DM_OPAQ); StandState.mistletoe field; kw_mistoe!/kw_newspred!/
  kw_dmauto!/kw_mistprt! recognize the 4 DM keywords (BC-only) and set active/newmod/
  dmalpha/prtmis. VALIDATED on YSM029-271: state correct (dmalpha=−0.5), keywords out of
  unrecognized_keywords, .sum data rows byte-identical (inert).
- **C1 remaining piece — DMINIT infection seeding (SCOPED, not yet done):** misdam.f — DM
  damage codes 30-34 (30 generic / 31 LP / 32 WL / 33 DF / 34 PP); `IMIST/DMRATE = severity`
  (capped 0-6). YSM029-271 = ONE initial infection (Pl, Damage1=31/Sev1=2 → DMR 2); the
  other 38 start clean (infect via spread). PREREQUISITE: jl's DB reader (fia_database.jl
  load_fia_stand!) does NOT currently load FVS_TreeInit Damage1/Severity1/Damage2/… — must
  add that (additive, default 0 when absent) before seeding `dmr[i]` + `dminf[i,ct,pool]`
  via DM_DMDMR. Then DMRNSD (dmrann.f) RNG seed. Careful change to the load-bearing DB
  reader — do next, validate the DMR-2 tree loads + seeding matches DMINIT.
- **C2 (coefficients + crown-width) — RESOLVED-MINIMAL 2026-08-13, mostly folded into C1/C4:**
  - dmblkd.f is dominated by **Shd1** (1496-element encoded spread-field TRAJECTORY table) +
    ShdPtr — that is the spatial-spread substrate consumed by dmfshd/dmtreg, so it moves to
    **C4** (extract with the spread core, where it's used). The genuinely-simple coefficients
    (DMDMR crown-third distribution, DMOPAQ species opacity) already landed in C1.
  - **BC crown-width = NO new model.** dmcw.f is just `DMTRCW(I) = CRWDTH(I)` — it returns
    FVS's existing per-tree crown width, which jl already computes (`t.crown_width`). No
    per-variant dmcw* for BC (BC isn't in the dmcw* list); the DM spread geometry reads the
    engine's crown width directly. Nothing to port here.
  ⇒ C2 has no standalone deliverable; its content is absorbed by C1 (coeffs) + C4 (Shd1).
- **NEXT: C3** — spatial substrate: dmshap.f (crown shape, 378), dmslop.f (slope, 69). ⚠ NOTE:
  bndist.f is NOT "between-tree distance" (a mislabel from the goal-doc survey) — `BNDIST` is the
  Binomial/Poisson/Negative-Binomial PDF family for neighbour counts (V≈M→Poisson, V<M→Binomial,
  V>M→NegBinom), called by DMNB; it belongs with the neighbour-density math in C4, not the geometry.
- **C4 — spread core (STARTED 2026-08-13).** DMTREG driver FULLY MAPPED (dmtreg.f:239-521):
  setup (DMOPTS/DMMTRX/DMFBRK/DMFSHD/DMNTRD/DMFINF) → per-ring DMNB neighbour build →
  per-species loop {DMFDNS density, DMSRC source vector, per-target-DMR-class k {DMTLST
  targets, per-ring m {DMAUTO autocorr source density, per-source-DMR-class n, per-target u
  {DMSAMP count, DMSLST select sources, per-source v: MISDGF dgf, DMSLOP slope offset,
  DMADLV spread-field accumulation over MESH bands (uses Shd1 shade), ATAN subtended-angle
  interception → NewSpr}}}} → self-intensification NewInt → DMOTHR apply} → DMCYCL life-hist
  → DMNDMR recompute DMR → BrkPnt→PBrkPt. RNG-heavy sampling (DMSAMP/DMSLST/DMSRC) ⇒ the
  spatial spread is a REALIZATION straddle; validate at the aggregate FVS_DM_Stnd_Sum_Metric
  (Mean_DMR / Inf_TPH / Mort_TPH trajectory), NOT bit-exact per-tree.
  ★ DETERMINISTIC GEOMETRY SUBSTRATE COMPLETE (ae259b4/b30fe49/0702c4e) — the per-tree geometry
  every spread loop reads, all engine-inert + precompile-clean:
    - **DMFBRK** — crown-third breakpoints in MESH units + DMCOM constants (MESH=2, FPM, MXHT=25,
      MXTHRX=7 rings, ORIGIN, TWOPIE, BPCNT=4).
    - **DMSHAP** — per-tree crown shape (5 shapes) via Fisher discriminant; 7 coeff tables (5×11)
      + BC species→group MAPBC, extracted programmatically.
    - **DMRDMX** — per-MESH-band frustum VOLUME/RADIUS branching on shape (dmsum.f:96-256);
      VALIDATED: cone frustum sum == analytic cone volume (523.6 ft³).
  ★ STRUCTURE CORRECTION (measured 2026-08-13): **DMFSHD is STOCHASTIC**, not part of the
  deterministic substrate — it simulates each tree's (x,y) on a 121×121 grid via DMRANN (Poisson)
  to build the shade/light field, so it belongs to the RNG layer (expect a realization straddle).
  And **Shd1/ShdPtr** (the INTEGER*2 1496-elt encoded seed-TRAJECTORY "black box", COMMON /DMMIST/)
  is consumed by **DMADLV** during spread accumulation — NOT by DMFSHD. So the Shd1 extract pairs
  with DMADLV, not the shade field.
  NEXT (the RNG-coupled spread core, validate at aggregate FVS_DM_Stnd_Sum_Metric, not per-tree):
  DMFSHD (grid shade sim) + DMNB/BNDIST (neighbour-count PDF) → DMSRC/DMSAMP/DMSLST (source sampling)
  + DMADLV (spread accumulation, needs the Shd1/ShdPtr extract) + DMAUTO (autocorr) → C5 DMCYCL/DMNDMR.
- **C4/C5 support layer LANDED + unit-validated (2026-08-13)** — the deterministic + RNG-support pieces the
  spread loop reads, all engine-inert:
    - Geometry: DMFBRK, DMSHAP, DMRDMX (validated: cone frustum = analytic vol).
    - Neighbour PDF/CDF: BNDIST+GAMMLN (validated vs Poisson/NegBinom/Binomial), DMNB annulus CDF (+ CrArea/Dstnce).
    - RNG: DMRANN (MINSTD, validated) + DMRNSD seed + DMSLOP slope offset.
    - DMR recompute: DMNDMR (stochastic crown-third rating).
    - Autocorrelation: DMAUTO + SF matrix (validated: SF=exp(-0.5·diff), density-preserving).
  DMMDMR = trivial per-tree DMR array copy (New→Old), folded in when wired — no standalone port.
- **NEXT — the coupled SPREAD BLOCK** (needs a coherent DMSPtr/DMSInd treelist-DM-index design first, built by
  DMFINF): DMFINF (index, 140 lines, det) → DMFDNS (target density, 93, det) / DMSRC (source vector, 136, det) /
  DMTLST (target list, 114, 1 RNG) → DMSAMP (109, 3 RNG) / DMSLST (160, 2 RNG) source sampling → DMADLV (349, spread
  accumulation, needs the Shd1/ShdPtr 1496-elt extract) → DMOTHR (apply, 114, det).
  ★ BACKBONE DESIGN (mapped 2026-08-13): DMFINF builds `Ptr(MAXSP, 0:6, {FST,LST})` (per species×DMR-class, the
  first/last positions) + `Index` (a treelist permutation) via **OPSORT** (base/opsort.f) — an UNSTABLE two-key
  QUICKSORT (primary ISP species, secondary DMRATE DMR, ascending). SAME class as RDPSRT: a jl STABLE sort would
  mis-order equal-(sp,DMR) ties. Replicate OPSORT's quicksort faithfully to keep the straddle minimal; but since the
  spread is already an RNG realization straddle and OPSORT ties only feed the RNG sampling, exact tie-order is NOT
  bit-critical (part of the accepted straddle). Everything downstream indexes trees through this Ptr/Index. Plus DMFSHD (231, 4 RNG, grid
  shade), DMMTRX (73, det, spread-matrix setup), DMNTRD (202, cycle-follow remap). Then C5 DMCYCL (593, det,
  life-history compartment advance). Then wire the DMTREG driver + C6 payoff (mistoe/mismrt + BHTG/DHTG). Validate
  the coupled block at the aggregate FVS_DM_Stnd_Sum_Metric trajectory (spatial draws = realization straddle).
- **★ SUPPORT LAYER COMPLETE (2026-08-13) — 17 routines landed + unit-validated, all engine-inert:**
  geometry DMFBRK/DMSHAP/DMRDMX; neighbour BNDIST+GAMMLN/DMNB; RNG DMRANN/DMSLOP; DMR DMNDMR;
  autocorr DMAUTO+SF; index opsort!/DMFINF/DMFDNS/DMSRC/DMTLST; sampling DMSAMP/DMSLST; trajectory
  Shd1/ShdPtr extract + DMBSHD decode. Each validated against a known quantity (analytic vol / known
  distributions / MINSTD / decode-vs-table / density sums).
- **★ DMADLV FULLY MAPPED (dmadlv.f, the hardest routine — port next):** walks each DMBSHD-decoded seed
  trajectory from an infected source. Per step: h = MshHt + CShd[k,m,ZZ] − ORIGIN (z-band); if inside the
  source crown (x≤Rad=DMRDMX[src,h,RADIUS]) → intensification `IFld[h] += VecWt·Op` (self); at the target
  distance (CShd[k,m,XX]==Dist) → spread `SFld[h] += Cnt·VecWt·Op` (to target); else en-route shading loss
  `VecWt −= VecWt·Shade[h]`. Op = DMOPQ2[sp] = 1−(1−DM_OPAQ[sp])^MESH (per-MESH-cell opacity, dmtreg.f:229).
  In the DMTREG driver call II=0 ⇒ the Shd/Shd0 edge-buffer distinction collapses to Shade[h]. Needs the
  DMFSHD `Shade` per-band field (stochastic grid sim, the one remaining input) + DMRDMX (done).
- **REMAINING = the ASSEMBLY + PAYOFF phase (a coherent block, validate together at the aggregate DMR
  trajectory / YSM mortality — NOT independently):** DMFSHD (shade grid, 231, 4 RNG) → DMADLV (accumulation)
  → DMOTHR (apply new spread/intensification, 114, det) → DMCYCL (life-history compartment advance, 593, det)
  → wire the DMTREG driver (the mapped loop) → DMMTRX/DMNTRD setup/remap → C6 payoff (mistoe/mismrt USEMRT +
  BHTG/DHTG growth-mult, wire BC into the DM dispatch). THIS closes the YSM gap (jl SDI→1586 vs oracle 926).
- **★ DMADLV LANDED (bf8ca54) — the spread/intensification accumulation CORE, hardest routine done.**
  Validated qualitatively (radius routes weight to SFld-spread vs IFld-intensification). ⇒ the ENTIRE NEWSPRED
  SPREAD MECHANISM is now ported (18 routines: full support/index/sampling/decode layer + DMADLV).
- **REMAINING = the final assembly + payoff (a coherent block, validate together at YSM DMR/mortality):**
    - DMFSHD (shade grid, 231, 4 RNG) — produces the per-band `shade[]` DMADLV consumes; simulates tree (x,y) on a
      121×121 grid via DMRANN (Poisson) → mean opacity per MESH band.
    - DMOTHR (114, det) — TRIVIAL per-species tuning scale of the NewSpr/NewInt accumulators: FacS=DMETUN·DMSTUN·F,
      FacI=DMETUN·DMITUN·F, F = 0.05·0.15·0.5·10·0.25·MESH³·TRAJWT = 7.5e-5; DMETUN/DMSTUN/DMITUN default 1.0. Wire
      with the driver (needs the NewSpr/NewInt/ISCT/IND1 driver-local state).
    - DMCYCL (593, det) — life-history compartment advance on the DMINF pools; the last large routine. MAPPED
      (dmcycl.f:390-583): per tree×crown-third, the transition flow is ImmLat=xImm·FProp2 (immature→latent),
      LatAct=xLat·FProp / SprAct=xSpr·FProp (latent/suppressed→active), ActSpr=xAct·BProp (active→suppressed);
      subtract from source, add to destination (xLat+=ImmLat; xSpr+=ActSpr; xAct+=LatAct+SprAct); then survival
      xImm/xLat/xSpr/xAct·=SpSurv; New-infection intake xImm+=New where New=(NewSpr+NewInt)/TVol; then the DMCAP
      capacity saturation New=xPrv+(DMCAP−xPrv)·(1−exp(−(xNow−xPrv)/DMCAP)) with xPrv/xNow = prev/new ACTIVE+SUPRSD.
      ★ ALL biocontrol (BC) terms (BCMORT/BCSUPP/HfLf/DMINF_BC/xImmBC…) are ZERO for base BC/YSM (no MISBCI) ⇒ drop
      them. Coefficients (RESOLVED, dmcycl.f:220-347): DMCYCL runs an ANNUAL loop over LastYr=IFINT(+Spin) years.
      Per crown-third, FProp(i,j)/BProp(i,j) = interpolate the per-species life-history rate curves DMLtRx (forward
      FvecX/FvecY, backward BvecX/BvecY; DMLtnp points) at the crown-third mid-height MESH band IHT; FProp2=1/DMFLWR(sp)
      (1/years-to-flower); SpSurv=DMSURV(sp); DMCAP(sp) default 3.0. ⇒ the port needs the DMLtRx life-history curve
      table + DMSURV/DMFLWR/DMCAP extracted from dminitbc.f (a bounded coeff extraction, like the DG tables) + a
      piecewise-linear interp helper. DMMTRX = no-op orchestrator (DMCW crown_width + DMSHAP + DMSUM/DMRDMX, all
      ported). Plus the event-monitor OPFIND/OPGET DMAUTO mid-cycle scheduling (dmcycl.f:294-329) + xOriginal save
      (:246) + DMNTRD cycle-follow remap. ⇒ NEWSPRED is now MAPPED END-TO-END; every routine understood.
    - Wire the DMTREG driver (the mapped loop) over these, then C6 payoff: mistoe/mismrt (USEMRT extra mortality) +
      BHTG/DHTG growth-mult, wire BC into the DM dispatch. THIS closes the YSM gap (jl SDI→1586 vs oracle 926).
- Multi-session; each chunk lands + validates before the next. Off-switch untouched (USER's).
