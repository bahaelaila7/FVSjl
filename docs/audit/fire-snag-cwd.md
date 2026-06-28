# FFE snags + CWD decay + additions (non-DDW parts) — bandaid audit

Files audited:
- `src/engine/fire/snag.jl`
- `src/engine/fire/fuel_decay.jl`
- `src/engine/fire/fuel_additions.jl`

FVS sources read (SN build dir unless noted): `fmsfall.f`, `fmsnag.f`, `fmsadd.f`, `fmsngdk.f`,
`fmscro.f` (fire/base), `fmcwd.f`, `fmcadd.f` (fire/base), `fmssum.f`, `fmvinit.f`, `fmcrbout.f`,
`fminit.f`, `fmsdit.f`, `fmmain.f`, `FMCOM.F77`.

## Summary

The decay/additions spine is faithful. `fmcwd!` matches `fmcwd.f:78-134` line-for-line (DKR table,
duff-first, PRDUFF=0.02 into the hard duff pool, the `ln(1−DKR)/ln(0.64)` hard→soft transfer).
`fmcadd_litterfall!`, `fmcadd_woody!`, and `_cwd2b_fall!` match `fmcadd.f:72-135` exactly.
The `snag_fall_density` core matches `fmsfall.f:128-178`. The "all ordinary-mortality snags created
HARD" decision is **faithful** (`fmvinit.f:131` sets `PSOFT(I)=0.0` for every SN species, so
`DENIH=DEND, DENIS=0`). `tsoft = (1.24·dbh+13.82)·DECAYX` matches `fmsngdk.f:80`. The dead-root decay
`(1−CRDCAY)^nyrs` is faithful (`fmcrbout.f:273`, with the climate `X=1.0` since `CLCWD` is commented
out; `CRDCAY=0.0425` at `fminit.f:918`), and the input-snag `(1−CRDCAY)^10` matches `fmsadd.f:313-318`.
The cone-taper `_cwd_cone_fractions` reproduces the `fmcwd.f:1000` block (BP breakpoints, RHRAT, BPH,
P1/P2 conic volumes).

Four GAPs flagged below, plus two out-of-scope notes.

---

## GAP 1 — Post-fire accelerated snag fall (PBSOFT/PBSMAL/PBTIME) not implemented

- **jl**: `snag_fall_density` (snag.jl:24) / `update_snags!` (snag.jl:155). The signature takes no
  `IYR`/`BURNYR`; it implements only the "normal conditions" `DFALLN` path.
- **FVS**: `fmsfall.f:98-120` computes `RSOFT`/`RSMAL` for the first `PBTIME` years after a burn, and
  `fmsnag.f:178-215` overrides the normal fall with `PBFRIS·DENIS`/`PBFRIH·DENIH` (taking the larger of
  normal vs post-burn). SN defaults are non-trivial: `fmvinit.f:1101-1104` set `PBSOFT=1.0`,
  `PBSMAL=0.9`, `PBSIZE=12.0`, `PBTIME=7.0` — i.e. after a fire, small (<12") and soft snags fall fast
  (90–100% gone within 7 yr). `BURNYR` is set by `FMBURN`.
- **Impact**: On fire scenarios (e.g. snt stand-4, where `fmburn.jl:94` creates fire-killed snags),
  the jl keeps post-fire snags standing at the slow normal rate. Standing-snag density / Stand-Dead
  carbon will run high for ~7 yr after any burn. Non-fire stands are unaffected.

## GAP 2 — Hard→soft DKTIME transition dropped; FVS_SnagSum hard/soft split is permanently "all hard"

- **jl**: `update_snags!` never moves `den_hard → den_soft`; `snag_decay_fraction` (snag.jl:50) is
  **defined but never called** (confirmed by grep). `add_snag!` seeds `den_soft=0`, and nothing ever
  raises it. `snag_summary` (snag.jl:228) therefore reports every snag as hard for all time.
  The docstring (snag.jl:226) claims the per-record split is "the current-state equivalent of the
  Fortran's DENIH/DENIS + HARD flag" — but the HARD flag is not modeled.
- **FVS**: `fmsnag.f:282-285` flips `HARD(I)=.FALSE.` once `(IYR−YRDEAD) ≥ DKTIME`. `fmssum.f:36-49`
  then books that record's `DENIH` into the **soft** totals (`IF (HARD(II)) THD+=DENIH ELSE TSF+=DENIH`).
  DKTIME for a 12" tree is 2/6/10 yr (snag classes 1/2/3), so snags flip within a cycle or two.
- **Impact**: The fall/carbon paths are *correct* (CWD1 uses the immutable initial `DENIH/DENIS`, and
  the jl's removal of the earlier erroneous `den_hard→den_soft`-at-fall was the right fix — see the
  snag.jl:178-184 comment). But by deleting the flip entirely, the separate **reporting** state was lost:
  the FVS_SnagSum table's hard vs soft columns diverge for any snag older than its DKTIME (all soft
  density shows up under hard). A per-record HARD flag flipping at DKTIME (independent of the density
  split) is the missing piece.

## GAP 3 — `fmscro!` ILIFE uses round() where FVS uses ceil()

- **jl**: `fmscro!` (fuel_additions.jl:76): `ilife = clamp(round(Int, min(tsoft, _fm_tfall(cls,sz))), 1, 60)`.
- **FVS**: `fmscro.f:126-131`: `ILIFE = INT(RLIFE)`; `IF (REAL(ILIFE) .LT. RLIFE .OR. ILIFE .LE. 0)
  ILIFE = ILIFE+1` — i.e. ceil (round **up**) for any non-integer positive value, floor of 1.
- **Impact**: For a fractional part < 0.5 the two disagree (e.g. `RLIFE≈2.3` on a ~15" fast-decay snag,
  DECAYX=0.07: FVS spreads the crown component over 3 yr, jl over 2 yr). Changes `annual = amt/ilife`
  and the CWD2B fall schedule (hence the Stand-Dead crown trajectory). Narrow but real; the jl is
  faithful when the min lands on an integer or a ≥.5 fraction.

## GAP 4 (low) — redcedar foliage TFALL hard-coded to 1 yr

- **jl**: `_fm_tfall` (fuel_additions.jl:54): `sz == 0 && return 1f0` for every species.
- **FVS**: `fmvinit.f:1017-1021`: `TFALL(I,0)=1.0` **except** redcedar (`I==2`) where `TFALL(2,0)=3.0`.
- **Impact**: Redcedar (SN sp index 2) foliage debris falls over 1 yr in the jl vs 3 yr in FVS,
  slightly front-loading its litter (cwd size-10) contribution. Only redcedar; all other species/sizes
  match the `_FM_TFALL1/3/4` tables (verified against `fmvinit.f:1023-1058`).

---

## Out-of-scope notes (DDW path — listed for completeness, not counted as in-scope flags)

### NOTE A (UNVERIFIED) — inventory `snapshot_ffe_oldcrown!` and a self-contradicting comment

`carbon.jl:388` / `summary.jl:101` call `snapshot_ffe_oldcrown!` at inventory, justified by an
output-matching comment (carbon.jl:386-387: "Without this the 1st cycle's crown-lift is skipped …
losing ~1.9 t/ac of fine down-wood"). This is the FFE **crown-lift / DDW** path, which the module
brief scopes out, but the function lives in `fuel_additions.jl` so it is noted here. Two concerns:
(1) the comment at carbon.jl:404 ("Both no-op on the first grow (ffe_old* unset ⇒ zero), matching
FVS's ICYC>1 gate") is **factually contradicted** by line 388, which *sets* `ffe_old*` at inventory —
so the first `compute_crown_lift!` is **not** a no-op. (2) `fmsdit.f:93` (`IF (ICYC.GT.1)`) makes the
first grow cycle's crown-lift exactly zero in FVS. Whether the jl's first *non-zero* crown-lift lands
in the correct FVS cycle slot (FVS cycle-2, faithful) or one cycle too early (bandaid) depends on the
FMOLDC/FMSDIT inter-cycle call ordering in the main driver, which I did not fully trace. Needs that
trace to settle; the prompt cites this as a previously-confirmed bandaid.

### NOTE B (UNVERIFIED) — bole-fall into down wood uses MERCH volume, FVS CWD1 uses FMSVL2 stem volume

`update_snags!` distributes `a = sn.bolevol` (merchantable cubic × V2T, snag.jl:188-200) down the cone.
FVS `CWD1` (`fmcwd.f:152-187`) distributes `TVOLI` from `FMSVL2(...,'D',.false.,.false.)`, which is the
snag **stem** volume, not necessarily the merch cubic. The jl's merch basis was validated for the
Stand-Dead snag *bole carbon* (memory note "FFE snag bole = MERCH"), but the down-wood *fall* addition
may use a different volume basis. This is the DDW addition (out of scope) and resolving FMSVL2's exact
volume requires reading the 200+-line `fmsvol.f`; flagged only as a lead.

### Minor — DZERO snag zeroing not replicated

`fmsnag.f:220-235` removes a snag record entirely once it drops below `DZERO=NZERO/50=0.0002`
(`fmvinit.f:125` NZERO=0.01). `update_snags!` lets densities asymptote toward zero without the cleanup.
Carbon-negligible; affects only snag-record bookkeeping.
