# Management & disturbance completeness tracker (C1-C5)

**Comprehensive** map of every management / disturbance / keyword-option in C1-C5,
built by classifying the **full 140-keyword master table** (`base/keywds.f`) — not a
hand-picked subset. The earlier version of this file missed ~18 keywords; this one is
exhaustive against keywds.f. Disturbance models that are whole extensions (FFE fire,
insects) are C7; econ is C8 — noted here for completeness but owned there.

> **DO NOT START until [`NATURAL_PROCESS_TODO.md`](NATURAL_PROCESS_TODO.md) is done**
> (it is: dynamics complete, .sum bit-exact 89/90, classification all 90). Same
> discipline: port the real FVS code; tests only catch bugs. C6-coupled items live in
> [`C6_DBS_TODO.md`](C6_DBS_TODO.md).

Legend: ✅ done · 🟡 partial · ⛔ unported · ⚪ N/A in SN · 🧊 C7/C8 extension

> **Audit 2026-06-24** (don't trust a ✅ without checking the code is *wired*, not just
> that a unit test passes): **MANAGED** was a **false ✅** — it sat in `KNOWN_NOOP` and
> `plot.managed` was never set, so the DGF planted growth term never fired. Found by
> tracing `plot.managed` from its DGF reader back to its (missing) setter; now ported +
> bit-exact vs Fortran (§5, `test_managed.jl`). Conversely **SPROUT/ESUCKR** was a stale
> **⛔** (actually done & bit-exact). Lesson: a keyword whose effect-variable is *read*
> in growth/mortality but only *set* in `KNOWN_NOOP` is a latent gap — grep the effect
> variable's **setters**, not just the keyword name.

## 1. CUTS — thinning / harvest methods (`cuts.f`, ICFLAG)

| keyword | semantics | status |
|---|---|---|
| THINBTA/THINATA/THINBBA/THINABA | from below/above to residual TPA/BA | ✅ |
| THINDBH | proportional DBH-class to residual TPA/BA | ✅ |
| THINAUTO | auto-thin to FULSTK on stocking trigger (recurring) | ✅ (±1-2 TPA) |
| THINPRSC | prescription thin — remove cut-code-marked (KUTKOD≥2) records at cuteff — snt01 stand 3 | ✅ (`_thinprsc!`; stand 3 bit-exact; nps>1 deferred) |
| THINSDI | thin to target SDI (Zeide summation + proportional CUTEFF) | ✅ |
| THINCC | thin to residual crown cover (CCCLS, forest-grown crown width) | ✅ |
| THINHT | thin a height class (label_325 on height) | ✅ |
| THINRDEN | relative-density thin (Curtis RD) | ✅ |
| THINRDSL | relative-density SDI-line (SILVAH RD) thin | ⚪ N/A in SN (RDCLS2 gated VARACD≠NE → no-op) |
| THINMIST | mistletoe (DMR) thin | ⚪ N/A in SN (no dwarf-mistletoe model → IDMR=0 → no-op) |
| THINPT / SETPTHIN | point (plot-specific) thin (per-point + PI/NONSTK) | ✅ |
| THINQFA | Q-factor diameter-dist thin (CUTQFA + 2-record) | ✅ |

## 2. CUTS — modifiers (`cuts.f`, IACTK 201-206 + setup keywords)

| keyword | effect | status |
|---|---|---|
| SPECPREF | per-species cut preference (RDPSRT order) | ✅ |
| SPLEAVE / LEAVESP | leave named species — exclude from every thin (cuts.f:1466) | ✅ **DONE** (kw_thin! icflag 206 → cuts! PASS 1 `_apply_spleave!` sets Control.leave_species: PRMS(1) <0 group / 0 reset-all / >0 species, PRMS(2)>0 = leave; `_cut_eligible` now takes the StandState and skips any LEAVESP species in both stocking + removal). Bit-exact vs Fortran (test_spleave.jl: THINBBA leaving SM ≠ no-SPLEAVE, matches Fortran) |
| CUTEFF | default cutting efficiency EFF (initre.f:5400) | ✅ **DONE (thinning side)** (kw_cuteff! → Control.cut_eff; THINPRSC/THINAUTO default their blank cuteff field to EFF instead of 1.0 — also fixed a wrong 0.98 stub default → 1.0 per grinit.f:172, and a latent THINPRSC blank-cuteff→0 bug). Bit-exact vs Fortran (test_cuteff.jl: THINAUTO eff=0.5 ≠ no-CUTEFF). NOTE: TOPKILL PRB does NOT default to EFF (verified — only HTGSTOP opt 38 has PRMS(2)=EFF; that narrow sub-case deferred, field mapping unverified) |
| TCONDMLT | tree-condition + special-status cut-priority weights (TCWT·IMC + SPCLWT·ISPECL, cuts.f:1074/1424) | ✅ **DONE (both terms)** (kw_thin! icflag 202 → cuts! PASS 1 sets Control.total_wt=PRM1 + special_wt=PRM2; RDPSRT cut key is the full cuts.f:1074 formula via `_cut_pref_wt`, IMC=clamp(mort_code,1,3); basdam code-55 → t.special/ISPECL). Bit-exact vs Fortran on BOTH weights (test_tcondmlt.jl: IMC=3 + TCWT, and code-55 ISPECL + SPCLWT, each cut first ≠ no-TCONDMLT). Default 0 ⇒ key unchanged |
| YARDLOSS | yarding-loss → scales removed merch/saw/bdft by (1−prlost) **and feeds the FFE down-wood/snag/crown fuel pools** | 🧊 **rolled into C7** (substantive effect is fuel-pool routing; standalone .sum effect nil; its `@test_broken` is the post-thin DGSCOR tail, not yardloss) |
| MINHARV | minimum-harvest thresholds — cancel a cut whose total removal is below ANY minimum (cuts.f:400/1556) | ✅ **DONE** (kw_thin! icflag 200 → cuts! PASS 1 sets ba_min/tcf_min/cf_min/scf_min/bf_min; after the methods, if BA/total/merch/sawlog-cubic/board-feet removed < the thresholds, restore the pre-thin TPA snapshot + return no removal). Bit-exact vs Fortran (test_minharv.jl: THINBBA + BAMIN=100 cancels the cut → identical to no-thin). Default 0 ⇒ no-op |
| TFIXAREA | total fixed plot area → small-tree expansion FP=1/TFPA (initre.f:816 / notre.f:45) | ✅ **DONE** (kw_tfixarea! sets s.plot.total_fixed_plot, which notre! already consumes — the field was wired into the expansion but never set). Bit-exact vs Fortran (test_tfixarea.jl: TFPA=0.02 changes TPA 536→786) |
| SPGROUP | species groups (vbase/initre.f:4726) referenced by a −N species field | ✅ (kw_spgroup! builds the table; the ISPCC<0 branch is wired into FIXDG/FIXHTG/FIXMORT/HTGSTOP/TOPKILL/MORTMULT/CRNMULT/TREESZCP/SPECPREF via sp_field_matches AND the species-filtering thins THINDBH/SDI/RDEN/CC/PT/QFA via _cut_eligible; bit-exact, test_spgroup.jl) |

## 3. Growth keyword multipliers / overrides (`dgdriv.f`/`htgf.f`/`regent.f`)

| keyword | effect | status |
|---|---|---|
| FIXDG | fix/scale diameter growth | ✅ (one-shot scaler, species×DBH window, scales tripled DG; bit-exact, test_fix_scalers.jl) |
| FIXHTG | fix/scale height growth | ✅ (one-shot scaler, species×DBH window, scales tripled HTG; bit-exact, test_fix_scalers.jl) |
| HTGSTOP / TOPKILL | scale height growth / top-kill (htgstp.f) | ✅ (act 110 HTG×PKIL + act 111 top-kill w/ NORMHT/ITRUNC/Behre/crown; deterministic + stochastic (htgstop_stoch) bit-exact through firing cycle, test_htgstp.jl) |
| BAIMULT | basal-area-increment multiplier (scales DDS) | ✅ (MULTS; bit-exact vs Fortran, test_multipliers.jl) |
| HTGMULT | height-growth multiplier | ✅ (MULTS; bit-exact vs Fortran) |
| CRNMULT | crown-ratio-change multiplier (sn/crown.f:319) | ✅ (active_crn_mult; scales the limited CR change over a DBH window, persistent; bit-exact within ±1 drift, test_crnmult.jl) |
| FIXCW | fix crown width | ✅ recognized no-op — **verified .sum-inert** (CRWDTH is output-only, never fed into CCF/growth; a live-Fortran FIXCW run is byte-identical to no FIXCW). In KNOWN_NOOP |
| REGDMULT / REGHMULT | regen diameter / height growth multiplier | ✅ (MULTS kinds 6/3; regent XRDGRO/XRHGRO; ±1 vs Fortran on regen cycles) |
| NOTRIPLE / NUMTRIP | tripling control (ICL4) | ✅ (NOTRIPLE→icl4=0, NUMTRIP n→icl4=n; bit-exact vs Fortran, test_tripling.jl) |

## 4. Mortality keyword overrides (`morts.f`/`fixmort.f`)

| keyword | effect | status |
|---|---|---|
| FIXMORT | keyword mortality rate override | ✅ normal path (replace/add/max/mult, DBH window, one-shot; bit-exact, test_fixmort.jl) **+ SIZE concentration** PRM(6)=10/20 (KBIG bottom-up/top-down): pools the window's kill into XMORE then re-imposes it whole-record on trees ranked by ∓(DBH+DG/bark) via the faithful `_rdpsrt!`, until XMORE is spent (morts.f:838). fixmort_big.key (pflag=10) bit-exact vs live Fortran on TPA/BA/SDI/QMD across 11 cycles. **+ KPOINT** point concentration PRM(6)=1 (point-by-point, morts.f:937) and **combined** PRM(6)=11/21 (size-within-point, morts.f:978) — DONE & bit-exact on the 11-point base stand (fixmort_kpoint/fixmort_kpbig, using t.plot_id=ITRE + points_inv=IPTINV). FIXMORT concentration now COMPLETE |
| MORTMSB / MATUREW | MSB mature-stand break-up mortality (msbmrt.f) | ⚪ EFFECTIVELY INERT for self-thinning/managed stands (verified): fires only when survivors EXCEED the 85% mature self-thinning line, but BAMAX/self-thinning hold the stand AT that line, so TMORE=0 — even a 30-cycle run to QMD 38 doesn't trigger it. Rare-trigger (overmature break-up only); deterministic if ported |
| MORTMULT | mortality-rate multiplier (background only + DBH window, morts.f:518/524) | ✅ (MULTS; DBH window D1≤DBH<D2 via active_mort_mult; bit-exact on bg-mortality cycles, windowed + windowless) |
| TREESZCP | per-species size cap (SIZCAP): DG bound + size-cap mortality + HT cap | ✅ (keyword + morts size-cap floor + htgf HT cap; nomort path bit-exact, see §SIZCAP) |
| SDIMAX | per-species max SDI (SDIDEF) + self-thinning percents PMSDIL/PMSDIU (initre.f:3072, option 89) | ✅ **DONE** (was another UNTRACKED gap, found 2026-06-24 — keyword unrecognized → ignored). `kw_sdimax!`: field 1 species (0/blank=all, −N=SPGROUP, code), field 2>0 ⇒ `sp_sdi_def[sp]` (preserved by `site_index_setup!`'s ≤0 guard, the MAXSDI flag), fields 5/6 ⇒ `pct_sdimax_mort_lo/hi` as fractions (PMSDIL≥10 / PMSDIU≤95, ÷100). Confirms FVSjl's fraction convention (the keyword stores percents). Bit-exact vs live Fortran (`test_sdimax.jl`, `sdimax.key`: all-species max SDI 300 caps the stand; ±Scribner board-foot noise) |
| BAMAX | user-pinned maximum basal area → SDImax self-thinning cap (initre.f:6800, option 66) | ✅ **DONE** (was a fully-UNTRACKED gap, found 2026-06-24 by the read-but-never-set audit: `control.ba_max` was consumed in `site_index.jl` but the BAMAX keyword was *unrecognized* → silently ignored). `kw_bamax!` sets `ba_max`; `site_index_setup!` derives `sp_sdi_def = BAMAX/(0.5454154·PMSDIU)` (sdical.f:208) so the SDImax mortality caps residual BA at BAMAX. **Also fixed two latent bugs in the never-before-exercised BAMAX→SDIDEF branch**: no PMSDIU 0.85 default (⇒ ÷0 ⇒ `sp_sdi_def=Inf` ⇒ no cap) and a stray `/100`. Bit-exact vs live Fortran (`test_bamax.jl`, `bamax.key`: a dense stand caps BA at 136 like Fortran; ±1 board-foot Scribner noise). The LMORT mortality-enable flag is implicit (FVSjl always runs mortality) |

## 5. Other stand management

| keyword | effect | status |
|---|---|---|
| PRUNE | pruning (option 108, act 249) | 🧊 .sum-INERT in SN (verified: no TPA/BA/TopHt/QMD change) — feeds ecopls.f pruned-log volume (C8 econ), the crown edit doesn't reach SN growth |
| FERTILIZ / FFERT | fertilizer growth response (ffin.f/ffert.f) | ✅ **DONE** (kw_fertiliz! → Control.fertilize_events; `fertilizer_growth!` in grow_cycle! after TRIPLE, grincr.f:564): 200-lb-N response boosts each tree's DDS by RDDS=exp(0.1108·lnD+0.003004·BAL/ln(D+1)) (cap 2.6) and HTG by 1.1626, for the IFLEN of the cycle's years within 10 yr of application, scaled by efficacy; carries over via ifert_date/ifert_eff. Bit-exact vs Fortran (test_fertiliz.jl: BA 126→138 at 2000). SN is outside the calibrated DF/GF range — Fortran warns but applies it, so we match |
| COMPRESS | record compression to a target (comprs.f act=250) | 🟡 **keyword DONE, algorithm tracked** (`kw_compress!` recognizes + schedules act 250 with target/PN1/date — no longer silently dropped; `cuts!` skips icflag 250 so records pass through uncompressed). The 1010-line IBM-SSP-eigen PCA clustering is scoped in [`COMPRESS_chunk_plan.md`](COMPRESS_chunk_plan.md) (helpers → EIGEN → classify → merge); a COMPRESS stand still diverges by the compression until that lands. `test_compress.jl` covers keyword recognition/scheduling only |
| ADDFILE | (option 22) redirect input to another unit | ⚪ NOT a tree-add — it just switches IREAD to a Fortran file UNIT (an include-file mechanism); no clean FVSjl mapping (FVSjl uses filenames, not units). 'ADDTREES' is not an SN keyword |
| MANAGED | managed-stand flag (DGF kplant term) | ✅ **DONE** (was a FALSE ✅ — found 2026-06-24: the keyword was in `KNOWN_NOOP` and `plot.managed` was never set, so the DGF planted term `dg_planted[sp]·kplant` (dgf.f:179/328) never fired). Now `kw_managed!` sets `plot.managed` per initre.f:10000 (immediate path: bare card / field-2≠0 ⇒ managed=1, field-2==0 ⇒ 0; dated OPNEW-act-82 path deferred). **Bit-exact vs live Fortran** on a loblolly-pine stand for BOTH managed and unmanaged (`test_managed.jl`, `managed.key`; the planted-pine boost self-thins the stand slightly more, matching Fortran). Non-planted species have `dg_planted=0` so MANAGED is correctly inert for them (which is why snt01 never exposed the gap) |
| MGMTID / RESETAGE / SETSITE | mgmt id / reset age / set site mid-run | 🟡 MGMTID read; ✅ RESETAGE; ✅ SETSITE (scheduled per-species site index + BAMAX/SDImax + dgcons! recompute, bit-exact vs Fortran; habitat param not wired) |

## 6. Volume / defect keywords (C5 — **.sum-affecting**, keyword-settable)

> NOTE: `DEFECT/BFDEFECT/MCDEFECT` set defect % from the KEY — so G1's defect IS
> reachable from a `.key` (not only DBS). These belong here, not just C6.

| keyword | effect | status |
|---|---|---|
| MCDEFECT | per-species CUBIC defect curve (CFDEFT) → reduces merch cubic | 🟡 **DONE** (kw_mcdefect! → Control.sp_cf_defect 9×MAXSP; FVSsn vols.f:294-332 SN branch: pulpwood part MCFV−SCFV reduced by ICDF%=NINT(ALGSLP(DBH,CFDEFT)·100) clamp[0,99], sawtimber untouched; ALGSLP segmented-linear over DBHCLS=[0,5..40]). Undated→immediate (affects cyc0). **Bit-exact vs live Fortran** (test_mcdefect.jl, fires −144..530 cuft/ac). Deferred: dated scheduling, per-tree DEFECT input, CFLA0/CFLA1 log-linear form model (default no-op) |
| BFDEFECT | per-species BOARD-FOOT defect curve (BFDEFT) → reduces board feet AND sawtimber cubic | 🟡 **DONE** (kw_bfdefect! → Control.sp_bf_defect; vols.f:419-432: BFV·(1−IBDF/100) AND SCFV·(1−IBDF/100), ≥99⇒both 0; then MCFV=PULPV+post-board SCFV couples it into merch cubic). Bit-exact vs Fortran incl. the coupled MCDEFECT+BFDEFECT case (test_mcdefect.jl) |
| DEFECT (per-tree) | per-tree CF/BF defect from the tree damage codes (basdam.f: agent 25=both, 26=cubic, 27=board; severity=%) | 🟡 **DONE** (treeinput.jl packs t.defect = CF·1e6+BF·1e4 from the damage pairs; compute_volumes! folds ICDF=DEFECT/1e6 / IBDF=DEFECT/1e4 mod 100 into the max with the CFDEFT/BFDEFT curves). Bit-exact vs Fortran (test_pertree_defect.jl, agent 26 sev 30 → merch cubic 1149→894). CFLA0/CFLA1/BFLA0/BFLA1 log-linear form model now ported via MCFDLN/BFFDLN (below; no Fortran oracle — FPE) |
| MCFDLN / BFFDLN | cubic/board log-linear form-model coefs CFLA0/CFLA1, BFLA0/BFLA1 (sdefln.f opt 39/40) | 🟡 **ported** (kw_mcfdln!/kw_bffdln! → Control.sp_cf_form0/1 + sp_bf_form0/1; compute_volumes! folds VOLCOR=exp(B0+B1·ln(V)) reduction into ICDF/IBDF, vols.f:303-310). Default 0/1 = no-op, bit-exact. ⚠ **NO Fortran oracle**: the SN build FPE-crashes when the form model is active (vols.f:306 ALOG before the TEMVOL==0 guard → log(0) on any MCFV==SCFV tree). FVSjl guards temvol>0 + is deterministic; test_mcfdln.jl pins no-op default + fires. NOT bit-exact-validated (DIVERGENCES.md) |
| VOLUME / BFVOLUME | per-species cubic / board-foot merch-standard overrides (volkey.f:9915/9905) | ✅ **DONE** (kw_volume!/kw_bfvolume! → per-stand Control.sp_* arrays via init_merch_standards! + scheduled apply_volume_overrides!; 0/+species/−SPGROUP). VOLUME DBHMIN gate bit-exact (test_volume_override.jl). **BFVOLUME bit-exact** incl. board feet via the BFPFLG=0 separate board call + Region-8 ≥10ft-product zeroing (test_bfvolume.jl, BFTOPD 9→11). FRMCLS/METHC/METHB ignored (no form-class/method selector in R8 Clark). LFIANVB=.FALSE. for SN |
| VOLEQNUM | per-species cubic volume-equation override (VEQNNC, initre.f:5061) | ✅ **DONE incl. board feet** (kw_voleqnum! → Control.voleqnum_overrides → apply_voleqnum_overrides! after VOLEQDEF; species via resolve_species alpha/FIA or −N group). Board feet kept on the default equation via the per-stand `sp_bf_vol_eq` snapshot + the BFPFLG=0 separate board-foot call (fvsvol.f:362). **Fully bit-exact vs Fortran** (test_voleqnum.jl, SM→AB eq: total cuft 1368→1419, board feet matches) |
| CFVOLEQU / BFVOLEQU | per-species cubic/board volume-equation (old keywords) | ⚪ DEPRECATED ("NO LONGER ACTIVE", errgro.f:514/524) — superseded by VOLEQNUM |
| FIAVBC | FIA volume/biomass calc switch | ⚪ OUT OF SCOPE — switches to the FIA National Volume Library; FVSjl has only the R8 Clark equations (the SN default LFIANVB=.FALSE.), so the FIA NVB path is unsupported. Recognized no-op (KNOWN_NOOP) so it doesn't error; a stand that sets FIAVBC would diverge (FVSjl stays on R8 Clark) |

## 7. Event monitor & activity scheduling (`evmon.f`/`opcycl.f`)

| keyword | effect | status |
|---|---|---|
| IF / THEN / ENDIF | conditional activity scheduling (snt01 stand 2) | ✅ event_monitor.jl (AST evaluator); stand 2 first 2 thins bit-exact; 3rd = class-boundary residual |
| COMPUTE | event-monitor variable assignment | ✅ (kw_compute! parses NAME=expr…END, evaluated each cycle in cuts! before IF conditions read them via compute_vars; ≡ direct ref, bit-exact lead stand, test_compute.jl) |
| TIMEINT | cycle length (period scaling) | ✅ uniform **and per-cycle** (field-1 = cycle N → `cycle_lengths[N+1]`, cumulated by `build_cycle_schedule!`; DDS·FINT/5 + HTG·FINT/5 + mortality^FINT + year/age from the IY schedule; snt01 bit-exact, TIMEINT-10 TPA ≤8 / volume ≤2% — calibrated-species sp33/65 non-5-yr period residual; test_timeint.jl). GROWTH (input-growth-data calibration keyword, ≠ per-cycle mult) still deferred |
| CYCLEAT | extra cycle-boundary year | ✅ **DONE** (opt 134, initre.f:13400 + fvs.f:116-135). `kw_cycleat!` collects de-duplicated calendar years; `build_cycle_schedule!` inserts each as a new boundary strictly inside the run (no end-extend / start-move), bumping `ncycle_eff`. Built the **non-uniform cycle schedule** (the FVS IY array) it needs: `cycle_year`/`cycle_period_at`/`current_cycle_year`, routing every per-cycle year derivation (summary, growth, mortality, cuts, fire, fert, econ, multiplier windows) through it — **bit-exact for uniform cycles** (snt01 unchanged). vs live Fortran: YEAR + PrdLen columns **bit-exact**; stand within the same TIMEINT non-5-yr residual (confirmed identical in a plain uniform-TIMEINT-3 control → period scaling, not CYCLEAT). test_cycleat.jl |
| ESTAB-block (TALLY/PLANT/NATURAL/SPROUT) | establishment scheduling | ✅ PLANT/NATURAL **and SPROUT** (ESUCKR stump-sprout regen — bit-exact vs live Fortran on `sprout.key`, 113 tests); ⛔ TALLY counts only |

## 8. Disturbance models (C7/C8 extensions — owned there)

> ★ **INSECTS/DISEASE ARE NOT PART OF THE SN VARIANT (2026-06-24).** The rebuilt SN binary links
> only the `ex*.f` no-op STUBS for every insect/disease model (exmpb/exdfb/exdftm/exmist/exbm/exrd);
> there is NO real model (mpbmza/dfbmza/… are absent from `bin/FVSsn_buildDir`). A `DFB`/`MPB`/… keyword
> in SN emits `FVS11 ERROR: REQUESTED EXTENSION IS NOT PART OF THIS PROGRAM` and leaves the simulation
> unchanged (verified: DFB on snt01 ⇒ byte-identical `.sum`). So these are **NOT unported omissions** —
> they are not in SN FVS at all (mostly western insects/diseases with no Southern hosts). The faithful
> SN behaviour is recognize-and-do-nothing, which FVSjl now does explicitly via
> `SN_UNSUPPORTED_EXTENSIONS` (southern.jl) ⇒ they are **out of scope for the SN port**, not a gap.

| keyword | model | status |
|---|---|---|
| FMIN … END | FFE fire (SIMFIRE/SALVAGE/fuels/snags/CWD/carbon) | 🧊 C7 (FFE IS in SN) |
| MPB / DFB / DFTM / WSBW / MISTOE / BRUST | mtn pine beetle / DF beetle / DF tussock moth / W spruce budworm / mistletoe / blister rust | ⚪ **NOT in SN** (ex* stub ⇒ "NOT PART OF THIS PROGRAM"); recognized-inert in `SN_UNSUPPORTED_EXTENSIONS` |
| RDIN / ANIN / RRIN | root-disease model (Western root disease / Annosus) input | ⚪ **NOT in SN** (exrd stub); recognized-inert |
| PRMFROST / CLIMATE | permafrost / climate-FVS modifiers | 🧊 |
| ECON / CHEAPO | economic analysis | 🧊 C8 (ANNUCST path exists) |

## Keyword-recognition audit — full keywds.f cross-check (2026-06-24)

Diffed the **142-keyword SN master list** (`keywds.f`) against everything FVSjl recognizes
(dispatch + `KNOWN_NOOP` + block sub-keywords). FVSjl recognizes 103 tokens; the rest are
**silently ignored** (fall through to the no-op `else`). This is exactly how MANAGED/BAMAX/
SDIMAX hid — so the unrecognized list is triaged here so each is a *visible* decision.

**⚠ Simulation-affecting, UNRECOGNIZED (silent `.sum` gaps, same class as MANAGED/BAMAX/SDIMAX):**
- `RANNSEED` (opt 61) — ✅ **DONE** (audit find #4). `kw_rannseed!` reseeds the main RNG via
  the existing `ranseed!` (forced odd; field-1==0 GETSED clock seed skipped as non-deterministic).
  Bit-exact vs live Fortran (`test_rannseed.jl`, seed 12345; ±DGSCOR/Scribner noise).
- `DGSTDEV` (opt 57) — ✅ **DONE** (audit find #5). `kw_dgstdev!` sets the new per-stand
  `control.dg_stddev_bound` (default 2.0); threaded into `dgscor!` + the OLDRN clamp in
  `diameter_growth!`/`small_tree_growth!` (was a hardcoded const). Bit-exact vs live Fortran
  (`test_dgstdev.jl`, DGSTDEV 0 ⇒ deterministic DG matches exactly; default 2.0 keeps snt01
  bit-exact).
- `SERLCORR` (opt 91) — ✅ **DONE** (audit find #7). `kw_serlcorr!` sets the new per-stand
  `control.dg_bjphi`/`dg_bjthet` (ARMA(1,1) AR/MA, defaults 0.74/0.42); `_stand_bjrho` recomputes
  the BJRHO autocorrelation series only when overridden (default path keeps the precomputed const
  ⇒ zero overhead + snt01 bit-exact), threaded into both `autcor` calls. Bit-exact vs live Fortran
  (`test_serlcorr.jl`, phi 0.50 / theta 0.30).
- `READCORD/READCORH/READCORR` + `REUSCORD/REUSCORH/REUSCORR` — ✅ **DONE** (NOT cross-run file
  persistence — earlier mischaracterization). Each READ* reads a block of MAXSP per-species
  correction terms (8F10.0 continuation lines, blank ⇒ 0) that modify a growth-model CONSTANT
  before the LSTART calibration: COR2 → `DGCON += ln(COR2)` (dgf.f:1168, large-tree DG); HCOR2 →
  `HTCON += ln(HCOR2)` (htgf.f:332, large-tree HT); RCOR2 → `RHCON = RCOR2` (regent.f:585, small-tree
  HT multiplier). REUSE* re-enables the prior terms without re-reading. `kw_readcor{d,h,r}!` +
  `read_species_corr!` + `kw_reuscor{d,h,r}!`; applied in `dgcons!` (COR2/HCOR2) and the REGENT con
  (RCOR2). Control gained `dg_cor2/htg_cor2/regh_cor2` (default 1) + `*_on` flags. Default off ⇒
  no-op (snt01 bit-exact). READCORD COR2=1.3 vs live Fortran: every structural column bit-exact,
  board-feet within ±Scribner Float32 noise; confirmed it changes the stand. `test_readcor.jl`.
- `NOCALIB` (opt 56) — ✅ **DONE** (audit find #6). Disables DG self-calibration per species
  (0/all, −N group, code). `control.dg_calib_sp` (LDGCAL) was declared-but-DEAD (defaulted
  all-`false`, never read) — flipped to all-`true` and now gates the COR fit in
  `calibrate_diameter_growth!` (a skipped species keeps `dg_cor`/`dg_cor_goal`=0). `kw_nocalib!`
  with SPDECD species decode. The SN LHTCAL side is naturally inert (FVSjl does no large-tree
  HT self-calibration — `htg_cor`=0 unless HCOR2). Bit-exact vs live Fortran (`test_nocalib.jl`,
  NOCALIB 0 ⇒ uncalibrated DG matches exactly; snt01 stays bit-exact with the all-true default).
- `CCADJ` (opt 444) — crown-competition-factor adjustment (sets `CCCOEF`, sstage.f:923
  `UPDATECCCOEF`). ✅ **Verified `.sum`-inert in SN** (like PRUNE/FIXCW): `CCCOEF` is read ONLY
  by `covolp.f` (the COVER canopy-cover report), `sstage.f` (the SSTAGE structure-stage code),
  and `evldx.f:430` (the cover event-monitor var) — never by the core growth/mortality/density.
  Empirically: live-Fortran snt01 `+CCADJ 0.5` is byte-identical to baseline (only the run
  timestamp differs). In KNOWN_NOOP; revisit only when COVER/SSTAGE output is ported (C6/output).
- `GROWTH` (opt 13, vbase/initre.f:2300) — 🟡 **RECOGNIZED + params captured** (`kw_growth!`,
  test_growth.jl). It sets how the INPUT tree records' growth fields are interpreted for the LSTART
  calibration: `IDG`/`IHTG` data-type codes (0 = none/increment [grinit default], 1/3 = the field is
  past DBH/HT → `PDBH`/`PHT`, 2 = increment) + measurement periods `FINT`/`FINTH`/`FINTM`. FVSjl now
  stores all 5 (`control.growth_*`); the **default (IDG/IHTG=0, periods=5) is the current bit-exact
  behaviour** (the DG field is the 5-yr increment) — verified a bare GROWTH is a no-op in live Fortran
  too (snt01 + GROWTH ⇒ byte-identical `.sum`). ⛔ **Deferred (the behavioral part):** the IDG=1/3
  past-DBH interpretation (intree.f:531-537 / cratet.f:170-185 / dgdriv.f:330) + non-default FINT
  scaling of the calibration — **WK3 past-DBH calibration territory** (sp33/65). These CHANGE the
  calibration, so they're not wired blind; they need a purpose-built past-DBH / non-5-yr-FINT scenario
  to validate bit-exact (the WK3 residual area). So GROWTH is no longer silently dropped, and the
  default/no-GROWTH case is faithful; only the non-default IDG/FINT behavior remains.
- `CYCLEAT` — ✅ **DONE** (see the keyword table above): explicit cycle-boundary years built on the
  new non-uniform IY schedule; bit-exact YEAR/PrdLen vs Fortran, stand within the TIMEINT residual.
  `SDICALC` — SDI
  method — ✅ **DONE** (audit find #9). The entanglement resolved cleanly: `zeide_sdi` (LZEIDE)
  is `true` after init (SN default; my first check read the bare *constructor*), and `mortality.jl`
  ALREADY consumes it — so `kw_sdicalc!` setting the SHARED `zeide_sdi` (+ `dbh_zeide`/`dbh_stage`)
  routes BOTH the reported SDI column and the SDImax mortality through the method. Made `stand_sdi`
  honor it (Zeide Σ vs the Reineke `SPROB·A+B·SDSQ` Taylor form, sdical.f:281-327) + the right
  threshold. **Bit-exact MULTI-CYCLE vs live Fortran** (`test_sdicalc.jl`, SDICALC→Reineke:
  TPA/BA/SDI/QMD match every cycle, ±Scribner board-feet); snt01 (default Zeide) unaffected. The
  earlier report-only attempt that diverged at cycle 1 was the tell that the mortality must follow
  — fixed by using the shared flag instead of a separate report flag.

  ~~method — **scoped, entangled with FVSjl's multi-path SDI:**~~ the fields exist + are read
  (`control.zeide_sdi`/`dbh_zeide`/`dbh_stage`/`sdi_method`, consumed by `mortality.jl`), so
  `kw_sdicalc!` setting them *looks* like a BAMAX-class wire-up — BUT `zeide_sdi` defaults
  `false` (Reineke) while `stand_sdi` (the `.sum` SDI column) hardcodes Zeide per-tree
  summation with the `dbh_sdi` threshold, and both are bit-exact for snt01. **Reconciliation
  DONE (the entanglement is understood, 2026-06-24):** SN default is `LZEIDE=.TRUE.`
  (grinit.f:129, Zeide). FVSjl is correct via TWO independent paths — `stand_sdi` does Zeide
  per-tree summation for the reported `.sum` column (matches LZEIDE=true), while `mortality.jl`
  uses the QMD path `sqrt(ΣD²/ΣTPA)` (the `zeide_sdi=false` branch), which is the SN SDImax
  computation independent of the report flag. So `zeide_sdi=false` is the *mortality* method,
  NOT a mis-set report flag. **Remaining port:** make `stand_sdi` honor the SDICALC report
  method (Zeide vs the Reineke `SDIC=SPROB*A+B*SDSQ` Taylor form, sdical.f:281-327) + the
  `dbh_zeide`/`dbh_stage` threshold (not the always-0 `dbh_sdi`), defaulting to today's Zeide/0
  behavior; `kw_sdicalc!` sets a *report* method flag + the thresholds. Validate the reported
  SDI column for both methods + a non-zero threshold vs Fortran. ⚠ **EMPIRICAL FINDING
  (2026-06-24, tried + reverted a report-only port):** SDICALC's method drives BOTH the
  reported SDI column AND the SDImax **mortality** — a `SDICALC`→Reineke stand matches Fortran's
  SDI column at cycle 0 but its **TPA diverges from cycle 1** (FT 486 vs report-only 507),
  because Fortran's SDImax self-thinning uses the SDICALC SDI method too. So the full port must
  ALSO route the mortality SDImax through the method (and resolve how the FVSjl `zeide_sdi`
  mortality flag — `false` for the SN Zeide default — maps to LZEIDE). Report-only is NOT enough;
  do both together + validate the multi-cycle TPA, not just the cycle-0 SDI column.
  `MGMTID` is read; `RESETAGE` ✅; `SETSITE` ✅ (below); `GROWTH` is **OPNEW-scheduled** (act
  443 / 120 / 444 / —) needing a per-cycle scheduled-activity handler (FVSjl has no non-cut
  activity dispatch in `grow_cycle!` yet — build it once, then plug each in):
  - `RESETAGE` (resage.f, act 443): ✅ **DONE** (audit find #8). Turned out NOT to need the
    scheduler — it is a pure function of the row year: `kw_resetage!` stores `age_reset_year`/
    `_age` (resolving a cycle-number date against INVYEAR), and `summary_row` rebases `age(Y) =
    age_reset_age + (Y − reset_year)` for `Y > reset_year` (the reset row keeps the old age,
    matching RESAGE running after DISPLY). Bit-exact vs live Fortran (`test_resetage.jl`,
    RESETAGE 2000 30 → AGE rebases 70→35→40…, MAI recomputed; snt01 unaffected).
  - `SETSITE` (act 120): ✅ **DONE** — scheduled mid-run site change. `kw_setsite!` schedules an
    act-120 `ScheduledActivity` (date + 6 params); `apply_setsite!` (called at the top of
    `grow_cycle!`, the RCON position) resets per-species `sp_site_index` (SITEAR, direct or % change,
    clamped ≥ 1) + optional BAMAX (→ `ba_max` + bamax-derived `sp_sdi_def`) + SDImax, then recomputes
    the DG constants via `dgcons!` (= RCON); the new site index also feeds HTG/small-tree growth
    directly. `species_selector` decodes field 4 (0/all, −group, index). vs live Fortran (SI 60→80 @
    2000): every structural column bit-exact, volume within ±1 cuft Float32 noise; confirmed it
    changes the projection. `test_setsite.jl`. ⚠ Habitat (param 1) not wired — SN growth keys off
    forest type, not the habitat code (documented gap; a non-zero habitat is ignored).
  (`BMIN` is NOT a simple gap — it is the WWPB insect-model input (exbm.f), an *extension*,
  belongs with insects below.)

> ⚠ **Audit status (2026-06-24):** the clean *immediate value-setter* wins are DONE
> (MANAGED/BAMAX/SDIMAX/RANNSEED/DGSTDEV — all bit-exact). Every keyword still on this list
> is now **focused-session work**, not a quick wire-up: the NOCALIB/SERLCORR/`*CORR`
> calibration cluster touches the precision-sensitive calibration core; CCADJ/GROWTH/
> CYCLEAT/SETSITE/RESETAGE are OPNEW-scheduled activities needing scheduler+handler
> integration (and none are exercised by a current test scenario); COMPRESS's algorithm +
> insects + C6 DBS are large subsystems. Port these with focused attention per chunk.

**Already triaged elsewhere in this doc** (no new action): `PRUNE` (.sum-inert), `YARDLOSS`
(C7 fuels), `MORTMSB` (effectively inert), `BFVOLEQU/CFVOLEQU` (deprecated), `THINMIST/
THINRDSL` (⚪ N/A in SN), `ADDFILE` (unit redirect), `DEFECT` (per-tree path handled in
treeinput.jl — but the keyword default-setter is unverified).

**C7/C8 extensions (owned there):** insects `MPB/DFB/DFTM/WSBW/BRUST/MISTOE`, root disease
`RDIN/ANIN/RRIN`, `CHEAPO` (econ), `CLIMATE/PRMFROST`, `ORGANON` (alt growth model).

**Output / report / control (C6 or output-only — no `.sum` math):** `DATABASE/DATASCRN/
DELOTAB` (DBS output), `SVS` (visualization), `COVER` (canopy report), the label/metadata
keywords `AGPLABEL/SPLABEL/SPCODES/STANDCN/MODTYPE/STRCLASS/LOCATE`, `OPEN/CLOSE` (file I/O),
`POINTREF/PTGROUP`, `ALSOTRY`, `CWEQN`.

> Next audit-driven ports (most-upstream first): `RANNSEED` + the `SERLCORR`/`*CORR`/`DGSTDEV`/
> `NOCALIB` calibration-&-correlation cluster (they govern the stochastic DGSCOR path that is
> already the dominant residual everywhere), then `CCADJ`, then `GROWTH`/`CYCLEAT`.

## Validation status — 3-way sweep vs live Fortran (2026-06-22)

The comprehensive 3-way sweep (162 scenarios × with/without management, vs live
Fortran) confirms the **cut logic of every ported method is correct**; the only
management-scenario residuals are *post-cut* tails, not thinning bugs:

| scenario | thin | finding |
|---|---|---|
| `s11`/`s28`/`s29_thinbta` | THINBTA | cut **bit-exact** (536→162 at the thin); post-thin ±4 TPA drift over 7 cycles = the **post-thin DGSCOR/serial-correlation tail** (cut re-ranks the stand → the stochastic DG/mortality responds slightly differently; the increment even flips sign cycle-to-cycle, so it does not propagate). Not a cut bug. |
| `s28_thindbh` | THINDBH | bit-identical (snt01 block 2). |
| `cut_thinprsc` | THINPRSC | Δ2 at the cut because the scenario uses `DESIGN …11.0` = **nps=11 plots** — the deferred multi-plot THINPRSC path (single-plot/snt01 stand 3 is bit-exact). + post-thin tail. |
| `cut_yardloss` | YARDLOSS | removes 0 merch → .sum-neutral; the ±9 TPA is the same post-THINDBH accretion tail, not yardloss. C7-coupled for the fuel pools. |
| `cut_thinsdi` | THINSDI | ✅ **ported** — bit-exact in TPA/BA every cycle (Zeide summation SDI + proportional CUTEFF); ±1 cuft tail only. |
| SPECPREF / IF-THEN | — | cut-preference + event-monitor blocks fire and cut correctly. |

**Conclusion:** ported thinning methods (THINBTA/ATA/BBA/ABA, THINDBH, single-plot
THINPRSC, SPECPREF, IF/THEN event monitor) are cut-exact. Remaining work is the
*unported* methods below + the multi-plot THINPRSC path; the post-thin numeric tail is
the same single-precision/serial-correlation floor seen in the natural-process runs.

## Triage: which ⛔ items are actually APPLIED in SN code (2026-06-23)

Grepped each keyword's effect-variable for READ references in `sn/`+`base/` (beyond
init/keyword-table). This separates real ports from set-but-not-read no-ops:

**Genuinely applied (real ports, each non-trivial):**
- `SIZCAP`/TREESZCP — ✅ DONE. The **TREESZCP** keyword (kw_treeszcp!, keyword_dispatch.jl)
  loads SIZCAP[is,1..4] immediately (no date): field 1 = species (0=all), 2 = cap DBH,
  3 = mortality rate, 4 = IDMFLG flag, 5 = HT cap (field order confirmed empirically vs
  live Fortran). The three effects: (a) DG bound (dg_bound, already present); (b) size-cap
  MORTALITY floor — ported in mortality.jl AFTER _varmrt!, before BAMAX, matching
  sn/morts.f:692: if (D+G)≥SIZCAP[is,1] & IFIX(SIZCAP[is,3])≠1 ⇒ killed=max(killed,
  P·SIZCAP[is,2]·FINT/5)≤P, where **G is OUTSIDE-bark, period-scaled (DG/BARK)·(FINT/5)**
  (the inside-bark diam_growth under-counts which trees reach the cap → too few killed);
  (c) HT cap — htgf.f:286-288 in height_growth!: if HT+HTG>SIZCAP[is,4] ⇒ HTG=max(SIZCAP[is,4]
  −HT, 0.1) (the 0.1 floor: trees already past the cap crawl, never shrink). Validated by
  test/integration/test_treeszcp.jl (3 scenarios vs Fortran, 106 asserts). Residuals: cap
  mid-cycle TPA/BA carry the regen response to cap-driven mortality (QMD bit-exact, endpoint
  matches); htcap TopHt drifts ≤4' as a declining-stand artifact (TPA/BA/QMD bit-exact).
- `FIXMORT` — ✅ normal path DONE (morts.f:1017). apply_fixmort! (keyword_dispatch.jl) overrides
  killed[] AFTER the BA-check (the last word on the kill), one-shot in the date's cycle, over a
  species×DBH window: IP 1 replace (P·rate), 2 add, 3 max, 4 multiply (kill·rate≤P), selected by
  PRM(5) (0/1/2/3), with Fortran's rate clamps. Needed a companion fix to mortality! ordering:
  **TPAMRT (the self-thinning line-reset, morts.f:772) is locked from the BA-check survivors
  BEFORE FIXMORT**, so the forced kill doesn't move next cycle's self-thinning line — without it
  the recovery ran TPA up to ~6% high. Bit-exact every cycle on 3 scenarios (replace, multiply,
  big-tree replace) vs live Fortran (test_fixmort.jl). ✅ **SIZE concentration DONE**
  (PRM(6)=10/20, morts.f:838-935): the window's mortality is pooled into XMORE per the IP, the
  in-window kills are zeroed where the IP replaces, then XMORE is re-imposed WHOLE-RECORD on the
  trees ranked by WORK3=∓(DBH+DG/bark) — KBIG=1 negates so the smallest grown trees go first
  (bottom up), KBIG=2 the largest (top down) — via the faithful `_rdpsrt!` (descending, FVS
  tie-break), killing each record fully (CREDIT+=tpa−killed) until CREDIT reaches XMORE (last
  record partial: killed+=XMORE−CREDIT). fixmort_big.key (pflag=10) is bit-exact vs live Fortran
  on TPA/BA/SDI/TopHt/QMD across all 11 cycles; the only diffs are ±1 in the volume cols (9-12)
  and per-acre growth (24-25) — the documented Float32 DGSCOR/volume noise, not the kill itself.
  ★ **KPOINT** point concentration (PRM(6)=1, morts.f:937) and **combined** size-within-point
  (PRM(6)=11/21, morts.f:978) now DONE & bit-exact (the base stand is 11-point, NOT single-plot —
  the earlier "no-op on IPTINV=1" note was wrong): the point walk uses `t.plot_id` (ITRE) +
  `points_inv` (IPTINV). Species groups (ISPCC<0) ✅ via SPGROUP + sp_field_matches.
- `FIXCW` — cwidth.f (crown-width override). ⚪ **OUTPUT-ONLY for the .sum** (verified): CRWDTH
  is referenced only by the calculator (cwidth.f), record bookkeeping that carries it along
  (comprs/tremov/triple), and OUTPUT consumers (sstage structure-class, svsnad SVS, evldx
  event-monitor var). It never feeds DGF/HTGF/MORTS/DENSE — those use crown RATIO, not width.
  So a FIXCW port changes no .sum growth number; defer until SVS/structure output is in scope.
- `HTGSTP` (HTGSTOP/TOPKILL) — ✅ DONE. htgstp! (keyword_dispatch.jl), called in grow_cycle!
  after TRIPLE/MORTS and before UPDATE (gradd.f:158). act 110 (HTGSTOP) scales HTG by PKIL;
  act 111 (TOPKILL) sets HT=H·(1−PKIL≤0.8), and for tall trees (H≥25, D≥6) whose Behre top
  diameter ≥4 marks a permanent broken top (NORMHT/ITRUNC) and cuts the crown ratio (ICR=−NEW).
  PKIL=BACHLO(AVEPRB,STDPBR), deterministic (=AVEPRB, no RNG) when STDPBR≤0; RANN escape when
  PRB<1; records walked in species-sorted IND1 order for RNG-exactness when stochastic. Needed a
  companion fix to crown_ratio_update! — the **negative-ICR bypass** (sn/crown.f:271): a crown
  already adjusted by topkill/pest models (ICR<0) is restored to +ICR and NOT recomputed that
  cycle; without it the top-killed trees' crown (hence DG/mortality) drifted and TPA ran ~10% high.
  Deterministic scenarios (HTGSTOP 0.5×, TOPKILL 0.5× >30') bit-exact every cycle vs live Fortran
  (test_htgstp.jl) AND a stochastic bare-plant scenario (STDPBR=0.2, htgstop_stoch) bit-exact through the firing cycle — the species-sorted IND1 RNG traversal is confirmed. IMC
  (management code) and ABIRTH (age) are set in Fortran but don't affect the .sum, so skipped.
- FIXDG/FIXHTG — ✅ DONE. grincr.f:451-525: DG/HTG·PRM(2) over a species×DBH window, applied
  in `apply_fix_scalers!` (keyword_dispatch.jl) after all growth / before MORTS. TWO things the
  earlier buggy attempt missed: (1) it is **ONE-SHOT** (OPDONE) — fires only in the cycle whose
  [start, start+period) range holds the keyword date, not every cycle (confirmed empirically:
  0.3× at 1995 drops QMD only in the 1995-cycle, then the gap persists ~constant, not runaway);
  (2) it must scale the **tripled** DG/HTG too — the stash dgU/dgL (htgU/htgL), matching FVS's
  DG(ITFN)/DG(ITFN+1). Reuses the GrowthMultiplier d1/d2 window. Bit-exact every cycle on 3
  scenarios (all/windowed DG, HTG) vs live Fortran (test_fix_scalers.jl). Species groups (ISPCC<0)
  ✅ via SPGROUP + sp_field_matches (test_spgroup.jl).

**Set-but-not-read in SN (0 application refs ⇒ likely NO-OP in SN, or external component):**
- CUTEFF, MINHARV, TCONDMLT — 0 refs.
- ⚠ CORRECTION: CRNMULT and TOPKILL were on this list but are NOT no-ops — both are applied in
  SN (CRNMULT scales the crown-ratio change in sn/crown.f:319; TOPKILL is the htgstp.f act 111
  top-kill). Both are now ✅ ported. The lesson: "0 application refs" from a coarse grep missed
  them because the grep keyed on the keyword name, not the COMMON variable (CRNMLT / IACT 111).
- SPLEAVE/LEAVESP — only `grinit.f:125 LEAVESP(I)=.FALSE.` (init), never checked in the cut logic.
- DEFECT/BFDEFECT/MCDEFECT — ★ CORRECTION (earlier note was WRONG): CFDEFT/BFDEFT ARE read, by
  **`bin/FVSsn_buildDir/vols.f`** — the REAL SN volume driver (the older R8 path that orchestrates
  the taper call via NATCRS and layers form+defect on top), NOT `vvolume/fvsvol.f` (the NVEL
  path). Verified empirically: a 50%-ish MCDEFECT slashes the live-Fortran merch-cubic `.sum`
  column. **MCDEFECT (cubic) now ported + bit-exact** (vols.f:294-332, ALGSLP over CFDEFT, SN
  pulpwood reduction; test_mcdefect.jl). REMAINING: BFDEFECT (vols.f:390-440 reduces BFV AND
  SCFV — needs the board-foot path), per-tree DEFECT input (digit-packed CF/BF/MC defect), and
  the CFLA0/CFLA1 / BFLA0/BFLA1 log-linear form model (default 0/1 = no-op, so latent until a
  species has non-default coefs — a separate gap that also affects un-keyworded stands).
- ⚠ CAVEAT: "0 refs" used a guessed effect-variable name; some may apply under a different
  COMMON name. Confirm empirically (does the keyword change the .sum?) before declaring no-op.

⇒ The cheap management wins (the 5 MULTS multipliers) are DONE. Every remaining item is a
focused chunk, not a quick port. Several listed-⛔ items are probably SN no-ops.

## Remaining work (current — 2026-06-23, after the growth/mortality surface + SPGROUP)

The entire growth/mortality keyword surface is ported & bit-exact: all CUTS methods (§1),
all MULTS multipliers, FIXDG/FIXHTG, FIXMORT, TREESZCP, HTGSTOP/TOPKILL, CRNMULT, SPGROUP
species groups. The DGSCOR cubic-volume drift is resolved as irreducible Float32 noise
(see DIVERGENCES.md §1). What is left, by kind of work:

**A. Genuinely-applied, .sum-affecting — real ports (ranked most-upstream → downstream):**
- ~~**Cycle calendar** — `TIMEINT`~~ 🟡 RESOLVED & IMPLEMENTED (uniform path). The hidden
  period-scaling mechanism was traced via a Fortran DEBUG dump: the DDS is FINT-independent, but
  dgdriv.f:325/715 applies SCALE=FINT/YR ⇒ DDS·(FINT/5) (YR=5 SN base); a 10-yr cycle ⇒ ×2 DG.
  Ported: diameter_growth! scales DDS by sfint/5; height by FINT/5; mortality^FINT; year/age step.
  Companion fixes: the morts grown-DBH G must NOT re-apply FINT/5 (diam_growth is already cycle-
  scaled); the calibration VMLT uses YR=5 not FINT; the .sum last-row year uses the real period.
  snt01 (period 5) bit-exact; TIMEINT-10 TPA ≤8, volume ≤2% (calibrated-species residual, like the
  DGSCOR tail). REMAINING: full YR-vs-FINT calibration split for bit-exactness; per-cycle GROWTH/
  CYCLEAT boundaries.
- ~~**Tripling control** — `NOTRIPLE` / `NUMTRIP`~~ ✅ DONE — wired through s.control.icl4
  (default 2; NOTRIPLE→0, NUMTRIP n→n). Was a real gap: FVSjl ignored NOTRIPLE (20 cols diverged)
  AND the bare-PLANT scenarios were silently passing for the wrong reason (FVSjl wrongly tripled
  them and the 3× averaging masked the regen DGSCOR tail at ±1). The fix exposed the true ±2 no-trip
  tail; the base+NOTRIPLE stand is bit-exact (test_tripling.jl).
- ~~ADDFILE~~ ⚪ unit-redirect include, not a tree-add (verified); no clean FVSjl mapping.
- ~~**COMPUTE**~~ ✅ DONE — kw_compute! stores NAME=expression defs (parsed with the event-monitor
  expression parser); cuts! evaluates them each cycle before the IF conditions, and _event_var resolves
  them; bit-exact (COMPUTE MYCYC=CYCLE ≡ direct CYCLE; firing thins match Fortran), test_compute.jl.
- **CYCLEAT / TIMEINT** scheduling boundaries (if not covered by the calendar item).
- ~~PRUNE~~ 🧊 .sum-inert in SN (verified) — C8 econ pruned-log volume only.
- **Volume overrides** — `VOLEQNUM`/`CFVOLEQU`/`BFVOLEQU` (🟡), `VOLUME`/`BFVOLUME`,
  `BFFDLN`/`MCFDLN`, `FIAVBC` (C5 volume side).
- **ESTAB TALLY** — tally-count regen (downstream, C4 regen). ~~SPROUT~~ ✅ DONE (ESUCKR, bit-exact).
- `COMPRESS` — only fires with the keyword or >~3000 records; comprs.f is a 762-line
  subsystem (incl. EIGEN!) — NOT least-dependent.

**B. Deferred sub-paths inside already-✅ items:**
- FIXMORT point/size concentration reallocation (PRM(6) — KBIG/KPOINT/combined, morts.f:838-1015) ✅ DONE.
- THINPRSC multi-plot (nps>1) path.
- ~~Species-group thin-method filter~~ ✅ DONE — `_cut_eligible` (and `_clsstk`/`_sdi_zeide`/
  `_rd_curtis`) take the SPGROUP table; the species-filtering thins (THINDBH/THINSDI/THINRDEN/
  THINCC/THINPT/THINQFA) thread `s.control.sp_groups`, so `THINDBH −1` resolves a group (bit-exact
  vs Fortran, test_spgroup.jl). REMAINING: group-**name** field refs (only −N numeric now).
- IF/THEN snt01 stand-2 3rd-thin class-boundary residual.

**C. Listed-⛔ but likely SN NO-OPS — confirm empirically (does it change the .sum?) first:**
- SPLEAVE/LEAVESP (only grinit init, never read in cut logic), CUTEFF / MINHARV / TCONDMLT
  (0 refs), MORTMSB (QMDMSB=999 inert), FIXCW (verified output-only),
  DEFECT/BFDEFECT/MCDEFECT (external NVEL lib; the keyword crashes this build — verify active).
  ⚠ "0 refs" mis-flagged CRNMULT+TOPKILL once; always confirm with a .sum diff.

**D. Out of scope here (C7/C8):** YARDLOSS (C7 fuel pools), FFE fire, insects, root disease
(RDIN/ANIN/RRIN), PRMFROST/CLIMATE, ECON/CHEAPO.
