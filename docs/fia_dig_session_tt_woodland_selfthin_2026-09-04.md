# Dig session — TT M331D woodland "self-thin count-straddle" residual (2026-09-04)

Closes the final residual of the TT M331D woodland sweep cap: the **32 non-MM stands** left after the
MM(sp14) UTVAR-routing fix (master d4ec8f47). The prior dig log
(`fia_dig_session_tt_woodland_2026-09-02.md`, on the sweep branch) re-characterised these 32 as a
"sparse-woodland self-thin / density-mortality timing count-straddle" (QMD tied 29/32; BA 16hi/12lo;
CCF 21hi/11lo; worst offender CN 1629352222290487, CCF 127%). This session TreeId-matched per-tree
traced them against the freshly-relinked oracle `/workspace/.ttwork/FVStt_g16`.

## Method
`FVS_TreeList` (TREELIDB) dumps from oracle and jl, matched by TreeId per cycle; compared per-tree
TPA (survival) and DBH/Ht (growth). Instrumented the oracle (`smhtgf.f` fort.87 dump of ISPC/tree/
D/HT/CR/TPCCF/ZRAND/HTGRTH) and jl (`_tt_smhtgf` call-site print of the same) to read the SMHTGF
inputs and the ZRAND deviate directly.

## MEASURED VERDICT — NOT a self-thin/mortality bug. It is the SMHTGF ZRAND draw-order straddle.

1. **Survival (mortality / self-thin) is BIT-IDENTICAL.** For the worst stand 1629352222290487 the
   per-tree TPA matches the oracle every cycle 2022–2072 (ΣTPA 17.6→16.8 both sides, every record).
   The dig-log "self-thin / density-mortality TIMING count-straddle" premise is **refuted** — the
   density-mortality kernel (`teton/mortality.jl` label-10 self-thin + BA-check) is already faithful.

2. **The divergence is small-tree HEIGHT growth of the height-missing seedling records**, propagated
   to DBH via SMDGF and thence to BA/SDI/CCF/QMD/volume. In 1629352222290487 two ESTAB seedlings
   (FIA 19 = subalpine fir/AF, FIA 93 = Engelmann spruce/ES; input DIA=0.1, HT=missing) drive the
   entire CCF divergence: cycle-0 identical (both d0.1/h1.0), then at cycle 1 the oracle grows the ES
   seedling h1.0→1.5 (+0.5 ft) while jl grows it h1.0→8.5 (+7.5 ft) — a ~15× height over-prediction;
   by 2072 oracle DBH 0.31 vs jl 3.8.

3. **The SMHTGF inputs are bit-identical; only the ZRAND ±2 BACHLO deviate differs.** SMHTGF is
   `HTGRTH = HTG1 + ZRAND·STDDEV` with `STDDEV ≈ HTG1`, so a +2 vs −0.5 draw swings the increment ~6×.
   Measured at cycle 1 (worst stand), with CR=95, TPCCF=25 and start-height 1.01 matching on both sides:

   | seedling            | oracle ZRAND | jl ZRAND |
   |---------------------|--------------|----------|
   | AF (sp9, TreeId 1019004) | −0.637   | −1.576   |
   | ES (sp8, TreeId 1093004) | −0.764   | **+1.587** |

   The oracle deviates (−0.637, −0.764) both re-appear in jl's stream but assigned to **different**
   tree records (jl gave −0.764 to a later record). This is a pure draw-order / stream-position
   reassignment of the same BACHLO deviate sequence — not a model, coefficient, or input error.

4. **Root of the misalignment (structural, not a local bug):**
   - FVS calls `TRIPLE` from `grincr.f` — it draws ZRAND **once per original tree during growth,
     then triples**, so the 3 copies share one deviate (oracle: all AF copies −0.637, all ES −0.764).
     jl triples at initialisation and pre-draws ZRAND **per record**, so its copies diverge and its
     draw COUNT differs (6 vs 2 for two tripled seedlings) — desyncing the shared RNG stream.
   - Even the FIRST small-tree draw is already off (jl ES = +1.587 vs oracle −0.764), so the stream
     is desynced **before** the small-tree phase (accumulated crown-dub / calibration / mortality
     BACHLO draws), i.e. the classic OLDRN stream-position drift.
   - FVS also runs a small-tree HEIGHT HCOR calibration for **all** species (regent.f `DO 90 ISPC=1,
     MAXSP`); jl ports it only for aspen/MM (`tt_regent_hcor_aspen_init!`). It is inert on the worst
     stands (they carry < NCALHT=5 sub-5" measured-HTG conifers) and, ported faithfully, would itself
     have to replicate the same calibration-time ZRAND draws — re-entering the very stream alignment
     that is the corner. Only 5/32 stands even reach the ≥5 threshold.

## Corner
Named primitive: **the tt/smhtgf.f ZRAND (BACHLO/OLDRN) per-tree random deviate — draw-order /
stream-position serial-correlation straddle** (same class as #206). The aggregate over the 32 is a
balanced straddle (QMD tied, BA 16/12, CCF 21/11 — the convexity of BA/CCF in DBH accounts for the
mild high-lean under a mean-zero deviate). Aligning it would require reworking jl's shared small-tree
ZRAND draw sequence (per-original-before-triple + auditing every prior BACHLO draw), a high-regression
change to the RNG path shared by every currently bit-exact-or-cornered TT conifer stand, with no
measured net improvement — so it is cornered, not fixed. No engine change; gate stays 339/11.

Cornered: `docs/fia_cornered_clusters.tsv` (M331D · structure_densephase, smhtgf-ZRAND) + the 32 CNs
in `docs/fia_cornered_stands.txt`.
