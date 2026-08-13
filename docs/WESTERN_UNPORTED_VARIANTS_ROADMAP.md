# Western Unported Variants — Port Roadmap

Scoping doc for the 10 unported **western** FVS variants (AK CA EC NC OC OP PN SO WC WS),
to be ported after the in-progress BC NEWSPRED (#196). All findings are **MEASURED from the
Fortran source** (variant dirs `/<v>/` + `bin/FVS<v>_buildDir/`), not inferred. No port code
was modified and the simulator was not run.

Repos: Fortran canonical `/workspace/ForestVegetationSimulator`; Julia port `/workspace/FVSjl`.

---

## Executive summary

**All 10 target variants share the same large-tree DG engine FVSjl already has:** the Wykoff
`LN(DDS)` diameter-growth model (`SUBROUTINE DGF(DIAM)`, verified byte-for-byte in every
`<v>/dgf.f` — same DDS structure with `DGLD/DGCR/DGBAL/DGDBAL/DGSITE/DGEL/DGCASP/DGSASP`
coefficient arrays as the ported CR/IE/EM/BM/CI). No dgf.f is byte-identical across variants —
**the differences are coefficient DATA, MAXSP, volume defaults, and a few special subsystems.**
The engine (Wykoff DDS, Reineke/Stage SDI self-thin, DUBSCR Weibull crown, REGENT small-tree,
shared NVEL volume) is already ported and reusable.

The port splits into **three clusters**:

1. **Westside R6 "Prognosis" family (7 variants): WC, PN, EC, CA, NC, SO, WS.** All share the
   exact same routine set as the already-ported **BM (Blue Mountains)** — `formcl.f + htdbh.f +
   sichg.f + ecocls.f + habtyp.f` + Wykoff-DDS `dgf` with aspect terms + DUBSCR crown +
   Reineke SDI + REGENT. Verified: `bm/formcl.f` vs `wc/formcl.f` differ **only in the national-
   forest coefficient tables** (same subroutine body). Cost = new coefficient blocks + one small
   variant-specific mortality/site routine each. **LOW–MED effort, batchable on one engine.**
2. **ORGANON-coupled (2 variants): OC, OP.** These carry the FVS-native Wykoff `dgf.f` too, but
   their primary integration **couples to the ORGANON growth model** (`orgspc.f` maps FVS species
   → ORGANON FIA codes, called from `DGDRIV`/`CRATET`; 45 ORGANON build refs vs 9–10 elsewhere).
   ORGANON is a **13,000-line separate growth subsystem** (`/organon` + `/vorganon`) not in FVSjl.
   **HIGH effort — the real cost driver.**
3. **Alaska (1 variant): AK.** Wykoff DDS **plus a unique PERMAFROST diameter-growth modifier**,
   a `SEAMRT` tolerance/percentile mortality, and R10 (Alaska) volume. Standalone, small
   (MAXSP=23). **MED effort.**

**Recommended order:** NC → WC → PN → {EC, CA, SO, WS} → AK → OC → OP. NC is the cheapest
(smallest species set, westside family, oracle already built at `/workspace/.ncwork` per memory);
WC is the group's structural anchor and porting it nearly gives PN for free (PN borrows WC's
`htgf/regent/cratet/dgdriv` byte-identically). OC/OP last because ORGANON is a whole new subsystem.

---

## Per-variant table

| Var | Region | MAXSP | DG model (routine) | Height / small-tree | Crown / SDI | Volume | Reuse anchor | Unique cost driver | Src lines¹ | Tier |
|-----|--------|-------|--------------------|---------------------|-------------|--------|--------------|--------------------|-----------|------|
| **NC** | Klamath Mtns / N. California | 12 | Wykoff LN(DDS) `nc/dgf.f` | REGENT + `htgr5.f` small-tree ht | DUBSCR / Reineke SDI | shared NVEL (R5/R6) | **BM** (westside set) + IE/EM Wykoff engine | `nwcmrt.f` mortality; `dunn.f` Dunning SI; smallest MAXSP | 11,016 | **LOW** |
| **WC** | West Cascades | 39 | Wykoff LN(DDS) `wc/dgf.f` | REGENT; `dgbnd.f` DG bound | DUBSCR / Reineke SDI | shared NVEL (R6) | **BM** (formcl/htdbh/sichg/ecocls identical structure) | group anchor; base `morts` | 10,227 | **MED** |
| **PN** | Pacific NW Coast | 39 | Wykoff LN(DDS) `pn/dgf.f` | **borrows WC** `htgf/regent/cratet/dgdriv` (byte-identical) | DUBSCR / Reineke SDI | shared NVEL | **WC** (near-clone) | own `dgf` coeffs + own `morts` only | 6,020² | **LOW** (after WC) |
| **EC** | East Cascades | 32 | Wykoff LN(DDS) `ec/dgf.f` | REGENT + `smhtgf.f` | DUBSCR / Reineke SDI | `bfvol` + NVEL | **BM/WC** (westside set) | `varmrt.f` mortality | 13,836 | **MED** |
| **CA** | Inland California / So. Cascades | 50 | Wykoff LN(DDS) `ca/dgf.f` | REGENT + `smhtgf.f` | DUBSCR / Reineke SDI | shared NVEL (R5) | **BM/WC** (westside set) | `varmrt.f`; `dunn.f` Dunning SI; large MAXSP | 11,984 | **MED** |
| **SO** | S. Central Oregon / Siskiyou | 33 | Wykoff LN(DDS) `so/dgf.f` | REGENT + `smhtgf.f`; `avht40/adjmai/maical` | DUBSCR / Reineke SDI | `bfvol` + NVEL | **BM/EC** (species "FROM EC") | `scomrt.f` mortality; `dunn.f` | 13,087 | **MED** |
| **WS** | Western Sierra Nevada | 43 | Wykoff LN(DDS) `ws/dgf.f` | REGENT + `smhtgf.f`; `dgbnd.f` | DUBSCR / Reineke SDI | shared NVEL (R5) | **BM/WC** (htdbh/sichg subset) | `varmrt.f`; `dunn.f`; no `formcl/ecocls` (lighter site) | 12,039 | **MED** |
| **AK** | Southeast Alaska | 23 | Wykoff LN(DDS) **+ PERMAFROST modifier** `ak/dgf.f` | REGENT; `htcalc/findag`; `esgent/esnspe` estab | DUBSCR / Reineke SDI (weighted SDIMAX) | **R10 (Alaska)** volume | none exact — Wykoff engine from CI/IE | **PERMAFROST DG modifier + `SEAMRT` mortality + R10 vol** | 11,905 | **MED** |
| **OC** | Oregon Coast / SW Oregon | 50 | Wykoff `oc/dgf.f` **+ ORGANON coupling** | REGENT / ORGANON | DUBSCR / Reineke; ORGANON crown | ORGANON vol (`getorgv`) + NVEL | none — needs new subsystem | **ORGANON growth model (SW-Oregon type), 13k lines** | 6,251² | **HIGH** |
| **OP** | Olympic Peninsula | 39 | Wykoff `op/dgf.f` **+ ORGANON coupling** | REGENT / ORGANON | DUBSCR / Reineke; ORGANON crown | ORGANON vol + NVEL | none — needs new subsystem | **ORGANON growth model (NW-Oregon type)** | 6,400² | **HIGH** |

¹ Line count of the `/<v>/` variant dir `.f` files only (shared base/common/volume excluded).
² OC/OP/PN dirs are small because they **borrow** most routines (PN from WC; OC/OP from ORGANON) at build time.

---

## Cluster 1 — Westside R6 "Prognosis" family: WC, PN, EC, CA, NC, SO, WS

**Evidence (presence matrix of the westside signature files):**

```
var  formcl htdbh  sichg  ecocls habtyp        <- shared westside routine set
bm    Y      Y      Y      Y      Y   (ANCHOR, already ported)
wc    Y      Y      Y      Y      Y
pn    Y      Y      Y      Y      Y
ec    Y      Y      Y      Y      Y
ca    Y      Y      Y      Y      Y
nc    Y      Y      Y      Y      Y
so    Y      Y      Y      Y      Y
ws    -      Y      Y      -      Y   (lighter: no formcl/ecocls)
```

**Structural proof:** `diff` of `bm/formcl.f` vs `wc/formcl.f` (ignoring comments/DATA) shows the
**same subroutine body**; the only differences are which national-forest form-class coefficient
tables are declared/loaded (BM: `MALHFC/OCHOFC/UMATFC/WLWHFC`; WC: `GIFPFC/MBSNFC/MTHDFC/…`).
Same story for `dgf.f` (identical DDS math, different coefficient blocks) and `crown.f`/`dubscr.f`.

**Shared, already-ported models this cluster needs:** Wykoff `LN(DDS)` DG, Reineke/Stage SDI
self-thinning, DUBSCR Weibull crown, REGENT small-tree/regeneration, shared **NVEL** volume
(the `vollib*/r5*/r6*` files present in every buildDir; FVSjl has western volume libs per memory).

**New per-variant work (small):** coefficient DATA blocks (dgf/formcl/htdbh/crown/regent) + one
variant-specific extra routine each:
- `varmrt.f` (variable mortality) — CA, EC, WS
- `nwcmrt.f` — NC · `scomrt.f` — SO · base `morts` — WC/PN
- `dunn.f` (Dunning site-index processing) — CA, NC, SO, WS
- `htgr5.f` (small-tree height) — NC · `smhtgf.f` — EC/CA/SO/WS

**PN special case:** `pn/` has no `dgdriv/htgf/regent/cratet/morts` of its own. Verified by
checksum that PN's buildDir copies of `htgf.f`, `regent.f`, `cratet.f`, `dgdriv.f` are
**byte-identical to `wc/`** (md5 match). PN = **WC clone** + its own `dgf.f` (coefficients) + own
`morts`. Once WC is ported, PN is a coefficient swap.

---

## Cluster 2 — ORGANON-coupled: OC, OP  (HIGH — the real cost driver)

`oc/orgspc.f` header: *"CONVERTS AN FVS SPECIES SEQUENCE NUMBER TO A VALID ORGANON SPECIES FIA
CODE … CALLED FROM SUBROUTINES CRATET AND DGDRIV."* OC uses the ORGANON **SW-Oregon** model type;
OP the NW-Oregon type. Build-integration depth: **45** ORGANON-referencing build files in
`FVSoc/op_buildDir` vs **9–10** in every other variant (those 9–10 are the shared
`getorgv`/tripling stubs that don't run ORGANON growth).

**Subsystem not in FVSjl:** `/organon` = 13,040 lines / 24 files (`org_intree`, `orgtab`,
`orgvol`, `orgtrip`, `statsorg`, `prepare`, …) + `/vorganon` = 1,331 lines. This is a whole
independent stand-growth model (diameter, height, crown, mortality, volume) with its own species
list and taper. OC/OP still carry a Wykoff `dgf.f` fallback, but a faithful port must reproduce
the ORGANON path. **Do these last.**

---

## Cluster 3 — Alaska: AK  (MED, standalone)

`ak/dgf.f` is Wykoff `LN(DDS)` **with a unique PERMAFROST diameter-growth modifier** (source:
`*** BEGIN PERMAFROST DIAMETER GROWTH MODIFIER COEFFICIENTS ***`, "ANNUAL DIAMETER GROWTH WITH
PERMAFROST MODIFIER", `PFCON` permafrost-intercept array + presence/absence factor). Mortality =
`SEAMRT` (distribute by percentile + species tolerance). Volume = **R10 (Southeast Alaska)**.
Establishment = `esgent/esnspe/esnutr/estock` (its own estab suite). No `formcl/htdbh/sichg/ecocls`
(not the westside site system). MAXSP=23. Engine reuse: the Wykoff DDS core from CI/IE; new work =
permafrost modifier + SEAMRT + R10 volume + AK estab.

---

## Effort tiers

| Tier | Variants | Why |
|------|----------|-----|
| **LOW** | NC, PN | NC = smallest MAXSP (12), westside family, oracle pre-built `/workspace/.ncwork`. PN = byte-identical WC clone. |
| **MED** | WC, EC, CA, SO, WS, AK | New coefficient blocks + 1 small routine each on the ported engine (WC/EC/CA/SO/WS); AK adds permafrost + SEAMRT + R10 vol but is standalone/small. |
| **HIGH** | OC, OP | Require the 14k-line ORGANON growth+volume subsystem — genuinely new. |

---

## Suggested port order

1. **NC** — cheapest, highest-leverage: MAXSP 12, westside Prognosis, oracle built. **Stands up
   and validates the westside engine reuse** off the BM anchor.
2. **WC** — the group's structural anchor (BM-identical routine bodies). Porting WC's coefficient
   blocks is the reusable template for the whole cluster.
3. **PN** — near-free after WC (borrows WC's `htgf/regent/cratet/dgdriv` verbatim; only `dgf`
   coeffs + `morts` are new).
4. **EC, CA, SO, WS** — batch on the same engine; each adds one mortality routine
   (`varmrt`/`scomrt`) + `dunn`/`smhtgf` coefficients. CA/WS are largest MAXSP (50/43).
5. **AK** — standalone MED; do after the westside batch since it shares no site system with them.
   Adds permafrost DG modifier, SEAMRT, R10 volume.
6. **OC**, then **OP** — HIGH; port the ORGANON subsystem once (shared `/organon`), then wire both.

---

## Models FVSjl already has vs genuinely new

| Already ported (reuse) | Genuinely new subsystem |
|------------------------|-------------------------|
| Wykoff `LN(DDS)` DG engine (CR/IE/EM/BM/CI/TT/UT) — **all 10 targets** | **ORGANON** growth+vol model (OC, OP) — 14k lines |
| Reineke/Stage SDI self-thinning; DUBSCR Weibull crown | **PERMAFROST** DG modifier (AK) |
| REGENT small-tree/regeneration; `smhtgf` potential-height | **R10 (Alaska)** volume library (AK) |
| Shared NVEL western volume libs (R5/R6, `bfvol/cubrds/logs`) | Minor: `varmrt/scomrt/nwcmrt/seamrt` mortality + `dunn` Dunning SI (small, per-variant) |
| Westside R6 site system (`formcl/htdbh/sichg/ecocls/habtyp`) via BM | — |

**ON (Ontario)** is eastern-Canada, not a western Prognosis variant (its own Penner-style model);
out of scope here and unrelated to this cluster.
