# PN (Pacific Northwest) variant port — audit

PN is a westside R6 Wykoff-DDS **near-clone of West Cascades (WC)** on the SAME shared engine: MAXSP=39,
Reineke/Stage SDI (LZEIDE=.FALSE.), DGSD=1.7, IFINT=10/IFINTH=5, RNG seed 55329. It differs from WC only
in the coefficient DATA + the species table (slot 6 = **SS/Sitka spruce** FIA 098, vs WC's blank). Per
doctrine, PN gets its OWN consts/functions so WC stays provably inert.

## Oracles
- `/workspace/.pnwork/FVSpn_clean` (relink bin/FVSpn_buildDir/*.o + isoc23 shim; `relink_pn.sh`).
- `/workspace/.pnwork/FVSpn_g16` (instrumentable; `build_g16.sh`). Canonical stand `pnt01` = the SAME
  tree data as wct01 (S248112) but forest **612 = Siuslaw → IFOR 2**, elev 7.0, habitat 40 (vs wct01's
  618/35.0/52), so site/DGCON inputs differ. `/workspace/.pnwork/pnt01.sum.save` reference.

## Validation — chunk 0/1/3 DGF (2026-08-13, DGCON bit-exact)

`pacificnorthwest.jl` (singleton + registration) + `species.jl` (39-species table, slot-6 SS) +
`diameter_growth.jl` (`pn_dgcons!` + `dgf!`). **pn/dgf.f has 20 GROUPS (WC had 19)** — SS is a NEW group
18 ⇒ WO→group 19, RW→group 20. Differences from `wc_dgcons!`: WO King's-SI XSITE transform is JSPC==19;
**NO JFOR (BLM→NF) remap** (pn/dgf.f:510 uses IFOR directly); DGFOR is 3-loc-class (20×3); DGDS col-2
nonzero for g7/g8/g9/g18. DEFAULT DDS + RA(g13) + RW(ISPC 17) branches byte-identical to WC.

> **PN DGCON 39/39 BIT-EXACT** vs FVSpn_g16 (instrumented ENTRY-DGCONS dump), worst |Δ|=4.2e-7. Validates
> every PN_* coefficient array + the DGCON formula. Harness `dgcon_validate.jl`. cyc1 LN(DDS) reference
> saved (`ref_dgf_pnt01_cyc1.txt`) for the end-to-end check once site/density land.

## REUSE MAP (from the coefficient-extraction agents) — the remaining chunks

| chunk | PN status | source |
|---|---|---|
| 2 site index | TODO — pn/forkod.f (JFOR=609,612,800,708,709,712→IFOR 1-6; 612→2), pn/habtyp.f PCOML + pn/ecocls.f ECOCLS (DIFFER from WC — extract pcoml/ecocls.csv), pn/sichg.f (in species CSV), pn/sitset.f FMSDI + site species. SDIDEF. | pn/*.f |
| density CCF | TODO — pn/ccfcal.f is **19 groups** (WC 16): INDCCF[39]/RD1/RD2/RD3 provided by agent (need pn_tree_ccf + stand_ccf/point_density PN branch). | pn/ccfcal.f |
| 4 HTG | TODO — check pn/htgf.f vs WC (findag/htcalc; htdbh differs). | pn/htgf.f |
| 5 crown | TODO — pn/crown.f is **17 groups** (WC 16; SS added): WEIBA/B0/B1/C0/C1 + crc0/crc1 provided by agent. | pn/crown.f |
| 6 regent | REUSE WC SMHGDG (vwc/smhgdg.f byte-identical) — but **NO DF Curtis→King** (VARACD≠'WC' ⇒ raw SITEAR). PN htdbh (htdbh_coeffs_pn.csv, 4 tables). | vwc/smhgdg.f + pn/htdbh.f |
| 7 mortality | REUSE WC ORGANON arrays (morts.f byte-identical) — but **XSITE1 = raw SITEAR(16)** (no King). | vwc/morts.f |
| 8 volume | REUSE shared R6_EQN (voleqdef.f byte-identical; forest '12' Siuslaw) + PN formcl (formcl_pn.csv) + merch (non-BLM = WC defaults). | shared + pn/formcl.f |

**Data files landed** (data/pacificnorthwest/): species_coefficients.csv (39), species_translation.csv,
htdbh_coeffs_pn.csv (234), formcl_pn.csv (1170), fire_species_props.csv (39), merch_specs.csv (39).
**Still to extract**: pcoml.csv + ecocls.csv (site chunk). **VARACD caveat** (agent): the DF Curtis→King
SI conversion in SMHGDG/mortality fires only for VARACD='WC' — PN must use RAW SITEAR(16). **Wiring TODO**:
setup dispatch (simulate.jl PN branch) + the shared-engine bark branches (_pn_bd/_pn_cal/_pn_dg/_pn_up —
same POWER wc_bratio, since PN bark is also imap=1) + stand_ccf/point_density/dub_missing_heights/
init_merch/setup_volume_equations PN branches. Then run pnt01 end-to-end and validate cyc1 LN(DDS) 27/27
+ the .sum (the WC end-to-end pattern; watch for the same COR-calibration/DDS→DG bark class that hit WC).
</content>
