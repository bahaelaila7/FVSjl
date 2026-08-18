# WPBR (White Pine Blister Rust) — chunk 0 goldens + oracle recipe

## Model identity (scope)
- Fortran: `wpbr/` (43 .f + BRCOM.F77). Keyword = **BRUST**, keywds.f option **75**.
- RNG: `brann.f` — MINSTD LCG `BRS1=DMOD(16807·BRS0, 2147483647)`, `SEL=BRS1/2^31`,
  seed **55329** (odd). ENTRY `BRNSED(LSET,SEED)` reseeds (AMOD even→+1). IDENTICAL
  to DFB DFBRAN / DFTM TMRANN — the 6-draw stream MATCHES `_G_TMRANN` bit-for-bit.
- Per-cycle driver: **BRTREG** (gradd.f:124 "begin simulation"). Setup: BRSETP
  (fvs.f:351, from MAIN, per-tree ground diam). Init: BRINIT (from INITRE).
  Other entry points (exbrus.f stub): BRDAM (dampro), BRCMPR (comprs), BRESTB
  (estab), BRTDEL (tredel), BRSOR (tresor), BRTRIP (triple), BRPR/BRROUT (output).
- Randomness (CALL BRANN) sites: BRECAN (expected cankers), BRCINI (random canker
  up/out/girdle), BRCGRO (mortality realization BRPB), BRCREM (prune/excise success),
  BRSTYP (stock-type assignment).
- Linkage: NO shipped binary runs WPBR (all link base/exbrus.o stub). Variant block
  data exists for base/NI, IE (brblkdie), CR (brblkdcr), SO (brblkdso) — all named
  `BLOCK DATA BRBLKD`, so a relink includes exactly ONE.
- Host maps (BRSPM): IE = WP@1, LM@13 (BRSPC WP,LM); SO = WP@1, SP@2 (WP,SP);
  CR = LP@11, PP@13 (LP,PP); base/NI = WP@1 (WP,SP). NBRSP=2. IBRDAM=36.

## RNG goldens (scratchpad/wpbr/driver_brann.f over pristine brann.f, gfortran-16)
Stream (seed 55329, 6 draws), Float32 IEEE hex:
    3EDDB57A 3F5AB21B 3F630A07 3F2734C9 3EF4FB2C 3F4AF646
BRNSED reseed:
    even 100 → forced-odd 101, draw1 = 3A4F3718
    odd 55, draw1 = 39E1AE10, draw2 = 3E70354E
    LSET=false reset-to-SS(=55), draw = 39E1AE10  (== odd-55 draw1)
All reproduced bit-exact by Julia wpbr_rand!/wpbr_seed! (test_wpbr.jl).

## Oracle relink (reusable for the dynamics chunks)
- `scratchpad/wpbr/build_ie_wpbr.sh` → `FVSie_wpbr` (IE .o set, exbrus.o swapped for
  real wpbr/*.o; brblkd/cr/so + wpbr/exbrus excluded; exppe stub added; isoc23 shim).
  LINK_OK, 39 wpbr objects.
- Turnkey: `scratchpad/wpbr/run_oracle.sh`. FVS CLI reads stdin line1 = keyword file,
  line2 = **separate tree-data file** (unit 02). Inline TREEDATA is IGNORED — the .tre
  MUST be a separate file. (This bit me first: inline records → 0 TPA.)
- VALIDATED: WPBR-off relink ≡ stock FVSie .sum BYTE-IDENTICAL (clean relink).
  WPBR-on (BRUST + RUSTINDX 0.05) ENGAGES: real canker mortality at 2030/2040
  (TPA 25→17, MOR 10→49, BA 45→35 by 2040) — the model is LIVE.

## Fidelity caveats for the dynamics chunk
- brin.f DEVFACT/STOCK index DFACT/PRPSTK/RESIST(KSP,·) by the RAW FVS species code,
  but those arrays are dimensioned (NBRSP=2,4) — a latent Fortran overrun. Only
  coincides with the BR host index for host FVS-codes ≤ NBRSP (WP@1, SP@2). Diverges
  for LM@13 / PP@13. The Julia port indexes by the in-bounds BR host index (matches
  the correct GROWRATE pattern). Inert-seam ⇒ no .sum effect NOW; audit when the
  dynamics chunk actually reads DFACT/PRPSTK/RESIST.
- RUSTINDX method 3/4 and CANKDATA read supplemental records — chunk-0 consumes the
  file-name/param records but does NOT open canker files or run BRCANK/BRIBES.
- OPNEW-scheduled activities (PRUNE/EXCISE/RIBES/DEVFACT/STOCK/PRNSPECS-idt>0/…) are
  recorded in w.activities but NOT acted on (no engine seam).

## Next sub-chunk (dependency-ordered)
1. BRSETP + BRINIT per-tree init (BRGD ground diameter, BRHTBC crown-base) — deterministic,
   dump-replay vs FVSie_wpbr. Needs the host-stand load recipe (separate .tre).
2. BRTARG/BRIBA/BRIBES/BRICAL — rust index (RI) + DFACT/RIAF + growth index (BGRI/BRGI)
   assignment. Mostly deterministic (RI from BA/ribes/age); the RIAF one-time draw.
3. BRECAN expected-canker generation (first BRANN draws in the cycle) — stream order.
4. BRCANK/BRCINI canker initialization (random up/out/girdle) — stream order.
5. BRCGRO canker growth → BRPB mortality + BRCSTA/BRTSTA status → top-kill/mortality.
6. BRTREG driver seam (gradd.f:124) wiring into simulate.jl (after MORTS/MISTOE, like
   the DFTM TMCOUP seam) + BRDAM damage codes + .sum-DELTA bit-exact-or-cornered.
