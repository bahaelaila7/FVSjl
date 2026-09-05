# NC (Klamath) LOC-705 Hoopa Region-7 BLMVOL volume bug — PAUSE + PLAN (2026-09-05)

Session: https://claude.ai/code/session_014RMPT9hw2GKinUa9hrorT7
Engine: wt/western-sweep = master d415384f + DVE-hardwood-vol fix 82020bc5.

## THE BUG (real, one-directional, deterministic, cycle-0) — NC is PAUSED on this
NC forest **IFOR 5 = HOOPA (LOCATION 705, Region 7)** — 57 stands in the population — has NO volume-routing
branch in jl, so it falls through to the Region-5 WO2W/DVE table (the fatter R5TAP taper) and mis-computes
volume: **TCuFt +16-22% OVER, BdFt -8-23% UNDER, at CYCLE 0 (structure bit-exact).** Deterministic, one-
directional, same TAXONOMY as the two volume fixes already landed this campaign (redwood broken-top cc24e470,
Siskiyou-R6 d415384f) — a forest/region VEQNNC routing gap on a crosswalk-exposed unvalidated forest.

Evidence (live FVSnc_clean A/B, cyc0):
- 7689975010901  (LOC705, red-alder/DF/tanoak): TCuFt jl 3160 vs O 2579 (+22.5%), BdFt 8446 vs 10958 (-23%). Structure (TPA97/BA163/TopHt7.9) bit-exact.
- 25058157010900 (LOC705, tanoak/bigleaf-maple): TCuFt 6579 vs 5646 (+16.5%), BdFt 27111 vs 29667 (-8.6%). Structure bit-exact.
- 301959964489998 (LOC705): flagged volume_persistent TCuFt +15.9%.

## ROOT CAUSE (traced to the NVEL dispatch)
NC's VEQNNC is forest/region-dependent (nc/sitset.f → VOLEQDEF(VAR='NC',IREGN=KODFOR/100,FORST)). LOC 705 →
IREGN=7, FORST='05'. The oracle VOLEQ table (FVSnc_clean .out on a 705 stand) is:
  all species → **B00BEHW<fia>**, DF → **B01BEHW202**.
`B`-prefix VOLEQ dispatches (volume/NVEL/volinit.f:366) to **BLMVOL → blmtap.f** (the BLM Behre's-hyperbola
stem-profile taper), NOT the R5 WO2W/DVE nor the R6 616BEHW Behre (nc_behre_vol). Form class for IFOR 5 =
**80** (formcl.f ELSE default: IFOR 5 is neither the SISKFC IFOR-4 nor the BLM712 IFOR-7 table). jl's
`compute_volumes_nc!` has branches for R5 (default table) and R6 (`isr6 = forest_idx==4`) but NONE for IFOR 5.

Note IFOR 7 (712 BLM Coos Bay) ALSO routes through BLMVOL (BLM712 form class) and jl's `isr6==4` misses it too
— but LOC 712 has 0 NC stands, so no current population impact (fix it for completeness alongside 705).

## THE FIX PLAN (port BLMVOL / blmtap.f — a new NVEL kernel, ~315 lines)
1. Port `nc_blmtap_dib(profile, D, H, D17, XLEN, htup)` = BLMTAP TLH=0 branch (Behre hyperbola):
   HBUTT=H-(XLEN+1.5); HTDIB=H-htup; A=BLMTHT(1..4,profile)·{1,D,H,D·H}; B=1-A; DIB=D17·(HTDIB/HBUTT)/(A·(HTDIB/HBUTT)+B).
   Tables BLMTHT(4,10) + BLMBA(6,9) are in blmtap.f:16-64. Port BLMTAPEQ (profile-index selection from VOLEQ)
   + the D17 / XLEN setup from the BLMVOL driver (read volume/NVEL/blmtap.f BLMVOL + blmtapeq).
2. Reuse the shared profile integration (`_fw2_tcubic` / `_fw2_merch_cuft`(VOL4+VOL7) / `_fw2_board`) with a
   `dibat = htup -> nc_blmtap_dib(...)` closure — the same pattern as `nc_wo2w_vol` (R5TAP) and the R6 path.
   BLM merch rules: MTOPP=6 IB (·bark), form class 80.
3. Add the IFOR-5 VOLEQ table `NC_R7_VOL_EQ` (B00BEHW<fia>, B01BEHW202 for DF) + an `isr7 = forest_idx==5`
   branch in `compute_volumes_nc!` (and `nc_snag_bole_cuft`) routing through nc_blmtap_vol; extend the R6
   check to include IFOR 7 (712) → BLMVOL with BLM712 form class.
4. Apply `r4_topkill` (CFTOPK/BFTOPK) for broken tops — METHC=6, same as the R5/R6/DVE paths.
5. VALIDATE per-tree vs an instrumented FVSnc_g16 (dvest BLMVOL dump, `--keywordfile` on a 705 stand) →
   cyc0 bit-exact TCuFt/BdFt on 7689975010901 / 25058157010900 / 301959964489998. Gate test_multicycle 339/11.

## STATUS
- DVE hardwood volume bug (the batch-1 cap) FIXED + validated → SHA **82020bc5** (branch fix-nc-dve-hardwood-vol).
- Dense-phase self-thinning count-straddle taxonomy (cursor 0-4000, 176 caps) cornered POP+SEED-verified
  (cyc0 honesty-gated 176/176 exact; seed-invariant SDImax knife-edge; pop two-sided). NC clusters added
  (263/261/M242 × count_divergence_UNVERIFIED). Synced to /workspace/FVSjl/docs/.
- Sweep PAUSED at cursor 6000 (batch 4000-6000, 238 caps) — of which ~dozens are the cornerable taxonomy and
  the LOC-705 stands are the REAL BUG above. **NC is NOT complete-to-7036: it is PAUSED on the 705 BLMVOL bug.**
