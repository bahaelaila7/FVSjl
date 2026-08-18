# DFB mortality path — handoff (branch dfb-mort, HEAD f9a06d43)

## DONE (bit-exact-or-cornered, committed f9a06d43 off 0d7c8d11)
DFBRAN + BACHLO + DFBMOD + DFBMRT + gated DFBDRV seam + DFBGO MANSTART/MANSCHED gate.
- Numeric core BIT-EXACT via dump-replay vs relinked FVSie_dfb (see below).
- End-to-end .sum-DELTA CORNERED by the documented IE cycle-1 growth/mortality straddle.
- Real FVS keyword is **DFB** (keywds.f opt 100), NOT "DFBEETLE" — dispatch fixed.

## Oracle relink recipe (IE)
`build_ie_dfb.sh` — builds the 20 real dfb/*.o for IE (all routines + dfblkdie, DROP
dfbincr.f which is a legacy duplicate of dfbin/dfbkey) + wsbwe/txnote.o, relinks
FVSie_dfb = /workspace/.iework/g16obj/*.o (minus exdfb.o) + the 20 dfb objs + IE *.c +
isoc23_shim. `relink_only.sh [suffix]` relinks from an already-built ieobj/ (so an
instrumented .o copied over its pristine name survives). Clean DFB-off ≡ stock FVSie_g16
(verified). Binaries were removed (8.6M each) — rerun build_ie_dfb.sh to regenerate.

## Mature-DF stand + run
iet01 (S248112, /workspace/.iework/ierun/iet01.tre, 30 recs, 8 DF, 4 DF>=9"): the
validated IE reference. Keyfiles in orun/: ictrl.key (control), idfb.key (DFB MANSTART +
MANSCHED 2 -> cycle-2 outbreak). Run: `printf 'idfb.key\niet01.tre\n' | FVSie_dfb`.
(orun/dfb5.key is the same on a 12-record synthetic subset — that subset hits a cornered
FVSjl growth regime, use iet01.)

## Golden dump (instrumented g16, MANSTART+MANSCHED, seed 56457) — in test_dfb.jl
DFBRAN stream: 3EE23A99 3E6A6620 3E56B516
DFBMOD (iet01 cyc2): BADF9=41B78CB8 BA9=42824F2E NUMYRS=4 -> DFKILL=412518D2 (10.31856)
DFBMRT SUMDBH=43A53859; DF>=9 TAMORT: 4031CDC0 40008AB8 401C9AF1 40456FE1
All reproduced BIT-EXACT in Julia (validate_core.jl).
BACHLO here is the UNIFORM branch (U=0.4419<=2/3) so the returned X=1.5*U carries no
transcendental — DFKILL is bit-exact regardless of log ULP.

## Cornering detail (measured)
FVSjl cyc2 BADF9/BA9 = 0.299 vs oracle 0.352 (DF>=9 basal area differs after cycle-1
growth) => DFKILL 8.76 vs 10.32 => .sum MOR delta +15 (FVSjl) vs +20 (oracle) at the
2000 row. NOT a DFB bug — the IE growth/volume straddle (RESULTS_WESTERN "forest-keyed
VOLEQDEF fixed-table gap"; CUFT ~0.75% cyc0, TPA straddle by cyc2). DFB math is exact on
identical inputs.

## Instrumentation files (scratchpad/dfb, for re-validation)
dfbran_instr.f (dumps S1/SEL), dfbmod_instr.f (DFKILL/BADF9/BA9/NUMYRS),
dfbmrt_instr.f (SUMDBH + per-record D/P/DCLAS/TAMORT/WK2), dfbgo_instr.f (LDFBON/LMIN/
NTODO/ISMETH/PROTBK gate). Copy over ieobj/<name>.o then relink_only.sh _instr.

## NEXT sub-chunks (NOT done)
1. RANSTART (ISMETH=2) stochastic inclusion: dfb_apply! already draws dfb_rand!<protbk;
   validate the DFBRAN inclusion-draw ORDER (it precedes BACHLO) + STOPROB/PROTBK path.
2. RANSCHED (IDBSCH=2) DFBSCH auto-scheduler: the DBEVNT probability loop + the
   OPNEW(2208) schedule + the ORSEED+1128 reset AFTER the loop (currently the reset is
   assumed unconditional, correct for MANSCHED). Port dfbsch.f, validate the schedule.
2b. OPNEW/OPFIND/OPDONE activity-2208/2209 seam if a turnkey exercises the regional
    scheduler (currently modelled directly via dfb_outbreak_due / mansched_years).
3. CUROUTBK / LINPRG in-progress-outbreak DFBMOD branch (PERDD/PREKLL/IYOUT) + DFBINV
   (pre-killed-DF count from treelist damage codes) — needed if a turnkey uses CUROUTBK.
4. DFBWIN windthrow (activity 2210, OKILL feed into DFBMRT) — WINDTHR keyword.
5. Variant sweep: exercise DFB on BC/BM/CI/CR/EC/PN/TT/UT (idfspc map in place; PN=16).
6. Tighter end-to-end: find/relink a stand where FVSjl IE control ≡ oracle control for a
   bit-exact (not cornered) .sum-DELTA, or accept the IE growth-straddle cornering.
