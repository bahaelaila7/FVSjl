# Sample external programs for FVS's out-of-tree hooks

Two FVS features delegate part of the simulation to an **external executable** that FVS
shells out to. FVSjl ports the FVS *side* of each contract bit-exactly (the file it writes
and the file it reads back), but the external program itself lives outside the model. These
two Python scripts are faithful, self-contained **stand-ins** for those programs: they honor
the exact on-disk I/O contract, with a simple documented decision rule you can replace with
a real model.

Both are validated end-to-end against the committed FVSjl paths (no oracle needed) — see
`test/…` and the round-trip below.

## 1. `addtrees_regen.py` — external regeneration model (ESTAB `ADDTREES`, esaddt.f, act 432)

FVS writes a one-cycle stand summary to `<stem>_<StandID>_<Year>_<Var>.es1`, runs
`SYSTEM("<cmd> <that.es1>")`, then reads the program's `.es2` activity block back
(`OPRDAT`) and schedules its activities (430=PLANT / 431=NATURAL → the validated regen
path). The script reads the `.es1` (field order per esaddt.f) and emits a `.es2`.

```
# FVS keyfile:
ADDTREES        1992         0         1
python3 /abs/path/examples/external_programs/addtrees_regen.py
```

The FVSjl `addtrees_bridge!` (`src/engine/establishment.jl`) invokes this exactly as FVS's
`CALL SYSTEM` would. Verified: bridge → `addtrees_regen.py` → `.es2` (`430 … 3 400 100 2
0.5 0`) → FVSjl schedules the PLANT and the regen appears in the `.sum`.

## 2. `ihvext_selector.py` — external harvest-priority selector (PPE MXHRVP `IHVEXT=1`)

With external selection on, MXHRVP does **not** compute per-stand harvest priority itself —
it reads an external optimizer's answer from `PPE_FFERdAccess.txt` (`sprdrd.f`/`sprdis.f`),
in `(A26, T30, F10.0)` columns. HVSEL then cuts stands in descending priority, and **never**
selects a stand whose supplied priority is ≤ 0 (off-road / inaccessible). The script reads a
plain stand table and writes that file.

```
python3 ihvext_selector.py stands.txt [out_dir]   # -> out_dir/PPE_FFERdAccess.txt
```

Verified: selector → `PPE_FFERdAccess.txt` → `ppe_read_rdaccess` yields the exact priority
map, and a ≤0 priority excludes that stand from selection.

## What is a "sample" here

The **I/O contract** (file names, formats, field order, terminators, sentinels) is faithful
and is what FVS/FVSjl actually depends on. The **decision rule** inside each script (which
species to plant, how to rank stands) is an illustrative placeholder — a production model
(e.g. the Blue Mountains estab model, or a transportation optimizer) drops in behind the same
I/O with no change to FVS.
