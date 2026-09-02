#!/usr/bin/env python3
"""
addtrees_regen.py — a sample external regeneration model for the FVS ADDTREES keyword.

FVS's ESTAB `ADDTREES ... / <cmd>` keyword (estb/esaddt.f, activity 432, IMETH=1)
delegates regeneration to an EXTERNAL program: FVS writes a one-cycle stand summary to
`<stem>_<StandID>_<Year>_<Var>.es1`, invokes `CALL SYSTEM("<cmd> <thatfile.es1>")`, then
reads the program's `<...>.es2` back and OPADD-schedules the activities it contains
(430/PLANT and 431/NATURAL route straight into FVS's already-validated regen path).

This script is a faithful, self-contained stand-in for that external program. It reads the
exact `.es1` written by esaddt.f and emits a valid `.es2` activity block. The *decision
rule* here is a deliberately simple, documented sample (a real model — e.g. the Blue
Mountains estab model esaddt.f was built for — would be far richer); what matters for the
FVS contract is that the I/O formats round-trip exactly, which they do.

--------------------------------------------------------------------------------------
.es1 INPUT (written by esaddt.f, one value per line; see the WRITE order in esaddt.f):
   1  StandID            (A30)   stand identifier (NPLT)
   2  Year               (I30)   planting/target year  (IPYR = PRMS(1)+cycle_end_year)
   3  ADDTREES field-2   (I30)   INT(PRMS(1)) — the keyword's year-offset field
   4  KODTYP             (I30)   habitat-type code
   5  ISLOP              (I30)   slope (percent)
   6  IASPEC             (I30)   aspect (degrees)
   7  ELEV*100           (F30.1) elevation (feet)
   8  bsdi               (F30.1) before-harvest SDI / GROSPC
   9  asdi               (F30.1) after-harvest  SDI / GROSPC
  10  ba                 (F30.1) before-harvest basal area / GROSPC
  11  aba                (F30.1) after-harvest  basal area / GROSPC
  12  stand TPA >1"      (F30.1)
  13..20  species TPA <1" (F30.1) in esaddt.f order:
          ABGR, ABLA, LAOC, PICO, PIMO, PIEN, PIPO, PSME

.es2 OUTPUT (read by esaddt.f -> OPRDAT; this is exactly what a native PLANT card injects):
   line 1        IKEEP            (I10)  1 = keep the .es2 file, else FVS deletes it
   line 2        StandID                 must equal NPLT (OPRDAT keys the block on it)
   line 3..N     IACTK IDT NPRMS  PRMS…  one scheduled activity per line, list-directed
   last line     End                     terminates this stand's block
   (430 = PLANT, params: species  TPA  survival%  age  height  0   -> NPRMS=6)
--------------------------------------------------------------------------------------

Usage (exactly how FVS calls it):
    python3 addtrees_regen.py  RUN_STANDID_2001_IE.es1
It writes RUN_STANDID_2001_IE.es2 next to the input.
"""
import sys
import os

# esaddt.f species order for the <1" TPA block (fields 13..20) and the FVS species index
# (ISP code) each maps to. Used only by the sample decision rule below.
SPP = [
    ("ABGR", 4),   # grand fir
    ("ABLA", 9),   # subalpine fir
    ("LAOC", 2),   # western larch
    ("PICO", 7),   # lodgepole pine
    ("PIMO", 1),   # western white pine
    ("PIEN", 8),   # Engelmann spruce
    ("PIPO", 10),  # ponderosa pine
    ("PSME", 3),   # Douglas-fir
]


def read_es1(path):
    """Parse the .es1 written by esaddt.f into a dict. One field per line, in WRITE order."""
    with open(path) as fh:
        lines = [ln.strip() for ln in fh.read().splitlines()]
    # tolerate a trailing blank line
    lines = [ln for ln in lines if ln != ""]
    f = {}
    f["standid"] = lines[0]
    f["year"]    = int(float(lines[1]))
    f["field2"]  = int(float(lines[2]))
    f["kodtyp"]  = int(float(lines[3]))
    f["slope"]   = int(float(lines[4]))
    f["aspect"]  = int(float(lines[5]))
    f["elev"]    = float(lines[6])
    f["bsdi"]    = float(lines[7])
    f["asdi"]    = float(lines[8])
    f["ba"]      = float(lines[9])
    f["aba"]     = float(lines[10])
    f["tpa_gt1"] = float(lines[11])
    f["spp_lt1"] = {name: float(lines[12 + i]) for i, (name, _isp) in enumerate(SPP)}
    return f


def decide_regen(f):
    """
    SAMPLE decision rule (replace with a real regen model as needed).

    Idea: an opening (low residual basal area) triggers planting; species is picked from
    the site (elevation as a crude proxy), density scales with how open the stand is, and
    if the stand still carries advance regen of a species we bias toward it. Everything
    here is deterministic so the .es2 is reproducible for validation.

    Returns a list of PLANT activities: (species_index, tpa, survival_pct, age, height).
    """
    acts = []
    # Only plant into a genuine opening.
    if f["aba"] >= 40.0:
        return acts

    openness = max(0.0, 1.0 - f["aba"] / 40.0)   # 0 (closed) .. 1 (bare)
    tpa = round(100.0 + 300.0 * openness)         # 100..400 TPA

    # Species: prefer whatever advance regen is already strongest on the plot; else pick by
    # elevation band (ponderosa low, Douglas-fir mid, subalpine fir high).
    strongest = max(f["spp_lt1"].items(), key=lambda kv: kv[1])
    if strongest[1] > 0.0:
        isp = dict(SPP)[strongest[0]]
    elif f["elev"] < 3500.0:
        isp = 10   # PIPO
    elif f["elev"] < 5500.0:
        isp = 3    # PSME
    else:
        isp = 9    # ABLA

    acts.append((isp, float(tpa), 100.0, 2.0, 0.5))
    return acts


def write_es2(es1_path, standid, year, acts, ikeep=1):
    es2_path = os.path.splitext(es1_path)[0] + ".es2"
    with open(es2_path, "w") as fh:
        fh.write(f"{ikeep:10d}\n")     # line 1: IKEEP (I10)
        fh.write(f"{standid}\n")       # line 2: stand id (OPRDAT block key)
        for (isp, tpa, surv, age, ht) in acts:
            # 430 = PLANT ; NPRMS=6 ; params: species tpa survival age height 0
            fh.write(f"430 {year} 6 {isp} {tpa:g} {surv:g} {age:g} {ht:g} 0\n")
        fh.write("End\n")
    return es2_path


def main(argv):
    if len(argv) != 2:
        sys.stderr.write("usage: addtrees_regen.py <file.es1>\n")
        return 2
    es1_path = argv[1]
    f = read_es1(es1_path)
    acts = decide_regen(f)
    out = write_es2(es1_path, f["standid"], f["year"], acts, ikeep=1)
    sys.stderr.write(f"addtrees_regen: {os.path.basename(es1_path)} -> "
                     f"{os.path.basename(out)} ({len(acts)} PLANT activity/-ies)\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
