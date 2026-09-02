#!/usr/bin/env python3
"""
ihvext_selector.py — a sample external harvest-priority selector for FVS PPE MXHRVP
(the IHVEXT=1 option, hvin.f/hvaloc.f/hvsel.f, read back by sprdrd.f/sprdis.f).

When MXHRVP runs with external selection (IHVEXT=1), FVS does NOT compute the per-stand
harvest priority itself — it defers to an EXTERNAL optimizer (typically a road-access /
transportation model) and reads that optimizer's answer from a fixed file,
`PPE_FFERdAccess.txt`. HVALOC then assigns those values to HVPRI, and HVSEL selects
stands in descending priority, with one hard rule: a stand whose supplied priority is
<= 0 is NEVER selected (it is off-road / inaccessible this cycle).

This script is a faithful, self-contained stand-in for that optimizer. The FVS-side
contract is entirely the OUTPUT file (the input the optimizer consumes is its own
business), so the OUTPUT format is reproduced exactly and the decision rule is a simple,
documented sample you can replace.

--------------------------------------------------------------------------------------
INPUT (this sample's own format — a plain text stand table; comments/blank lines ignored):
      <StandID>  <standing_value>  <road_dist_mi>  [<accessible 0/1>]
  e.g.
      STANDA   120.0   0.4   1
      STANDB    80.0   2.1   1
      STANDC    55.0   9.0   0      # too far / no road -> priority forced <= 0

OUTPUT: PPE_FFERdAccess.txt, in the EXACT format sprdrd.f reads — `(A26, T30, F10.0)`:
      cols  1-26  StandID (left-justified, trimmed on read)
      cols 27-29  blank
      cols 30-39  priority value (F10.0 field; read as a float)
  a line whose id field contains `-999` terminates the file; a blank id line is skipped;
  a blank value field is the missing-value sentinel (-99999.0 on read).
--------------------------------------------------------------------------------------

Usage:
    python3 ihvext_selector.py  stands.txt              # writes ./PPE_FFERdAccess.txt
    python3 ihvext_selector.py  stands.txt  out_dir/    # writes out_dir/PPE_FFERdAccess.txt
"""
import sys
import os


def read_stands(path):
    """Parse the sample stand table. Returns list of (id, value, road_dist, accessible)."""
    rows = []
    with open(path) as fh:
        for raw in fh:
            line = raw.split("#", 1)[0].strip()   # strip inline comments
            if not line:
                continue
            parts = line.split()
            sid = parts[0]
            value = float(parts[1]) if len(parts) > 1 else 0.0
            dist  = float(parts[2]) if len(parts) > 2 else 0.0
            acc   = int(parts[3]) if len(parts) > 3 else 1
            rows.append((sid, value, dist, acc))
    return rows


def priority(value, road_dist, accessible):
    """
    SAMPLE priority rule (replace with a real transportation model as needed).

    Deliver more of the standing value per unit of hauling: priority falls off with road
    distance. An inaccessible stand gets a non-positive priority so HVSEL will never pick
    it (mirrors the real IHVEXT 'priority <= 0 => skip' rule). Deterministic for
    reproducible validation.
    """
    if not accessible:
        return 0.0
    return value / (1.0 + road_dist)


def write_rdaccess(rows, out_path):
    with open(out_path, "w") as fh:
        for (sid, value, dist, acc) in rows:
            p = priority(value, dist, acc)
            # (A26, T30, F10.0): id in cols 1-26, blank 27-29, value in cols 30-39.
            fh.write(f"{sid[:26]:<26}   {p:>10.1f}\n")
        fh.write("-999\n")   # terminator (sprdrd.f)
    return out_path


def main(argv):
    if len(argv) < 2:
        sys.stderr.write("usage: ihvext_selector.py <stands.txt> [out_dir]\n")
        return 2
    in_path = argv[1]
    out_dir = argv[2] if len(argv) > 2 else "."
    out_path = os.path.join(out_dir, "PPE_FFERdAccess.txt")
    rows = read_stands(in_path)
    write_rdaccess(rows, out_path)
    sys.stderr.write(f"ihvext_selector: {len(rows)} stands -> {out_path}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
