#!/usr/bin/env python3
"""Build the TM_ID -> (byte offset, row count) index into the TreeMap tree table.

The tree table is sorted ascending by TM_ID, so each plot's rows are contiguous.
Output: a compact binary index (uint32 count, then per plot: uint32 tm_id,
uint64 offset, uint32 nrows) enabling O(1) random access to any plot's tree list.

Usage:
    build_index.py <Tree_Table.csv> <out_index.bin>
"""
import struct
import sys
import time


def main(csv_path: str, out_path: str) -> None:
    t0 = time.time()
    idx: dict[int, list[int]] = {}  # tm_id -> [start_offset, nrows]
    with open(csv_path, "rb") as f:
        f.readline()  # skip header
        off = f.tell()
        line = f.readline()
        while line:
            tm = int(line[: line.index(b",")])
            e = idx.get(tm)
            if e is None:
                idx[tm] = [off, 1]
            else:
                e[1] += 1
            off = f.tell()
            line = f.readline()

    keys = sorted(idx)
    with open(out_path, "wb") as o:
        o.write(struct.pack("<I", len(keys)))
        for k in keys:
            so, nr = idx[k]
            o.write(struct.pack("<IQI", k, so, nr))

    total = sum(v[1] for v in idx.values())
    print(f"indexed {len(keys):,} plots, {total:,} tree rows in {time.time()-t0:.1f}s")
    print(f"tm_id range {keys[0]} .. {keys[-1]}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
