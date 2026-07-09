#!/usr/bin/env python3
"""Generate bigBed test fixtures + a JSON oracle of range-query results, using
pybigtools (pip install pybigtools). Run in an env that has it; the .bb files and
oracle.json are committed so the C bigbed decoder test does not need pybigtools.

Oracle format: {"queries": [{"chrom","start","end","records":[[start,end,rest...]]}]}
`records(chrom,start,end)` returns every feature OVERLAPPING [start,end)."""
import json
import sys

import pybigtools

OUT = sys.argv[1] if len(sys.argv) > 1 else "."

def build(path, chroms, rows):
    b = pybigtools.open(path, "w")
    b.write(chroms, iter(sorted(rows, key=lambda r: (r[0], r[1]))))
    b.close()

# fixture 1: small, multi-transcript, 2 chroms, name = transcript group (field 4)
chroms = {"chr1": 1_000_000, "chr2": 500_000}
rows = [
    ("chr1", 100, 200, "txA\t0\t+"),
    ("chr1", 300, 450, "txA\t0\t+"),
    ("chr1", 900, 1000, "txB\t0\t-"),
    ("chr1", 5000, 5200, "txC\t0\t+"),
    ("chr2", 100, 250, "txD\t0\t-"),
]
build(f"{OUT}/small.bb", chroms, rows)

# fixture 2: many records so the R-tree has multiple leaves + multiple data blocks
chroms2 = {"chr1": 10_000_000}
rows2 = [("chr1", i * 1000, i * 1000 + 500, f"t{i}\t{i % 1000}\t{'+' if i % 2 else '-'}")
         for i in range(2000)]
build(f"{OUT}/many.bb", chroms2, rows2)

# oracle: a spread of range queries (edge overlaps, spanning, empty, cross-chrom)
queries = [
    ("chr1", 120, 350), ("chr1", 0, 100), ("chr1", 200, 300), ("chr1", 950, 5100),
    ("chr2", 0, 1000), ("chr2", 300, 400), ("chr1", 999999, 1000000),
]
oracle = {"files": {}}
for fn in ("small.bb", "many.bb"):
    b = pybigtools.open(f"{OUT}/{fn}", "r")
    qs = queries + [("chr1", 1_500_000, 1_520_000), ("chr1", 4_000_000, 4_002_000)]
    res = []
    for (c, s, e) in qs:
        if c not in b.chroms():
            continue
        try:
            recs = [list(r) for r in b.records(c, s, e)]
        except Exception:
            recs = []
        res.append({"chrom": c, "start": s, "end": e, "records": recs})
    oracle["files"][fn] = {"chroms": b.chroms(), "queries": res}

json.dump(oracle, open(f"{OUT}/oracle.json", "w"), indent=1)
print("wrote small.bb, many.bb, oracle.json to", OUT)
