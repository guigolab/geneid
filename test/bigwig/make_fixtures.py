#!/usr/bin/env python3
"""Generate bigWig test fixtures + a JSON oracle of range-query results, using
pybigtools (pip install pybigtools). Run in an env that has it; the .bw files and
oracle.json are committed so the C bigwig decoder test does not need pybigtools.

Ground truth = pybigtools `values(chrom, start, end, fillna=0)` (per-base signal,
uncovered bases = 0). Stored run-length-encoded (non-zero runs) to keep the
oracle small; check.py expands both the oracle and the C decoder's intervals to
per-base arrays over the query window and compares them."""
import json
import sys

import pybigtools

OUT = sys.argv[1] if len(sys.argv) > 1 else "."


def build(path, chroms, rows):
    b = pybigtools.open(path, "w")
    b.write(chroms, iter(sorted(rows, key=lambda r: (r[0], r[1]))))
    b.close()


# fixture 1: small, 2 chroms, gaps + adjacent runs (bedGraph sections)
chroms = {"chr1": 1_000_000, "chr2": 500_000}
rows = [
    ("chr1", 100, 200, 1.0),
    ("chr1", 200, 300, 2.5),    # adjacent, different value
    ("chr1", 900, 1000, 7.0),   # gap before
    ("chr1", 5000, 5200, 3.0),
    ("chr2", 100, 250, 4.0),
]
build(f"{OUT}/small.bw", chroms, rows)

# fixture 2: many contiguous runs -> multiple (zlib) data blocks
chroms2 = {"chr1": 10_000_000}
rows2 = [("chr1", i * 1000, i * 1000 + 600, float((i % 50) + 1)) for i in range(2000)]
build(f"{OUT}/many.bw", chroms2, rows2)

# fixture 3: many chroms (one section each) -> a MULTI-LEVEL R-tree, so the query
# recurses through internal (non-leaf) nodes, not just a single leaf root.
chroms3 = {f"c{i:04d}": 1000 for i in range(400)}
rows3 = [(f"c{i:04d}", 100, 300, float(i % 7 + 1)) for i in range(400)]
build(f"{OUT}/deep.bw", chroms3, rows3)


def rle(vals, base):
    """Non-zero runs of a per-base array -> [[start, end, value], ...]."""
    runs = []
    j, n = 0, len(vals)
    while j < n:
        v = float(vals[j])
        if v != v or v == 0.0:            # NaN or 0 -> gap
            j += 1
            continue
        k = j
        while k < n and float(vals[k]) == v:
            k += 1
        runs.append([base + j, base + k, v])
        j = k
    return runs


queries = [
    ("chr1", 0, 100), ("chr1", 150, 250), ("chr1", 200, 300), ("chr1", 199, 201),
    ("chr1", 850, 1050), ("chr1", 4900, 5300), ("chr2", 0, 400), ("chr2", 200, 260),
    ("chr1", 1_500_000, 1_520_000),   # empty region
]
# queries into the multi-level tree: several chroms across the index (hit early,
# middle, late leaves), an edge overlap, and an empty window
deep_queries = [("c0000", 0, 400), ("c0007", 150, 250), ("c0100", 299, 301),
                ("c0200", 0, 100), ("c0399", 100, 300), ("c0250", 400, 500)]

oracle = {"files": {}}
for fn in ("small.bw", "many.bw", "deep.bw"):
    b = pybigtools.open(f"{OUT}/{fn}", "r")
    if fn == "deep.bw":
        qs = deep_queries
    else:
        qs = queries + [("chr1", 1_499_000, 1_502_000), ("chr1", 4_000_000, 4_003_000)]
    res = []
    for (c, s, e) in qs:
        if c not in b.chroms():
            continue
        vals = b.values(c, s, e, fillna=0)
        res.append({"chrom": c, "start": s, "end": e, "runs": rle(vals, s)})
    oracle["files"][fn] = {"chroms": b.chroms(), "queries": res}

json.dump(oracle, open(f"{OUT}/oracle.json", "w"), indent=1)
print("wrote small.bw, many.bw, oracle.json to", OUT)
