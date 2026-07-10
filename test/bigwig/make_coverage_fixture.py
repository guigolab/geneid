#!/usr/bin/env python3
"""Regenerate the stranded bigWig coverage fixtures for the `rnaseq_bw` regression
case from the text expression GFF, using pybigtools (pip install pybigtools).

The GFF is dirRNAseq2geneid.sh step-3 output: `seqname . . start end depth strand . .`
(1-based inclusive coords, depth in col 6). Each strand becomes one bigWig of
0-based half-open [start-1, end) intervals -- exactly the coverage the text -S
path carries, so geneid's bigWig -S output is byte-identical to the text -S golden.

Usage: make_coverage_fixture.py [geneid-repo-root]   (default: ../.. from here)"""
import os
import sys

import pybigtools

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "..")
GFF = os.path.join(ROOT, "samples", "ENCFF001.1.MORC.stranded.expression.shuffled.gff")
CHROM, SIZE = "chr21", 46709983   # length of samples/chr21.fa

plus, minus = [], []
for ln in open(GFF):
    f = ln.split("\t")
    if len(f) < 7:
        continue
    s, e, v, strand = int(f[3]) - 1, int(f[4]), float(f[5]), f[6]
    (plus if strand == "+" else minus).append((CHROM, s, e, v))

for name, rows in (("plus", plus), ("minus", minus)):
    rows.sort(key=lambda r: (r[0], r[1]))
    for a, b in zip(rows, rows[1:]):
        assert a[2] <= b[1], f"overlapping coverage intervals: {a} {b}"  # bedGraph is a partition
    out = os.path.join(ROOT, "samples", f"ENCFF001.1.MORC.stranded.{name}.bw")
    w = pybigtools.open(out, "w")
    w.write({CHROM: SIZE}, iter(rows))
    w.close()
    print(f"wrote {out} ({len(rows)} intervals)")
