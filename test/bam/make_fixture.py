#!/usr/bin/env python3
"""Generate a small indexed BAM fixture + a JSON oracle of per-region coverage
for the bamcov reader test. Needs samtools on PATH (to sort+index); the oracle
is pure-stdlib. The .bam/.bai and oracle.json are committed so the C test needs
only a WITH_HTSLIB toolchain + python3 (not samtools) to run.

Coverage rule mirrors src/bamcov.c: CIGAR M/=/X add +1 per reference base; N and
D advance the reference without adding depth (reads split at introns/deletions).
Oracle stores non-zero coverage runs [start,end,depth] per query window."""
import json
import os
import random
import subprocess
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(__file__) or "."
CHROM, SIZE = "chr1", 200_000
random.seed(20260711)

# Reproducible reads over [100000,110000): plain 120M plus spliced 60M<N>60M so
# the split-at-N path (N advances ref, adds no depth) is exercised.
reads = []   # (pos0, [(op,len),...])
for _ in range(160):
    p = random.randint(100_000, 109_880)
    reads.append((p, [("M", 120)]))
for _ in range(24):
    p = random.randint(100_000, 106_000)
    gap = random.randint(400, 1500)
    reads.append((p, [("M", 60), ("N", gap), ("M", 60)]))


def covered(reads):
    depth = {}
    for pos, ops in reads:
        ref = pos
        for op, ln in ops:
            if op in ("M", "=", "X"):
                for x in range(ref, ref + ln):
                    depth[x] = depth.get(x, 0) + 1
            if op in ("M", "=", "X", "N", "D"):
                ref += ln
    return depth


depth = covered(reads)


def runs(lo, hi):
    out = []
    x = lo
    while x < hi:
        v = depth.get(x, 0)
        if v == 0:
            x += 1
            continue
        j = x
        while j < hi and depth.get(j, 0) == v:
            j += 1
        out.append([x, j, v])
        x = j
    return out


# SAM -> sorted BAM + index via samtools
def cig(ops):
    return "".join(f"{ln}{op}" for op, ln in ops)


sam = os.path.join(OUT, "_synth.sam")
with open(sam, "w") as f:
    f.write("@HD\tVN:1.6\tSO:unsorted\n")
    f.write(f"@SQ\tSN:{CHROM}\tLN:{SIZE}\n")
    for i, (pos, ops) in enumerate(reads):
        qlen = sum(ln for op, ln in ops if op in ("M", "=", "X", "I", "S"))
        f.write(f"r{i}\t0\t{CHROM}\t{pos+1}\t60\t{cig(ops)}\t*\t0\t0\t{'A'*qlen}\t{'I'*qlen}\n")

bam = os.path.join(OUT, "synth.bam")
subprocess.run(f"samtools sort -o {bam} {sam}", shell=True, check=True)
subprocess.run(f"samtools index {bam}", shell=True, check=True)
os.remove(sam)

queries = [(CHROM, 100_000, 110_000), (CHROM, 99_900, 100_200),
           (CHROM, 104_500, 105_500), (CHROM, 109_800, 110_100),
           (CHROM, 150_000, 151_000), ("chrX", 0, 1000)]
oracle = {"queries": []}
for c, s, e in queries:
    oracle["queries"].append({"chrom": c, "start": s, "end": e,
                              "runs": runs(s, e) if c == CHROM else []})
json.dump(oracle, open(os.path.join(OUT, "oracle.json"), "w"), indent=1)
print(f"wrote synth.bam (+.bai), oracle.json; {len(reads)} reads, max depth {max(depth.values())}")


# --- junction fixture: XS-tagged spliced reads (only these become introns) ----
# Each entry: (junction_start0, junction_end0, strand, n_reads). Two junctions
# share coords but differ in strand (distinct); one read carries two junctions;
# a few reads have NO XS tag and must be ignored.
juncs = [(102_000, 103_000, "+", 5),
         (102_000, 103_000, "-", 2),   # same coords, other strand
         (105_000, 105_500, "-", 3),
         (108_000, 108_400, "+", 1)]
jreads = []   # (pos0, cigar, xs_or_None)
for js, je, strand, k in juncs:
    for _ in range(k):
        jreads.append((js - 40, f"40M{je-js}N40M", strand))
# a read spanning TWO junctions (both '+') -> (102000,103000) and (104000,104500),
# plus 2 unspliced reads with no XS tag (must be ignored for junctions)
jreads.append((101_960, "40M1000N1000M500N40M", "+"))
jreads.append((103_500, "80M", None))
jreads.append((106_100, "80M", None))

jsam = os.path.join(OUT, "_junc.sam")
with open(jsam, "w") as f:
    f.write("@HD\tVN:1.6\tSO:unsorted\n")
    f.write(f"@SQ\tSN:{CHROM}\tLN:{SIZE}\n")
    for i, (pos, c, xs) in enumerate(jreads):
        mlen = sum(int(x) for x in __import__("re").findall(r"(\d+)[M=XIS]", c))
        tag = f"\tXS:A:{xs}" if xs else ""
        f.write(f"j{i}\t0\t{CHROM}\t{pos+1}\t60\t{c}\t*\t0\t0\t{'A'*mlen}\t{'I'*mlen}{tag}\n")

jbam = os.path.join(OUT, "synth_junc.bam")
subprocess.run(f"samtools sort -o {jbam} {jsam}", shell=True, check=True)
subprocess.run(f"samtools index {jbam}", shell=True, check=True)
os.remove(jsam)


# Build the junction oracle directly from the read set (authoritative), so the
# two-junction read and the no-XS reads are accounted for exactly.
def all_junctions(reads):
    import re
    out = []
    for pos, c, xs in reads:
        if xs not in ("+", "-"):
            continue
        ref = pos
        for ln, op in re.findall(r"(\d+)([MIDNSHP=X])", c):
            ln = int(ln)
            if op == "N":
                out.append((ref, ref + ln, xs))
            if op in ("M", "=", "X", "D", "N"):
                ref += ln
    return out


occ = all_junctions(jreads)
jqueries = [(CHROM, 101_000, 106_000), (CHROM, 102_000, 102_001),
            (CHROM, 107_000, 109_000), (CHROM, 120_000, 121_000), ("chrX", 0, 1000)]
joracle = {"queries": []}
for c, s, e in jqueries:
    agg = {}
    if c == CHROM:
        for js, je, st in occ:
            if s <= js < e:
                agg[(js, je, st)] = agg.get((js, je, st), 0) + 1
    joracle["queries"].append({"chrom": c, "start": s, "end": e,
                               "junctions": sorted([list(k) + [v] for k, v in agg.items()])})
json.dump(joracle, open(os.path.join(OUT, "oracle_junc.json"), "w"), indent=1)
print(f"wrote synth_junc.bam (+.bai), oracle_junc.json; {len(jreads)} junction reads")
