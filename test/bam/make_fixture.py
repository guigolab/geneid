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


# --- junction fixture: spliced reads exercising bamJunctionQuery's strand
# resolution -- XS (genomic strand, used directly), minimap2 ts (relative to the
# read, flipped by orientation), and no tag (emitted as '.'). Each read is
# (pos0, cigar, flag, tag) with tag = ("XS",s) | ("ts",s) | None. Motif inference
# is NOT tested here (it lives above bamJunctionQuery, in ReadIntronsBam).
def spliced(js, je, flank=40):
    return js - flank, f"{flank}M{je-js}N{flank}M"

jreads = []
# XS-tagged: two junctions share coords but differ in strand (distinct)
for js, je, s, k in [(102_000, 103_000, "+", 5), (102_000, 103_000, "-", 2),
                     (105_000, 105_500, "-", 3), (108_000, 108_400, "+", 1)]:
    for _ in range(k):
        p, c = spliced(js, je)
        jreads.append((p, c, 0, ("XS", s)))
# a read spanning TWO junctions -> (102000,103000) and (104000,104500), XS +
jreads.append((101_960, "40M1000N1000M500N40M", 0, ("XS", "+")))
# minimap2 ts: forward read (flag 0) -> genomic == ts; reverse read (flag 16) -> flipped
p, c = spliced(110_000, 110_500); jreads.append((p, c, 0, ("ts", "+")))    # genomic +
p, c = spliced(111_000, 111_400); jreads.append((p, c, 16, ("ts", "+")))   # genomic -
# no tag -> emitted with strand '.'
p, c = spliced(112_000, 112_300); jreads.append((p, c, 0, None))

import re


def read_mlen(c):
    return sum(int(x) for x in re.findall(r"(\d+)[M=XIS]", c))


def genomic_strand(flag, tag):
    if tag is None:
        return "."
    kind, s = tag
    if kind == "XS":
        return s
    return ("-" if s == "+" else "+") if (flag & 16) else s   # ts: flip if reverse


jsam = os.path.join(OUT, "_junc.sam")
with open(jsam, "w") as f:
    f.write("@HD\tVN:1.6\tSO:unsorted\n")
    f.write(f"@SQ\tSN:{CHROM}\tLN:{SIZE}\n")
    for i, (pos, c, flag, tag) in enumerate(jreads):
        m = read_mlen(c)
        t = f"\t{tag[0]}:A:{tag[1]}" if tag else ""
        f.write(f"j{i}\t{flag}\t{CHROM}\t{pos+1}\t60\t{c}\t*\t0\t0\t{'A'*m}\t{'I'*m}{t}\n")

jbam = os.path.join(OUT, "synth_junc.bam")
subprocess.run(f"samtools sort -o {jbam} {jsam}", shell=True, check=True)
subprocess.run(f"samtools index {jbam}", shell=True, check=True)
os.remove(jsam)


# Authoritative oracle: every N junction, strand resolved as bamJunctionQuery does.
def all_junctions(reads):
    out = []
    for pos, c, flag, tag in reads:
        st = genomic_strand(flag, tag)
        ref = pos
        for ln, op in re.findall(r"(\d+)([MIDNSHP=X])", c):
            ln = int(ln)
            if op == "N":
                out.append((ref, ref + ln, st))
            if op in ("M", "=", "X", "D", "N"):
                ref += ln
    return out


occ = all_junctions(jreads)
jqueries = [(CHROM, 101_000, 106_000), (CHROM, 102_000, 102_001),
            (CHROM, 107_000, 109_000), (CHROM, 109_000, 113_000),  # ts + / ts - / no-tag '.'
            (CHROM, 120_000, 121_000), ("chrX", 0, 1000)]
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


# --- stranded fixture: reads with controlled FLAGs to exercise readTxnStrand ---
# (pos0, flag): read1/read2 x forward/reverse x SE, non-overlapping 100M reads.
sreads = [(200_000, 83),   # PAIRED|PROPER|READ1|REVERSE
          (200_500, 67),   # PAIRED|PROPER|READ1 (forward)
          (201_000, 147),  # PAIRED|PROPER|READ2|REVERSE
          (201_500, 131),  # PAIRED|PROPER|READ2 (forward)
          (202_000, 0)]    # single-end, forward
RLEN = 100

ssam = os.path.join(OUT, "_str.sam")
with open(ssam, "w") as f:
    f.write("@HD\tVN:1.6\tSO:unsorted\n")
    f.write(f"@SQ\tSN:{CHROM}\tLN:{SIZE}\n")
    for i, (pos, flag) in enumerate(sreads):
        f.write(f"s{i}\t{flag}\t{CHROM}\t{pos+1}\t60\t{RLEN}M\t=\t{pos+1}\t0\t{'A'*RLEN}\t{'I'*RLEN}\n")

sbam = os.path.join(OUT, "synth_str.bam")
subprocess.run(f"samtools sort -o {sbam} {ssam}", shell=True, check=True)
subprocess.run(f"samtools index {sbam}", shell=True, check=True)
os.remove(ssam)


def txn(flag, libmode):
    rev = bool(flag & 16)
    read2 = bool(flag & 1) and bool(flag & 128)
    if libmode == "fr":
        fwd = rev if read2 else (not rev)
    else:  # rf (dUTP)
        fwd = (not rev) if read2 else rev
    return "+" if fwd else "-"


# query [199000,203000) covers all 5 reads; expected runs per (want, libmode)
soracle = {"queries": []}
for want, lib in [("+", "rf"), ("-", "rf"), ("+", "fr"), ("-", "fr"), (".", "none")]:
    runs = []
    for pos, flag in sreads:
        if want == "." or txn(flag, lib) == want:
            runs.append([pos, pos + RLEN, 1])
    soracle["queries"].append({"chrom": CHROM, "start": 199_000, "end": 203_000,
                               "want": want, "libmode": lib, "runs": sorted(runs)})
json.dump(soracle, open(os.path.join(OUT, "oracle_strand.json"), "w"), indent=1)
print(f"wrote synth_str.bam (+.bai), oracle_strand.json; {len(sreads)} stranded reads")
