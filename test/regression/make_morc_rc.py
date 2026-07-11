#!/usr/bin/env python3
"""Regenerate the reverse-strand MORC3 fixture for the `morc_rc_rnaseq` case.

Reverse-complements the chr21 MORC3 window so the (forward '+') MORC3 gene lands
on the '-' strand of a standalone `morc_rc` sequence, and maps its RNA-seq
evidence to match: intron junctions and stranded coverage get reversed
coordinates and a flipped strand (plus<->minus). Running geneid on this then
predicts the same 17-CDS/2-UTR gene on '-', exercising the reverse-strand
coordinate paths (sr[] genomic<->RSequence flip, evidence on '-') end to end.

Needs pybigtools + chr21.fa (run.sh derives chr21.fa from the tracked .gz). The
outputs (samples/morc_rc.fa, .introns.gff, .plus.bw, .minus.bw) are committed."""
import os
import pybigtools

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
W0, W1 = 36_315_000, 36_380_000                 # genomic window, 1-based inclusive
Lw = W1 - W0 + 1
CH = "morc_rc"
COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
FLIP = {"+": "-", "-": "+"}


def rc_pos(p):                                  # genomic 1-based -> rc-fasta 1-based
    return Lw - (p - W0 + 1) + 1


def s(name):
    return os.path.join(ROOT, "samples", name)


chrom = "".join(l.strip() for l in open(s("chr21.fa")) if not l.startswith(">")).upper()
rc = "".join(COMP.get(c, "N") for c in reversed(chrom[W0 - 1:W1]))
with open(s("morc_rc.fa"), "w") as f:
    f.write(f">{CH}\n")
    for i in range(0, len(rc), 60):
        f.write(rc[i:i + 60] + "\n")

# introns: reverse coords + flip strand, then SORT by acceptor (begin) as
# ReadExonsGFF requires.
introns = []
for ln in open(s("ENCFF001.1.MORC.introns.gff")):
    c = ln.split("\t")
    b, e, cnt, st = int(c[3]), int(c[4]), c[5], c[6]
    introns.append((rc_pos(e), rc_pos(b), cnt, FLIP[st]))   # begin=rc(end), end=rc(begin)
introns.sort()
with open(s("morc_rc.introns.gff"), "w") as o:
    for b, e, cnt, st in introns:
        o.write(f"{CH}\t.\tIntron\t{b}\t{e}\t{cnt}\t{st}\t.\t.\n")

# stranded coverage: reverse coords + swap plus/minus
plus, minus = [], []
for ln in open(s("ENCFF001.1.MORC.stranded.expression.shuffled.gff")):
    c = ln.split("\t")
    if len(c) < 7:
        continue
    b, e, v, st = int(c[3]), int(c[4]), float(c[5]), c[6]
    rows = plus if FLIP[st] == "+" else minus
    rows.append((CH, rc_pos(e) - 1, rc_pos(b), v))          # 0-based half-open, reversed
for name, rows in (("plus", plus), ("minus", minus)):
    rows.sort(key=lambda r: (r[0], r[1]))
    for a, nxt in zip(rows, rows[1:]):
        assert a[2] <= nxt[1], (a, nxt)                     # non-overlapping
    w = pybigtools.open(s(f"morc_rc.{name}.bw"), "w")
    w.write({CH: Lw}, iter(rows))
    w.close()

print(f"wrote morc_rc.fa ({Lw} bp), morc_rc.introns.gff ({len(introns)}), "
      f"morc_rc.plus.bw ({len(plus)}), morc_rc.minus.bw ({len(minus)})")
