"""Splice-site and start-codon window extraction (replaces the legacy SSgff step).

For each gene model this pulls the fixed-length sequence windows the profile
estimator (``stats.sites``) consumes: a donor window per intron (5' splice), an
acceptor window per intron (3' splice), and one start-codon window per gene.
Windows are in transcription orientation with the anchor base at a fixed
position, matching the geometry of the reference ``*.canonical.*.tbl`` files:

    position 31 (1-based) is the anchor base:
      donor    -> last coding base of the upstream exon   (intron GT at 32-33)
      acceptor -> first coding base of the downstream exon (intron AG at 29-30)
      start    -> the A of the ATG start codon             (ATG at 31-33)

Each window is ``BEFORE`` (30) bases upstream + the anchor + ``AFTER`` (29) bases
downstream = 60 bp. Sites too close to a sequence end to yield a full window are
skipped (the estimator requires equal-length windows).
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

from ..core.seq import revcomp
from .base import GeneModel

BEFORE = 30  # bases upstream of the anchor (transcription orientation)
AFTER = 29  # bases downstream of the anchor
WINDOW = BEFORE + 1 + AFTER  # 60


def _window(chrom: str, anchor: int, strand: str) -> str | None:
    """Extract the 60 bp transcription-oriented window with ``anchor`` (1-based
    genomic) at profile position 31, or ``None`` if it would run off the end."""
    if strand == "+":
        lo, hi = anchor - 1 - BEFORE, anchor + AFTER
    else:
        lo, hi = anchor - 1 - AFTER, anchor + BEFORE
    if lo < 0 or hi > len(chrom):
        return None
    seq = chrom[lo:hi]
    if strand == "-":
        seq = revcomp(seq)
    return seq.upper()


def donor_windows(model: GeneModel, chrom: str) -> list[str]:
    """One donor (5' splice) window per intron, transcription-oriented."""
    out = []
    for s, e in model.introns():
        anchor = (s - 1) if model.strand == "+" else (e + 1)
        w = _window(chrom, anchor, model.strand)
        if w is not None:
            out.append(w)
    return out


def acceptor_windows(model: GeneModel, chrom: str) -> list[str]:
    """One acceptor (3' splice) window per intron, transcription-oriented."""
    out = []
    for s, e in model.introns():
        anchor = (e + 1) if model.strand == "+" else (s - 1)
        w = _window(chrom, anchor, model.strand)
        if w is not None:
            out.append(w)
    return out


def start_window(model: GeneModel, chrom: str) -> str | None:
    """The start-codon window (anchor = A of the ATG), or ``None`` if truncated."""
    anchor = model.start if model.strand == "+" else model.end
    return _window(chrom, anchor, model.strand)


def is_canonical_donor(seq: str) -> bool:
    return seq[31:33] == "GT"  # profile positions 32-33


def is_canonical_acceptor(seq: str) -> bool:
    return seq[28:30] == "AG"  # profile positions 29-30


def is_canonical_start(seq: str) -> bool:
    return seq[30:33] == "ATG"  # profile positions 31-33


@dataclass
class SiteSets:
    """Canonical vs non-canonical site windows, split for separate estimation."""

    donor: list[str]
    acceptor: list[str]
    start: list[str]
    noncanonical_donor: list[str]
    noncanonical_acceptor: list[str]
    noncanonical_start: list[str]


def collect_sites(models: Iterable[GeneModel], genome: Mapping[str, str]) -> SiteSets:
    """Extract every donor/acceptor/start window across ``models`` and partition
    them into canonical (GT-AG / ATG) and non-canonical bags."""
    s = SiteSets([], [], [], [], [], [])
    for m in models:
        chrom = genome[m.seqid]
        for w in donor_windows(m, chrom):
            (s.donor if is_canonical_donor(w) else s.noncanonical_donor).append(w)
        for w in acceptor_windows(m, chrom):
            (s.acceptor if is_canonical_acceptor(w) else s.noncanonical_acceptor).append(w)
        w = start_window(m, chrom)
        if w is not None:
            (s.start if is_canonical_start(w) else s.noncanonical_start).append(w)
    return s
