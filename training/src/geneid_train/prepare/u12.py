"""U12-type intron models for training the minor-spliceosome profiles.

A single genome almost never has enough U12-type introns to train their splice
profiles (they are <1% of introns, and the branch point is unlabelled). This
module instead consumes a cross-species set of *known* U12 introns — the IAOD
(Intron Annotation and Orthology Database, Moyer/Larue/Roy/Padgett 2020, the
modern successor to U12DB) — as full intron sequences, and trains the U12 donor
and acceptor profiles from the pooled set.

Two things differ from the U2 site profiles (``stats.sites``):

- **Subtype split.** geneid keeps separate U12 profiles for the GT-AG and AT-AC
  subtypes (donor + acceptor + branch trio each), so introns are partitioned by
  their terminal dinucleotides.
- **No invariant-dinucleotide clamp.** U2 donors/acceptors mask the invariant
  GT/AG (``mask_invariant_dinuc``); U12 terminal dinucleotides are *more
  degenerate* (real GT-AG U12 introns also appear as AT-AC and rarer variants),
  so the U12 profiles are left unclamped to score the true, looser distribution.

The IAOD fastas are intron-only (donor at the 5' end, acceptor at the 3' end);
the U12 5'SS (``GTATCCTT``/``ATATCCTT``), branch and PPT motifs are all intronic,
so intron-only windows are sufficient for these profiles.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from ..stats.sites import Matrix, log_ratio, position_matrix, submatrix


@dataclass
class U12Intron:
    """One U12-type intron: full intronic sequence plus IAOD provenance."""

    seq: str
    species: str = ""
    chrom: str = ""
    strand: str = "+"
    start: int = 0
    end: int = 0

    @property
    def subtype(self) -> str:
        """``gtag`` / ``atac`` by terminal dinucleotides, else ``other``."""
        d, a = self.seq[:2], self.seq[-2:]
        if d == "GT" and a == "AG":
            return "gtag"
        if d == "AT" and a == "AC":
            return "atac"
        return "other"


def parse_iaod_fasta(path) -> list[U12Intron]:
    """Parse an IAOD ``*_U12.fasta`` file. Headers look like
    ``>Homo sapiens|chrom|strand|start|end|length|...|intron#``."""
    out: list[U12Intron] = []
    header: str | None = None
    seq_parts: list[str] = []

    def flush() -> None:
        if header is None:
            return
        f = header[1:].split("|")
        seq = "".join(seq_parts).upper()
        out.append(
            U12Intron(
                seq=seq,
                species=f[0] if len(f) > 0 else "",
                chrom=f[1] if len(f) > 1 else "",
                strand=f[2] if len(f) > 2 else "+",
                start=int(f[3]) if len(f) > 3 and f[3].isdigit() else 0,
                end=int(f[4]) if len(f) > 4 and f[4].isdigit() else 0,
            )
        )

    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush()
            header, seq_parts = line, []
        elif line:
            seq_parts.append(line)
    flush()
    return out


def load_u12_introns(paths: Iterable[str]) -> list[U12Intron]:
    """Load and pool U12 introns from several IAOD fastas (e.g. many genomes)."""
    introns: list[U12Intron] = []
    for p in paths:
        introns.extend(parse_iaod_fasta(p))
    return introns


def by_subtype(introns: Iterable[U12Intron]) -> dict[str, list[U12Intron]]:
    groups: dict[str, list[U12Intron]] = {"gtag": [], "atac": [], "other": []}
    for it in introns:
        groups[it.subtype].append(it)
    return groups


def donor_windows(introns: Iterable[U12Intron], length: int) -> list[str]:
    """The first ``length`` intronic bases (the U12 5' splice site region)."""
    return [it.seq[:length] for it in introns if len(it.seq) >= length]


def acceptor_windows(introns: Iterable[U12Intron], length: int) -> list[str]:
    """The last ``length`` intronic bases (PPT + acceptor region)."""
    return [it.seq[-length:] for it in introns if len(it.seq) >= length]


def consensus(seqs: Sequence[str], n: int, *, from_end: bool = False) -> str:
    """Most-frequent base per column over the first (or last) ``n`` positions —
    used to confirm the trained set really is U12 (donor -> GTATCCTT / ATATCCTT)."""
    cols: list[Counter] = [Counter() for _ in range(n)]
    for s in seqs:
        w = s[-n:] if from_end else s[:n]
        for i, ch in enumerate(w):
            cols[i][ch] += 1
    return "".join(c.most_common(1)[0][0] if c else "N" for c in cols)


def train_u12_profile(
    seqs: Sequence[str],
    background: Matrix,
    *,
    order: int,
    start: int,
    end: int,
) -> Matrix:
    """Train an *unclamped* U12 site profile: per-position order-k log-ratio of the
    site vs ``background`` over window ``[start, end]``. Unlike the U2 profiles no
    invariant-dinucleotide masking is applied — the terminal dinucleotides are
    part of the (degenerate) U12 signal, not a fixed anchor.
    """
    site = position_matrix(seqs, order=order)
    return submatrix(log_ratio(site, background), start, end)
