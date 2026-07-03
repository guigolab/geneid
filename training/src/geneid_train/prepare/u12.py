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


def markov_background(seqs: Iterable[str], order: int) -> Matrix:
    """A position-independent order-k background: the pooled P(base | prefix) over
    all ``seqs``, keyed at position 1 so :func:`log_ratio` broadcasts it across
    every profile position. This is the null the U12 profiles are scored against
    (generic intronic composition), analogous to the genome background used for
    the U2 profiles."""
    from collections import defaultdict

    counts: dict[str, int] = defaultdict(int)
    prefix_counts: dict[str, int] = defaultdict(int)
    for s in seqs:
        s = s.upper()
        for i in range(len(s) - order):
            oligo = s[i : i + order + 1]
            if set(oligo) - set("ACGT"):
                continue
            counts[oligo] += 1
            prefix_counts[oligo[:order]] += 1
    from itertools import product

    out: Matrix = {}
    for pre in ("".join(p) for p in product("ACGT", repeat=order)):
        denom = prefix_counts.get(pre, 0)
        for base in "ACGT":
            out[(1, pre + base)] = counts.get(pre + base, 0) / denom if denom else 0.0
    return out


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


# --- branch-point location ---------------------------------------------------
#
# The IAOD fastas don't mark the branch adenosine, and it sits at a variable
# distance from the 3' acceptor, so it can't be read off a fixed position. But
# the U12 branch motif is strongly conserved, so we seed with a cross-taxa U12
# branch PWM (e.g. geneid's U12_Branch_point_profile), score every candidate
# window in the acceptor-upstream region, and take the best — a bootstrap that
# both locates the branch and yields aligned windows to retrain the profile.


@dataclass
class BranchHit:
    """A located branch point: score, distance of the anchor base from the
    acceptor (last intron base), and the aligned window used for retraining."""

    score: float
    distance: int
    window: str


def score_window(window: str, pwm: Matrix, order: int) -> float:
    """Sum of the order-k log-ratio values of ``window`` against a profile PWM
    (positions 1..len keyed by the ``order+1``-mer starting at each position).
    Windows with a non-ACGT base score ``-inf`` so they are never chosen."""
    length = max(p for p, _ in pwm)
    total = 0.0
    for pos in range(1, length + 1):
        oligo = window[pos - 1 : pos - 1 + order + 1]
        if len(oligo) < order + 1 or set(oligo) - set("ACGT"):
            return float("-inf")
        total += pwm.get((pos, oligo), 0.0)
    return total


def locate_branch(
    seq: str,
    pwm: Matrix,
    *,
    order: int,
    offset: int,
    acc_context: int = 50,
    min_dist: int = 7,
) -> BranchHit | None:
    """Find the best-scoring branch window in the acceptor-upstream region of an
    intron. ``offset`` is the profile position of the branch anchor base; the
    returned ``distance`` is that base's distance from the last intron base.
    Scans windows whose anchor lies in ``[min_dist, acc_context]`` from the
    acceptor. Returns ``None`` if the intron is too short."""
    length = max(p for p, _ in pwm)
    win_len = length + order  # bases needed for `length` order-k positions
    n = len(seq)
    best: BranchHit | None = None
    # window start s (0-based); anchor base is at s + (offset-1)
    lo = n - acc_context
    hi = n - min_dist - (length - offset)
    for s in range(max(0, lo), min(hi, n - win_len) + 1):
        window = seq[s : s + win_len]
        sc = score_window(window, pwm, order)
        if sc == float("-inf"):
            continue
        anchor = s + (offset - 1)
        dist = n - anchor
        if best is None or sc > best.score:
            best = BranchHit(sc, dist, window)
    return best


def locate_branches(
    introns: Iterable[U12Intron], pwm: Matrix, *, order: int, offset: int,
    acc_context: int = 50, min_dist: int = 7,
) -> list[BranchHit]:
    """Locate the branch point in each intron; skips introns with no valid hit."""
    hits = []
    for it in introns:
        h = locate_branch(
            it.seq, pwm, order=order, offset=offset,
            acc_context=acc_context, min_dist=min_dist,
        )
        if h is not None:
            hits.append(h)
    return hits
