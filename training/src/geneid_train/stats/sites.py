"""Position weight arrays / order-k Markov site matrices.

Replaces the legacy AWK chain (Getkmatrix / logratio_kmatrix / information /
submatrix). A site profile is built in stages, each of which is validated against
the reference xgXerMont run:

1. ``position_matrix`` — per-position order-k conditional probabilities P(next|prev)
   from a set of aligned site windows  (reference: ``*.canonical.<site>.di-matrix``)
2. ``log_ratio`` — natural-log ratio of the site matrix vs a background model
   (reference: ``*-log.di-matrix``)
3. information-weighting + boundary selection -> the profile rows written into the
   parameter file  (reference: ``*-log-info.di-matrix`` == the ``*_profile`` block)

Stage 3 is not yet implemented. Matrices are dicts keyed by ``(position, oligo)``
with 1-based positions, matching the geneid ``.param`` profile layout.

geneid profile "order" k means oligos of length k+1 (order 1 = dinucleotides).
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from itertools import product
from pathlib import Path

Matrix = dict[tuple[int, str], float]

_ACGT = frozenset("ACGT")

MASK = -9999.0  # geneid sentinel for a forbidden PWM cell

# Dirichlet pseudocount per k+1-tuple, matching legacy Getkmatrix.awk. Every
# oligo gets +PCOUNT and every prefix +4*PCOUNT, so unobserved oligos still get a
# defined (non-zero) probability and every position has a complete matrix.
PCOUNT = 0.25


def position_matrix(seqs: Iterable[str], order: int = 1, pcount: float = PCOUNT) -> Matrix:
    """Per-position order-k conditional probabilities from aligned site windows.

    For each position ``i`` (1-based) the value for oligo ``prefix+base`` is the
    smoothed conditional P(base | prefix):
    ``(pcount + count(oligo)) / (4*pcount + count(prefix))``. Positions run
    1..L-order and the matrix is complete (all 4**(order+1) oligos per position).

    A whole window containing any non-ACGT character is skipped (matching the
    legacy estimator), not just the offending column.
    """
    seqs = [s.upper() for s in seqs]
    if not seqs:
        return {}
    length = min(len(s) for s in seqs)
    counts: dict[tuple[int, str], int] = defaultdict(int)
    prefix_counts: dict[tuple[int, str], int] = defaultdict(int)
    for s in seqs:
        window = s[:length]
        if set(window) - _ACGT:
            continue
        for i in range(1, length - order + 1):
            counts[(i, window[i - 1 : i - 1 + order + 1])] += 1
            prefix_counts[(i, window[i - 1 : i - 1 + order])] += 1

    out: Matrix = {}
    denom_pseudo = 4 * pcount
    for pos in range(1, length - order + 1):
        for prefix in ("".join(p) for p in product("ACGT", repeat=order)):
            denom = denom_pseudo + prefix_counts.get((pos, prefix), 0)
            for base in "ACGT":
                oligo = prefix + base
                num = pcount + counts.get((pos, oligo), 0)
                out[(pos, oligo)] = num / denom if denom else 0.0
    return out


def log_ratio_zero_order(site: Matrix, background: Matrix) -> Matrix:
    """Log-ratio for order-0 (single-nucleotide) profiles, matching legacy
    ``logratio_zero_order.awk``.

    Order-0 site matrices are built with no pseudocounts (see ``position_matrix``
    with ``pcount=0``), so a base that never occurs at a position has an exact raw
    frequency of 0 -- and a base that occurs at every site has frequency exactly
    1. Those two cases are special-cased to -9999 / 0 rather than computing
    ``log(0/bg)`` or ``log(1/bg)``; this is what naturally masks an invariant
    position (e.g. the A of ATG) without a separate explicit masking step.
    """
    out: Matrix = {}
    for key, fg in site.items():
        bg = background.get(key) or background.get((1, key[1]))
        if fg == 0:
            out[key] = MASK
        elif fg == 1:
            out[key] = 0.0
        elif bg:
            out[key] = math.log(fg / bg)
    return out


def log_ratio(site: Matrix, background: Matrix) -> Matrix:
    """Natural-log ratio of a site matrix against a background matrix, per cell.

    Background may be keyed per-position (matching keys) or position-independent
    (keyed by position 1); the latter is broadcast across positions.
    """
    pos_independent = all(p == 1 for p, _ in background) and background
    out: Matrix = {}
    for (pos, oligo), fg in site.items():
        bg = background.get((pos, oligo))
        if bg is None and pos_independent:
            bg = background.get((1, oligo))
        if bg:
            out[(pos, oligo)] = math.log(fg / bg)
    return out


def submatrix(matrix: Matrix, start: int, end: int) -> Matrix:
    """Restrict to positions [start, end] and renumber them to 1..(end-start+1),
    matching the legacy ``submatrix.awk``."""
    return {
        (pos - start + 1, oligo): v for (pos, oligo), v in matrix.items() if start <= pos <= end
    }


def mask_invariant_dinuc(matrix: Matrix, st: int, nd: int, rd: int, anchor: str) -> Matrix:
    """Apply the invariant splice-dinucleotide constraint to an order-1 profile
    (replaces preparedimatrix{donor,acceptor}4parameter.awk).

    ``anchor`` is the invariant dinucleotide ("GT" for donors, "AG" for acceptors)
    sitting at profile positions ``nd``/``nd+1``. Across the three order-1 columns
    overlapping it:

    - ``st`` (=nd-1): cell kept as ``0`` if its oligo ends in anchor[0], else masked
    - ``nd``:         cell kept as ``0`` if its oligo == anchor, else masked
    - ``rd`` (=nd+1): cell keeps its log value if its oligo starts with anchor[1],
      else masked

    Every other cell keeps its log-ratio value.
    """
    a, b = anchor[0], anchor[1]
    out: Matrix = {}
    for (pos, oligo), v in matrix.items():
        if pos == st:
            out[(pos, oligo)] = 0.0 if oligo.endswith(a) else MASK
        elif pos == nd:
            out[(pos, oligo)] = 0.0 if oligo == anchor else MASK
        elif pos == rd:
            out[(pos, oligo)] = v if oligo.startswith(b) else MASK
        else:
            out[(pos, oligo)] = v
    return out


# --- boundary selection (replaces frequency.awk / information.awk / BitScoreGraph) ---

# Per-site defaults from the legacy driver (geneidTRAINer1_3.pl): the anchor
# ``offset`` seed, the information-content threshold (bits) above which a position
# joins the profile window, and the raw-position clip range the info graph is
# restricted to before selection.
SITE_DEFAULTS: dict[str, dict] = {
    "donor": {"offset": 30, "info_thresh": 0.15, "clip": (25, 38)},
    "acceptor": {"offset": 30, "info_thresh": 0.04, "clip": (2, 33)},
    "start": {"offset": 30, "info_thresh": 0.15, "clip": (25, 37)},
    "branch": {"offset": 32, "info_thresh": 0.30, "clip": (28, 41)},
}


def frequency(seqs: Iterable[str]) -> dict[tuple[int, str], float]:
    """Per-position single-nucleotide frequencies (replaces ``frequency.awk`` k=1).

    Returns ``{(pos, base): count/total}`` with 1-based positions. Unlike
    :func:`position_matrix`, a non-ACGT character masks only its own column at
    that position (not the whole window), and the denominator is the per-position
    count of valid ACGT observations. No pseudocounts.
    """
    seqs = [s.upper() for s in seqs]
    if not seqs:
        return {}
    length = min(len(s) for s in seqs)
    counts: dict[tuple[int, str], int] = defaultdict(int)
    totals: dict[int, int] = defaultdict(int)
    for s in seqs:
        for i in range(1, length + 1):
            base = s[i - 1]
            if base in _ACGT:
                counts[(i, base)] += 1
                totals[i] += 1
    out: dict[tuple[int, str], float] = {}
    for pos in range(1, length + 1):
        tot = totals.get(pos, 0)
        for base in "ACGT":
            out[(pos, base)] = counts.get((pos, base), 0) / tot if tot else 0.0
    return out


def info_content(
    site: dict[tuple[int, str], float], background: dict[tuple[int, str], float]
) -> dict[int, float]:
    """Per-position information content in bits (replaces ``information.awk`` k=1).

    ``info[pos] = sum_base p*log2(p/bg)`` over the four nucleotides, where ``p`` is
    the site frequency and ``bg`` the background frequency; a term is dropped when
    either ``p`` or ``bg`` is zero. This is the relative entropy of the site vs
    background distribution at each position — the quantity the boundary selector
    thresholds on.
    """
    out: dict[int, float] = defaultdict(float)
    positions = {pos for pos, _ in site}
    for pos in positions:
        for base in "ACGT":
            p = site.get((pos, base), 0.0)
            bg = background.get((pos, base), 0.0)
            if p > 0 and bg > 0:
                out[pos] += p * math.log(p / bg) / math.log(2)
    return dict(out)


@dataclass(frozen=True)
class SiteWindow:
    """A selected profile window plus the invariant-dinucleotide anchor columns.

    ``start``/``end`` are raw (pre-submatrix) positions passed to :func:`submatrix`;
    ``length`` is the resulting profile length. ``st``/``nd``/``rd`` are the
    profile-local (1-based, post-renumber) columns overlapping the invariant
    splice/start motif, ready for :func:`mask_invariant_dinuc`; they are ``None``
    for order-0 profiles, which carry no dinucleotide anchor.
    """

    start: int
    end: int
    offset: int
    length: int
    st: int | None
    nd: int | None
    rd: int | None


def select_window(
    info: dict[int, float],
    *,
    site: str,
    order: int,
    offset: int | None = None,
    info_thresh: float | None = None,
    clip: tuple[int, int] | None = None,
) -> SiteWindow:
    """Select the profile window and anchor columns from per-position info content.

    Reproduces the legacy ``BitScoreGraph`` + post-processing in the Perl driver:
    seed the candidate set with the anchor's flanking positions (``offset-1`` and
    ``offset+1``), add every clipped position whose info content exceeds
    ``info_thresh``, then take ``start = min`` (floored at 1) and ``end = max``.
    The offset is then renumbered into profile-local coordinates and the
    site-specific anchor columns are derived from it.

    ``site`` is one of ``donor``/``acceptor``/``start``/``branch``; the remaining
    keyword args default to :data:`SITE_DEFAULTS` for that site.
    """
    d = SITE_DEFAULTS.get(site, {})
    if offset is None:
        offset = d["offset"]
    if info_thresh is None:
        info_thresh = d["info_thresh"]
    if clip is None:
        clip = d["clip"]
    lo, hi = clip

    candidates = {offset - 1, offset + 1}
    for pos, bits in info.items():
        if lo <= pos <= hi and bits > info_thresh:
            candidates.add(pos)
    start = max(1, min(candidates))
    end = max(candidates)

    # renumber the anchor offset into the submatrix's 1-based coordinates
    off = offset - order
    end = end - order
    off = off - start + 1

    st = nd = rd = None
    if site == "donor" and order >= 1:
        st, nd, rd = off + 2, off + 3, off + 4
    elif site == "acceptor" and order >= 1:
        st, nd, rd = off - 1, off, off + 1
    elif site == "start" and order >= 2:
        st, nd, rd = off - 2, off - 1, off

    return SiteWindow(start, end, off, end - start + 1, st, nd, rd)


def _fmt(value: float) -> str:
    """Format a matrix value as geneid does: 6 significant figures, integers bare."""
    if value == int(value):
        return str(int(value))
    return f"{value:g}"


def profile_header(window: SiteWindow, order: int, cutoff: float = -7.0) -> list[str]:
    """The geneid profile header ``len offset cutoff order [a b]`` for a window.

    Order>=1 profiles carry the trailing ``0 1`` pair (the dinucleotide-context
    flags geneid expects); order-0 profiles omit it.
    """
    fields = [window.length, window.offset, cutoff, order]
    if order >= 1:
        fields += [0, 1]
    return [_fmt(float(f)) for f in fields]


def format_profile(matrix: Matrix, header: list[str]) -> list[str]:
    """Render a site profile as geneid param lines: a header line, the standard
    comment, then ``pos oligo value`` rows sorted by position then oligo."""
    lines = [" ".join(header), "# Transition probabilities at every position"]
    for (pos, oligo), v in sorted(matrix.items()):
        lines.append(f"{pos} {oligo} {_fmt(v)}")
    return lines


def read_matrix(path: str | Path) -> Matrix:
    """Read a geneid ``.di-matrix`` / profile-style file.

    Accepts either column order: ``<pos> <oligo> <value>`` (the profile / order>=1
    ``.di-matrix`` convention) or ``<oligo> <pos> <value>`` (the legacy order-0
    ``*.order-0-matrix`` files use base-first columns)."""
    out: Matrix = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            if parts[0].isdigit():
                pos, oligo = int(parts[0]), parts[1]
            elif parts[1].isdigit():
                oligo, pos = parts[0], int(parts[1])
            else:
                continue
            out[(pos, oligo)] = float(parts[2])
    return out
