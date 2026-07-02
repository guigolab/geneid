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
