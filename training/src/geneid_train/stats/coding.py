"""Coding-potential Markov models (replaces geneidCEGMA.pm + deriveCodingPotential).

geneid scores an ORF with two log-ratio matrices comparing a codon-phase-specific
model of coding sequence against a phase-agnostic model of introns:

- ``Markov_Initial_probability_matrix``  — order-4 *frequencies* of 5-mers, one
  distribution per codon frame (used to seed the first bases of an exon)
- ``Markov_Transition_probability_matrix`` — order-5 Markov transitions
  P(base | preceding 5 bases), one per codon frame (used to extend)

Each is ``log(P_coding(oligo|frame) / P_intron(oligo|frame=0))`` (natural log).
The coding model cycles through frames 0/1/2 as it walks a CDS (which starts
in-frame at the ATG); the intron model uses a single frame 0.

Faithfully reproduces the legacy estimator, including its quirks: the counting
loop stops one base early (``i < len - order - 1``); a non-ACGT position is
skipped *without* advancing the frame; and the FREQ per-frame denominator is
taken one step after the count (a harmless off-by-one that averages out).
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from itertools import product

_ACGT = frozenset("ACGT")

# per-frame pseudocounts from the legacy driver
CODING_PSEUDO = 0.25
INTRON_PSEUDO = 10.0

# Model = {(oligo, frame): probability_or_logratio}, oligo = (order+1)-mer.
Model = dict[tuple[str, int], float]


def _kmers(order: int) -> list[str]:
    return ["".join(p) for p in product("ACGT", repeat=order)]


def _counts(seqs: Iterable[str], order: int, pseudo: float, nframes: int) -> dict:
    """Per-frame (prefix, nt) counts, initialised to ``pseudo``. ``nframes`` is 1
    (phase-agnostic, intron) or 3 (phase-specific, coding). Returns
    ``(table, totals)`` where table[(prefix, nt, fr)] is a count and totals[fr]
    the running per-frame position count used by the FREQ normaliser."""
    frame = nframes - 1  # legacy "frame" arg: 0 or 2
    table: dict[tuple[str, str, int], float] = {}
    for pre in _kmers(order):
        for nt in "ACGT":
            for fr in range(frame + 1):
                table[(pre, nt, fr)] = pseudo
    totals = [0, 0, 0]
    for seq in seqs:
        s = seq.upper()
        fr = 0
        for i in range(0, len(s) - order - 1):
            pre = s[i : i + order]
            nt = s[i + order]
            if _ACGT.issuperset(pre) and nt in _ACGT:
                table[(pre, nt, fr)] += 1
                fr = 0 if frame == 0 else (fr + 1) % 3
                totals[fr] += 1
        # (frame is not advanced on skipped positions — matches the legacy `next`)
    return table, totals


def initial_model(
    seqs: Iterable[str], order: int, pseudo: float, nframes: int
) -> Model:
    """FREQ model: P(oligo | frame) = count / (per-frame total position count)."""
    seqs = list(seqs)
    table, totals = _counts(seqs, order, pseudo, nframes)
    out: Model = {}
    for (pre, nt, fr), c in table.items():
        denom = totals[fr]
        out[(pre + nt, fr)] = c / denom if denom else 0.0
    return out


def transition_model(
    seqs: Iterable[str], order: int, pseudo: float, nframes: int
) -> Model:
    """MM model: P(nt | prefix, frame) = count / sum_nt count over the prefix."""
    seqs = list(seqs)
    table, _ = _counts(seqs, order, pseudo, nframes)
    frame = nframes - 1
    out: Model = {}
    prefixes = _kmers(order)
    for pre in prefixes:
        for fr in range(frame + 1):
            denom = sum(table[(pre, nt, fr)] for nt in "ACGT")
            for nt in "ACGT":
                out[(pre + nt, fr)] = table[(pre, nt, fr)] / denom if denom else 0.0
    return out


def coding_log_ratio(coding: Model, intron: Model) -> Model:
    """log(P_coding(oligo|frame) / P_intron(oligo|frame 0)) per coding-model cell."""
    out: Model = {}
    for (oligo, fr), pc in coding.items():
        pi = intron.get((oligo, 0))
        if pc > 0 and pi:
            out[(oligo, fr)] = math.log(pc / pi)
    return out


def choose_orders(coding_bases: int, noncoding_bases: int) -> tuple[int, int]:
    """Return ``(initial_order, transition_order)`` for the training-set size.

    The reference driver always trains order-5 transitions / order-4 initials
    when it derives a coding potential at all; the size test only gates *whether*
    to derive one. We mirror that (5/4) whenever there is enough data, and
    fall back to a smaller model on sparse sets to avoid overfitting.
    """
    enough = (
        (coding_bases > 400_000 and noncoding_bases > 100_000)
        or (coding_bases > 375_000 and noncoding_bases > 150_000)
        or (noncoding_bases > 35_000 and coding_bases > 25 * noncoding_bases)
    )
    return (4, 5) if enough else (3, 4)


def derive_coding_potential(
    cds_seqs: Sequence[str], intron_seqs: Sequence[str], transition_order: int = 5
) -> tuple[Model, Model, int, int]:
    """Build the initial and transition coding-potential log-ratio matrices.

    Returns ``(initial_logs, transition_logs, coding_bases, noncoding_bases)``.
    The initial model is order ``transition_order - 1``; both compare a
    3-frame coding model to a 1-frame intron model.
    """
    coding_bases = sum(len(s) for s in cds_seqs)
    noncoding_bases = sum(len(s) for s in intron_seqs)
    init_order = transition_order - 1

    coding_init = initial_model(cds_seqs, init_order, CODING_PSEUDO, 3)
    intron_init = initial_model(intron_seqs, init_order, INTRON_PSEUDO, 1)
    coding_trans = transition_model(cds_seqs, transition_order, CODING_PSEUDO, 3)
    intron_trans = transition_model(intron_seqs, transition_order, INTRON_PSEUDO, 1)

    initial_logs = coding_log_ratio(coding_init, intron_init)
    transition_logs = coding_log_ratio(coding_trans, intron_trans)
    return initial_logs, transition_logs, coding_bases, noncoding_bases


def format_markov_matrix(model: Model) -> list[str]:
    """Render a log-ratio model as geneid param lines: ``oligo index frame value``,
    index = lexicographic rank of the oligo, ordered oligo then frame."""
    by_oligo: dict[str, list[int]] = {}
    for oligo, fr in model:
        by_oligo.setdefault(oligo, []).append(fr)
    lines = []
    for idx, oligo in enumerate(sorted(by_oligo)):
        for fr in sorted(by_oligo[oligo]):
            lines.append(f"{oligo} {idx} {fr} {model[(oligo, fr)]:.3f}")
    return lines
