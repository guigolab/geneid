"""Genome-wide background sequence model (replaces the legacy getBackground).

The site profiles are log-ratios of a site model against a background model of
generic genomic composition. The legacy trainer estimates that background from
``numseqs`` (100000) random ``k``-mers (k=62) drawn from the genome, then builds
per-position single-nucleotide frequencies and an order-1 dinucleotide matrix
from them. Reproduced here; the sampling is seeded so a given genome yields a
reproducible background (the legacy used an unseeded RNG, so exact byte-repro of
its background was never possible).
"""

from __future__ import annotations

import random
from collections.abc import Iterable

from .sites import Matrix, frequency, position_matrix

DEFAULT_K = 62
DEFAULT_N = 100_000


def sample_kmers(
    genome_seqs: Iterable[str], k: int = DEFAULT_K, n: int = DEFAULT_N, seed: int = 0
) -> list[str]:
    """Draw ``n`` random ``k``-mers from the genome, rejecting any with a non-ACGT
    base. Sequences are concatenated (matching the legacy) only up to the length
    needed to supply the samples, then sampled uniformly with a seeded RNG."""
    need = n * k + 1
    pieces: list[str] = []
    total = 0
    for s in genome_seqs:
        pieces.append(s)
        total += len(s)
        if total >= need:
            break
    seq = "".join(pieces).upper()
    span = len(seq) - k
    if span <= 0:
        raise ValueError(f"genome too short ({len(seq)} bp) for {k}-mer background")
    rng = random.Random(seed)
    out: list[str] = []
    while len(out) < n:
        r = rng.randint(0, span)
        kmer = seq[r : r + k]
        if all(c in "ACGT" for c in kmer):
            out.append(kmer)
    return out


def background_models(kmers: list[str]) -> tuple[Matrix, Matrix]:
    """Build the ``(order-0 frequency, order-1 dinucleotide)`` background models
    from a set of equal-length background sequences, using the same estimators as
    the site profiles (``frequency`` for info content, ``position_matrix`` for the
    log-ratio denominator)."""
    freq = frequency(kmers)
    dimatrix = position_matrix(kmers, order=1)
    return freq, dimatrix


def from_genome(
    genome_seqs: Iterable[str], k: int = DEFAULT_K, n: int = DEFAULT_N, seed: int = 0
) -> tuple[Matrix, Matrix]:
    """Convenience: sample a background from the genome and build both models."""
    return background_models(sample_kmers(genome_seqs, k=k, n=n, seed=seed))
