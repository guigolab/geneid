"""Gene-model length statistics for the parameter file's GenAmic rules.

Reproduces the intron/intergenic distance bounds the legacy ``WriteStatsFile``
injects into the gene model: the intron range gates intron-spanning intragenic
connections, the intergenic range gates gene-to-gene connections.
"""

from __future__ import annotations

import math
from collections.abc import Sequence


def _mean_sd(values: Sequence[float]) -> tuple[float, float]:
    """Mean and *population* standard deviation (matches geneidCEGMA::average)."""
    n = len(values)
    mean = sum(values) / n
    var = sum((v - mean) ** 2 for v in values) / n
    return mean, math.sqrt(var)


def intron_range(
    intron_lengths: Sequence[int],
    *,
    short_cap: int = 40,
    long_cap: int = 100_000,
) -> tuple[float, float]:
    """Return ``(min_intron, max_intron)`` for the gene model.

    ``min = 0.75 * shortest_intron`` capped at ``short_cap`` (40); ``max =
    mean + 3*sd`` capped at ``long_cap`` (100000) — the legacy WriteStatsFile
    heuristic.
    """
    shortest = min(intron_lengths)
    lo = shortest * 0.75
    if lo > short_cap:
        lo = float(short_cap)
    mean, sd = _mean_sd(list(intron_lengths))
    hi = mean + 3 * sd
    if hi > long_cap:
        hi = float(long_cap)
    return lo, hi


def format_range(lo: float, hi: float | str) -> str:
    """Render a ``min:max`` distance token (``Infinity`` passes through)."""
    def one(x: float | str) -> str:
        if isinstance(x, str):
            return x
        return str(int(x)) if x == int(x) else f"{x:.10g}"

    return f"{one(lo)}:{one(hi)}"
