"""Gene-model length statistics for the parameter file's GenAmic rules.

Reproduces the intron/intergenic distance bounds the legacy ``WriteStatsFile``
injects into the gene model: the intron range gates intron-spanning intragenic
connections, the intergenic range gates gene-to-gene connections.
"""

from __future__ import annotations

from collections.abc import Sequence


def _percentile(values: Sequence[float], q: float) -> float:
    """The ``q`` quantile (0..1) by linear interpolation between order statistics."""
    s = sorted(values)
    n = len(s)
    if n == 1:
        return float(s[0])
    idx = q * (n - 1)
    lo = int(idx)
    frac = idx - lo
    if lo + 1 < n:
        return s[lo] + frac * (s[lo + 1] - s[lo])
    return float(s[lo])


def intron_range(
    intron_lengths: Sequence[int],
    *,
    short_cap: int = 40,
    long_cap: int = 100_000,
    max_quantile: float = 0.999,
) -> tuple[float, float]:
    """Return ``(min_intron, max_intron)`` for the gene model.

    ``min = 0.75 * shortest_intron`` capped at ``short_cap`` (40). ``max`` is the
    ``max_quantile`` (default p99.9) of the intron lengths, capped at ``long_cap``
    (100000). Intron lengths are strongly right-skewed, so the legacy
    ``mean + 3*sd`` heuristic clips a real long tail (e.g. ~1.8% of introns on
    xgXerMont) — and geneid cannot span an intron longer than this max, so a too-low
    value fragments long-intron genes into separate models. A high percentile spans
    essentially all real introns while staying under the safety cap.
    """
    shortest = min(intron_lengths)
    lo = shortest * 0.75
    if lo > short_cap:
        lo = float(short_cap)
    hi = _percentile(list(intron_lengths), max_quantile)
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
