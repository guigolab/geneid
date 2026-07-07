"""Gene-model length statistics for the parameter file's GenAmic rules.

Reproduces the intron/intergenic distance bounds the legacy ``WriteStatsFile``
injects into the gene model: the intron range gates intron-spanning intragenic
connections, the intergenic range gates gene-to-gene connections.
"""

from __future__ import annotations

import math
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
    max_intron: float | None = None,
) -> tuple[float, float]:
    """Return ``(min_intron, max_intron)`` for the gene model.

    ``min = 0.75 * shortest_intron`` capped at ``short_cap`` (40). ``max`` is the
    ``max_quantile`` (default p99.9) of the intron lengths, capped at ``long_cap``
    (100000). Intron lengths are strongly right-skewed, so the legacy
    ``mean + 3*sd`` heuristic clips a real long tail (e.g. ~1.8% of introns on
    xgXerMont) — and geneid cannot span an intron longer than this max, so a too-low
    value fragments long-intron genes into separate models. A high percentile spans
    essentially all real introns while staying under the safety cap.

    ``max_intron``, when given, overrides the percentile/``long_cap`` computation
    and is used directly as the max. That is the intended setting alongside the
    soft intron-length penalty (``Intron_length_model``): the hard max becomes a
    generous *safety bound* (e.g. 500 kb for human, which admits all but the ~9-35
    longest introns) while the smooth penalty — not a cliff — grades the long ones.
    """
    shortest = min(intron_lengths)
    lo = shortest * 0.75
    if lo > short_cap:
        lo = float(short_cap)
    if max_intron is not None:
        return lo, float(max_intron)
    hi = _percentile(list(intron_lengths), max_quantile)
    if hi > long_cap:
        hi = float(long_cap)
    return lo, hi


def intron_length_model(intron_lengths: Sequence[int]) -> tuple[float, float]:
    """Fit a log-normal to intron lengths; return ``(mu, sigma)`` in log space.

    These are the maximum-likelihood log-normal parameters — ``mu`` = mean and
    ``sigma`` = population standard deviation of ``ln(length)``. They are emitted
    into the param's ``Intron_length_model`` section, where geneid uses them for a
    smooth, length-dependent intron score penalty: a *soft* replacement for the
    hard gene-model ``max`` gate (see ``intron_range``). Introns near the typical
    length sit near the distribution mode and are essentially unpenalised; the
    penalty grows only in the long tail. Intron lengths are strongly right-skewed,
    so a log-normal fits far better than a normal on the raw lengths.

    The paired penalty *weight* (``Intron_length_score_weight``) defaults to 0 in
    geneid, so emitting this model is backward-compatible until the weight is
    turned on (by the optimizer or by hand).
    """
    logs = [math.log(n) for n in intron_lengths if n > 0]
    n = len(logs)
    if n == 0:
        return 0.0, 0.0
    mu = sum(logs) / n
    var = sum((x - mu) ** 2 for x in logs) / n
    return mu, math.sqrt(var)


def format_range(lo: float, hi: float | str) -> str:
    """Render a ``min:max`` distance token (``Infinity`` passes through)."""
    def one(x: float | str) -> str:
        if isinstance(x, str):
            return x
        return str(int(x)) if x == int(x) else f"{x:.10g}"

    return f"{one(lo)}:{one(hi)}"
