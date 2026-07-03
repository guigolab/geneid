"""U2 branch-point discovery by a self-supervised EM motif finder.

The U2 branch point was historically found by hand with MEME. This module
automates that: it fits a two-component mixture (a branch-motif PWM vs an
intronic background) to the acceptor-upstream region of a genome's *own* introns
by EM, so the result is species-appropriate — sharp for conserved (fungal)
branches, still usable when degenerate (mammalian). Nothing is human-tuned: the
motif is *seeded* only from the universal U2-snRNA-complementary consensus (the
invariant branch adenosine) and everything else is learned from the target genome.

This is deliberately *not* a port of BPP (Zhang et al. 2017), whose fixed human
offset windows and human-trained weights would misfire on short fungal introns and
non-human composition. We keep only the idea — unsupervised EM motif discovery.

The model is OOPS (one occurrence per sequence): each intron has exactly one branch
in its 3' region. The EM also yields, per intron, the most likely branch position,
whose distance distribution to the 3' splice site sets geneid's branch-distance
knobs (``opt_dist``/``min_dist``/``acc_context``) directly from data rather than by
search. Emission of the geneid ``Branch_point_profile`` + the accuracy-gated toggle
live elsewhere; this module is the discovery + distance estimation.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

from ..prepare.base import GeneModel

_ACGT = ("A", "C", "G", "T")
_ACGT_SET = frozenset(_ACGT)

# Motif geometry. Width 7 with the branch adenosine at index 5 matches the yeast
# UACUAAC / mammalian yUnAy consensus (branch A is the 6th base). AG at the 3' end
# of the intron is excluded from the search window (2 nt).
BRANCH_WIDTH = 7
BRANCH_A = 5
_AG = 2  # terminal AG excluded from the branch search window
# search the acceptor-upstream region; capped so short introns still fit (adaptive)
MAX_WINDOW = 45


@dataclass
class BranchModel:
    """A discovered branch model: the per-position motif PWM (``pwm[m][base]``),
    the intronic background composition, and the motif geometry."""

    pwm: list[dict[str, float]]
    background: dict[str, float]
    width: int
    anchor: int
    require_anchor: str | None = "A"


def branch_windows(
    models: Iterable[GeneModel],
    genome: Mapping[str, str],
    *,
    max_window: int = MAX_WINDOW,
    width: int = BRANCH_WIDTH,
) -> list[str]:
    """The acceptor-upstream search window of every intron: up to ``max_window``
    intronic bases ending just before the terminal AG (adaptive — short introns
    contribute whatever they have). Windows shorter than ``width`` or containing a
    non-ACGT base are dropped."""
    out: list[str] = []
    for m in models:
        if m.seqid not in genome:
            continue
        for iseq in m.intron_seqs(genome):
            n = len(iseq)
            if n < width + _AG:
                continue
            w = iseq[max(0, n - _AG - max_window) : n - _AG]
            if len(w) >= width and set(w) <= _ACGT_SET:
                out.append(w)
    return out


def _background(windows: Sequence[str]) -> dict[str, float]:
    counts = dict.fromkeys(_ACGT, 0.0)
    for w in windows:
        for b in w:
            counts[b] += 1
    total = sum(counts.values()) or 1.0
    return {b: counts[b] / total for b in _ACGT}


def _seed_pwm(width: int, anchor: int, background: dict[str, float]) -> list[dict[str, float]]:
    """Seed every position at the background composition except the anchor, which
    is set strongly to the invariant branch adenosine (the only universal,
    non-species-specific prior). EM learns the flanking preferences from the data."""
    pwm = [dict(background) for _ in range(width)]
    pwm[anchor] = {"A": 0.94, "C": 0.02, "G": 0.02, "T": 0.02}
    return pwm


def _normalize(counts: dict[str, float]) -> dict[str, float]:
    total = sum(counts.values()) or 1.0
    return {b: counts[b] / total for b in _ACGT}


def _position_odds(window: str, start: int, pwm: list[dict[str, float]],
                   background: dict[str, float], width: int) -> float:
    """Likelihood ratio of the motif starting at ``start`` vs all-background."""
    odds = 1.0
    for m in range(width):
        b = window[start + m]
        odds *= pwm[m][b] / background[b]
    return odds


def _anchored_starts(window: str, width: int, anchor: int, require: str | None) -> list[int]:
    """Candidate motif starts. The branch adenosine is the invariant 2'-OH
    nucleophile, so by default only positions with that base at the anchor slot are
    considered — without this constraint EM drifts to whatever k-mer is most
    over-represented (in T/pyrimidine-rich intron 3' ends, a spurious T/G motif)."""
    starts = range(0, len(window) - width + 1)
    if require is None:
        return list(starts)
    return [j for j in starts if window[j + anchor] == require]


_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def fit_branch_em(
    windows: Sequence[str],
    *,
    width: int = BRANCH_WIDTH,
    anchor: int = BRANCH_A,
    require_anchor: str | None = "A",
    max_iter: int = 40,
    pseudo: float = 0.1,
    tol: float = 1e-4,
    max_fit: int = 5000,
    seed: int = 0,
) -> BranchModel:
    """Fit the branch motif by EM (OOPS). Each iteration: E-step computes, per
    window, the posterior that the branch sits at each *anchor-consistent* start
    (odds vs background, normalised over those starts); M-step re-estimates the PWM
    from the posterior-weighted base counts. The background is fixed at the observed
    composition. Converges when the PWM stops moving (max abs change < ``tol``).

    ``require_anchor='A'`` (default) restricts candidates to A-anchored positions —
    the biological invariant that keeps EM on the real branch rather than the most
    frequent k-mer. Set ``None`` to discover an unconstrained motif. At most
    ``max_fit`` windows (sampled) are used to fit the PWM; that is ample for a
    7-column model and bounds the runtime.

    Encodes windows to int arrays and precomputes the anchor-consistent starts once
    (they don't change across iterations) so the hot loop is plain arithmetic."""
    usable = [w for w in windows if len(w) >= width and set(w) <= _ACGT_SET]
    if not usable:
        raise ValueError("no usable branch windows")
    background = _background(usable)

    fit = usable
    if len(usable) > max_fit:
        import random

        fit = random.Random(seed).sample(usable, max_fit)

    # encode once: each window -> (int-coded bases, list of anchor-consistent starts)
    a_idx = _IDX.get(require_anchor) if require_anchor is not None else None
    encoded: list[tuple[list[int], list[int]]] = []
    for w in fit:
        codes = [_IDX[b] for b in w]
        starts = [
            j for j in range(len(w) - width + 1)
            if a_idx is None or codes[j + anchor] == a_idx
        ]
        if starts:
            encoded.append((codes, starts))

    bg = [background[b] for b in _ACGT]
    pwm = [[p[b] for b in _ACGT] for p in _seed_pwm(width, anchor, background)]

    for _ in range(max_iter):
        ratio = [[pwm[m][b] / bg[b] for b in range(4)] for m in range(width)]
        counts = [[pseudo] * 4 for _ in range(width)]
        for codes, starts in encoded:
            odds = []
            for j in starts:
                o = 1.0
                for m in range(width):
                    o *= ratio[m][codes[j + m]]
                odds.append(o)
            z = sum(odds)
            if z <= 0:
                continue
            for j, o in zip(starts, odds):
                gamma = o / z
                for m in range(width):
                    counts[m][codes[j + m]] += gamma
        new_pwm = []
        for col in counts:
            s = sum(col) or 1.0
            new_pwm.append([c / s for c in col])
        delta = max(abs(new_pwm[m][b] - pwm[m][b]) for m in range(width) for b in range(4))
        pwm = new_pwm
        if delta < tol:
            break

    pwm_dicts = [dict(zip(_ACGT, col)) for col in pwm]
    return BranchModel(pwm_dicts, background, width, anchor, require_anchor)


def locate_branch(window: str, model: BranchModel) -> tuple[int, float] | None:
    """Best (A-anchored) branch start position in ``window`` and its log-odds
    score, or ``None`` if no anchor-consistent position exists."""
    import math

    starts = _anchored_starts(window, model.width, model.anchor, model.require_anchor)
    best_j, best_odds = None, -1.0
    for j in starts:
        o = _position_odds(window, j, model.pwm, model.background, model.width)
        if best_j is None or o > best_odds:
            best_j, best_odds = j, o
    if best_j is None:
        return None
    return best_j, math.log(best_odds) if best_odds > 0 else float("-inf")


def branch_distances(windows: Sequence[str], model: BranchModel) -> list[int]:
    """Distance (nt) from the located branch adenosine to the 3' splice site, per
    window. The window ends ``_AG`` nt before the 3'SS, so for a branch starting at
    ``j`` the branch-A distance is ``len(window) + _AG - (j + anchor)``."""
    dists = []
    for w in windows:
        if len(w) < model.width or not set(w) <= _ACGT_SET:
            continue
        hit = locate_branch(w, model)
        if hit is None:
            continue
        j, _ = hit
        dists.append(len(w) + _AG - (j + model.anchor))
    return dists


@dataclass
class BranchDistances:
    """geneid branch-distance knobs estimated from the data: minimum, optimal
    (penalty-free) and the acceptor-context scan span."""

    acc_context: int
    min_dist: int
    opt_dist: int


def distance_knobs(distances: Sequence[int], *, offset: int = BRANCH_A) -> BranchDistances:
    """Set the branch-distance knobs from the observed branch→3'SS distribution:
    ``opt_dist`` at the median, ``min_dist`` at the 5th percentile, ``acc_context``
    spanning the 95th (so the scan reaches essentially every real branch). geneid
    requires ``acc_context - offset - opt_dist > 0``, so ``acc_context`` is floored
    accordingly."""
    d = sorted(distances)
    n = len(d)
    if n == 0:
        raise ValueError("no branch distances")

    def pct(p: float) -> int:
        return d[min(n - 1, max(0, int(p * n)))]

    opt = pct(0.50)
    lo = max(1, pct(0.05))
    hi = max(pct(0.95), opt)
    acc = max(hi + 5, offset + opt + 1)
    return BranchDistances(acc_context=acc, min_dist=lo, opt_dist=opt)


# --- emit the geneid Branch_point_profile ------------------------------------
#
# geneid's U2 branch profile header is `len offset cutoff order a b acc_context
# min_dist opt_dist pen_scale` (readparam ReadProfile). We emit an order-0 log-odds
# PWM: offset = the branch-A index (so geneid's PositionBP lands on the branch A),
# order 0, and the distance knobs come from distance_knobs() — set from the data,
# not searched. The score contribution is governed separately by the param scalar
# Branch_point_score_weight (see geneid BuildAcceptors.c); default 0 = the branch
# is scored and reported (bp_score/bp_pos) without affecting splice-site selection.

BRANCH_CUTOFF = -20.0  # permissive: let the branch always be located/reported
BRANCH_PEN_SCALE = 6  # quadratic distance-penalty scale (geneid default)


def _fmt(v: float) -> str:
    return str(int(v)) if v == int(v) else f"{v:g}"


def branch_profile_lines(
    model: BranchModel,
    knobs: BranchDistances,
    *,
    cutoff: float = BRANCH_CUTOFF,
    pen_scale: int = BRANCH_PEN_SCALE,
) -> list[str]:
    """geneid ``Branch_point_profile`` data lines: the header then per-position
    ``pos base log-odds`` rows (order-0 log-ratio of the motif PWM vs background)."""
    import math

    header = [
        model.width, model.anchor, cutoff, 0, 0, 1,
        knobs.acc_context, knobs.min_dist, knobs.opt_dist, pen_scale,
    ]
    lines = [" ".join(_fmt(float(x)) for x in header),
             "# Transition probabilities at every position"]
    for pos in range(model.width):
        for base in _ACGT:
            lo = math.log(model.pwm[pos][base] / model.background[base])
            lines.append(f"{pos + 1} {base} {_fmt(round(lo, 6))}")
    return lines


def branch_profile_section(model: BranchModel, knobs: BranchDistances, **kw) -> str:
    """The full ``Branch_point_profile`` section text, ready to splice in before
    ``Acceptor_profile`` (via :meth:`Param.insert_text_before`)."""
    return "Branch_point_profile\n" + "\n".join(branch_profile_lines(model, knobs, **kw)) + "\n"


def branch_weight_scalar(weight: float) -> str:
    """The optional ``Branch_point_score_weight`` scalar section (0 = report-only)."""
    return f"Branch_point_score_weight\n{_fmt(weight)}\n"
