"""Parameter optimisation: grid-search the exon/site weights (replaces the
legacy OptimizeParameter loop).

For each ``(eWF, oWF)`` grid point the exon-weight / exon-factor / site-factor
fields of the parameter file are set, geneid is run on a held-out evaluation set,
and the prediction is scored with :mod:`geneid_train.evaluate`. The point with
the best exon-level ``SNSP`` (tie-broken by gene ``SNSPg``, nucleotide ``CC``,
then fewest missing/wrong exons — the legacy ``sorteval`` order) wins, and its
weights are written into the optimised parameter file.

Two search strategies are available:

- :func:`optimize` — the legacy uniform grid over a single ``(eWF, oWF)`` applied
  to all exon types (fast, good for a coarse starting point).
- :func:`coordinate_descent` — refines the four exon types (First/Internal/
  Terminal/Single) *independently*: the Exon_weights / Exon_factor / Site_factor
  columns exist to be tuned per type, so this walks the 8-dimensional space one
  coordinate at a time. Seed it with :func:`uniform_point` from a coarse grid.

- :func:`global_optimize` — Latin-hypercube exploration + compass-search
  refinement over the 8-D per-type weight box (:class:`SearchSpace`). When the
  search space names branch profiles, four branch-distance axes per profile
  (acc_context / min_dist / opt_dist / pen_scale) are appended, so the U12 (or
  U2) branch knobs are tuned by the same search.
"""

from __future__ import annotations

import itertools
import random
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field, replace
from pathlib import Path

from .core.param import Param
from .evaluate import Accuracy, evaluate_files

# the four exon types tuned independently (columns 0-3 of Exon_weights /
# Exon_factor / Site_factor; column 4, where present, is UTR and is left alone)
EXON_TYPES = ("First", "Internal", "Terminal", "Single")

# default grid bounds from the legacy driver (init, final, step)
DEFAULT_EWF = (-4.5, -2.5, 0.5)
DEFAULT_OWF = (0.25, 0.50, 0.05)

_TYPED_EXONS = ("First", "Internal", "Terminal", "Single")


def _frange(init: float, final: float, step: float) -> list[float]:
    vals, v = [], init
    while v <= final + 1e-9:
        vals.append(round(v, 6))
        v += step
    return vals


def _set_all(param: Param, keyword: str, value: float, n: int = 4) -> None:
    """Set the first ``n`` fields of every occurrence of a vector section (one per
    isochore) to ``value``, preserving trailing fields (e.g. the single-exon
    column)."""
    count = sum(1 for k in param.keywords() if k == keyword)
    for i in range(count):
        old = param.vector(keyword, index=i)
        k = min(n, len(old))
        param.set_scalar(keyword, " ".join([f"{value:g}"] * k + old[k:]), index=i)


def apply_weights(param: Param, ewf: float, owf: float) -> None:
    """Set Exon_weights=[ewf]*4, Exon_factor=[owf]*4, Site_factor=[1-owf]*4 across
    all isochores."""
    _set_all(param, "Exon_weights", ewf)
    _set_all(param, "Exon_factor", owf)
    _set_all(param, "Site_factor", 1 - owf)


def run_geneid(geneid_bin: str, param_path: str, fasta: str) -> str:
    """Run ``geneid -GP param fasta`` and return the typed-CDS-exon prediction
    lines (First/Internal/Terminal/Single), dropping comments and generic exons."""
    out = subprocess.run(
        [geneid_bin, "-GP", param_path, fasta],
        capture_output=True, text=True, check=True,
    ).stdout
    kept = []
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) >= 9 and f[2] in _TYPED_EXONS:
            kept.append(line)
    return "\n".join(kept) + "\n"


def accuracy_key(a: Accuracy) -> tuple:
    """Ranking key: exon SNSP desc, then gene SNSPg desc, nucleotide CC desc,
    fewest missing then wrong exons (the legacy ``sorteval`` order)."""
    return (-a.snsp, -a.snspg, -a.cc, a.ra_me, a.ra_we)


@dataclass
class GridResult:
    ewf: float
    owf: float
    accuracy: Accuracy

    def sort_key(self) -> tuple:
        return accuracy_key(self.accuracy)


def _score_point(
    base_text: str, ewf: float, owf: float, geneid_bin: str, fasta: str, gff: str, workdir: Path
) -> GridResult:
    param = Param.from_text(base_text)
    apply_weights(param, ewf, owf)
    ptmp = workdir / f"param_{ewf}_{owf}.tmp"
    param.write(ptmp)
    pred = workdir / f"pred_{ewf}_{owf}.gff"
    pred.write_text(run_geneid(geneid_bin, str(ptmp), fasta))
    acc = evaluate_files(str(pred), gff)
    return GridResult(ewf, owf, acc)


def optimize(
    base_param_text: str,
    eval_fasta: str,
    eval_gff: str,
    *,
    geneid_bin: str = "geneid",
    ewf_grid: tuple[float, float, float] = DEFAULT_EWF,
    owf_grid: tuple[float, float, float] = DEFAULT_OWF,
    workers: int = 4,
) -> tuple[str, list[GridResult]]:
    """Grid-search eWF×oWF; return the optimised param text and all results
    (best first). geneid runs in parallel across grid points."""
    bin_path = shutil.which(geneid_bin) or geneid_bin
    ewfs = _frange(*ewf_grid)
    owfs = _frange(*owf_grid)
    results: list[GridResult] = []
    with tempfile.TemporaryDirectory() as td:
        workdir = Path(td)
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futs = [
                pool.submit(
                    _score_point, base_param_text, e, o, bin_path, eval_fasta, eval_gff, workdir
                )
                for e in ewfs
                for o in owfs
            ]
            results = [f.result() for f in futs]
    results.sort(key=lambda r: r.sort_key())
    best = results[0]
    param = Param.from_text(base_param_text)
    apply_weights(param, best.ewf, best.owf)
    return param.to_text(), results


# --- per-exon-type coordinate-descent optimisation ---------------------------


@dataclass(frozen=True)
class WeightPoint:
    """Per-exon-type exon weights (``ewf``) and factors (``owf``), ordered
    First, Internal, Terminal, Single. ``Exon_factor`` takes ``owf`` and
    ``Site_factor`` takes ``1 - owf`` per type, mirroring the uniform search."""

    ewf: tuple[float, float, float, float]
    owf: tuple[float, float, float, float]

    def with_value(self, axis: str, idx: int, value: float) -> WeightPoint:
        vals = list(getattr(self, axis))
        vals[idx] = value
        return replace(self, **{axis: tuple(vals)})


def uniform_point(ewf: float, owf: float) -> WeightPoint:
    """A WeightPoint with every exon type set to the same ``ewf``/``owf`` — the
    starting point equivalent to the uniform grid's result."""
    return WeightPoint((ewf,) * 4, (owf,) * 4)


def _set_columns(param: Param, keyword: str, values: tuple[float, ...]) -> None:
    """Set columns 0..len(values)-1 of every occurrence of a vector section,
    preserving trailing columns (e.g. the UTR column)."""
    count = sum(1 for k in param.keywords() if k == keyword)
    for i in range(count):
        old = param.vector(keyword, index=i)
        new = [f"{v:g}" for v in values] + old[len(values):]
        param.set_scalar(keyword, " ".join(new), index=i)


def apply_weight_point(param: Param, point: WeightPoint) -> None:
    """Write a per-type WeightPoint into Exon_weights / Exon_factor / Site_factor."""
    _set_columns(param, "Exon_weights", point.ewf)
    _set_columns(param, "Exon_factor", point.owf)
    _set_columns(param, "Site_factor", tuple(round(1 - o, 6) for o in point.owf))


# --- branch-point distance knobs (U12 / U2 branch profiles) ------------------
#
# A branch-point profile header is ``len offset cutoff order a b acc_context
# min_dist opt_dist pen_scale`` (readparam.c ReadProfile). The last four fields
# are the distance knobs geneid uses to score the branch: it scans the window
# [acc-acc_context, acc-min_dist] upstream of the acceptor for the best branch,
# subtracting a quadratic penalty ``pen_scale * ((|d-opt_dist| / (acc_context -
# offset - opt_dist))**2)`` (BuildAcceptors.c ComputeU2BranchProfile). The legacy
# optimiser only gridded acc_context and min_dist; opt_dist and pen_scale were
# left at defaults — this exposes all four to the same search.


@dataclass(frozen=True)
class BranchKnobs:
    """The four branch-point distance parameters of a branch profile header."""

    acc_context: int
    min_dist: int
    opt_dist: int
    pen_scale: float


def set_branch_knobs(param: Param, profile: str, knobs: BranchKnobs) -> None:
    """Rewrite the distance knobs (fields 6-9) of every occurrence of a branch
    ``profile`` header, preserving len/offset/cutoff/order/a/b. ``opt_dist`` is
    clamped so geneid's penalty denominator ``acc_context - offset - opt_dist``
    stays positive."""
    count = sum(1 for k in param.keywords() if k == profile)
    for i in range(count):
        block = param._find(profile, i)
        di = block.first_data_index()
        if di is None:
            continue
        fields = block.raw[di].split()
        # ensure the 6 leading fields exist (a=0, b=1 defaults if a bare header)
        while len(fields) < 6:
            fields.append("0" if len(fields) == 4 else "1")
        offset = int(float(fields[1]))
        opt = min(knobs.opt_dist, knobs.acc_context - offset - 1)
        head = fields[:6]
        tail = [str(knobs.acc_context), str(knobs.min_dist), str(opt), f"{knobs.pen_scale:g}"]
        ending = "\n" if block.raw[di].endswith("\n") else ""
        block.raw[di] = " ".join(head + tail) + ending


@dataclass
class CDResult:
    point: WeightPoint
    accuracy: Accuracy


def _score_weight_point(
    base_text: str, point: WeightPoint, geneid_bin: str,
    fasta: str, gff: str, workdir: Path, tag: str,
    branch: tuple[tuple[str, BranchKnobs], ...] = (),
) -> Accuracy:
    param = Param.from_text(base_text)
    apply_weight_point(param, point)
    for name, knobs in branch:
        set_branch_knobs(param, name, knobs)
    ptmp = workdir / f"p_{tag}.param"
    param.write(ptmp)
    pred = workdir / f"pred_{tag}.gff"
    pred.write_text(run_geneid(geneid_bin, str(ptmp), fasta))
    return evaluate_files(str(pred), gff)


def coordinate_descent(
    base_param_text: str,
    eval_fasta: str,
    eval_gff: str,
    *,
    geneid_bin: str = "geneid",
    init: WeightPoint | None = None,
    ewf_values: list[float] | None = None,
    owf_values: list[float] | None = None,
    workers: int = 4,
    max_rounds: int = 3,
) -> tuple[str, CDResult, list[CDResult]]:
    """Refine per-exon-type weights by coordinate descent, maximising held-out
    exon SNSP. Each of the 8 coordinates (4 eWF + 4 oWF) is line-searched in turn
    (its candidate values run in parallel) holding the others fixed; rounds repeat
    until no coordinate improves or ``max_rounds`` is reached.

    Returns ``(optimised_param_text, best CDResult, history of improvements)``.
    Seed with :func:`uniform_point` from a coarse :func:`optimize` for a good start.
    """
    bin_path = shutil.which(geneid_bin) or geneid_bin
    ewf_values = ewf_values if ewf_values is not None else _frange(*DEFAULT_EWF)
    owf_values = owf_values if owf_values is not None else _frange(*DEFAULT_OWF)
    point = init if init is not None else uniform_point(-4.0, 0.30)

    counter = itertools.count()
    with tempfile.TemporaryDirectory() as td:
        workdir = Path(td)

        def score(p: WeightPoint) -> Accuracy:
            return _score_weight_point(
                base_param_text, p, bin_path, eval_fasta, eval_gff, workdir,
                str(next(counter)),
            )

        best_acc = score(point)
        history = [CDResult(point, best_acc)]
        for _ in range(max_rounds):
            improved = False
            for axis, values in (("ewf", ewf_values), ("owf", owf_values)):
                for idx in range(4):
                    trials = [
                        point.with_value(axis, idx, v)
                        for v in values
                        if v != getattr(point, axis)[idx]
                    ]
                    if not trials:
                        continue
                    with ThreadPoolExecutor(max_workers=workers) as pool:
                        accs = list(pool.map(score, trials))
                    candidates = [(point, best_acc)] + list(zip(trials, accs))
                    candidates.sort(key=lambda c: accuracy_key(c[1]))
                    best_point, best_of = candidates[0]
                    if best_point != point:
                        point, best_acc = best_point, best_of
                        improved = True
                        history.append(CDResult(point, best_acc))
            if not improved:
                break

    param = Param.from_text(base_param_text)
    apply_weight_point(param, point)
    return param.to_text(), CDResult(point, best_acc), history


# --- global search + local refinement over the per-type weight box -----------


# default bounds for the four branch-distance axes (acc_context, min_dist,
# opt_dist, pen_scale). Chosen so the penalty denominator (acc_context - offset -
# opt_dist, offset ~9) stays positive across the whole box: min acc_context 40 >
# 9 + max opt_dist 25. Legacy gridded only acc_context 40-70 and min_dist 7-9.
BRANCH_BOUNDS: tuple[tuple[float, float], ...] = (
    (40.0, 70.0),  # acc_context
    (5.0, 12.0),   # min_dist
    (10.0, 25.0),  # opt_dist
    (2.0, 10.0),   # pen_scale
)


@dataclass
class SearchSpace:
    """The box the global search explores. The first 8 axes are the per-type eWF
    and oWF weights (decoded to a :class:`WeightPoint`). When ``branch_profiles``
    is non-empty, four more axes per profile — acc_context, min_dist, opt_dist,
    pen_scale — are appended and decoded to :class:`BranchKnobs`, so the U12 (or
    U2) branch-distance knobs are tuned by the same search."""

    ewf_bounds: tuple[float, float] = (-6.0, 0.0)
    owf_bounds: tuple[float, float] = (0.10, 0.70)
    branch_profiles: tuple[str, ...] = ()
    branch_bounds: tuple[tuple[float, float], ...] = BRANCH_BOUNDS

    def bounds(self) -> list[tuple[float, float]]:
        base = [self.ewf_bounds] * 4 + [self.owf_bounds] * 4
        return base + list(self.branch_bounds) * len(self.branch_profiles)

    def clip(self, v: list[float]) -> list[float]:
        return [min(hi, max(lo, x)) for x, (lo, hi) in zip(v, self.bounds())]

    def decode(self, v: list[float]) -> WeightPoint:
        return WeightPoint(tuple(v[:4]), tuple(round(x, 6) for x in v[4:8]))

    def decode_branch(self, v: list[float]) -> tuple[tuple[str, BranchKnobs], ...]:
        """Decode the branch axes (v[8:]) into ``(profile, BranchKnobs)`` pairs;
        empty when no branch profiles are being tuned."""
        out = []
        for j, name in enumerate(self.branch_profiles):
            a, m, o, p = v[8 + 4 * j : 12 + 4 * j]
            out.append((name, BranchKnobs(round(a), round(m), round(o), round(p, 3))))
        return tuple(out)


def latin_hypercube(bounds: list[tuple[float, float]], n: int, seed: int = 0) -> list[list[float]]:
    """``n`` Latin-hypercube samples over the box (one stratified draw per axis)."""
    rng = random.Random(seed)
    dims = len(bounds)
    pts = [[0.0] * dims for _ in range(n)]
    for d, (lo, hi) in enumerate(bounds):
        strata = list(range(n))
        rng.shuffle(strata)
        for i in range(n):
            u = (strata[i] + rng.random()) / n
            pts[i][d] = lo + u * (hi - lo)
    return pts


@dataclass
class SearchResult:
    point: WeightPoint
    accuracy: Accuracy
    n_evaluations: int = 0
    history: list[CDResult] = field(default_factory=list)
    branch: tuple[tuple[str, BranchKnobs], ...] = ()


def global_optimize(
    base_param_text: str,
    eval_fasta: str,
    eval_gff: str,
    *,
    geneid_bin: str = "geneid",
    space: SearchSpace | None = None,
    n_samples: int = 32,
    step: tuple[float, float] = (1.0, 0.1),
    min_step: tuple[float, float] = (0.125, 0.0125),
    branch_step: float = 2.0,
    branch_min_step: float = 0.5,
    workers: int = 4,
    seed: int = 0,
    max_evals: int = 400,
) -> tuple[str, SearchResult]:
    """Global Latin-hypercube exploration followed by compass-search refinement.

    ``n_samples`` points are drawn over :class:`SearchSpace` and scored in
    parallel; the best seeds a pattern search that probes each axis at ``±step``
    (eWF, oWF steps; ``branch_step`` for any branch-distance axes), moving to the
    best neighbour and halving the step when a full sweep fails to improve, until
    every step drops below its ``min_step`` or the ``max_evals`` geneid-run budget
    is spent. Ranking is held-out exon SNSP (:func:`accuracy_key`). Returns the
    optimised param text and a :class:`SearchResult` (best weights, branch knobs,
    accuracy, eval count).
    """
    bin_path = shutil.which(geneid_bin) or geneid_bin
    space = space or SearchSpace()
    bounds = space.bounds()
    counter = itertools.count()
    evals = 0

    with tempfile.TemporaryDirectory() as td:
        workdir = Path(td)

        def score(v: list[float]) -> Accuracy:
            cv = space.clip(v)
            return _score_weight_point(
                base_param_text, space.decode(cv), bin_path,
                eval_fasta, eval_gff, workdir, str(next(counter)),
                branch=space.decode_branch(cv),
            )

        def score_many(vs: list[list[float]]) -> list[Accuracy]:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                return list(pool.map(score, vs))

        # --- global: Latin-hypercube sampling ---
        samples = latin_hypercube(bounds, n_samples, seed)
        accs = score_many(samples)
        evals += len(samples)
        best_v, best_acc = min(zip(samples, accs), key=lambda p: accuracy_key(p[1]))
        best_v = list(best_v)
        history = [CDResult(space.decode(space.clip(best_v)), best_acc)]

        # --- local: compass (pattern) search with step halving ---
        n_branch = len(bounds) - 8
        cur_step = [step[0]] * 4 + [step[1]] * 4 + [branch_step] * n_branch
        min_s = [min_step[0]] * 4 + [min_step[1]] * 4 + [branch_min_step] * n_branch
        while evals < max_evals and any(cur_step[d] >= min_s[d] for d in range(len(bounds))):
            neighbours = []
            for d in range(len(bounds)):
                if cur_step[d] < min_s[d]:
                    continue
                for sign in (+1, -1):
                    nb = list(best_v)
                    nb[d] = nb[d] + sign * cur_step[d]
                    neighbours.append(space.clip(nb))
            n_accs = score_many(neighbours)
            evals += len(neighbours)
            cand_v, cand_acc = min(
                zip(neighbours, n_accs), key=lambda p: accuracy_key(p[1])
            )
            if accuracy_key(cand_acc) < accuracy_key(best_acc):
                best_v, best_acc = list(cand_v), cand_acc
                history.append(CDResult(space.decode(best_v), best_acc))
            else:
                cur_step = [s / 2 for s in cur_step]  # contract and refine

    best_v = space.clip(best_v)
    point = space.decode(best_v)
    branch = space.decode_branch(best_v)
    param = Param.from_text(base_param_text)
    apply_weight_point(param, point)
    for name, knobs in branch:
        set_branch_knobs(param, name, knobs)
    return param.to_text(), SearchResult(point, best_acc, evals, history, branch)
