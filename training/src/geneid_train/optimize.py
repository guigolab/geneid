"""Parameter optimisation: grid-search the exon/site weights (replaces the
legacy OptimizeParameter loop).

For each ``(eWF, oWF)`` grid point the exon-weight / exon-factor / site-factor
fields of the parameter file are set, geneid is run on a held-out evaluation set,
and the prediction is scored with :mod:`geneid_train.evaluate`. The point with
the best exon-level ``SNSP`` (tie-broken by gene ``SNSPg``, nucleotide ``CC``,
then fewest missing/wrong exons — the legacy ``sorteval`` order) wins, and its
weights are written into the optimised parameter file.

Only the default eWF×oWF grid is implemented (matching the reference run); the
branch/U12 acceptor-context + min-branch axes are left for the U12 work.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from .core.param import Param
from .evaluate import Accuracy, evaluate_files

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


@dataclass
class GridResult:
    ewf: float
    owf: float
    accuracy: Accuracy

    def sort_key(self) -> tuple:
        a = self.accuracy
        # SNSP desc, SNSPg desc, CC desc, raME asc, raWE asc  (legacy sorteval)
        return (-a.snsp, -a.snspg, -a.cc, a.ra_me, a.ra_we)


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
