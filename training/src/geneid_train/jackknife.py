"""Leave-group-out cross-validation (replaces the legacy runJacknife).

Estimates a parameter file's real accuracy without a separate held-out set: the
training loci are split into ~10 groups; for each group a fresh parameter file is
trained on all the *other* loci (carrying the optimised exon/site weights) and
scored on the held-out group with geneid + :mod:`geneid_train.evaluate`. Counts
from every fold accumulate into one cross-validated SN/SP estimate.

Retraining every fold is the whole point (it exposes over-fitting), so a fixed
background model is reused across folds to keep them comparable and fast.
"""

from __future__ import annotations

import subprocess
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path

from .core.param import Param
from .evaluate import Accuracy, Locus, Totals, accumulate, finalize, read_prediction_gff
from .optimize import apply_weights, run_geneid
from .prepare.base import GeneModel
from .stats.background import from_genome
from .stats.sites import Matrix
from .train import train


def make_folds(names: Sequence[str], k: int = 10) -> list[set[str]]:
    """Split sorted locus names into ~``k`` contiguous groups of size
    ``len//k + 1`` (matching the legacy group sizing)."""
    ordered = sorted(names)
    n = len(ordered)
    size = n // k + 1
    return [set(ordered[i : i + size]) for i in range(0, n, size)]


def _write_fasta(records: Mapping[str, str], names: set[str], path: Path) -> int:
    written = 0
    with open(path, "w") as fh:
        for name in sorted(names):
            seq = records.get(name)
            if seq is None:
                continue
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
            written += 1
    return written


def jackknife(
    models: Sequence[GeneModel],
    genome: Mapping[str, str],
    species: str,
    locus_fasta: Mapping[str, str],
    annotations: Sequence[Locus],
    *,
    geneid_bin: str = "geneid",
    folds: int = 10,
    weights: tuple[float, float] | None = None,
    background: tuple[Matrix, Matrix] | None = None,
    seed: int = 0,
) -> Accuracy:
    """Run leave-group-out CV and return the aggregated accuracy.

    ``locus_fasta`` maps each locus name to its (flanked) prediction sequence;
    ``annotations`` are the per-locus real exon sets (from
    :func:`evaluate.read_annotation_gff`). ``weights`` = ``(eWF, oWF)`` applies the
    optimised weights to each fold's param before scoring.
    """
    ann_by_name = {locus.name: locus for locus in annotations}
    # one background for every fold, so folds differ only by their training genes
    if background is None:
        background = from_genome(genome.values(), seed=seed)

    fold_groups = make_folds([m.gene_id for m in models], folds)
    totals = Totals()
    with tempfile.TemporaryDirectory() as td:
        workdir = Path(td)
        for i, held in enumerate(fold_groups):
            train_models = [m for m in models if m.gene_id not in held]
            if not train_models:
                continue
            param = Param.from_text(
                train(train_models, genome, species, background=background)
            )
            if weights is not None:
                apply_weights(param, *weights)
            ppath = workdir / f"fold{i}.param"
            param.write(ppath)

            fapath = workdir / f"fold{i}.fa"
            if not _write_fasta(locus_fasta, held, fapath):
                continue
            try:
                pred_text = run_geneid(geneid_bin, str(ppath), str(fapath))
            except subprocess.CalledProcessError:
                continue
            predpath = workdir / f"fold{i}.pred.gff"
            predpath.write_text(pred_text)
            pred = read_prediction_gff(predpath)

            for name in held:
                locus = ann_by_name.get(name)
                if locus is None:
                    continue
                accumulate(pred.get(name, []), locus.exons, locus.length, totals)
    return finalize(totals)
