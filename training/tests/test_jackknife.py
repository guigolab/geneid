import os
from pathlib import Path

import pytest

from geneid_train.jackknife import _write_fasta, make_folds

from .conftest import ref_dir


def test_make_folds_partitions_without_overlap():
    names = [f"g{i:02d}" for i in range(23)]
    folds = make_folds(names, k=10)
    # size = 23//10 + 1 = 3  -> ceil(23/3) = 8 folds
    assert all(len(f) <= 3 for f in folds)
    union: set[str] = set()
    for f in folds:
        assert not (union & f)  # disjoint
        union |= f
    assert union == set(names)  # complete coverage


def test_make_folds_small_set():
    folds = make_folds(["b", "a", "c"], k=10)  # size = 0+1 = 1 -> 3 folds
    assert [sorted(f) for f in folds] == [["a"], ["b"], ["c"]]


def test_write_fasta_selects_and_wraps(tmp_path):
    records = {"x": "ACGT" * 20, "y": "TTTT", "z": "GG"}
    path = tmp_path / "out.fa"
    n = _write_fasta(records, {"x", "z"}, path)
    assert n == 2
    text = path.read_text()
    assert ">x\n" in text and ">z\n" in text and ">y\n" not in text
    # x (80 bp) wraps at 60
    assert "ACGT" * 15 + "\n" in text


# ---- end-to-end cross-validation with a real geneid binary ------------------

REF = ref_dir()
GENEID = os.environ.get("GENEID_BIN") or "/Users/talioto/repositories/geneid_fresh/bin/geneid"
GENOME = (
    REF.parents[2] / "xgXerMont_curated.no_mt.scrubbed.fa.gz" if REF else None
)


@pytest.mark.integration
@pytest.mark.skipif(
    REF is None or not (GENOME and GENOME.exists()) or not Path(GENEID).exists(),
    reason="needs GENEID_TRAIN_REFDIR, the genome, and a geneid binary",
)
def test_jackknife_runs_and_aggregates():
    from geneid_train.core.fasta import read_fasta, read_fasta_subset
    from geneid_train.core.gff import read_gff3
    from geneid_train.evaluate import read_annotation_gff
    from geneid_train.jackknife import jackknife
    from geneid_train.prepare.base import (
        build_models,
        collapse_isoforms,
        filter_complete,
        filter_non_overlapping,
    )

    sp = "Xerocrassa_montserratensis"
    gff3 = REF.parent / "get_candidates" / "good_candidates.1trans.gff3"
    recs = read_gff3(gff3)
    models = collapse_isoforms(build_models(recs), recs)
    train_loci = {
        ln.strip()[1:]
        for ln in open(REF / f"{sp}.train.gp.fa")
        if ln.startswith(">")
    }
    want = {f"SUPER_{i}" for i in range(1, 9)}
    models = [m for m in models if m.gene_id in train_loci and m.seqid in want]
    genome = read_fasta_subset(str(GENOME), {m.seqid for m in models})
    models = filter_non_overlapping(filter_complete(models, genome))

    locus_fasta = read_fasta(REF / f"{sp}.train.gp.fa")
    annots = read_annotation_gff(REF / f"{sp}.train.gp.gff")

    acc = jackknife(
        models, genome, sp, locus_fasta, annots,
        geneid_bin=GENEID, folds=2, weights=(-4.5, 0.35), seed=1,
    )
    # each fold retrained and predicted -> non-empty aggregated counts
    assert acc.totals.tp > 0
    assert acc.totals.exr > 0 and acc.totals.exp > 0
    assert 0.0 <= acc.snsp <= 1.0
