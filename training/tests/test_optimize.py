import os
from pathlib import Path

import pytest

from geneid_train.core.param import Param
from geneid_train.evaluate import Accuracy
from geneid_train.optimize import GridResult, _frange, apply_weights, optimize

from .conftest import ref_dir

_MINI_PARAM = (
    "number_of_isochores\n1\n"
    "Exon_weights\n-4 -4 -4 -4 0\n"
    "Exon_factor\n0.4 0.4 0.4 0.4\n"
    "Site_factor\n0.6 0.6 0.6 0.6 0.6\n"
)


def test_frange_inclusive_of_final():
    assert _frange(-4.5, -2.5, 0.5) == [-4.5, -4.0, -3.5, -3.0, -2.5]
    assert _frange(0.25, 0.50, 0.05) == [0.25, 0.30, 0.35, 0.40, 0.45, 0.50]


def test_apply_weights_sets_leading_and_preserves_trailing():
    p = Param.from_text(_MINI_PARAM)
    apply_weights(p, ewf=-3.5, owf=0.3)
    assert p.vector("Exon_weights") == ["-3.5", "-3.5", "-3.5", "-3.5", "0"]
    assert p.vector("Exon_factor") == ["0.3", "0.3", "0.3", "0.3"]
    # site factor = 1 - owf, trailing 0.6 preserved
    assert p.vector("Site_factor") == ["0.7", "0.7", "0.7", "0.7", "0.6"]


def _acc(snsp, snspg=0.0, cc=0.0, ra_me=0.0, ra_we=0.0):
    return Accuracy(
        sn=0, sp=0, cc=cc, sne=snsp, spe=snsp, snsp=snsp,
        sng=snspg, spg=snspg, snspg=snspg, ra_me=ra_me, ra_we=ra_we,
    )


def test_grid_result_sort_prefers_higher_snsp_then_tiebreaks():
    a = GridResult(-4.5, 0.3, _acc(0.80))
    b = GridResult(-4.0, 0.3, _acc(0.70))
    c = GridResult(-3.5, 0.3, _acc(0.80, snspg=0.5))  # ties a on SNSP, wins on SNSPg
    results = sorted([a, b, c], key=lambda r: r.sort_key())
    assert results[0] is c
    assert results[1] is a
    assert results[2] is b


# ---- end-to-end with a real geneid binary -----------------------------------

REF = ref_dir()
GENEID = os.environ.get("GENEID_BIN") or (
    "/Users/talioto/repositories/geneid_fresh/bin/geneid"
)


@pytest.mark.integration
@pytest.mark.skipif(
    REF is None or not Path(GENEID).exists(),
    reason="needs GENEID_TRAIN_REFDIR and a geneid binary (GENEID_BIN)",
)
def test_optimize_runs_grid_and_selects_best():
    sp = "Xerocrassa_montserratensis"
    base = (REF / f"{sp}.geneid.param").read_text()
    fasta = str(REF / f"{sp}.eval.gp.fa")
    gff = str(REF / f"{sp}.eval.gp.gff")
    # a tiny 2x2 grid keeps the test quick
    opt_text, results = optimize(
        base, fasta, gff, geneid_bin=GENEID,
        ewf_grid=(-4.5, -4.0, 0.5), owf_grid=(0.30, 0.35, 0.05), workers=2,
    )
    assert len(results) == 4
    # results are sorted best-first by SNSP
    snsps = [r.accuracy.snsp for r in results]
    assert snsps == sorted(snsps, reverse=True)
    # the optimised param carries the winning weights
    best = results[0]
    p = Param.from_text(opt_text)
    assert p.vector("Exon_weights")[0] == f"{best.ewf:g}"
    assert p.vector("Exon_factor")[0] == f"{best.owf:g}"
