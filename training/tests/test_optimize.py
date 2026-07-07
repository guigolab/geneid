import os
from pathlib import Path

import pytest

from geneid_train.core.param import Param
from geneid_train.evaluate import Accuracy
from geneid_train.optimize import (
    BRANCH_BOUNDS,
    BranchKnobs,
    GridResult,
    SearchSpace,
    WeightPoint,
    _frange,
    accuracy_key,
    apply_weight_point,
    apply_weights,
    latin_hypercube,
    optimize,
    set_branch_knobs,
    set_intron_length_weight,
    uniform_point,
)

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


def test_accuracy_key_orders_by_snsp_then_snspg():
    assert accuracy_key(_acc(0.80)) < accuracy_key(_acc(0.70))
    assert accuracy_key(_acc(0.80, snspg=0.5)) < accuracy_key(_acc(0.80, snspg=0.4))


def test_uniform_point_sets_all_types_equal():
    p = uniform_point(-3.5, 0.3)
    assert p.ewf == (-3.5, -3.5, -3.5, -3.5)
    assert p.owf == (0.3, 0.3, 0.3, 0.3)


def test_weight_point_with_value_is_immutable_single_change():
    p = uniform_point(-4.0, 0.3)
    q = p.with_value("ewf", 2, -2.0)  # Terminal only
    assert q.ewf == (-4.0, -4.0, -2.0, -4.0)
    assert p.ewf == (-4.0, -4.0, -4.0, -4.0)  # original untouched
    assert q.owf == p.owf


def test_latin_hypercube_stratified_and_in_bounds():
    bounds = [(-6.0, 0.0), (0.1, 0.7)]
    pts = latin_hypercube(bounds, n=10, seed=3)
    assert len(pts) == 10
    for d, (lo, hi) in enumerate(bounds):
        col = sorted(p[d] for p in pts)
        assert all(lo <= x <= hi for x in col)
        # one sample per stratum: each of the 10 equal bins holds exactly one point
        width = (hi - lo) / 10
        bins = {min(9, int((x - lo) / width)) for x in col}
        assert len(bins) == 10


def test_latin_hypercube_seed_reproducible():
    b = [(-6.0, 0.0), (0.1, 0.7)]
    assert latin_hypercube(b, 8, seed=5) == latin_hypercube(b, 8, seed=5)
    assert latin_hypercube(b, 8, seed=5) != latin_hypercube(b, 8, seed=6)


def test_searchspace_decode_and_clip():
    space = SearchSpace(ewf_bounds=(-6.0, 0.0), owf_bounds=(0.1, 0.7))
    assert space.bounds() == [(-6.0, 0.0)] * 4 + [(0.1, 0.7)] * 4
    # out-of-box values are clipped per axis
    v = space.clip([-9, 1, -3, -4, 0.9, 0.0, 0.3, 0.4])
    assert v == [-6.0, 0.0, -3, -4, 0.7, 0.1, 0.3, 0.4]
    point = space.decode(v)
    assert point.ewf == (-6.0, 0.0, -3, -4)
    assert point.owf == (0.7, 0.1, 0.3, 0.4)


def test_set_branch_knobs_preserves_head_and_clamps_opt_dist():
    p = Param.from_text("U12_Branch_point_profile\n12 9 2.5 2 0 1 50 7 17 6\n1 AAA 0\n")
    set_branch_knobs(p, "U12_Branch_point_profile", BranchKnobs(45, 8, 20, 5.0))
    hdr = p.vector("U12_Branch_point_profile")
    # len/offset/cutoff/order/a/b preserved; knobs replaced
    assert hdr[:6] == ["12", "9", "2.5", "2", "0", "1"]
    assert hdr[6:] == ["45", "8", "20", "5"]
    # opt_dist is clamped so acc_context - offset - opt_dist stays positive:
    # acc_context 40, offset 9 -> opt_dist must be < 31
    set_branch_knobs(p, "U12_Branch_point_profile", BranchKnobs(40, 8, 99, 5.0))
    assert p.vector("U12_Branch_point_profile")[8] == "30"


def test_searchspace_branch_axes_appended_and_decoded():
    space = SearchSpace(branch_profiles=("U12_Branch_point_profile",))
    b = space.bounds()
    assert len(b) == 12  # 8 weights + 4 branch axes
    assert b[8:] == list(BRANCH_BOUNDS)
    v = [-4] * 4 + [0.3] * 4 + [45.4, 7.6, 18.2, 6.0]
    assert space.decode(v).ewf == (-4, -4, -4, -4)  # branch dims don't leak into owf
    (name, kn), = space.decode_branch(v)
    assert name == "U12_Branch_point_profile"
    assert (kn.acc_context, kn.min_dist, kn.opt_dist) == (45, 8, 18)  # rounded to int


def test_set_intron_length_weight_sets_scalar():
    p = Param.from_text("Intron_length_model\n7.2 1.5\nIntron_length_score_weight\n0\n")
    set_intron_length_weight(p, 0.75)
    assert p.scalar("Intron_length_score_weight") == "0.75"


def test_set_intron_length_weight_noop_when_absent():
    # a param without the section is left untouched (older params predate the feature)
    p = Param.from_text(_MINI_PARAM)
    set_intron_length_weight(p, 1.5)  # must not raise
    assert "Intron_length_score_weight" not in p.keywords()


def test_searchspace_intron_length_axis_appended_and_decoded():
    space = SearchSpace(tune_intron_length=True)
    b = space.bounds()
    assert len(b) == 9  # 8 weights + 1 lambda axis
    assert b[8] == space.intron_length_bounds
    v = [-4] * 4 + [0.3] * 4 + [1.234]
    assert space.decode(v).owf == (0.3, 0.3, 0.3, 0.3)  # lambda doesn't leak into owf
    assert space.decode_intron_length(v) == 1.234
    assert SearchSpace().decode_intron_length(v) is None  # off by default


def test_searchspace_branch_and_intron_length_axes_order():
    space = SearchSpace(branch_profiles=("U12_Branch_point_profile",), tune_intron_length=True)
    assert len(space.bounds()) == 13  # 8 weights + 4 branch + 1 lambda
    v = [-4] * 4 + [0.3] * 4 + [45.0, 8.0, 18.0, 6.0] + [0.9]
    (name, kn), = space.decode_branch(v)  # branch axes still read from v[8:12]
    assert kn.acc_context == 45
    assert space.decode_intron_length(v) == 0.9  # lambda is the last axis


def test_apply_weight_point_per_type_columns_preserve_utr():
    p = Param.from_text(_MINI_PARAM)
    apply_weight_point(p, WeightPoint((-1, -2, -3, -4), (0.1, 0.2, 0.3, 0.4)))
    assert p.vector("Exon_weights") == ["-1", "-2", "-3", "-4", "0"]  # UTR 0 kept
    assert p.vector("Exon_factor") == ["0.1", "0.2", "0.3", "0.4"]
    # Site_factor = 1 - owf per type, trailing UTR 0.6 preserved
    assert p.vector("Site_factor") == ["0.9", "0.8", "0.7", "0.6", "0.6"]


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


@pytest.mark.integration
@pytest.mark.skipif(
    REF is None or not Path(GENEID).exists(),
    reason="needs GENEID_TRAIN_REFDIR and a geneid binary (GENEID_BIN)",
)
def test_coordinate_descent_beats_or_matches_uniform_seed():
    from geneid_train.optimize import coordinate_descent

    sp = "Xerocrassa_montserratensis"
    base = (REF / f"{sp}.geneid.param").read_text()
    fasta = str(REF / f"{sp}.eval.gp.fa")
    gff = str(REF / f"{sp}.eval.gp.gff")
    seed = uniform_point(-4.5, 0.35)
    opt_text, best, history = coordinate_descent(
        base, fasta, gff, geneid_bin=GENEID, init=seed,
        ewf_values=[-4.5], owf_values=[0.30, 0.35], workers=2, max_rounds=1,
    )
    # descent only ever accepts improvements over the seed, so best >= seed
    assert best.accuracy.snsp >= history[0].accuracy.snsp
    assert len(best.point.ewf) == 4 and len(best.point.owf) == 4
    # optimised param carries the per-type owf (First may differ from the rest)
    p = Param.from_text(opt_text)
    assert p.vector("Exon_factor") == [f"{o:g}" for o in best.point.owf]


@pytest.mark.integration
@pytest.mark.skipif(
    REF is None or not Path(GENEID).exists(),
    reason="needs GENEID_TRAIN_REFDIR and a geneid binary (GENEID_BIN)",
)
def test_global_optimize_returns_valid_result():
    from geneid_train.optimize import SearchSpace, global_optimize

    sp = "Xerocrassa_montserratensis"
    base = (REF / f"{sp}.geneid.param").read_text()
    fasta = str(REF / f"{sp}.eval.gp.fa")
    gff = str(REF / f"{sp}.eval.gp.gff")
    # tiny sample + small eval budget keep the test quick
    opt_text, res = global_optimize(
        base, fasta, gff, geneid_bin=GENEID, space=SearchSpace((-5.0, -2.0), (0.2, 0.5)),
        n_samples=4, step=(1.0, 0.1), min_step=(0.5, 0.05), workers=4, seed=1, max_evals=24,
    )
    # LHS(4) + at least one compass sweep; a sweep adds 2*ndim at once so the
    # budget is a soft cap (it stops checking at max_evals, not mid-sweep)
    assert res.n_evaluations >= 4
    assert 0.0 <= res.accuracy.snsp <= 1.0
    assert len(res.point.ewf) == 4 and len(res.point.owf) == 4
    # every weight stays inside the search box
    assert all(-5.0 <= e <= -2.0 for e in res.point.ewf)
    assert all(0.2 <= o <= 0.5 for o in res.point.owf)
    p = Param.from_text(opt_text)
    assert p.vector("Exon_weights")[:4] == [f"{e:g}" for e in res.point.ewf]
