import math

import pytest

from geneid_train.stats.genemodel import (
    format_range,
    intron_length_model,
    intron_range,
)

from .conftest import ref_dir


def test_intron_range_caps_short_at_40():
    # shortest intron 1000 -> 0.75*1000 = 750, capped at 40
    lo, hi = intron_range([1000, 1000, 1000, 1000])
    assert lo == 40.0


def test_intron_range_short_below_cap():
    lo, _ = intron_range([20, 100, 200])  # 0.75*20 = 15 < 40
    assert lo == 15.0


def test_format_range():
    assert format_range(40.0, "Infinity") == "40:Infinity"
    assert format_range(24.75, 25394.023) == "24.75:25394.023"


def test_intron_length_model_recovers_lognormal_params():
    # Draw ln(length) from a fixed grid so mu/sigma are exactly the mean/pop-sd of
    # those logs; the fit must recover them.
    log_vals = [4.0, 5.0, 6.0, 7.0, 8.0]
    lengths = [round(math.exp(v)) for v in log_vals]
    mu, sigma = intron_length_model(lengths)
    logs = [math.log(n) for n in lengths]
    exp_mu = sum(logs) / len(logs)
    exp_sigma = math.sqrt(sum((x - exp_mu) ** 2 for x in logs) / len(logs))
    assert mu == pytest.approx(exp_mu)
    assert sigma == pytest.approx(exp_sigma)


def test_intron_range_max_intron_override():
    lens = [1000] * 999 + [90000]
    lo, hi = intron_range(lens, max_intron=500_000)
    assert hi == 500000.0  # used directly, bypassing the p99.9/long_cap path
    lo2, hi2 = intron_range(lens)  # default is the (much lower) p99.9 estimate
    assert hi2 < hi
    assert lo == lo2  # min is unaffected by the override


def test_intron_length_model_constant_lengths_zero_sigma():
    mu, sigma = intron_length_model([2000, 2000, 2000])
    assert mu == pytest.approx(math.log(2000))
    assert sigma == pytest.approx(0.0)


def test_intron_length_model_ignores_nonpositive_and_empty():
    assert intron_length_model([]) == (0.0, 0.0)
    # non-positive lengths (shouldn't occur) are filtered, not crash on log(0)
    assert intron_length_model([0, -5]) == (0.0, 0.0)


REF = ref_dir()


def test_intron_range_max_is_p999_skew_robust():
    # a right-skewed set: bulk small + a long tail. mean+3sd would sit far above
    # the bulk; p99.9 tracks the actual tail and excludes ~0.1%.
    lens = [1000] * 999 + [90000]
    lo, hi = intron_range(lens)
    over = sum(1 for x in lens if x > hi)
    assert over <= 1  # at most the top ~0.1% excluded
    assert hi <= 100_000


@pytest.mark.integration
@pytest.mark.skipif(REF is None, reason="set GENEID_TRAIN_REFDIR to the train_geneid dir")
def test_intron_range_reference_spans_long_tail():
    name = "Xerocrassa_montserratensis.train.intron.tbl"
    lens = [len(ln.split("\t")[1].strip()) for ln in open(REF / name)]
    lo, hi = intron_range(lens)
    assert lo == 24.75
    # the p99.9 max spans the skewed long tail -> well above the legacy mean+3sd
    # (~25394) that clipped ~1.8% of real introns, and under the 100 kb safety cap
    assert hi > 25394
    assert hi <= 100_000
    assert sum(1 for x in lens if x > hi) <= len(lens) // 1000 + 1
