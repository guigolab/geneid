import pytest

from geneid_train.stats.genemodel import format_range, intron_range

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
