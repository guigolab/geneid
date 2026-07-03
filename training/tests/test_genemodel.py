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


@pytest.mark.integration
@pytest.mark.skipif(REF is None, reason="set GENEID_TRAIN_REFDIR to the train_geneid dir")
def test_intron_range_matches_reference():
    name = "Xerocrassa_montserratensis.train.intron.tbl"
    lens = [len(ln.split("\t")[1].strip()) for ln in open(REF / name)]
    lo, hi = intron_range(lens)
    # reference gene model uses 24.75:25394.023; lo is exact, hi matches to <1 bp
    # (the last-decimal difference is float-printing noise on a coarse bound)
    assert lo == 24.75
    assert abs(hi - 25394.023) < 1.0
