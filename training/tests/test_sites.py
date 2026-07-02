import math

import pytest

from geneid_train.stats.sites import log_ratio, position_matrix, read_matrix

from .conftest import ref_dir


def test_position_matrix_order1_conditionals():
    m = position_matrix(["ACGT", "ACGA"], order=1, pcount=0.0)  # unsmoothed for exact fractions
    assert m[(1, "AC")] == 1.0  # A always followed by C at pos 1
    assert m[(2, "CG")] == 1.0
    assert m[(3, "GT")] == 0.5  # G -> T or A at pos 3
    assert m[(3, "GA")] == 0.5
    # positions run 1..L-order
    assert max(p for p, _ in m) == 3


def test_position_matrix_pseudocounts():
    # single AC: P(AC) = (0.25+1)/(1.0+1) = 0.625; unobserved P(AA) = 0.25/2 = 0.125
    m = position_matrix(["AC"], order=1, pcount=0.25)
    assert m[(1, "AC")] == 0.625
    assert m[(1, "AA")] == 0.125
    # complete matrix: all 16 dinucleotides present at position 1
    assert len([k for k in m if k[0] == 1]) == 16


def test_position_matrix_skips_whole_window_with_non_acgt():
    m = position_matrix(["ANGT", "ACGT"], order=1, pcount=0.0)
    # ANGT dropped entirely -> only ACGT contributes
    assert m[(1, "AC")] == 1.0
    assert m[(1, "AA")] == 0.0


def test_log_ratio_natural_log():
    out = log_ratio({(1, "AA"): 0.5}, {(1, "AA"): 0.25})
    assert out[(1, "AA")] == pytest.approx(math.log(2.0))


def test_log_ratio_broadcasts_position_independent_background():
    site = {(2, "AA"): 0.4}
    bg = {(1, "AA"): 0.2}  # position-independent background
    out = log_ratio(site, bg)
    assert out[(2, "AA")] == pytest.approx(math.log(2.0))


# ---- validation against the real xgXerMont reference run -------------------

REF = ref_dir()
SP = "Xerocrassa_montserratensis"


@pytest.mark.integration
@pytest.mark.skipif(REF is None, reason="set GENEID_TRAIN_REFDIR to the train_geneid dir")
def test_donor_position_matrix_matches_reference():
    seqs = [ln.split("\t")[1].strip() for ln in open(REF / f"{SP}.canonical.donor.tbl")]
    ours = position_matrix(seqs, order=1)
    ref = read_matrix(REF / f"{SP}.canonical.donor.di-matrix")
    shared = set(ours) & set(ref)
    assert len(shared) > 100
    worst = max(abs(ours[k] - ref[k]) for k in shared)
    assert worst < 1e-3, f"max deviation {worst}"


@pytest.mark.integration
@pytest.mark.skipif(REF is None, reason="set GENEID_TRAIN_REFDIR to the train_geneid dir")
def test_donor_log_ratio_matches_reference():
    site = read_matrix(REF / f"{SP}.canonical.donor.di-matrix")
    bg = read_matrix(REF / f"{SP}_background.info.di-matrix")
    ours = log_ratio(site, bg)
    ref = read_matrix(REF / f"{SP}.canonical.donor-log.di-matrix")
    shared = set(ours) & set(ref)
    assert len(shared) > 100
    worst = max(abs(ours[k] - ref[k]) for k in shared)
    assert worst < 1e-3, f"max deviation {worst}"
