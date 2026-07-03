import math

import pytest

from geneid_train.stats.sites import (
    frequency,
    info_content,
    log_ratio,
    log_ratio_zero_order,
    mask_invariant_dinuc,
    position_matrix,
    read_matrix,
    select_window,
    submatrix,
)

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


def test_log_ratio_zero_order_natural_masking():
    # freq==0 -> -9999 (never observed); freq==1 -> 0 (invariant); else log ratio
    site = {(4, "A"): 1.0, (4, "C"): 0.0, (7, "G"): 0.4}
    bg = {(4, "A"): 0.3, (4, "C"): 0.2, (7, "G"): 0.2}
    out = log_ratio_zero_order(site, bg)
    assert out[(4, "A")] == 0.0
    assert out[(4, "C")] == -9999.0
    assert out[(7, "G")] == pytest.approx(math.log(0.4 / 0.2))


@pytest.mark.integration
@pytest.mark.skipif(REF is None, reason="set GENEID_TRAIN_REFDIR to the train_geneid dir")
def test_start_order0_masking_matches_reference_shape():
    """Order-0 (Start) profile: no explicit masking step -- the ATG anchor falls
    out naturally from exact 0/1 raw frequencies (no pseudocounts). The reference
    log-ratio's 0 / -9999 cells are exactly where our raw frequency is 1 / 0,
    independent of the background model used."""
    seqs = [ln.split()[1].strip() for ln in open(REF / f"{SP}.canonical.start.tbl")]
    site = position_matrix(seqs, order=0, pcount=0.0)
    ref = read_matrix(REF / f"{SP}.canonical.start-log.order-0-matrix")
    zero_cells = {k for k, v in ref.items() if v == 0.0}
    masked_cells = {k for k, v in ref.items() if v == -9999.0}
    assert zero_cells == {k for k, v in site.items() if v == 1.0}
    assert masked_cells == {k for k, v in site.items() if v == 0.0}


def test_mask_invariant_dinuc_donor():
    # profile positions 3/4/5 carry the invariant GT
    m = {(3, "AG"): 1.5, (3, "AA"): 1.5, (4, "GT"): 2.0, (4, "AA"): 2.0,
         (5, "TA"): 0.7, (5, "AA"): 0.7, (6, "AA"): -0.3}
    out = mask_invariant_dinuc(m, st=3, nd=4, rd=5, anchor="GT")
    assert out[(3, "AG")] == 0.0 and out[(3, "AA")] == -9999.0  # ends in G -> 0
    assert out[(4, "GT")] == 0.0 and out[(4, "AA")] == -9999.0  # == GT -> 0
    assert out[(5, "TA")] == 0.7 and out[(5, "AA")] == -9999.0  # starts T -> keep
    assert out[(6, "AA")] == -0.3  # untouched


# ---- boundary selection (info content -> profile window + anchor) -----------


def read_freq(path):
    """Parse a legacy ``.freq`` file (``base pos count freq``) -> {(pos, base): freq}."""
    out = {}
    for ln in open(path):
        p = ln.split()
        if len(p) >= 4:
            out[(int(p[1]), p[0])] = float(p[3])
    return out


def read_info(path):
    """Parse a legacy tempinfolog file (``pos bits``) -> {pos: bits}."""
    out = {}
    for ln in open(path):
        p = ln.split()
        if len(p) == 2 and p[0].lstrip("-").isdigit():
            out[int(p[0])] = float(p[1])
    return out


def test_frequency_per_position_denominator():
    # N masks only its own column: pos1 sees A twice + one N -> total 2
    f = frequency(["AC", "AG", "NC"])
    assert f[(1, "A")] == 1.0  # 2/2 valid
    assert f[(2, "C")] == pytest.approx(2 / 3)
    assert f[(2, "G")] == pytest.approx(1 / 3)


def test_info_content_relative_entropy_bits():
    # one position, site all-A vs uniform background -> log2(4) = 2 bits
    site = {(1, b): (1.0 if b == "A" else 0.0) for b in "ACGT"}
    bg = {(1, b): 0.25 for b in "ACGT"}
    info = info_content(site, bg)
    assert info[1] == pytest.approx(2.0)


def test_select_window_donor_from_info():
    # peaks at 30..37; seed 29/31. -> start 29, end 37-order=36, anchor 3/4/5
    info = {p: (1.0 if 30 <= p <= 37 else 0.01) for p in range(25, 39)}
    w = select_window(info, site="donor", order=1)
    assert (w.start, w.end) == (29, 36)
    assert (w.st, w.nd, w.rd) == (3, 4, 5)


def test_select_window_acceptor_from_info():
    # everything above 0.04 from pos 2..32 -> start 2, end 31, anchor 27/28/29
    info = {p: 0.5 for p in range(2, 33)}
    w = select_window(info, site="acceptor", order=1)
    assert (w.start, w.end) == (2, 31)
    assert (w.st, w.nd, w.rd) == (27, 28, 29)


@pytest.mark.integration
@pytest.mark.skipif(REF is None, reason="set GENEID_TRAIN_REFDIR to the train_geneid dir")
def test_frequency_matches_reference():
    seqs = [ln.split("\t")[1].strip() for ln in open(REF / f"{SP}.canonical.donor.tbl")]
    ours = frequency(seqs)
    ref = read_freq(REF / f"{SP}.canonical.donor.freq")
    worst = max(abs(ours[k] - ref[k]) for k in ref)
    assert worst < 1e-3, f"max deviation {worst}"


@pytest.mark.integration
@pytest.mark.skipif(REF is None, reason="set GENEID_TRAIN_REFDIR to the train_geneid dir")
def test_info_content_matches_reference():
    site = read_freq(REF / f"{SP}.canonical.donor.freq")
    bg = read_freq(REF / f"{SP}_background.info.freq")
    ours = info_content(site, bg)
    ref = read_info(REF / f"{SP}.canonical.donor-{SP}_background.info")
    shared = set(ours) & set(ref)
    assert len(shared) > 5
    worst = max(abs(ours[k] - ref[k]) for k in shared)
    assert worst < 1e-4, f"max deviation {worst}"


@pytest.mark.integration
@pytest.mark.skipif(REF is None, reason="set GENEID_TRAIN_REFDIR to the train_geneid dir")
def test_select_window_matches_reference_donor_and_acceptor():
    dinfo = read_info(REF / f"{SP}.canonical.donor-{SP}_background.info")
    dw = select_window(dinfo, site="donor", order=1)
    assert (dw.start, dw.end, dw.st, dw.nd, dw.rd) == (29, 36, 3, 4, 5)
    ainfo = read_info(REF / f"{SP}.canonical.acceptor-{SP}_background.info")
    aw = select_window(ainfo, site="acceptor", order=1)
    assert (aw.start, aw.end, aw.st, aw.nd, aw.rd) == (2, 31, 27, 28, 29)


@pytest.mark.integration
@pytest.mark.skipif(REF is None, reason="set GENEID_TRAIN_REFDIR to the train_geneid dir")
def test_full_donor_profile_matches_param():
    """End-to-end: donor site sequences -> the trained param's Donor_profile block,
    with the window and GT-anchor columns selected automatically from info content."""
    seqs = [ln.split("\t")[1].strip() for ln in open(REF / f"{SP}.canonical.donor.tbl")]
    site = position_matrix(seqs, order=1)
    bg = read_matrix(REF / f"{SP}_background.info.di-matrix")
    info = info_content(
        read_freq(REF / f"{SP}.canonical.donor.freq"),
        read_freq(REF / f"{SP}_background.info.freq"),
    )
    w = select_window(info, site="donor", order=1)
    prof = mask_invariant_dinuc(
        submatrix(log_ratio(site, bg), w.start, w.end), w.st, w.nd, w.rd, "GT"
    )
    ref = read_matrix(REF / f"{SP}.canonical.donor-log-info.di-matrix")  # == param Donor_profile
    assert set(prof) == set(ref)
    worst = max(abs(prof[k] - ref[k]) for k in ref)
    assert worst < 1e-3, f"max deviation {worst}"


@pytest.mark.integration
@pytest.mark.skipif(REF is None, reason="set GENEID_TRAIN_REFDIR to the train_geneid dir")
def test_full_acceptor_profile_matches_param():
    """Same pipeline as the donor test, with the AG anchor and window selected
    automatically from acceptor info content."""
    seqs = [ln.split("\t")[1].strip() for ln in open(REF / f"{SP}.canonical.acceptor.tbl")]
    site = position_matrix(seqs, order=1)
    bg = read_matrix(REF / f"{SP}_background.info.di-matrix")
    info = info_content(
        read_freq(REF / f"{SP}.canonical.acceptor.freq"),
        read_freq(REF / f"{SP}_background.info.freq"),
    )
    w = select_window(info, site="acceptor", order=1)
    prof = mask_invariant_dinuc(
        submatrix(log_ratio(site, bg), w.start, w.end), w.st, w.nd, w.rd, "AG"
    )
    ref = read_matrix(REF / f"{SP}.canonical.acceptor-log-info.di-matrix")
    assert set(prof) == set(ref)
    worst = max(abs(prof[k] - ref[k]) for k in ref)
    assert worst < 1e-3, f"max deviation {worst}"
