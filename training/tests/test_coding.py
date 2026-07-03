import math

import pytest

from geneid_train.stats.coding import (
    choose_orders,
    coding_log_ratio,
    derive_coding_potential,
    format_markov_matrix,
    initial_model,
    transition_model,
)

from .conftest import ref_dir


def test_transition_model_conditional_probabilities():
    # ACGTACGT (order 1, single frame): A->C, C->G, G->T, T->A observed
    m = transition_model(["ACGTACGT"], order=1, pseudo=0.0, nframes=1)
    assert m[("AC", 0)] == 1.0
    assert m[("CG", 0)] == 1.0
    assert m[("GT", 0)] == 1.0
    # unobserved prefix falls to 0 (no counts, no pseudo)
    assert m[("AA", 0)] == 0.0


def test_transition_model_frames_sum_to_one():
    seqs = ["ACGTACGTACGTACGTACGT"]
    m = transition_model(seqs, order=1, pseudo=0.25, nframes=3)
    for fr in range(3):
        for pre in "ACGT":
            s = sum(m[(pre + nt, fr)] for nt in "ACGT")
            assert s == pytest.approx(1.0)


def test_initial_model_is_a_frequency_distribution():
    seqs = ["ACGTACGTACGTACGT"]
    # without pseudocounts the per-frame oligo frequencies form a distribution
    m = initial_model(seqs, order=1, pseudo=0.0, nframes=1)
    total = sum(v for (o, fr), v in m.items() if fr == 0)
    assert total == pytest.approx(1.0)


def test_coding_log_ratio_natural_log():
    coding = {("AAAAA", 0): 0.5, ("AAAAA", 1): 0.4}
    intron = {("AAAAA", 0): 0.25}
    out = coding_log_ratio(coding, intron)
    assert out[("AAAAA", 0)] == pytest.approx(math.log(2.0))
    assert out[("AAAAA", 1)] == pytest.approx(math.log(1.6))


def test_choose_orders():
    assert choose_orders(500_000, 200_000) == (4, 5)
    assert choose_orders(1_300_000, 30_000_000) == (4, 5)
    assert choose_orders(1_000, 1_000) == (3, 4)


def test_format_markov_matrix_indexing():
    model = {("AAAA", 0): 1.234, ("AAAA", 1): -0.5, ("AAAC", 0): 0.1}
    lines = format_markov_matrix(model)
    assert lines[0] == "AAAA 0 0 1.234"
    assert lines[1] == "AAAA 0 1 -0.500"
    assert lines[2] == "AAAC 1 0 0.100"


# ---- validation against the trained xgXerMont param -------------------------

REF = ref_dir()


def _read_param_markov(path, header):
    """Parse a Markov matrix block (oligo index frame value) from a param file."""
    out = {}
    active = False
    for ln in open(path):
        s = ln.strip()
        if s == header:
            active = True
            continue
        if active:
            p = s.split()
            if len(p) == 4 and p[0].isalpha():
                out[(p[0], int(p[2]))] = float(p[3])
            elif out:  # first non-data line after the block ends it
                break
    return out


@pytest.mark.integration
@pytest.mark.skipif(REF is None, reason="set GENEID_TRAIN_REFDIR to the train_geneid dir")
def test_coding_potential_matches_param():
    def read_seqs(name):
        return [ln.split("\t")[1].strip() for ln in open(REF / name)]

    cds = read_seqs("Xerocrassa_montserratensis.train.cds.tbl")
    introns = read_seqs("Xerocrassa_montserratensis.train.intron.tbl")
    init, trans, cbases, nbases = derive_coding_potential(cds, introns, transition_order=5)

    assert choose_orders(cbases, nbases) == (4, 5)  # this training set warrants 5/4

    param = REF / "Xerocrassa_montserratensis.geneid.param"
    ref_init = _read_param_markov(param, "Markov_Initial_probability_matrix")
    ref_trans = _read_param_markov(param, "Markov_Transition_probability_matrix")
    assert len(ref_init) == 3072 and len(ref_trans) == 12288

    # both matrices reproduce the param exactly at its 3-decimal precision
    worst_init = max(abs(round(init[k], 3) - ref_init[k]) for k in ref_init)
    worst_trans = max(abs(round(trans[k], 3) - ref_trans[k]) for k in ref_trans)
    assert worst_init == 0.0, f"initial max dev {worst_init}"
    assert worst_trans == 0.0, f"transition max dev {worst_trans}"
