"""Tests for emitting the U12 profile trio into a geneid .param (param/u12.py)."""

import os
import subprocess
import tempfile
from pathlib import Path

import pytest

from geneid_train.core.param import Param, Profile
from geneid_train.param.u12 import (
    DONOR_SIDE,
    donor_loglik,
    load_bundled_u12,
    load_u12_donor_model,
    overlay_profile,
    partition_sections,
    register_shift,
)

HERE = Path(__file__).resolve().parent
REF_U12 = HERE.parents[1] / "param" / "human3isoU12.param"
GENEID = os.environ.get("GENEID_BIN") or "/Users/talioto/repositories/geneid_fresh/bin/geneid"


def test_insert_text_before_places_section():
    p = Param.from_text("Acceptor_profile\n1 A 0\nDonor_profile\n1 A 0\n")
    p.insert_text_before("Acceptor_profile", "U12_Branch_point_profile\n1 AA 0\n")
    kws = p.keywords()
    assert kws == ["U12_Branch_point_profile", "Acceptor_profile", "Donor_profile"]
    # round-trips as text with the new section spliced in front
    assert p.to_text().startswith("U12_Branch_point_profile\n1 AA 0\n")


def test_register_shift_recovers_known_offset():
    # a reference PWM at positions 3..5 and an IAOD PWM (same shape) at positions
    # 1..3 must register with shift +2 (iaod pos j -> ref pos j+2)
    ref = {(3, "AA"): 2.0, (4, "AA"): -1.0, (5, "AA"): 3.0}
    iaod = {(1, "AA"): 2.0, (2, "AA"): -1.0, (3, "AA"): 3.0}
    assert register_shift(iaod, ref) == 2


def test_overlay_profile_keeps_geometry_and_neutralises_gaps():
    # reference geometry: dimension 4, order 1 (dinucleotides), offset 1
    ref = Profile("Donor_profile", [4.0, 1.0, -7.0, 1.0, 0.0, 1.0], [])
    iaod = {(1, "AC"): 5.5, (1, "AA"): -2.0}  # only one position of signal
    lines = overlay_profile(ref, iaod, shift=1)  # iaod pos1 -> ref pos2
    assert lines[0] == "4 1 -7 1 0 1"
    assert lines[1].startswith("#")
    rows = {(int(p), o): float(v) for p, o, v in (ln.split() for ln in lines[2:])}
    # every position 1..4 x every dinucleotide present (complete matrix)
    assert len(rows) == 4 * 16
    # the IAOD value landed at reference position 2
    assert rows[(2, "AC")] == 5.5
    assert rows[(2, "AA")] == -2.0
    # positions outside the IAOD window are neutral
    assert rows[(1, "AC")] == 0.0
    assert rows[(4, "GT")] == 0.0


def test_partition_routes_by_profile_name():
    text = (
        "U12_Branch_point_profile\n12 9 2 2\n1 AAA 0\n"
        "U12gtag_Acceptor_profile\n13 9 1 2\n1 AAA 0\n"
        "U12gtag_Donor_profile\n11 1 0 2\n1 AAA 0\n"
    )
    sec = partition_sections(text)
    assert "U12_Branch_point_profile" in sec.acceptor_side
    assert "U12gtag_Acceptor_profile" in sec.acceptor_side
    assert "U12gtag_Donor_profile" in sec.donor_side
    assert "U12gtag_Donor_profile" not in sec.acceptor_side
    assert "U12gtag_Donor_profile" in DONOR_SIDE  # routing set drives partition


def test_assemble_emits_u12_score_thresholds_before_exon_weights():
    from geneid_train.param.assemble import assemble_param

    stub = ["1 A 0"]
    txt = assemble_param(
        species="X", start_profile=stub, acceptor_profile=stub, donor_profile=stub,
        markov_order=1, markov_initial=stub, markov_transition=stub,
        intron_range="40:80000", intergenic_range="200:Infinity",
        u12=load_bundled_u12(),  # defaults: splice 9, exon 8
    )
    lines = txt.splitlines()
    assert lines[lines.index("U12_Splice_Score_Threshold") + 1] == "9"
    assert lines[lines.index("U12_Exon_Score_Threshold") + 1] == "8"
    # geneid reads these gates in the optional-scalar block before Exon_weights
    assert txt.index("U12_Splice_Score_Threshold") < txt.index("Exon_weights")


def test_no_u12_thresholds_when_u12_absent():
    from geneid_train.param.assemble import assemble_param

    stub = ["1 A 0"]
    txt = assemble_param(
        species="X", start_profile=stub, acceptor_profile=stub, donor_profile=stub,
        markov_order=1, markov_initial=stub, markov_transition=stub,
        intron_range="40:80000", intergenic_range="200:Infinity",
    )
    assert "U12_Splice_Score_Threshold" not in txt


def test_bundled_u12_loads_and_has_full_trio():
    sec = load_bundled_u12()
    acc_names = ("U12_Branch_point_profile", "U12gtag_Acceptor_profile", "U12atac_Acceptor_profile")
    for name in acc_names:
        assert name in sec.acceptor_side
    for name in ("U12gtag_Donor_profile", "U12atac_Donor_profile"):
        assert name in sec.donor_side
    # the branch header carries the four distance knobs geneid needs to tune
    branch_hdr = [
        ln for ln in sec.acceptor_side.splitlines()
    ]
    idx = branch_hdr.index("U12_Branch_point_profile")
    fields = branch_hdr[idx + 1].split()
    assert len(fields) == 10  # len offset cutoff order a b acc_context min_dist opt_dist pen_scale


def test_donor_loglik_sums_position_logs():
    import math

    freq = {(1, "A"): 0.5, (2, "C"): 0.25}
    # exact match: log(0.5+eps) + log(0.25+eps)
    got = donor_loglik("AC", freq, 2, pseudo=0.0)
    assert abs(got - (math.log(0.5) + math.log(0.25))) < 1e-9


def test_bundled_u12_donor_model_is_gtatcctt():
    freq, length, floor = load_u12_donor_model()
    assert length >= 8
    # the U12 5' consensus GTATCCTT must be the per-position argmax of the PWM
    consensus = "".join(
        max("ACGT", key=lambda b: freq.get((pos, b), 0.0)) for pos in range(1, 9)
    )
    assert consensus == "GTATCCTT"
    assert floor < 0  # a log-likelihood floor


# ---- end-to-end: geneid actually enables U12 from the bundled profiles -------


@pytest.mark.integration
@pytest.mark.skipif(
    not Path(GENEID).exists() or not REF_U12.exists(),
    reason="needs a geneid binary (GENEID_BIN) and a reference param",
)
def test_geneid_enables_u12_from_bundle():
    # splice the bundled U12 trio into a plain single-isochore param and confirm
    # geneid loads the branch profile and scores U12 sites (i.e. U12 is switched on)
    # use a non-U12 single-isochore param as the host
    host = Path("/Users/talioto/repositories/geneid_fresh/param/dros.param")
    if not host.exists():
        pytest.skip("no host param")
    sec = load_bundled_u12()
    p = Param.read(host)
    p.insert_text_before("Acceptor_profile", sec.acceptor_side)
    p.insert_text_before("Donor_profile", sec.donor_side)
    with tempfile.TemporaryDirectory() as td:
        pm = Path(td) / "u12.param"
        p.write(pm)
        fa = Path(td) / "seq.fa"
        import random

        random.seed(0)
        s = "".join(random.choice("ACGT") for _ in range(4000))
        fa.write_text(">s\n" + "\n".join(s[i : i + 60] for i in range(0, len(s), 60)) + "\n")
        out = subprocess.run(
            [GENEID, "-v", "-U", "-P", str(pm), str(fa)],
            capture_output=True, text=True,
        )
    combined = out.stdout + out.stderr
    assert "Wrong format" not in combined
    assert "BranchPoint" in combined  # branch profile loaded
    assert "U12gtag" in combined  # U12gtag sites scored -> class active
