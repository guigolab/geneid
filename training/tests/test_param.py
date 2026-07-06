"""Param read/write round-trip and typed-accessor tests.

The round-trip is the Phase 1 correctness oracle: parsing then re-serializing a
geneid ``.param`` must be byte-identical, since the frozen format is a contract
with the compiled ``geneid`` binary.
"""

from pathlib import Path

import pytest

from geneid_train.core.param import Param

from .conftest import real_param_files


def test_roundtrip_byte_identical_mini(mini_param_path: Path):
    text = mini_param_path.read_text()
    assert Param.from_text(text).to_text() == text


def test_scalar_and_vector_accessors(mini_param_path: Path):
    p = Param.read(mini_param_path)
    assert p.scalar("NO_SCORE") == "0"
    assert p.num_isochores == 1
    assert p.vector("Absolute_cutoff_exons") == ["-15", "-15", "-15", "-15"]
    assert p.scalar("maximum_number_of_donors_per_acceptor_site") == "7"


def test_u12_detection_and_profiles(mini_param_path: Path):
    p = Param.read(mini_param_path)
    assert p.has("U12_Splice_Score_Threshold")
    assert p.scalar("U12_Exon_weight") == "-3"
    names = p.profile_names()
    assert names == ["Start_profile", "Donor_profile", "Stop_profile"]
    start = p.profile("Start_profile")
    assert start.length == 3
    assert start.order == 0
    assert len(start.rows) == 12
    assert start.rows[0] == (1, "A", 0.1)


def test_set_scalar_is_surgical(mini_param_path: Path):
    text = mini_param_path.read_text()
    p = Param.from_text(text)
    p.set_scalar("NO_SCORE", -5)
    out = p.to_text()
    # exactly one line changed
    before = text.splitlines()
    after = out.splitlines()
    assert len(before) == len(after)
    diff = [(a, b) for a, b in zip(before, after, strict=True) if a != b]
    assert diff == [("0", "-5")]
    # and the change parses back
    assert Param.from_text(out).scalar("NO_SCORE") == "-5"


def test_missing_section_raises(mini_param_path: Path):
    p = Param.read(mini_param_path)
    with pytest.raises(KeyError):
        p.scalar("does_not_exist")


@pytest.mark.integration
@pytest.mark.parametrize("path", real_param_files(), ids=lambda p: p.name)
def test_roundtrip_byte_identical_real_params(path: Path):
    text = path.read_text()
    assert Param.from_text(text).to_text() == text
