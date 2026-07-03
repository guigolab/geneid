import glob
import os

import pytest

from geneid_train.prepare.u12 import (
    acceptor_windows,
    by_subtype,
    consensus,
    donor_windows,
    load_u12_introns,
    parse_iaod_fasta,
    train_u12_profile,
)
from geneid_train.stats.sites import MASK, position_matrix

_FIXTURE = (
    ">Homo sapiens|1|+|100|210|110|2|5\n"
    "GTATCCTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACAG\n"
    ">Homo sapiens|1|-|300|410|110|1|6\n"
    "GTATCCTTTTGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTTTTCAG\n"
    ">Zea mays|3|+|1|60|59|1|1\n"
    "ATATCCTTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTATAC\n"
    ">Odd sp|4|+|1|20|19|1|1\n"
    "GCACGTACGTACGTACGTGG\n"
)


def _write(tmp_path):
    p = tmp_path / "X_U12.fasta"
    p.write_text(_FIXTURE)
    return p


def test_parse_iaod_fasta_headers_and_seq(tmp_path):
    introns = parse_iaod_fasta(_write(tmp_path))
    assert len(introns) == 4
    a = introns[0]
    assert a.species == "Homo sapiens" and a.chrom == "1" and a.strand == "+"
    assert a.start == 100 and a.end == 210
    assert a.seq.startswith("GTATCCTT") and a.seq.endswith("AG")


def test_subtype_classification():
    from geneid_train.prepare.u12 import U12Intron

    assert U12Intron("GT" + "A" * 10 + "AG").subtype == "gtag"
    assert U12Intron("AT" + "A" * 10 + "AC").subtype == "atac"
    assert U12Intron("GC" + "A" * 10 + "GG").subtype == "other"


def test_by_subtype_and_windows(tmp_path):
    introns = load_u12_introns([_write(tmp_path)])
    groups = by_subtype(introns)
    assert len(groups["gtag"]) == 2  # two GT-AG
    assert len(groups["atac"]) == 1
    assert len(groups["other"]) == 1
    dw = donor_windows(groups["gtag"], 8)
    assert dw == ["GTATCCTT", "GTATCCTT"]
    aw = acceptor_windows(groups["gtag"], 3)
    assert all(w.endswith("AG") for w in aw)


def test_consensus_first_and_last():
    assert consensus(["GTATCCTTAA", "GTATCCTTAA"], 8) == "GTATCCTT"
    assert consensus(["CCCCAG", "TTTTAG"], 2, from_end=True) == "AG"


def test_train_u12_profile_is_unclamped():
    # a U12 donor-like set; the terminal dinucleotide must NOT be masked to -9999
    seqs = ["GTATCCTTAC", "GTATCCTTAG", "GTGTCCTTAC", "GTATCCTTAA"]
    bg = position_matrix(["ACGT" * 10, "TGCA" * 10, "GATC" * 10], order=1)
    prof = train_u12_profile(seqs, bg, order=1, start=1, end=6)
    assert prof  # non-empty
    assert all(v != MASK for v in prof.values())  # nothing clamped (unlike U2)
    # positions renumbered to 1..(end-start+1)
    assert max(p for p, _ in prof) == 6


# ---- real IAOD data (opt-in) ------------------------------------------------

U12DIR = os.environ.get("GENEID_TRAIN_U12DIR")


@pytest.mark.integration
@pytest.mark.skipif(not U12DIR, reason="set GENEID_TRAIN_U12DIR to a dir of *_U12.fasta files")
def test_iaod_pooled_consensus_is_u12():
    introns = load_u12_introns(sorted(glob.glob(f"{U12DIR}/*_U12.fasta")))
    assert len(introns) > 5000
    groups = by_subtype(introns)
    assert len(groups["gtag"]) > 3000 and len(groups["atac"]) > 500
    # the canonical U12 5' splice sites must emerge from the pooled set
    assert consensus(donor_windows(groups["gtag"], 8), 8) == "GTATCCTT"
    assert consensus(donor_windows(groups["atac"], 8), 8) == "ATATCCTT"
