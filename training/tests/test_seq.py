from geneid_train.core.seq import (
    STANDARD,
    GeneticCode,
    has_internal_stop,
    is_complete_cds,
    revcomp,
)


def test_revcomp():
    assert revcomp("ATGC") == "GCAT"
    assert revcomp("aaccggtt") == "aaccggtt"[::-1].translate(str.maketrans("acgt", "tgca"))


def test_translate_standard():
    # ATG AAA TTT TAA -> M K F *
    assert STANDARD.translate("ATGAAATTTTAA") == "MKF*"


def test_is_complete_cds():
    assert is_complete_cds("ATGAAAAAATTTTAA")  # M K K F *
    assert not is_complete_cds("ATGAAATTT")  # no terminal stop
    assert not is_complete_cds("AAAAAATTTTAA")  # no start
    assert not is_complete_cds("ATGTAAAAATAA")  # internal stop
    assert not is_complete_cds("ATGAAATT")  # not multiple of 3


def test_has_internal_stop_ignores_terminal():
    assert not has_internal_stop("ATGAAATAA")
    assert has_internal_stop("ATGTAAAAATAA")


def test_alternative_genetic_code_ciliate():
    # ciliate: TAA/TAG code for Q, only TGA is a stop
    codons = dict(STANDARD.codons)
    codons["TAA"] = "Q"
    codons["TAG"] = "Q"
    ciliate = GeneticCode(codons=codons)
    assert ciliate.stops == frozenset({"TGA"})
    # a CDS with TAA internally is complete under the ciliate code, not the standard one
    cds = "ATGTAAAAATGA"  # M Q K *  (ciliate)
    assert is_complete_cds(cds, ciliate)
    assert not is_complete_cds(cds, STANDARD)
